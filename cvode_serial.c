/* Aug 2010, Timothy Bourke (INRIA) */

#include <caml/mlvalues.h>
#include <caml/bigarray.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h> /* definition of type realtype */

/* linear solvers */
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <cvode/cvode_lapack.h>
#endif

#include <stdio.h>

#define VARIANT_LMM_ADAMS 0
#define VARIANT_LMM_BDF   1

#define RECORD_BANDRANGE_MUPPER 0
#define RECORD_BANDRANGE_MLOWER 1

#define RECORD_SPRANGE_PRETYPE 0
#define RECORD_SPRANGE_MAXL    1

/* untagged: */
#define VARIANT_LINEAR_SOLVER_DENSE	    0
#define VARIANT_LINEAR_SOLVER_LAPACKDENSE   1
#define VARIANT_LINEAR_SOLVER_DIAG	    2
/* tagged: */
#define VARIANT_LINEAR_SOLVER_BAND	    0
#define VARIANT_LINEAR_SOLVER_LAPACKBAND    1
#define VARIANT_LINEAR_SOLVER_SPGMR	    2
#define VARIANT_LINEAR_SOLVER_SPBCG	    3
#define VARIANT_LINEAR_SOLVER_SPTFQMR	    4

#define VARIANT_SOLVER_RESULT_CONTINUE		0
#define VARIANT_SOLVER_RESULT_ROOTSFOUND	1
#define VARIANT_SOLVER_RESULT_STOPTIMEREACHED	2

#define RECORD_INTEGRATOR_STATS_STEPS			0
#define RECORD_INTEGRATOR_STATS_RHS_EVALS		1
#define RECORD_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS	2
#define RECORD_INTEGRATOR_STATS_ERROR_TEST_FAILURES	3
#define RECORD_INTEGRATOR_STATS_LAST_INTERNAL_ORDER	4
#define RECORD_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER	5
#define RECORD_INTEGRATOR_STATS_INITIAL_STEP_SIZE	6
#define RECORD_INTEGRATOR_STATS_LAST_STEP_SIZE		7
#define RECORD_INTEGRATOR_STATS_NEXT_STEP_SIZE		8
#define RECORD_INTEGRATOR_STATS_INTERNAL_TIME		9

#define RECORD_ERROR_DETAILS_ERROR_CODE	    0
#define RECORD_ERROR_DETAILS_MODULE_NAME    1
#define RECORD_ERROR_DETAILS_FUNCTION_NAME  2
#define RECORD_ERROR_DETAILS_ERROR_MESSAGE  3

/* XXX:
 * - realtype must equal double
 * - we assume (in BIGARRAY_INT) that an int is 32-bits
 *   (this should be configured per platform.)
 */

/*
 * TODO:
 * - call Gc.full_major () from the f and roots routines to see if we get
 *   any segmentation fault problems.
 */

/* TODO:
 * Possible extensions:
 * - add a closure for returning the Jacobian.
 * - allow T0 to be configured
 *
 */


#define T0    RCONST(0.0)      /* initial time */
#define BIGARRAY_FLOAT (CAML_BA_FLOAT64 | CAML_BA_C_LAYOUT)
#define BIGARRAY_INT (CAML_BA_INT32 | CAML_BA_C_LAYOUT)
#define MAX_ERRMSG_LEN 256

static void check_flag(const char *call, int flag, void *to_free);

// TODO: Is there any risk that the Ocaml GC will try to free the two
//	 closures? Do we have to somehow record that we're using them,
//	 and then release them again in the free routine?
//	 SEE: ml_cvode_data_alloc and finalize
struct ml_cvode_data {
    void *cvode_mem;
    long int neq;
    intnat num_roots;
    value *closure_f;
    value *closure_roots;
    value *closure_errh;
    FILE *err_file;
};

typedef struct ml_cvode_data* ml_cvode_data_p;

static void finalize(value vdata)
{
    ml_cvode_data_p data = (ml_cvode_data_p)Data_custom_val(vdata);

    // TODO:
    // The Ocaml Manual (18.9.1) says:
    // ``Note: the finalize, compare, hash, serialize and deserialize
    // functions attached to custom block descriptors must never trigger a
    // garbage collection. Within these functions, do not call any of the
    // Caml allocation functions, and do not perform a callback into Caml
    // code. Do not use CAMLparam to register the parameters to these
    // functions, and do not use CAMLreturn to return the result.''
    //
    // But, obviously, we're calling two caml functions. We need to find out
    // if this is ok.
    caml_remove_generational_global_root(data->closure_f);
    data->closure_f = NULL;

    caml_remove_generational_global_root(data->closure_roots);
    data->closure_roots = NULL;

    if (data->closure_errh != NULL) {
	caml_remove_generational_global_root(data->closure_errh);
	data->closure_errh = NULL;
    }

    if (data->cvode_mem != NULL) {
	CVodeFree(&(data->cvode_mem));
    }

    if (data->err_file != NULL) {
	fclose(data->err_file);
	data->err_file = NULL;
    }
}

// TODO:
// The Ocaml Manual (18.9.3) says:
// ``The contents of custom blocks are not scanned by the garbage collector,
// and must therefore not contain any pointer inside the Caml heap. In other
// terms, never store a Caml value in a custom block, and do not use Field,
// Store_field nor caml_modify to access the data part of a custom block.
// Conversely, any C data structure (not containing heap pointers) can be
// stored in a custom block.''
//
// But, obviously, we're storing two closure values in the struct. We need
// to find out if and when this is ok.
//
static value ml_cvode_data_alloc(mlsize_t approx_size)
{
    return caml_alloc_final(sizeof(struct ml_cvode_data), &finalize,
			    approx_size, 10);
}

static ml_cvode_data_p cvode_data_from_ml(value vdata)
{
    ml_cvode_data_p data = (ml_cvode_data_p)Data_custom_val(vdata);
    if (data->cvode_mem == NULL) {
	caml_failwith("This session has been freed");
    }

    return data;
}

static int check_exception(value r)
{
    CAMLparam0();
    static value *recoverable_failure = NULL;

    if (!Is_exception_result(r)) return 0;

    if (recoverable_failure == NULL) {
	recoverable_failure =
	    caml_named_value("cvode_RecoverableFailure");
    }

    value exn = Extract_exception(r);
    CAMLreturn((Field(exn, 0) == *recoverable_failure) ? 1 : -1);
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();

    value *closure_f = ((ml_cvode_data_p)user_data)->closure_f;

    intnat y_l = NV_LENGTH_S(y);
    intnat ydot_l = NV_LENGTH_S(ydot);

    value y_ba = caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(y), &y_l);
    value ydot_ba = caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(ydot), &ydot_l);

    // XXX: the data payloads inside y_ba and ydot_ba are only valid during
    //      this call, afterward that memory goes back to cvode. These
    //      bigarrays must not be retained by closure_f! If it wants a
    //      permanent copy, then it has to make it manually.
    //
    //      Eventually y_ba and ydot_ba will be reclaimed by the ocaml gc,
    //      which should not, however, free the attached payload.
    value r = caml_callback3_exn(*closure_f, caml_copy_double(t), y_ba, ydot_ba);

    Caml_ba_array_val(y_ba)->dim[0] = 0;
    Caml_ba_array_val(ydot_ba)->dim[0] = 0;

    CAMLreturn(check_exception(r));
}

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    CAMLparam0();

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    intnat y_l = NV_LENGTH_S(y);
    value y_ba = caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(y), &y_l);

    value gout_ba = caml_ba_alloc(BIGARRAY_FLOAT, 1, gout, &(data->num_roots));

    // XXX: see notes for f()
    value r = caml_callback3_exn(*(data->closure_roots), caml_copy_double(t),
				 y_ba, gout_ba);

    Caml_ba_array_val(y_ba)->dim[0] = 0;
    Caml_ba_array_val(gout_ba)->dim[0] = 0;

    CAMLreturn(check_exception(r));
}

static mlsize_t approx_size_cvode_mem(void *cvode_mem)
{
    mlsize_t used = 0;
    long int lenrw = 0;
    long int leniw = 0;
    int flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);

    if (flag == CV_SUCCESS) {
    	used = lenrw * sizeof(realtype) + leniw * sizeof(long int);
    }

    return used;
}

static void set_linear_solver(void *cvode_mem, value ls, int n)
{
    printf("set_linear_solver: 0\n"); // XXX
    int flag;

    if (Is_block(ls)) {
	printf("Tag_val(ls)=%d\n", Tag_val(ls)); // XXX
	int field0 = Field(ls, 0); /* mupper, pretype */
	int field1 = Field(ls, 1); /* mlower, maxl */

	switch (Tag_val(ls)) {
	case VARIANT_LINEAR_SOLVER_BAND:
	    flag = CVBand(cvode_mem, n, field0, field1);
	    check_flag("CVBand", flag, NULL);
	    break;

#if SUNDIALS_BLAS_LAPACK == 1
	case VARIANT_LINEAR_SOLVER_LAPACKBAND:
	    field0 = Field(Field(ls, 1), 0);
	    field1 = Field(Field(ls, 1), 1);
	    flag = CVLapackBand(cvode_mem, n, field0, field1);
	    check_flag("CVLapackBand", flag, NULL);
	    break;
#endif

	case VARIANT_LINEAR_SOLVER_SPGMR:
	    field0 = Field(Field(ls, 1), 0);
	    field1    = Field(Field(ls, 1), 1);
	    flag = CVSpgmr(cvode_mem, field0, field1);
	    check_flag("CVSpgmr", flag, NULL);
	    break;

	case VARIANT_LINEAR_SOLVER_SPBCG:
	    field0 = Field(Field(ls, 1), 0);
	    field1    = Field(Field(ls, 1), 1);
	    flag = CVSpbcg(cvode_mem, field0, field1);
	    check_flag("CVSpbcg", flag, NULL);
	    break;

	case VARIANT_LINEAR_SOLVER_SPTFQMR:
	    field0 = Field(Field(ls, 1), 0);
	    field1    = Field(Field(ls, 1), 1);
	    flag = CVSptfqmr(cvode_mem, field0, field1);
	    check_flag("CVSPtfqmr", flag, NULL);
	    break;

	default:
	    caml_failwith("Illegal linear solver block value.");
	    break;
	}

    } else {
	printf("Int_val(ls)=%d\n", Int_val(ls)); // XXX
	switch (Int_val(ls)) {
	case VARIANT_LINEAR_SOLVER_DENSE:
	    flag = CVDense(cvode_mem, n);
	    check_flag("CVDense", flag, NULL);
	    break;

#if SUNDIALS_BLAS_LAPACK == 1
	case VARIANT_LINEAR_SOLVER_LAPACKDENSE:
	    flag = CVLapackDense(cvode_mem, n);
	    check_flag("CVLapackDense", flag, NULL);
	    break;
#endif

	case VARIANT_LINEAR_SOLVER_DIAG:
	    flag = CVDiag(cvode_mem);
	    check_flag("CVDiag", flag, NULL);
	    break;

	default:
	    caml_failwith("Illegal linear solver value.");
	    break;
	}
    }

    printf("set_linear_solver: 0\n"); // XXX
}
 
CAMLprim value c_init(value lmm, value iter, value initial, value num_roots)
{
    CAMLparam4(lmm, iter, initial, num_roots);

    int flag;

    int lmm_c;
    switch (Int_val(lmm)) {
    case VARIANT_LMM_ADAMS:
	lmm_c = CV_ADAMS;
	break;

    case VARIANT_LMM_BDF:
	lmm_c = CV_BDF;
	break;

    default:
	caml_failwith("Illegal lmm value.");
	break;
    }

    int iter_c;
    if (Is_block(iter)) {
	iter_c = CV_NEWTON;
    } else {
	iter_c = CV_FUNCTIONAL;
    }

    int initial_l = Caml_ba_array_val(initial)->dim[0];
    realtype *initial_d = Caml_ba_data_val(initial);
    N_Vector initial_nv = N_VMake_Serial(initial_l, initial_d);

    void *cvode_mem = CVodeCreate(lmm_c, iter_c);

    value vdata = ml_cvode_data_alloc(approx_size_cvode_mem(cvode_mem));
    ml_cvode_data_p data = (ml_cvode_data_p)Data_custom_val(vdata);

    data->cvode_mem = cvode_mem;
    data->neq = initial_l;
    data->err_file = NULL;
    data->closure_f = caml_named_value("cvode_serial_callback_f");
    data->closure_roots = caml_named_value("cvode_serial_callback_roots");
    data->closure_errh = NULL;
    // TODO: check if these two calls are necessary and ok:
    caml_register_generational_global_root(data->closure_f);
    caml_register_generational_global_root(data->closure_roots);

    data->num_roots = Int_val(num_roots);

    if (data->cvode_mem == NULL) {
	free(data);
	caml_failwith("CVodeCreate returned NULL");
	CAMLreturn0;
    }

    flag = CVodeInit(data->cvode_mem, f, T0, initial_nv);
    N_VDestroy(initial_nv);
    check_flag("CVodeInit", flag, data);

    flag = CVodeRootInit(data->cvode_mem, data->num_roots, roots);
    check_flag("CVodeRootInit", flag, data);

    CVodeSetUserData(data->cvode_mem, (void *)data);

    // setup linear solvers (if necessary)
    if (iter_c == CV_NEWTON) {
	set_linear_solver(data->cvode_mem, Field(iter, 0), data->neq);
    }

    // default tolerances
    N_Vector abstol = N_VNew_Serial(data->neq); 
    int i;
    for (i=0; i < data->neq; ++i) {
	NV_Ith_S(abstol, i) = RCONST(1.0e-8);
    }
    flag = CVodeSVtolerances(data->cvode_mem, RCONST(1.0e-4), abstol);
    check_flag("CVodeSVtolerances", flag, data);
    N_VDestroy_Serial(abstol);

    CAMLreturn(vdata);
}

CAMLprim value c_neqs(value vdata)
{
    CAMLparam1(vdata);
    ml_cvode_data_p data = cvode_data_from_ml(vdata);
    CAMLreturn(Val_int(data->neq));
}

CAMLprim value c_nroots(value vdata)
{
    CAMLparam1(vdata);
    ml_cvode_data_p data = cvode_data_from_ml(vdata);
    CAMLreturn(Val_int(data->num_roots));
}

CAMLprim value c_set_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int atol_l = Caml_ba_array_val(abstol)->dim[0];
    realtype *atol_d = Caml_ba_data_val(abstol);
    N_Vector atol_nv = N_VMake_Serial(atol_l, atol_d);

    int flag = CVodeSVtolerances(data->cvode_mem, Double_val(reltol), atol_nv);
    N_VDestroy(atol_nv);
    check_flag("CVodeSVtolerances", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_reinit(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int y0_l = Caml_ba_array_val(y0)->dim[0];
    realtype *y0_d = Caml_ba_data_val(y0);
    N_Vector y0_nv = N_VMake_Serial(y0_l, y0_d);

    int flag = CVodeReInit(data->cvode_mem, Double_val(t0), y0_nv);
    N_VDestroy(y0_nv);
    check_flag("CVodeReInit", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_get_roots(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = Caml_ba_data_val(roots);

    if (roots_l < data->num_roots) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = CVodeGetRootInfo(data->cvode_mem, roots_d);
    check_flag("CVodeGetRootInfo", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_free(value vdata)
{
    CAMLparam1(vdata);
    finalize(vdata);
    Store_field(vdata, 1, (value)NULL);
    CAMLreturn0;
}

static void check_flag(const char *call, int flag, void *to_free)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	|| flag == CV_ROOT_RETURN
	|| flag == CV_TSTOP_RETURN) return;

    if (to_free != NULL) free(to_free);

    switch (flag) {
    case CV_ILL_INPUT:
	caml_raise_constant(*caml_named_value("cvode_IllInput"));
	break;

    case CV_TOO_CLOSE:
	caml_raise_constant(*caml_named_value("cvode_TooClose"));
	break;

    case CV_TOO_MUCH_WORK:
	caml_raise_constant(*caml_named_value("cvode_TooMuchWork"));
	break;

    case CV_TOO_MUCH_ACC:
	caml_raise_constant(*caml_named_value("cvode_TooMuchAccuracy"));
	break;

    case CV_ERR_FAILURE:
	caml_raise_constant(*caml_named_value("cvode_ErrFailure"));
	break;

    case CV_CONV_FAILURE:
	caml_raise_constant(*caml_named_value("cvode_ConvergenceFailure"));
	break;

    case CV_LINIT_FAIL:
	caml_raise_constant(*caml_named_value("cvode_LinearInitFailure"));
	break;

    case CV_LSETUP_FAIL:
	caml_raise_constant(*caml_named_value("cvode_LinearSetupFailure"));
	break;

    case CV_LSOLVE_FAIL:
	caml_raise_constant(*caml_named_value("cvode_LinearSolveFailure"));
	break;

    case CV_RHSFUNC_FAIL:
	caml_raise_constant(*caml_named_value("cvode_RhsFuncFailure"));
	break;

    case CV_FIRST_RHSFUNC_ERR:
	caml_raise_constant(*caml_named_value("cvode_FirstRhsFuncError"));
	break;

    case CV_REPTD_RHSFUNC_ERR:
	caml_raise_constant(*caml_named_value("cvode_RepeatedRhsFuncError"));
	break;

    case CV_UNREC_RHSFUNC_ERR:
	caml_raise_constant(*caml_named_value("cvode_UnrecoverableRhsFuncError"));
	break;

    case CV_RTFUNC_FAIL:
	caml_raise_constant(*caml_named_value("cvode_RootFuncFailure"));
	break;

    default:
	/* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 CVodeGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

static value solver(value vdata, value nextt, value y, int onestep)
{
    CAMLparam2(vdata, nextt);

    realtype t = 0.0;
    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int leny = Bigarray_val(y)->dim[0];

    N_Vector y_nv = N_VMake_Serial(leny, Caml_ba_data_val(y));

    // TODO:
    // The payload of y (a big array) must not be shifted by the Ocaml GC
    // during this function call, even though Caml will be reentered
    // through the callback f. Is this guaranteed?
    int flag = CVode(data->cvode_mem, Double_val(nextt), y_nv, &t,
		     onestep ? CV_ONE_STEP : CV_NORMAL);
    N_VDestroy(y_nv);
    check_flag("CVode", flag, NULL);

    value r = caml_alloc_tuple(2);
    Store_field(r, 0, caml_copy_double(t));

    switch (flag) {
    case CV_ROOT_RETURN:
	Store_field(r, 1, Val_int(VARIANT_SOLVER_RESULT_ROOTSFOUND));
	break;

    case CV_TSTOP_RETURN:
	Store_field(r, 1, Val_int(VARIANT_SOLVER_RESULT_STOPTIMEREACHED));
	break;

    default:
	Store_field(r, 1, Val_int(VARIANT_SOLVER_RESULT_CONTINUE));
    }

    CAMLreturn(r);
}

CAMLprim value c_advance(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value c_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 1));
}

CAMLprim value c_integrator_stats(value vdata)
{
    CAMLparam1(vdata);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);
    int flag;

    long int nsteps;
    long int nfevals;    
    long int nlinsetups;
    long int netfails;

    int qlast;
    int qcur;	 

    realtype hinused;
    realtype hlast;
    realtype hcur;
    realtype tcur;

    flag = CVodeGetIntegratorStats(data->cvode_mem,
	&nsteps,
	&nfevals,    
	&nlinsetups,
	&netfails,
	&qlast,
	&qcur,	 
	&hinused,
	&hlast,
	&hcur,
	&tcur
    ); 
    check_flag("CVodeGetIntegratorStats", flag, NULL);

    value r = caml_alloc_tuple(10);
    Store_field(r, RECORD_INTEGRATOR_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_INTEGRATOR_STATS_RHS_EVALS, Val_long(nfevals));
    Store_field(r, RECORD_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS, Val_long(nlinsetups));
    Store_field(r, RECORD_INTEGRATOR_STATS_ERROR_TEST_FAILURES, Val_long(netfails));

    Store_field(r, RECORD_INTEGRATOR_STATS_LAST_INTERNAL_ORDER, Val_int(qlast));
    Store_field(r, RECORD_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER, Val_int(qcur));

    Store_field(r, RECORD_INTEGRATOR_STATS_INITIAL_STEP_SIZE, caml_copy_double(hinused));
    Store_field(r, RECORD_INTEGRATOR_STATS_LAST_STEP_SIZE, caml_copy_double(hlast));
    Store_field(r, RECORD_INTEGRATOR_STATS_NEXT_STEP_SIZE, caml_copy_double(hcur));
    Store_field(r, RECORD_INTEGRATOR_STATS_INTERNAL_TIME, caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim value c_last_step_size(value vdata)
{
    CAMLparam1(vdata);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);
    int flag;
    realtype hlast;

    flag = CVodeGetLastStep(data->cvode_mem, &hlast);
    check_flag("CVodeGetLastStep", flag, NULL);

    CAMLreturn(caml_copy_double(hlast));
}

CAMLprim value c_next_step_size(value vdata)
{
    CAMLparam1(vdata);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);
    int flag;
    realtype hcur;

    flag = CVodeGetCurrentStep(data->cvode_mem, &hcur);
    check_flag("CVodeGetCurrentStep", flag, NULL);

    CAMLreturn(caml_copy_double(hcur));
}

/* optional input functions */

CAMLprim value c_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    if (data->err_file != NULL) {
	fclose(data->err_file);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    data->err_file = fopen(String_val(vpath), mode);
    if (data->err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = CVodeSetErrFile(data->cvode_mem, data->err_file);
    check_flag("CVodeSetErrFile", flag, NULL);

    CAMLreturn0;
}

static void errh(
    int error_code,
    const char *module,
    const char *func,
    char *msg,
    void *eh_data)
{
    CAMLparam0();
    value *closure_errh = ((ml_cvode_data_p)eh_data)->closure_errh;

    value a = caml_alloc_tuple(4);
    Store_field(a, RECORD_ERROR_DETAILS_ERROR_CODE, Val_int(error_code));
    Store_field(a, RECORD_ERROR_DETAILS_MODULE_NAME, caml_copy_string(module));
    Store_field(a, RECORD_ERROR_DETAILS_FUNCTION_NAME, caml_copy_string(func));
    Store_field(a, RECORD_ERROR_DETAILS_ERROR_MESSAGE, caml_copy_string(msg));

    caml_callback(*closure_errh, a);

    CAMLreturn0;
}

// external set_error_handler : session -> (error_details -> unit) -> unit 
CAMLprim value c_set_error_handler(value vdata)
{
    CAMLparam1(vdata);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    if (data->closure_errh != NULL) {
	caml_remove_generational_global_root(data->closure_errh);
    }
    data->closure_errh = caml_named_value("cvode_serial_callback_errh");
    caml_register_generational_global_root(data->closure_errh);

    int flag = CVodeSetErrHandlerFn(data->cvode_mem, errh, (void *)data);
    check_flag("CVodeSetErrHandlerFn", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_ord(value vdata, value maxord)
{
    CAMLparam2(vdata, maxord);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMaxOrd(data->cvode_mem, Int_val(maxord));
    check_flag("CVodeSetMaxOrd", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_num_steps(value vdata, value mxsteps)
{
    CAMLparam2(vdata, mxsteps);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMaxNumSteps(data->cvode_mem, Long_val(mxsteps));
    check_flag("CVodeSetMaxNumSteps", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_hnil_warns(value vdata, value mxhnil)
{
    CAMLparam2(vdata, mxhnil);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMaxHnilWarns(data->cvode_mem, Int_val(mxhnil));
    check_flag("CVodeSetMaxHnilWarns", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_stability_limit_detection(value vdata, value stldet)
{
    CAMLparam2(vdata, stldet);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetStabLimDet(data->cvode_mem, Bool_val(stldet));
    check_flag("CVodeSetStabLimDet", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_initial_step_size(value vdata, value hin)
{
    CAMLparam2(vdata, hin);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetInitStep(data->cvode_mem, Double_val(hin));
    check_flag("CVodeSetInitStep", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_min_abs_step_size(value vdata, value hmin)
{
    CAMLparam2(vdata, hmin);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMinStep(data->cvode_mem, Double_val(hmin));
    check_flag("CVodeSetMinStep", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_abs_step_size(value vdata, value hmax)
{
    CAMLparam2(vdata, hmax);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMaxStep(data->cvode_mem, Double_val(hmax));
    check_flag("CVodeSetMaxStep", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_stop_time(value vdata, value tstop)
{
    CAMLparam2(vdata, tstop);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetStopTime(data->cvode_mem, Double_val(tstop));
    check_flag("CVodeSetStopTime", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_error_test_failures(value vdata, value maxnef)
{
    CAMLparam2(vdata, maxnef);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMaxErrTestFails(data->cvode_mem, Int_val(maxnef));
    check_flag("CVodeSetMaxErrTestFails", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_nonlinear_iterations(value vdata, value maxcor)
{
    CAMLparam2(vdata, maxcor);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMaxNonlinIters(data->cvode_mem, Int_val(maxcor));
    check_flag("CVodeSetMaxNonlinIters", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_convergence_failures(value vdata, value maxncf)
{
    CAMLparam2(vdata, maxncf);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetMaxConvFails(data->cvode_mem, Int_val(maxncf));
    check_flag("CVodeSetMaxConvFails", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_nonlinear_convergence_coeffficient(value vdata, value nlscoef)
{
    CAMLparam2(vdata, nlscoef);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetNonlinConvCoef(data->cvode_mem, Double_val(nlscoef));
    check_flag("CVodeSetNonlinConvCoef", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_nonlinear_iteration_type(value vdata, value iter)
{
    CAMLparam2(vdata, iter);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int iter_c;
    if (Is_block(iter)) {
	iter_c = CV_NEWTON;
    } else {
	iter_c = CV_FUNCTIONAL;
    }

    int flag = CVodeSetIterType(data->cvode_mem, iter_c);
    check_flag("CVodeSetIterType", flag, NULL);

    if (iter_c == CV_NEWTON) {
	set_linear_solver(data->cvode_mem, Field(iter, 0), data->neq);
    }

    CAMLreturn0;
}

CAMLprim value c_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = Caml_ba_data_val(rootdirs);

    if (rootdirs_l < data->num_roots) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = CVodeSetRootDirection(data->cvode_mem, rootdirs_d);
    check_flag("CVodeSetRootDirection", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_disable_inactive_root_warnings(value vdata)
{
    CAMLparam1(vdata);

    ml_cvode_data_p data = cvode_data_from_ml(vdata);

    int flag = CVodeSetNoInactiveRootWarn(data->cvode_mem);
    check_flag("CVodeSetNoInactiveRootWarn", flag, NULL);

    CAMLreturn0;
}

