/* Aug 2010, Timothy Bourke (INRIA) */

#include "caml/mlvalues.h"
#include "caml/bigarray.h"
#include "caml/memory.h"
#include "caml/callback.h"
#include "caml/custom.h"
#include "caml/fail.h"

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

/*
 * Possible extensions:
 * - allow different solvers to be specified (in init)
 * - add a closure for returning the Jacobian.
 * - add handling for flag == CV_TSTOP_RETURN to solver
 * - allow T0 to be configured
 * - add an interface for CVodeSVtolerances(cvode_mem, reltol, abstol);
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
    intnat num_roots;
    value *closure_f;
    value *closure_roots;
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
    caml_remove_generational_global_root(data->closure_roots);

    if (data->cvode_mem != NULL) {
	CVodeFree(&(data->cvode_mem));
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
    static value *recoverable_failure = NULL;

    if (!Is_exception_result(r)) return 0;

    if (recoverable_failure == NULL) {
	recoverable_failure =
	    caml_named_value("cvode_RecoverableFailure");
    }

    value exn = Extract_exception(r);
    return (Field(exn, 0) == *recoverable_failure) ? 1 : -1;
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
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

    return check_exception(r);
}

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    intnat y_l = NV_LENGTH_S(y);
    value y_ba = caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(y), &y_l);

    value gout_ba = caml_ba_alloc(BIGARRAY_FLOAT, 1, gout, &(data->num_roots));

    // XXX: see notes for f()
    value r = caml_callback3_exn(*(data->closure_roots), caml_copy_double(t),
				 y_ba, gout_ba);

    Caml_ba_array_val(y_ba)->dim[0] = 0;
    Caml_ba_array_val(gout_ba)->dim[0] = 0;

    return check_exception(r);
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
    data->closure_f = caml_named_value("cvode_serial_callback_f");
    data->closure_roots = caml_named_value("cvode_serial_callback_roots");
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
	set_linear_solver(data->cvode_mem, Field(iter, 0), initial_l);
    }

    // default tolerances
    N_Vector abstol = N_VNew_Serial(initial_l); 
    int i;
    for (i=0; i < initial_l; ++i) {
	NV_Ith_S(abstol, i) = RCONST(1.0e-8);
    }
    flag = CVodeSVtolerances(data->cvode_mem, RCONST(1.0e-4), abstol);
    check_flag("CVodeSVtolerances", flag, data);
    N_VDestroy_Serial(abstol);

    CAMLreturn(vdata);
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
    Store_field(r, 1, flag == CV_ROOT_RETURN ? Val_true : Val_false);

    CAMLreturn(r);
}

CAMLprim value c_advance(value vdata, value nextt, value y)
{
    CAMLparam0();
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value c_step(value vdata, value nextt, value y)
{
    CAMLparam0();
    CAMLreturn(solver(vdata, nextt, y, 1));
}

