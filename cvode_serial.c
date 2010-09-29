/* Aug 2010, Timothy Bourke (INRIA)
 *
 * Trickier definitions for Sundials interface.
 *
 */

/*
 * TODO:
 * - call Gc.full_major () from the f and roots routines to see if we get
 *   any segmentation fault problems.
 * - see notes throughout program related to garbage collection.
 */

#include "cvode_serial.h"

#ifdef RESTRICT_INTERNAL_PRECISION
#ifdef __GNUC__
#include <fpu_control.h>
#endif
#endif

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

static const char *callback_ocaml_names[] = {
    "cvode_serial_callback_rhsfn",
    "cvode_serial_callback_rootsfn",
    "cvode_serial_callback_errorhandler",
    "cvode_serial_callback_errorweight",
    "cvode_serial_callback_jacfn",
    "cvode_serial_callback_bandjacfn",
    "cvode_serial_callback_presetupfn",
    "cvode_serial_callback_presolvefn",
    "cvode_serial_callback_jactimesfn",
};

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0
#define NVECTORIZE_BA(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDBA(nv) N_VDestroy(nv)

static void finalize_closure(value** closure_field)
{
    if (*closure_field != NULL) {
	caml_remove_generational_global_root(*closure_field);
	*closure_field = NULL;
    }
}

static void finalize(value vdata)
{
    CVODE_DATA_FROM_ML(data, vdata);

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
    finalize_closure(&(data->closure_rhsfn));
    finalize_closure(&(data->closure_rootsfn));
    finalize_closure(&(data->closure_errh));
    finalize_closure(&(data->closure_errw));

    finalize_closure(&(data->closure_jacfn));
    finalize_closure(&(data->closure_bandjacfn));
    finalize_closure(&(data->closure_presetupfn));
    finalize_closure(&(data->closure_presolvefn));
    finalize_closure(&(data->closure_jactimesfn));

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
value ml_cvode_data_alloc(void* cvode_mem)
{
    CAMLparam0();
    CAMLlocal1(r);

    mlsize_t approx_size = 0;
    long int lenrw, leniw;
    if (CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw) == CV_SUCCESS) {
    	approx_size = lenrw * sizeof(realtype) + leniw * sizeof(long int);
    }

    r = caml_alloc_final(sizeof(struct ml_cvode_data),
			 &finalize, approx_size, 10);

    ml_cvode_data_p d = CVODE_DATA(r);
    d->cvode_mem = cvode_mem;

    d->err_file = NULL;
    d->closure_rhsfn = NULL;
    d->closure_rootsfn = NULL;
    d->closure_errh = NULL;
    d->closure_errw = NULL;
    d->closure_jacfn = NULL;
    d->closure_bandjacfn = NULL;
    d->closure_presetupfn = NULL;
    d->closure_presolvefn = NULL;
    d->closure_jactimesfn = NULL;
    
    CAMLreturn(r);
}

void ml_cvode_check_flag(const char *call, int flag, void *to_free)
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

    case CV_BAD_K:
	caml_raise_constant(*caml_named_value("cvode_BadK"));
	break;

    case CV_BAD_T:
	caml_raise_constant(*caml_named_value("cvode_BadT"));
	break;

    case CV_BAD_DKY:
	caml_raise_constant(*caml_named_value("cvode_BadDky"));
	break;

    default:
	/* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 CVodeGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

static int precond_type(value vptype)
{
    CAMLparam1(vptype);

    int ptype;
    switch (Int_val(vptype)) {
    case VARIANT_PRECONDITIONING_TYPE_PRECNONE:
	ptype = PREC_NONE;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECLEFT:
	ptype = PREC_LEFT;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECRIGHT:
	ptype = PREC_RIGHT;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECBOTH:
	ptype = PREC_BOTH;
	break;
    }

    CAMLreturn(ptype);
}

static int check_exception(value r)
{
    CAMLparam0();
    CAMLlocal1(exn);

    static value *recoverable_failure = NULL;

    if (!Is_exception_result(r)) return 0;

    if (recoverable_failure == NULL) {
	recoverable_failure =
	    caml_named_value("cvode_RecoverableFailure");
    }

    exn = Extract_exception(r);
    CAMLreturn((Field(exn, 0) == *recoverable_failure) ? 1 : -1);
}

/* callbacks */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(y_ba, ydot_ba, r);

    value *closure_rhsfn = ((ml_cvode_data_p)user_data)->closure_rhsfn;

    y_ba = WRAP_NVECTOR(y);
    ydot_ba = WRAP_NVECTOR(ydot);

    // TODO: the data payloads inside y_ba and ydot_ba are only valid
    //	     during this call, afterward that memory goes back to cvode.
    //	     These bigarrays must not be retained by closure_rhsfn! If
    //	     it wants a permanent copy, then it has to make it manually.
    //
    //       Eventually y_ba and ydot_ba will be reclaimed by the ocaml gc,
    //       which should not, however, free the attached payload.
    r = caml_callback3_exn(*closure_rhsfn, caml_copy_double(t), y_ba, ydot_ba);

    RELINQUISH_WRAPPEDNV(y_ba);
    RELINQUISH_WRAPPEDNV(ydot_ba);

    CAMLreturn(check_exception(r));
}

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(y_ba, gout_ba, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    y_ba = WRAP_NVECTOR(y);
    gout_ba = caml_ba_alloc(BIGARRAY_FLOAT, 1, gout, &(data->num_roots));

    // TODO: see notes for f()
    r = caml_callback3_exn(*(data->closure_rootsfn), caml_copy_double(t),
				 y_ba, gout_ba);

    RELINQUISH_WRAPPEDNV(y_ba);
    Caml_ba_array_val(gout_ba)->dim[0] = 0;

    CAMLreturn(check_exception(r));
}

static void errh(
	int error_code,
	const char *module,
	const char *func,
	char *msg,
	void *eh_data)
{
    CAMLparam0();
    CAMLlocal1(a);

    value *closure_errh = ((ml_cvode_data_p)eh_data)->closure_errh;

    a = caml_alloc_tuple(4);
    Store_field(a, RECORD_ERROR_DETAILS_ERROR_CODE, Val_int(error_code));
    Store_field(a, RECORD_ERROR_DETAILS_MODULE_NAME, caml_copy_string(module));
    Store_field(a, RECORD_ERROR_DETAILS_FUNCTION_NAME, caml_copy_string(func));
    Store_field(a, RECORD_ERROR_DETAILS_ERROR_MESSAGE, caml_copy_string(msg));

    caml_callback(*closure_errh, a);

    CAMLreturn0;
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(y_ba, ewt_ba, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    y_ba = WRAP_NVECTOR(y);
    ewt_ba = WRAP_NVECTOR(ewt);

    // TODO: see notes for f()
    r = caml_callback2_exn(*(data->closure_errw), y_ba, ewt_ba);

    RELINQUISH_WRAPPEDNV(y_ba);
    RELINQUISH_WRAPPEDNV(ewt_ba);

    CAMLreturn(check_exception(r));
}

CAMLprim value c_enable_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
 
    int flag = CVodeSetErrHandlerFn(data->cvode_mem, errh, (void *)data);
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim value c_wf_tolerances(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
 
    int flag = CVodeWFtolerances(data->cvode_mem, errw);
    CHECK_FLAG("CVodeWFtolerances", flag);

    CAMLreturn0;
}

static value make_jac_arg(realtype t, N_Vector y, N_Vector fy, value tmp)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(4);
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_Y, WRAP_NVECTOR(y));
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_FY, WRAP_NVECTOR(fy));
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

static value make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(3);
    Store_field(r, 0, WRAP_NVECTOR(tmp1));
    Store_field(r, 1, WRAP_NVECTOR(tmp2));
    Store_field(r, 2, WRAP_NVECTOR(tmp3));
    CAMLreturn(r);
}

static relinquish_jac_arg(value arg, int triple)
{
    CAMLparam0();
    CAMLlocal1(tmp);

    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_JACOBIAN_ARG_JAC_Y));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_JACOBIAN_ARG_JAC_FY));

    tmp = Field(arg, RECORD_JACOBIAN_ARG_JAC_TMP);

    if (triple) {
	RELINQUISH_WRAPPEDNV(Field(tmp, 0));
	RELINQUISH_WRAPPEDNV(Field(tmp, 1));
	RELINQUISH_WRAPPEDNV(Field(tmp, 2));
    } else {
	RELINQUISH_WRAPPEDNV(tmp);
    }

    CAMLreturn0;
}

static int jacfn(
	int n,
	realtype t,
	N_Vector y,
	N_Vector fy,	     
	DlsMat Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal3(arg, r, matrix);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));

    matrix = caml_alloc(1, Abstract_tag);
    Store_field(matrix, 0, (value)Jac);

    r = caml_callback2_exn(*(data->closure_jacfn), arg, matrix);
    relinquish_jac_arg(arg, 1);
    // note: matrix is also invalid after the callback

    CAMLreturn(check_exception(r));
}

static int bandjacfn(
	int N,
	int mupper,
	int mlower, 	 
	realtype t,
	N_Vector y,
	N_Vector fy, 	 
	DlsMat Jac,
	void *user_data, 	 
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal1(r);
    CAMLlocalN(args, 4);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    args[0] = Val_int(mupper);
    args[1] = Val_int(mlower);
    args[2] = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));
    args[3] = caml_alloc(1, Abstract_tag);
    Store_field(args[3], 0, (value)Jac);

    r = caml_callbackN_exn(*(data->closure_bandjacfn), 4, args);

    relinquish_jac_arg(args[2], 1);
    // note: matrix is also invalid after the callback

    CAMLreturn(check_exception(r));
}

static int presetupfn(
    realtype t,
    N_Vector y,
    N_Vector fy,
    booleantype jok,
    booleantype *jcurPtr,
    realtype gamma,
    void *user_data,
    N_Vector tmp1,
    N_Vector tmp2,
    N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal2(arg, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));

    r = caml_callback3_exn(*(data->closure_presetupfn),
	    arg, Val_bool(jok), caml_copy_double(gamma));

    relinquish_jac_arg(arg, 1);

    if (!Is_exception_result(r)) {
	*jcurPtr = Bool_val(r);
    }

    CAMLreturn(check_exception(r));
}

static value make_spils_solve_arg(
	N_Vector r,
	realtype gamma,
	realtype delta,
	int lr)

{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(4);
    Store_field(v, RECORD_SPILS_SOLVE_ARG_RHS, WRAP_NVECTOR(r));
    Store_field(v, RECORD_SPILS_SOLVE_ARG_GAMMA, caml_copy_double(gamma));
    Store_field(v, RECORD_SPILS_SOLVE_ARG_DELTA, caml_copy_double(delta));
    Store_field(v, RECORD_SPILS_SOLVE_ARG_LEFT, lr == 1 ? Val_true : Val_false);

    CAMLreturn(v);
}

static relinquish_spils_solve_arg(value arg)
{
    CAMLparam0();
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_SPILS_SOLVE_ARG_RHS));
    CAMLreturn0;
}

static int presolvefn(
	realtype t,
	N_Vector y,
	N_Vector fy,
	N_Vector r,
	N_Vector z,
	realtype gamma,
	realtype delta,
	int lr,
	void *user_data,
	N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal4(arg, solvearg, zv, rv);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, WRAP_NVECTOR(tmp));
    solvearg = make_spils_solve_arg(r, gamma, delta, lr);
    zv = WRAP_NVECTOR(z);

    rv = caml_callback3_exn(*(data->closure_presolvefn), arg, solvearg, zv);

    relinquish_jac_arg(arg, 0);
    relinquish_spils_solve_arg(solvearg);
    RELINQUISH_WRAPPEDNV(zv);

    CAMLreturn(check_exception(rv));
}

static int jactimesfn(
    N_Vector v,
    N_Vector Jv,
    realtype t,
    N_Vector y,
    N_Vector fy,
    void *user_data,
    N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal4(arg, varg, jvarg, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, WRAP_NVECTOR(tmp));
    varg = WRAP_NVECTOR(v);
    jvarg = WRAP_NVECTOR(Jv);

    r = caml_callback3_exn(*(data->closure_jactimesfn), arg, varg, jvarg);

    relinquish_jac_arg(arg, 0);
    RELINQUISH_WRAPPEDNV(varg);
    RELINQUISH_WRAPPEDNV(jvarg);

    CAMLreturn(check_exception(r));
}

CAMLprim value c_dls_enable_dense_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVDlsSetDenseJacFn(data->cvode_mem, jacfn);
    CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim value c_dls_enable_band_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVDlsSetBandJacFn(data->cvode_mem, bandjacfn);
    CHECK_FLAG("CVDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim value c_enable_preconditioner(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVSpilsSetPreconditioner(data->cvode_mem,
	    presetupfn, presolvefn);
    CHECK_FLAG("CVSpilsSetPreconditioner", flag);
    CAMLreturn0;
}

CAMLprim value c_enable_jac_times_vec_fn(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVSpilsSetJacTimesVecFn(data->cvode_mem, jactimesfn);
    CHECK_FLAG("CVSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

/* basic interface */

static void set_linear_solver(void *cvode_mem, value ls, int n)
{
    int flag;

    if (Is_block(ls)) {
	int field0 = Field(Field(ls, 0), 0); /* mupper, pretype */
	int field1 = Field(Field(ls, 0), 1); /* mlower, maxl */
	value sprange, bandrange;

	switch (Tag_val(ls)) {
	case VARIANT_LINEAR_SOLVER_BAND:
	    flag = CVBand(cvode_mem, n, Int_val(field0), Int_val(field1));
	    CHECK_FLAG("CVBand", flag);
	    break;

#if SUNDIALS_BLAS_LAPACK == 1
	case VARIANT_LINEAR_SOLVER_LAPACKBAND:
	    flag = CVLapackBand(cvode_mem, n, Int_val(field0), Int_val(field1));
	    CHECK_FLAG("CVLapackBand", flag);
	    break;
#endif

	case VARIANT_LINEAR_SOLVER_SPGMR:
	    flag = CVSpgmr(cvode_mem, precond_type(field0), Int_val(field1));
	    CHECK_FLAG("CVSpgmr", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_SPBCG:
	    flag = CVSpbcg(cvode_mem, precond_type(field0), Int_val(field1));
	    CHECK_FLAG("CVSpbcg", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_SPTFQMR:
	    flag = CVSptfqmr(cvode_mem, precond_type(field0), Int_val(field1));
	    CHECK_FLAG("CVSPtfqmr", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_BANDED_SPGMR:
	    sprange = Field(ls, 0);
	    bandrange = Field(ls, 1);

	    flag = CVSpgmr(cvode_mem, precond_type(Field(sprange, 0)),
				      Int_val(Field(sprange, 1)));
	    CHECK_FLAG("CVSpgmr", flag);

	    flag = CVBandPrecInit(cvode_mem, n, Int_val(Field(bandrange, 0)),
						Int_val(Field(bandrange, 1)));
	    CHECK_FLAG("CVBandPrecInit", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_BANDED_SPBCG:
	    sprange = Field(ls, 0);
	    bandrange = Field(ls, 1);

	    flag = CVSpbcg(cvode_mem, precond_type(Field(sprange, 0)),
				      Int_val(Field(sprange, 1)));
	    CHECK_FLAG("CVSpbcg", flag);

	    flag = CVBandPrecInit(cvode_mem, n, Int_val(Field(bandrange, 0)),
						Int_val(Field(bandrange, 1)));
	    CHECK_FLAG("CVBandPrecInit", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_BANDED_SPTFQMR:
	    sprange = Field(ls, 0);
	    bandrange = Field(ls, 1);

	    flag = CVSptfqmr(cvode_mem, precond_type(Field(sprange, 0)),
				        Int_val(Field(sprange, 1)));
	    CHECK_FLAG("CVSptfqmr", flag);

	    flag = CVBandPrecInit(cvode_mem, n, Int_val(Field(bandrange, 0)),
						Int_val(Field(bandrange, 1)));
	    CHECK_FLAG("CVBandPrecInit", flag);
	    break;

	default:
	    caml_failwith("Illegal linear solver block value.");
	    break;
	}

    } else {
	switch (Int_val(ls)) {
	case VARIANT_LINEAR_SOLVER_DENSE:
	    flag = CVDense(cvode_mem, n);
	    CHECK_FLAG("CVDense", flag);
	    break;

#if SUNDIALS_BLAS_LAPACK == 1
	case VARIANT_LINEAR_SOLVER_LAPACKDENSE:
	    flag = CVLapackDense(cvode_mem, n);
	    CHECK_FLAG("CVLapackDense", flag);
	    break;
#endif

	case VARIANT_LINEAR_SOLVER_DIAG:
	    flag = CVDiag(cvode_mem);
	    CHECK_FLAG("CVDiag", flag);
	    break;

	default:
	    caml_failwith("Illegal linear solver value.");
	    break;
	}
    }
}
 
CAMLprim value c_init(value lmm, value iter, value initial, value num_roots,
		      value t0)
{
    CAMLparam4(lmm, iter, initial, num_roots);
    CAMLlocal1(vdata);

    if (sizeof(int) != 4) {
	caml_failwith("The library assumes that an int (in C) has 32-bits.");
    }

#ifdef RESTRICT_INTERNAL_PRECISION
#ifdef __GNUC__
    fpu_control_t fpu_cw;
    _FPU_GETCW(fpu_cw);
    fpu_cw = (fpu_cw & ~_FPU_EXTENDED & ~_FPU_SINGLE) | _FPU_DOUBLE;
    _FPU_SETCW(fpu_cw);
#endif
#endif

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

    N_Vector initial_nv = NVECTORIZE_BA(initial);

    void *cvode_mem = CVodeCreate(lmm_c, iter_c);

    vdata = ml_cvode_data_alloc(cvode_mem);
    CVODE_DATA_FROM_ML(data, vdata);

    data->neq = Caml_ba_array_val(initial)->dim[0];
    data->num_roots = Int_val(num_roots);

    if (data->cvode_mem == NULL) {
	free(data);
	caml_failwith("CVodeCreate returned NULL");
	CAMLreturn0;
    }

    flag = CVodeInit(data->cvode_mem, f, Double_val(t0), initial_nv);
    RELINQUISH_NVECTORIZEDBA(initial_nv);
    ml_cvode_check_flag("CVodeInit", flag, data);

    if (data->num_roots > 0) {
	flag = CVodeRootInit(data->cvode_mem, data->num_roots, roots);
	ml_cvode_check_flag("CVodeRootInit", flag, data);
    }

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
    ml_cvode_check_flag("CVodeSVtolerances", flag, data);
    N_VDestroy_Serial(abstol);

    CAMLreturn(vdata);
}

CAMLprim value c_neqs(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    CAMLreturn(Val_int(data->neq));
}

CAMLprim value c_nroots(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    CAMLreturn(Val_int(data->num_roots));
}

CAMLprim value c_sv_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    CVODE_DATA_FROM_ML(data, vdata);

    N_Vector atol_nv = NVECTORIZE_BA(abstol);

    int flag = CVodeSVtolerances(data->cvode_mem, Double_val(reltol), atol_nv);
    RELINQUISH_NVECTORIZEDBA(atol_nv);
    CHECK_FLAG("CVodeSVtolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    CVODE_DATA_FROM_ML(data, vdata);

    int flag = CVodeSStolerances(data->cvode_mem,
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("CVodeSStolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_re_init(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    CVODE_DATA_FROM_ML(data, vdata);

    N_Vector y0_nv = NVECTORIZE_BA(y0);
    int flag = CVodeReInit(data->cvode_mem, Double_val(t0), y0_nv);
    RELINQUISH_NVECTORIZEDBA(y0_nv);
    CHECK_FLAG("CVodeReInit", flag);

    CAMLreturn0;
}

CAMLprim value c_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    CVODE_DATA_FROM_ML(data, vdata);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = Caml_ba_data_val(roots);

    if (roots_l < data->num_roots) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = CVodeGetRootInfo(data->cvode_mem, roots_d);
    CHECK_FLAG("CVodeGetRootInfo", flag);

    CAMLreturn0;
}

CAMLprim value c_free(value vdata)
{
    CAMLparam1(vdata);
    finalize(vdata);
    Store_field(vdata, 1, (value)NULL);
    CAMLreturn0;
}

static value solver(value vdata, value nextt, value y, int onestep)
{
    CAMLparam2(vdata, nextt);
    CAMLlocal1(r);

    realtype t = 0.0;
    CVODE_DATA_FROM_ML(data, vdata);

    int leny = Bigarray_val(y)->dim[0];

    N_Vector y_nv = N_VMake_Serial(leny, Caml_ba_data_val(y));

    // TODO:
    // The payload of y (a big array) must not be shifted by the Ocaml GC
    // during this function call, even though Caml will be reentered
    // through the callback f. Is this guaranteed?
    int flag = CVode(data->cvode_mem, Double_val(nextt), y_nv, &t,
		     onestep ? CV_ONE_STEP : CV_NORMAL);
    N_VDestroy(y_nv);
    CHECK_FLAG("CVode", flag);

    r = caml_alloc_tuple(2);
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

CAMLprim value c_normal(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value c_one_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 1));
}

CAMLprim value c_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    CVODE_DATA_FROM_ML(data, vdata);
    N_Vector y_nv = NVECTORIZE_BA(vy);

    int flag = CVodeGetDky(data->cvode_mem, Double_val(vt), Int_val(vk), y_nv);
    CHECK_FLAG("CVodeGetDky", flag);
    RELINQUISH_NVECTORIZEDBA(y_nv);
    
    CAMLreturn0;
}

CAMLprim value c_get_integrator_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    CVODE_DATA_FROM_ML(data, vdata);
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
    CHECK_FLAG("CVodeGetIntegratorStats", flag);

    r = caml_alloc_tuple(10);
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

CAMLprim value c_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    CVODE_DATA_FROM_ML(data, vdata);

    if (data->err_file != NULL) {
	fclose(data->err_file);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    data->err_file = fopen(String_val(vpath), mode);
    if (data->err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = CVodeSetErrFile(data->cvode_mem, data->err_file);
    CHECK_FLAG("CVodeSetErrFile", flag);

    CAMLreturn0;
}

CAMLprim value c_register_handler(value vdata, value handler)
{
    CAMLparam2(vdata, handler);

    CVODE_DATA_FROM_ML(data, vdata);
    const char* ocaml_name = callback_ocaml_names[Int_val(handler)];
    value **handler_field;

    switch (Int_val(handler)) {
    case VARIANT_HANDLER_RHSFN:
	handler_field = &(data->closure_rhsfn);
	break;

    case VARIANT_HANDLER_ROOTSFN:
	handler_field = &(data->closure_rootsfn);
	break;

    case VARIANT_HANDLER_ERRORHANDLER:
	handler_field = &(data->closure_errh);
	break;

    case VARIANT_HANDLER_ERRORWEIGHT:
	handler_field = &(data->closure_errw);
	break;

    case VARIANT_HANDLER_JACFN:
	handler_field = &(data->closure_jacfn);
	break;

    case VARIANT_HANDLER_BANDJACFN:
	handler_field = &(data->closure_bandjacfn);
	break;

    case VARIANT_HANDLER_PRESETUPFN:
	handler_field = &(data->closure_presetupfn);
	break;

    case VARIANT_HANDLER_PRESOLVEFN:
	handler_field = &(data->closure_presolvefn);
	break;

    case VARIANT_HANDLER_JACTIMESFN:
	handler_field = &(data->closure_jactimesfn);
	break;

    default:
	break;
    }

    if ((*handler_field) != NULL) {
	caml_remove_generational_global_root(*handler_field);
    }
    (*handler_field) = caml_named_value(ocaml_name);
    // TODO: check if this call is necessary and ok:
    caml_register_generational_global_root(*handler_field);

    CAMLreturn0;
}

CAMLprim value c_set_iter_type(value vdata, value iter)
{
    CAMLparam2(vdata, iter);

    CVODE_DATA_FROM_ML(data, vdata);

    int iter_c;
    if (Is_block(iter)) {
	iter_c = CV_NEWTON;
    } else {
	iter_c = CV_FUNCTIONAL;
    }

    int flag = CVodeSetIterType(data->cvode_mem, iter_c);
    CHECK_FLAG("CVodeSetIterType", flag);

    if (iter_c == CV_NEWTON) {
	set_linear_solver(data->cvode_mem, Field(iter, 0), data->neq);
    }

    CAMLreturn0;
}

CAMLprim value c_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    CVODE_DATA_FROM_ML(data, vdata);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = Caml_ba_data_val(rootdirs);

    if (rootdirs_l < data->num_roots) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = CVodeSetRootDirection(data->cvode_mem, rootdirs_d);
    CHECK_FLAG("CVodeSetRootDirection", flag);

    CAMLreturn0;
}

CAMLprim value c_get_err_weights(value vcvode_mem, value verrws)
{
    CAMLparam2(vcvode_mem, verrws);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    N_Vector errws_nv = NVECTORIZE_BA(verrws);

    int flag = CVodeGetErrWeights(cvode_mem, errws_nv);
    RELINQUISH_NVECTORIZEDBA(errws_nv);
    CHECK_FLAG("CVodeGetErrWeights", flag);

    CAMLreturn0;
}

CAMLprim value c_get_est_local_errors(value vcvode_mem, value vele)
{
    CAMLparam2(vcvode_mem, vele);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    N_Vector ele_nv = NVECTORIZE_BA(vele);

    int flag = CVodeGetEstLocalErrors(cvode_mem, ele_nv);
    RELINQUISH_NVECTORIZEDBA(ele_nv);
    CHECK_FLAG("CVodeGetEstLocalErrors", flag);

    CAMLreturn0;
}

CAMLprim value c_set_prec_type(value vcvode_mem, value vptype)
{
    CAMLparam2(vcvode_mem, vptype);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVSpilsSetPrecType(cvode_mem, precond_type(vptype));
    CHECK_FLAG("CVSpilsSetPrecType", flag);

    CAMLreturn0;
}

CAMLprim value c_vmax_norm(value u)
{
    CAMLparam0();
    CAMLlocal1(r);

    N_Vector u_nv = NVECTORIZE_BA(u);
    r = caml_copy_double(N_VMaxNorm(u_nv));
    RELINQUISH_NVECTORIZEDBA(u_nv);

    CAMLreturn(r);
}

