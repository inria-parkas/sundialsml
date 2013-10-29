/***********************************************************************
 *                                                                     *
 *     OCaml interface to Sundials (serial) CVODE and IDA solvers      *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

/* The parts of the Sundials interface that distinguish between Serial
   NVectors (handled by Bigarrays) and generic NVectors (handled by a
   wrapper type). */

#include <cvode/cvode.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

/* linear solvers */
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#include <cvode/cvode_bandpre.h>

#include "cvode_ml.h"
#include "nvector_ml.h"

#ifdef RESTRICT_INTERNAL_PRECISION
#ifdef __GNUC__
#include <fpu_control.h>
#endif
#endif

// Call with CVODE_ML_BIGARRAYS to compile for the Serial NVector to
// Bigarray interface code.

#ifdef CVODE_ML_BIGARRAYS

#define CVTYPE(fname) c_ba_cvode_ ## fname
#include <nvector/nvector_serial.h>

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0

#define NVECTORIZE_VAL(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDVAL(nv) N_VDestroy(nv)

#else

#define CVTYPE(fname) c_nvec_cvode_ ## fname
#include <sundials/sundials_nvector.h>

#define WRAP_NVECTOR(v) NVEC_DATA(v)
#define RELINQUISH_WRAPPEDNV(v) {}

#define NVECTORIZE_VAL(v) NVEC_VAL(v)
#define RELINQUISH_NVECTORIZEDVAL(nv) {}

#endif

#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define CVTYPESTR(fname) QUOTE(CVTYPE(fname))

/* callbacks */

#define CAML_FN(name)					\
    static value *name;					\
    if (name == NULL)					\
	name = caml_named_value (CVTYPESTR (name));

static void errh(
	int error_code,
	const char *module,
	const char *func,
	char *msg,
	void *eh_data)
{
    CAMLparam0();
    CAMLlocal1(a);
    value *backref = eh_data;

    CAML_FN (call_errh);

    a = caml_alloc_tuple(4);
    Store_field(a, RECORD_CVODE_ERROR_DETAILS_ERROR_CODE,
                Val_int(error_code));
    Store_field(a, RECORD_CVODE_ERROR_DETAILS_MODULE_NAME,
                caml_copy_string(module));
    Store_field(a, RECORD_CVODE_ERROR_DETAILS_FUNCTION_NAME,
                caml_copy_string(func));
    Store_field(a, RECORD_CVODE_ERROR_DETAILS_ERROR_MESSAGE,
                caml_copy_string(msg));

    caml_callback2(*call_errh, *backref, a);

    CAMLreturn0;
}

CAMLprim void CVTYPE(set_err_handler_fn)(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = CVodeSetErrHandlerFn(CVODE_MEM_FROM_ML(vdata), errh,
				    CVODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(clear_err_handler_fn)(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSetErrHandlerFn(CVODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn0;
}

#define CHECK_RECOVERABLE      1
#define DONT_CHECK_RECOVERABLE 0
static int check_exception(value session, value r, int check_recoverable)
{
    CAMLparam2(session, r);
    CAMLlocal1(exn);

    static value *recoverable_failure = NULL;

    if (!Is_exception_result(r)) return 0;

    r = Extract_exception(r);

    if (check_recoverable) {
	if (recoverable_failure == NULL) {
	    recoverable_failure =
		caml_named_value("cvode_RecoverableFailure");
	}
	if (Field(r, 0) == *recoverable_failure)
	    CAMLreturnT (int, 1);
    }

    /* Unrecoverable error.  Save the exception and return -1.  */
    exn = caml_alloc_small (1,0);
    Field (exn,0) = r;
    Store_field (session, RECORD_CVODE_SESSION_EXN_TEMP, exn);
    CAMLreturnT (int, -1);
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_rhsfn);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = WRAP_NVECTOR(y);
    args[3] = WRAP_NVECTOR(ydot);

    // the data payloads inside y_d (args[2]) and ydot_d (args[3]) are only
    // valid during this call, afterward that memory goes back to cvode.
    // These bigarrays must not be retained by closure_rhsfn! If it wants a
    // permanent copy, then it has to make it manually.
    r = Int_val (caml_callbackN(*call_rhsfn,
                                sizeof (args) / sizeof (*args),
                                args));

    RELINQUISH_WRAPPEDNV(args[2]);
    RELINQUISH_WRAPPEDNV(args[3]);

    CAMLreturnT(int, r);
}

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, r);
    CAMLlocalN(args, 3);

    value *backref = user_data;
    intnat nroots;

    /* The length of gout is only available at the nroots field of the session
     * structure, so a dereference of the backref is unavoidable.  Therefore,
     * we do all of the setup here and directly call the user-supplied OCaml
     * function without going through an OCaml trampoline.  */
    WEAK_DEREF (session, *backref);

    nroots = CVODE_NROOTS_FROM_ML (session);

    args[0] = caml_copy_double (t);
    args[1] = WRAP_NVECTOR (y);
    args[2] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    r = caml_callbackN_exn (CVODE_ROOTSFN_FROM_ML (session),
			    sizeof (args) / sizeof (*args),
			    args);

    RELINQUISH_WRAPPEDNV (args[1]);

    CAMLreturnT(int, check_exception(session, r, CHECK_RECOVERABLE));
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    int r;
    value *backref = user_data;
    CAML_FN (call_errw);

    args[0] = *backref;
    args[1] = WRAP_NVECTOR (y);
    args[2] = WRAP_NVECTOR (ewt);

    r = Int_val (caml_callbackN (*call_errw,
				 sizeof (args) / sizeof (*args),
				 args));

    RELINQUISH_WRAPPEDNV(args[1]);
    RELINQUISH_WRAPPEDNV(args[2]);

    CAMLreturnT (int, r);
}

static value make_jac_arg(realtype t, N_Vector y, N_Vector fy, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(4);
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_Y, WRAP_NVECTOR(y));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_FY, WRAP_NVECTOR(fy));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_TMP, tmp);

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

#define TRIPLE 1
#define SINGLE 0
static void relinquish_jac_arg(value arg, int triple)
{
    CAMLparam1(arg);
    CAMLlocal1(tmp);

    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_CVODE_JACOBIAN_ARG_JAC_Y));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_CVODE_JACOBIAN_ARG_JAC_FY));

    tmp = Field(arg, RECORD_CVODE_JACOBIAN_ARG_JAC_TMP);

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
	long int n,
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
    CAMLlocalN (args, 3);
    int r;
    value *backref = user_data;
    CAML_FN (call_jacfn);

    args[0] = *backref;
    args[1] = make_jac_arg (t, y, fy, make_triple_tmp (tmp1, tmp2, tmp3));
    args[2] = caml_alloc_final (2, NULL, 0, 1);
    Store_field (args[2], 1, (value)Jac);

    r = Int_val (caml_callbackN (*call_jacfn,
				 sizeof (args) / sizeof (*args),
				 args));

    relinquish_jac_arg(args[1], TRIPLE);
    // note: matrix is also invalid after the callback

    CAMLreturnT(int, r);
}

static int bandjacfn(
	long int N,
	long int mupper,
	long int mlower, 	 
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
    CAMLlocalN(args, 5);
    int r;
    value *backref = user_data;
    CAML_FN (call_bandjacfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));
    args[2] = Val_int(mupper);
    args[3] = Val_int(mlower);
    args[4] = caml_alloc_final(2, NULL, 0, 1);
    Store_field (args[4], 1, (value)Jac);

    r = Int_val (caml_callbackN(*call_bandjacfn,
                                sizeof (args) / sizeof (*args),
                                args));

    relinquish_jac_arg(args[1], TRIPLE);
    // note: args[4] is also invalid after the callback

    CAMLreturnT(int, r);
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
    CAMLlocal2(session, r);
    CAMLlocalN(args, 3);
    value *backref = user_data;

    WEAK_DEREF (session, *backref);

    args[0] = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));
    args[1] = Val_bool(jok);
    args[2] = caml_copy_double(gamma);

    r = caml_callbackN_exn(CVODE_PRESETUPFN_FROM_ML (session),
                           sizeof (args) / sizeof (*args),
                           args);

    relinquish_jac_arg(args[0], TRIPLE);

    if (!Is_exception_result(r)) {
	*jcurPtr = Bool_val(r);
    }

    CAMLreturnT(int, check_exception(session, r, CHECK_RECOVERABLE));
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
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_RHS, WRAP_NVECTOR(r));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_GAMMA,
                caml_copy_double(gamma));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_DELTA,
                caml_copy_double(delta));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_LEFT,
                lr == 1 ? Val_true : Val_false);

    CAMLreturn(v);
}

static CAMLprim void relinquish_spils_solve_arg(value arg)
{
    CAMLparam1(arg);
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_CVODE_SPILS_SOLVE_ARG_RHS));
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
    CAMLlocal1(rv);
    CAMLlocalN(args, 4);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_presolvefn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, fy, WRAP_NVECTOR(tmp));
    args[2] = make_spils_solve_arg(r, gamma, delta, lr);
    args[3] = WRAP_NVECTOR(z);

    retcode = Int_val (caml_callbackN(*call_presolvefn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    relinquish_jac_arg(args[1], SINGLE);
    relinquish_spils_solve_arg(args[2]);
    RELINQUISH_WRAPPEDNV(args[3]);

    CAMLreturnT(int, retcode);
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
    CAMLlocal1(r);
    CAMLlocalN(args, 4);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_jactimesfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, fy, WRAP_NVECTOR(tmp));
    args[2] = WRAP_NVECTOR(v);
    args[3] = WRAP_NVECTOR(Jv);

    retcode = Int_val (caml_callbackN(*call_jactimesfn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    relinquish_jac_arg(args[1], SINGLE);
    RELINQUISH_WRAPPEDNV(args[2]);
    RELINQUISH_WRAPPEDNV(args[3]);

    CAMLreturnT(int, retcode);
}

CAMLprim void CVTYPE(wf_tolerances)(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = CVodeWFtolerances(CVODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("CVodeWFtolerances", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_set_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetDenseJacFn(CVODE_MEM_FROM_ML(vdata), jacfn);
    CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_clear_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetDenseJacFn(CVODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_set_band_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetBandJacFn(CVODE_MEM_FROM_ML(vdata), bandjacfn);
    CHECK_FLAG("CVDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_clear_band_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetBandJacFn(CVODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("CVDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(set_preconditioner)(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVSpilsSetPreconditioner(CVODE_MEM_FROM_ML(vdata),
					presetupfn, presolvefn);
    CHECK_FLAG("CVSpilsSetPreconditioner", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(set_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVSpilsSetJacTimesVecFn(CVODE_MEM_FROM_ML(vdata), jactimesfn);
    CHECK_FLAG("CVSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(clear_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVSpilsSetJacTimesVecFn(CVODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("CVSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

/* basic interface */

/* CVodeCreate() + CVodeInit().  */
CAMLprim value CVTYPE(init)(value weakref, value lmm, value iter, value initial,
			    value t0)
{
    CAMLparam5(weakref, lmm, iter, initial, t0);
    CAMLlocal1(r);

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
    case VARIANT_CVODE_LMM_ADAMS:
	lmm_c = CV_ADAMS;
	break;

    case VARIANT_CVODE_LMM_BDF:
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

    void *cvode_mem = CVodeCreate(lmm_c, iter_c);
    if (cvode_mem == NULL)
	caml_failwith("CVodeCreate returned NULL");

    N_Vector initial_nv = NVECTORIZE_VAL(initial);
    flag = CVodeInit(cvode_mem, f, Double_val(t0), initial_nv);
    RELINQUISH_NVECTORIZEDVAL(initial_nv);
    if (flag != CV_SUCCESS) {
	CVodeFree (cvode_mem);
	CHECK_FLAG("CVodeInit", flag);
    }

    value *backref;
    backref = malloc (sizeof (*backref));
    if (backref == NULL) {
	CVodeFree (cvode_mem);
	caml_failwith ("Out of memory");
    }
    *backref = weakref;
    caml_register_generational_global_root (backref);
    CVodeSetUserData (cvode_mem, backref);

    r = caml_alloc_tuple (3);
    Store_field (r, 0, (value)cvode_mem);
    Store_field (r, 1, (value)backref);
    Store_field (r, 2, Val_long (0)); /* no err_file = NULL */

    CAMLreturn(r);
}

/* Set the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim void CVTYPE(root_init) (value vdata, value vnroots)
{
    CAMLparam2 (vdata, vnroots);
    void *cvode_mem = CVODE_MEM_FROM_ML (vdata);
    int nroots = Int_val (vnroots);
    int flag = CVodeRootInit (cvode_mem, nroots, roots);
    CHECK_FLAG ("CVodeRootInit", flag);
    Store_field (vdata, RECORD_CVODE_SESSION_NROOTS, vnroots);
    CAMLreturn0;
}

CAMLprim void CVTYPE(sv_tolerances)(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVECTORIZE_VAL(abstol);

    int flag = CVodeSVtolerances(CVODE_MEM_FROM_ML(vdata),
				 Double_val(reltol), atol_nv);
    RELINQUISH_NVECTORIZEDVAL(atol_nv);
    CHECK_FLAG("CVodeSVtolerances", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(reinit)(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    N_Vector y0_nv = NVECTORIZE_VAL(y0);
    int flag = CVodeReInit(CVODE_MEM_FROM_ML(vdata), Double_val(t0), y0_nv);
    RELINQUISH_NVECTORIZEDVAL(y0_nv);
    CHECK_FLAG("CVodeReInit", flag);

    CAMLreturn0;
}

static value solver(value vdata, value nextt, value vy, int onestep)
{
    CAMLparam3(vdata, nextt, vy);
    CAMLlocal1(ret);
    realtype tret;
    int flag;
    N_Vector y;
    enum cvode_solver_result_tag result;

    y = NVECTORIZE_VAL (vy);
    // Caml_ba_data_val(y) must not be shifted by the OCaml GC during this
    // function call, which calls Caml through the callback f.  Is this
    // guaranteed?
    flag = CVode (CVODE_MEM_FROM_ML (vdata), Double_val (nextt), y, &tret,
		  onestep ? CV_ONE_STEP : CV_NORMAL);
    RELINQUISH_NVECTORIZEDVAL (y);

    switch (flag) {
    case CV_SUCCESS:
	result = VARIANT_CVODE_SOLVER_RESULT_CONTINUE;
	break;

    case CV_ROOT_RETURN:
	result = VARIANT_CVODE_SOLVER_RESULT_ROOTSFOUND;
	break;

    case CV_TSTOP_RETURN:
	result = VARIANT_CVODE_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_CVODE_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_CVODE_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	CHECK_FLAG ("CVode", flag);
    }

    /* Hmm...should this go in the production code or not?  */
    if (Is_block (Field (vdata, RECORD_CVODE_SESSION_EXN_TEMP)))
	abort ();

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}

CAMLprim value CVTYPE(solve_normal)(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value CVTYPE(solve_one_step)(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 1));
}

CAMLprim void CVTYPE(get_dky)(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVECTORIZE_VAL(vy);

    int flag = CVodeGetDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
			   Int_val(vk), y_nv);
    CHECK_FLAG("CVodeGetDky", flag);
    RELINQUISH_NVECTORIZEDVAL(y_nv);
    
    CAMLreturn0;
}

CAMLprim void CVTYPE(get_err_weights)(value vcvode_mem, value verrws)
{
    CAMLparam2(vcvode_mem, verrws);

    N_Vector errws_nv = NVECTORIZE_VAL(verrws);

    int flag = CVodeGetErrWeights(CVODE_MEM_FROM_ML(vcvode_mem), errws_nv);
    RELINQUISH_NVECTORIZEDVAL(errws_nv);
    CHECK_FLAG("CVodeGetErrWeights", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(get_est_local_errors)(value vcvode_mem, value vele)
{
    CAMLparam2(vcvode_mem, vele);

    N_Vector ele_nv = NVECTORIZE_VAL(vele);

    int flag = CVodeGetEstLocalErrors(CVODE_MEM_FROM_ML(vcvode_mem), ele_nv);
    RELINQUISH_NVECTORIZEDVAL(ele_nv);
    CHECK_FLAG("CVodeGetEstLocalErrors", flag);

    CAMLreturn0;
}

