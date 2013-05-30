/***********************************************************************
 *                                                                     *
 *              Ocaml interface to Sundials CVODE solver               *
 *                                                                     *
 *           Timothy Bourke (INRIA) and Marc Pouzet (LIENS)            *
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under the terms of the GNU Library General Public License, with    *
 *  the special exception on linking described in file LICENSE.        *
 *                                                                     *
 ***********************************************************************/

/* The parts of the Sundials IDA interface that distinguish between Serial
   NVectors (handled by Bigarrays) and generic NVectors (handled by a wrapper
   type). */

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/fail.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/bigarray.h>

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <ida/ida_band.h>
#include <ida/ida_spgmr.h>

#include "ida_ml.h"
#include "nvector_ml.h"

// Call with IDA_ML_BIGARRAYS to compile for the Serial NVector to
// Bigarray interface code.

#ifdef IDA_ML_BIGARRAYS

#define IDATYPE(fname) c_ba_ida_ ## fname
#include <nvector/nvector_serial.h>

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0

#define NVECTORIZE_VAL(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDVAL(nv) N_VDestroy(nv)

#else

#define IDATYPE(fname) c_nvec_ida_ ## fname
#include <sundials/sundials_nvector.h>

#define WRAP_NVECTOR(v) NVEC_DATA(v)
#define RELINQUISH_WRAPPEDNV(v) {}

#define NVECTORIZE_VAL(v) NVEC_VAL(v)
#define RELINQUISH_NVECTORIZEDVAL(nv) {}

#endif

#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define IDATYPESTR(fname) QUOTE(IDATYPE(fname))

#define Val_none (Val_int(0))

static void errh(
	int error_code,
	const char *module,
	const char *func,
	char *msg,
	void *eh_data)
{
    CAMLparam0();
    CAMLlocal1(a);

    static value *ida_ml_errh;
    if (ida_ml_errh == NULL)
	ida_ml_errh = caml_named_value(IDATYPESTR(ida_ml_errh));

    a = caml_alloc_tuple(RECORD_IDA_ERROR_DETAILS_SIZE);
    Store_field(a, RECORD_IDA_ERROR_DETAILS_ERROR_CODE,
		Val_int(error_code));
    Store_field(a, RECORD_IDA_ERROR_DETAILS_MODULE_NAME,
		caml_copy_string(module));
    Store_field(a, RECORD_IDA_ERROR_DETAILS_FUNCTION_NAME,
		caml_copy_string(func));
    Store_field(a, RECORD_IDA_ERROR_DETAILS_ERROR_MESSAGE,
		caml_copy_string(msg));

    caml_callback2(*ida_ml_errh, Val_long((long int)eh_data), a);

    CAMLreturn0;
}

CAMLprim void IDATYPE(set_err_handler_fn)(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), errh, NULL);
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(clear_err_handler_fn)(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn0;
}

/* To be called on the return value of a user-supplied callback function.
 * Saves any exceptions before returning an appropriate return code.
 *
 * IDA requires the following callbacks from the user (some are optional):
 *  - the residual function
 *  - Jacobian function
 *  - the roots function
 * In IDA's C interface, the residual and Jacobian functions must return 0 for
 * success, 1 for recoverable error (meaning the integrator should retry with
 * different parameters), and -1 for unrecoverable error (meaning the
 * integrator should abort).  The roots function has a similar protocol except
 * that all nonzero exit code aborts the integrator.
 *
 * In the OCaml interface, these functions should instead return unit for
 * success, raise RecoverableFailure for recoverable error, and raise any other
 * exception for unrecoverable error.  The following function, callback_return,
 * implements this interface.  If the returned value is not an exception, it
 * simply returns 0.  It returns 1 if RecoverableFailure was raised and
 * check_recoverable is true.  Otherwise, stores the exception in the session
 * structure and returns -1.  In the last case, the stored exception should be
 * extracted and re-raised after execution leaves the integrator.
 */
enum cbr_flag {			/* Use this enum for the last parameter. */
    DONT_CHECK_RECOVERABLE = 0,
    CHECK_RECOVERABLE = 1
};
static int callback_return (value *session, value r,
			    enum cbr_flag check_recoverable)
{
    CAMLparam1(r);
    CAMLlocal1(exn);

    static value *recoverable_failure = NULL;

    /* The OCaml function may have recursively called the solver on the same
     * session instance, which would overwrite the user data.  We need to be
     * conservative and reinstate the user data here.
     *
     * FIXME: note there's no IDAGetUserData() so we can't rely on the
     * recursive call to re-establish the previous value of the user data
     * (right now the field is always reset to NULL).  However, perhaps we
     * could cast the mem pointer to IDAMem and access the user data field
     * directly.  I don't know if that kind of usage for IDAMem is supposed
     * to be legitimate, though.  */
    IDASetUserData (IDA_MEM_FROM_ML (*session), session);

    if (!Is_exception_result(r)) return 0;

    r = Extract_exception(r);

    if (check_recoverable) {
	if (recoverable_failure == NULL) {
	    recoverable_failure =
		caml_named_value("ida_RecoverableFailure");
	}
	if (Field(r, 0) == *recoverable_failure)
	    CAMLreturnT (int, 1);
    }

    /* Unrecoverable error.  Save the exception and return -1.  */
    exn = caml_alloc_small (1,0);
    Field (exn,0) = r;
    Store_field (*session, RECORD_IDA_SESSION_EXN_TEMP, exn);
    CAMLreturnT (int, -1);
}

static int resfn (realtype t, N_Vector y, N_Vector yp,
		  N_Vector resval, void *user_data)
{
    CAMLparam0 ();
    CAMLlocal1 (r);
    CAMLlocalN (args, 4);
    value *session = (value *)user_data;

    args[0] = caml_copy_double(t);
    args[1] = WRAP_NVECTOR (y);
    args[2] = WRAP_NVECTOR (yp);
    args[3] = WRAP_NVECTOR (resval);

    r = caml_callbackN_exn (Field (*session, RECORD_IDA_SESSION_RESFN),
			    sizeof (args) / sizeof (*args),
			    args);

    RELINQUISH_WRAPPEDNV (args[1]);
    RELINQUISH_WRAPPEDNV (args[2]);
    RELINQUISH_WRAPPEDNV (args[3]);

    CAMLreturnT (int, callback_return (session, r, CHECK_RECOVERABLE));
}

static value make_jac_arg(realtype t, realtype coef, N_Vector y, N_Vector yp,
			  N_Vector res, value tmp)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_IDA_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_COEF, caml_copy_double(coef));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_Y, WRAP_NVECTOR(y));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_YP, WRAP_NVECTOR(yp));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_RES, WRAP_NVECTOR(res));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_TMP, tmp);

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

static value make_double_tmp(N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, WRAP_NVECTOR(tmp1));
    Store_field(r, 1, WRAP_NVECTOR(tmp2));
    CAMLreturn(r);
}

static void relinquish_jac_arg(value arg)
{
    CAMLparam0();
    CAMLlocal1(tmp);
    int i;

    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_Y));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_YP));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_RES));

    tmp = Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_TMP);

    for (i = 0; i < Wosize_val (tmp); ++i)
	RELINQUISH_WRAPPEDNV(Field(tmp, i));

    CAMLreturn0;
}

static int jacfn (long int neq, realtype t, realtype coef,
		  N_Vector y, N_Vector yp, N_Vector res,
		  DlsMat jac, void *user_data,
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocal3 (arg, vjac, r);
    value *session = (value *)user_data;

    arg = make_jac_arg (t, coef, y, yp, res,
			make_triple_tmp (tmp1, tmp2, tmp3));
    vjac = caml_alloc_final (2, NULL, 0, 1);
    Store_field (vjac, 1, (value)jac);

    r = caml_callback2_exn (Field (*session, RECORD_IDA_SESSION_JACFN),
			    arg, vjac);

    relinquish_jac_arg (arg);

    CAMLreturn (callback_return (session, r, CHECK_RECOVERABLE));
}

static int bandjacfn (long int neq, long int mupper, long int mlower,
		      realtype t, realtype coef, N_Vector y, N_Vector yp,
		      N_Vector res, DlsMat jac, void *user_data,
		      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    value *session = (value *)user_data;
    caml_failwith ("to be implemented");
    IDASetUserData (IDA_MEM_FROM_ML (*session), user_data);
}

static int rootsfn (realtype t, N_Vector y, N_Vector yp,
		    realtype *gout, void *user_data)
{
    CAMLparam0 ();
    CAMLlocalN (args, 4);
    CAMLlocal4 (r, vy, vyp, roots);
    value *session = (value *)user_data;
    intnat nroots = Field (*session, RECORD_IDA_SESSION_NEQS);

    args[0] = caml_copy_double (t);
    args[1] = WRAP_NVECTOR (y);
    args[2] = WRAP_NVECTOR (yp);
    args[3] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    r = caml_callbackN_exn (Field (*session, RECORD_IDA_SESSION_ROOTSFN),
			    sizeof (args) / sizeof (*args),
			    args);

    RELINQUISH_WRAPPEDNV (args[1]);
    RELINQUISH_WRAPPEDNV (args[2]);

    CAMLreturn (callback_return (session, r, DONT_CHECK_RECOVERABLE));
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(y_d, ewt_d, r);

    value *session = user_data;

    y_d = WRAP_NVECTOR(y);
    ewt_d = WRAP_NVECTOR(ewt);

    r = caml_callback2_exn(IDA_ERRW_FROM_ML (*session), y_d, ewt_d);

    RELINQUISH_WRAPPEDNV(y_d);
    RELINQUISH_WRAPPEDNV(ewt_d);

    IDASetUserData (IDA_MEM_FROM_ML (*session), session);

    CAMLreturn(callback_return(session, r, DONT_CHECK_RECOVERABLE));
}

static int presetupfn(
    realtype t,
    N_Vector y,
    N_Vector yp,
    N_Vector res,
    realtype cj,
    void *user_data,
    N_Vector tmp1,
    N_Vector tmp2,
    N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal2(arg, r);
    value *session = user_data;

#if 0
    static value *cvode_ml_presetupfn;
    if (cvode_ml_presetupfn == NULL)
	cvode_ml_presetupfn = caml_named_value(CVTYPESTR(cvode_ml_presetupfn));
#endif

    arg = make_jac_arg(t, cj, y, yp, res, make_triple_tmp(tmp1, tmp2, tmp3));
    r = caml_callback_exn(Field (IDA_MEM_FROM_ML (*session),
				 RECORD_IDA_SESSION_PRESETUPFN),
			  arg);

    relinquish_jac_arg(arg);

    CAMLreturn(callback_return(session, r, CHECK_RECOVERABLE));
}

static int presolvefn(
	realtype t,
	N_Vector y,
	N_Vector yp,
	N_Vector res,
	N_Vector r,
	N_Vector z,
	realtype cj,
	realtype delta,
	void *user_data,
	N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal1(rv);
    CAMLlocalN(args, 4);
    value *session = user_data;

    args[0] = make_jac_arg(t, cj, y, yp, res, WRAP_NVECTOR(tmp));
    args[1] = WRAP_NVECTOR(r);
    args[2] = WRAP_NVECTOR(z);
    args[3] = caml_copy_double (delta);

    rv = caml_callbackN_exn(Field (IDA_MEM_FROM_ML (*session),
				   RECORD_IDA_SESSION_PRESOLVEFN),
			    sizeof (args) / sizeof (*args),
			    args);

    relinquish_jac_arg(args[0]);
    RELINQUISH_WRAPPEDNV(args[1]);
    RELINQUISH_WRAPPEDNV(args[2]);

    CAMLreturn(callback_return(session, rv, CHECK_RECOVERABLE));
}

static int jactimesfn(
    realtype t,
    N_Vector y,
    N_Vector yp,
    N_Vector res,
    N_Vector v,
    N_Vector Jv,
    realtype cj,
    void *user_data,
    N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal1(r);
    CAMLlocalN(args, 4);
    value *session = user_data;

    args[0] = make_jac_arg(t, cj, y, yp, res, make_double_tmp (tmp1, tmp2));
    args[1] = WRAP_NVECTOR(v);
    args[2] = WRAP_NVECTOR(Jv);

    r = caml_callbackN_exn(Field (IDA_MEM_FROM_ML (*session),
				  RECORD_IDA_SESSION_JACTIMESFN),
			   sizeof (args) / sizeof (*args),
			   args);

    relinquish_jac_arg(args[1]);
    RELINQUISH_WRAPPEDNV(args[2]);
    RELINQUISH_WRAPPEDNV(args[3]);

    CAMLreturn(callback_return(session, r, DONT_CHECK_RECOVERABLE));
}



CAMLprim void IDATYPE(wf_tolerances)(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDAWFtolerances(IDA_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("IDAWFtolerances", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_set_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetDenseJacFn(IDA_MEM_FROM_ML(vdata), jacfn);
    CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_clear_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetDenseJacFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_set_band_jac_fn)(value vdata, value fbandjacfn)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vdata), bandjacfn);
    CHECK_FLAG("IDADlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_clear_band_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDADlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(set_preconditioner)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDASpilsSetPreconditioner(IDA_MEM_FROM_ML(vdata),
					 presetupfn, presolvefn);
    CHECK_FLAG("IDASpilsSetPreconditioner", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(set_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDASpilsSetJacTimesVecFn(IDA_MEM_FROM_ML(vdata), jactimesfn);
    CHECK_FLAG("IDASpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(clear_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDASpilsSetJacTimesVecFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDASpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}


/* IDAInit + IDARootInit + linear solver setup.  The residual and root
 * functions are set to C stubs; the actual functions should be assigned to the
 * ida_session record to be created on the OCaml side.  */
CAMLprim value IDATYPE (init) (value linsolver, value vy, value vyp,
			       value vneqs, value vnroots, value vt0)
{
    CAMLparam5 (linsolver, vy, vyp, vneqs, vnroots);
    CAMLxparam1 (vt0);
    CAMLlocal1 (r);
    int flag;
    N_Vector y, yp;
    int neqs = Int_val (vneqs);
    int nroots = Int_val (vnroots);

    void *ida_mem = IDACreate ();
    if (ida_mem == NULL) {
	caml_failwith ("IDACreate failed");
    }

    y = NVECTORIZE_VAL (vy);
    yp = NVECTORIZE_VAL (vyp);
    flag = IDAInit (ida_mem, resfn, Double_val (vt0), y, yp);
    RELINQUISH_NVECTORIZEDVAL (y);
    RELINQUISH_NVECTORIZEDVAL (yp);
    CHECK_FLAG ("IDAInit", flag);

    flag = IDARootInit (ida_mem, nroots, rootsfn);
    CHECK_FLAG ("IDARootInit", flag);

    ida_ml_set_linear_solver (ida_mem, linsolver, neqs);

    flag = IDASStolerances (ida_mem, RCONST (1.0e-4), RCONST (1.0e-8));
    CHECK_FLAG ("IDASStolerances", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, (value)ida_mem);
    Store_field(r, 1, (value)NULL); // no err_file = NULL

    CAMLreturn (r);
}

CAMLprim value IDATYPE(init_bytecode) (value args[], int n)
{
    return IDATYPE(init)(args[0], args[1], args[2], args[3], args[4], args[5]);
}

CAMLprim void IDATYPE(sv_tolerances) (value ida_mem, value vrtol, value vavtol)
{
    CAMLparam3 (ida_mem, vrtol, vavtol);
    N_Vector avtol;
    int flag;

    avtol = NVECTORIZE_VAL (vavtol);
    flag = IDASVtolerances (IDA_MEM_FROM_ML (ida_mem),
			    Double_val (vrtol), avtol);
    RELINQUISH_NVECTORIZEDVAL (avtol);
    CHECK_FLAG ("IDASVtolerances", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(reinit)(value vdata, value t0, value y0, value yp0)
{
    CAMLparam4(vdata, t0, y0, yp0);

    N_Vector y0_nv = NVECTORIZE_VAL(y0);
    N_Vector yp0_nv = NVECTORIZE_VAL(yp0);
    int flag = IDAReInit(IDA_MEM_FROM_ML(vdata), Double_val(t0), y0_nv, yp0_nv);
    RELINQUISH_NVECTORIZEDVAL(yp0_nv);
    RELINQUISH_NVECTORIZEDVAL(y0_nv);
    CHECK_FLAG("IDAReInit", flag);

    CAMLreturn0;
}

static value solve (value vdata, value nextt, value vy, value vyp, int onestep)
{
    CAMLparam4 (vdata, nextt, vy, vyp);
    CAMLlocal1 (ret);
    void *ida_mem = IDA_MEM_FROM_ML (vdata);
    realtype tret;
    int flag;
    N_Vector y, yp;
    enum ida_solver_result_tag result;

    IDASetUserData (ida_mem, &vdata);

    y = NVECTORIZE_VAL (vy);
    yp = NVECTORIZE_VAL (vyp);
    flag = IDASolve (ida_mem, Double_val (nextt), &tret, y, yp,
	             onestep ? IDA_ONE_STEP : IDA_NORMAL);
    RELINQUISH_NVECTORIZEDVAL (y);
    RELINQUISH_NVECTORIZEDVAL (yp);

    /* For precaution; if we screw up somewhere, we'll promptly segfault
     * rather than follow a dangling pointer.  */
    IDASetUserData (ida_mem, NULL);

    switch (flag) {
    case IDA_SUCCESS:
	result = VARIANT_IDA_SOLVER_RESULT_CONTINUE;
	break;

    case IDA_ROOT_RETURN:
	result = VARIANT_IDA_SOLVER_RESULT_ROOTSFOUND;
	break;

    case IDA_TSTOP_RETURN:
	result = VARIANT_IDA_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    case IDA_RES_FAIL:
	ret = Field (vdata, RECORD_IDA_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_IDA_SESSION_EXN_TEMP, Val_none);
	    /* FIXME: In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	/*FALLTHROUGH*/
    default:
	CHECK_FLAG ("IDASolve", flag);
    }

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}


CAMLprim value IDATYPE(normal)(value vdata, value nextt, value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 0));
}

CAMLprim value IDATYPE(one_step)(value vdata, value nextt, value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 1));
}

CAMLprim void IDATYPE(get_dky)(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVECTORIZE_VAL(vy);

    int flag = IDAGetDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
			 Int_val(vk), y_nv);
    CHECK_FLAG("IDAGetDky", flag);
    RELINQUISH_NVECTORIZEDVAL(y_nv);

    CAMLreturn0;
}

CAMLprim void IDATYPE(get_err_weights)(value vida_mem, value verrws)
{
    CAMLparam2(vida_mem, verrws);

    N_Vector errws_nv = NVECTORIZE_VAL(verrws);

    int flag = IDAGetErrWeights(IDA_MEM_FROM_ML(vida_mem), errws_nv);
    RELINQUISH_NVECTORIZEDVAL(errws_nv);
    CHECK_FLAG("IDAGetErrWeights", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(get_est_local_errors)(value vida_mem, value vele)
{
    CAMLparam2(vida_mem, vele);

    N_Vector ele_nv = NVECTORIZE_VAL(vele);

    int flag = IDAGetEstLocalErrors(IDA_MEM_FROM_ML(vida_mem), ele_nv);
    RELINQUISH_NVECTORIZEDVAL(ele_nv);
    CHECK_FLAG("IDAGetEstLocalErrors", flag);

    CAMLreturn0;
}

static void calc_ic (void *ida_mem, value *session, int icopt, realtype tout1)
{
    CAMLparam0 ();
    CAMLlocal1 (exn);
    int flag;

    IDASetUserData (ida_mem, session);

    flag = IDACalcIC (ida_mem, icopt, tout1);

    /* For precaution; if we screw up somewhere, we'll promptly segfault
     * instead of accessing a dangling pointer.  */
    IDASetUserData (ida_mem, NULL);

    if (flag < 0) {
	switch (flag) {
	case IDA_RES_FAIL: case IDA_LSETUP_FAIL: case IDA_LSOLVE_FAIL:
	    /* If an exception is saved in the session, grab it and re-raise. */
	    exn = Field (*session, RECORD_IDA_SESSION_EXN_TEMP);
	    if (Is_block (exn)) {
		Store_field (*session, RECORD_IDA_SESSION_EXN_TEMP, Val_none);
		/* FIXME: In bytecode, caml_raise() duplicates some parts of
		 * the stacktrace.  This does not seem to happen in native code
		 * execution.  */
		caml_raise (Field (exn, 0));
	    }
	    /* Otherwise, raise the generic exception Ida.ResFuncFailure  */
	    /*FALLTHROUGH*/
	default:
	    CHECK_FLAG ("IDACalcIC", flag);
	}
    }
    CAMLreturn0;
}

CAMLprim void IDATYPE(calc_ic_y)(value vida_mem, value tout1)
{
    CAMLparam2 (vida_mem, tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    calc_ic (ida_mem, &vida_mem, IDA_Y_INIT, Double_val (tout1));

    CAMLreturn0;
}

CAMLprim void IDATYPE(calc_ic_ya_ydp)(value vida_mem, value vid, value tout1)
{
    CAMLparam2 (vida_mem, tout1);
    int flag;
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    N_Vector id = NVECTORIZE_VAL (vid);
    flag = IDASetId (ida_mem, id);
    RELINQUISH_NVECTORIZEDVAL (id);
    CHECK_FLAG ("IDASetId", flag);

    calc_ic (ida_mem, &vida_mem, IDA_YA_YDP_INIT, Double_val (tout1));
    CAMLreturn0;
}

CAMLprim void IDATYPE(set_constraints) (value vida_mem, value vconstraints)
{
    CAMLparam2(vida_mem, vconstraints);
    int flag;

    N_Vector constraints = NVECTORIZE_VAL (vconstraints);
    flag = IDASetConstraints (IDA_MEM_FROM_ML (vida_mem), constraints);
    RELINQUISH_NVECTORIZEDVAL (constraints);
    CHECK_FLAG ("IDASetConstraints", flag);

    CAMLreturn0;
}
