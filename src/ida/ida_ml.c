/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include "../config.h"

/* Sundials IDA interface functions that do not involve NVectors. */

#include <caml/memory.h>
#include <caml/fail.h>
#include <caml/alloc.h>
#include <caml/callback.h>
#include <caml/bigarray.h>

#ifdef SUNDIALSML_WITHSENS
/* IDAS (with sensitivity) */

#include <idas/idas.h>
#include <idas/idas_impl.h>

/* linear solvers */
#if   400 <= SUNDIALS_LIB_VERSION
#include <idas/idas_ls.h>
#elif 300 <= SUNDIALS_LIB_VERSION
#include <idas/idas_direct.h>
#include <idas/idas_spils.h>
#else
#include <idas/idas_dense.h>
#include <idas/idas_band.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_spbcgs.h>
#endif

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <idas/idas_lapack.h>
#endif

#else  /* IDA (without sensitivity) */

#include <ida/ida.h>
#include <ida/ida_impl.h>

/* linear solvers */
#if   400 <= SUNDIALS_LIB_VERSION
#include <ida/ida_ls.h>
#elif 300 <= SUNDIALS_LIB_VERSION
#include <ida/ida_direct.h>
#include <ida/ida_spils.h>
#else
#include <ida/ida_dense.h>
#include <ida/ida_band.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_sptfqmr.h>
#include <ida/ida_spbcgs.h>
#endif

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <ida/ida_lapack.h>
#endif

#endif

#include <stdio.h>

#define MAX_ERRMSG_LEN 256

#include "ida_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"
#include "../lsolvers/sundials_nonlinearsolver_ml.h"
#include "../lsolvers/sundials_linearsolver_ml.h"


CAMLprim value sunml_ida_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (IDA, exns);
    CAMLreturn (Val_unit);
}

int sunml_ida_translate_exception(value session, value exn,
			    recoverability recoverable)
{
    CAMLparam2(session, exn);
    CAMLlocal1(bucket);

    exn = Extract_exception(exn);

    if (recoverable && Field(exn, 0) == SUNDIALS_EXN_TAG (RecoverableFailure))
	CAMLreturnT (int, 1);

    /* Unrecoverable error.  Save the exception and return -1.  */
    bucket = caml_alloc_small (1,0);
    Field (bucket, 0) = exn;
    Store_field (session, RECORD_IDA_SESSION_EXN_TEMP, bucket);
    CAMLreturnT (int, -1);
}

/* callbacks */

static void errh(
	int error_code,
	const char *module,
	const char *func,
	char *msg,
	void *eh_data)
{
    CAMLparam0();
    CAMLlocal2(session, a);
    value *backref = eh_data;

    a = caml_alloc_tuple(RECORD_SUNDIALS_ERROR_DETAILS_SIZE);
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_CODE,
		Val_int(error_code));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_MODULE_NAME,
		caml_copy_string(module));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_FUNCTION_NAME,
		caml_copy_string(func));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_MESSAGE,
		caml_copy_string(msg));

    WEAK_DEREF (session, *backref);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (Field(session, RECORD_IDA_SESSION_ERRH), a);
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined error handler");

    CAMLreturn0;
}

CAMLprim value sunml_ida_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), errh,
				  IDA_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

static int resfn (realtype t, N_Vector y, N_Vector yp,
		  N_Vector resval, void *user_data)
{
    CAMLparam0 ();
    CAMLlocalN (args, 4);
    CAMLlocal2 (session, cb);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK (y);
    args[2] = NVEC_BACKLINK (yp);
    args[3] = NVEC_BACKLINK (resval);

    WEAK_DEREF (session, *(value*)user_data);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (IDA_RESFN_FROM_ML (session), 4, args);

    CAMLreturnT (int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

value sunml_ida_make_jac_arg(realtype t, realtype coef, N_Vector y, N_Vector yp,
		       N_Vector res, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_IDA_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_COEF, caml_copy_double(coef));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_YP, NVEC_BACKLINK(yp));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_RES, NVEC_BACKLINK(res));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

value sunml_ida_make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(3);
    Store_field(r, 0, NVEC_BACKLINK(tmp1));
    Store_field(r, 1, NVEC_BACKLINK(tmp2));
    Store_field(r, 2, NVEC_BACKLINK(tmp3));
    CAMLreturn(r);
}

value sunml_ida_make_double_tmp(N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, NVEC_BACKLINK(tmp1));
    Store_field(r, 1, NVEC_BACKLINK(tmp2));
    CAMLreturn(r);
}

#if 300 <= SUNDIALS_LIB_VERSION
static int jacfn (realtype t,
		  realtype coef,
		  N_Vector y,
		  N_Vector yp,
		  N_Vector res,
		  SUNMatrix jac,
		  void *user_data,
		  N_Vector tmp1,
		  N_Vector tmp2,
		  N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 2);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_ida_make_jac_arg (t, coef, y, yp, res,
				sunml_ida_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = MAT_BACKLINK(jac);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#else

/* Dense and band Jacobians only work with serial NVectors.  */
static int jacfn (long int neq,
		  realtype t,
		  realtype coef,
		  N_Vector y,
		  N_Vector yp,
		  N_Vector res,
		  DlsMat jac,
		  void *user_data,
		  N_Vector tmp1,
		  N_Vector tmp2,
		  N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(jac));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_ida_make_jac_arg (t, coef, y, yp, res,
				sunml_ida_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bandjacfn (long int neq,
		      long int mupper,
		      long int mlower,
		      realtype t,
		      realtype coef,
		      N_Vector y,
		      N_Vector yp,
		      N_Vector res,
		      DlsMat jac,
		      void *user_data,
		      N_Vector tmp1,
		      N_Vector tmp2,
		      N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(jac));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_ida_make_jac_arg (t, coef, y, yp, res,
				sunml_ida_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#endif

static int rootsfn (realtype t, N_Vector y, N_Vector yp,
		    realtype *gout, void *user_data)
{
    CAMLparam0 ();
    CAMLlocal1 (session);
    CAMLlocalN (args, 4);
    intnat nroots;

    WEAK_DEREF (session, *(value*)user_data);

    nroots = IDA_NROOTS_FROM_ML (session);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);
    args[2] = NVEC_BACKLINK (yp);
    args[3] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (IDA_ROOTSFN_FROM_ML (session), 4, args);

    CAMLreturnT (int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    value *backref = user_data;

    WEAK_DEREF (session, *backref);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (Field (session, RECORD_IDA_SESSION_ERRW),
				  NVEC_BACKLINK (y), NVEC_BACKLINK (ewt));
    if (Is_exception_result (r)) {
	r = Extract_exception (r);
	if (Field (r, 0) != SUNDIALS_EXN_TAG (NonPositiveEwt))
	    sunml_warn_discarded_exn (r, "user-defined error weight fun");
	CAMLreturnT (int, -1);
    }

    CAMLreturnT (int, 0);
}

static int precsetupfn(realtype t,
		       N_Vector y,
		       N_Vector yp,
		       N_Vector res,
		       realtype cj,
		       void *user_data
#if SUNDIALS_LIB_VERSION < 300
		       ,
		       N_Vector tmp1,
		       N_Vector tmp2,
		       N_Vector tmp3
#endif
		       )
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDA_SPILS_PRECFNS_PREC_SETUP_FN);
    cb = Field (cb, 0);

    arg = sunml_ida_make_jac_arg(t, cj, y, yp, res, Val_unit);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, arg);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int precsolvefn(
	realtype t,
	N_Vector y,
	N_Vector yp,
	N_Vector res,
	N_Vector rvec,
	N_Vector z,
	realtype cj,
	realtype delta,
	void *user_data
#if SUNDIALS_LIB_VERSION < 300
	,
	N_Vector tmp
#endif
	)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal2(session, cb);

    args[0] = sunml_ida_make_jac_arg(t, cj, y, yp, res, Val_unit);
    args[1] = NVEC_BACKLINK (rvec);
    args[2] = NVEC_BACKLINK (z);
    args[3] = caml_copy_double (delta);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDA_SPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    CAMLreturnT (int, CHECK_EXCEPTION (session, r, RECOVERABLE));
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
    CAMLlocalN(args, 3);
    CAMLlocal2(session, cb);

    args[0] = sunml_ida_make_jac_arg (t, cj, y, yp, res,
				sunml_ida_make_double_tmp (tmp1, tmp2));
    args[1] = NVEC_BACKLINK (v);
    args[2] = NVEC_BACKLINK (Jv);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    CAMLreturnT (int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 300
static int jacsetupfn(realtype t,
		      N_Vector y,
		      N_Vector yp,
		      N_Vector res,
		      realtype cj,
		      void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = sunml_ida_make_jac_arg(t, cj, y, yp, res, Val_unit);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, arg);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_ida_dls_dense (value vida_mem, value vneqs, value vset_jac)
{
    CAMLparam3(vida_mem, vneqs, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDADense (ida_mem, neqs);
    CHECK_FLAG ("IDADense", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetDenseJacFn(ida_mem, jacfn);
	CHECK_DLS_FLAG("IDADlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_dls_lapack_dense (value vida_mem, value vneqs,
				       value vset_jac)
{
    CAMLparam3 (vida_mem, vneqs, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDALapackDense (ida_mem, neqs);
    CHECK_FLAG ("IDALapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetDenseJacFn (ida_mem, jacfn);
	CHECK_DLS_FLAG("IDADlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_dls_band (value vida_mem, value vneqs,
			       value mupper, value mlower, value vset_jac)
{
    CAMLparam5(vida_mem, vneqs, mupper, mlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDABand (ida_mem, neqs, Long_val (mupper), Long_val (mlower));
    CHECK_FLAG ("IDABand", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetBandJacFn(ida_mem, bandjacfn);
	CHECK_DLS_FLAG("IDADlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_dls_lapack_band (value vida_mem, value vneqs,
				      value mupper, value mlower,
				      value vset_jac)
{
    CAMLparam5(vida_mem, vneqs, mupper, mlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDALapackBand (ida_mem, neqs, Long_val (mupper), Long_val (mlower));
    CHECK_FLAG ("IDALapackBand", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetBandJacFn(ida_mem, bandjacfn);
	CHECK_DLS_FLAG("IDADlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_linear_solver (value vida_mem, value vlsolv,
					      value vojmat, value vhasjac)
{
    CAMLparam4(vida_mem, vlsolv, vojmat, vhasjac);
#if 400 <= SUNDIALS_LIB_VERSION
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = IDASetLinearSolver(ida_mem, lsolv, jmat);
    CHECK_LS_FLAG ("IDASetLinearSolver", flag);
    if (Bool_val (vhasjac)) {
	flag = IDASetJacFn(ida_mem, jacfn);
	CHECK_LS_FLAG("IDASetJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_dls_set_linear_solver (value vida_mem, value vlsolv,
					    value vjmat, value vhasjac)
{
    CAMLparam4(vida_mem, vlsolv, vjmat, vhasjac);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = IDADlsSetLinearSolver(ida_mem, lsolv, jmat);
    CHECK_DLS_FLAG ("IDADlsSetLinearSolver", flag);
    if (Bool_val (vhasjac)) {
	flag = IDADlsSetJacFn(ida_mem, jacfn);
	CHECK_DLS_FLAG("IDADlsSetJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_set_linear_solver (value vida_mem, value vlsolv)
{
    CAMLparam2(vida_mem, vlsolv);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);

    int flag = IDASpilsSetLinearSolver(ida_mem, lsolv);
    CHECK_SPILS_FLAG ("IDASpilsSetLinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_set_preconditioner (value vsession,
					       value vset_presetup)
{
    CAMLparam2 (vsession, vset_presetup);
    void *mem = IDA_MEM_FROM_ML (vsession);
#if 400 <= SUNDIALS_LIB_VERSION
    IDALsPrecSetupFn setup = Bool_val (vset_presetup) ? precsetupfn : NULL;
    int flag = IDASetPreconditioner (mem, setup, precsolvefn);
    CHECK_LS_FLAG ("IDASetPreconditioner", flag);
#else
    IDASpilsPrecSetupFn setup = Bool_val (vset_presetup) ? precsetupfn : NULL;
    int flag = IDASpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_SPILS_FLAG ("IDASpilsSetPreconditioner", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_set_jac_times(value vdata, value vhas_setup,
						      value vhas_times)
{
    CAMLparam3(vdata, vhas_setup, vhas_times);
#if 400 <= SUNDIALS_LIB_VERSION
    IDALsJacTimesSetupFn setup = Bool_val (vhas_setup) ? jacsetupfn : NULL;
    IDALsJacTimesVecFn   times = Bool_val (vhas_times) ? jactimesfn : NULL;

    int flag = IDASetJacTimes(IDA_MEM_FROM_ML(vdata), setup, times);
    CHECK_LS_FLAG("IDASetJacTimes", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    IDASpilsJacTimesSetupFn setup = Bool_val (vhas_setup) ? jacsetupfn : NULL;
    IDASpilsJacTimesVecFn   times = Bool_val (vhas_times) ? jactimesfn : NULL;

    int flag = IDASpilsSetJacTimes(IDA_MEM_FROM_ML(vdata), setup, times);
    CHECK_SPILS_FLAG("IDASpilsSetJacTimes", flag);
#else
    IDASpilsJacTimesVecFn jac = Bool_val (vhas_times) ? jactimesfn : NULL;
    int flag = IDASpilsSetJacTimesVecFn(IDA_MEM_FROM_ML(vdata), jac);
    CHECK_SPILS_FLAG("IDASpilsSetJacTimesVecFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_spgmr (value vida_mem, value vmaxl)
{
    CAMLparam2 (vida_mem, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASpgmr (ida_mem, Int_val (vmaxl));
    CHECK_FLAG ("IDASpgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_spbcgs (value vida_mem, value vmaxl)
{
    CAMLparam2 (vida_mem, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASpbcg (ida_mem, Int_val (vmaxl));
    CHECK_FLAG ("IDASpbcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_sptfqmr (value vida_mem, value vmaxl)
{
    CAMLparam2 (vida_mem, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASptfqmr (ida_mem, Int_val (vmaxl));
    CHECK_FLAG ("IDASptfqmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_nonlinear_solver(value vida_mem, value vnlsolv)
{
    CAMLparam2(vida_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *ida_mem = IDA_MEM_FROM_ML(vida_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    int flag = IDASetNonlinearSolver(ida_mem, nlsolv);
    CHECK_FLAG ("IDASetNonlinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_wf_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDAWFtolerances(IDA_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("IDAWFtolerances", flag);

    CAMLreturn (Val_unit);
}

/* Sets the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim value sunml_ida_root_init (value vida_mem, value vnroots)
{
    CAMLparam2 (vida_mem, vnroots);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int nroots = Int_val (vnroots);
    int flag;
    flag = IDARootInit (ida_mem, nroots, rootsfn);
    CHECK_FLAG ("IDARootInit", flag);
    Store_field (vida_mem, RECORD_IDA_SESSION_NROOTS, vnroots);
    CAMLreturn (Val_unit);
}

/* IDACreate + IDAInit + IDASetUser.  The vy and vyp vectors must have the same
 * length (not checked).  Returns an IDAMem (= void*) along with a global root
 * wrapping weakref that is attached to the IDAMem.  Before doing the remaining
 * initialization, these pointers have to be wrapped in a session so that if
 * initialization fails the GC frees IDAMem and the weak global root.  Doing
 * that cleanup in C would be ugly, e.g. we wouldn't be able to reuse
 * c_ida_set_linear_solver() so we have to duplicate it or hack some ad-hoc
 * extensions to that function.  */
CAMLprim value sunml_ida_init (value weakref, value vt0, value vy, value vyp)
{
    CAMLparam4 (weakref, vy, vyp, vt0);
    CAMLlocal2 (r, vida_mem);
    int flag;
    N_Vector y, yp;
    void *ida_mem = NULL;
    value *backref = NULL;

    ida_mem = IDACreate ();
    if (ida_mem == NULL)
	caml_failwith ("IDACreate failed");

    vida_mem = caml_alloc_final(1, NULL, 1, 5);
    IDA_MEM(vida_mem) = ida_mem;

    y = NVEC_VAL (vy);
    yp = NVEC_VAL (vyp);
    flag = IDAInit (ida_mem, resfn, Double_val (vt0), y, yp);
    if (flag != IDA_SUCCESS) {
	IDAFree (&ida_mem);
	CHECK_FLAG ("IDAInit", flag);
    }

    backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	IDAFree (&ida_mem);
	caml_failwith ("Out of memory");
    }
    IDASetUserData (ida_mem, backref);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, vida_mem);
    Store_field(r, 1, (value)backref);

    CAMLreturn (r);
}

CAMLprim value sunml_ida_sv_tolerances (value ida_mem, value vrtol, value vavtol)
{
    CAMLparam3 (ida_mem, vrtol, vavtol);
    N_Vector avtol;
    int flag;

    avtol = NVEC_VAL (vavtol);
    flag = IDASVtolerances (IDA_MEM_FROM_ML (ida_mem),
			    Double_val (vrtol), avtol);
    CHECK_FLAG ("IDASVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_reinit(value vdata, value t0, value y0, value yp0)
{
    CAMLparam4(vdata, t0, y0, yp0);

    N_Vector y0_nv = NVEC_VAL(y0);
    N_Vector yp0_nv = NVEC_VAL(yp0);
    int flag = IDAReInit(IDA_MEM_FROM_ML(vdata), Double_val(t0),
			 y0_nv, yp0_nv);
    CHECK_FLAG("IDAReInit", flag);

    CAMLreturn (Val_unit);
}

static value solve (value vdata, value nextt, value vy, value vyp, int onestep)
{
    CAMLparam4 (vdata, nextt, vy, vyp);
    CAMLlocal1 (ret);
    void *ida_mem = IDA_MEM_FROM_ML (vdata);
    realtype tret;
    int flag;
    N_Vector y, yp;
    enum ida_solver_result_tag result = -1;

    y = NVEC_VAL (vy);
    yp = NVEC_VAL (vyp);
    flag = IDASolve (ida_mem, Double_val (nextt), &tret, y, yp,
	             onestep ? IDA_ONE_STEP : IDA_NORMAL);

    switch (flag) {
    case IDA_SUCCESS:
	result = VARIANT_IDA_SOLVER_RESULT_SUCCESS;
	break;

    case IDA_ROOT_RETURN:
	result = VARIANT_IDA_SOLVER_RESULT_ROOTSFOUND;
	break;

    case IDA_TSTOP_RETURN:
	result = VARIANT_IDA_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_IDA_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_IDA_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	CHECK_FLAG ("IDASolve", flag);
    }

    assert (Field (vdata, RECORD_IDA_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}


CAMLprim value sunml_ida_solve_normal (value vdata, value nextt,
				   value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 0));
}

CAMLprim value sunml_ida_solve_one_step (value vdata, value nextt,
				     value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 1));
}

CAMLprim value sunml_ida_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVEC_VAL(vy);

    int flag = IDAGetDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
			 Int_val(vk), y_nv);
    CHECK_FLAG("IDAGetDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_get_err_weights(value vida_mem, value verrws)
{
    CAMLparam2(vida_mem, verrws);

    N_Vector errws_nv = NVEC_VAL(verrws);

    int flag = IDAGetErrWeights(IDA_MEM_FROM_ML(vida_mem), errws_nv);
    CHECK_FLAG("IDAGetErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_get_est_local_errors(value vida_mem, value vele)
{
    CAMLparam2(vida_mem, vele);

    N_Vector ele_nv = NVEC_VAL(vele);

    int flag = IDAGetEstLocalErrors(IDA_MEM_FROM_ML(vida_mem), ele_nv);
    CHECK_FLAG("IDAGetEstLocalErrors", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value sunml_ida_set_id (value vida_mem, value vid)
{
    CAMLparam2(vida_mem, vid);
    N_Vector id;

    id = NVEC_VAL (vid);
    int flag = IDASetId (IDA_MEM_FROM_ML(vida_mem), id);
    CHECK_FLAG("IDASetId", flag);

    CAMLreturn (Val_unit);
}


static void calc_ic (void *ida_mem, value session, int icopt, realtype tout1,
		     value vy, value vyp)
{
    CAMLparam3 (session, vy, vyp);
    CAMLlocal1 (exn);
    int flag;
    N_Vector y, yp;

    flag = IDACalcIC (ida_mem, icopt, tout1);

    if (flag < 0) {
	/* If an exception is saved in the session, grab it and re-raise. */
	exn = Field (session, RECORD_IDA_SESSION_EXN_TEMP);
	if (Is_block (exn)) {
	    Store_field (session, RECORD_IDA_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (exn, 0));
	}
	/* Otherwise, raise a generic exception like Ida.ResFuncFailure */
	CHECK_FLAG ("IDACalcIC", flag);
    }

    /* Retrieve the calculated initial conditions if y,yp are given.  */
    y  = (Is_block (vy))  ? NVEC_VAL (Field (vy, 0))  : NULL;
    yp = (Is_block (vyp)) ? NVEC_VAL (Field (vyp, 0)) : NULL;
    if (y != NULL || yp != NULL) {
	flag = IDAGetConsistentIC (ida_mem, y, yp);
	CHECK_FLAG ("IDAGetConsistentIC", flag);
    }
    CAMLreturn0;
}

CAMLprim value sunml_ida_calc_ic_y(value vida_mem, value vy, value tout1)
{
    CAMLparam3 (vida_mem, vy, tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    calc_ic (ida_mem, vida_mem, IDA_Y_INIT, Double_val (tout1), vy, Val_none);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_calc_ic_ya_ydp(value vida_mem, value y, value yp,
				    value tout1)
{
    CAMLparam4 (vida_mem, y, yp, tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    calc_ic (ida_mem, vida_mem, IDA_YA_YDP_INIT, Double_val (tout1), y, yp);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_constraints (value vida_mem, value vconstraints)
{
    CAMLparam2(vida_mem, vconstraints);
    int flag;

    N_Vector constraints = NVEC_VAL (vconstraints);
    flag = IDASetConstraints (IDA_MEM_FROM_ML (vida_mem), constraints);
    CHECK_FLAG ("IDASetConstraints", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_clear_constraints (value vida_mem)
{
    CAMLparam1(vida_mem);
    int flag;

    flag = IDASetConstraints (IDA_MEM_FROM_ML (vida_mem), NULL);
    CHECK_FLAG ("IDASetConstraints", flag);

    CAMLreturn (Val_unit);
}

void sunml_ida_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == IDA_SUCCESS
	|| flag == IDA_ROOT_RETURN
	|| flag == IDA_TSTOP_RETURN) return;

    switch (flag) {
    case IDA_ILL_INPUT:
	caml_raise_constant(IDA_EXN(IllInput));

    case IDA_CONV_FAIL:
	caml_raise_constant(IDA_EXN(ConvergenceFailure));

    case IDA_TOO_MUCH_WORK:
	caml_raise_constant(IDA_EXN(TooMuchWork));

    case IDA_TOO_MUCH_ACC:
	caml_raise_constant(IDA_EXN(TooMuchAccuracy));

    case IDA_LINIT_FAIL:
	caml_raise_constant(IDA_EXN(LinearInitFailure));

    case IDA_LSETUP_FAIL:
	caml_raise_constant(IDA_EXN(LinearSetupFailure));

    case IDA_LSOLVE_FAIL:
	caml_raise_constant(IDA_EXN(LinearSolveFailure));

    case IDA_BAD_EWT:
	caml_raise_constant(IDA_EXN(BadEwt));

    case IDA_NO_RECOVERY:
	caml_raise_constant(IDA_EXN(NoRecovery));

    case IDA_RES_FAIL:
	caml_raise_constant(IDA_EXN(ResFuncFailure));

    case IDA_FIRST_RES_FAIL:
	caml_raise_constant(IDA_EXN(FirstResFuncFailure));

    case IDA_REP_RES_ERR:
	caml_raise_constant(IDA_EXN(RepeatedResFuncFailure));

    case IDA_RTFUNC_FAIL:
	caml_raise_constant(IDA_EXN(RootFuncFailure));

    case IDA_CONSTR_FAIL:
	caml_raise_constant(IDA_EXN(ConstraintFailure));

    case IDA_BAD_K:
	caml_raise_constant(IDA_EXN(BadK));

    case IDA_BAD_T:
	caml_raise_constant(IDA_EXN(BadT));

#if SUNDIALS_LIB_VERSION >= 400
    case IDA_NLS_INIT_FAIL:
	caml_raise_constant(IDA_EXN(NonlinearInitFailure));

    case IDA_NLS_SETUP_FAIL:
	caml_raise_constant(IDA_EXN(NonlinearSetupFailure));

    case IDA_VECTOROP_ERR:
	caml_raise_constant(IDA_EXN(VectorOpErr));
#endif

    default:
	/* e.g. IDA_MEM_NULL, IDA_MEM_FAIL */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 IDAGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

#if 400 <= SUNDIALS_LIB_VERSION
void sunml_ida_check_ls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == IDALS_SUCCESS) return;

    switch (flag) {
	case IDALS_ILL_INPUT:
	    caml_raise_constant(IDA_EXN(IllInput));

	case IDALS_MEM_FAIL:
	    caml_raise_out_of_memory();

	case IDALS_SUNMAT_FAIL:
	case IDALS_SUNLS_FAIL:
	case IDALS_JACFUNC_UNRECVR:
	case IDALS_JACFUNC_RECVR:
	default:
	    /* e.g. IDALS_MEM_NULL, IDALS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    IDAGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#else
void sunml_ida_check_dls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == IDADLS_SUCCESS) return;

    switch (flag) {
	case IDADLS_ILL_INPUT:
	    caml_raise_constant(IDA_EXN(IllInput));

	case IDADLS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if SUNDIALS_LIB_VERSION >= 300
	case IDADLS_SUNMAT_FAIL:
#endif
	case IDADLS_JACFUNC_UNRECVR:
	case IDADLS_JACFUNC_RECVR:
	default:
	    /* e.g. IDADLS_MEM_NULL, IDADLS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    IDADlsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

void sunml_ida_check_spils_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == IDASPILS_SUCCESS) return;

    switch (flag) {
	case IDASPILS_ILL_INPUT:
	    caml_raise_constant(IDA_EXN(IllInput));

	case IDASPILS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if SUNDIALS_LIB_VERSION >+ 300
	case IDASPILS_SUNLS_FAIL:
#endif
	default:
	    /* e.g. IDASPILS_MEM_NULL, IDASPILS_PMEM_NULL, IDASPILS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     IDASpilsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#endif

CAMLprim value sunml_ida_session_finalize(value vdata)
{
    if (IDA_MEM_FROM_ML(vdata) != NULL) {
	void *ida_mem = IDA_MEM_FROM_ML(vdata);
	value *backref = IDA_BACKREF_FROM_ML(vdata);
	IDAFree(&ida_mem);
	sunml_sundials_free_value(backref);
    }

    return Val_unit;
}
 
CAMLprim value sunml_ida_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDASStolerances(IDA_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("IDASStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < IDA_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = IDAGetRootInfo(IDA_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("IDAGetRootInfo", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_get_integrator_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

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

    flag = IDAGetIntegratorStats(IDA_MEM_FROM_ML(vdata),
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
    CHECK_FLAG("IDAGetIntegratorStats", flag);

    r = caml_alloc_tuple(RECORD_IDA_INTEGRATOR_STATS_SIZE);
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_RES_EVALS, Val_long(nfevals));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS, Val_long(nlinsetups));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_ERROR_TEST_FAILURES, Val_long(netfails));

    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_LAST_INTERNAL_ORDER, Val_int(qlast));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER, Val_int(qcur));

    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_INITIAL_STEP_SIZE, caml_copy_double(hinused));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_LAST_STEP_SIZE, caml_copy_double(hlast));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_NEXT_STEP_SIZE, caml_copy_double(hcur));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_INTERNAL_TIME, caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim value sunml_ida_set_error_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

    int flag = IDASetErrFile(IDA_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("IDASetErrFile", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < IDA_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = IDASetRootDirection(IDA_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("IDASetRootDirection", flag);

    CAMLreturn (Val_unit);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

CAMLprim value sunml_ida_get_work_space(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = IDAGetWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrw, &leniw);
    CHECK_FLAG("IDAGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
}


CAMLprim value sunml_ida_get_num_steps(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumSteps(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_ida_get_num_res_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumResEvals(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumResEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_ida_get_num_lin_solv_setups(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumLinSolvSetups(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_ida_get_num_err_test_fails(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumErrTestFails(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_ida_get_last_order(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    int v;

    flag = IDAGetLastOrder(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetLastOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value sunml_ida_get_current_order(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    int v;

    flag = IDAGetCurrentOrder(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetCurrentOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value sunml_ida_get_actual_init_step(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    realtype v;

    flag = IDAGetActualInitStep(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_ida_get_last_step(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1 (tmp);

    int flag;
    realtype v;

    flag = IDAGetLastStep(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetLastStep", flag);

    tmp = caml_copy_double(v);
    CAMLreturn(tmp);
}

CAMLprim value sunml_ida_get_current_step(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    realtype v;

    flag = IDAGetCurrentStep(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_ida_get_current_time(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    realtype v;

    flag = IDAGetCurrentTime(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetCurrentTime", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_ida_set_nonlin_conv_coef_ic(value vida_mem, value vcoef)
{
    CAMLparam2(vida_mem, vcoef);

    int flag = IDASetNonlinConvCoefIC(IDA_MEM_FROM_ML(vida_mem),
				      Double_val(vcoef));
    CHECK_FLAG("IDASetNonlinConvCoefIC", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_num_steps_ic(value vida_mem, value vmaxnh)
{
    CAMLparam2(vida_mem, vmaxnh);

    int flag = IDASetMaxNumStepsIC(IDA_MEM_FROM_ML(vida_mem), Int_val(vmaxnh));
    CHECK_FLAG("IDASetMaxNumStepsIC", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_num_jacs_ic(value vida_mem, value vmaxnj)
{
    CAMLparam2(vida_mem, vmaxnj);

    int flag = IDASetMaxNumJacsIC(IDA_MEM_FROM_ML(vida_mem), Int_val(vmaxnj));
    CHECK_FLAG("IDASetMaxNumJacsIC", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_num_iters_ic(value vida_mem, value vmaxnit)
{
    CAMLparam2(vida_mem, vmaxnit);

    int flag = IDASetMaxNumItersIC(IDA_MEM_FROM_ML(vida_mem), Int_val(vmaxnit));
    CHECK_FLAG("IDASetMaxNumItersIC", flag);

    CAMLreturn (Val_unit);
}

#if SUNDIALS_LIB_VERSION == 270
// work around a missing prototype in 2.7.0
int IDASetMaxBacksIC(void *, int);
#endif

CAMLprim value sunml_ida_set_max_backs_ic(value vida_mem, value vmaxbacks)
{
    CAMLparam2(vida_mem, vmaxbacks);

#if SUNDIALS_LIB_VERSION >= 270
    int flag = IDASetMaxBacksIC(IDA_MEM_FROM_ML(vida_mem), Int_val(vmaxbacks));
    CHECK_FLAG("IDASetMaxBacksIC", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_line_search_ic(value vida_mem, value vls)
{
    CAMLparam2(vida_mem, vls);

    int flag = IDASetLineSearchOffIC(IDA_MEM_FROM_ML(vida_mem), !Bool_val(vls));
    CHECK_FLAG("IDASetLineSearchOffIC", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_step_tolerance_ic(value vida_mem, value vsteptol)
{
    CAMLparam2(vida_mem, vsteptol);

    int flag = IDASetStepToleranceIC(IDA_MEM_FROM_ML(vida_mem),
				     Double_val(vsteptol));
    CHECK_FLAG("IDASetStepToleranceIC", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_get_num_backtrack_ops (value vida_mem)
{
    CAMLparam1 (vida_mem);
    int flag;
    long nbo;
    flag = IDAGetNumBacktrackOps (IDA_MEM_FROM_ML (vida_mem), &nbo);
    CHECK_FLAG ("IDAGetNumBcktrackOps", flag);
    CAMLreturn (Val_int (nbo));
}

CAMLprim value sunml_ida_set_max_ord(value vida_mem, value maxord)
{
    CAMLparam2(vida_mem, maxord);


    int flag = IDASetMaxOrd(IDA_MEM_FROM_ML(vida_mem), Int_val(maxord));
    CHECK_FLAG("IDASetMaxOrd", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_num_steps(value vida_mem, value mxsteps)
{
    CAMLparam2(vida_mem, mxsteps);


    int flag = IDASetMaxNumSteps(IDA_MEM_FROM_ML(vida_mem), Long_val(mxsteps));
    CHECK_FLAG("IDASetMaxNumSteps", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_init_step(value vida_mem, value hin)
{
    CAMLparam2(vida_mem, hin);


    int flag = IDASetInitStep(IDA_MEM_FROM_ML(vida_mem), Double_val(hin));
    CHECK_FLAG("IDASetInitStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_step(value vida_mem, value hmax)
{
    CAMLparam2(vida_mem, hmax);


    int flag = IDASetMaxStep(IDA_MEM_FROM_ML(vida_mem), Double_val(hmax));
    CHECK_FLAG("IDASetMaxStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_stop_time(value vida_mem, value tstop)
{
    CAMLparam2(vida_mem, tstop);


    int flag = IDASetStopTime(IDA_MEM_FROM_ML(vida_mem), Double_val(tstop));
    CHECK_FLAG("IDASetStopTime", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_err_test_fails(value vida_mem, value maxnef)
{
    CAMLparam2(vida_mem, maxnef);


    int flag = IDASetMaxErrTestFails(IDA_MEM_FROM_ML(vida_mem), Int_val(maxnef));
    CHECK_FLAG("IDASetMaxErrTestFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_nonlin_iters(value vida_mem, value maxcor)
{
    CAMLparam2(vida_mem, maxcor);


    int flag = IDASetMaxNonlinIters(IDA_MEM_FROM_ML(vida_mem), Int_val(maxcor));
    CHECK_FLAG("IDASetMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_max_conv_fails(value vida_mem, value maxncf)
{
    CAMLparam2(vida_mem, maxncf);


    int flag = IDASetMaxConvFails(IDA_MEM_FROM_ML(vida_mem), Int_val(maxncf));
    CHECK_FLAG("IDASetMaxConvFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_nonlin_conv_coef(value vida_mem, value nlscoef)
{
    CAMLparam2(vida_mem, nlscoef);


    int flag = IDASetNonlinConvCoef(IDA_MEM_FROM_ML(vida_mem), Double_val(nlscoef));
    CHECK_FLAG("IDASetNonlinConvCoef", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value sunml_ida_set_no_inactive_root_warn(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag = IDASetNoInactiveRootWarn(IDA_MEM_FROM_ML(vida_mem));
    CHECK_FLAG("IDASetNoInactiveRootWarn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_suppress_alg (value vida_mem, value vb)
{
    CAMLparam2(vida_mem, vb);

    int flag = IDASetSuppressAlg(IDA_MEM_FROM_ML(vida_mem),
				 Bool_val (vb));
    CHECK_FLAG("IDASetSuppressAlg", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_set_gs_type(value vida_mem, value vgstype)
{
    CAMLparam2(vida_mem, vgstype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = IDASpilsSetGSType(IDA_MEM_FROM_ML(vida_mem),
				 sunml_lsolver_gs_type(vgstype));
    CHECK_SPILS_FLAG("IDASpilsSetGSType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_set_max_restarts(value vida_mem, value vmaxr)
{
    CAMLparam2(vida_mem, vmaxr);
#if SUNDIALS_LIB_VERSION < 300
    int flag = IDASpilsSetMaxRestarts(IDA_MEM_FROM_ML(vida_mem),
				      Int_val(vmaxr));
    CHECK_SPILS_FLAG("IDASpilsSetMaxRestarts", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_eps_lin(value vida_mem, value eplifac)
{
    CAMLparam2(vida_mem, eplifac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDASetEpsLin(IDA_MEM_FROM_ML(vida_mem), Double_val(eplifac));
    CHECK_LS_FLAG("IDASetEpsLin", flag);
#else
    int flag = IDASpilsSetEpsLin(IDA_MEM_FROM_ML(vida_mem),
				 Double_val(eplifac));
    CHECK_SPILS_FLAG("IDASpilsSetEpsLin", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_set_increment_factor(value vida_mem, value dqincfac)
{
    CAMLparam2(vida_mem, dqincfac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDASetIncrementFactor(IDA_MEM_FROM_ML(vida_mem),
				     Double_val(dqincfac));
    CHECK_LS_FLAG("IDASetIncrementFactor", flag);
#else
    int flag = IDASpilsSetIncrementFactor(IDA_MEM_FROM_ML(vida_mem),
					  Double_val(dqincfac));
    CHECK_SPILS_FLAG("IDASpilsSetIncrementFactor", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_spils_set_maxl(value vida_mem, value maxl)
{
    CAMLparam2(vida_mem, maxl);
#if SUNDIALS_LIB_VERSION < 300
    int flag = IDASpilsSetMaxl(IDA_MEM_FROM_ML(vida_mem), Int_val(maxl));
    CHECK_SPILS_FLAG("IDASpilsSetMaxl", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}


/* statistic accessor functions */

CAMLprim value sunml_ida_get_tol_scale_factor(value vida_mem)
{
    CAMLparam1(vida_mem);

    realtype r;
    int flag = IDAGetTolScaleFactor(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetTolScaleFactor", flag);

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_ida_get_num_nonlin_solv_iters(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDAGetNumNonlinSolvIters(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_get_num_nonlin_solv_conv_fails(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDAGetNumNonlinSolvConvFails(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_get_nonlin_solv_stats(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(ret);

    long int nniters, nncfails;
    int flag = IDAGetNonlinSolvStats(IDA_MEM_FROM_ML(vida_mem),
				     &nniters, &nncfails);
    CHECK_FLAG("IDAGetNonlinSolvStats", flag);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, Val_long (nniters));
    Store_field (ret, 1, Val_long (nncfails));

    CAMLreturn(ret);
}

CAMLprim value sunml_ida_get_num_g_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDAGetNumGEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetNumGEvals", flag);

    CAMLreturn(Val_long(r));
}


CAMLprim value sunml_ida_dls_get_work_space(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetLinWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrwLS, &leniwLS);
    CHECK_LS_FLAG("IDAGetLinWorkSpace", flag);
#else
    int flag = IDADlsGetWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrwLS, &leniwLS);
    CHECK_DLS_FLAG("IDADlsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));

    CAMLreturn(r);
}


CAMLprim value sunml_ida_get_num_jac_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumJacEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumJacEvals", flag);
#else
    int flag = IDADlsGetNumJacEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_DLS_FLAG("IDADlsGetNumJacEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_get_num_lin_res_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumResEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumResEvals", flag);
#else
    int flag = IDADlsGetNumResEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_DLS_FLAG("IDADlsGetNumResEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

/* spils functions */

CAMLprim value sunml_ida_spils_get_num_lin_iters(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumLinIters(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumLinIters", flag);
#else
    int flag = IDASpilsGetNumLinIters(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_SPILS_FLAG("IDASpilsGetNumLinIters", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_spils_get_num_lin_conv_fails(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumLinConvFails(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumLinConvFails", flag);
#else
    int flag = IDASpilsGetNumConvFails(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_SPILS_FLAG("IDASpilsGetNumConvFails", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_spils_get_work_space(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = IDAGetLinWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrw, &leniw);
    CHECK_LS_FLAG("IDAGetLinWorkSpace", flag);
#else
    flag = IDASpilsGetWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrw, &leniw);
    CHECK_SPILS_FLAG("IDASpilsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_ida_spils_get_num_prec_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumPrecEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumPrecEvals", flag);
#else
    int flag = IDASpilsGetNumPrecEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_SPILS_FLAG("IDASpilsGetNumPrecEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_spils_get_num_prec_solves(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumPrecSolves(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumPrecSolves", flag);
#else
    int flag = IDASpilsGetNumPrecSolves(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_SPILS_FLAG("IDASpilsGetNumPrecSolves", flag);
#endif

    CAMLreturn(Val_long(r));
}

#if SUNDIALS_LIB_VERSION <= 301
SUNDIALS_EXPORT int IDASpilsGetNumJTSetupEvals(void *ida_mem,
					       long int *njtsetups);
#endif

CAMLprim value sunml_ida_spils_get_num_jtsetup_evals(value vida_mem)
{
    CAMLparam1(vida_mem);
    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumJTSetupEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumJTSetupEvals", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    int flag = IDASpilsGetNumJTSetupEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_SPILS_FLAG("IDASpilsGetNumJTSetupEvals", flag);
#else
    r = 0;
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_spils_get_num_jtimes_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumJtimesEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumJtimesEvals", flag);
#else
    int flag = IDASpilsGetNumJtimesEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_SPILS_FLAG("IDASpilsGetNumJtimesEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_ida_spils_get_num_lin_res_evals (value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDAGetNumLinResEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_LS_FLAG("IDAGetNumLinResEvals", flag);
#else
    int flag = IDASpilsGetNumResEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_SPILS_FLAG("IDASpilsGetNumResEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

