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

#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

// When we compile with sensitivity (CVODES), we are obliged to use the
// cvodes/cvodes_* header files. In fact, nearly everything functions
// correctly if the cvode/cvode_* header files are used instead (the
// function prototypes are identical), except for the function
// sunml_cvode_set_alternate which relies on the internal representation of
// CVodeMem (for cv_lsolve, etc.). In any case, it seems a better idea to
// use the appropriate header files even if this introduces a minor
// complication in the build system.

#ifdef SUNDIALSML_WITHSENS
/* CVODES (with sensitivity) */

#include <cvodes/cvodes.h>

/* linear solvers */
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_bandpre.h>

#if   400 <= SUNDIALS_LIB_VERSION
#include <cvodes/cvodes_ls.h>
#elif 300 <= SUNDIALS_LIB_VERSION
#include <cvodes/cvodes_direct.h>
#include <cvodes/cvodes_spils.h>
#else
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#endif

#if SUNDIALS_LIB_VERSION < 410
#include <cvodes/cvodes_impl.h>
#endif

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <cvodes/cvodes_lapack.h>
#endif

#else
/* CVODE (without sensitivity) */

#include <cvode/cvode.h>

/* linear solvers */
#include <cvode/cvode_diag.h>
#include <cvode/cvode_bandpre.h>

#if SUNDIALS_LIB_VERSION < 410
#include <cvode/cvode_impl.h>
#endif

#if   400 <= SUNDIALS_LIB_VERSION
#include <cvode/cvode_ls.h>
#elif 300 <= SUNDIALS_LIB_VERSION
#include <cvode/cvode_direct.h>
#include <cvode/cvode_spils.h>
#else
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#endif

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <cvode/cvode_lapack.h>
#endif

#endif

#include "../lsolvers/sundials_matrix_ml.h"
#include "../lsolvers/sundials_linearsolver_ml.h"
#include "../lsolvers/sundials_nonlinearsolver_ml.h"
#include "cvode_ml.h"
#include "../nvectors/nvector_ml.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256


CAMLprim value sunml_cvode_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (CVODE, exns);
    CAMLreturn (Val_unit);
}


/* callbacks */

static void errh(int error_code,
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
    value r = caml_callback_exn (Field(session, RECORD_CVODE_SESSION_ERRH), a);
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined error handler");

    CAMLreturn0;
}

CAMLprim value sunml_cvode_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = CVodeSetErrHandlerFn(CVODE_MEM_FROM_ML(vdata), errh,
				    CVODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSetErrHandlerFn(CVODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

#if defined SUNDIALS_BUILD_WITH_MONITORING && !SUNDIALSML_WITHSENS
int monitorfn(void *cvode_mem, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);

    WEAK_DEREF (session, *(value*)user_data);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(Field(session, RECORD_CVODE_SESSION_MONITORFN),
				session);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

CAMLprim value sunml_cvode_set_monitor_fn(value vdata, value vset)
{
    CAMLparam2(vdata, vset);
#if defined SUNDIALS_BUILD_WITH_MONITORING && !SUNDIALSML_WITHSENS
    int flag = CVodeSetMonitorFn(CVODE_MEM_FROM_ML(vdata),
				 Bool_val(vset) ? monitorfn : NULL);
    CHECK_FLAG("CVodeSetMonitorFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_monitor_frequency(value vdata, value vfreq)
{
    CAMLparam2(vdata, vfreq);
#if defined SUNDIALS_BUILD_WITH_MONITORING && !SUNDIALSML_WITHSENS
    int flag = CVodeSetMonitorFrequency(CVODE_MEM_FROM_ML(vdata),
					Int_val(vfreq));
    CHECK_FLAG("CVodeSetMonitorFrequency", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

int sunml_cvode_translate_exception (value session, value exn,
			       recoverability recoverable)
{
    CAMLparam2(session, exn);
    CAMLlocal1(bucket);

    if (recoverable && Field(exn, 0) == SUNDIALS_EXN_TAG (RecoverableFailure))
	CAMLreturnT (int, 1);

    /* Unrecoverable error.  Save the exception and return -1.  */
    bucket = caml_alloc_small (1,0);
    Field (bucket, 0) = exn;
    Store_field (session, RECORD_CVODE_SESSION_EXN_TEMP, bucket);
    CAMLreturnT (int, -1);
}

static int rhsfn(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);

    WEAK_DEREF (session, *(value*)user_data);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(Field(session, RECORD_CVODE_SESSION_RHSFN),
				 caml_copy_double(t),
				 NVEC_BACKLINK(y),
				 NVEC_BACKLINK(ydot));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if 580 <= SUNDIALS_LIB_VERSION
static int nlsrhsfn(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_CVODE_SESSION_NLS_RHSFN),
				 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int roots(sunrealtype t, N_Vector y, sunrealtype *gout, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    value *backref = user_data;
    intnat nroots;

    WEAK_DEREF (session, *backref);

    nroots = CVODE_NROOTS_FROM_ML (session);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);
    args[2] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(session, RECORD_CVODE_SESSION_ROOTSFN),
				  3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, UNRECOVERABLE));
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocal1 (session);
    value *backref = user_data;

    WEAK_DEREF (session, *backref);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (Field(session, RECORD_CVODE_SESSION_ERRW),
				  NVEC_BACKLINK (y), NVEC_BACKLINK (ewt));
    if (Is_exception_result (r)) {
	r = Extract_exception (r);
	if (Field (r, 0) != SUNDIALS_EXN_TAG (NonPositiveEwt))
	    sunml_warn_discarded_exn (r, "user-defined error weight fun");
	CAMLreturnT (int, -1);
    }

    CAMLreturnT (int, 0);
}

value sunml_cvode_make_jac_arg(sunrealtype t, N_Vector y, N_Vector fy, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_CVODE_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_FY, NVEC_BACKLINK(fy));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

value sunml_cvode_make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(3);
    Store_field(r, 0, NVEC_BACKLINK(tmp1));
    Store_field(r, 1, NVEC_BACKLINK(tmp2));
    Store_field(r, 2, NVEC_BACKLINK(tmp3));
    CAMLreturn(r);
}

#if 300 <= SUNDIALS_LIB_VERSION
static int jacfn(
	sunrealtype t,
	N_Vector y,
	N_Vector fy,
	SUNMatrix Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 2);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_cvode_make_jac_arg (t, y, fy,
				  sunml_cvode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = MAT_BACKLINK(Jac);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#else

/* Dense and band Jacobians only work with serial NVectors.  */
static int jacfn(
	long int n,
	sunrealtype t,
	N_Vector y,
	N_Vector fy,	     
	DlsMat Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(Jac));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_cvode_make_jac_arg (t, y, fy,
				  sunml_cvode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bandjacfn(
	long int N,
	long int mupper,
	long int mlower,
	sunrealtype t,
	N_Vector y,
	N_Vector fy,
	DlsMat Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(Jac));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_cvode_make_jac_arg(t, y, fy,
				 sunml_cvode_make_triple_tmp(tmp1, tmp2, tmp3));
    args[1] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static int linsysfn(
	sunrealtype t,
	N_Vector y,
	N_Vector fy,
	SUNMatrix M,
	booleantype jok,
	booleantype *jcur,
	sunrealtype gamma,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 4);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1);

    args[0] = sunml_cvode_make_jac_arg (t, y, fy,
				  sunml_cvode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = MAT_BACKLINK(M);
    args[2] = Val_bool(jok);
    args[3] = caml_copy_double(gamma);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    /* Update jcur; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcur = Bool_val (r);
	CAMLreturnT(int, 0);
    }
    r = Extract_exception (r);

    CAMLreturnT(int, sunml_cvode_translate_exception (session, r, RECOVERABLE));
}
#endif

static int precsetupfn(sunrealtype t,
		       N_Vector y,
		       N_Vector fy,
		       booleantype jok,
		       booleantype *jcurPtr,
		       sunrealtype gamma,
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
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = sunml_cvode_make_jac_arg(t, y, fy, Val_unit);
    args[1] = Val_bool(jok);
    args[2] = caml_copy_double(gamma);

    cb = CVODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODE_SPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    /* Update jcurPtr; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcurPtr = Bool_val (r);
	CAMLreturnT(int, 0);
    }
    r = Extract_exception (r);

    CAMLreturnT(int, sunml_cvode_translate_exception (session, r, RECOVERABLE));
}

static value make_spils_solve_arg(
	N_Vector r,
	sunrealtype gamma,
	sunrealtype delta,
	int lr)

{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(RECORD_CVODE_SPILS_SOLVE_ARG_SIZE);
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_RHS, NVEC_BACKLINK(r));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_GAMMA,
                caml_copy_double(gamma));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_DELTA,
                caml_copy_double(delta));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_LEFT, Val_bool (lr == 1));

    CAMLreturn(v);
}

static int precsolvefn(
	sunrealtype t,
	N_Vector y,
	N_Vector fy,
	N_Vector rvec,
	N_Vector z,
	sunrealtype gamma,
	sunrealtype delta,
	int lr,
	void *user_data
#if SUNDIALS_LIB_VERSION < 300
	,
	N_Vector tmp
#endif
	)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = sunml_cvode_make_jac_arg(t, y, fy, Val_unit);
    args[1] = make_spils_solve_arg(rvec, gamma, delta, lr);
    args[2] = NVEC_BACKLINK(z);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_PRECFNS_FROM_ML(session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODE_SPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int jactimesfn(N_Vector v,
		      N_Vector Jv,
		      sunrealtype t,
		      N_Vector y,
		      N_Vector fy,
		      void *user_data,
		      N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = sunml_cvode_make_jac_arg(t, y, fy, NVEC_BACKLINK(tmp));
    args[1] = NVEC_BACKLINK(v);
    args[2] = NVEC_BACKLINK(Jv);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if 300 <= SUNDIALS_LIB_VERSION
static int jacsetupfn(sunrealtype t,
		      N_Vector y,
		      N_Vector fy,
		      void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = sunml_cvode_make_jac_arg(t, y, fy, Val_unit);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, arg);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

#if 530 <= SUNDIALS_LIB_VERSION
static int jactimesrhsfn(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
int projfn(sunrealtype t, N_Vector ycur, N_Vector corr,
           sunrealtype epsProj, N_Vector err, void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, verr);
    CAMLlocalN(args, 5);

    WEAK_DEREF (session, *(value*)user_data);

    if (err == NULL) {
	verr = Val_none;
    } else {
	Store_some(verr, NVEC_BACKLINK(err));
    }

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(ycur);
    args[2] = NVEC_BACKLINK(corr);
    args[3] = caml_copy_double(epsProj);
    args[4] = verr;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_CVODE_SESSION_PROJFN),
				 5, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_cvode_dls_dense (value vcvode_mem, value vneqs, value vset_jac)
{
    CAMLparam3(vcvode_mem, vneqs, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    long neqs = Long_val(vneqs);
    int flag;

    flag = CVDense (cvode_mem, neqs);
    CHECK_FLAG ("CVDense", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetDenseJacFn(cvode_mem, jacfn);
	CHECK_DLS_FLAG("CVDlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_dls_lapack_dense (value vcvode_mem, value vneqs,
					 value vset_jac)
{
    CAMLparam3 (vcvode_mem, vneqs, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = CVLapackDense (cvode_mem, neqs);
    CHECK_FLAG ("CVLapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetDenseJacFn (cvode_mem, jacfn);
	CHECK_DLS_FLAG("CVDlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_dls_band (value vcvode_mem, value vneqs,
				 value vmupper, value vmlower,
				 value vset_jac)
{
    CAMLparam5(vcvode_mem, vneqs, vmupper, vmlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = CVBand (cvode_mem, neqs, Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("CVBand", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetBandJacFn(cvode_mem, bandjacfn);
	CHECK_DLS_FLAG("CVDlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_dls_lapack_band (value vcvode_mem, value vneqs,
					value vmupper, value vmlower,
					value vset_jac)
{
    CAMLparam5(vcvode_mem, vneqs, vmupper, vmlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    long neqs = Long_val(vneqs);
    int flag;

    flag = CVLapackBand (cvode_mem, neqs,
			 Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("CVLapackBand", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetBandJacFn(cvode_mem, bandjacfn);
	CHECK_DLS_FLAG("CVDlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_linear_solver (value vcvode_mem, value vlsolv,
					      value vojmat, value vhasjac,
					      value vhaslsf)
{
    CAMLparam5(vcvode_mem, vlsolv, vojmat, vhasjac, vhaslsf);
#if 400 <= SUNDIALS_LIB_VERSION
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = CVodeSetLinearSolver(cvode_mem, lsolv, jmat);
    CHECK_LS_FLAG ("CVodeSetLinearSolver", flag);
    flag = CVodeSetJacFn(cvode_mem, Bool_val(vhasjac) ? jacfn : NULL);
    CHECK_LS_FLAG("CVodeSetJacFn", flag);

#if 500 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vhaslsf)) {
	CVodeSetLinSysFn(cvode_mem, linsysfn);
	CHECK_LS_FLAG("CVodeSetLinSysFn", flag);
    }
#endif

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_dls_set_linear_solver (value vcvode_mem, value vlsolv,
					      value vjmat, value vhasjac)
{
    CAMLparam4(vcvode_mem, vlsolv, vjmat, vhasjac);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = CVDlsSetLinearSolver(cvode_mem, lsolv, jmat);
    CHECK_DLS_FLAG ("CVDlsSetLinearSolver", flag);
    flag = CVDlsSetJacFn(cvode_mem, Bool_val(vhasjac) ? jacfn : NULL);
    CHECK_DLS_FLAG("CVDlsSetJacFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_linear_solver (value vcvode_mem, value vlsolv)
{
    CAMLparam2(vcvode_mem, vlsolv);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    int flag;

    flag = CVSpilsSetLinearSolver(cvode_mem, lsolv);
    CHECK_SPILS_FLAG ("CVSpilsSetLinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_preconditioner (value vsession,
					       value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = CVODE_MEM_FROM_ML (vsession);
#if 400 <= SUNDIALS_LIB_VERSION
    CVLsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = CVodeSetPreconditioner (mem, setup, precsolvefn);
    CHECK_LS_FLAG ("CVodeSetPreconditioner", flag);
#else
    CVSpilsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = CVSpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_SPILS_FLAG ("CVSpilsSetPreconditioner", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_banded_preconditioner (value vsession,
						      value vneqs,
						      value vmupper,
						      value vmlower)
{
    CAMLparam3 (vsession, vmupper, vmlower);
    long neqs = Index_val (vneqs);
    int flag = CVBandPrecInit (CVODE_MEM_FROM_ML (vsession), neqs,
			       Index_val (vmupper), Index_val (vmlower));
    CHECK_FLAG ("CVBandPrecInit", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_jac_times(value vdata, value vhas_setup,
					 value vhas_times)
{
    CAMLparam3(vdata, vhas_setup, vhas_times);
#if 400 <= SUNDIALS_LIB_VERSION
    CVLsJacTimesSetupFn setup = Bool_val (vhas_setup) ? jacsetupfn : NULL;
    CVLsJacTimesVecFn   times = Bool_val (vhas_times) ? jactimesfn : NULL;

    int flag = CVodeSetJacTimes(CVODE_MEM_FROM_ML(vdata), setup, times);
    CHECK_LS_FLAG("CVodeSetJacTimes", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    CVSpilsJacTimesSetupFn setup = Bool_val (vhas_setup) ? jacsetupfn : NULL;
    CVSpilsJacTimesVecFn   times = Bool_val (vhas_times) ? jactimesfn : NULL;

    int flag = CVSpilsSetJacTimes(CVODE_MEM_FROM_ML(vdata), setup, times);
    CHECK_SPILS_FLAG("CVSpilsSetJacTimes", flag);
#else
    CVSpilsJacTimesVecFn times = Bool_val (vhas_times) ? jactimesfn : NULL;
    int flag = CVSpilsSetJacTimesVecFn(CVODE_MEM_FROM_ML(vdata), times);
    CHECK_SPILS_FLAG("CVSpilsSetJacTimesVecFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_jac_times_rhsfn(value vdata, value vhas_rhsfn)
{
    CAMLparam2(vdata, vhas_rhsfn);
#if 530 <= SUNDIALS_LIB_VERSION
    CVRhsFn rhsfn = Bool_val (vhas_rhsfn) ? jactimesrhsfn : NULL;

    int flag = CVodeSetJacTimesRhsFn(CVODE_MEM_FROM_ML(vdata), rhsfn);
    CHECK_LS_FLAG("CVodeSetJacTimesRhsFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

#if 400 <= SUNDIALS_LIB_VERSION

// hack to work around lack of CVodeGetUserData
typedef struct {
  sunrealtype cv_uround;
  CVRhsFn cv_f;
  void *cv_user_data;
  //...
} *StartOf_CVodeMem;

static value sunml_cvode_session_to_value(void *cvode_mem)
{
    value session;
    // void *user_data = CVodeGetUserData(cvode_mem);
    void *user_data = ((StartOf_CVodeMem)cvode_mem)->cv_user_data;

    WEAK_DEREF (session, *(value*)user_data);
    return session;
}

static void* sunml_cvode_session_from_value(value vcvode_mem)
{
    return (CVODE_MEM_FROM_ML(vcvode_mem));
}
#endif

CAMLprim value sunml_cvode_set_nonlinear_solver(value vcvode_mem, value vnlsolv)
{
    CAMLparam2(vcvode_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_cvode_session_to_value,
				   sunml_cvode_session_from_value);

    int flag = CVodeSetNonlinearSolver(cvode_mem, nlsolv);
    CHECK_FLAG ("CVodeSetNonlinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_wf_tolerances (value vdata)
{
    CAMLparam1(vdata);
 
    int flag = CVodeWFtolerances(CVODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("CVodeWFtolerances", flag);

    CAMLreturn (Val_unit);
}

/* basic interface */

/* CVodeCreate() + CVodeInit().  */
CAMLprim value sunml_cvode_init(value weakref, value lmm, value iter, value initial,
			        value t0)
{
    CAMLparam5(weakref, lmm, iter, initial, t0);
    CAMLlocal2(r, vcvode_mem);

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

#if 400 <= SUNDIALS_LIB_VERSION
    void *cvode_mem = CVodeCreate(lmm_c);
#else
    void *cvode_mem = CVodeCreate(lmm_c,
				  Bool_val(iter) ? CV_NEWTON : CV_FUNCTIONAL);
#endif
    if (cvode_mem == NULL)
	caml_failwith("CVodeCreate returned NULL");

    vcvode_mem = sunml_wrap_session_pointer(cvode_mem);

    N_Vector initial_nv = NVEC_VAL(initial);
    flag = CVodeInit(cvode_mem, rhsfn, Double_val(t0), initial_nv);
    if (flag != CV_SUCCESS) {
	CVodeFree (&cvode_mem);
	CHECK_FLAG("CVodeInit", flag);
    }

    value *backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	CVodeFree (&cvode_mem);
	caml_raise_out_of_memory();
    }
    CVodeSetUserData (cvode_mem, backref);

    r = caml_alloc_tuple (2);
    Store_field (r, 0, vcvode_mem);
    Store_field (r, 1, (value)backref);

    CAMLreturn(r);
}

/* Set the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim value sunml_cvode_root_init (value vdata, value vnroots)
{
    CAMLparam2 (vdata, vnroots);
    void *cvode_mem = CVODE_MEM_FROM_ML (vdata);
    int nroots = Int_val (vnroots);
    int flag = CVodeRootInit (cvode_mem, nroots, roots);
    CHECK_FLAG ("CVodeRootInit", flag);
    Store_field (vdata, RECORD_CVODE_SESSION_NROOTS, vnroots);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_sv_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = CVodeSVtolerances(CVODE_MEM_FROM_ML(vdata),
				 Double_val(reltol), atol_nv);
    CHECK_FLAG("CVodeSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_reinit(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    N_Vector y0_nv = NVEC_VAL(y0);
    int flag = CVodeReInit(CVODE_MEM_FROM_ML(vdata), Double_val(t0), y0_nv);
    CHECK_FLAG("CVodeReInit", flag);

    CAMLreturn (Val_unit);
}

static value solver(value vdata, value nextt, value vy, int onestep)
{
    CAMLparam3(vdata, nextt, vy);
    CAMLlocal1(ret);
    sunrealtype tret;
    int flag;
    N_Vector y;
    enum cvode_solver_result_tag result = -1;

    y = NVEC_VAL (vy);
    // Caml_ba_data_val(y) must not be shifted by the OCaml GC during this
    // function call, which calls Caml through the callback f.  Is this
    // guaranteed?
    flag = CVode (CVODE_MEM_FROM_ML (vdata), Double_val (nextt), y, &tret,
		  onestep ? CV_ONE_STEP : CV_NORMAL);

    switch (flag) {
    case CV_SUCCESS:
	result = VARIANT_CVODE_SOLVER_RESULT_SUCCESS;
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
	sunml_cvode_check_flag("CVode", flag, CVODE_MEM_FROM_ML(vdata));
    }

    assert (Field (vdata, RECORD_CVODE_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}

CAMLprim value sunml_cvode_solve_normal(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value sunml_cvode_solve_one_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 1));
}

CAMLprim value sunml_cvode_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVEC_VAL(vy);

    int flag = CVodeGetDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
			   Int_val(vk), y_nv);
    CHECK_FLAG("CVodeGetDky", flag);
    
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_get_err_weights(value vcvode_mem, value verrws)
{
    CAMLparam2(vcvode_mem, verrws);

    N_Vector errws_nv = NVEC_VAL(verrws);

    int flag = CVodeGetErrWeights(CVODE_MEM_FROM_ML(vcvode_mem), errws_nv);
    CHECK_FLAG("CVodeGetErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_get_est_local_errors(value vcvode_mem, value vele)
{
    CAMLparam2(vcvode_mem, vele);

    N_Vector ele_nv = NVEC_VAL(vele);

    int flag = CVodeGetEstLocalErrors(CVODE_MEM_FROM_ML(vcvode_mem), ele_nv);
    CHECK_FLAG("CVodeGetEstLocalErrors", flag);

    CAMLreturn (Val_unit);
}

value sunml_cvode_last_lin_exception(void *cvode_mem)
{
#if 400 <= SUNDIALS_LIB_VERSION
    long int lsflag;

    if (cvode_mem != NULL
	  && CVodeGetLastLinFlag(cvode_mem, &lsflag) == CVLS_SUCCESS)
	return sunml_lsolver_exception_from_flag(lsflag);
#endif
    return Val_none;
}

void sunml_cvode_check_flag(const char *call, int flag, void *cvode_mem)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	    || flag == CV_ROOT_RETURN
	    || flag == CV_TSTOP_RETURN) return;

    switch (flag) {
	case CV_ILL_INPUT:
	    caml_raise_constant(CVODE_EXN(IllInput));

	case CV_TOO_CLOSE:
	    caml_raise_constant(CVODE_EXN(TooClose));

	case CV_TOO_MUCH_WORK:
	    caml_raise_constant(CVODE_EXN(TooMuchWork));

	case CV_TOO_MUCH_ACC:
	    caml_raise_constant(CVODE_EXN(TooMuchAccuracy));

	case CV_ERR_FAILURE:
	    caml_raise_constant(CVODE_EXN(ErrFailure));

	case CV_CONV_FAILURE:
	    caml_raise_constant(CVODE_EXN(ConvergenceFailure));

	case CV_LINIT_FAIL:
	    caml_raise_constant(CVODE_EXN(LinearInitFailure));

	case CV_LSETUP_FAIL:
	    caml_raise_with_arg(CVODE_EXN(LinearSetupFailure),
				sunml_cvode_last_lin_exception(cvode_mem));

	case CV_LSOLVE_FAIL:
	    caml_raise_with_arg(CVODE_EXN(LinearSolveFailure),
				sunml_cvode_last_lin_exception(cvode_mem));

	case CV_RHSFUNC_FAIL:
	    caml_raise_constant(CVODE_EXN(RhsFuncFailure));

	case CV_FIRST_RHSFUNC_ERR:
	    caml_raise_constant(CVODE_EXN(FirstRhsFuncFailure));

	case CV_REPTD_RHSFUNC_ERR:
	    caml_raise_constant(CVODE_EXN(RepeatedRhsFuncFailure));

	case CV_UNREC_RHSFUNC_ERR:
	    caml_raise_constant(CVODE_EXN(UnrecoverableRhsFuncFailure));

	case CV_RTFUNC_FAIL:
	    caml_raise_constant(CVODE_EXN(RootFuncFailure));

#if 320 <= SUNDIALS_LIB_VERSION
	case CV_CONSTR_FAIL:
	    caml_raise_constant(CVODE_EXN(ConstraintFailure));
#endif

	case CV_BAD_K:
	    caml_raise_constant(CVODE_EXN(BadK));

	case CV_BAD_T:
	    caml_raise_constant(CVODE_EXN(BadT));

#if 400 <= SUNDIALS_LIB_VERSION
	case CV_NLS_INIT_FAIL:
	    caml_raise_constant(CVODE_EXN(NonlinearInitFailure));

	case CV_NLS_SETUP_FAIL:
	    caml_raise_constant(CVODE_EXN(NonlinearSetupFailure));

	case CV_VECTOROP_ERR:
	    caml_raise_constant(CVODE_EXN(VectorOpErr));
#endif
#if 500 <= SUNDIALS_LIB_VERSION
	case CV_NLS_FAIL:
	    caml_raise_constant(CVODE_EXN(NonlinearFailure));
#endif
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
	case CV_PROJFUNC_FAIL:
	    caml_raise_constant(CVODE_EXN(ProjFuncFailure));

	case CV_REPTD_PROJFUNC_ERR:
	    caml_raise_constant(CVODE_EXN(RepeatedProjFuncError));

	case CV_PROJ_MEM_NULL:
	    caml_raise_constant(CVODE_EXN(ProjectionNotEnabled));
#endif

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

#if 400 <= SUNDIALS_LIB_VERSION
void sunml_cvode_check_ls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CVLS_SUCCESS) return;

    switch (flag) {
	case CVLS_ILL_INPUT:
	    caml_raise_constant(CVODE_EXN(IllInput));

	case CVLS_MEM_FAIL:
	    caml_raise_out_of_memory();

	case CVLS_SUNMAT_FAIL:
	case CVLS_SUNLS_FAIL:
	case CVLS_JACFUNC_UNRECVR:
	case CVLS_JACFUNC_RECVR:
	default:
	    /* e.g. CVLS_MEM_NULL, CVLS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     CVodeGetLinReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#else
void sunml_cvode_check_dls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CVDLS_SUCCESS) return;

    switch (flag) {
	case CVDLS_ILL_INPUT:
	    caml_raise_constant(CVODE_EXN(IllInput));

	case CVDLS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if 300 <= SUNDIALS_LIB_VERSION
	case CVDLS_SUNMAT_FAIL:
#endif
	case CVDLS_JACFUNC_UNRECVR:
	case CVDLS_JACFUNC_RECVR:
	default:
	    /* e.g. CVDLS_MEM_NULL, CVDLS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVDlsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

void sunml_cvode_check_spils_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CVSPILS_SUCCESS) return;

    switch (flag) {
	case CVSPILS_ILL_INPUT:
	    caml_raise_constant(CVODE_EXN(IllInput));

	case CVSPILS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if 300 <= SUNDIALS_LIB_VERSION
	case CVSPILS_SUNLS_FAIL:
#endif
	default:
	    /* e.g. CVSPILS_MEM_NULL, CVSPILS_PMEM_NULL, CVSPILS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     CVSpilsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#endif

/* basic interface */

CAMLprim value sunml_cvode_session_finalize(value vdata)
{
    if (CVODE_MEM_FROM_ML(vdata) != NULL) {
	void *cvode_mem = CVODE_MEM_FROM_ML(vdata);
	value *backref = CVODE_BACKREF_FROM_ML(vdata);
	CVodeFree(&cvode_mem);
	sunml_sundials_free_value(backref);
    }

    return Val_unit;
}

CAMLprim value sunml_cvode_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("CVodeSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < CVODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = CVodeGetRootInfo(CVODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("CVodeGetRootInfo", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_get_integrator_stats(value vdata)
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

    sunrealtype hinused;
    sunrealtype hlast;
    sunrealtype hcur;
    sunrealtype tcur;

    flag = CVodeGetIntegratorStats(CVODE_MEM_FROM_ML(vdata),
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

    r = caml_alloc_tuple(RECORD_CVODE_INTEGRATOR_STATS_SIZE);
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_RHS_EVALS, Val_long(nfevals));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS, Val_long(nlinsetups));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_ERROR_TEST_FAILURES, Val_long(netfails));

    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_LAST_INTERNAL_ORDER, Val_int(qlast));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER, Val_int(qcur));

    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_INITIAL_STEP_SIZE, caml_copy_double(hinused));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_LAST_STEP_SIZE, caml_copy_double(hlast));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_NEXT_STEP_SIZE, caml_copy_double(hcur));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_INTERNAL_TIME, caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_get_linear_solver_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS

    int flag;

    long int njevals;
    long int nfevalsLS;
    long int nliters;
    long int nlcfails;
    long int npevals;
    long int npsolves;
    long int njtsetup;
    long int njtimes;

    flag = CVodeGetLinSolveStats(CVODE_MEM_FROM_ML(vdata),
				 &njevals,
				 &nfevalsLS,
				 &nliters,
				 &nlcfails,
				 &npevals,
				 &npsolves,
				 &njtsetup,
				 &njtimes); 
    CHECK_FLAG("CVodeGetLinSolveStats", flag);

    r = caml_alloc_tuple(RECORD_CVODE_LINEAR_SOLVER_STATS_SIZE);
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_JAC_EVALS,
		   Val_long(njevals));
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_LIN_RHS_EVALS,
		   Val_long(nfevalsLS));
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_LIN_ITERS,
		   Val_long(nliters));
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_LIN_CONV_FAILS,
		   Val_long(nlcfails));
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_PREC_EVALS,
		   Val_long(npevals));
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_PREC_SOLVES,
		   Val_long(npsolves));
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_JTSETUP_EVALS,
		   Val_long(njtsetup));
    Store_field(r, RECORD_CVODE_LINEAR_SOLVER_STATS_JTIMES_EVALS,
		   Val_long(njtimes));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_cvode_set_error_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

    int flag = CVodeSetErrFile(CVODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("CVodeSetErrFile", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_functional (value vdata)
{
    CAMLparam1 (vdata);
#if SUNDIALS_LIB_VERSION < 400
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vdata), CV_FUNCTIONAL);
    CHECK_FLAG ("CVodeSetIterType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_newton (value vdata)
{
    CAMLparam1 (vdata);
#if SUNDIALS_LIB_VERSION < 400
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vdata), CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < CVODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = CVodeSetRootDirection(CVODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("CVodeSetRootDirection", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_prec_type(value vcvode_mem, value vptype)
{
    CAMLparam2(vcvode_mem, vptype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = CVSpilsSetPrecType(CVODE_MEM_FROM_ML(vcvode_mem),
				  sunml_lsolver_precond_type(vptype));
    CHECK_SPILS_FLAG("CVSpilsSetPrecType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_diag (value vcvode_mem)
{
    CAMLparam1 (vcvode_mem);
    int flag;

    flag = CVDiag(CVODE_MEM_FROM_ML (vcvode_mem));
    CHECK_FLAG("CVDiag", flag);
    CAMLreturn (Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

CAMLprim value sunml_cvode_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = CVodeGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_FLAG("CVodeGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_get_num_steps(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumSteps(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvode_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvode_get_num_lin_solv_setups(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumLinSolvSetups(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvode_get_num_err_test_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumErrTestFails(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvode_get_last_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    int v;

    flag = CVodeGetLastOrder(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetLastOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value sunml_cvode_get_current_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    int v;

    flag = CVodeGetCurrentOrder(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value sunml_cvode_get_current_state(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(vnv);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector nv;
    int flag = CVodeGetCurrentState(CVODE_MEM_FROM_ML(vcvode_mem), &nv);
    CHECK_FLAG("CVodeGetCurrentState", flag);

    vnv = NVEC_BACKLINK(nv);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim value sunml_cvode_get_nonlin_system_data(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(vnv);
#if 540 <= SUNDIALS_LIB_VERSION
    sunrealtype tn, gamma, rl1;
    N_Vector ypred, yn, fn, zn1;
    void *user_data;

    int flag = CVodeGetNonlinearSystemData(CVODE_MEM_FROM_ML(vcvode_mem),
		    &tn, &ypred, &yn, &fn, &gamma, &rl1, &zn1, &user_data);
    CHECK_FLAG("CVodeGetNonlinearSystemData", flag);

    vnv = caml_alloc_tuple(RECORD_CVODE_NONLIN_SYSTEM_DATA_SIZE);
    Store_field(vnv, RECORD_CVODE_NONLIN_SYSTEM_DATA_TN,
	    caml_copy_double(tn));
    Store_field(vnv, RECORD_CVODE_NONLIN_SYSTEM_DATA_YPRED,
	    NVEC_BACKLINK(ypred));
    Store_field(vnv, RECORD_CVODE_NONLIN_SYSTEM_DATA_YN,
	    NVEC_BACKLINK(yn));
    Store_field(vnv, RECORD_CVODE_NONLIN_SYSTEM_DATA_FN,
	    NVEC_BACKLINK(fn));
    Store_field(vnv, RECORD_CVODE_NONLIN_SYSTEM_DATA_GAMMA,
	    caml_copy_double(gamma));
    Store_field(vnv, RECORD_CVODE_NONLIN_SYSTEM_DATA_RL1,
	    caml_copy_double(rl1));
    Store_field(vnv, RECORD_CVODE_NONLIN_SYSTEM_DATA_ZN1,
	    NVEC_BACKLINK(zn1));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim void sunml_cvode_compute_state(value vcvode_mem,
					 value vycor, value vyz)
{
    CAMLparam3(vcvode_mem, vycor, vyz);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = CVodeComputeState(CVODE_MEM_FROM_ML(vcvode_mem),
				 NVEC_VAL(vycor), NVEC_VAL(vyz));
    CHECK_FLAG("CVodeComputeState", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim value sunml_cvode_get_current_gamma(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    double gamma;

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetCurrentGamma(CVODE_MEM_FROM_ML(vcvode_mem), &gamma);
    CHECK_FLAG("CVodeGetCurrentGamma", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(gamma));
}

CAMLprim value sunml_cvode_get_current_gamma_unboxed(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    double gamma;

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetCurrentGamma(CVODE_MEM_FROM_ML(vcvode_mem), &gamma);
    CHECK_FLAG("CVodeGetCurrentGamma", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturnT(double, gamma);
}

CAMLprim value sunml_cvode_get_actual_init_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    sunrealtype v;

    flag = CVodeGetActualInitStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_cvode_get_last_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    sunrealtype v;

    flag = CVodeGetLastStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetLastStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_cvode_get_current_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    sunrealtype v;

    flag = CVodeGetCurrentStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_cvode_get_current_time(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    sunrealtype v;

    flag = CVodeGetCurrentTime(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentTime", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_cvode_set_max_ord(value vcvode_mem, value maxord)
{
    CAMLparam2(vcvode_mem, maxord);

    int flag = CVodeSetMaxOrd(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxord));
    CHECK_FLAG("CVodeSetMaxOrd", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_max_num_steps(value vcvode_mem, value mxsteps)
{
    CAMLparam2(vcvode_mem, mxsteps);


    int flag = CVodeSetMaxNumSteps(CVODE_MEM_FROM_ML(vcvode_mem), Long_val(mxsteps));
    CHECK_FLAG("CVodeSetMaxNumSteps", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_max_hnil_warns(value vcvode_mem, value mxhnil)
{
    CAMLparam2(vcvode_mem, mxhnil);


    int flag = CVodeSetMaxHnilWarns(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(mxhnil));
    CHECK_FLAG("CVodeSetMaxHnilWarns", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_stab_lim_det(value vcvode_mem, value stldet)
{
    CAMLparam2(vcvode_mem, stldet);


    int flag = CVodeSetStabLimDet(CVODE_MEM_FROM_ML(vcvode_mem), Bool_val(stldet));
    CHECK_FLAG("CVodeSetStabLimDet", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_init_step(value vcvode_mem, value hin)
{
    CAMLparam2(vcvode_mem, hin);


    int flag = CVodeSetInitStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hin));
    CHECK_FLAG("CVodeSetInitStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_min_step(value vcvode_mem, value hmin)
{
    CAMLparam2(vcvode_mem, hmin);


    int flag = CVodeSetMinStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hmin));
    CHECK_FLAG("CVodeSetMinStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_max_step(value vcvode_mem, value hmax)
{
    CAMLparam2(vcvode_mem, hmax);


    int flag = CVodeSetMaxStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hmax));
    CHECK_FLAG("CVodeSetMaxStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_stop_time(value vcvode_mem, value tstop)
{
    CAMLparam2(vcvode_mem, tstop);


    int flag = CVodeSetStopTime(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(tstop));
    CHECK_FLAG("CVodeSetStopTime", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_max_err_test_fails(value vcvode_mem, value maxnef)
{
    CAMLparam2(vcvode_mem, maxnef);


    int flag = CVodeSetMaxErrTestFails(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxnef));
    CHECK_FLAG("CVodeSetMaxErrTestFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_max_nonlin_iters(value vcvode_mem, value maxcor)
{
    CAMLparam2(vcvode_mem, maxcor);


    int flag = CVodeSetMaxNonlinIters(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxcor));
    CHECK_FLAG("CVodeSetMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_max_conv_fails(value vcvode_mem, value maxncf)
{
    CAMLparam2(vcvode_mem, maxncf);


    int flag = CVodeSetMaxConvFails(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxncf));
    CHECK_FLAG("CVodeSetMaxConvFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_nonlin_conv_coef(value vcvode_mem, value nlscoef)
{
    CAMLparam2(vcvode_mem, nlscoef);


    int flag = CVodeSetNonlinConvCoef(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(nlscoef));
    CHECK_FLAG("CVodeSetNonlinConvCoef", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_constraints (value vcvode_mem, value vconstraints)
{
    CAMLparam2(vcvode_mem, vconstraints);
#if 320 <= SUNDIALS_LIB_VERSION
    int flag;

    N_Vector constraints = NVEC_VAL (vconstraints);
    flag = CVodeSetConstraints (CVODE_MEM_FROM_ML (vcvode_mem), constraints);
    CHECK_FLAG ("CVodeSetConstraints", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_clear_constraints (value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
#if 320 <= SUNDIALS_LIB_VERSION
    int flag;

    flag = CVodeSetConstraints (CVODE_MEM_FROM_ML (vcvode_mem), NULL);
    CHECK_FLAG ("CVodeSetConstraints", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_proj_fn(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeSetProjFn(CVODE_MEM_FROM_ML (vcvode_mem), projfn);
    CHECK_FLAG ("CVodeSetProjFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_nls_rhs_fn(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
#if 580 <= SUNDIALS_LIB_VERSION
    int flag;

    flag = CVodeSetNlsRhsFn(CVODE_MEM_FROM_ML (vcvode_mem), nlsrhsfn);
    CHECK_FLAG ("CVodeSetNlsRhsFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_proj_err_est(value vcvode_mem, value vonoff)
{
    CAMLparam2(vcvode_mem, vonoff);
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeSetProjErrEst(CVODE_MEM_FROM_ML (vcvode_mem), Bool_val(vonoff));
    CHECK_FLAG ("CVodeSetProjErrEst", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_proj_frequency(value vcvode_mem, value vfreq)
{
    CAMLparam2(vcvode_mem, vfreq);
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeSetProjFrequency(CVODE_MEM_FROM_ML (vcvode_mem), Int_val(vfreq));
    CHECK_FLAG ("CVodeSetProjFrequency", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_max_num_proj_fails(value vcvode_mem, value vmfails)
{
    CAMLparam2(vcvode_mem, vmfails);
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeSetMaxNumProjFails(CVODE_MEM_FROM_ML (vcvode_mem),
				   Int_val(vmfails));
    CHECK_FLAG ("CVodeSetMaxNumProjFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_eps_proj(value vcvode_mem, value veps)
{
    CAMLparam2(vcvode_mem, veps);
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeSetEpsProj(CVODE_MEM_FROM_ML (vcvode_mem), Double_val(veps));
    CHECK_FLAG ("CVodeSetEpsProj", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_proj_fail_eta(value vcvode_mem, value veps)
{
    CAMLparam2(vcvode_mem, veps);
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeSetProjFailEta(CVODE_MEM_FROM_ML (vcvode_mem), Double_val(veps));
    CHECK_FLAG ("CVodeSetProjFailEta", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_no_inactive_root_warn(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag = CVodeSetNoInactiveRootWarn(CVODE_MEM_FROM_ML(vcvode_mem));
    CHECK_FLAG("CVodeSetNoInactiveRootWarn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_gs_type(value vcvode_mem, value vgstype)
{
    CAMLparam2(vcvode_mem, vgstype);
#if SUNDIALS_LIB_VERSION < 300

    int flag = CVSpilsSetGSType(CVODE_MEM_FROM_ML(vcvode_mem),
				sunml_lsolver_gs_type(vgstype));
    CHECK_SPILS_FLAG("CVSpilsSetGSType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_jac_eval_frequency(value vcvode_mem,
						  value vfreq)
{
    CAMLparam2(vcvode_mem, vfreq);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetJacEvalFrequency(CVODE_MEM_FROM_ML(vcvode_mem),
					Long_val(vfreq));
    CHECK_LS_FLAG("CVodeSetJacEvalFrequency", flag);
#elif 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetMaxStepsBetweenJac(CVODE_MEM_FROM_ML(vcvode_mem),
					  Long_val(vfreq));
    CHECK_LS_FLAG("CVodeSetMaxStepsBetweenJac", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_lsetup_frequency(value vcvode_mem, value varg)
{
    CAMLparam2(vcvode_mem, varg);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetLSetupFrequency(CVODE_MEM_FROM_ML(vcvode_mem),
				       Int_val(varg));
    CHECK_FLAG("CVodeSetLSetupFrequency", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_linear_solution_scaling(value vcvode_mem,
						       value vonoff)
{
    CAMLparam2(vcvode_mem, vonoff);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetLinearSolutionScaling(CVODE_MEM_FROM_ML(vcvode_mem),
					     Bool_val(vonoff));
    CHECK_FLAG("CVodeSetLinearSolutionScaling", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_eps_lin(value vcvode_mem, value eplifac)
{
    CAMLparam2(vcvode_mem, eplifac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetEpsLin(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(eplifac));
    CHECK_LS_FLAG("CVodeSetEpsLin", flag);
#else
    int flag = CVSpilsSetEpsLin(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(eplifac));
    CHECK_SPILS_FLAG("CVSpilsSetEpsLin", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_set_ls_norm_factor(value vcvode_mem, value vfac)
{
    CAMLparam2(vcvode_mem, vfac);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetLSNormFactor(CVODE_MEM_FROM_ML(vcvode_mem),
				    Double_val(vfac));
    CHECK_LS_FLAG("CVodeSetLSNormFactor", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_maxl(value vcvode_mem, value maxl)
{
    CAMLparam2(vcvode_mem, maxl);
#if SUNDIALS_LIB_VERSION < 300
    int flag = CVSpilsSetMaxl(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxl));
    CHECK_SPILS_FLAG("CVSpilsSetMaxl", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* statistic accessor functions */

CAMLprim value sunml_cvode_get_num_stab_lim_order_reds(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumStabLimOrderReds(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumStabLimOrderReds", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_tol_scale_factor(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    sunrealtype r;
    int flag = CVodeGetTolScaleFactor(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetTolScaleFactor", flag);

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_cvode_get_num_nonlin_solv_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_num_nonlin_solv_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_nonlin_solv_stats(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int nniters, nncfails;
    int flag = CVodeGetNonlinSolvStats(CVODE_MEM_FROM_ML(vcvode_mem),
				       &nniters, &nncfails);
    CHECK_FLAG("CVodeGetNonlinSolvStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nniters));
    Store_field(r, 1, Val_long(nncfails));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_get_num_g_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumGEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumGEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_num_proj_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    long int nproj;
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeGetNumProjEvals(CVODE_MEM_FROM_ML (vcvode_mem), &nproj);
    CHECK_FLAG ("CVodeGetNumProjEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(nproj));
}

CAMLprim value sunml_cvode_get_num_proj_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    long int npfails;
#if 530 <= SUNDIALS_LIB_VERSION && !SUNDIALSML_WITHSENS
    int flag;

    flag = CVodeGetNumProjFails(CVODE_MEM_FROM_ML (vcvode_mem), &npfails);
    CHECK_FLAG ("CVodeGetNumProjFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(npfails));
}

CAMLprim value sunml_cvode_dls_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetLinWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem),
				    &lenrwLS, &leniwLS);
    CHECK_LS_FLAG("CVodeGetLinWorkSpace", flag);
#else
    int flag = CVDlsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem),
				 &lenrwLS, &leniwLS);
    CHECK_DLS_FLAG("CVDlsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumJacEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumJacEvals", flag);
#else
    int flag = CVDlsGetNumJacEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_DLS_FLAG("CVDlsGetNumJacEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_dls_get_num_lin_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumLinRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumLinRhsEvals", flag);
#else
    int flag = CVDlsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_DLS_FLAG("CVDlsGetNumRhsEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_diag_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = CVDiagGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrwLS, &leniwLS);
    CHECK_FLAG("CVDiagGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_diag_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDiagGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVDiagGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_bandprec_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwBP;
    long int leniwBP;

    int flag = CVBandPrecGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrwBP, &leniwBP);
    CHECK_FLAG("CVBandPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrwBP));
    Store_field(r, 1, Val_long(leniwBP));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_bandprec_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVBandPrecGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVBandPrecGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

/* spils functions */

CAMLprim value sunml_cvode_spils_spgmr (value vcvode_mem, value vmaxl, value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVSpgmr (cvode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_spbcgs (value vcvode_mem, value vmaxl, value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVSpbcg (cvode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpbcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_sptfqmr (value vcvode_mem, value vmaxl,
				      value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVSptfqmr (cvode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSptfqmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_get_num_lin_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumLinIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumLinIters", flag);
#else
    int flag = CVSpilsGetNumLinIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumLinIters", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_num_lin_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumLinConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumLinConvFails", flag);
#else
    int flag = CVSpilsGetNumConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumConvFails", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_spils_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = CVodeGetLinWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_LS_FLAG("CVodeGetLinWorkSpace", flag);
#else
    flag = CVSpilsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_SPILS_FLAG("CVSpilsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_get_num_prec_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumPrecEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumPrecEvals", flag);
#else
    int flag = CVSpilsGetNumPrecEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumPrecEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_num_prec_solves(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumPrecSolves(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumPrecSolves", flag);
#else
    int flag = CVSpilsGetNumPrecSolves(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumPrecSolves", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_num_jtsetup_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumJTSetupEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumJTSetupEvals", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    int flag = CVSpilsGetNumJTSetupEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumJTSetupEvals", flag);
#else
    r = 0;
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_num_jtimes_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumJtimesEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumJtimesEvals", flag);
#else
    int flag = CVSpilsGetNumJtimesEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumJtimesEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_get_num_lin_rhs_evals (value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeGetNumLinRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_LS_FLAG("CVodeGetNumLinRhsEvals", flag);
#else
    int flag = CVSpilsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumRhsEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

