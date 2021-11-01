/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2015 Institut National de Recherche en Informatique et   *
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

/* ARKODE */

#include <arkode/arkode.h>

#if 400 <= SUNDIALS_LIB_VERSION
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_mristep.h>
#endif

/* linear solvers */
#include <arkode/arkode_bandpre.h>
#if   400 <= SUNDIALS_LIB_VERSION
#include <arkode/arkode_ls.h>
#elif 300 <= SUNDIALS_LIB_VERSION
#include <arkode/arkode_direct.h>
#include <arkode/arkode_spils.h>
#else
#include <arkode/arkode_dense.h>
#include <arkode/arkode_band.h>
#include <arkode/arkode_spgmr.h>
#include <arkode/arkode_spfgmr.h>
#include <arkode/arkode_spbcgs.h>
#include <arkode/arkode_sptfqmr.h>
#include <arkode/arkode_pcg.h>
#endif

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <arkode/arkode_lapack.h>
#endif

#include "../lsolvers/sundials_matrix_ml.h"
#include "../lsolvers/sundials_linearsolver_ml.h"
#include "../lsolvers/sundials_nonlinearsolver_ml.h"
#include "../sundials/sundials_ml.h"
#include "arkode_ml.h"
#include "../nvectors/nvector_ml.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

#if 580 <= SUNDIALS_LIB_VERSION
#define ISTEPPER(v) (*(MRIStepInnerStepper*)Data_custom_val(v))
#define ISTEPPER_FROM_ML(v) (ISTEPPER(Field((v), RECORD_ARKODE_MRI_ISTEPPER_RAWPTR)))
#endif

// must correspond with the constructors of Arkode.interpolant_type
static int ark_interpolant_types[] = {
    ARK_INTERP_HERMITE,
    ARK_INTERP_LAGRANGE,
};

CAMLprim value sunml_arkode_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (ARKODE, exns);
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
    value r = caml_callback_exn (Field(session, RECORD_ARKODE_SESSION_ERRH), a);
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined error handler");

    CAMLreturn0;
}

CAMLprim value sunml_arkode_ark_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), errh,
				      ARKODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("ARKStepSetErrHandlerFn", flag);
#else
    int flag = ARKodeSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), errh,
				     ARKODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("ARKodeSetErrHandlerFn", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("ARKStepSetErrHandlerFn", flag);
#else
    int flag = ARKodeSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("ARKodeSetErrHandlerFn", flag);
#endif

    CAMLreturn (Val_unit);
}

int sunml_arkode_translate_exception (value session, value exn,
			        recoverability recoverable)
{
    CAMLparam2(session, exn);
    CAMLlocal1(bucket);

    if (recoverable && Field(exn, 0) == SUNDIALS_EXN_TAG (RecoverableFailure))
	CAMLreturnT (int, 1);

    /* Unrecoverable error.  Save the exception and return -1.  */
    bucket = caml_alloc_small (1,0);
    Field (bucket, 0) = exn;
    Store_field (session, RECORD_ARKODE_SESSION_EXN_TEMP, bucket);
    CAMLreturnT (int, -1);
}

static int istepper_translate_exception(value exn, recoverability recoverable)
{
    CAMLparam1(exn);

    if ((recoverable & 0x1)
	    && Field(exn, 0) == SUNDIALS_EXN_TAG (RecoverableFailure))
	CAMLreturnT (int, 1);

    if (recoverable) sunml_warn_discarded_exn (exn, "nonlinear solver");

    /* Unrecoverable error -1. Unfortunately we lose the exception. */
    CAMLreturnT (int, -1);
}

static int rhsfn1(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_ARKODE_SESSION_RHSFN1),
				 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int rhsfn2(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_ARKODE_SESSION_RHSFN2),
				 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if 580 <= SUNDIALS_LIB_VERSION
static int nlsrhsfn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_ARKODE_SESSION_NLS_RHSFN),
				 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    value *backref = user_data;
    intnat nroots;

    WEAK_DEREF (session, *backref);

    nroots = ARKODE_NROOTS_FROM_ML (session);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);
    args[2] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(session, RECORD_ARKODE_SESSION_ROOTSFN),
				  3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, UNRECOVERABLE));
}

#if 500 <= SUNDIALS_LIB_VERSION
static int preinnerfn(realtype t, N_Vector *f, int num_vecs, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 2);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = Field(session, RECORD_ARKODE_SESSION_PREINNERARRAY);

    if (Wosize_val(args[1]) != num_vecs) {
	args[1] = caml_alloc_tuple(num_vecs);
    }
    sunml_nvectors_into_array(num_vecs, args[1], f);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session,
				    RECORD_ARKODE_SESSION_PREINNERFN),
				 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static int postinnerfn(realtype t, N_Vector y, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 2);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session,
				    RECORD_ARKODE_SESSION_POSTINNERFN),
				 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static int stagepredictfn(realtype t, N_Vector zpred, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 2);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(zpred);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session,
				    RECORD_ARKODE_SESSION_STAGEPREDICTFN),
				 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocal1 (session);
    value *backref = user_data;

    WEAK_DEREF (session, *backref);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (Field(session, RECORD_ARKODE_SESSION_ERRW),
				  NVEC_BACKLINK (y), NVEC_BACKLINK (ewt));
    if (Is_exception_result (r)) {
	r = Extract_exception (r);
	if (Field (r, 0) != SUNDIALS_EXN_TAG (NonPositiveEwt))
	    sunml_warn_discarded_exn (r, "user-defined error weight fun");
	CAMLreturnT (int, -1);
    }

    CAMLreturnT (int, 0);
}

static int resw(N_Vector y, N_Vector rwt, void *user_data)
{
    CAMLparam0();
    CAMLlocal1 (session);
    value *backref = user_data;

    WEAK_DEREF (session, *backref);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (Field(session, RECORD_ARKODE_SESSION_RESW),
				  NVEC_BACKLINK (y), NVEC_BACKLINK (rwt));
    if (Is_exception_result (r)) {
	r = Extract_exception (r);
	if (Field (r, 0) != SUNDIALS_EXN_TAG (NonPositiveEwt))
	    sunml_warn_discarded_exn (r,
		    "user-defined residual weight fun");
	CAMLreturnT (int, -1);
    }

    CAMLreturnT (int, 0);
}

static int adaptfn(N_Vector y, realtype t,
		   realtype h1, realtype h2, realtype h3, 
		   realtype e1, realtype e2, realtype e3,
		   int q, int p, 
		   realtype *hnew, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    value *backref = user_data;
    WEAK_DEREF (session, *backref);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);
    args[2] = caml_alloc_tuple(RECORD_ARKODE_ADAPTIVITY_ARGS_SIZE);

    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_H1, caml_copy_double(h1));
    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_H2, caml_copy_double(h2));
    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_H3, caml_copy_double(h3));
    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_E1, caml_copy_double(e1));
    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_E2, caml_copy_double(e2));
    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_E3, caml_copy_double(e3));
    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_Q,  Val_int(q));
    Store_field(args[2],RECORD_ARKODE_ADAPTIVITY_ARGS_P,  Val_int(p));

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(session, RECORD_ARKODE_SESSION_ADAPTFN),
				  3, args);

    /* Update hnew; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*hnew = Double_val (r);
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}

static int stabfn(N_Vector y, realtype t, realtype *hstab, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 2);

    value *backref = user_data;
    WEAK_DEREF (session, *backref);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(session, RECORD_ARKODE_SESSION_STABFN),
				  2, args);

    /* Update hstab; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*hstab = Double_val (r);
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}

static int resizefn(N_Vector y, N_Vector ytemplate, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 2);

    value *backref = user_data;
    WEAK_DEREF (session, *backref);

    args[0] = NVEC_BACKLINK (y);
    args[1] = NVEC_BACKLINK (ytemplate);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_ARKODE_SESSION_RESIZEFN),
				 2, args);

    CAMLreturnT(int, Is_exception_result(r));
}

#if 270 <= SUNDIALS_LIB_VERSION
static int poststepfn(realtype t, N_Vector y, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 2);

    value *backref = user_data;
    WEAK_DEREF (session, *backref);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_ARKODE_SESSION_POSTSTEPFN),
				 2, args);

    CAMLreturnT(int, Is_exception_result(r));
}
#endif

value sunml_arkode_make_jac_arg(realtype t, N_Vector y, N_Vector fy, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_ARKODE_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_ARKODE_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_ARKODE_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_ARKODE_JACOBIAN_ARG_JAC_FY, NVEC_BACKLINK(fy));
    Store_field(r, RECORD_ARKODE_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

value sunml_arkode_make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
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

static int jacfn(realtype t,
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

    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_arkode_make_jac_arg (t, y, fy,
				  sunml_arkode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = MAT_BACKLINK(Jac);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#else // SUNDIALS_LIB_VERSION < 300

/* Dense and band Jacobians only work with serial NVectors.  */
static int jacfn(long int n,
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
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);

    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(Jac));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_arkode_make_jac_arg (t, y, fy,
				  sunml_arkode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bandjacfn(long int N,
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
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(Jac));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_arkode_make_jac_arg(t, y, fy,
				 sunml_arkode_make_triple_tmp(tmp1, tmp2, tmp3));
    args[1] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#endif

#if 500 <= SUNDIALS_LIB_VERSION
static int linsysfn(
	realtype t,
	N_Vector y,
	N_Vector fy,
	SUNMatrix A,
	SUNMatrix M,
	booleantype jok,
	booleantype *jcur,
	realtype gamma,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 5);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);

    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1);

    args[0] = sunml_arkode_make_jac_arg (t, y, fy,
				  sunml_arkode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = MAT_BACKLINK(A);
    args[2] = M == NULL ? Val_none : Some_val(MAT_BACKLINK(M));
    args[3] = Val_bool(jok);
    args[4] = caml_copy_double(gamma);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 5, args);

    /* Update jcur; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcur = Bool_val (r);
	CAMLreturnT(int, 0);
    }
    r = Extract_exception (r);

    CAMLreturnT(int, sunml_arkode_translate_exception (session, r, RECOVERABLE));
}
#endif

static int precsetupfn(realtype t,
		       N_Vector y,
		       N_Vector fy,
		       booleantype jok,
		       booleantype *jcurPtr,
		       realtype gamma,
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

    args[0] = sunml_arkode_make_jac_arg(t, y, fy, Val_unit);
    args[1] = Val_bool(jok);
    args[2] = caml_copy_double(gamma);

    cb = ARKODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_SPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    /* Update jcurPtr; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcurPtr = Bool_val (r);
	CAMLreturnT(int, 0);
    }
    r = Extract_exception (r);

    CAMLreturnT(int, sunml_arkode_translate_exception (session, r, RECOVERABLE));
}

static value make_spils_solve_arg(N_Vector r,
				  realtype gamma,
				  realtype delta,
				  int lr)
{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(RECORD_ARKODE_SPILS_SOLVE_ARG_SIZE);
    Store_field(v, RECORD_ARKODE_SPILS_SOLVE_ARG_RHS, NVEC_BACKLINK(r));
    Store_field(v, RECORD_ARKODE_SPILS_SOLVE_ARG_GAMMA,
                caml_copy_double(gamma));
    Store_field(v, RECORD_ARKODE_SPILS_SOLVE_ARG_DELTA,
                caml_copy_double(delta));
    Store_field(v, RECORD_ARKODE_SPILS_SOLVE_ARG_LEFT, Val_bool (lr == 1));

    CAMLreturn(v);
}

static int precsolvefn(realtype t,
		       N_Vector y,
		       N_Vector fy,
		       N_Vector rvec,
		       N_Vector z,
		       realtype gamma,
		       realtype delta,
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

    args[0] = sunml_arkode_make_jac_arg(t, y, fy, Val_unit);
    args[1] = make_spils_solve_arg(rvec, gamma, delta, lr);
    args[2] = NVEC_BACKLINK(z);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_LS_PRECFNS_FROM_ML(session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_SPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int jactimesfn(N_Vector v,
		      N_Vector Jv,
		      realtype t,
		      N_Vector y,
		      N_Vector fy,
		      void *user_data,
		      N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = sunml_arkode_make_jac_arg(t, y, fy, NVEC_BACKLINK(tmp));
    args[1] = NVEC_BACKLINK(v);
    args[2] = NVEC_BACKLINK(Jv);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if 300 <= SUNDIALS_LIB_VERSION
static int jacsetupfn(realtype t,
		      N_Vector y,
		      N_Vector fy,
		      void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = sunml_arkode_make_jac_arg(t, y, fy, Val_unit);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, arg);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

#if 530 <= SUNDIALS_LIB_VERSION
static int jactimesrhsfn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_arkode_dls_dense (value varkode_mem,
				       value vneqs,
				       value vset_jac)
{
    CAMLparam3(varkode_mem, vneqs, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val(vneqs);
    int flag;

    flag = ARKDense (arkode_mem, neqs);
    CHECK_FLAG ("ARKDense", flag);
    if (Bool_val (vset_jac)) {
	flag = ARKDlsSetDenseJacFn(ARKODE_MEM_FROM_ML(varkode_mem), jacfn);
	CHECK_FLAG("ARKDlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_lapack_dense (value varkode_mem,
					      value vneqs,
					      value vset_jac)
{
    CAMLparam3 (varkode_mem, vneqs, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = ARKLapackDense (arkode_mem, neqs);
    CHECK_FLAG ("ARKLapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = ARKDlsSetDenseJacFn (arkode_mem, jacfn);
	CHECK_FLAG("ARKDlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_band (value varkode_mem, value vneqs,
				      value vmupper, value vmlower,
				      value vset_jac)
{
    CAMLparam5(varkode_mem, vneqs, vmupper, vmlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = ARKBand (arkode_mem, neqs, Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("ARKBand", flag);
    if (Bool_val (vset_jac)) {
	flag = ARKDlsSetBandJacFn(arkode_mem, bandjacfn);
	CHECK_FLAG("ARKDlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_lapack_band (value varkode_mem, value vneqs,
					     value vmupper, value vmlower,
					     value vset_jac)
{
    CAMLparam5(varkode_mem, vneqs, vmupper, vmlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val(vneqs);
    int flag;

    flag = ARKLapackBand (arkode_mem, neqs,
			  Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("ARKLapackBand", flag);
    if (Bool_val (vset_jac)) {
	flag = ARKDlsSetBandJacFn(arkode_mem, bandjacfn);
	CHECK_FLAG("ARKDlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_linear_solver (value varkode_mem,
						   value vlsolv,
						   value vojmat,
						   value vhasjac,
						   value vhaslsf)
{
    CAMLparam5(varkode_mem, vlsolv, vojmat, vhasjac, vhaslsf);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = ARKStepSetLinearSolver(arkode_mem, lsolv, jmat);
    CHECK_LS_FLAG ("ARKStepSetLinearSolver", flag);
    flag = ARKStepSetJacFn(arkode_mem, Bool_val(vhasjac) ? jacfn : NULL);
    CHECK_LS_FLAG("ARKStepSetJacFn", flag);

#if 500 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vhaslsf)) {
	ARKStepSetLinSysFn(arkode_mem, linsysfn);
	CHECK_LS_FLAG("ARKStepSetLinSysFn", flag);
    }
#endif

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_linear_solver (value varkode_mem,
						   value vlsolv,
						   value vojmat,
						   value vhasjac,
						   value vhaslsf)
{
    CAMLparam5(varkode_mem, vlsolv, vojmat, vhasjac, vhaslsf);
#if 540 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = MRIStepSetLinearSolver(arkode_mem, lsolv, jmat);
    CHECK_LS_FLAG ("MRIStepSetLinearSolver", flag);
    flag = MRIStepSetJacFn(arkode_mem, Bool_val(vhasjac) ? jacfn : NULL);
    CHECK_LS_FLAG("MRIStepSetJacFn", flag);

    if (Bool_val(vhaslsf)) {
	MRIStepSetLinSysFn(arkode_mem, linsysfn);
	CHECK_LS_FLAG("MRIStepSetLinSysFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_set_linear_solver (value varkode_mem, value vlsolv,
					       value vjmat, value vhasjac)
{
    CAMLparam4(varkode_mem, vlsolv, vjmat, vhasjac);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = ARKDlsSetLinearSolver(arkode_mem, lsolv, jmat);
    CHECK_FLAG ("ARKDlsSetLinearSolver", flag);
    flag = ARKDlsSetJacFn(arkode_mem, Bool_val(vhasjac) ? jacfn : NULL);
    CHECK_FLAG("ARKDlsSetJacFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_linear_solver (value varkode_mem,
						     value vlsolv)
{
    CAMLparam2(varkode_mem, vlsolv);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    int flag;

    flag = ARKSpilsSetLinearSolver(arkode_mem, lsolv);
    CHECK_FLAG ("ARKSpilsSetLinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_preconditioner (value vsession,
						    value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = ARKODE_MEM_FROM_ML (vsession);
#if 400 <= SUNDIALS_LIB_VERSION
    ARKLsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = ARKStepSetPreconditioner (mem, setup, precsolvefn);
    CHECK_FLAG ("ARKStepSetPreconditioner", flag);
#else
    ARKSpilsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = ARKSpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_FLAG ("ARKSpilsSetPreconditioner", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_banded_preconditioner (value vsession,
						       value vneqs,
						       value vmupper,
						       value vmlower)
{
    CAMLparam4 (vsession, vneqs, vmupper, vmlower);
    long neqs = Index_val (vneqs);
    int flag = ARKBandPrecInit (ARKODE_MEM_FROM_ML (vsession), neqs,
			        Index_val (vmupper), Index_val (vmlower));
    CHECK_FLAG ("ARKBandPrecInit", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_jac_times(value vdata, value vhas_setup,
					      value vhas_times)
{
    CAMLparam3(vdata, vhas_setup, vhas_times);
#if 400 <= SUNDIALS_LIB_VERSION
    ARKLsJacTimesSetupFn setup = Bool_val (vhas_setup) ? jacsetupfn : NULL;
    ARKLsJacTimesVecFn   times = Bool_val (vhas_times) ? jactimesfn : NULL;

    int flag = ARKStepSetJacTimes(ARKODE_MEM_FROM_ML(vdata), setup, times);
    CHECK_LS_FLAG("ARKStepSetJacTimes", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    ARKSpilsJacTimesSetupFn setup = Bool_val (vhas_setup) ? jacsetupfn : NULL;
    ARKSpilsJacTimesVecFn   times = Bool_val (vhas_times) ? jactimesfn : NULL;

    int flag = ARKSpilsSetJacTimes(ARKODE_MEM_FROM_ML(vdata), setup, times);
    CHECK_FLAG("ARKSpilsSetJacTimes", flag);
#else
    ARKSpilsJacTimesVecFn jac = Bool_val (vhas_times) ? jactimesfn : NULL;
    int flag = ARKSpilsSetJacTimesVecFn(ARKODE_MEM_FROM_ML(vdata), jac);
    CHECK_FLAG("ARKSpilsSetJacTimesVecFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_jac_times_rhsfn(value vdata,
						    value vhas_rhsfn)
{
    CAMLparam2(vdata, vhas_rhsfn);
#if 530 <= SUNDIALS_LIB_VERSION
    ARKRhsFn rhsfn = Bool_val (vhas_rhsfn) ? jactimesrhsfn : NULL;

    int flag = ARKStepSetJacTimesRhsFn(ARKODE_MEM_FROM_ML(vdata), rhsfn);
    CHECK_LS_FLAG("ARKStepSetJacTimesRhsFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

#if 400 <= SUNDIALS_LIB_VERSION
// hack to work around lack of CVodeGetUserData
typedef struct {
  realtype uround;
  void *user_data;
  //...
} *StartOf_ARKodeMem;

static value sunml_arkode_ark_session_to_value(void *arkode_mem)
{
    value session;
    // void *user_data = ARKStepGetUserData(arkode_mem);
    void *user_data = ((StartOf_ARKodeMem)arkode_mem)->user_data;

    WEAK_DEREF (session, *(value*)user_data);
    return session;
}

static value sunml_arkode_mri_session_to_value(void *arkode_mem)
{
    value session;
    // void *user_data = MRIStepGetUserData(arkode_mem);
    void *user_data = ((StartOf_ARKodeMem)arkode_mem)->user_data;

    WEAK_DEREF (session, *(value*)user_data);
    return session;
}

static void* sunml_arkode_session_from_value(value varkode_mem)
{
    return (ARKODE_MEM_FROM_ML(varkode_mem));
}
#endif

CAMLprim value sunml_arkode_ark_set_nonlinear_solver(value varkode_mem,
						     value vnlsolv)
{
    CAMLparam2(varkode_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_arkode_ark_session_to_value,
				   sunml_arkode_session_from_value);

    int flag = ARKStepSetNonlinearSolver(arkode_mem, nlsolv);
    CHECK_FLAG ("ARKStepSetNonlinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_nonlinear_solver(value varkode_mem,
						     value vnlsolv)
{
    CAMLparam2(varkode_mem, vnlsolv);
#if 540 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_arkode_mri_session_to_value,
				   sunml_arkode_session_from_value);

    int flag = MRIStepSetNonlinearSolver(arkode_mem, nlsolv);
    CHECK_FLAG ("MRIStepSetNonlinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_wf_tolerances (value vdata)
{
    CAMLparam1(vdata);
 
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepWFtolerances(ARKODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("ARKStepWFtolerances", flag);
#else
    int flag = ARKodeWFtolerances(ARKODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("ARKodeWFtolerances", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_wf_tolerances (value vdata)
{
    CAMLparam1(vdata);
 
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepWFtolerances(ARKODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("MRIStepWFtolerances", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_wf_tolerances (value vdata)
{
    CAMLparam1(vdata);
 
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepWFtolerances(ARKODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("ERKStepWFtolerances", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

/* Mass matrix routines. */

#if 300 <= SUNDIALS_LIB_VERSION

static int massfn(realtype t,
		  SUNMatrix M,
		  void *user_data,
		  N_Vector tmp1,
		  N_Vector tmp2,
		  N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 3);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);

    cb = ARKODE_MASS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = caml_copy_double(t);
    args[1] = sunml_arkode_make_triple_tmp (tmp1, tmp2, tmp3);
    args[2] = MAT_BACKLINK(M);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#else // SUNDIALS_LIB_VERSION < 300

static int massfn(long int n,
		  realtype t,
		  DlsMat M,
		  void *user_data,
		  N_Vector tmp1,
		  N_Vector tmp2,
		  N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 3);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);

    cb = ARKODE_MASS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(M));
	Store_field(cb, 1, dmat);
    }

    args[0] = caml_copy_double(t);
    args[1] = sunml_arkode_make_triple_tmp (tmp1, tmp2, tmp3);
    args[2] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bandmassfn(long int N,
		      long int mupper,
		      long int mlower,
		      realtype t,
		      DlsMat M,
		      void *user_data, 	 
		      N_Vector tmp1,
		      N_Vector tmp2,
		      N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_MASS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(M));
	Store_field(cb, 1, bmat);
    }

    args[0] = caml_alloc_tuple(RECORD_ARKODE_BANDRANGE_SIZE);
    Store_field(args[0], RECORD_ARKODE_BANDRANGE_MUPPER, Val_long(mupper));
    Store_field(args[0], RECORD_ARKODE_BANDRANGE_MLOWER, Val_long(mlower));
    args[1] = caml_copy_double(t);
    args[2] = sunml_arkode_make_triple_tmp(tmp1, tmp2, tmp3);
    args[3] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#endif

static int masstimesfn(N_Vector v,
		       N_Vector Mv,
		       realtype t,
		       void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(v);
    args[2] = NVEC_BACKLINK(Mv);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_MASS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if 300 <= SUNDIALS_LIB_VERSION
static int masssetupfn(realtype t,
		       void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_MASS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, caml_copy_double(t));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int massprecsetupfn(realtype t,
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

    WEAK_DEREF (session, *(value*)user_data);

    cb = ARKODE_MASS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_SPILS_MASS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, caml_copy_double(t));

    /* Update jcurPtr; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	CAMLreturnT(int, 0);
    }

    r = Extract_exception (r);
    CAMLreturnT(int, sunml_arkode_translate_exception (session, r, RECOVERABLE));
}

static value make_spils_mass_solve_arg(N_Vector r,
				       realtype delta,
				       int lr)
{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_SIZE);
    Store_field(v, RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_RHS, NVEC_BACKLINK(r));
    Store_field(v, RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_DELTA,
                caml_copy_double(delta));
    Store_field(v, RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_LEFT, Val_bool (lr == 1));

    CAMLreturn(v);
}

static int massprecsolvefn(realtype t,
			   N_Vector rvec,
			   N_Vector z,
			   realtype delta,
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

    args[0] = caml_copy_double(t);
    args[1] = make_spils_mass_solve_arg(rvec, delta, lr);
    args[2] = NVEC_BACKLINK(z);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_MASS_PRECFNS_FROM_ML(session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_SPILS_MASS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_arkode_dls_mass_dense (value varkode_mem, value vneqs)
{
    CAMLparam2(varkode_mem, vneqs);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val(vneqs);
    int flag;

    flag = ARKMassDense (arkode_mem, neqs, massfn);
    CHECK_FLAG ("ARKMassDense", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_mass_lapack_dense (value varkode_mem, value vneqs)
{
    CAMLparam2 (varkode_mem, vneqs);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = ARKMassLapackDense (arkode_mem, neqs, massfn);
    CHECK_FLAG ("ARKMassLapackDense", flag);
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_mass_band (value varkode_mem, value vneqs,
				       value vmupper, value vmlower)
{
    CAMLparam4(varkode_mem, vneqs, vmupper, vmlower);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = ARKMassBand (arkode_mem, neqs,
				    Long_val (vmupper),
				    Long_val (vmlower),
				    bandmassfn);
    CHECK_FLAG ("ARKMassBand", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_mass_lapack_band (value varkode_mem, value vneqs,
					      value vmupper, value vmlower)
{
    CAMLparam4(varkode_mem, vneqs, vmupper, vmlower);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    long neqs = Long_val(vneqs);
    int flag;

    flag = ARKMassLapackBand (arkode_mem, neqs,
			                  Long_val (vmupper),
					  Long_val (vmlower),
					  bandmassfn);
    CHECK_FLAG ("ARKMassLapackBand", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_mass_linear_solver (value varkode_mem,
							value vlsolv,
					                value vojmat,
							value vtimedep)
{
    CAMLparam4(varkode_mem, vlsolv, vojmat, vtimedep);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = ARKStepSetMassLinearSolver(arkode_mem, lsolv, jmat,
				      Bool_val(vtimedep));
    CHECK_LS_FLAG ("ARKStepSetMassLinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_dls_set_mass_linear_solver (value varkode_mem,
				    value vlsolv, value vjmat, value vtime_dep)
{
    CAMLparam4(varkode_mem, vlsolv, vjmat, vtime_dep);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = ARKDlsSetMassLinearSolver(arkode_mem, lsolv, jmat,
				     Bool_val(vtime_dep));
    CHECK_FLAG ("ARKDlsSetMassLinearSolver", flag);
    flag = ARKDlsSetMassFn(arkode_mem, massfn);
    CHECK_FLAG("ARKDlsSetMassFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_mass_linear_solver (value varkode_mem,
							  value vlsolv,
							  value vhassetup,
							  value vtime_dep)
{
    CAMLparam4(varkode_mem, vlsolv, vhassetup, vtime_dep);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    int flag;

    flag = ARKSpilsSetMassLinearSolver(arkode_mem, lsolv, Bool_val(vtime_dep));
    CHECK_FLAG ("ARKSpilsSetMassLinearSolver", flag);
    flag = ARKSpilsSetMassTimes(
		    arkode_mem,
		    Bool_val(vhassetup) ? masssetupfn : NULL,
		    masstimesfn,
		    (void *)Field(arkode_mem, RECORD_ARKODE_SESSION_BACKREF));
    CHECK_FLAG ("ARKSpilsSetMassTimes", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_mass_preconditioner(value vsession,
						        value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = ARKODE_MEM_FROM_ML (vsession);
#if 400 <= SUNDIALS_LIB_VERSION
    ARKLsMassPrecSetupFn setup
	= Bool_val (vset_precsetup) ? massprecsetupfn : NULL;
    int flag = ARKStepSetMassPreconditioner (mem, setup, massprecsolvefn);
    CHECK_FLAG ("ARKStepSetMassPreconditioner", flag);
#else
    ARKSpilsMassPrecSetupFn setup
	= Bool_val (vset_precsetup) ? massprecsetupfn : NULL;
    int flag = ARKSpilsSetMassPreconditioner (mem, setup, massprecsolvefn);
    CHECK_FLAG ("ARKSpilsSetMassPreconditioner", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_mass_fn(value varkode_mem)
{
    CAMLparam1(varkode_mem);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMassFn(ARKODE_MEM_FROM_ML (varkode_mem), massfn);
    CHECK_LS_FLAG("ARKStepSetMassFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_mass_times(value varkode_mem,
					       value vhassetup)
{
    CAMLparam2(varkode_mem, vhassetup);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;
    flag = ARKStepSetMassTimes(
		    arkode_mem,
		    Bool_val(vhassetup) ? masssetupfn : NULL,
		    masstimesfn,
		    (void *)Field(arkode_mem, RECORD_ARKODE_SESSION_BACKREF));
    CHECK_FLAG ("ARKStepSetMassTimes", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;
    flag = ARKSpilsSetMassTimes(
		    arkode_mem,
		    Bool_val(vhassetup) ? masssetupfn : NULL,
		    masstimesfn,
		    (void *)Field(arkode_mem, RECORD_ARKODE_SESSION_BACKREF));
    CHECK_FLAG ("ARKSpilsSetMassTimes", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/** ARKode generic interface: ARKStep, ERKStep, MRIStep **/

static value sunml_arkode_last_lin_exception(void *arkode_mem)
{
#if 400 <= SUNDIALS_LIB_VERSION
    long int lsflag;

    // ARKStepGetLastLinFlag and MRIStepGetLastLinFlag share the same
    // underlying implementation
    if (arkode_mem != NULL
	  && ARKStepGetLastLinFlag(arkode_mem, &lsflag) == ARKLS_SUCCESS)
	return sunml_lsolver_exception_from_flag(lsflag);
#endif
    return Val_none;
}

static value sunml_arkode_ark_last_mass_exception(void *arkode_mem)
{
#if 400 <= SUNDIALS_LIB_VERSION
    long int mlsflag;

    if (arkode_mem != NULL
	  && ARKStepGetLastMassFlag(arkode_mem, &mlsflag) == ARKLS_SUCCESS)
	return sunml_lsolver_exception_from_flag(mlsflag);
#endif
    return Val_none;
}

#if 400 <= SUNDIALS_LIB_VERSION
static value sunml_stepper_exception_from_flag(int flag)
{
    CAMLparam0();
    CAMLlocal1(vro);

#if 400 <= SUNDIALS_LIB_VERSION
    switch (flag) {
	case ARK_ILL_INPUT:
	    Store_some(vro, ARKODE_EXN(IllInput));
	    break;

	case ARK_TOO_CLOSE:
	    Store_some(vro, ARKODE_EXN(TooClose));
	    break;

	case ARK_TOO_MUCH_WORK:
	    Store_some(vro, ARKODE_EXN(TooMuchWork));
	    break;

	case ARK_TOO_MUCH_ACC:
	    Store_some(vro, ARKODE_EXN(TooMuchAccuracy));
	    break;

	case ARK_ERR_FAILURE:
	    Store_some(vro, ARKODE_EXN(ErrFailure));
	    break;

	case ARK_CONV_FAILURE:
	    Store_some(vro, ARKODE_EXN(ConvergenceFailure));
	    break;

	case ARK_LINIT_FAIL:
	    Store_some(vro, ARKODE_EXN(LinearInitFailure));
	    break;

	case ARK_RHSFUNC_FAIL:
	    Store_some(vro, ARKODE_EXN(RhsFuncFailure));
	    break;

	case ARK_FIRST_RHSFUNC_ERR:
	    Store_some(vro, ARKODE_EXN(FirstRhsFuncFailure));
	    break;

	case ARK_REPTD_RHSFUNC_ERR:
	    Store_some(vro, ARKODE_EXN(RepeatedRhsFuncFailure));
	    break;

	case ARK_UNREC_RHSFUNC_ERR:
	    Store_some(vro, ARKODE_EXN(UnrecoverableRhsFuncFailure));
	    break;

	case ARK_RTFUNC_FAIL:
	    Store_some(vro, ARKODE_EXN(RootFuncFailure));
	    break;

	case ARK_POSTPROCESS_FAIL:
	    Store_some(vro, ARKODE_EXN(PostprocStepFailure));
	    break;

	case ARK_BAD_K:
	    Store_some(vro, ARKODE_EXN(BadK));
	    break;

	case ARK_BAD_T:
	    Store_some(vro, ARKODE_EXN(BadT));
	    break;

	case ARK_VECTOROP_ERR:
	    Store_some(vro, ARKODE_EXN(VectorOpErr));
	    break;

	default:
	    vro = Val_none;
    }
#endif

    CAMLreturn(vro);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static value sunml_arkode_mri_innerstep_exception(void *arkode_mem)
{
    CAMLparam0();
    CAMLlocal1(vr);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag;

    if (arkode_mem != NULL
	  && MRIStepGetLastInnerStepFlag(arkode_mem, &flag) == ARK_SUCCESS)
    {
	vr = sunml_stepper_exception_from_flag(flag);
    } else {
	vr = Val_none;
    }
#endif
    CAMLreturn(vr);
}
#endif

void sunml_arkode_check_flag(const char *call, int flag, void *arkode_mem)
{
    CAMLparam0();
    CAMLlocal1(va);
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == ARK_SUCCESS
	    || flag == ARK_ROOT_RETURN
	    || flag == ARK_TSTOP_RETURN) return;

    switch (flag) {
	case ARK_ILL_INPUT:
	    caml_raise_constant(ARKODE_EXN(IllInput));

	case ARK_TOO_CLOSE:
	    caml_raise_constant(ARKODE_EXN(TooClose));

	case ARK_TOO_MUCH_WORK:
	    caml_raise_constant(ARKODE_EXN(TooMuchWork));

	case ARK_TOO_MUCH_ACC:
	    caml_raise_constant(ARKODE_EXN(TooMuchAccuracy));

	case ARK_ERR_FAILURE:
	    caml_raise_constant(ARKODE_EXN(ErrFailure));

	case ARK_CONV_FAILURE:
	    caml_raise_constant(ARKODE_EXN(ConvergenceFailure));

	case ARK_LINIT_FAIL:
	    caml_raise_constant(ARKODE_EXN(LinearInitFailure));

	case ARK_LSETUP_FAIL:
	    // NB: Only ARKStep or MRIStep should return this error
	    caml_raise_with_arg(ARKODE_EXN(LinearSetupFailure),
				sunml_arkode_last_lin_exception(arkode_mem));

	case ARK_LSOLVE_FAIL:
	    // NB: Only ARKStep or MRIStep should return this error
	    caml_raise_with_arg(ARKODE_EXN(LinearSolveFailure),
				sunml_arkode_last_lin_exception(arkode_mem));

	case ARK_MASSINIT_FAIL:
	    caml_raise_constant(ARKODE_EXN(MassInitFailure));

	case ARK_MASSSETUP_FAIL:
	    // NB: Only ARKStep should return this error
	    caml_raise_with_arg(ARKODE_EXN(MassSetupFailure),
				sunml_arkode_ark_last_mass_exception(arkode_mem));

	case ARK_MASSSOLVE_FAIL:
	    // NB: Only ARKStep should return this error
	    caml_raise_with_arg(ARKODE_EXN(MassSolveFailure),
				sunml_arkode_ark_last_mass_exception(arkode_mem));

	case ARK_MASSMULT_FAIL:
	    caml_raise_constant(ARKODE_EXN(MassMultFailure));

	case ARK_RHSFUNC_FAIL:
	    caml_raise_constant(ARKODE_EXN(RhsFuncFailure));

	case ARK_FIRST_RHSFUNC_ERR:
	    caml_raise_constant(ARKODE_EXN(FirstRhsFuncFailure));

	case ARK_REPTD_RHSFUNC_ERR:
	    caml_raise_constant(ARKODE_EXN(RepeatedRhsFuncFailure));

	case ARK_UNREC_RHSFUNC_ERR:
	    caml_raise_constant(ARKODE_EXN(UnrecoverableRhsFuncFailure));

	case ARK_RTFUNC_FAIL:
	    caml_raise_constant(ARKODE_EXN(RootFuncFailure));

#if 270 <= SUNDIALS_LIB_VERSION
	case ARK_POSTPROCESS_FAIL:
	    caml_raise_constant(ARKODE_EXN(PostprocStepFailure));
#endif

	case ARK_BAD_K:
	    caml_raise_constant(ARKODE_EXN(BadK));

	case ARK_BAD_T:
	    caml_raise_constant(ARKODE_EXN(BadT));

#if 400 <= SUNDIALS_LIB_VERSION
	case ARK_NLS_INIT_FAIL:
	    caml_raise_constant(ARKODE_EXN(NonlinearInitFailure));

	case ARK_NLS_SETUP_FAIL:
	    caml_raise_constant(ARKODE_EXN(NonlinearSetupFailure));

	case ARK_NLS_SETUP_RECVR:
	    caml_raise_constant(ARKODE_EXN(NonlinearSetupRecoverable));

	case ARK_NLS_OP_ERR:
	    caml_raise_constant(ARKODE_EXN(NonlinearOperationError));

	case ARK_VECTOROP_ERR:
	    caml_raise_constant(ARKODE_EXN(VectorOpErr));
	
	case ARK_INNERSTEP_FAIL:
	    // NB: Only MRIStep should return this error
	    va = sunml_arkode_mri_innerstep_exception(arkode_mem);
	    caml_raise_with_arg(ARKODE_EXN_TAG(InnerStepFail), va);
#endif

	default:
#if 400 <= SUNDIALS_LIB_VERSION
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     ARKStepGetReturnFlagName(flag));
	    // ARKStepGetReturnFlagName, ERKStepGetReturnFlagName, and
	    // MRIStepGetReturnFlagName share the same underlying
	    // implementation.
#else
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     ARKodeGetReturnFlagName(flag));
#endif
	    caml_failwith(exmsg);
    }

    CAMLreturn0;
}

/** ARKStep basic interface **/

/* ARKStepCreate() */
CAMLprim value sunml_arkode_ark_init(value weakref, value hasfi, value hasfe,
			             value y0, value t0)
{
    CAMLparam5(weakref, hasfi, hasfe, y0, t0);
    CAMLlocal2(r, varkode_mem);

    value *backref;

#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector nv_y0 = NVEC_VAL(y0);
    void *arkode_mem = ARKStepCreate(Bool_val(hasfe)  ? rhsfn2 : NULL,
				     Bool_val(hasfi)  ? rhsfn1 : NULL,
				     Double_val(t0),
				     nv_y0);

    if (arkode_mem == NULL)
	caml_failwith("ARKStepCreate returned NULL");

    varkode_mem = caml_alloc_final(1, NULL, 1, 5);
    ARKODE_MEM(varkode_mem) = arkode_mem;

    backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	ARKStepFree (&arkode_mem);
	caml_raise_out_of_memory();
    }
    ARKStepSetUserData (arkode_mem, backref);
#else
    int flag;
    void *arkode_mem = ARKodeCreate();
    if (arkode_mem == NULL)
	caml_failwith("ARKodeCreate returned NULL");

    varkode_mem = caml_alloc_final(1, NULL, 1, 5);
    ARKODE_MEM(varkode_mem) = arkode_mem;

    N_Vector nv_y0 = NVEC_VAL(y0);
    flag = ARKodeInit(arkode_mem,
		      Bool_val(hasfe)  ? rhsfn2 : NULL,
		      Bool_val(hasfi)  ? rhsfn1 : NULL,
		      Double_val(t0),
		      nv_y0);
    if (flag != ARK_SUCCESS) {
	ARKodeFree (&arkode_mem);
	CHECK_FLAG("ARKodeInit", flag);
    }

    backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	ARKodeFree (&arkode_mem);
	caml_raise_out_of_memory();
    }
    ARKodeSetUserData (arkode_mem, backref);
#endif

    r = caml_alloc_tuple (2);
    Store_field (r, 0, varkode_mem);
    Store_field (r, 1, (value)backref);

    CAMLreturn(r);
}

/* Set the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim value sunml_arkode_ark_root_init (value vdata, value vnroots)
{
    CAMLparam2 (vdata, vnroots);
    void *arkode_mem = ARKODE_MEM_FROM_ML (vdata);
    int nroots = Int_val (vnroots);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepRootInit (arkode_mem, nroots, roots);
    CHECK_FLAG ("ARKStepRootInit", flag);
#else
    int flag = ARKodeRootInit (arkode_mem, nroots, roots);
    CHECK_FLAG ("ARKodeRootInit", flag);
#endif
    Store_field (vdata, RECORD_ARKODE_SESSION_NROOTS, vnroots);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_sv_tolerances(value vdata, value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSVtolerances(ARKODE_MEM_FROM_ML(vdata),
			 	   Double_val(reltol), atol_nv);
    CHECK_FLAG("ARKStepSVtolerances", flag);
#else
    int flag = ARKodeSVtolerances(ARKODE_MEM_FROM_ML(vdata),
				  Double_val(reltol), atol_nv);
    CHECK_FLAG("ARKodeSVtolerances", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_reinit(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    ARKRhsFn fe = NULL;
    ARKRhsFn fi = NULL;

    switch (ARKODE_PROBLEM(vdata)) {
    case VARIANT_ARKODE_PROBLEM_TYPE_IMPLICIT_ONLY:
	fi = rhsfn1;
	break;
    case VARIANT_ARKODE_PROBLEM_TYPE_EXPLICIT_ONLY:
	fe = rhsfn2;
	break;
    case VARIANT_ARKODE_PROBLEM_TYPE_IMPLICIT_AND_EXPLICIT:
	fe = rhsfn2;
	fi = rhsfn1;
	break;
    }

    N_Vector y0_nv = NVEC_VAL(y0);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepReInit(ARKODE_MEM_FROM_ML(vdata),
			     fe, fi, Double_val(t0), y0_nv);
    CHECK_FLAG("ARKStepReInit", flag);
#else
    int flag = ARKodeReInit(ARKODE_MEM_FROM_ML(vdata),
			    fe, fi, Double_val(t0), y0_nv);
    CHECK_FLAG("ARKodeReInit", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_reset(value vdata, value vt, value vy)
{
    CAMLparam3(vdata, vt, vy);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepReset(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			    NVEC_VAL(vy));
    CHECK_FLAG("ARKStepReset", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

static value ark_solver(value vdata, value nextt, value vy, int onestep)
{
    CAMLparam3(vdata, nextt, vy);
    CAMLlocal1(ret);
    realtype tret;
    int flag;
    N_Vector y;
    enum arkode_solver_result_tag result = -1;
    const char* call;

    y = NVEC_VAL (vy);
    // Caml_ba_data_val(y) must not be shifted by the OCaml GC during this
    // function call, which calls Caml through the callback f.  Is this
    // guaranteed?
#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepEvolve(ARKODE_MEM_FROM_ML (vdata), Double_val (nextt),
			 y, &tret, onestep ? ARK_ONE_STEP : ARK_NORMAL);
    call = "ARKStepEvolve";
#else
    flag = ARKode (ARKODE_MEM_FROM_ML (vdata), Double_val (nextt),
		   y, &tret, onestep ? ARK_ONE_STEP : ARK_NORMAL);
    call = "ARKode";
#endif

    switch (flag) {
    case ARK_SUCCESS:
	result = VARIANT_ARKODE_SOLVER_RESULT_SUCCESS;
	break;

    case ARK_ROOT_RETURN:
	result = VARIANT_ARKODE_SOLVER_RESULT_ROOTSFOUND;
	break;

    case ARK_TSTOP_RETURN:
	result = VARIANT_ARKODE_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	sunml_arkode_check_flag(call, flag, ARKODE_MEM_FROM_ML(vdata));
    }

    assert (Field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}

CAMLprim value sunml_arkode_ark_solve_normal(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(ark_solver(vdata, nextt, y, 0));
}

CAMLprim value sunml_arkode_ark_solve_one_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(ark_solver(vdata, nextt, y, 1));
}

CAMLprim value sunml_arkode_ark_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVEC_VAL(vy);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetDky(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			     Int_val(vk), y_nv);
    CHECK_FLAG("ARKStepGetDky", flag);
#else
    int flag = ARKodeGetDky(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			    Int_val(vk), y_nv);
    CHECK_FLAG("ARKodeGetDky", flag);
#endif
    
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_err_weights(value varkode_mem, value verrws)
{
    CAMLparam2(varkode_mem, verrws);

    N_Vector errws_nv = NVEC_VAL(verrws);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetErrWeights(ARKODE_MEM_FROM_ML(varkode_mem), errws_nv);
    CHECK_FLAG("ARKStepGetErrWeights", flag);
#else
    int flag = ARKodeGetErrWeights(ARKODE_MEM_FROM_ML(varkode_mem), errws_nv);
    CHECK_FLAG("ARKodeGetErrWeights", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_res_weights(value varkode_mem, value vresws)
{
    CAMLparam2(varkode_mem, vresws);

#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector resws_nv = NVEC_VAL(vresws);
    int flag = ARKStepGetResWeights(ARKODE_MEM_FROM_ML(varkode_mem), resws_nv);
    CHECK_FLAG("ARKStepGetResWeights", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    N_Vector resws_nv = NVEC_VAL(vresws);
    int flag = ARKodeGetResWeights(ARKODE_MEM_FROM_ML(varkode_mem), resws_nv);
    CHECK_FLAG("ARKodeGetResWeights", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_est_local_errors(value varkode_mem,
						     value vele)
{
    CAMLparam2(varkode_mem, vele);

    N_Vector ele_nv = NVEC_VAL(vele);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetEstLocalErrors(ARKODE_MEM_FROM_ML(varkode_mem), ele_nv);
    CHECK_FLAG("ARKStepGetEstLocalErrors", flag);
#else
    int flag = ARKodeGetEstLocalErrors(ARKODE_MEM_FROM_ML(varkode_mem), ele_nv);
    CHECK_FLAG("ARKodeGetEstLocalErrors", flag);
#endif

    CAMLreturn (Val_unit);
}

#if 400 <= SUNDIALS_LIB_VERSION
void sunml_arkode_check_ls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == ARKLS_SUCCESS) return;

    switch (flag) {
	case ARKLS_ILL_INPUT:
	    caml_raise_constant(ARKODE_EXN(IllInput));

	case ARKLS_MEM_FAIL:
	    caml_raise_out_of_memory();

	case ARKLS_SUNMAT_FAIL:
	case ARKLS_SUNLS_FAIL:
	case ARKLS_JACFUNC_UNRECVR:
	case ARKLS_JACFUNC_RECVR:
	case ARKLS_MASSFUNC_UNRECVR:
	case ARKLS_MASSFUNC_RECVR:
	default:
	    /* e.g. ARKLS_MEM_NULL, ARKLS_LMEM_NULL */
	    // ARKStepGetLinReturnFlagName and MRIStepGetLinReturnFlagName
	    // share the same underlying implementation
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    ARKStepGetLinReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#else
void sunml_arkode_check_dls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == ARKDLS_SUCCESS) return;

    switch (flag) {
	case ARKDLS_ILL_INPUT:
	    caml_raise_constant(ARKODE_EXN(IllInput));

	case ARKDLS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if 300 <= SUNDIALS_LIB_VERSION
	case ARKDLS_SUNMAT_FAIL:
#endif
	case ARKDLS_JACFUNC_UNRECVR:
	case ARKDLS_JACFUNC_RECVR:
	case ARKDLS_MASSFUNC_UNRECVR:
	case ARKDLS_MASSFUNC_RECVR:
	default:
	    /* e.g. ARKDLS_MEM_NULL, ARKDLS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    ARKDlsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

void sunml_arkode_check_spils_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == ARKSPILS_SUCCESS) return;

    switch (flag) {
	case ARKSPILS_ILL_INPUT:
	    caml_raise_constant(ARKODE_EXN(IllInput));

	case ARKSPILS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if 300 <= SUNDIALS_LIB_VERSION
	case ARKSPILS_SUNLS_FAIL:
#endif
	default:
	    /* e.g. ARKSPILS_MEM_NULL, ARKSPILS_PMEM_NULL, ARKSPILS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     ARKSpilsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#endif

CAMLprim value sunml_arkode_ark_session_finalize(value vdata)
{
    if (ARKODE_MEM_FROM_ML(vdata) != NULL) {
	void *arkode_mem = ARKODE_MEM_FROM_ML(vdata);
	value *backref = ARKODE_BACKREF_FROM_ML(vdata);
#if 400 <= SUNDIALS_LIB_VERSION
	ARKStepFree(&arkode_mem);
#else
	ARKodeFree(&arkode_mem);
#endif
	sunml_sundials_free_value(backref);
    }

    return Val_unit;
}

CAMLprim value sunml_arkode_ark_ss_tolerances(value vdata, value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSStolerances(ARKODE_MEM_FROM_ML(vdata),
				   Double_val(reltol),
				   Double_val(abstol));
    CHECK_FLAG("ARKStepSStolerances", flag);
#else
    int flag = ARKodeSStolerances(ARKODE_MEM_FROM_ML(vdata),
				  Double_val(reltol),
				  Double_val(abstol));
    CHECK_FLAG("ARKodeSStolerances", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetRootInfo(ARKODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("ARKStepGetRootInfo", flag);
#else
    int flag = ARKodeGetRootInfo(ARKODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("ARKodeGetRootInfo", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_current_state(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vnv);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector nv;
    int flag = ARKStepGetCurrentState(ARKODE_MEM_FROM_ML(varkode_mem), &nv);
    CHECK_FLAG("ARKStepGetCurrentState", flag);

    vnv = NVEC_BACKLINK(nv);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim value sunml_arkode_ark_get_nonlin_system_data(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vnv);
#if 540 <= SUNDIALS_LIB_VERSION
    realtype tcur, gamma;
    N_Vector zpred, zi, Fi, sdata;
    void *user_data;

    int flag = ARKStepGetNonlinearSystemData(ARKODE_MEM_FROM_ML(varkode_mem),
		    &tcur, &zpred, &zi, &Fi, &gamma, &sdata, &user_data);
    CHECK_FLAG("ARKStepGetNonlinearSystemData", flag);

    vnv = caml_alloc_tuple(RECORD_ARKODE_NONLIN_SYSTEM_DATA_SIZE);
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_TCUR,
	    caml_copy_double(tcur));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_ZPRED,
	    NVEC_BACKLINK(zpred));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_ZI,
	    NVEC_BACKLINK(zi));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_FI,
	    NVEC_BACKLINK(Fi));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_GAMMA,
	    caml_copy_double(gamma));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_SDATA,
	    NVEC_BACKLINK(sdata));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim value sunml_arkode_ark_compute_state(value varkode_mem,
					      value vzcor, value vz)
{
    CAMLparam3(varkode_mem, vzcor, vz);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepComputeState(ARKODE_MEM_FROM_ML(varkode_mem),
				   NVEC_VAL(vzcor), NVEC_VAL(vz));
    CHECK_FLAG("ARKStepComputeState", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim value sunml_arkode_ark_get_current_gamma(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    double gamma;

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetCurrentGamma(ARKODE_MEM_FROM_ML(varkode_mem), &gamma);
    CHECK_FLAG("ARKStepGetCurrentGamma", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(gamma));
}

CAMLprim value sunml_arkode_ark_get_timestepper_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;

    long int expsteps;
    long int accsteps;
    long int step_attempts;
    long int nfe_evals;
    long int nfi_evals;
    long int nlinsetups;
    long int netfails;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetTimestepperStats(ARKODE_MEM_FROM_ML(vdata),
                                      &expsteps,
                                      &accsteps,
                                      &step_attempts,
                                      &nfe_evals,
                                      &nfi_evals,
                                      &nlinsetups,
                                      &netfails);
    CHECK_FLAG("ARKStepGetTimestepperStats", flag);
#else
    long int nsteps;
    realtype hinused;
    realtype hlast;
    realtype hcur;
    realtype tcur;

    flag = ARKodeGetIntegratorStats(ARKODE_MEM_FROM_ML(vdata),
				    &nsteps,
				    &expsteps,
				    &accsteps,
				    &step_attempts,
				    &nfe_evals,
				    &nfi_evals,
				    &nlinsetups,
				    &netfails,
				    &hinused,
				    &hlast,
				    &hcur,
				    &tcur);
    CHECK_FLAG("ARKodeGetIntegratorStats", flag);
#endif

    r = caml_alloc_tuple(RECORD_ARKODE_ARK_TIMESTEPPER_STATS_SIZE);
    Store_field(r, RECORD_ARKODE_ARK_TIMESTEPPER_STATS_EXP_STEPS,
		   Val_long(expsteps));
    Store_field(r, RECORD_ARKODE_ARK_TIMESTEPPER_STATS_ACC_STEPS,
		   Val_long(accsteps));
    Store_field(r, RECORD_ARKODE_ARK_TIMESTEPPER_STATS_STEP_ATTEMPTS,
		   Val_long(step_attempts));
    Store_field(r, RECORD_ARKODE_ARK_TIMESTEPPER_STATS_NUM_NFE_EVALS,
		   Val_long(nfe_evals));
    Store_field(r, RECORD_ARKODE_ARK_TIMESTEPPER_STATS_NUM_NFI_EVALS,
		   Val_long(nfi_evals));
    Store_field(r, RECORD_ARKODE_ARK_TIMESTEPPER_STATS_LIN_SOLV_SETUPS,
		   Val_long(nlinsetups));
    Store_field(r, RECORD_ARKODE_ARK_TIMESTEPPER_STATS_NUM_ERR_TEST_FAILS,
		   Val_long(netfails));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_get_step_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;

    long int nsteps;
    realtype hinused;
    realtype hlast;
    realtype hcur;
    realtype tcur;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetStepStats(ARKODE_MEM_FROM_ML(vdata),
                               &nsteps,
                               &hinused,
                               &hlast,
                               &hcur,
                               &tcur);
    CHECK_FLAG("ARKStepGetStepStats", flag);
#else
    long int expsteps;
    long int accsteps;
    long int step_attempts;
    long int nfe_evals;
    long int nfi_evals;
    long int nlinsetups;
    long int netfails;

    flag = ARKodeGetIntegratorStats(ARKODE_MEM_FROM_ML(vdata),
				    &nsteps,
				    &expsteps,
				    &accsteps,
				    &step_attempts,
				    &nfe_evals,
				    &nfi_evals,
				    &nlinsetups,
				    &netfails,
				    &hinused,
				    &hlast,
				    &hcur,
				    &tcur);
    CHECK_FLAG("ARKodeGetIntegratorStats", flag);
#endif

    r = caml_alloc_tuple(RECORD_ARKODE_STEP_STATS_SIZE);
    Store_field(r, RECORD_ARKODE_STEP_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_ARKODE_STEP_STATS_ACTUAL_INIT_STEP,
						    caml_copy_double(hinused));
    Store_field(r, RECORD_ARKODE_STEP_STATS_LAST_STEP,
						    caml_copy_double(hlast));
    Store_field(r, RECORD_ARKODE_STEP_STATS_CURRENT_STEP,
						    caml_copy_double(hcur));
    Store_field(r, RECORD_ARKODE_STEP_STATS_CURRENT_TIME,
						    caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_set_error_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetErrFile(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ARKStepSetErrFile", flag);
#else
    int flag = ARKodeSetErrFile(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ARKodeSetErrFile", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_diagnostics(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ARKStepSetDiagnostics", flag);
#else
    int flag = ARKodeSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ARKodeSetDiagnostics", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_clear_diagnostics(value vdata)
{
    CAMLparam1(vdata);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("ARKStepSetDiagnostics", flag);
#else
    int flag = ARKodeSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("ARKodeSetDiagnostics", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetRootDirection(ARKODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("ARKStepSetRootDirection", flag);
#else
    int flag = ARKodeSetRootDirection(ARKODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("ARKodeSetRootDirection", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_resize(value varkode_mem,
				       value vhasfn, value vhscale,
				       value vt0, value vynew)
{
    CAMLparam5(varkode_mem, vhasfn, vhscale, vt0, vynew);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepResize(
		 arkode_mem,
		 NVEC_VAL(vynew),
		 Double_val(vhscale),
		 Double_val(vt0),
		 Bool_val(vhasfn) ? resizefn : NULL,
		 Bool_val(vhasfn) ? ARKODE_BACKREF_FROM_ML(varkode_mem) : NULL);
    CHECK_FLAG("ARKStepResize", flag);
#else
    int flag = ARKodeResize(
		 arkode_mem,
		 NVEC_VAL(vynew),
		 Double_val(vhscale),
		 Double_val(vt0),
		 Bool_val(vhasfn) ? resizefn : NULL,
		 Bool_val(vhasfn) ? ARKODE_BACKREF_FROM_ML(varkode_mem) : NULL);
    CHECK_FLAG("ARKodeResize", flag);
#endif

    CAMLreturn (Val_unit);
}

#if 400 <= SUNDIALS_LIB_VERSION
static value val_butcher_table(ARKodeButcherTable bt)
{
    CAMLparam0();
    CAMLlocal5(vbt, va, vc, vb, vd);
    CAMLlocal2(vobt, vod);
    int i, j;
    intnat q, p, stages;
    realtype **A, *c, *b, *d;

    if (bt == NULL) {
	vobt = Val_none;
    } else {
	q = bt->q;
	p = bt->p;
	stages = bt->stages;

	// Allocate OCaml arrays
	va = sunml_sundials_realarray2_create(stages, stages);
	A = ARRAY2_ACOLS(va);
	vc = caml_ba_alloc(BIGARRAY_FLOAT, 1, NULL, &stages);
	c = REAL_ARRAY(vc);
	vb = caml_ba_alloc(BIGARRAY_FLOAT, 1, NULL, &stages);
	b = REAL_ARRAY(vb);
	vod = Val_none;

	// Copy contents
	for (i = 0; i < stages; i++) {
	    c[i] = bt->c[i];
	    b[i] = bt->b[i];
	    for (j = 0; j < stages; j++) {
		A[i][j] = bt->A[i][j];
	    }
	}

	// Allocate and copy embedding if necessary
	if (bt->d != NULL) {
	    vd = caml_ba_alloc(BIGARRAY_FLOAT, 1, NULL, &stages);
	    d = REAL_ARRAY(vd);
	    Store_some(vod, vd);

	    for (i = 0; i < stages; i++) {
		d[i] = bt->d[i];
	    }
	}

	vbt = caml_alloc_tuple(RECORD_ARKODE_BUTCHER_TABLE_SIZE);
	Store_field(vbt, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER, Val_int(q));
	Store_field(vbt, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER, Val_int(p));
	Store_field(vbt, RECORD_ARKODE_BUTCHER_TABLE_STAGES, Val_int(stages));
	Store_field(vbt, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES, va);
	Store_field(vbt, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES, vc);
	Store_field(vbt, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS, vb);
	Store_field(vbt, RECORD_ARKODE_BUTCHER_TABLE_BEMBED, vod);
	Store_some(vobt, vbt);
    }
    CAMLreturn(vobt);
}
#endif

CAMLprim value sunml_arkode_ark_get_current_butcher_tables(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal3(vobti, vobte, vr);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag;
    ARKodeButcherTable bti, bte;
    flag = ARKStepGetCurrentButcherTables(ARKODE_MEM_FROM_ML(varkode_mem),
					  &bti, &bte);
    CHECK_FLAG("ARKStepGetCurrentButcherTables", flag);

    vobti = val_butcher_table(bti);
    vobte = val_butcher_table(bte);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, vobti);
    Store_field(vr, 1, vobte);
#elif 270 <= SUNDIALS_LIB_VERSION
    CAMLlocal2(vbti, vbte);
    CAMLlocal5(vai, vci, vbi, vdi, vodi);
    CAMLlocal5(vae, vce, vbe, vde, vode);
    int s, q, p, i, j;
    realtype *Ai, *Ae, *ai, *ae;

    vai = sunml_sundials_realarray2_create(ARK_S_MAX, ARK_S_MAX);
    vae = sunml_sundials_realarray2_create(ARK_S_MAX, ARK_S_MAX);
    vci = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    vce = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    vbi = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    vbe = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    vdi = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    Store_some(vodi, vdi);
    vde = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    Store_some(vode, vde);

    Ai = calloc(ARK_S_MAX * ARK_S_MAX, sizeof(realtype));
    if (Ai == NULL) caml_raise_out_of_memory();
    Ae = calloc(ARK_S_MAX * ARK_S_MAX, sizeof(realtype));
    if (Ae ==NULL) {
	free(Ai);
	caml_raise_out_of_memory();
    }

    ARKodeGetCurrentButcherTables(ARKODE_MEM_FROM_ML(varkode_mem),
				  &s, &q, &p,
				  Ai, Ae,
				  REAL_ARRAY(vci), REAL_ARRAY(vce),
				  REAL_ARRAY(vbi), REAL_ARRAY(vbe),
				  REAL_ARRAY(vdi), REAL_ARRAY(vde));
    // No need to check flag: can only fail if arkode_mem == NULL
    // (must not risk an exception before freeing Ai and Ae)

    // NB: translate to column-major from row-major
    ai = ARRAY2_DATA(vai);
    ae = ARRAY2_DATA(vae);
    for (i=0; i < ARK_S_MAX; i++) {
	for (j=0; j < ARK_S_MAX; j++) {
	    ARK_A(ai,j,i) = ARK_A(Ai,i,j);
	    ARK_A(ae,j,i) = ARK_A(Ae,i,j);
	}
    }
    free(Ai);
    free(Ae);

    vbti = caml_alloc_tuple(RECORD_ARKODE_BUTCHER_TABLE_SIZE);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER, Val_int(q));
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER, Val_int(p));
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_STAGES, Val_int(s));
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES, vai);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES, vci);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS, vbi);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_BEMBED, vodi);
    Store_some(vobti, vbti);

    vbte = caml_alloc_tuple(RECORD_ARKODE_BUTCHER_TABLE_SIZE);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER, Val_int(q));
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER, Val_int(p));
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_STAGES, Val_int(s));
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES, vae);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES, vce);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS, vbe);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_BEMBED, vode);
    Store_some(vobte, vbte);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, vobti);
    Store_field(vr, 1, vobte);
#else // SUNDIALS_LIB_VERSION < 270
    CAMLlocal2(vbti, vbte);
    CAMLlocal5(vai, vci, vbi, vdi, vodi);
    CAMLlocal1(vae);
    int s, q, p, i, j;
    realtype *Ai, *Ae, *ai, *ae;

    vai = sunml_sundials_realarray2_create(ARK_S_MAX, ARK_S_MAX);
    vae = sunml_sundials_realarray2_create(ARK_S_MAX, ARK_S_MAX);
    vci = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    vbi = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    vdi = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    Store_some(vodi, vdi);

    Ai = calloc(ARK_S_MAX * ARK_S_MAX, sizeof(realtype));
    if (Ai == NULL) caml_raise_out_of_memory();
    Ae = calloc(ARK_S_MAX * ARK_S_MAX, sizeof(realtype));
    if (Ae ==NULL) {
	free(Ai);
	caml_raise_out_of_memory();
    }

    ARKodeGetCurrentButcherTables(ARKODE_MEM_FROM_ML(varkode_mem),
				  &s, &q, &p,
				  Ai, Ae,
				  REAL_ARRAY(vci),
				  REAL_ARRAY(vbi),
				  REAL_ARRAY(vdi));
    // No need to check return flag: can only fail if arkode_mem == NULL
    // (must not risk an exception before freeing Ai and Ae)

    // NB: translate to column-major from row-major
    ai = ARRAY2_ACOLS(vai);
    ae = ARRAY2_ACOLS(vae);
    for (i=0; i < ARK_S_MAX; i++) {
	for (j=0; j < ARK_S_MAX; j++) {
	    ARK_A(ai,j,i) = ARK_A(Ai,i,j);
	    ARK_A(ae,j,i) = ARK_A(Ae,i,j);
	}
    }
    free(Ai);
    free(Ae);

    vbti = caml_alloc_tuple(RECORD_ARKODE_BUTCHER_TABLE_SIZE);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER, Val_int(q));
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER, Val_int(p));
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_STAGES, Val_int(stages));
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES, vai);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES, vci);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS, vbi);
    Store_field(vbti, RECORD_ARKODE_BUTCHER_TABLE_BEMBED, vodi);
    Store_some(vobti, vbti);

    vbte = caml_alloc_tuple(RECORD_ARKODE_BUTCHER_TABLE_SIZE);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER, Val_int(q));
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER, Val_int(p));
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_STAGES, Val_int(stages));
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES, vai);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES, vci);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS, vbi);
    Store_field(vbte, RECORD_ARKODE_BUTCHER_TABLE_BEMBED, vodi);
    Store_some(vobte, vbte);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, vobti);
    Store_field(vr, 1, vobte);
#endif

    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_butcher_table_load_erk(value vmethod)
{
    CAMLparam1(vmethod);
    CAMLlocal2(vobt, vbt);

#if 400 <= SUNDIALS_LIB_VERSION
    ARKodeButcherTable bt;
    bt = ARKodeButcherTable_LoadERK(Int_val(vmethod));

    vobt = val_butcher_table(bt);
    if (vobt == Val_none) caml_failwith("Failed to load ERK table");

    ARKodeButcherTable_Free(bt);
    vbt = Some_val(vobt);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vbt);
}

CAMLprim value sunml_arkode_butcher_table_load_dirk(value vmethod)
{
    CAMLparam1(vmethod);
    CAMLlocal2(vobt, vbt);

#if 400 <= SUNDIALS_LIB_VERSION
    ARKodeButcherTable bt;
    bt = ARKodeButcherTable_LoadDIRK(Int_val(vmethod));

    vobt = val_butcher_table(bt);
    if (vobt == Val_none) caml_failwith("Failed to load DIRK table");

    ARKodeButcherTable_Free(bt);
    vbt = Some_val(vobt);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vbt);
}

#if 400 <= SUNDIALS_LIB_VERSION
static ARKodeButcherTable butcher_table_val(value vob)
{
    CAMLparam0();
    ARKodeButcherTable r = NULL;
    CAMLlocal2(vb, vd);

    if (vob != Val_none) {
	vb = Some_val(vob);
	vd = Field(vb, RECORD_ARKODE_BUTCHER_TABLE_BEMBED);
	r = ARKodeButcherTable_Create(
		Int_val(Field(vb, RECORD_ARKODE_BUTCHER_TABLE_STAGES)),
		(vd == Val_none) ? 0
		  : Int_val(Field(vb, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER)),
		Int_val(Field(vb, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER)),
		REAL_ARRAY(Field(vb, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES)),
		ARRAY2_DATA(Field(vb, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES)),
		REAL_ARRAY(Field(vb, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS)),
		(vd == Val_none) ? NULL : REAL_ARRAY(Some_val(vd)));
    }

    CAMLreturnT(ARKodeButcherTable, r);
}
#endif

CAMLprim value sunml_arkode_butcher_table_write(value vbt, value volog)
{
    CAMLparam2(vbt, volog);
#if 400 <= SUNDIALS_LIB_VERSION
    FILE *vlog = (volog == Val_none) ? NULL : ML_CFILE(Some_val(volog));

    ARKodeButcherTable bt = butcher_table_val(vbt);
    ARKodeButcherTable_Write(bt, vlog);
    ARKodeButcherTable_Free(bt);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_butcher_table_check_order(value volog, value vbt)
{
    CAMLparam2(volog, vbt);
    CAMLlocal2(vr, vop);
#if 400 <= SUNDIALS_LIB_VERSION
    FILE *log = (volog == Val_none) ? NULL : ML_CFILE(Some_val(volog));
    int flag, q, p;

    ARKodeButcherTable bt = butcher_table_val(vbt);
    flag = ARKodeButcherTable_CheckOrder(bt, &q, &p, log);
    ARKodeButcherTable_Free(bt);

    if (p == 0) {
	vop = Val_none;
    } else {
	Store_some(vop, Val_int(p));
    }

    vr = caml_alloc_tuple(4);
    Store_field(vr, 0, Val_int(q));
    Store_field(vr, 1, vop);
    Store_field(vr, 2, Val_bool(flag < 0));
    Store_field(vr, 3, Val_bool(flag > 0));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_butcher_table_check_ark_order(value volog,
							  value vbt1,
							  value vbt2)
{
    CAMLparam3(volog, vbt1, vbt2);
    CAMLlocal2(vr, vop);
#if 400 <= SUNDIALS_LIB_VERSION
    FILE *log = (volog == Val_none) ? NULL : ML_CFILE(Some_val(volog));
    int flag, q, p;

    ARKodeButcherTable bt1 = butcher_table_val(vbt1);
    ARKodeButcherTable bt2 = butcher_table_val(vbt2);
    flag = ARKodeButcherTable_CheckARKOrder(bt1, bt2, &q, &p, log);
    ARKodeButcherTable_Free(bt2);
    ARKodeButcherTable_Free(bt1);

    if (flag < 0)
	caml_failwith("Unexpected error in ARKodeButcherTable_CheckARKOrder");

    if (p == 0) {
	vop = Val_none;
    } else {
	Store_some(vop, Val_int(p));
    }

    vr = caml_alloc_tuple(3);
    Store_field(vr, 0, Val_int(q));
    Store_field(vr, 1, vop);
    Store_field(vr, 3, Val_bool(flag > 0));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_ark_set_tables(value varkode_mem,
					   value vq, value vp,
					   value vobi, value vobe)
{
    CAMLparam5(varkode_mem, vq, vp, vobi, vobe);
    int flag;

#if   400 <= SUNDIALS_LIB_VERSION
    ARKodeButcherTable bi = butcher_table_val(vobi);
    ARKodeButcherTable be = butcher_table_val(vobe);

    flag = ARKStepSetTables(ARKODE_MEM_FROM_ML(varkode_mem),
			    Int_val(vq), Int_val(vp), bi, be);
    CHECK_FLAG("ARKStepSetTables", flag);

    if (bi != NULL) ARKodeButcherTable_Free(bi);
    if (be != NULL) ARKodeButcherTable_Free(be);

#else // SUNDIALS_LIB_VERSION < 400
    CAMLlocal4(vbi, vbe, vbembedi, vbembede);
    realtype *Ai = NULL, *Ae = NULL;
    realtype *ai, *ae;
    int s, p, q, i, j;

    if ((vobi != Val_none) && (vobe != Val_none)) {
	vbi = Some_val(vobi);
	vbe = Some_val(vobe);
	vbembedi = Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_BEMBED);
	vbembede = Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_BEMBED);

	s = Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGES));
	q = Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER));
	p = Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER));

	// NB: translate to row-major from column-major
	Ai = calloc(s * s, sizeof(realtype));
	if (Ai == NULL) caml_raise_out_of_memory();
	Ae = calloc(s * s, sizeof(realtype));
	if (Ae ==NULL) {
	    free(Ai);
	    caml_raise_out_of_memory();
	}
	ai = ARRAY2_DATA(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES));
	ae = ARRAY2_DATA(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES));
	for (i=0; i < s; i++) {
	    for (j=0; j < s; j++) {
		ARK_A(Ai,j,i) = ARK_A(ai,i,j);
		ARK_A(Ae,j,i) = ARK_A(ae,i,j);
	    }
	}

#if 270 <= SUNDIALS_LIB_VERSION
	flag = ARKodeSetARKTables(ARKODE_MEM_FROM_ML(varkode_mem),
		s, q, p,
		REAL_ARRAY(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES)),
		REAL_ARRAY(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES)),
		Ai, Ae,
		REAL_ARRAY(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS)),
		REAL_ARRAY(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS)),
		((vbembedi == Val_none) ? NULL : REAL_ARRAY(Some_val(vbembedi))),
		((vbembede == Val_none) ? NULL : REAL_ARRAY(Some_val(vbembede))));
#else // SUNDIALS_LIB_VERSION < 270
	flag = ARKodeSetARKTables(ARKODE_MEM_FROM_ML(varkode_mem),
		s, q, p
		REAL_ARRAY(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES)),
		Ai, Ae,
		REAL_ARRAY(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS)),
		((vbembedi == Val_none) ? NULL : REAL_ARRAY(Some_val(vbembedi))));
#endif
	free(Ai);
	free(Ae);
	CHECK_FLAG("ARKodeSetARKTables", flag);

    } else if ((vobi != Val_none) && (vobe == Val_none)) {
	vbi= Some_val(vobi);
	vbembedi = Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_BEMBED);
	s = Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGES));

	// NB: translate to row-major from column-major
	Ai = calloc(s * s, sizeof(realtype));
	if (Ai == NULL) caml_raise_out_of_memory();
	ai = ARRAY2_DATA(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES));
	for (i=0; i < s; i++) {
	    for (j=0; j < s; j++) {
		ARK_A(Ai,j,i) = ARK_A(ai,i,j);
	    }
	}

	flag = ARKodeSetIRKTable(ARKODE_MEM_FROM_ML(varkode_mem),
	    Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGES)),
	    Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER)),
	    Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER)),
	    REAL_ARRAY(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES)),
	    Ai,
	    REAL_ARRAY(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS)),
	    ((vbembedi == Val_none) ? NULL : REAL_ARRAY(Some_val(vbembedi))));
	free(Ai);
	CHECK_FLAG("ARKodeSetIRKTable", flag);

    } else if ((vobi == Val_none) && (vobe != Val_none)) {
	vbe = Some_val(vobe);
	vbembede = Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_BEMBED);
	s = Int_val(Field(vbi, RECORD_ARKODE_BUTCHER_TABLE_STAGES));

	// NB: translate to row-major from column-major
	Ae = calloc(s * s, sizeof(realtype));
	if (Ae == NULL) caml_raise_out_of_memory();
	ae = ARRAY2_DATA(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES));
	for (i=0; i < s; i++) {
	    for (j=0; j < s; j++) {
		ARK_A(Ae,j,i) = ARK_A(ae,i,j);
	    }
	}

	flag = ARKodeSetERKTable(ARKODE_MEM_FROM_ML(varkode_mem),
	    Int_val(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_STAGES)),
	    Int_val(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER)),
	    Int_val(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING_ORDER)),
	    REAL_ARRAY(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES)),
	    Ae,
	    REAL_ARRAY(Field(vbe, RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS)),
	    ((vbembede == Val_none) ? NULL : REAL_ARRAY(Some_val(vbembede))));
	free(Ae);
	CHECK_FLAG("ARKodeSetERKTable", flag);
    }

#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_table_num(value varkode_mem, value vnums)
{
    CAMLparam2(varkode_mem, vnums);
    int itable = Int_val(Field(vnums, 0));
    int etable = Int_val(Field(vnums, 1));
    int flag;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepSetTableNum(ARKODE_MEM_FROM_ML(varkode_mem), itable, etable);
    CHECK_FLAG("ARKStepSetTableNum", flag);
#else
    if (itable < 0) {
	flag = ARKodeSetERKTableNum(ARKODE_MEM_FROM_ML(varkode_mem), etable);
	CHECK_FLAG("ARKodeSetERKTableNum", flag);
    } else if (etable < 0) {
	flag = ARKodeSetIRKTableNum(ARKODE_MEM_FROM_ML(varkode_mem), itable);
	CHECK_FLAG("ARKodeSetIRKTableNum", flag);
    } else {
	flag = ARKodeSetARKTableNum(ARKODE_MEM_FROM_ML(varkode_mem),
				    itable, etable);
	CHECK_FLAG("ARKodeSetARKTableNum", flag);
    }
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_adaptivity_method(value varkode_mem, value vmeth)
{
    CAMLparam2(varkode_mem, vmeth);
    CAMLlocal2(vks, vorder);

    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    if (Tag_val(vmeth) == VARIANT_ARKODE_ADAPTIVITY_METHOD_ADAPTIVITYFN) {
#if 400 <= SUNDIALS_LIB_VERSION
	flag = ARKStepSetAdaptivityFn(arkode_mem,
				      adaptfn,
				      ARKODE_BACKREF_FROM_ML(varkode_mem));
	CHECK_FLAG("ARKStepSetAdaptivityFn", flag);
#else
	flag = ARKodeSetAdaptivityFn(arkode_mem,
				     adaptfn,
				     ARKODE_BACKREF_FROM_ML(varkode_mem));
	CHECK_FLAG("ARKodeSetAdaptivityFn", flag);
#endif

    } else {
	realtype adapt_params[3] = { 0 };

	vks = Field(Field(vmeth, 0), RECORD_ARKODE_ADAPTIVITY_PARAMS_KS);
	vorder = Field(Field(vmeth, 0),
			RECORD_ARKODE_ADAPTIVITY_PARAMS_METHOD_ORDER);

	if (vks != Val_none) {
	    adapt_params[0] = Double_val(Field(Some_val(vks), 0));
	    adapt_params[1] = Double_val(Field(Some_val(vks), 1));
	    adapt_params[2] = Double_val(Field(Some_val(vks), 2));
	}

#if 400 <= SUNDIALS_LIB_VERSION
	flag = ARKStepSetAdaptivityMethod(arkode_mem, 
					  Tag_val(vmeth), 
					  vks == Val_none,
					  Bool_val(vorder),
					  vks != Val_none ? adapt_params : NULL);
	CHECK_FLAG("ARKStepSetAdaptivityMethod", flag);
#else
	flag = ARKodeSetAdaptivityMethod(arkode_mem, 
					 Tag_val(vmeth), 
					 vks == Val_none,
					 Bool_val(vorder),
					 vks != Val_none ? adapt_params : NULL);
	CHECK_FLAG("ARKodeSetAdaptivityMethod", flag);
#endif
    }

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_stability_fn(value varkode_mem, value vhasf)
{
    CAMLparam2(varkode_mem, vhasf);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetStabilityFn(arkode_mem,
				     Bool_val(vhasf) ? stabfn : NULL,
				     ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetStabilityFn", flag);
#else
    int flag = ARKodeSetStabilityFn(arkode_mem,
				    Bool_val(vhasf) ? stabfn : NULL,
				    ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetStabilityFn", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_imex(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetImEx(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetImEx", flag);
#else
    int flag = ARKodeSetImEx(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetImEx", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_implicit(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetImplicit(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetImplicit", flag);
#else
    int flag = ARKodeSetImplicit(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetImplicit", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_explicit(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetExplicit(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetExplicit", flag);
#else
    int flag = ARKodeSetExplicit(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetExplicit", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_fixed_point(value varkode_mem, value vfpm)
{
    CAMLparam2(varkode_mem, vfpm);
#if SUNDIALS_LIB_VERSION < 400
    int flag = ARKodeSetFixedPoint(ARKODE_MEM_FROM_ML(varkode_mem),
				   Long_val(vfpm));
    CHECK_FLAG("ARKodeSetFixedPoint", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_newton(value varkode_mem)
{
    CAMLparam1(varkode_mem);
#if SUNDIALS_LIB_VERSION < 400
    int flag = ARKodeSetNewton(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetNewton", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_linear(value varkode_mem, value vtimedepend)
{
    CAMLparam2(varkode_mem, vtimedepend);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetLinear(ARKODE_MEM_FROM_ML(varkode_mem),
			        Bool_val(vtimedepend));
    CHECK_FLAG("ARKStepSetLinear", flag);
#else
    int flag = ARKodeSetLinear(ARKODE_MEM_FROM_ML(varkode_mem),
			       Bool_val(vtimedepend));
    CHECK_FLAG("ARKodeSetLinear", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_nonlinear(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetNonlinear(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetNonlinear", flag);
#else
    int flag = ARKodeSetNonlinear(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetNonlinear", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_predictor_method(value varkode_mem, value vmethod)
{
    CAMLparam2(varkode_mem, vmethod);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetPredictorMethod(ARKODE_MEM_FROM_ML(varkode_mem),
				         Int_val(vmethod));
    CHECK_FLAG("ARKStepSetPredictorMethod", flag);
#else
    int method = Int_val(vmethod);
    if (method == 5)
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    int flag = ARKodeSetPredictorMethod(ARKODE_MEM_FROM_ML(varkode_mem),
					method);
    CHECK_FLAG("ARKodeSetPredictorMethod", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_postprocess_step_fn(value varkode_mem,
							value vhasf)
{
    CAMLparam2(varkode_mem, vhasf);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKStepSetPostprocessStepFn(arkode_mem,
					   Bool_val(vhasf) ? poststepfn : NULL);
    CHECK_FLAG("ARKStepSetPostprocessStepFn", flag);
#elif 270 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKodeSetPostprocessStepFn(arkode_mem,
					  Bool_val(vhasf) ? poststepfn : NULL);
    CHECK_FLAG("ARKodeSetPostprocessStepFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Boiler plate definitions for ARKStep interface.
 */

CAMLprim value sunml_arkode_ark_resv_tolerance(value vdata, value abstol)
{
    CAMLparam2(vdata, abstol);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepResVtolerance(ARKODE_MEM_FROM_ML(vdata), NVEC_VAL(abstol));
    CHECK_FLAG("ARKStepResVtolerance", flag);
#else
    int flag = ARKodeResVtolerance(ARKODE_MEM_FROM_ML(vdata), NVEC_VAL(abstol));
    CHECK_FLAG("ARKodeResVtolerance", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_ress_tolerance(value vdata, value abstol)
{
    CAMLparam2(vdata, abstol);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepResStolerance(ARKODE_MEM_FROM_ML(vdata),
				    Double_val(abstol));
    CHECK_FLAG("ARKStepResStolerance", flag);
#else
    int flag = ARKodeResStolerance(ARKODE_MEM_FROM_ML(vdata),
				   Double_val(abstol));
    CHECK_FLAG("ARKodeResStolerance", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_resf_tolerance(value vdata)
{
    CAMLparam1(vdata);
 
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepResFtolerance(ARKODE_MEM_FROM_ML(vdata), resw);
    CHECK_FLAG("ARKStepResFtolerance", flag);
#else
    int flag = ARKodeResFtolerance(ARKODE_MEM_FROM_ML(vdata), resw);
    CHECK_FLAG("ARKodeResFtolerance", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem), &lenrw, &leniw);
    CHECK_FLAG("ARKStepGetWorkSpace", flag);
#else
    flag = ARKodeGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem), &lenrw, &leniw);
    CHECK_FLAG("ARKodeGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_get_num_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetNumSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetNumSteps", flag);
#else
    flag = ARKodeGetNumSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumSteps", flag);
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_ark_get_num_acc_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetNumAccSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetNumAccSteps", flag);
#else
    flag = ARKodeGetNumAccSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumAccSteps", flag);
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_ark_get_num_exp_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetNumExpSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetNumExpSteps", flag);
#else
    flag = ARKodeGetNumExpSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumExpSteps", flag);
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_ark_get_num_step_attempts(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetNumStepAttempts(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetNumStepAttempts", flag);
#else
    flag = ARKodeGetNumStepAttempts(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumStepAttempts", flag);
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_ark_get_num_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int nfe;
    long int nfi;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &nfe, &nfi);
    CHECK_FLAG("ARKStepGetNumRhsEvals", flag);
#else
    flag = ARKodeGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &nfe, &nfi);
    CHECK_FLAG("ARKodeGetNumRhsEvals", flag);
#endif

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nfe));
    Store_field(r, 1, Val_long(nfi));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_get_num_lin_solv_setups(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetNumLinSolvSetups(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumLinSolvSetups", flag);
#else
    flag = ARKodeGetNumLinSolvSetups(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumLinSolvSetups", flag);
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_ark_get_num_err_test_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetNumErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetNumErrTestFails", flag);
#else
    flag = ARKodeGetNumErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumErrTestFails", flag);
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_ark_get_actual_init_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetActualInitStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetActualInitStep", flag);
#else
    flag = ARKodeGetActualInitStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetActualInitStep", flag);
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_ark_get_last_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetLastStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetLastStep", flag);
#else
    flag = ARKodeGetLastStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetLastStep", flag);
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_ark_get_current_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetCurrentStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetCurrentStep", flag);
#else
    flag = ARKodeGetCurrentStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetCurrentStep", flag);
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_ark_get_current_time(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetCurrentTime(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKStepGetCurrentTime", flag);
#else
    flag = ARKodeGetCurrentTime(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetCurrentTime", flag);
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_ark_set_max_num_steps(value varkode_mem, value mxsteps)
{
    CAMLparam2(varkode_mem, mxsteps);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxNumSteps(ARKODE_MEM_FROM_ML(varkode_mem),
				     Long_val(mxsteps));
    CHECK_FLAG("ARKStepSetMaxNumSteps", flag);
#else
    int flag = ARKodeSetMaxNumSteps(ARKODE_MEM_FROM_ML(varkode_mem),
				    Long_val(mxsteps));
    CHECK_FLAG("ARKodeSetMaxNumSteps", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_hnil_warns(value varkode_mem, value mxhnil)
{
    CAMLparam2(varkode_mem, mxhnil);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxHnilWarns(ARKODE_MEM_FROM_ML(varkode_mem),
				      Int_val(mxhnil));
    CHECK_FLAG("ARKStepSetMaxHnilWarns", flag);
#else
    int flag = ARKodeSetMaxHnilWarns(ARKODE_MEM_FROM_ML(varkode_mem),
				     Int_val(mxhnil));
    CHECK_FLAG("ARKodeSetMaxHnilWarns", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_init_step(value varkode_mem, value hin)
{
    CAMLparam2(varkode_mem, hin);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetInitStep(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(hin));
    CHECK_FLAG("ARKStepSetInitStep", flag);
#else
    int flag = ARKodeSetInitStep(ARKODE_MEM_FROM_ML(varkode_mem),
				 Double_val(hin));
    CHECK_FLAG("ARKodeSetInitStep", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_min_step(value varkode_mem, value hmin)
{
    CAMLparam2(varkode_mem, hmin);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMinStep(ARKODE_MEM_FROM_ML(varkode_mem),
		 		 Double_val(hmin));
    CHECK_FLAG("ARKStepSetMinStep", flag);
#else
    int flag = ARKodeSetMinStep(ARKODE_MEM_FROM_ML(varkode_mem),
				Double_val(hmin));
    CHECK_FLAG("ARKodeSetMinStep", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_step(value varkode_mem, value hmax)
{
    CAMLparam2(varkode_mem, hmax);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxStep(ARKODE_MEM_FROM_ML(varkode_mem),
		 		 Double_val(hmax));
    CHECK_FLAG("ARKStepSetMaxStep", flag);
#else
    int flag = ARKodeSetMaxStep(ARKODE_MEM_FROM_ML(varkode_mem),
				Double_val(hmax));
    CHECK_FLAG("ARKodeSetMaxStep", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_stop_time(value varkode_mem, value tstop)
{
    CAMLparam2(varkode_mem, tstop);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetStopTime(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(tstop));
    CHECK_FLAG("ARKStepSetStopTime", flag);
#else
    int flag = ARKodeSetStopTime(ARKODE_MEM_FROM_ML(varkode_mem),
				 Double_val(tstop));
    CHECK_FLAG("ARKodeSetStopTime", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_err_test_fails(value varkode_mem, value maxnef)
{
    CAMLparam2(varkode_mem, maxnef);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem),
					 Int_val(maxnef));
    CHECK_FLAG("ARKStepSetMaxErrTestFails", flag);
#else
    int flag = ARKodeSetMaxErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem),
					Int_val(maxnef));
    CHECK_FLAG("ARKodeSetMaxErrTestFails", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_nonlin_iters(value varkode_mem, value maxcor)
{
    CAMLparam2(varkode_mem, maxcor);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxNonlinIters(ARKODE_MEM_FROM_ML(varkode_mem),
				        Int_val(maxcor));
    CHECK_FLAG("ARKStepSetMaxNonlinIters", flag);
#else
    int flag = ARKodeSetMaxNonlinIters(ARKODE_MEM_FROM_ML(varkode_mem),
				       Int_val(maxcor));
    CHECK_FLAG("ARKodeSetMaxNonlinIters", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_conv_fails(value varkode_mem, value maxncf)
{
    CAMLparam2(varkode_mem, maxncf);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxConvFails(ARKODE_MEM_FROM_ML(varkode_mem),
				      Int_val(maxncf));
    CHECK_FLAG("ARKStepSetMaxConvFails", flag);
#else
    int flag = ARKodeSetMaxConvFails(ARKODE_MEM_FROM_ML(varkode_mem),
				     Int_val(maxncf));
    CHECK_FLAG("ARKodeSetMaxConvFails", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_stage_predict_fn(value varkode_mem,
						     value vset)
{
    CAMLparam2(varkode_mem, vset);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetStagePredictFn(ARKODE_MEM_FROM_ML(varkode_mem),
				    Bool_val(vset) ? stagepredictfn : NULL);
    CHECK_FLAG("ARKStepSetStagePredictFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_nonlin_conv_coef(value varkode_mem, value nlscoef)
{
    CAMLparam2(varkode_mem, nlscoef);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetNonlinConvCoef(ARKODE_MEM_FROM_ML(varkode_mem),
				        Double_val(nlscoef));
    CHECK_FLAG("ARKStepSetNonlinConvCoef", flag);
#else
    int flag = ARKodeSetNonlinConvCoef(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(nlscoef));
    CHECK_FLAG("ARKodeSetNonlinConvCoef", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_constraints(value varkode_mem, value vnv)
{
    CAMLparam2(varkode_mem, vnv);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetConstraints(ARKODE_MEM_FROM_ML(varkode_mem),
				     NVEC_VAL(vnv));
    CHECK_FLAG("ARKStepSetConstraints", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_nls_rhs_fn(value varkode_mem)
{
    CAMLparam1(varkode_mem);
#if 580 <= SUNDIALS_LIB_VERSION
    int flag;

    flag = ARKStepSetNlsRhsFn(ARKODE_MEM_FROM_ML (varkode_mem), nlsrhsfn);
    CHECK_FLAG ("ARKStepSetNlsRhsFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_no_inactive_root_warn(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetNoInactiveRootWarn(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetNoInactiveRootWarn", flag);
#else
    int flag = ARKodeSetNoInactiveRootWarn(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetNoInactiveRootWarn", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_cfl_fraction(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetCFLFraction(ARKODE_MEM_FROM_ML(varkode_mem),
				    Double_val(varg));
    CHECK_FLAG("ARKStepSetCFLFraction", flag);
#else
    int flag = ARKodeSetCFLFraction(ARKODE_MEM_FROM_ML(varkode_mem),
				    Double_val(varg));
    CHECK_FLAG("ARKodeSetCFLFraction", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_defaults(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetDefaults(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetDefaults", flag);
#else
    int flag = ARKodeSetDefaults(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetDefaults", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_delta_gamma_max(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetDeltaGammaMax(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ARKStepSetDeltaGammaMax", flag);
#else
    int flag = ARKodeSetDeltaGammaMax(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ARKodeSetDeltaGammaMax", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_error_bias(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetErrorBias(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKStepSetErrorBias", flag);
#else
    int flag = ARKodeSetErrorBias(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKodeSetErrorBias", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_fixed_step(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetFixedStep(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKStepSetFixedStep", flag);
#else
    int flag = ARKodeSetFixedStep(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKodeSetFixedStep", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_num_constr_fails(value varkode_mem, 
							 value vmaxfails)
{
    CAMLparam2(varkode_mem, vmaxfails);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxNumConstrFails(ARKODE_MEM_FROM_ML(varkode_mem),
					   Int_val(vmaxfails));
    CHECK_FLAG("ARKStepSetMaxNumConstrFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_fixed_step_bounds(value varkode_mem,
						      value vlb, value vub)
{
    CAMLparam3(varkode_mem, vlb, vub);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetFixedStepBounds(ARKODE_MEM_FROM_ML(varkode_mem),
				        Double_val(vlb), Double_val(vub));
    CHECK_FLAG("ARKStepSetFixedStepBounds", flag);
#else
    int flag = ARKodeSetFixedStepBounds(ARKODE_MEM_FROM_ML(varkode_mem),
				        Double_val(vlb), Double_val(vub));
    CHECK_FLAG("ARKodeSetFixedStepBounds", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_cfail_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxCFailGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKStepSetMaxCFailGrowth", flag);
#else
    int flag = ARKodeSetMaxCFailGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxCFailGrowth", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_efail_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxEFailGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKStepSetMaxEFailGrowth", flag);
#else
    int flag = ARKodeSetMaxEFailGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxEFailGrowth", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_first_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxFirstGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKStepSetMaxFirstGrowth", flag);
#else
    int flag = ARKodeSetMaxFirstGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxFirstGrowth", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_max_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				   Double_val(varg));
    CHECK_FLAG("ARKStepSetMaxGrowth", flag);
#else
    int flag = ARKodeSetMaxGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxGrowth", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_min_reduction(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 530 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMinReduction(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ARKStepSetMinReduction", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_lsetup_frequency(value varkode_mem,
						     value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetLSetupFrequency(ARKODE_MEM_FROM_ML(varkode_mem),
				         Int_val(varg));
    CHECK_FLAG("ARKStepSetLSetupFrequency", flag);
#elif 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxStepsBetweenLSet(ARKODE_MEM_FROM_ML(varkode_mem),
				             Int_val(varg));
    CHECK_FLAG("ARKStepSetMaxStepsBetweenLSet", flag);
#else
    int flag = ARKodeSetMaxStepsBetweenLSet(ARKODE_MEM_FROM_ML(varkode_mem),
				            Int_val(varg));
    CHECK_FLAG("ARKodeSetMaxStepsBetweenLSet", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_nonlin_crdown(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetNonlinCRDown(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ARKStepSetNonlinCRDown", flag);
#else
    int flag = ARKodeSetNonlinCRDown(ARKODE_MEM_FROM_ML(varkode_mem),
				     Double_val(varg));
    CHECK_FLAG("ARKodeSetNonlinCRDown", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_nonlin_rdiv(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetNonlinRDiv(ARKODE_MEM_FROM_ML(varkode_mem),
				    Double_val(varg));
    CHECK_FLAG("ARKStepSetNonlinRDiv", flag);
#else
    int flag = ARKodeSetNonlinRDiv(ARKODE_MEM_FROM_ML(varkode_mem),
				   Double_val(varg));
    CHECK_FLAG("ARKodeSetNonlinRDiv", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_optimal_params(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetOptimalParams(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKStepSetOptimalParams", flag);
#else
    int flag = ARKodeSetOptimalParams(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetOptimalParams", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_order(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetOrder(ARKODE_MEM_FROM_ML(varkode_mem), Int_val(varg));
    CHECK_FLAG("ARKStepSetOrder", flag);
#else
    int flag = ARKodeSetOrder(ARKODE_MEM_FROM_ML(varkode_mem), Int_val(varg));
    CHECK_FLAG("ARKodeSetOrder", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_interpolant_type(value varkode_mem,
						     value vinterptype)
{
    CAMLparam2(varkode_mem, vinterptype);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetInterpolantType(ARKODE_MEM_FROM_ML(varkode_mem),
		ark_interpolant_types[Int_val(vinterptype)]);
    CHECK_FLAG("ARKStepSetInterpolantType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_interpolant_degree(value varkode_mem,
						       value vinterpdegree)
{
    CAMLparam2(varkode_mem, vinterpdegree);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetInterpolantDegree(ARKODE_MEM_FROM_ML(varkode_mem),
					   Int_val(vinterpdegree));
    CHECK_FLAG("ARKStepSetInterpolantDegree", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_safety_factor(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetSafetyFactor(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ARKStepSetSafetyFactor", flag);
#else
    int flag = ARKodeSetSafetyFactor(ARKODE_MEM_FROM_ML(varkode_mem),
				     Double_val(varg));
    CHECK_FLAG("ARKodeSetSafetyFactor", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_small_num_efails(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetSmallNumEFails(ARKODE_MEM_FROM_ML(varkode_mem),
				        Int_val(varg));
    CHECK_FLAG("ARKStepSetSmallNumEFails", flag);
#else
    int flag = ARKodeSetSmallNumEFails(ARKODE_MEM_FROM_ML(varkode_mem),
				       Int_val(varg));
    CHECK_FLAG("ARKodeSetSmallNumEFails", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_gs_type(value varkode_mem, value vgstype)
{
    CAMLparam2(varkode_mem, vgstype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = ARKSpilsSetGSType(ARKODE_MEM_FROM_ML(varkode_mem),
				 sunml_lsolver_gs_type(vgstype));
    CHECK_FLAG("ARKSpilsSetGSType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_mass_gs_type(value varkode_mem, value vgstype)
{
    CAMLparam2(varkode_mem, vgstype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = ARKSpilsSetMassGSType(ARKODE_MEM_FROM_ML(varkode_mem),
				     sunml_lsolver_gs_type(vgstype));
    CHECK_FLAG("ARKSpilsSetMassGSType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_jac_eval_frequency(value varkode_mem,
						       value vevalfreq)
{
    CAMLparam2(varkode_mem, vevalfreq);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetJacEvalFrequency(ARKODE_MEM_FROM_ML(varkode_mem),
					  Long_val(vevalfreq));
    CHECK_LS_FLAG("ARKStepSetJacEvalFrequency", flag);
#elif 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMaxStepsBetweenJac(ARKODE_MEM_FROM_ML(varkode_mem),
					    Long_val(vevalfreq));
    CHECK_LS_FLAG("ARKStepSetMaxStepsBetweenJac", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_linear_solution_scaling(value varkode_mem,
							    value vonoff)
{
    CAMLparam2(varkode_mem, vonoff);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetLinearSolutionScaling(ARKODE_MEM_FROM_ML(varkode_mem),
					       Bool_val(vonoff));
    CHECK_FLAG("ARKStepSetLinearSolutionScaling", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_eps_lin(value varkode_mem, value eplifac)
{
    CAMLparam2(varkode_mem, eplifac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetEpsLin(ARKODE_MEM_FROM_ML(varkode_mem),
				Double_val(eplifac));
    CHECK_FLAG("ARKStepSetEpsLin", flag);
#else
    int flag = ARKSpilsSetEpsLin(ARKODE_MEM_FROM_ML(varkode_mem),
				 Double_val(eplifac));
    CHECK_FLAG("ARKSpilsSetEpsLin", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_ls_norm_factor(value varkode_mem, value vfac)
{
    CAMLparam2(varkode_mem, vfac);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetLSNormFactor(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(vfac));
    CHECK_FLAG("ARKStepSetLSNormFactor", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_mass_eps_lin(value varkode_mem, value eplifac)
{
    CAMLparam2(varkode_mem, eplifac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMassEpsLin(ARKODE_MEM_FROM_ML(varkode_mem),
				    Double_val(eplifac));
    CHECK_FLAG("ARKStepSetMassEpsLin", flag);
#else
    int flag = ARKSpilsSetMassEpsLin(ARKODE_MEM_FROM_ML(varkode_mem),
				     Double_val(eplifac));
    CHECK_FLAG("ARKSpilsSetMassEpsLin", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_set_mass_ls_norm_factor(value varkode_mem,
							value vfac)
{
    CAMLparam2(varkode_mem, vfac);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepSetMassLSNormFactor(ARKODE_MEM_FROM_ML(varkode_mem),
				          Double_val(vfac));
    CHECK_FLAG("ARKStepSetMassLSNormFactor", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_maxl(value varkode_mem, value maxl)
{
    CAMLparam2(varkode_mem, maxl);
#if SUNDIALS_LIB_VERSION < 300
    int flag = ARKSpilsSetMaxl(ARKODE_MEM_FROM_ML(varkode_mem), Int_val(maxl));
    CHECK_FLAG("ARKSpilsSetMaxl", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_mass_maxl(value varkode_mem, value maxl)
{
    CAMLparam2(varkode_mem, maxl);
#if SUNDIALS_LIB_VERSION < 300
    int flag = ARKSpilsSetMassMaxl(ARKODE_MEM_FROM_ML(varkode_mem),
				   Int_val(maxl));
    CHECK_FLAG("ARKSpilsSetMassMaxl", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_prec_type(value varkode_mem, value vptype)
{
    CAMLparam2(varkode_mem, vptype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = ARKSpilsSetPrecType(ARKODE_MEM_FROM_ML(varkode_mem),
				   sunml_lsolver_precond_type(vptype));
    CHECK_FLAG("ARKSpilsSetPrecType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_mass_prec_type(value varkode_mem,
						 value vptype)
{
    CAMLparam2(varkode_mem, vptype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = ARKSpilsSetMassPrecType(ARKODE_MEM_FROM_ML(varkode_mem),
				       sunml_lsolver_precond_type(vptype));
    CHECK_FLAG("ARKSpilsSetMassPrecType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* statistic accessor functions */

CAMLprim value sunml_arkode_ark_get_tol_scale_factor(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    realtype r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetTolScaleFactor(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetTolScaleFactor", flag);
#else
    int flag = ARKodeGetTolScaleFactor(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetTolScaleFactor", flag);
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_arkode_ark_get_num_nonlin_solv_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumNonlinSolvIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumNonlinSolvIters", flag);
#else
    int flag = ARKodeGetNumNonlinSolvIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetNumNonlinSolvIters", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_nonlin_solv_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumNonlinSolvConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumNonlinSolvConvFails", flag);
#else
    int flag = ARKodeGetNumNonlinSolvConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetNumNonlinSolvConvFails", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_nonlin_solv_stats(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    long int nniters, nncfails;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNonlinSolvStats(ARKODE_MEM_FROM_ML(varkode_mem),
				         &nniters, &nncfails);
    CHECK_FLAG("ARKStepGetNonlinSolvStats", flag);
#else
    int flag = ARKodeGetNonlinSolvStats(ARKODE_MEM_FROM_ML(varkode_mem),
				       &nniters, &nncfails);
    CHECK_FLAG("ARKodeGetNonlinSolvStats", flag);
#endif

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nniters));
    Store_field(r, 1, Val_long(nncfails));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_get_num_g_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumGEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumGEvals", flag);
#else
    int flag = ARKodeGetNumGEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetNumGEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_constr_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumConstrFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumConstrFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_lin_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetLinWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				      &lenrwLS, &leniwLS);
    CHECK_FLAG("ARKStepGetLinWorkSpace", flag);
#else
    int flag = ARKDlsGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				  &lenrwLS, &leniwLS);
    CHECK_FLAG("ARKDlsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_dls_get_mass_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetMassWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				       &lenrwLS, &leniwLS);
    CHECK_FLAG("ARKStepGetMassWorkSpace", flag);
#else
    int flag = ARKDlsGetMassWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				      &lenrwLS, &leniwLS);
    CHECK_FLAG("ARKDlsGetMassWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_get_num_jac_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumJacEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumJacEvals", flag);
#else
    int flag = ARKDlsGetNumJacEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumJacEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_mass_setups(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassSetups(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassSetups", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    int flag = ARKDlsGetNumMassSetups(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumMassSetups", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_mass_mult_setups(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassMultSetups(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassMultSetups", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_mass_solves(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassSolves", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    int flag = ARKDlsGetNumMassSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumMassSolves", flag);
#else
    int flag = ARKodeGetNumMassSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetNumMassSolves", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_dls_get_num_mass_mult(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassMult(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassMult", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    int flag = ARKDlsGetNumMassMult(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumMassMult", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_dls_get_num_lin_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumLinRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumLinRhsEvals", flag);
#else
    int flag = ARKDlsGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumRhsEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_bandprec_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    long int lenrwBP;
    long int leniwBP;

    int flag = ARKBandPrecGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				       &lenrwBP, &leniwBP);
    CHECK_FLAG("ARKBandPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrwBP));
    Store_field(r, 1, Val_long(leniwBP));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_bandprec_get_num_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKBandPrecGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKBandPrecGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

/* spils functions */

CAMLprim value sunml_arkode_spils_spgmr (value varkode_mem,
				     value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKSpgmr (arkode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("ARKSpgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_spbcgs (value varkode_mem,
				      value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKSpbcg (arkode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("ARKSpbcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_sptfqmr (value varkode_mem, value vmaxl,
				      value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKSptfqmr (arkode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("ARKSptfqmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_spfgmr (value varkode_mem, value vmaxl,
				      value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKSpfgmr (arkode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("ARKSpfgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_pcg (value varkode_mem, value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKPcg (arkode_mem, sunml_lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("ARKPcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_num_lin_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumLinIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumLinIters", flag);
#else
    int flag = ARKSpilsGetNumLinIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumLinIters", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_lin_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumLinConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumLinConvFails", flag);
#else
    int flag = ARKSpilsGetNumConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumConvFails", flag);
#endif

    CAMLreturn(Val_long(r));
}

#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION <= 301
int ARKSpilsGetNumMTSetups(void *arkode_mem, long int *nmtsetups);
#endif

CAMLprim value sunml_arkode_ark_get_num_mtsetups(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if   400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMTSetups(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMTSetups", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    int flag = ARKSpilsGetNumMTSetups(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMTSetups", flag);
#else
    r = 0;
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_mass_mult(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassMult(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassMult", flag);
#elif 263 <= SUNDIALS_LIB_VERSION
    int flag = ARKSpilsGetNumMtimesEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMtimesEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetLinWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				  &lenrw, &leniw);
    CHECK_FLAG("ARKStepGetLinWorkSpace", flag);
#else
    flag = ARKSpilsGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				&lenrw, &leniw);
    CHECK_FLAG("ARKSpilsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_get_num_prec_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumPrecEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumPrecEvals", flag);
#else
    int flag = ARKSpilsGetNumPrecEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumPrecEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_prec_solves(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumPrecSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumPrecSolves", flag);
#else
    int flag = ARKSpilsGetNumPrecSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumPrecSolves", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_jtsetup_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if   400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumJTSetupEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumJTSetupEvals", flag);
#elif 300 <= SUNDIALS_LIB_VERSION
    int flag = ARKSpilsGetNumJTSetupEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumJTSetupEvals", flag);
#else
    r = 0;
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_jtimes_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumJtimesEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumJtimesEvals", flag);
#else
    int flag = ARKSpilsGetNumJtimesEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumJtimesEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_lin_rhs_evals (value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumLinRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumLinRhsEvals", flag);
#else
    int flag = ARKSpilsGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumRhsEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

/* spils mass functions */

CAMLprim value sunml_arkode_spils_mass_spgmr (value varkode_mem,
					  value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSpgmr (arkode_mem, sunml_lsolver_precond_type (vtype),
				     Int_val (vmaxl),
				     masstimesfn,
				     ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG ("ARKMassSpgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_mass_spbcgs (value varkode_mem,
					   value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSpbcg (arkode_mem, sunml_lsolver_precond_type (vtype),
				     Int_val (vmaxl),
				     masstimesfn,
				     ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG ("ARKMassSpbcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_mass_sptfqmr (value varkode_mem,
					    value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSptfqmr (arkode_mem, sunml_lsolver_precond_type (vtype),
				       Int_val (vmaxl),
				       masstimesfn,
				       ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG ("ARKMassSptfqmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_mass_spfgmr (value varkode_mem,
					   value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSpfgmr (arkode_mem, sunml_lsolver_precond_type (vtype),
				      Int_val (vmaxl),
				      masstimesfn,
				      ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG ("ARKMassSpfgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_mass_pcg (value varkode_mem,
					value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassPcg (arkode_mem, sunml_lsolver_precond_type (vtype),
				   Int_val (vmaxl),
				   masstimesfn,
				   ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG ("ARKMassPcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ark_get_num_mass_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassIters", flag);
#else
    int flag = ARKSpilsGetNumMassIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMassIters", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_mass_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassConvFails", flag);
#else
    int flag = ARKSpilsGetNumMassConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMassConvFails", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_mass_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

#if 400 <= SUNDIALS_LIB_VERSION
    flag = ARKStepGetMassWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				    &lenrw, &leniw);
    CHECK_FLAG("ARKStepGetMassWorkSpace", flag);
#else
    flag = ARKSpilsGetMassWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				    &lenrw, &leniw);
    CHECK_FLAG("ARKSpilsGetMassWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_ark_get_num_mass_prec_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassPrecEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassPrecEvals", flag);
#else
    int flag = ARKSpilsGetNumMassPrecEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMassPrecEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_get_num_mass_prec_solves(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepGetNumMassPrecSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKStepGetNumMassPrecSolves", flag);
#else
    int flag = ARKSpilsGetNumMassPrecSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMassPrecSolves", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_ark_write_parameters(value varkode_mem, value vlog)
{
    CAMLparam2(varkode_mem, vlog);
#if 410 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepWriteParameters(ARKODE_MEM_FROM_ML(varkode_mem),
				      ML_CFILE(vlog));
    CHECK_FLAG("ARKStepWriteParameters", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_unit);
}

CAMLprim value sunml_arkode_ark_write_butcher(value varkode_mem, value vlog)
{
    CAMLparam2(varkode_mem, vlog);
#if 410 <= SUNDIALS_LIB_VERSION
    int flag = ARKStepWriteButcher(ARKODE_MEM_FROM_ML(varkode_mem),
				   ML_CFILE(vlog));
    CHECK_FLAG("ARKStepWriteButcher", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * ERKStep basic interface
 */

/* ERKStepCreate() */
CAMLprim value sunml_arkode_erk_init(value weakref, value y0, value t0)
{
    CAMLparam3(weakref, y0, t0);
    CAMLlocal2(r, varkode_mem);
#if 400 <= SUNDIALS_LIB_VERSION
    value *backref;

    void *arkode_mem = ERKStepCreate(rhsfn1, Double_val(t0), NVEC_VAL(y0));

    if (arkode_mem == NULL)
	caml_failwith("ERKStepCreate returned NULL");

    varkode_mem = caml_alloc_final(1, NULL, 1, 5);
    ARKODE_MEM(varkode_mem) = arkode_mem;

    backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	ERKStepFree (&arkode_mem);
	caml_raise_out_of_memory();
    }
    ERKStepSetUserData (arkode_mem, backref);

    r = caml_alloc_tuple (2);
    Store_field (r, 0, varkode_mem);
    Store_field (r, 1, (value)backref);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

/* Set the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim value sunml_arkode_erk_root_init (value vdata, value vnroots)
{
    CAMLparam2 (vdata, vnroots);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (vdata);
    int nroots = Int_val (vnroots);
    int flag = ERKStepRootInit (arkode_mem, nroots, roots);
    CHECK_FLAG ("ERKStepRootInit", flag);
    Store_field (vdata, RECORD_ARKODE_SESSION_NROOTS, vnroots);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_sv_tolerances(value vdata, value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = ERKStepSVtolerances(ARKODE_MEM_FROM_ML(vdata),
			 	   Double_val(reltol), atol_nv);
    CHECK_FLAG("ERKStepSVtolerances", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_reinit(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepReInit(ARKODE_MEM_FROM_ML(vdata),
			     rhsfn1, Double_val(t0), NVEC_VAL(y0));
    CHECK_FLAG("ERKStepReInit", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_reset(value vdata, value vt, value vy)
{
    CAMLparam3(vdata, vt, vy);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepReset(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			    NVEC_VAL(vy));
    CHECK_FLAG("ERKStepReset", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

static value erk_solver(value vdata, value nextt, value vy, int onestep)
{
    CAMLparam3(vdata, nextt, vy);
    CAMLlocal1(ret);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype tret;
    int flag;
    N_Vector y;
    enum arkode_solver_result_tag result = -1;

    y = NVEC_VAL (vy);
    // Caml_ba_data_val(y) must not be shifted by the OCaml GC during this
    // function call, which calls Caml through the callback f.  Is this
    // guaranteed?
    flag = ERKStepEvolve(ARKODE_MEM_FROM_ML (vdata), Double_val (nextt),
			 y, &tret, onestep ? ARK_ONE_STEP : ARK_NORMAL);

    switch (flag) {
    case ARK_SUCCESS:
	result = VARIANT_ARKODE_SOLVER_RESULT_SUCCESS;
	break;

    case ARK_ROOT_RETURN:
	result = VARIANT_ARKODE_SOLVER_RESULT_ROOTSFOUND;
	break;

    case ARK_TSTOP_RETURN:
	result = VARIANT_ARKODE_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	sunml_arkode_check_flag("ERKStepEvolve", flag,
				ARKODE_MEM_FROM_ML(vdata));
    }

    assert (Field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (ret);
}

CAMLprim value sunml_arkode_erk_solve_normal(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(erk_solver(vdata, nextt, y, 0));
}

CAMLprim value sunml_arkode_erk_solve_one_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(erk_solver(vdata, nextt, y, 1));
}

CAMLprim value sunml_arkode_erk_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetDky(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			     Int_val(vk), NVEC_VAL(vy));
    CHECK_FLAG("ERKStepGetDky", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_session_finalize(value vdata)
{
#if 400 <= SUNDIALS_LIB_VERSION
    if (ARKODE_MEM_FROM_ML(vdata) != NULL) {
	void *arkode_mem = ARKODE_MEM_FROM_ML(vdata);
	value *backref = ARKODE_BACKREF_FROM_ML(vdata);
	ERKStepFree(&arkode_mem);
	sunml_sundials_free_value(backref);
    }
#endif
    return Val_unit;
}

CAMLprim value sunml_arkode_erk_ss_tolerances(value vdata, value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSStolerances(ARKODE_MEM_FROM_ML(vdata),
				   Double_val(reltol),
				   Double_val(abstol));
    CHECK_FLAG("ERKStepSStolerances", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_get_timestepper_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);
#if 400 <= SUNDIALS_LIB_VERSION

    int flag;

    long int expsteps;
    long int accsteps;
    long int step_attempts;
    long int nf_evals;
    long int netfails;

    flag = ERKStepGetTimestepperStats(ARKODE_MEM_FROM_ML(vdata),
                                      &expsteps,
                                      &accsteps,
                                      &step_attempts,
                                      &nf_evals,
                                      &netfails);
    CHECK_FLAG("ERKStepGetTimestepperStats", flag);

    r = caml_alloc_tuple(RECORD_ARKODE_ERK_TIMESTEPPER_STATS_SIZE);
    Store_field(r, RECORD_ARKODE_ERK_TIMESTEPPER_STATS_EXP_STEPS,
		   Val_long(expsteps));
    Store_field(r, RECORD_ARKODE_ERK_TIMESTEPPER_STATS_ACC_STEPS,
		   Val_long(accsteps));
    Store_field(r, RECORD_ARKODE_ERK_TIMESTEPPER_STATS_STEP_ATTEMPTS,
		   Val_long(step_attempts));
    Store_field(r, RECORD_ARKODE_ERK_TIMESTEPPER_STATS_NUM_NF_EVALS,
		   Val_long(nf_evals));
    Store_field(r, RECORD_ARKODE_ERK_TIMESTEPPER_STATS_NUM_ERR_TEST_FAILS,
		   Val_long(netfails));

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_arkode_erk_get_step_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);
#if 400 <= SUNDIALS_LIB_VERSION
    long int nsteps;
    realtype hinused;
    realtype hlast;
    realtype hcur;
    realtype tcur;

    int flag = ERKStepGetStepStats(ARKODE_MEM_FROM_ML(vdata),
				   &nsteps,
				   &hinused,
				   &hlast,
				   &hcur,
				   &tcur);
    CHECK_FLAG("ERKStepGetStepStats", flag);

    r = caml_alloc_tuple(RECORD_ARKODE_STEP_STATS_SIZE);
    Store_field(r, RECORD_ARKODE_STEP_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_ARKODE_STEP_STATS_ACTUAL_INIT_STEP,
						    caml_copy_double(hinused));
    Store_field(r, RECORD_ARKODE_STEP_STATS_LAST_STEP,
						    caml_copy_double(hlast));
    Store_field(r, RECORD_ARKODE_STEP_STATS_CURRENT_STEP,
						    caml_copy_double(hcur));
    Store_field(r, RECORD_ARKODE_STEP_STATS_CURRENT_TIME,
						    caml_copy_double(tcur));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_arkode_erk_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);
#if 400 <= SUNDIALS_LIB_VERSION
    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = ERKStepSetRootDirection(ARKODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("ERKStepSetRootDirection", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_resize(value varkode_mem,
				       value vhasfn, value vhscale,
				       value vt0, value vynew)
{
    CAMLparam5(varkode_mem, vhasfn, vhscale, vt0, vynew);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ERKStepResize(
		 arkode_mem,
		 NVEC_VAL(vynew),
		 Double_val(vhscale),
		 Double_val(vt0),
		 Bool_val(vhasfn) ? resizefn : NULL,
		 Bool_val(vhasfn) ? ARKODE_BACKREF_FROM_ML(varkode_mem) : NULL);
    CHECK_FLAG("ERKStepResize", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_get_current_butcher_table(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vr);
#if 400 <= SUNDIALS_LIB_VERSION
    ARKodeButcherTable bt;
    int flag = ERKStepGetCurrentButcherTable(ARKODE_MEM_FROM_ML(varkode_mem),
					     &bt);
    CHECK_FLAG("ERKStepGetCurrentButcherTable", flag);

    vr = val_butcher_table(bt);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_erk_set_table(value varkode_mem, value vob)
{
    CAMLparam2(varkode_mem, vob);
#if 400 <= SUNDIALS_LIB_VERSION
    ARKodeButcherTable bt = butcher_table_val(vob);

    int flag = ERKStepSetTable(ARKODE_MEM_FROM_ML(varkode_mem), bt);
    CHECK_FLAG("ERKStepSetTable", flag);

    if (bt != NULL) ARKodeButcherTable_Free(bt);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_table_num(value varkode_mem, value vnum)
{
    CAMLparam2(varkode_mem, vnum);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetTableNum(ARKODE_MEM_FROM_ML(varkode_mem),
				  Int_val(vnum));
    CHECK_FLAG("ERKStepSetTableNum", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_adaptivity_method(value varkode_mem,
						      value vmeth)
{
    CAMLparam2(varkode_mem, vmeth);
    CAMLlocal2(vks, vorder);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    if (Tag_val(vmeth) == VARIANT_ARKODE_ADAPTIVITY_METHOD_ADAPTIVITYFN) {
	flag = ERKStepSetAdaptivityFn(arkode_mem,
				      adaptfn,
				      ARKODE_BACKREF_FROM_ML(varkode_mem));
	CHECK_FLAG("ERKStepSetAdaptivityFn", flag);

    } else {
	realtype adapt_params[3] = { 0 };

	vks = Field(Field(vmeth, 0), RECORD_ARKODE_ADAPTIVITY_PARAMS_KS);
	vorder = Field(Field(vmeth, 0),
			RECORD_ARKODE_ADAPTIVITY_PARAMS_METHOD_ORDER);

	if (vks != Val_none) {
	    adapt_params[0] = Double_val(Field(Some_val(vks), 0));
	    adapt_params[1] = Double_val(Field(Some_val(vks), 1));
	    adapt_params[2] = Double_val(Field(Some_val(vks), 2));
	}

	flag = ERKStepSetAdaptivityMethod(arkode_mem, 
					  Tag_val(vmeth), 
					  vks == Val_none,
					  Bool_val(vorder),
					  vks != Val_none ? adapt_params : NULL);
	CHECK_FLAG("ERKStepSetAdaptivityMethod", flag);
    }

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_stability_fn(value varkode_mem, value vhasf)
{
    CAMLparam2(varkode_mem, vhasf);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ERKStepSetStabilityFn(arkode_mem,
				     Bool_val(vhasf) ? stabfn : NULL,
				     ARKODE_BACKREF_FROM_ML(varkode_mem));
    CHECK_FLAG("ERKStepSetStabilityFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_postprocess_step_fn(value varkode_mem,
							value vhasf)
{
    CAMLparam2(varkode_mem, vhasf);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ERKStepSetPostprocessStepFn(arkode_mem,
					   Bool_val(vhasf) ? poststepfn : NULL);
    CHECK_FLAG("ERKStepSetPostprocessStepFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Boiler plate definitions for ERKStep interface.
 */

/* main solver optional output functions */

CAMLprim value sunml_arkode_erk_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);
#if 400 <= SUNDIALS_LIB_VERSION

    int flag;
    long int lenrw;
    long int leniw;

    flag = ERKStepGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem), &lenrw, &leniw);
    CHECK_FLAG("ERKStepGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_arkode_erk_get_num_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int v = 0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetNumSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetNumSteps", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_erk_get_actual_init_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype v = 0.0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetActualInitStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetActualInitStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_erk_get_last_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype v = 0.0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetLastStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetLastStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_erk_get_current_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype v = 0.0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetCurrentStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetCurrentStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_erk_get_current_time(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype v = 0.0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetCurrentTime(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetCurrentTime", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_erk_get_err_weights(value varkode_mem, value verrws)
{
    CAMLparam2(varkode_mem, verrws);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector errws_nv = NVEC_VAL(verrws);
    int flag = ERKStepGetErrWeights(ARKODE_MEM_FROM_ML(varkode_mem), errws_nv);
    CHECK_FLAG("ERKStepGetErrWeights", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_get_tol_scale_factor(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype r = 0.0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetTolScaleFactor(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ERKStepGetTolScaleFactor", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_arkode_erk_get_num_exp_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int v = 0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetNumExpSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetNumExpSteps", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_erk_get_num_acc_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int v = 0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetNumAccSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetNumAccSteps", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_erk_get_num_step_attempts(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int v = 0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetNumStepAttempts(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetNumStepAttempts", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_erk_get_num_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

#if 400 <= SUNDIALS_LIB_VERSION
    long int nf;

    int flag = ERKStepGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &nf);
    CHECK_FLAG("ERKStepGetNumRhsEvals", flag);

    r = Val_long(nf);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_erk_get_num_err_test_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int v = 0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetNumErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ERKStepGetNumErrTestFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_erk_get_est_local_errors(value varkode_mem,
						     value vele)
{
    CAMLparam2(varkode_mem, vele);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetEstLocalErrors(ARKODE_MEM_FROM_ML(varkode_mem),
					NVEC_VAL(vele));
    CHECK_FLAG("ERKStepGetEstLocalErrors", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* optional inputs for ERKStep */

CAMLprim value sunml_arkode_erk_set_defaults(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetDefaults(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ERKStepSetDefaults", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}


CAMLprim value sunml_arkode_erk_set_diagnostics(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ERKStepSetDiagnostics", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_clear_diagnostics(value vdata)
{
    CAMLparam1(vdata);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("ERKStepSetDiagnostics", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_error_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetErrFile(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ERKStepSetErrFile", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), errh,
				      ARKODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("ERKStepSetErrHandlerFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("ERKStepSetErrHandlerFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_fixed_step(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetFixedStep(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ERKStepSetFixedStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_num_constr_fails(value varkode_mem, 
							 value vmaxfails)
{
    CAMLparam2(varkode_mem, vmaxfails);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxNumConstrFails(ARKODE_MEM_FROM_ML(varkode_mem),
					   Int_val(vmaxfails));
    CHECK_FLAG("ERKStepSetMaxNumConstrFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_init_step(value varkode_mem, value hin)
{
    CAMLparam2(varkode_mem, hin);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetInitStep(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(hin));
    CHECK_FLAG("ERKStepSetInitStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_hnil_warns(value varkode_mem, value mxhnil)
{
    CAMLparam2(varkode_mem, mxhnil);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxHnilWarns(ARKODE_MEM_FROM_ML(varkode_mem),
				      Int_val(mxhnil));
    CHECK_FLAG("ERKStepSetMaxHnilWarns", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_num_steps(value varkode_mem, value mxsteps)
{
    CAMLparam2(varkode_mem, mxsteps);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxNumSteps(ARKODE_MEM_FROM_ML(varkode_mem),
				     Long_val(mxsteps));
    CHECK_FLAG("ERKStepSetMaxNumSteps", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_step(value varkode_mem, value hmax)
{
    CAMLparam2(varkode_mem, hmax);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxStep(ARKODE_MEM_FROM_ML(varkode_mem),
		 		 Double_val(hmax));
    CHECK_FLAG("ERKStepSetMaxStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_min_step(value varkode_mem, value hmin)
{
    CAMLparam2(varkode_mem, hmin);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMinStep(ARKODE_MEM_FROM_ML(varkode_mem),
		 		 Double_val(hmin));
    CHECK_FLAG("ERKStepSetMinStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_stop_time(value varkode_mem, value tstop)
{
    CAMLparam2(varkode_mem, tstop);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetStopTime(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(tstop));
    CHECK_FLAG("ERKStepSetStopTime", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_err_test_fails(value varkode_mem, value maxnef)
{
    CAMLparam2(varkode_mem, maxnef);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem),
					 Int_val(maxnef));
    CHECK_FLAG("ERKStepSetMaxErrTestFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

/* Optional inputs for time-step adaptivity */

CAMLprim value sunml_arkode_erk_set_cfl_fraction(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetCFLFraction(ARKODE_MEM_FROM_ML(varkode_mem),
				    Double_val(varg));
    CHECK_FLAG("ERKStepSetCFLFraction", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_error_bias(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetErrorBias(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ERKStepSetErrorBias", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_fixed_step_bounds(value varkode_mem,
						      value vlb, value vub)
{
    CAMLparam3(varkode_mem, vlb, vub);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetFixedStepBounds(ARKODE_MEM_FROM_ML(varkode_mem),
				        Double_val(vlb), Double_val(vub));
    CHECK_FLAG("ERKStepSetFixedStepBounds", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_efail_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxEFailGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ERKStepSetMaxEFailGrowth", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_first_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxFirstGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ERKStepSetMaxFirstGrowth", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_max_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMaxGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				   Double_val(varg));
    CHECK_FLAG("ERKStepSetMaxGrowth", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_min_reduction(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 530 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetMinReduction(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ERKStepSetMinReduction", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_safety_factor(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetSafetyFactor(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ERKStepSetSafetyFactor", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_small_num_efails(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetSmallNumEFails(ARKODE_MEM_FROM_ML(varkode_mem),
				        Int_val(varg));
    CHECK_FLAG("ERKStepSetSmallNumEFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_order(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetOrder(ARKODE_MEM_FROM_ML(varkode_mem), Int_val(varg));
    CHECK_FLAG("ERKStepSetOrder", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_interpolant_type(value varkode_mem,
						     value vinterptype)
{
    CAMLparam2(varkode_mem, vinterptype);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetInterpolantType(ARKODE_MEM_FROM_ML(varkode_mem),
		ark_interpolant_types[Int_val(vinterptype)]);
    CHECK_FLAG("ERKStepSetInterpolantType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_interpolant_degree(value varkode_mem,
						       value vinterpdegree)
{
    CAMLparam2(varkode_mem, vinterpdegree);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetInterpolantDegree(ARKODE_MEM_FROM_ML(varkode_mem),
					   Int_val(vinterpdegree));
    CHECK_FLAG("ERKStepSetInterpolantDegree", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);
#if 400 <= SUNDIALS_LIB_VERSION
    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = ERKStepGetRootInfo(ARKODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("ERKStepGetRootInfo", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_get_num_g_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetNumGEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ERKStepGetNumGEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_erk_get_num_constr_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepGetNumConstrFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ERKStepGetNumConstrFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_erk_set_constraints(value varkode_mem, value vnv)
{
    CAMLparam2(varkode_mem, vnv);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetConstraints(ARKODE_MEM_FROM_ML(varkode_mem),
				     NVEC_VAL(vnv));
    CHECK_FLAG("ERKStepSetConstraints", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_set_no_inactive_root_warn(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepSetNoInactiveRootWarn(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ERKStepSetNoInactiveRootWarn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_erk_write_parameters(value varkode_mem, value vlog)
{
    CAMLparam2(varkode_mem, vlog);
#if 410 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepWriteParameters(ARKODE_MEM_FROM_ML(varkode_mem),
				      ML_CFILE(vlog));
    CHECK_FLAG("ERKStepWriteParameters", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_unit);
}

CAMLprim value sunml_arkode_erk_write_butcher(value varkode_mem, value vlog)
{
    CAMLparam2(varkode_mem, vlog);
#if 410 <= SUNDIALS_LIB_VERSION
    int flag = ERKStepWriteButcher(ARKODE_MEM_FROM_ML(varkode_mem),
				   ML_CFILE(vlog));
    CHECK_FLAG("ERKStepWriteButcher", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * MRIStep basic interface
 */

/* MRIStepCreate() */
CAMLprim value sunml_arkode_mri_init(value weakref, value vistepper,
				     value y0, value t0)
{
    CAMLparam4(weakref, vistepper, y0, t0);
    CAMLlocal3(r, varkode_mem, vistepper_val);
#if 500 <= SUNDIALS_LIB_VERSION
    value *backref;
    void *arkode_mem;

    vistepper_val = Field(vistepper, RECORD_ARKODE_MRI_ISTEPPER_VAL);

    if (Is_block(vistepper_val)
	&& Tag_val(vistepper_val) == VARIANT_ARKODE_MRI_ISTEPPER_ARKSTEP)
    {
	arkode_mem =
	    MRIStepCreate(rhsfn1, Double_val(t0), NVEC_VAL(y0),
			  MRISTEP_ARKSTEP,
			  ARKODE_MEM_FROM_ML(Field(vistepper_val, 0)));
    } else {
#if 580 <= SUNDIALS_LIB_VERSION
	arkode_mem =
	    MRIStepCreate(rhsfn1, Double_val(t0), NVEC_VAL(y0),
			  MRISTEP_CUSTOM,
			  ISTEPPER_FROM_ML(vistepper));
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    }

    if (arkode_mem == NULL)
	caml_failwith("MRIStepCreate returned NULL");

    varkode_mem = caml_alloc_final(1, NULL, 1, 5);
    ARKODE_MEM(varkode_mem) = arkode_mem;

    backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	MRIStepFree (&arkode_mem);
	caml_raise_out_of_memory();
    }
    MRIStepSetUserData (arkode_mem, backref);

    r = caml_alloc_tuple (2);
    Store_field (r, 0, varkode_mem);
    Store_field (r, 1, (value)backref);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

/* Set the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim value sunml_arkode_mri_root_init (value vdata, value vnroots)
{
    CAMLparam2 (vdata, vnroots);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (vdata);
    int nroots = Int_val (vnroots);
    int flag = MRIStepRootInit (arkode_mem, nroots, roots);
    CHECK_FLAG ("MRIStepRootInit", flag);
    Store_field (vdata, RECORD_ARKODE_SESSION_NROOTS, vnroots);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_reinit(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);
#if 500 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepReInit(ARKODE_MEM_FROM_ML(vdata),
			     rhsfn1, Double_val(t0), NVEC_VAL(y0));
    CHECK_FLAG("MRIStepReInit", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_reset(value vdata, value vt, value vy)
{
    CAMLparam3(vdata, vt, vy);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepReset(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			    NVEC_VAL(vy));
    CHECK_FLAG("MRIStepReset", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

static value mri_solver(value vdata, value nextt, value vy, int onestep)
{
    CAMLparam3(vdata, nextt, vy);
    CAMLlocal1(ret);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype tret;
    int flag;
    N_Vector y;
    enum arkode_solver_result_tag result = -1;

    y = NVEC_VAL (vy);
    // Caml_ba_data_val(y) must not be shifted by the OCaml GC during this
    // function call, which calls Caml through the callback f.  Is this
    // guaranteed?
    flag = MRIStepEvolve(ARKODE_MEM_FROM_ML (vdata), Double_val (nextt),
			 y, &tret, onestep ? ARK_ONE_STEP : ARK_NORMAL);

    switch (flag) {
    case ARK_SUCCESS:
	result = VARIANT_ARKODE_SOLVER_RESULT_SUCCESS;
	break;

    case ARK_ROOT_RETURN:
	result = VARIANT_ARKODE_SOLVER_RESULT_ROOTSFOUND;
	break;

    case ARK_TSTOP_RETURN:
	result = VARIANT_ARKODE_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	sunml_arkode_check_flag("MRIStepEvolve", flag,
				ARKODE_MEM_FROM_ML(vdata));
    }

    assert (Field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (ret);
}

CAMLprim value sunml_arkode_mri_solve_normal(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(mri_solver(vdata, nextt, y, 0));
}

CAMLprim value sunml_arkode_mri_solve_one_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(mri_solver(vdata, nextt, y, 1));
}

CAMLprim value sunml_arkode_mri_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetDky(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			     Int_val(vk), NVEC_VAL(vy));
    CHECK_FLAG("MRIStepGetDky", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_session_finalize(value vdata)
{
#if 400 <= SUNDIALS_LIB_VERSION
    if (ARKODE_MEM_FROM_ML(vdata) != NULL) {
	void *arkode_mem = ARKODE_MEM_FROM_ML(vdata);
	value *backref = ARKODE_BACKREF_FROM_ML(vdata);
	MRIStepFree(&arkode_mem);
	sunml_sundials_free_value(backref);
    }
#endif
    return Val_unit;
}

CAMLprim value sunml_arkode_mri_ss_tolerances(value vdata, value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSStolerances(ARKODE_MEM_FROM_ML(vdata),
				   Double_val(reltol),
				   Double_val(abstol));
    CHECK_FLAG("MRIStepSStolerances", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_sv_tolerances(value vdata, value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSVtolerances(ARKODE_MEM_FROM_ML(vdata),
			 	   Double_val(reltol), atol_nv);
    CHECK_FLAG("MRIStepSVtolerances", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);
#if 400 <= SUNDIALS_LIB_VERSION
    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = MRIStepSetRootDirection(ARKODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("MRIStepSetRootDirection", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_resize(value varkode_mem,
				       value vhasfn, value vt0, value vynew)
{
    CAMLparam4(varkode_mem, vhasfn, vt0, vynew);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = MRIStepResize(
		 arkode_mem,
		 NVEC_VAL(vynew),
		 Double_val(vt0),
		 Bool_val(vhasfn) ? resizefn : NULL,
		 Bool_val(vhasfn) ? ARKODE_BACKREF_FROM_ML(varkode_mem) : NULL);
    CHECK_FLAG("MRIStepResize", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_table(value varkode_mem,
					  value vq, value vobt)
{
    CAMLparam3(varkode_mem, vq, vobt);
#if 500 <= SUNDIALS_LIB_VERSION
    ARKodeButcherTable bt = butcher_table_val(vobt);

    int flag = MRIStepSetTable(ARKODE_MEM_FROM_ML(varkode_mem),
			       Int_val(vq), bt);
    CHECK_FLAG("MRIStepSetTable", flag);

    if (bt != NULL) ARKodeButcherTable_Free(bt);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_table_num(value varkode_mem, value vt)
{
    CAMLparam2(varkode_mem, vt);
#if 500 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetTableNum(ARKODE_MEM_FROM_ML(varkode_mem), Int_val(vt));
    CHECK_FLAG("MRIStepSetTableNum", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_postprocess_step_fn(value varkode_mem,
							value vhasf)
{
    CAMLparam2(varkode_mem, vhasf);
#if 400 <= SUNDIALS_LIB_VERSION
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = MRIStepSetPostprocessStepFn(arkode_mem,
					   Bool_val(vhasf) ? poststepfn : NULL);
    CHECK_FLAG("MRIStepSetPostprocessStepFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_nls_rhs_fn(value varkode_mem)
{
    CAMLparam1(varkode_mem);
#if 580 <= SUNDIALS_LIB_VERSION
    int flag;

    flag = MRIStepSetNlsRhsFn(ARKODE_MEM_FROM_ML (varkode_mem), nlsrhsfn);
    CHECK_FLAG ("MRIStepSetNlsRhsFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* MRIStep: StepCoupling */

#define ML_MRI_COUPLING(v) (*(MRIStepCoupling *)Data_custom_val(v))

#if 540 <= SUNDIALS_LIB_VERSION
static void finalize_mri_coupling(value vcptr)
{

    MRIStepCoupling MRIC = ML_MRI_COUPLING(vcptr);
    int i;

    if (MRIC != NULL) {
	for (i = 0; i < MRIC->nmat; i++) {
	    free(MRIC->G[i]);
	    MRIC->G[i] = NULL;
	}
	free(MRIC->G);
	MRIC->G = NULL;
	free(MRIC);
    }
}
#endif

CAMLprim value sunml_arkode_mri_coupling_make(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal1(vcptr);
#if 540 <= SUNDIALS_LIB_VERSION
    CAMLlocal3(vg, vgi, vc);
    int nmat, stages, q, p;
    int i, j;
    MRIStepCoupling MRIC;

    nmat   = Int_val(Field(vargs, 0));
    stages = Int_val(Field(vargs, 1));
    q      = Int_val(Field(vargs, 2));
    p      = Int_val(Field(vargs, 3));
    vg     = Field(vargs, 4);
    vc     = Field(vargs, 5);

    // adapted from MRIStepCoupling_Alloc

    MRIC = calloc(1, sizeof(*MRIC));
    if (MRIC == NULL) caml_raise_out_of_memory();

    MRIC->nmat = nmat;
    MRIC->stages = stages;
    MRIC->q = q;
    MRIC->p = p;

    MRIC->G = (realtype ***) calloc( nmat, sizeof(realtype**) );
    if (MRIC->G == NULL) {
	MRIStepCoupling_Free(MRIC);
	caml_raise_out_of_memory();
    }

    // allocate arrays in C
    for (i = 0; i < nmat; i++) {
	MRIC->G[i] = (realtype **) calloc( stages, sizeof(realtype*) );
	if (MRIC->G[i] == NULL) {
	    MRIStepCoupling_Free(MRIC);
	    caml_raise_out_of_memory();
	}
    }

    // link nested elements to underlying big arrays
    for (i = 0; i < nmat; i++) {
	vgi = Field(vg, i);
	for (j = 0; j < stages; j++) {
	    MRIC->G[i][j] = REAL_ARRAY(Field(vgi, j));
	}
    }

    MRIC->c = REAL_ARRAY(vc);

    // Create a cptr wrapper
    vcptr = caml_alloc_final(1, &finalize_mri_coupling, 1, 20);
    ML_MRI_COUPLING(vcptr) = MRIC;
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vcptr);
}
 
#if 540 <= SUNDIALS_LIB_VERSION
// Create an MRIStep.Coupling.t from a C value (ownership transferred to OCaml)
static value sunml_arkode_mri_coupling_wrap(MRIStepCoupling MRIC)
{
    CAMLparam0();
    CAMLlocal5(vr, vcptr, vg, vgi, vc);
    int i, j;

    // wrap as bigarrays
    vc = caml_ba_alloc_dims(BIGARRAY_FLOAT | CAML_BA_MANAGED, 1,
			    MRIC->c, MRIC->stages);

    vg = caml_alloc_tuple(MRIC->nmat);
    for (i = 0; i < MRIC->nmat; i++) {
	vgi = caml_alloc_tuple(MRIC->stages);	
	Store_field(vg, i, vgi);
	for (j = 0; j < MRIC->stages; j++) {
	    Store_field(vgi, j,
		caml_ba_alloc_dims(BIGARRAY_FLOAT | CAML_BA_MANAGED, 1,
				   MRIC->G[i][j], MRIC->stages));
	}
    }

    vcptr = caml_alloc_final(1, &finalize_mri_coupling, 1, 20);
    ML_MRI_COUPLING(vcptr) = MRIC;

    // create the record
    vr = caml_alloc_tuple(RECORD_ARKODE_MRI_COUPLING_SIZE);
    Store_field(vr, RECORD_ARKODE_MRI_COUPLING_CPTR, vcptr);
    Store_field(vr, RECORD_ARKODE_MRI_COUPLING_NMAT, Val_int(MRIC->nmat));
    Store_field(vr, RECORD_ARKODE_MRI_COUPLING_STAGES, Val_int(MRIC->stages));
    Store_field(vr, RECORD_ARKODE_MRI_COUPLING_METHOD_ORDER, Val_int(MRIC->q));
    Store_field(vr, RECORD_ARKODE_MRI_COUPLING_EMBEDDING, Val_int(MRIC->p));
    Store_field(vr, RECORD_ARKODE_MRI_COUPLING_MATRICES, vg);
    Store_field(vr, RECORD_ARKODE_MRI_COUPLING_ABSCISSAE, vc);

    CAMLreturn(vr);
}
#endif

CAMLprim value sunml_arkode_mri_coupling_load_table(value vtable)
{
    CAMLparam1(vtable);
    CAMLlocal1(vr);
#if 540 <= SUNDIALS_LIB_VERSION
    int imethod;
    MRIStepCoupling MRIC = NULL;

    switch (Int_val(vtable)) {
    case VARIANT_ARKODE_MRI_COUPLING_MIS_KW3:
	imethod = MIS_KW3;
	break;

    case VARIANT_ARKODE_MRI_COUPLING_GARK_ERK45a:
	imethod = MRI_GARK_ERK45a;
	break;

    case VARIANT_ARKODE_MRI_COUPLING_GARK_IRK21a:
	imethod = MRI_GARK_IRK21a;
	break;

    case VARIANT_ARKODE_MRI_COUPLING_GARK_ESDIRK34a:
	imethod = MRI_GARK_ESDIRK34a;
	break;

    default:
	caml_invalid_argument("unexpected error in MRIStep.Coupling.load_table");
    }

    MRIC = MRIStepCoupling_LoadTable(imethod);
    if (MRIC == NULL) caml_failwith("MRIStepCoupling_LoadTable");

    vr = sunml_arkode_mri_coupling_wrap(MRIC);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_mri_coupling_copy(value vt)
{
    CAMLparam1(vt);
    CAMLlocal1(vr);
#if 540 <= SUNDIALS_LIB_VERSION
    MRIStepCoupling MRIC =
	ML_MRI_COUPLING(Field(vt, RECORD_ARKODE_MRI_COUPLING_CPTR));

    MRIC = MRIStepCoupling_Copy(MRIC);
    if (MRIC == NULL) caml_failwith("MRIStepCoupling_Copy");

    vr = sunml_arkode_mri_coupling_wrap(MRIC);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_mri_coupling_mistomri(value vbt, value vq, value vp)
{
    CAMLparam3(vbt, vq, vp);
    CAMLlocal1(vr);
#if 540 <= SUNDIALS_LIB_VERSION
    ARKodeButcherTable bt = butcher_table_val(vbt);

    MRIStepCoupling MRIC = MRIStepCoupling_MIStoMRI(bt, Int_val(vq), Int_val(vp));
    if (bt != NULL) ARKodeButcherTable_Free(bt);
    if (MRIC == NULL) caml_failwith("MRIStepCoupling_MIStoMRI");

    vr = sunml_arkode_mri_coupling_wrap(MRIC);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_mri_coupling_space(value vt)
{
    CAMLparam1(vt);
    CAMLlocal1(vr);
#if 540 <= SUNDIALS_LIB_VERSION
    sundials_ml_index liw, lrw;

    MRIStepCoupling_Space(
	ML_MRI_COUPLING(Field(vt, RECORD_ARKODE_MRI_COUPLING_CPTR)),
	&liw, &lrw);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_index(lrw));
    Store_field(vr, 1, Val_index(liw));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_mri_coupling_write(value vt, value vlog)
{
    CAMLparam2(vt, vlog);
#if 540 <= SUNDIALS_LIB_VERSION
    MRIStepCoupling_Write(
	ML_MRI_COUPLING(Field(vt, RECORD_ARKODE_MRI_COUPLING_CPTR)),
	ML_CFILE(vlog));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_arkode_mri_set_coupling(value varkode_mem, value vt)
{
    CAMLparam2(varkode_mem, vt);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetCoupling(
		 ARKODE_MEM_FROM_ML (varkode_mem),
		 ML_MRI_COUPLING(Field(vt, RECORD_ARKODE_MRI_COUPLING_CPTR)));
    CHECK_FLAG("MRIStepSetCoupling", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_arkode_mri_get_current_coupling(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vr);
#if 540 <= SUNDIALS_LIB_VERSION
    MRIStepCoupling MRIC = NULL;

    int flag = MRIStepGetCurrentCoupling(ARKODE_MEM_FROM_ML (varkode_mem), &MRIC);
    CHECK_FLAG("MRIStepGetCurrentCoupling", flag);

    MRIC = MRIStepCoupling_Copy(MRIC);
    if (MRIC == NULL) caml_failwith("MRIStepCoupling_Copy");

    vr = sunml_arkode_mri_coupling_wrap(MRIC);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_mri_write_coupling(value varkode_mem, value vlog)
{
    CAMLparam2(varkode_mem, vlog);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepWriteCoupling(ARKODE_MEM_FROM_ML (varkode_mem),
				    ML_CFILE(vlog));
    CHECK_FLAG("MRIStepWriteCoupling", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

/* MRIStep: Linear and Nonlinear Solvers */

CAMLprim value sunml_arkode_mri_get_lin_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);
#if 540 <= SUNDIALS_LIB_VERSION
    long int lenrwLS;
    long int leniwLS;

    int flag = MRIStepGetLinWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				      &lenrwLS, &leniwLS);
    CHECK_FLAG("MRIStepGetLinWorkSpace", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_arkode_mri_get_num_jac_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumJacEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumJacEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_num_lin_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumLinRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumLinRhsEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_set_jac_times(value vdata, value vhas_setup,
					      value vhas_times)
{
    CAMLparam3(vdata, vhas_setup, vhas_times);
#if 540 <= SUNDIALS_LIB_VERSION
    ARKLsJacTimesSetupFn setup = Bool_val (vhas_setup) ? jacsetupfn : NULL;
    ARKLsJacTimesVecFn   times = Bool_val (vhas_times) ? jactimesfn : NULL;

    int flag = MRIStepSetJacTimes(ARKODE_MEM_FROM_ML(vdata), setup, times);
    CHECK_LS_FLAG("MRIStepSetJacTimes", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_jac_times_rhsfn(value vdata,
						    value vhas_rhsfn)
{
    CAMLparam2(vdata, vhas_rhsfn);
#if 540 <= SUNDIALS_LIB_VERSION
    ARKRhsFn rhsfn = Bool_val (vhas_rhsfn) ? jactimesrhsfn : NULL;

    int flag = MRIStepSetJacTimesRhsFn(ARKODE_MEM_FROM_ML(vdata), rhsfn);
    CHECK_LS_FLAG("MRIStepSetJacTimesRhsFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_preconditioner (value vsession,
						    value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = ARKODE_MEM_FROM_ML (vsession);
#if 540 <= SUNDIALS_LIB_VERSION
    ARKLsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = MRIStepSetPreconditioner (mem, setup, precsolvefn);
    CHECK_FLAG ("MRIStepSetPreconditioner", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_jac_eval_frequency(value varkode_mem,
						       value vevalfreq)
{
    CAMLparam2(varkode_mem, vevalfreq);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetJacEvalFrequency(ARKODE_MEM_FROM_ML(varkode_mem),
					  Long_val(vevalfreq));
    CHECK_LS_FLAG("MRIStepSetJacEvalFrequency", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_linear_solution_scaling(value varkode_mem,
							    value vonoff)
{
    CAMLparam2(varkode_mem, vonoff);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetLinearSolutionScaling(ARKODE_MEM_FROM_ML(varkode_mem),
					       Bool_val(vonoff));
    CHECK_FLAG("MRIStepSetLinearSolutionScaling", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_eps_lin(value varkode_mem, value eplifac)
{
    CAMLparam2(varkode_mem, eplifac);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetEpsLin(ARKODE_MEM_FROM_ML(varkode_mem),
				Double_val(eplifac));
    CHECK_FLAG("MRIStepSetEpsLin", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_ls_norm_factor(value varkode_mem, value vfac)
{
    CAMLparam2(varkode_mem, vfac);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetLSNormFactor(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(vfac));
    CHECK_FLAG("MRIStepSetLSNormFactor", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_get_num_lin_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumLinIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumLinIters", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_num_lin_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumLinConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumLinConvFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_num_prec_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumPrecEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumPrecEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_num_prec_solves(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumPrecSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumPrecSolves", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_num_jtsetup_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumJTSetupEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumJTSetupEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_num_jtimes_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumJtimesEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumJtimesEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_set_linear(value varkode_mem, value vtimedepend)
{
    CAMLparam2(varkode_mem, vtimedepend);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetLinear(ARKODE_MEM_FROM_ML(varkode_mem),
			        Bool_val(vtimedepend));
    CHECK_FLAG("MRIStepSetLinear", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_nonlinear(value varkode_mem)
{
    CAMLparam1(varkode_mem);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetNonlinear(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("MRIStepSetNonlinear", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_stage_predict_fn(value varkode_mem,
						     value vset)
{
    CAMLparam2(varkode_mem, vset);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetStagePredictFn(ARKODE_MEM_FROM_ML(varkode_mem),
				    Bool_val(vset) ? stagepredictfn : NULL);
    CHECK_FLAG("MRIStepSetStagePredictFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * MRIStep Inner Stepper interface.
 */

#if 580 <= SUNDIALS_LIB_VERSION
static void finalize_istepper(value vistepper_cptr)
{
    MRIStepInnerStepper stepper = ISTEPPER(vistepper_cptr);

    value *pvcallbacks = NULL;
    MRIStepInnerStepper_GetContent(stepper, (void **)&pvcallbacks);
    if (pvcallbacks != NULL) sunml_sundials_free_value(pvcallbacks);

    MRIStepInnerStepper_Free(&stepper);
}

static void finalize_sundials_istepper(value vistepper_cptr)
{
    MRIStepInnerStepper stepper = ISTEPPER(vistepper_cptr);
    MRIStepInnerStepper_Free(&stepper);
}

static int istepper_evolvefn(MRIStepInnerStepper stepper,
			     realtype t0,
			     realtype tout,
			     N_Vector v)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal1(callbacks);

    value *pvcallbacks = NULL;
    MRIStepInnerStepper_GetContent(stepper, (void **)&pvcallbacks);

    args[0] = caml_copy_double(t0);
    args[1] = caml_copy_double(tout);
    args[2] = NVEC_BACKLINK(v);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(*pvcallbacks, RECORD_ARKODE_MRI_ISTEPPER_CALLBACKS_EVOLVE_FN),
	3, args);
    if (!Is_exception_result (r)) CAMLreturnT(int, 0);

    r = Extract_exception (r);
    CAMLreturnT(int, istepper_translate_exception (r, RECOVERABLE));
}

static int istepper_fullrhsfn(MRIStepInnerStepper stepper,
			      realtype t,
			      N_Vector v,
			      N_Vector f,
			      int mode)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal2(callbacks, vmode);

    value *pvcallbacks = NULL;
    MRIStepInnerStepper_GetContent(stepper, (void **)&pvcallbacks);

    switch (mode) {
    case ARK_FULLRHS_START:
	vmode = Val_int(VARIANT_ARKODE_MRI_ISTEPPER_FULL_RHSFN_MODE_START);
	break;

    case ARK_FULLRHS_END:
	vmode = Val_int(VARIANT_ARKODE_MRI_ISTEPPER_FULL_RHSFN_MODE_END);
	break;

    case ARK_FULLRHS_OTHER:
    default:
	vmode = Val_int(VARIANT_ARKODE_MRI_ISTEPPER_FULL_RHSFN_MODE_OTHER);
	break;
    }

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(v);
    args[2] = NVEC_BACKLINK(f);
    args[3] = vmode;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(*pvcallbacks, RECORD_ARKODE_MRI_ISTEPPER_CALLBACKS_FULL_RHS_FN),
	4, args);
    if (!Is_exception_result (r)) CAMLreturnT(int, 0);

    r = Extract_exception (r);
    CAMLreturnT(int, istepper_translate_exception (r, RECOVERABLE));
}

static int istepper_resetfn(MRIStepInnerStepper stepper,
			    realtype tR,
			    N_Vector vR)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal1(callbacks);

    value *pvcallbacks = NULL;
    MRIStepInnerStepper_GetContent(stepper, (void **)&pvcallbacks);

    args[0] = caml_copy_double(tR);
    args[1] = NVEC_BACKLINK(vR);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(*pvcallbacks, RECORD_ARKODE_MRI_ISTEPPER_CALLBACKS_RESET_FN),
	2, args);
    if (!Is_exception_result (r)) CAMLreturnT(int, 0);

    r = Extract_exception (r);
    CAMLreturnT(int, istepper_translate_exception (r, RECOVERABLE));
}
#endif

CAMLprim value sunml_arkode_mri_istepper_from_arkstep(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);
#if 580 <= SUNDIALS_LIB_VERSION
    int flag;
    MRIStepInnerStepper stepper;

    flag = ARKStepCreateMRIStepInnerStepper(ARKODE_MEM_FROM_ML(varkode_mem),
					    &stepper);
    CHECK_FLAG("ARKStepCreateMRIStepInnerStepper", flag);

    r = caml_alloc_final(1, &finalize_sundials_istepper, 0, 1);
    ISTEPPER(r) = stepper;
#else
    r = Val_unit;
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_arkode_mri_istepper_create(value vcallbacks,
						value vhasresetfn)
{
    CAMLparam2(vcallbacks, vhasresetfn);
    CAMLlocal1(vistepper);
#if 580 <= SUNDIALS_LIB_VERSION
    int flag;
    MRIStepInnerStepper stepper;

    flag = MRIStepInnerStepper_Create(&stepper);
    CHECK_FLAG("MRIStepInnerStepper_Create", flag);

    flag = MRIStepInnerStepper_SetEvolveFn(stepper, istepper_evolvefn);
    if (flag != ARK_SUCCESS) {
	MRIStepInnerStepper_Free(&stepper);
	sunml_arkode_check_flag("MRIStepInnerStepper_SetEvolveFn", flag, NULL);
    }

    flag = MRIStepInnerStepper_SetFullRhsFn(stepper, istepper_fullrhsfn);
    if (flag != ARK_SUCCESS) {
	MRIStepInnerStepper_Free(&stepper);
	sunml_arkode_check_flag("MRIStepInnerStepper_SetFullRhsFn", flag, NULL);
    }

    if (Bool_val(vhasresetfn))
    {
	flag = MRIStepInnerStepper_SetResetFn(stepper, istepper_resetfn);
	if (flag != ARK_SUCCESS) {
	    MRIStepInnerStepper_Free(&stepper);
	    sunml_arkode_check_flag("MRIStepInnerStepper_SetResetFn", flag, NULL);
	}
    }

    value *pvcallbacks = sunml_sundials_malloc_value(vcallbacks);
    if (pvcallbacks == NULL) caml_raise_out_of_memory();

    flag = MRIStepInnerStepper_SetContent(stepper, (void *)pvcallbacks);
    if (flag != ARK_SUCCESS) {
	MRIStepInnerStepper_Free(&stepper);
	sunml_sundials_free_value(pvcallbacks);
	sunml_arkode_check_flag("MRIStepInnerStepper_SetContent", flag, NULL);
    }

    vistepper = caml_alloc_final(1, &finalize_istepper, 0, 1);
    ISTEPPER(vistepper) = stepper;
#endif
    CAMLreturn(vistepper);
}

CAMLprim value sunml_arkode_mri_istepper_add_forcing(value visteppercptr,
						     value vt, value vff)
{
    CAMLparam3(visteppercptr, vt, vff);
#if 580 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepInnerStepper_AddForcing(ISTEPPER(visteppercptr),
		    Double_val(vt), NVEC_VAL(vff));
    CHECK_FLAG("MRIStepInnerStepper_AddForcing", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_arkode_mri_istepper_get_forcing_data(value visteppercptr)
{
    CAMLparam1(visteppercptr);
    CAMLlocal1(vr);
#if 580 <= SUNDIALS_LIB_VERSION
    realtype tshift;
    realtype tscale;
    N_Vector *forcing;
    int nforcing;

    int flag = MRIStepInnerStepper_GetForcingData(ISTEPPER(visteppercptr),
		    &tshift, &tscale, &forcing, &nforcing);
    CHECK_FLAG("MRIStepInnerStepper_GetForcingData", flag);

    if (forcing == NULL) caml_invalid_argument("no forcing data");

    vr = caml_alloc_tuple(RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_SIZE);
    Store_field(vr, RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_TSHIFT,
		caml_copy_double(tshift));
    Store_field(vr, RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_TSCALE,
		caml_copy_double(tscale));
    Store_field(vr, RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_FORCING,
		sunml_wrap_to_nvector_table(nforcing, forcing));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Boiler plate definitions for MRIStep interface.
 */

/* main solver optional output functions */

CAMLprim value sunml_arkode_mri_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);
#if 400 <= SUNDIALS_LIB_VERSION

    int flag;
    long int lenrw;
    long int leniw;

    flag = MRIStepGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem), &lenrw, &leniw);
    CHECK_FLAG("MRIStepGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_arkode_mri_get_num_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vr);

#if 500 <= SUNDIALS_LIB_VERSION
    long int ns;

    int flag = MRIStepGetNumSteps(ARKODE_MEM_FROM_ML(varkode_mem), &ns);
    CHECK_FLAG("MRIStepGetNumSteps", flag);

    vr = Val_long(ns);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vr);
}

CAMLprim value sunml_arkode_mri_get_last_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype v = 0.0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetLastStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("MRIStepGetLastStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_mri_get_num_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vr);

#if 500 <= SUNDIALS_LIB_VERSION
    long int ne;

    int flag = MRIStepGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &ne);
    CHECK_FLAG("MRIStepGetNumRhsEvals", flag);

    vr = Val_long(ne);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vr);
}

/* optional inputs for MRIStep */

CAMLprim value sunml_arkode_mri_set_defaults(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetDefaults(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("MRIStepSetDefaults", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_interpolant_type(value varkode_mem,
						     value vinterptype)
{
    CAMLparam2(varkode_mem, vinterptype);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetInterpolantType(ARKODE_MEM_FROM_ML(varkode_mem),
		ark_interpolant_types[Int_val(vinterptype)]);
    CHECK_FLAG("MRIStepSetInterpolantType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_interpolant_degree(value varkode_mem,
						       value vinterpdegree)
{
    CAMLparam2(varkode_mem, vinterpdegree);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetInterpolantDegree(ARKODE_MEM_FROM_ML(varkode_mem),
					   Int_val(vinterpdegree));
    CHECK_FLAG("MRIStepSetInterpolantDegree", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_diagnostics(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("MRIStepSetDiagnostics", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_clear_diagnostics(value vdata)
{
    CAMLparam1(vdata);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("MRIStepSetDiagnostics", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_error_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetErrFile(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("MRIStepSetErrFile", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), errh,
				      ARKODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("MRIStepSetErrHandlerFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("MRIStepSetErrHandlerFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_fixed_step(value varkode_mem, value vslow)
{
    CAMLparam2(varkode_mem, vslow);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetFixedStep(ARKODE_MEM_FROM_ML(varkode_mem),
				   Double_val(vslow));
    CHECK_FLAG("MRIStepSetFixedStep", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_max_hnil_warns(value varkode_mem, value mxhnil)
{
    CAMLparam2(varkode_mem, mxhnil);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetMaxHnilWarns(ARKODE_MEM_FROM_ML(varkode_mem),
				      Int_val(mxhnil));
    CHECK_FLAG("MRIStepSetMaxHnilWarns", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_max_num_steps(value varkode_mem, value mxsteps)
{
    CAMLparam2(varkode_mem, mxsteps);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetMaxNumSteps(ARKODE_MEM_FROM_ML(varkode_mem),
				     Long_val(mxsteps));
    CHECK_FLAG("MRIStepSetMaxNumSteps", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_stop_time(value varkode_mem, value tstop)
{
    CAMLparam2(varkode_mem, tstop);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetStopTime(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(tstop));
    CHECK_FLAG("MRIStepSetStopTime", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_pre_inner_fn(value varkode_mem, value vset)
{
    CAMLparam2(varkode_mem, vset);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetPreInnerFn(ARKODE_MEM_FROM_ML(varkode_mem),
				    Bool_val(vset) ? preinnerfn : NULL);
    CHECK_FLAG("MRIStepSetPreInnerFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_post_inner_fn(value varkode_mem, value vset)
{
    CAMLparam2(varkode_mem, vset);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetPostInnerFn(ARKODE_MEM_FROM_ML(varkode_mem),
			 	     Bool_val(vset) ? postinnerfn : NULL);
    CHECK_FLAG("MRIStepSetPostInnerFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_get_current_time(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype v = 0.0;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetCurrentTime(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("MRIStepGetCurrentTime", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_mri_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);
#if 400 <= SUNDIALS_LIB_VERSION
    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = MRIStepGetRootInfo(ARKODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("MRIStepGetRootInfo", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_get_num_g_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumGEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumGEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_current_state(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vnv);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector nv;
    int flag = MRIStepGetCurrentState(ARKODE_MEM_FROM_ML(varkode_mem), &nv);
    CHECK_FLAG("MRIStepGetCurrentState", flag);

    vnv = NVEC_BACKLINK(nv);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim value sunml_arkode_mri_get_nonlin_system_data(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(vnv);
#if 540 <= SUNDIALS_LIB_VERSION
    realtype tcur, gamma;
    N_Vector zpred, zi, Fi, sdata;
    void *user_data;

    int flag = MRIStepGetNonlinearSystemData(ARKODE_MEM_FROM_ML(varkode_mem),
		    &tcur, &zpred, &zi, &Fi, &gamma, &sdata, &user_data);
    CHECK_FLAG("MRIStepGetNonlinearSystemData", flag);

    vnv = caml_alloc_tuple(RECORD_ARKODE_NONLIN_SYSTEM_DATA_SIZE);
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_TCUR,
	    caml_copy_double(tcur));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_ZPRED,
	    NVEC_BACKLINK(zpred));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_ZI,
	    NVEC_BACKLINK(zi));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_FI,
	    NVEC_BACKLINK(Fi));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_GAMMA,
	    caml_copy_double(gamma));
    Store_field(vnv, RECORD_ARKODE_NONLIN_SYSTEM_DATA_SDATA,
	    NVEC_BACKLINK(sdata));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim value sunml_arkode_mri_compute_state(value varkode_mem,
					      value vzcor, value vz)
{
    CAMLparam3(varkode_mem, vzcor, vz);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepComputeState(ARKODE_MEM_FROM_ML(varkode_mem),
				   NVEC_VAL(vzcor), NVEC_VAL(vz));
    CHECK_FLAG("MRIStepComputeState", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim value sunml_arkode_mri_set_no_inactive_root_warn(value varkode_mem)
{
    CAMLparam1(varkode_mem);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetNoInactiveRootWarn(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("MRIStepSetNoInactiveRootWarn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_write_parameters(value varkode_mem, value vlog)
{
    CAMLparam2(varkode_mem, vlog);
#if 410 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepWriteParameters(ARKODE_MEM_FROM_ML(varkode_mem),
				      ML_CFILE(vlog));
    CHECK_FLAG("MRIStepWriteParameters", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_unit);
}

CAMLprim value sunml_arkode_mri_set_nonlin_conv_coef(value varkode_mem, value nlscoef)
{
    CAMLparam2(varkode_mem, nlscoef);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetNonlinConvCoef(ARKODE_MEM_FROM_ML(varkode_mem),
				        Double_val(nlscoef));
    CHECK_FLAG("MRIStepSetNonlinConvCoef", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_nonlin_crdown(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetNonlinCRDown(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("MRIStepSetNonlinCRDown", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_nonlin_rdiv(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetNonlinRDiv(ARKODE_MEM_FROM_ML(varkode_mem),
				    Double_val(varg));
    CHECK_FLAG("MRIStepSetNonlinRDiv", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_delta_gamma_max(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetDeltaGammaMax(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("MRIStepSetDeltaGammaMax", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_lsetup_frequency(value varkode_mem,
						     value varg)
{
    CAMLparam2(varkode_mem, varg);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetLSetupFrequency(ARKODE_MEM_FROM_ML(varkode_mem),
				         Int_val(varg));
    CHECK_FLAG("MRIStepSetLSetupFrequency", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_predictor_method(value varkode_mem, value vmethod)
{
    CAMLparam2(varkode_mem, vmethod);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetPredictorMethod(ARKODE_MEM_FROM_ML(varkode_mem),
				         Int_val(vmethod));
    CHECK_FLAG("MRIStepSetPredictorMethod", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_set_max_nonlin_iters(value varkode_mem, value maxcor)
{
    CAMLparam2(varkode_mem, maxcor);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepSetMaxNonlinIters(ARKODE_MEM_FROM_ML(varkode_mem),
				        Int_val(maxcor));
    CHECK_FLAG("MRIStepSetMaxNonlinIters", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_get_current_gamma(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    double gamma;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetCurrentGamma(ARKODE_MEM_FROM_ML(varkode_mem), &gamma);
    CHECK_FLAG("MRIStepGetCurrentGamma", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(caml_copy_double(gamma));
}

CAMLprim value sunml_arkode_mri_get_tol_scale_factor(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    realtype r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetTolScaleFactor(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetTolScaleFactor", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_arkode_mri_get_err_weights(value varkode_mem, value verrws)
{
    CAMLparam2(varkode_mem, verrws);
#if 540 <= SUNDIALS_LIB_VERSION
    N_Vector errws_nv = NVEC_VAL(verrws);
    int flag = MRIStepGetErrWeights(ARKODE_MEM_FROM_ML(varkode_mem), errws_nv);
    CHECK_FLAG("MRIStepGetErrWeights", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mri_get_num_lin_solv_setups(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int v;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumLinSolvSetups(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumLinSolvSetups", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_mri_get_num_nonlin_solv_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumNonlinSolvIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumNonlinSolvIters", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_num_nonlin_solv_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = MRIStepGetNumNonlinSolvConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("MRIStepGetNumNonlinSolvConvFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_mri_get_nonlin_solv_stats(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);
#if 540 <= SUNDIALS_LIB_VERSION
    long int nniters, nncfails;
    int flag = MRIStepGetNonlinSolvStats(ARKODE_MEM_FROM_ML(varkode_mem),
				         &nniters, &nncfails);
    CHECK_FLAG("MRIStepGetNonlinSolvStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nniters));
    Store_field(r, 1, Val_long(nncfails));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

