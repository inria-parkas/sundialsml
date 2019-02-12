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

/* linear solvers */
#include <arkode/arkode_direct.h>
#include <arkode/arkode_spils.h>
#include <arkode/arkode_bandpre.h>
#include <arkode/arkode_impl.h>

#if SUNDIALS_LIB_VERSION < 300
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
#include "../sundials/sundials_ml.h"
#include "arkode_ml.h"
#include "../nvectors/nvector_ml.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

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

CAMLprim value sunml_arkode_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = ARKodeSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), errh,
				     ARKODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("ARKodeSetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = ARKodeSetErrHandlerFn(ARKODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("ARKodeSetErrHandlerFn", flag);

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

static int irhsfn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_ARKODE_SESSION_IRHSFN),
				 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int erhsfn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_ARKODE_SESSION_ERHSFN),
				 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

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

#if SUNDIALS_LIB_VERSION >= 270
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

#if SUNDIALS_LIB_VERSION >= 300

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

#if SUNDIALS_LIB_VERSION >= 300
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

static int linit(ARKodeMem ark_mem)
{
    CAMLparam0();
    CAMLlocal2 (session, cb);

    WEAK_DEREF (session, *(value*)ark_mem->ark_user_data);

    cb = ARKODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_ALTERNATE_CALLBACKS_LINIT);
    cb = Some_val(cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, session);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

static int lsetup(ARKodeMem ark_mem, int convfail,
		  N_Vector ypred, N_Vector fpred, booleantype *jcurPtr,
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    WEAK_DEREF (session, *(value*)ark_mem->ark_user_data);

    switch (convfail) {
    case ARK_NO_FAILURES:
	convfail = VARIANT_ARKODE_ALTERNATE_NO_FAILURES;
	break;

    case ARK_FAIL_BAD_J:
	convfail = VARIANT_ARKODE_ALTERNATE_FAIL_BAD_J;
	break;

    case ARK_FAIL_OTHER:
	convfail = VARIANT_ARKODE_ALTERNATE_FAIL_OTHER;
	break;
    }
    args = caml_alloc_tuple (RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_SIZE);
    Store_field (args, RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_CONV_FAIL,
		 Val_int (convfail));
    Store_field (args, RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_Y,
		 NVEC_BACKLINK(ypred));
    Store_field (args, RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_RHS,
		 NVEC_BACKLINK(fpred));
    Store_field (args, RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_TMP,
		 sunml_arkode_make_triple_tmp(tmp1, tmp2, tmp3));

    cb = ARKODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_ALTERNATE_CALLBACKS_LSETUP);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(cb, session, args);

    /* Update jcurPtr; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcurPtr = Bool_val(r);
	CAMLreturnT (int, 0);
    }
    r = Extract_exception (r);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}


#if SUNDIALS_LIB_VERSION >= 300
static int lsolve(ARKodeMem ark_mem, N_Vector b, N_Vector ycur, N_Vector fcur)
#else
static int lsolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight, N_Vector ycur,
		  N_Vector fcur)
#endif
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_ARKODE_ALTERNATE_LSOLVE_ARGS_SIZE);
    Store_field (args, RECORD_ARKODE_ALTERNATE_LSOLVE_ARGS_Y,
		 NVEC_BACKLINK (ycur));
    Store_field (args, RECORD_ARKODE_ALTERNATE_LSOLVE_ARGS_RHS,
		 NVEC_BACKLINK (fcur));

    WEAK_DEREF (session, *(value*)ark_mem->ark_user_data);

    cb = Field (session, RECORD_ARKODE_SESSION_LS_CALLBACKS);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_ALTERNATE_CALLBACKS_LSOLVE);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (cb, session, args, NVEC_BACKLINK (b));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

CAMLprim value sunml_arkode_get_gamma(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    r = caml_alloc_small (2 * Double_wosize, Double_array_tag);
    Store_double_field (r, 0, arkode_mem->ark_gamma);
    Store_double_field (r, 1, arkode_mem->ark_gammap);
    CAMLreturn(r);
}

CAMLprim value sunml_arkode_set_alternate (value varkode_mem, value vhas_init,
				       value vhas_setup)
{
    CAMLparam3(varkode_mem, vhas_init, vhas_setup);
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    arkode_mem->ark_linit  = Bool_val(vhas_init)  ? linit : NULL;
    arkode_mem->ark_lsetup  = Bool_val(vhas_setup) ? lsetup : NULL;
#if SUNDIALS_LIB_VERSION < 300
    arkode_mem->ark_setupNonNull = Bool_val(vhas_setup);
#endif
    arkode_mem->ark_lsolve = lsolve;
    arkode_mem->ark_lmem   = NULL;
    arkode_mem->ark_msolve_type = 3; // custom

    CAMLreturn (Val_unit);
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_arkode_dls_dense (value varkode_mem,
				   value vneqs, value vset_jac)
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

CAMLprim value sunml_arkode_dls_lapack_dense (value varkode_mem, value vneqs,
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

CAMLprim value sunml_arkode_dls_set_linear_solver (value varkode_mem, value vlsolv,
					       value vjmat, value vhasjac)
{
    CAMLparam4(varkode_mem, vlsolv, vjmat, vhasjac);
#if SUNDIALS_LIB_VERSION >= 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = ARKDlsSetLinearSolver(arkode_mem, lsolv, jmat);
    CHECK_FLAG ("ARKDlsSetLinearSolver", flag);
    if (Bool_val (vhasjac)) {
	flag = ARKDlsSetJacFn(arkode_mem, jacfn);
	CHECK_FLAG("ARKDlsSetJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_linear_solver (value varkode_mem,
						 value vlsolv)
{
    CAMLparam2(varkode_mem, vlsolv);
#if SUNDIALS_LIB_VERSION >= 300
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

CAMLprim value sunml_arkode_spils_set_preconditioner (value vsession,
						  value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = ARKODE_MEM_FROM_ML (vsession);
    ARKSpilsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = ARKSpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_FLAG ("ARKSpilsSetPreconditioner", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_banded_preconditioner (value vsession,
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

CAMLprim value sunml_arkode_spils_set_jac_times(value vdata, value vhas_setup,
							 value vhas_times)
{
    CAMLparam3(vdata, vhas_setup, vhas_times);
#if SUNDIALS_LIB_VERSION >= 300
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

CAMLprim value sunml_arkode_wf_tolerances (value vdata)
{
    CAMLparam1(vdata);
 
    int flag = ARKodeWFtolerances(ARKODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("ARKodeWFtolerances", flag);

    CAMLreturn (Val_unit);
}

/* Mass matrix routines. */

#if SUNDIALS_LIB_VERSION >= 300

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

#if SUNDIALS_LIB_VERSION >= 300
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

static int minit(ARKodeMem ark_mem)
{
    CAMLparam0();
    CAMLlocal2 (session, cb);

    WEAK_DEREF (session, *(value*)ark_mem->ark_user_data);

    cb = ARKODE_MASS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_ALTERNATE_MASS_CALLBACKS_MINIT);
    cb = Some_val(cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, session);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

static int msetup(ARKodeMem ark_mem, 
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)ark_mem->ark_user_data);

    cb = ARKODE_MASS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_ALTERNATE_MASS_CALLBACKS_MSETUP);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(cb, session,
				 sunml_arkode_make_triple_tmp(tmp1, tmp2, tmp3));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 300
static int msolve(ARKodeMem ark_mem, N_Vector b)
#else
static int msolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight)
#endif
{
    CAMLparam0();
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)ark_mem->ark_user_data);

    cb = Field (session, RECORD_ARKODE_SESSION_MASS_CALLBACKS);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_ALTERNATE_MASS_CALLBACKS_MSOLVE);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (cb, session, NVEC_BACKLINK (b));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

CAMLprim value sunml_arkode_set_mass_alternate (value varkode_mem, value vhas_init,
					    value vhas_setup)
{
    CAMLparam3(varkode_mem, vhas_init, vhas_setup);
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    arkode_mem->ark_minit  = Bool_val(vhas_init)  ? minit : NULL;
    arkode_mem->ark_msetup  = Bool_val(vhas_setup) ? msetup : NULL;
#if SUNDIALS_LIB_VERSION < 300
    arkode_mem->ark_MassSetupNonNull = Bool_val(vhas_setup);
#endif
    arkode_mem->ark_msolve = msolve;
    arkode_mem->ark_mass_mem  = NULL;
    arkode_mem->ark_msolve_type = 4; // custom

    CAMLreturn (Val_unit);
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

CAMLprim value sunml_arkode_dls_set_mass_linear_solver (value varkode_mem,
				    value vlsolv, value vjmat, value vtime_dep)
{
    CAMLparam4(varkode_mem, vlsolv, vjmat, vtime_dep);
#if SUNDIALS_LIB_VERSION >= 300
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
#if SUNDIALS_LIB_VERSION >= 300
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

CAMLprim value sunml_arkode_spils_set_mass_preconditioner (value vsession,
						       value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = ARKODE_MEM_FROM_ML (vsession);
    ARKSpilsMassPrecSetupFn setup
	= Bool_val (vset_precsetup) ? massprecsetupfn : NULL;
    int flag = ARKSpilsSetMassPreconditioner (mem, setup, massprecsolvefn);
    CHECK_FLAG ("ARKSpilsSetMassPreconditioner", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_mass_times(value varkode_mem,
					     value vhassetup)
{
    CAMLparam2(varkode_mem, vhassetup);
#if SUNDIALS_LIB_VERSION >= 300
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

/* basic interface */

/* ARKodeCreate() + ARKodeInit().  */
CAMLprim value sunml_arkode_init(value weakref, value hasfi, value hasfe,
			     value y0, value t0)
{
    CAMLparam5(weakref, hasfi, hasfe, y0, t0);
    CAMLlocal2(r, varkode_mem);

    int flag;

    void *arkode_mem = ARKodeCreate();
    if (arkode_mem == NULL)
	caml_failwith("ARKodeCreate returned NULL");

    varkode_mem = caml_alloc_final(1, NULL, 1, 5);
    ARKODE_MEM(varkode_mem) = arkode_mem;

    N_Vector nv_y0 = NVEC_VAL(y0);
    flag = ARKodeInit(arkode_mem,
		      Bool_val(hasfe)  ? erhsfn : NULL,
		      Bool_val(hasfi)  ? irhsfn : NULL,
		      Double_val(t0),
		      nv_y0);
    if (flag != ARK_SUCCESS) {
	ARKodeFree (&arkode_mem);
	CHECK_FLAG("ARKodeInit", flag);
    }

    value *backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	ARKodeFree (&arkode_mem);
	caml_raise_out_of_memory();
    }
    ARKodeSetUserData (arkode_mem, backref);

    r = caml_alloc_tuple (2);
    Store_field (r, 0, varkode_mem);
    Store_field (r, 1, (value)backref);

    CAMLreturn(r);
}

/* Set the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim value sunml_arkode_root_init (value vdata, value vnroots)
{
    CAMLparam2 (vdata, vnroots);
    void *arkode_mem = ARKODE_MEM_FROM_ML (vdata);
    int nroots = Int_val (vnroots);
    int flag = ARKodeRootInit (arkode_mem, nroots, roots);
    CHECK_FLAG ("ARKodeRootInit", flag);
    Store_field (vdata, RECORD_ARKODE_SESSION_NROOTS, vnroots);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_sv_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = ARKodeSVtolerances(ARKODE_MEM_FROM_ML(vdata),
				  Double_val(reltol), atol_nv);
    CHECK_FLAG("ARKodeSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_reinit(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    ARKRhsFn fe = NULL;
    ARKRhsFn fi = NULL;

    switch (ARKODE_PROBLEM(vdata)) {
    case VARIANT_ARKODE_PROBLEM_TYPE_IMPLICIT_ONLY:
	fi = irhsfn;
	break;
    case VARIANT_ARKODE_PROBLEM_TYPE_EXPLICIT_ONLY:
	fe = erhsfn;
	break;
    case VARIANT_ARKODE_PROBLEM_TYPE_IMPLICIT_AND_EXPLICIT:
	fe = erhsfn;
	fi = irhsfn;
	break;
    }

    N_Vector y0_nv = NVEC_VAL(y0);
    int flag = ARKodeReInit(ARKODE_MEM_FROM_ML(vdata),
			    fe,
			    fi,
			    Double_val(t0),
			    y0_nv);
    CHECK_FLAG("ARKodeReInit", flag);

    CAMLreturn (Val_unit);
}

static value solver(value vdata, value nextt, value vy, int onestep)
{
    CAMLparam3(vdata, nextt, vy);
    CAMLlocal1(ret);
    realtype tret;
    int flag;
    N_Vector y;
    enum arkode_solver_result_tag result = -1;

    y = NVEC_VAL (vy);
    // Caml_ba_data_val(y) must not be shifted by the OCaml GC during this
    // function call, which calls Caml through the callback f.  Is this
    // guaranteed?
    flag = ARKode (ARKODE_MEM_FROM_ML (vdata), Double_val (nextt), y, &tret,
		  onestep ? ARK_ONE_STEP : ARK_NORMAL);

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
	CHECK_FLAG ("ARKode", flag);
    }

    assert (Field (vdata, RECORD_ARKODE_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}

CAMLprim value sunml_arkode_solve_normal(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value sunml_arkode_solve_one_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 1));
}

CAMLprim value sunml_arkode_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVEC_VAL(vy);

    int flag = ARKodeGetDky(ARKODE_MEM_FROM_ML(vdata), Double_val(vt),
			    Int_val(vk), y_nv);
    CHECK_FLAG("ARKodeGetDky", flag);
    
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_get_err_weights(value varkode_mem, value verrws)
{
    CAMLparam2(varkode_mem, verrws);

    N_Vector errws_nv = NVEC_VAL(verrws);

    int flag = ARKodeGetErrWeights(ARKODE_MEM_FROM_ML(varkode_mem), errws_nv);
    CHECK_FLAG("ARKodeGetErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_get_est_local_errors(value varkode_mem, value vele)
{
    CAMLparam2(varkode_mem, vele);

    N_Vector ele_nv = NVEC_VAL(vele);

    int flag = ARKodeGetEstLocalErrors(ARKODE_MEM_FROM_ML(varkode_mem), ele_nv);
    CHECK_FLAG("ARKodeGetEstLocalErrors", flag);

    CAMLreturn (Val_unit);
}


void sunml_arkode_check_flag(const char *call, int flag)
{
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
	    caml_raise_constant(ARKODE_EXN(LinearSetupFailure));

	case ARK_LSOLVE_FAIL:
	    caml_raise_constant(ARKODE_EXN(LinearSolveFailure));

	case ARK_MASSINIT_FAIL:
	    caml_raise_constant(ARKODE_EXN(MassInitFailure));

	case ARK_MASSSETUP_FAIL:
	    caml_raise_constant(ARKODE_EXN(MassSetupFailure));

	case ARK_MASSSOLVE_FAIL:
	    caml_raise_constant(ARKODE_EXN(MassSolveFailure));

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

#if SUNDIALS_LIB_VERSION >= 270
	case ARK_POSTPROCESS_FAIL:
	    caml_raise_constant(ARKODE_EXN(PostprocStepFailure));
#endif

	case ARK_BAD_K:
	    caml_raise_constant(ARKODE_EXN(BadK));

	case ARK_BAD_T:
	    caml_raise_constant(ARKODE_EXN(BadT));

	default:
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    ARKodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

#if SUNDIALS_LIB_VERSION >= 400
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
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    ARKDlsGetReturnFlagName(flag));
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

#if SUNDIALS_LIB_VERSION >= 300
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

#if SUNDIALS_LIB_VERSION >+ 300
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

CAMLprim value sunml_arkode_session_finalize(value vdata)
{
    if (ARKODE_MEM_FROM_ML(vdata) != NULL) {
	void *arkode_mem = ARKODE_MEM_FROM_ML(vdata);
	value *backref = ARKODE_BACKREF_FROM_ML(vdata);
	ARKodeFree(&arkode_mem);
	sunml_sundials_free_value(backref);
    }

    return Val_unit;
}

CAMLprim value sunml_arkode_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = ARKodeSStolerances(ARKODE_MEM_FROM_ML(vdata),
				  Double_val(reltol),
				  Double_val(abstol));
    CHECK_FLAG("ARKodeSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = ARKodeGetRootInfo(ARKODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("ARKodeGetRootInfo", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_get_integrator_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;

    long int nsteps;
    long int expsteps;
    long int accsteps;
    long int step_attempt;
    long int nfe_evals;
    long int nfi_evals;
    long int nlinsetups;
    long int netfails;

    realtype hinused;
    realtype hlast;
    realtype hcur;
    realtype tcur;

    flag = ARKodeGetIntegratorStats(ARKODE_MEM_FROM_ML(vdata),
				    &nsteps,
				    &expsteps,
				    &accsteps,
				    &step_attempt,
				    &nfe_evals,
				    &nfi_evals,
				    &nlinsetups,
				    &netfails,
				    &hinused,
				    &hlast,
				    &hcur,
				    &tcur);
    CHECK_FLAG("ARKodeGetIntegratorStats", flag);

    r = caml_alloc_tuple(RECORD_ARKODE_INTEGRATOR_STATS_SIZE);
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_STEPS,Val_long(nsteps));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_EXP_STEPS,
							Val_long(expsteps));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_ACC_STEPS,
							Val_long(accsteps));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_STEP_ATTEMPTS,
							Val_long(step_attempt));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_NUM_NFE_EVALS,
							Val_long(nfe_evals));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_NUM_NFI_EVALS,
							Val_long(nfi_evals));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_LIN_SOLV_SETUPS,
							Val_long(nlinsetups));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_NUM_ERR_TEST_FAILS,
							Val_long(netfails));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_ACTUAL_INIT_STEP,
						    caml_copy_double(hinused));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_LAST_STEP,
						    caml_copy_double(hlast));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_CURRENT_STEP,
						    caml_copy_double(hcur));
    Store_field(r, RECORD_ARKODE_INTEGRATOR_STATS_CURRENT_TIME,
						    caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_set_error_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

    int flag = ARKodeSetErrFile(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ARKodeSetErrFile", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_diagnostics(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

    int flag = ARKodeSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("ARKodeSetDiagnostics", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_clear_diagnostics(value vdata)
{
    CAMLparam1(vdata);

    int flag = ARKodeSetDiagnostics(ARKODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("ARKodeSetDiagnostics", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < ARKODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = ARKodeSetRootDirection(ARKODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("ARKodeSetRootDirection", flag);

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

CAMLprim value sunml_arkode_resize(value varkode_mem,
			       value vhasfn, value vhscale,
			       value vt0, value vynew)
{
    CAMLparam5(varkode_mem, vhasfn, vhscale, vt0, vynew);
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKodeResize(
		    arkode_mem,
		    NVEC_VAL(vynew),
		    Double_val(vhscale),
		    Double_val(vt0),
		    Bool_val(vhasfn) ? resizefn : NULL,
		    Bool_val(vhasfn) ? arkode_mem->ark_user_data : NULL);
    CHECK_FLAG("ARKodeResize", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_get_current_butcher_tables(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal4(ai, ae, ci, ce);
    CAMLlocal4(bi, be, b2i, b2e);
    CAMLlocal4(rkm, rtci, rtce, r);

    int flag;
    int s;
    int q;
    int p;

    ai  = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX * ARK_S_MAX);
    ae  = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX * ARK_S_MAX);
    ci  = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    bi  = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    b2i = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);

#if SUNDIALS_LIB_VERSION >= 270
    ce  = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    be  = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    b2e = caml_ba_alloc_dims (BIGARRAY_FLOAT, 1, NULL, ARK_S_MAX);
    flag = ARKodeGetCurrentButcherTables(ARKODE_MEM_FROM_ML(varkode_mem),
					 &s, &q, &p,
					 REAL_ARRAY(ai),
					 REAL_ARRAY(ae),
					 REAL_ARRAY(ci),
					 REAL_ARRAY(ce),
					 REAL_ARRAY(bi),
					 REAL_ARRAY(be),
					 REAL_ARRAY(b2i),
					 REAL_ARRAY(b2e));
    CHECK_FLAG("ARKodeGetCurrentButcherTables", flag);

    rtci = caml_alloc_tuple(RECORD_ARKODE_RK_TIMESCOEFS_SIZE);
    Store_field(rtci, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES, ci);
    Store_field(rtci, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS, bi);
    Store_some(Field(rtci, RECORD_ARKODE_RK_TIMESCOEFS_BEMBED), b2i);

    rtce = caml_alloc_tuple(RECORD_ARKODE_RK_TIMESCOEFS_SIZE);
    Store_field(rtce, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES, ce);
    Store_field(rtce, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS, be);
    Store_some(Field(rtce, RECORD_ARKODE_RK_TIMESCOEFS_BEMBED), b2e);
#else
    flag = ARKodeGetCurrentButcherTables(ARKODE_MEM_FROM_ML(varkode_mem),
					 &s, &q, &p,
					 REAL_ARRAY(ai),
					 REAL_ARRAY(ae),
					 REAL_ARRAY(ci),
					 REAL_ARRAY(bi),
					 REAL_ARRAY(b2i));
    CHECK_FLAG("ARKodeGetCurrentButcherTables", flag);

    rtci = caml_alloc_tuple(RECORD_ARKODE_RK_TIMESCOEFS_SIZE);
    Store_field(rtci, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES, ci);
    Store_field(rtci, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS, bi);
    Store_some(Field(rtci, RECORD_ARKODE_RK_TIMESCOEFS_BEMBED), b2i);

    rtce = rtci;
#endif

    rkm = caml_alloc_tuple(RECORD_ARKODE_RK_METHOD_SIZE);
    Store_field(rkm, RECORD_ARKODE_RK_METHOD_STAGES, Val_int(s));
    Store_field(rkm, RECORD_ARKODE_RK_METHOD_GLOBAL_ORDER, Val_int(q));
    Store_field(rkm, RECORD_ARKODE_RK_METHOD_GLOBAL_EMBEDDED_ORDER,Val_int(p));

    r = caml_alloc_tuple(5);
    Store_field(r, 0, rkm);
    Store_field(r, 1, ai);
    Store_field(r, 2, ae);
    Store_field(r, 3, rtci);
    Store_field(r, 4, rtce);

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_set_erk_table(value varkode_mem,
				      value vrkm, value vae, value vtc)
{
    CAMLparam4(varkode_mem, vrkm, vae, vtc);
    CAMLlocal1(vtcb2);

    vtcb2 = Field(vtc, RECORD_ARKODE_RK_TIMESCOEFS_BEMBED);

    int flag =
	ARKodeSetERKTable(ARKODE_MEM_FROM_ML(varkode_mem),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_STAGES)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_ORDER)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_EMBEDDED_ORDER)),
	    REAL_ARRAY(Field(vtc, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES)),
	    REAL_ARRAY(vae),
	    REAL_ARRAY(Field(vtc, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS)),
	    ((vtcb2 == Val_none) ? NULL : REAL_ARRAY(Some_val(vtcb2))));
    CHECK_FLAG("ARKodeSetERKTable", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_irk_table(value varkode_mem,
				      value vrkm, value vai, value vtc)
{
    CAMLparam4(varkode_mem, vrkm, vai, vtc);
    CAMLlocal1(vtcb2);

    vtcb2 = Field(vtc, RECORD_ARKODE_RK_TIMESCOEFS_BEMBED);

    int flag =
	ARKodeSetIRKTable(ARKODE_MEM_FROM_ML(varkode_mem),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_STAGES)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_ORDER)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_EMBEDDED_ORDER)),
	    REAL_ARRAY(Field(vtc, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES)),
	    REAL_ARRAY(vai),
	    REAL_ARRAY(Field(vtc, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS)),
	    ((vtcb2 == Val_none) ? NULL : REAL_ARRAY(Some_val(vtcb2))));
    CHECK_FLAG("ARKodeSetIRKTable", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_ark_tables(value varkode_mem,
				       value vrkm, value vai, value vae,
				       value vtcs)
{
    CAMLparam5(varkode_mem, vrkm, vai, vae, vtcs);
    CAMLlocal4(vtcsi, vtcse, vb2i, vb2e);
    int flag;

    vtcsi = Field(vtcs, 0);
    vtcse = Field(vtcs, 1);
    vb2i = Field(vtcsi, RECORD_ARKODE_RK_TIMESCOEFS_BEMBED);
    vb2e = Field(vtcse, RECORD_ARKODE_RK_TIMESCOEFS_BEMBED);

#if SUNDIALS_LIB_VERSION >= 270
    flag = ARKodeSetARKTables(ARKODE_MEM_FROM_ML(varkode_mem),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_STAGES)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_ORDER)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_EMBEDDED_ORDER)),
	    REAL_ARRAY(Field(vtcsi, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES)),
	    REAL_ARRAY(Field(vtcse, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES)),
	    REAL_ARRAY(vai),
	    REAL_ARRAY(vae),
	    REAL_ARRAY(Field(vtcsi, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS)),
	    REAL_ARRAY(Field(vtcse, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS)),
	    ((vb2i == Val_none) ? NULL : REAL_ARRAY(Some_val(vb2i))),
	    ((vb2e == Val_none) ? NULL : REAL_ARRAY(Some_val(vb2e))));
#else
    flag = ARKodeSetARKTables(ARKODE_MEM_FROM_ML(varkode_mem),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_STAGES)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_ORDER)),
	    Int_val(Field(vrkm, RECORD_ARKODE_RK_METHOD_GLOBAL_EMBEDDED_ORDER)),
	    REAL_ARRAY(Field(vtcsi, RECORD_ARKODE_RK_TIMESCOEFS_STAGE_TIMES)),
	    REAL_ARRAY(vai),
	    REAL_ARRAY(vae),
	    REAL_ARRAY(Field(vtcsi, RECORD_ARKODE_RK_TIMESCOEFS_COEFFICIENTS)),
	    REAL_ARRAY(Some_val(vb2i)));
#endif
    CHECK_FLAG("ARKodeSetARKTables", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_erk_table_num(value varkode_mem, value vnum)
{
    CAMLparam2(varkode_mem, vnum);

    int flag = ARKodeSetERKTableNum(ARKODE_MEM_FROM_ML(varkode_mem),
				    Int_val(vnum));
    CHECK_FLAG("ARKodeSetERKTableNum", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_irk_table_num(value varkode_mem, value vnum)
{
    CAMLparam2(varkode_mem, vnum);

    int flag = ARKodeSetIRKTableNum(ARKODE_MEM_FROM_ML(varkode_mem),
				    Int_val(vnum));
    CHECK_FLAG("ARKodeSetIRKTableNum", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_ark_table_num(value varkode_mem, value vnums)
{
    CAMLparam2(varkode_mem, vnums);

    int flag = ARKodeSetARKTableNum(ARKODE_MEM_FROM_ML(varkode_mem),
				    Int_val(Field(vnums, 0)),
				    Int_val(Field(vnums, 1)));
    CHECK_FLAG("ARKodeSetARKTableNum", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_adaptivity_method(value varkode_mem, value vmeth)
{
    CAMLparam2(varkode_mem, vmeth);
    CAMLlocal2(vks, vorder);

    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    if (Tag_val(vmeth) == VARIANT_ARKODE_ADAPTIVITY_METHOD_ADAPTIVITYFN) {
	flag = ARKodeSetAdaptivityFn(arkode_mem,
				     adaptfn,
				     arkode_mem->ark_user_data);
	CHECK_FLAG("ARKodeSetAdaptivityFn", flag);

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

	flag = ARKodeSetAdaptivityMethod(arkode_mem, 
					 Tag_val(vmeth), 
					 vks == Val_none,
					 Bool_val(vorder),
					 vks != Val_none ? adapt_params : NULL);
	CHECK_FLAG("ARKodeSetAdaptivityMethod", flag);
    }

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_stability_fn(value varkode_mem, value vhasf)
{
    CAMLparam2(varkode_mem, vhasf);
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKodeSetStabilityFn(arkode_mem,
				    Bool_val(vhasf) ? stabfn : NULL,
				    arkode_mem->ark_user_data);
    CHECK_FLAG("ARKodeSetStabilityFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_imex(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetImEx(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetImEx", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_implicit(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetImplicit(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetImplicit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_explicit(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetExplicit(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetExplicit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_fixed_point(value varkode_mem, value vfpm)
{
    CAMLparam2(varkode_mem, vfpm);

    int flag = ARKodeSetFixedPoint(ARKODE_MEM_FROM_ML(varkode_mem),
				   Long_val(vfpm));
    CHECK_FLAG("ARKodeSetFixedPoint", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_newton(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetNewton(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetNewton", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_linear(value varkode_mem, value vtimedepend)
{
    CAMLparam2(varkode_mem, vtimedepend);

    int flag = ARKodeSetLinear(ARKODE_MEM_FROM_ML(varkode_mem),
			       Bool_val(vtimedepend));
    CHECK_FLAG("ARKodeSetLinear", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_nonlinear(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetNonlinear(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetNonlinear", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_predictor_method(value varkode_mem, value vmethod)
{
    CAMLparam2(varkode_mem, vmethod);

    int flag = ARKodeSetPredictorMethod(ARKODE_MEM_FROM_ML(varkode_mem),
				        Int_val(vmethod));
    CHECK_FLAG("ARKodeSetPredictorMethod", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_postprocess_step_fn(value varkode_mem, value vhasf)
{
    CAMLparam2(varkode_mem, vhasf);
#if SUNDIALS_LIB_VERSION >= 270
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKodeSetPostprocessStepFn(arkode_mem,
					  Bool_val(vhasf) ? poststepfn : NULL);
    CHECK_FLAG("ARKodeSetPostprocessStepFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Boiler plate definitions for Arkode interface.
 */

CAMLprim value sunml_arkode_resv_tolerance(value vdata, value abstol)
{
    CAMLparam2(vdata, abstol);

    int flag = ARKodeResVtolerance(ARKODE_MEM_FROM_ML(vdata), NVEC_VAL(abstol));
    CHECK_FLAG("ARKodeResVtolerance", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_ress_tolerance(value vdata, value abstol)
{
    CAMLparam2(vdata, abstol);

    int flag = ARKodeResStolerance(ARKODE_MEM_FROM_ML(vdata),
				   Double_val(abstol));
    CHECK_FLAG("ARKodeResStolerance", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_resf_tolerance(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = ARKodeResFtolerance(ARKODE_MEM_FROM_ML(vdata), resw);
    CHECK_FLAG("ARKodeResFtolerance", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = ARKodeGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem), &lenrw, &leniw);
    CHECK_FLAG("ARKodeGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_get_num_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

    flag = ARKodeGetNumSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_get_num_acc_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

    flag = ARKodeGetNumAccSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumAccSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_get_num_exp_steps(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

    flag = ARKodeGetNumExpSteps(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumExpSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_get_num_step_attempts(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

    flag = ARKodeGetNumStepAttempts(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumStepAttempts", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_get_num_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int nfe;
    long int nfi;

    flag = ARKodeGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &nfe, &nfi);
    CHECK_FLAG("ARKodeGetNumRhsEvals", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nfe));
    Store_field(r, 1, Val_long(nfi));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_get_num_lin_solv_setups(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

    flag = ARKodeGetNumLinSolvSetups(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_get_num_err_test_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

    flag = ARKodeGetNumErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_arkode_get_actual_init_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

    flag = ARKodeGetActualInitStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_get_last_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

    flag = ARKodeGetLastStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetLastStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_get_current_step(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

    flag = ARKodeGetCurrentStep(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_get_current_time(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    realtype v;

    flag = ARKodeGetCurrentTime(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKodeGetCurrentTime", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_arkode_set_max_num_steps(value varkode_mem, value mxsteps)
{
    CAMLparam2(varkode_mem, mxsteps);

    int flag = ARKodeSetMaxNumSteps(ARKODE_MEM_FROM_ML(varkode_mem),
				    Long_val(mxsteps));
    CHECK_FLAG("ARKodeSetMaxNumSteps", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_hnil_warns(value varkode_mem, value mxhnil)
{
    CAMLparam2(varkode_mem, mxhnil);

    int flag = ARKodeSetMaxHnilWarns(ARKODE_MEM_FROM_ML(varkode_mem),
				     Int_val(mxhnil));
    CHECK_FLAG("ARKodeSetMaxHnilWarns", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_init_step(value varkode_mem, value hin)
{
    CAMLparam2(varkode_mem, hin);

    int flag = ARKodeSetInitStep(ARKODE_MEM_FROM_ML(varkode_mem),
				 Double_val(hin));
    CHECK_FLAG("ARKodeSetInitStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_min_step(value varkode_mem, value hmin)
{
    CAMLparam2(varkode_mem, hmin);

    int flag = ARKodeSetMinStep(ARKODE_MEM_FROM_ML(varkode_mem),
				Double_val(hmin));
    CHECK_FLAG("ARKodeSetMinStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_step(value varkode_mem, value hmax)
{
    CAMLparam2(varkode_mem, hmax);

    int flag = ARKodeSetMaxStep(ARKODE_MEM_FROM_ML(varkode_mem),
				Double_val(hmax));
    CHECK_FLAG("ARKodeSetMaxStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_stop_time(value varkode_mem, value tstop)
{
    CAMLparam2(varkode_mem, tstop);

    int flag = ARKodeSetStopTime(ARKODE_MEM_FROM_ML(varkode_mem),
				 Double_val(tstop));
    CHECK_FLAG("ARKodeSetStopTime", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_err_test_fails(value varkode_mem, value maxnef)
{
    CAMLparam2(varkode_mem, maxnef);

    int flag = ARKodeSetMaxErrTestFails(ARKODE_MEM_FROM_ML(varkode_mem),
					Int_val(maxnef));
    CHECK_FLAG("ARKodeSetMaxErrTestFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_nonlin_iters(value varkode_mem, value maxcor)
{
    CAMLparam2(varkode_mem, maxcor);

    int flag = ARKodeSetMaxNonlinIters(ARKODE_MEM_FROM_ML(varkode_mem),
				       Int_val(maxcor));
    CHECK_FLAG("ARKodeSetMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_conv_fails(value varkode_mem, value maxncf)
{
    CAMLparam2(varkode_mem, maxncf);

    int flag = ARKodeSetMaxConvFails(ARKODE_MEM_FROM_ML(varkode_mem),
				     Int_val(maxncf));
    CHECK_FLAG("ARKodeSetMaxConvFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_nonlin_conv_coef(value varkode_mem, value nlscoef)
{
    CAMLparam2(varkode_mem, nlscoef);

    int flag = ARKodeSetNonlinConvCoef(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(nlscoef));
    CHECK_FLAG("ARKodeSetNonlinConvCoef", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_no_inactive_root_warn(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetNoInactiveRootWarn(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetNoInactiveRootWarn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_cfl_fraction(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetCFLFraction(ARKODE_MEM_FROM_ML(varkode_mem),
				    Double_val(varg));
    CHECK_FLAG("ARKodeSetCFLFraction", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_defaults(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetDefaults(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetDefaults", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_delta_gamma_max(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetDeltaGammaMax(ARKODE_MEM_FROM_ML(varkode_mem),
				      Double_val(varg));
    CHECK_FLAG("ARKodeSetDeltaGammaMax", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_dense_order(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetDenseOrder(ARKODE_MEM_FROM_ML(varkode_mem),
				   Int_val(varg));
    CHECK_FLAG("ARKodeSetDenseOrder", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_error_bias(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetErrorBias(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKodeSetErrorBias", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_fixed_step(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetFixedStep(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKodeSetFixedStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_fixed_step_bounds(value varkode_mem,
					      value vlb, value vub)
{
    CAMLparam3(varkode_mem, vlb, vub);

    int flag = ARKodeSetFixedStepBounds(ARKODE_MEM_FROM_ML(varkode_mem),
				        Double_val(vlb), Double_val(vub));
    CHECK_FLAG("ARKodeSetFixedStepBounds", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_cfail_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetMaxCFailGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxCFailGrowth", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_efail_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetMaxEFailGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxEFailGrowth", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_first_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetMaxFirstGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				       Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxFirstGrowth", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_growth(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetMaxGrowth(ARKODE_MEM_FROM_ML(varkode_mem),
				  Double_val(varg));
    CHECK_FLAG("ARKodeSetMaxGrowth", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_max_steps_between_lset(value varkode_mem,
						   value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetMaxStepsBetweenLSet(ARKODE_MEM_FROM_ML(varkode_mem),
				            Int_val(varg));
    CHECK_FLAG("ARKodeSetMaxStepsBetweenLSet", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_nonlin_crdown(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetNonlinCRDown(ARKODE_MEM_FROM_ML(varkode_mem),
				     Double_val(varg));
    CHECK_FLAG("ARKodeSetNonlinCRDown", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_nonlin_rdiv(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetNonlinRDiv(ARKODE_MEM_FROM_ML(varkode_mem),
				   Double_val(varg));
    CHECK_FLAG("ARKodeSetNonlinRDiv", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_optimal_params(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag = ARKodeSetOptimalParams(ARKODE_MEM_FROM_ML(varkode_mem));
    CHECK_FLAG("ARKodeSetOptimalParams", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_order(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetOrder(ARKODE_MEM_FROM_ML(varkode_mem), Int_val(varg));
    CHECK_FLAG("ARKodeSetOrder", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_safety_factor(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetSafetyFactor(ARKODE_MEM_FROM_ML(varkode_mem),
				     Double_val(varg));
    CHECK_FLAG("ARKodeSetSafetyFactor", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_set_small_num_efails(value varkode_mem, value varg)
{
    CAMLparam2(varkode_mem, varg);

    int flag = ARKodeSetSmallNumEFails(ARKODE_MEM_FROM_ML(varkode_mem),
				       Int_val(varg));
    CHECK_FLAG("ARKodeSetSmallNumEFails", flag);

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

CAMLprim value sunml_arkode_spils_set_eps_lin(value varkode_mem, value eplifac)
{
    CAMLparam2(varkode_mem, eplifac);

    int flag = ARKSpilsSetEpsLin(ARKODE_MEM_FROM_ML(varkode_mem),
				 Double_val(eplifac));
    CHECK_FLAG("ARKSpilsSetEpsLin", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_set_mass_eps_lin(value varkode_mem, value eplifac)
{
    CAMLparam2(varkode_mem, eplifac);

    int flag = ARKSpilsSetMassEpsLin(ARKODE_MEM_FROM_ML(varkode_mem),
				     Double_val(eplifac));
    CHECK_FLAG("ARKSpilsSetMassEpsLin", flag);

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

/* statistic accessor functions */

CAMLprim value sunml_arkode_get_tol_scale_factor(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    realtype r;
    int flag = ARKodeGetTolScaleFactor(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetTolScaleFactor", flag);

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_arkode_get_num_nonlin_solv_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKodeGetNumNonlinSolvIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_get_num_nonlin_solv_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKodeGetNumNonlinSolvConvFails(ARKODE_MEM_FROM_ML(varkode_mem),
					       &r);
    CHECK_FLAG("ARKodeGetNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_get_nonlin_solv_stats(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    long int nniters, nncfails;
    int flag = ARKodeGetNonlinSolvStats(ARKODE_MEM_FROM_ML(varkode_mem),
				       &nniters, &nncfails);
    CHECK_FLAG("ARKodeGetNonlinSolvStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nniters));
    Store_field(r, 1, Val_long(nncfails));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_get_num_g_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKodeGetNumGEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKodeGetNumGEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_dls_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = ARKDlsGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				  &lenrwLS, &leniwLS);
    CHECK_FLAG("ARKDlsGetWorkSpace", flag);

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

    int flag = ARKDlsGetMassWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				      &lenrwLS, &leniwLS);
    CHECK_FLAG("ARKDlsGetMassWorkSpace", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_dls_get_num_jac_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKDlsGetNumJacEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_dls_get_num_mass_setups(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if SUNDIALS_LIB_VERSION >= 300
    int flag = ARKDlsGetNumMassSetups(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumMassSetups", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_dls_get_num_mass_solves(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
#if SUNDIALS_LIB_VERSION >= 300
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
#if SUNDIALS_LIB_VERSION >= 300
    int flag = ARKDlsGetNumMassMult(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumMassMult", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_dls_get_num_rhs_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKDlsGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKDlsGetNumRhsEvals", flag);

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

CAMLprim value sunml_arkode_spils_get_num_lin_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumLinIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION <= 301
int ARKSpilsGetNumMTSetups(void *arkode_mem, long int *nmtsetups);
#endif

CAMLprim value sunml_arkode_spils_get_num_mtsetup_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if SUNDIALS_LIB_VERSION >= 300
    int flag = ARKSpilsGetNumMTSetups(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMTSetups", flag);
#else
    r = 0;
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_mtimes_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;

#if SUNDIALS_LIB_VERSION >= 263
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

    flag = ARKSpilsGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				&lenrw, &leniw);
    CHECK_FLAG("ARKSpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_spils_get_num_prec_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumPrecEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_prec_solves(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumPrecSolves(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_jtsetup_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r;
#if SUNDIALS_LIB_VERSION >= 300
    int flag = ARKSpilsGetNumJTSetupEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumJTSetupEvals", flag);
#else
    r = 0;
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_jtimes_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumJtimesEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_rhs_evals (value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumRhsEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

/* spils mass functions */

CAMLprim value sunml_arkode_spils_mass_spgmr (value varkode_mem,
					  value vmaxl, value vtype)
{
    CAMLparam3 (varkode_mem, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSpgmr (arkode_mem, sunml_lsolver_precond_type (vtype),
				     Int_val (vmaxl),
				     masstimesfn,
				     arkode_mem->ark_user_data);
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
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSpbcg (arkode_mem, sunml_lsolver_precond_type (vtype),
				     Int_val (vmaxl),
				     masstimesfn,
				     arkode_mem->ark_user_data);
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
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSptfqmr (arkode_mem, sunml_lsolver_precond_type (vtype),
				       Int_val (vmaxl),
				       masstimesfn,
				       arkode_mem->ark_user_data);
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
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSpfgmr (arkode_mem, sunml_lsolver_precond_type (vtype),
				      Int_val (vmaxl),
				      masstimesfn,
				      arkode_mem->ark_user_data);
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
    ARKodeMem arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassPcg (arkode_mem, sunml_lsolver_precond_type (vtype),
				   Int_val (vmaxl),
				   masstimesfn,
				   arkode_mem->ark_user_data);
    CHECK_FLAG ("ARKMassPcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_spils_get_num_mass_iters(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumMassIters(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMassIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_mass_conv_fails(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumMassConvFails(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMassConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_mass_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = ARKSpilsGetMassWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				    &lenrw, &leniw);
    CHECK_FLAG("ARKSpilsGetMassWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_arkode_spils_get_num_mass_prec_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumMassPrecEvals(ARKODE_MEM_FROM_ML(varkode_mem), &r);
    CHECK_FLAG("ARKSpilsGetNumMassPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_arkode_spils_get_num_mass_prec_solves(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    long int r;
    int flag = ARKSpilsGetNumMassPrecSolves(ARKODE_MEM_FROM_ML(varkode_mem),
					    &r);
    CHECK_FLAG("ARKSpilsGetNumMassPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

