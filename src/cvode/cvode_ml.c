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
#include <cvodes/cvodes_direct.h>
#include <cvodes/cvodes_spils.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_bandpre.h>

#if SUNDIALS_LIB_VERSION < 300
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#endif

#include <cvodes/cvodes_impl.h>

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <cvodes/cvodes_lapack.h>
#endif

#else
/* CVODE (without sensitivity) */

#include <cvode/cvode.h>

/* linear solvers */
#include <cvode/cvode_direct.h>
#include <cvode/cvode_spils.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_bandpre.h>

#if SUNDIALS_LIB_VERSION < 300
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#endif
#include <cvode/cvode_impl.h>

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <cvode/cvode_lapack.h>
#endif

#endif

#include "../lsolvers/sundials_matrix_ml.h"
#include "../lsolvers/sundials_linearsolver_ml.h"
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
	sundials_ml_warn_discarded_exn (Extract_exception (r),
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

int cvode_translate_exception (value session, value exn,
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

static int rhsfn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal1(session);
    CAMLlocalN(args, 3);

    WEAK_DEREF (session, *(value*)user_data);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(Field(session, RECORD_CVODE_SESSION_RHSFN),
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
	    sundials_ml_warn_discarded_exn (r, "user-defined error weight fun");
	CAMLreturnT (int, -1);
    }

    CAMLreturnT (int, 0);
}

value cvode_make_jac_arg(realtype t, N_Vector y, N_Vector fy, value tmp)
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

value cvode_make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
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
static int jacfn(
	realtype t,
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

    args[0] = cvode_make_jac_arg (t, y, fy,
				  cvode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = MAT_BACKLINK(Jac);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#else // SUNDIALS_LIB_VERSION < 300

/* Dense and band Jacobians only work with serial NVectors.  */
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
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, c_matrix_dense_wrap(Jac));
	Store_field(cb, 1, dmat);
    }

    args[0] = cvode_make_jac_arg (t, y, fy,
				  cvode_make_triple_tmp (tmp1, tmp2, tmp3));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
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
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, c_matrix_band_wrap(Jac));
	Store_field(cb, 1, bmat);
    }

    args[0] = cvode_make_jac_arg(t, y, fy,
				 cvode_make_triple_tmp(tmp1, tmp2, tmp3));
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

    args[0] = cvode_make_jac_arg(t, y, fy, Val_unit);
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

    CAMLreturnT(int, cvode_translate_exception (session, r, RECOVERABLE));
}

static value make_spils_solve_arg(
	N_Vector r,
	realtype gamma,
	realtype delta,
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
	realtype t,
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

    args[0] = cvode_make_jac_arg(t, y, fy, Val_unit);
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
		      realtype t,
		      N_Vector y,
		      N_Vector fy,
		      void *user_data,
		      N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = cvode_make_jac_arg(t, y, fy, NVEC_BACKLINK(tmp));
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

#if SUNDIALS_LIB_VERSION >= 300
static int jacsetupfn(realtype t,
		      N_Vector y,
		      N_Vector fy,
		      void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = cvode_make_jac_arg(t, y, fy, Val_unit);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, arg);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int linit(CVodeMem cv_mem)
{
    CAMLparam0();
    CAMLlocal2 (session, cb);

    WEAK_DEREF (session, *(value*)cv_mem->cv_user_data);

    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODE_ALTERNATE_CALLBACKS_LINIT);
    cb = Some_val(cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(cb, session);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

static int lsetup(CVodeMem cv_mem, int convfail,
		  N_Vector ypred, N_Vector fpred, booleantype *jcurPtr,
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    WEAK_DEREF (session, *(value*)cv_mem->cv_user_data);

    switch (convfail) {
    case CV_NO_FAILURES:
	convfail = VARIANT_CVODE_ALTERNATE_NO_FAILURES;
	break;

    case CV_FAIL_BAD_J:
	convfail = VARIANT_CVODE_ALTERNATE_FAIL_BAD_J;
	break;

    case CV_FAIL_OTHER:
	convfail = VARIANT_CVODE_ALTERNATE_FAIL_OTHER;
	break;
    }
    args = caml_alloc_tuple (RECORD_CVODE_ALTERNATE_LSETUP_ARGS_SIZE);
    Store_field (args, RECORD_CVODE_ALTERNATE_LSETUP_ARGS_CONV_FAIL,
		 Val_int (convfail));
    Store_field (args, RECORD_CVODE_ALTERNATE_LSETUP_ARGS_Y,
		 NVEC_BACKLINK(ypred));
    Store_field (args, RECORD_CVODE_ALTERNATE_LSETUP_ARGS_RHS,
		 NVEC_BACKLINK(fpred));
    Store_field (args, RECORD_CVODE_ALTERNATE_LSETUP_ARGS_TMP,
		 cvode_make_triple_tmp(tmp1, tmp2, tmp3));

    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODE_ALTERNATE_CALLBACKS_LSETUP);
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

static int lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector ycur,
		  N_Vector fcur)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_CVODE_ALTERNATE_LSOLVE_ARGS_SIZE);
    Store_field (args, RECORD_CVODE_ALTERNATE_LSOLVE_ARGS_Y,
		 NVEC_BACKLINK (ycur));
    Store_field (args, RECORD_CVODE_ALTERNATE_LSOLVE_ARGS_RHS,
		 NVEC_BACKLINK (fcur));
    Store_field (args, RECORD_CVODE_ALTERNATE_LSOLVE_ARGS_EWT,
		 NVEC_BACKLINK (weight));

    WEAK_DEREF (session, *(value*)cv_mem->cv_user_data);

    cb = Field (session, RECORD_CVODE_SESSION_LS_CALLBACKS);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODE_ALTERNATE_CALLBACKS_LSOLVE);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (cb, session, args, NVEC_BACKLINK (b));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

CAMLprim value sunml_cvode_get_gamma(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);
    CVodeMem cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    r = caml_alloc_small (2 * Double_wosize, Double_array_tag);
    Store_double_field (r, 0, cvode_mem->cv_gamma);
    Store_double_field (r, 1, cvode_mem->cv_gammap);
    CAMLreturn(r);
}

CAMLprim value sunml_cvode_set_alternate (value vcvode_mem, value vhas_init,
				      value vhas_setup)
{
    CAMLparam3(vcvode_mem, vhas_init, vhas_setup);
    CVodeMem cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    cvode_mem->cv_linit  = Bool_val(vhas_init)  ? linit : NULL;
    cvode_mem->cv_lsetup  = Bool_val(vhas_setup) ? lsetup : NULL;
#if SUNDIALS_LIB_VERSION < 300
    cvode_mem->cv_setupNonNull = Bool_val(vhas_setup);
#endif
    cvode_mem->cv_lsolve = lsolve;
    cvode_mem->cv_lmem   = NULL;

    CAMLreturn (Val_unit);
}

#ifdef SUNDIALSML_WITHSENS
CAMLprim value sunml_cvode_adj_set_alternate (value vparent, value vwhich,
					  value vhas_init, value vhas_setup)
{
    CAMLparam4(vparent, vwhich, vhas_init, vhas_setup);
    CVodeMem cvode_mem = CVodeGetAdjCVodeBmem(CVODE_MEM_FROM_ML (vparent),
					      Int_val(vwhich));

    cvode_mem->cv_linit  = Bool_val(vhas_init)  ? linit : NULL;
    cvode_mem->cv_lsetup  = Bool_val(vhas_setup) ? lsetup : NULL;
#if SUNDIALS_LIB_VERSION < 300
    cvode_mem->cv_setupNonNull = Bool_val(vhas_setup);
#endif
    cvode_mem->cv_lsolve = lsolve;
    cvode_mem->cv_lmem   = NULL;

    CAMLreturn (Val_unit);
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

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
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

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
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

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
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

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
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

CAMLprim value sunml_cvode_dls_set_linear_solver (value vcvode_mem, value vlsolv,
					      value vjmat, value vhasjac)
{
    CAMLparam4(vcvode_mem, vlsolv, vjmat, vhasjac);
#if SUNDIALS_LIB_VERSION >= 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVDlsSetLinearSolver(cvode_mem, lsolv, jmat);
    CHECK_DLS_FLAG ("CVDlsSetLinearSolver", flag);
    if (Bool_val (vhasjac)) {
	flag = CVDlsSetJacFn(cvode_mem, jacfn);
	CHECK_DLS_FLAG("CVDlsSetJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_linear_solver (value vcvode_mem, value vlsolv)
{
    CAMLparam2(vcvode_mem, vlsolv);
#if SUNDIALS_LIB_VERSION >= 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpilsSetLinearSolver(cvode_mem, lsolv);
    CHECK_SPILS_FLAG ("CVSpilsSetLinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_preconditioner (value vsession,
						 value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = CVODE_MEM_FROM_ML (vsession);
    CVSpilsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = CVSpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_SPILS_FLAG ("CVSpilsSetPreconditioner", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_banded_preconditioner (value vsession,
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

CAMLprim value sunml_cvode_spils_set_jac_times(value vdata, value vhas_setup,
							value vhas_times)
{
    CAMLparam3(vdata, vhas_setup, vhas_times);
#if SUNDIALS_LIB_VERSION >= 300
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

    int iter_c;
    if (Is_block(iter)) {
	iter_c = CV_NEWTON;
    } else {
	iter_c = CV_FUNCTIONAL;
    }

    void *cvode_mem = CVodeCreate(lmm_c, iter_c);
    if (cvode_mem == NULL)
	caml_failwith("CVodeCreate returned NULL");

    vcvode_mem = caml_alloc_final(1, NULL, 1, 5);
    CVODE_MEM(vcvode_mem) = cvode_mem;

    N_Vector initial_nv = NVEC_VAL(initial);
    flag = CVodeInit(cvode_mem, rhsfn, Double_val(t0), initial_nv);
    if (flag != CV_SUCCESS) {
	CVodeFree (&cvode_mem);
	CHECK_FLAG("CVodeInit", flag);
    }

    value *backref = c_sundials_malloc_value(weakref);
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
    realtype tret;
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
	CHECK_FLAG ("CVode", flag);
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


void sunml_cvode_check_flag(const char *call, int flag)
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
	    caml_raise_constant(CVODE_EXN(LinearSetupFailure));

	case CV_LSOLVE_FAIL:
	    caml_raise_constant(CVODE_EXN(LinearSolveFailure));

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

	case CV_BAD_K:
	    caml_raise_constant(CVODE_EXN(BadK));

	case CV_BAD_T:
	    caml_raise_constant(CVODE_EXN(BadT));

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

void sunml_cvode_check_dls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CVDLS_SUCCESS) return;

    switch (flag) {
	case CVDLS_ILL_INPUT:
	    caml_raise_constant(CVODE_EXN(IllInput));

	case CVDLS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if SUNDIALS_LIB_VERSION >= 300
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

#if SUNDIALS_LIB_VERSION >+ 300
	case CVSPILS_SUNLS_FAIL:
#endif
	default:
	    /* e.g. CVSPILS_MEM_NULL, CVSPILS_PMEM_NULL, CVSPILS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     CVSpilsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

/* basic interface */

CAMLprim value sunml_cvode_session_finalize(value vdata)
{
    if (CVODE_MEM_FROM_ML(vdata) != NULL) {
	void *cvode_mem = CVODE_MEM_FROM_ML(vdata);
	value *backref = CVODE_BACKREF_FROM_ML(vdata);
	CVodeFree(&cvode_mem);
	c_sundials_free_value(backref);
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

    realtype hinused;
    realtype hlast;
    realtype hcur;
    realtype tcur;

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
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vdata), CV_FUNCTIONAL);
    CHECK_FLAG ("CVodeSetIterType", flag);
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
				  lsolver_precond_type(vptype));
    CHECK_SPILS_FLAG("CVSpilsSetPrecType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_diag (value vcvode_mem)
{
    CAMLparam1 (vcvode_mem);
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vcvode_mem), CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
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

CAMLprim value sunml_cvode_get_actual_init_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetActualInitStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_cvode_get_last_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetLastStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetLastStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_cvode_get_current_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetCurrentStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_cvode_get_current_time(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

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
				lsolver_gs_type(vgstype));
    CHECK_SPILS_FLAG("CVSpilsSetGSType", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_set_eps_lin(value vcvode_mem, value eplifac)
{
    CAMLparam2(vcvode_mem, eplifac);

    int flag = CVSpilsSetEpsLin(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(eplifac));
    CHECK_SPILS_FLAG("CVSpilsSetEpsLin", flag);

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

    realtype r;
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

CAMLprim value sunml_cvode_dls_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = CVDlsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrwLS, &leniwLS);
    CHECK_DLS_FLAG("CVDlsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_dls_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDlsGetNumJacEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_DLS_FLAG("CVDlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_dls_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDlsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_DLS_FLAG("CVDlsGetNumRhsEvals", flag);

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

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpgmr (cvode_mem, lsolver_precond_type (vtype), Int_val (vmaxl));
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

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpbcg (cvode_mem, lsolver_precond_type (vtype), Int_val (vmaxl));
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

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSptfqmr (cvode_mem, lsolver_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSptfqmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_spils_get_num_lin_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumLinIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_spils_get_num_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_spils_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = CVSpilsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_SPILS_FLAG("CVSpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_cvode_spils_get_num_prec_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_spils_get_num_prec_solves(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecSolves(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_spils_get_num_jtsetup_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    long int r;
#if SUNDIALS_LIB_VERSION >= 300
    int flag = CVSpilsGetNumJTSetupEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumJTSetupEvals", flag);
#else
    r = 0;
#endif
    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_spils_get_num_jtimes_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumJtimesEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_cvode_spils_get_num_rhs_evals (value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_SPILS_FLAG("CVSpilsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

