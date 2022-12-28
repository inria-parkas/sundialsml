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
#include <cvodes/cvodes.h>
#include <sundials/sundials_band.h>

#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

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

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <cvodes/cvodes_lapack.h>
#endif

#include "../lsolvers/sundials_matrix_ml.h"
#include "../lsolvers/sundials_linearsolver_ml.h"
#include "../lsolvers/sundials_nonlinearsolver_ml.h"
#include "../sundials/sundials_ml.h"
#include "../cvode/cvode_ml.h"
#include "cvodes_ml.h"
#include "../nvectors/nvector_ml.h"

CAMLprim value sunml_cvodes_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (CVODES, exns);
    CAMLreturn (Val_unit);
}

/* Interface with nvectors */

CAMLprim value sunml_cvodes_alloc_nvector_array(value vn)
{
    CAMLparam1(vn);
    value r;
    int n = Int_val(vn);

    if (n > 0) {
	r = caml_alloc_tuple(n);
    } else {
	r = Val_int(0);
    }

    CAMLreturn(r);
}

static value make_double_tmp(N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, NVEC_BACKLINK(tmp1));
    Store_field(r, 1, NVEC_BACKLINK(tmp2));
    CAMLreturn(r);
}


/* Callbacks */

static int quadrhsfn(sunrealtype t, N_Vector y, N_Vector yQdot, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal2(session, cb);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(yQdot);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_SENSEXT_FROM_ML (session);
    cb = CVODES_QUADRHSFN_FROM_EXT (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int sensrhsfn(int ns, sunrealtype t, N_Vector y, N_Vector ydot,
		     N_Vector *ys, N_Vector *ysdot, void *user_data,
		     N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal3(args, session, sensext);

    WEAK_DEREF (session, *(value*)user_data);
    sensext = CVODE_SENSEXT_FROM_ML(session);

    args = caml_alloc_tuple (RECORD_CVODES_SENSRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_YP, NVEC_BACKLINK (ydot));
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_TMP,
		 make_double_tmp (tmp1, tmp2));

    sunml_nvectors_into_array(ns, CVODES_SENSARRAY1_FROM_EXT(sensext), ys);
    sunml_nvectors_into_array(ns, CVODES_SENSARRAY2_FROM_EXT(sensext), ysdot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(CVODES_SENSRHSFN_FROM_EXT(sensext),
				 args,
				 CVODES_SENSARRAY1_FROM_EXT(sensext),
				 CVODES_SENSARRAY2_FROM_EXT(sensext));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int sensrhsfn1(int ns, sunrealtype t, N_Vector y, N_Vector ydot,
		      int is, N_Vector ys, N_Vector ysdot, void *user_data,
		      N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal4(args, sens, session, cb);

    WEAK_DEREF (session, *(value*)user_data);

    args = caml_alloc_tuple (RECORD_CVODES_SENSRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_YP, NVEC_BACKLINK (ydot));
    Store_field (args, RECORD_CVODES_SENSRHSFN_ARGS_TMP,
		 make_double_tmp (tmp1, tmp2));

    cb = CVODE_SENSEXT_FROM_ML (session);
    cb = CVODES_SENSRHSFN1_FROM_EXT (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, Val_int (is));
    if (! Is_exception_result (r))
	r = caml_callback3_exn (r, args, NVEC_BACKLINK (ys),
				NVEC_BACKLINK (ysdot));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int quadsensrhsfn(int ns, sunrealtype t, N_Vector y, N_Vector *ys,
		         N_Vector yqdot, N_Vector *yqsdot, void *user_data,
		         N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal3(args, session, sensext);

    WEAK_DEREF (session, *(value*)user_data);
    sensext = CVODE_SENSEXT_FROM_ML(session);

    args = caml_alloc_tuple (RECORD_CVODES_QUADSENSRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_QUADSENSRHSFN_ARGS_T,
		 caml_copy_double (t));
    Store_field (args, RECORD_CVODES_QUADSENSRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_QUADSENSRHSFN_ARGS_SENS,
		 CVODES_SENSARRAY1_FROM_EXT(sensext));
    Store_field (args, RECORD_CVODES_QUADSENSRHSFN_ARGS_YQP,
		 NVEC_BACKLINK(yqdot));
    Store_field (args, RECORD_CVODES_QUADSENSRHSFN_ARGS_TMP,
		 make_double_tmp (tmp1, tmp2));

    sunml_nvectors_into_array(ns, CVODES_SENSARRAY1_FROM_EXT(sensext), ys);
    sunml_nvectors_into_array(ns, CVODES_SENSARRAY2_FROM_EXT(sensext), yqsdot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(CVODES_QUADSENSRHSFN_FROM_EXT(sensext), args,
				 CVODES_SENSARRAY2_FROM_EXT(sensext));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int brhsfn(sunrealtype t, N_Vector y, N_Vector yb, N_Vector ybdot,
		  void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_CVODES_ADJ_BRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_YB, NVEC_BACKLINK (yb));

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_SENSEXT_FROM_ML (session);
    cb = CVODES_BRHSFN_FROM_EXT (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (cb, args, NVEC_BACKLINK (ybdot));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int brhsfn_sens(sunrealtype t, N_Vector y, N_Vector *ys, N_Vector yb,
		       N_Vector ybdot, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, sensext);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    sensext = CVODE_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(sensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));

    args = caml_alloc_tuple (RECORD_CVODES_ADJ_BRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_YB, NVEC_BACKLINK (yb));

    sunml_nvectors_into_array(ns, CVODES_BSENSARRAY_FROM_EXT(sensext), ys);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(CVODES_BRHSFN_SENS_FROM_EXT(sensext),
				 args, CVODES_BSENSARRAY_FROM_EXT(sensext),
				 NVEC_BACKLINK (ybdot));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bquadrhsfn(sunrealtype t, N_Vector y, N_Vector yb, N_Vector qbdot,
		      void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_T,
		 caml_copy_double(t));
    Store_field (args, RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_Y, NVEC_BACKLINK(y));
    Store_field (args, RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_YB,
		 NVEC_BACKLINK (yb));

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_SENSEXT_FROM_ML (session);
    cb = CVODES_BQUADRHSFN_FROM_EXT(cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (cb, args, NVEC_BACKLINK(qbdot));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bquadrhsfn_sens(sunrealtype t, N_Vector y, N_Vector *ys, N_Vector yb,
			   N_Vector qbdot, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, sensext);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    sensext = CVODE_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(sensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));

    args = caml_alloc_tuple (RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_T,
		 caml_copy_double (t));
    Store_field (args, RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_Y,
		 NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_ADJ_BQUADRHSFN_ARGS_YB,
		 NVEC_BACKLINK (yb));

    sunml_nvectors_into_array(ns, CVODES_BSENSARRAY_FROM_EXT(sensext), ys);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (CVODES_BQUADRHSFN_SENS_FROM_EXT(sensext),
				  args, CVODES_BSENSARRAY_FROM_EXT(sensext),
				  NVEC_BACKLINK (qbdot));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

value sunml_cvodes_make_jac_arg(sunrealtype t, N_Vector y, N_Vector yb,
			  N_Vector fyb, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_CVODES_ADJ_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_YB, NVEC_BACKLINK(yb));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_FYB, NVEC_BACKLINK(fyb));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

static value make_spils_solve_arg(
	N_Vector rvecb,
	sunrealtype gammab,
	sunrealtype deltab,
	int lrb)

{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_SIZE);
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_RVEC, NVEC_BACKLINK(rvecb));
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_GAMMA,
                caml_copy_double(gammab));
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_DELTA,
                caml_copy_double(deltab));
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_LR, Val_bool (lrb == 1));

    CAMLreturn(v);
}

static int bprecsolvefn(
	sunrealtype t,
	N_Vector y,
	N_Vector yb,
	N_Vector fyb,
	N_Vector rvecb,
	N_Vector zvecb,
	sunrealtype gammab,
	sunrealtype deltab,
	int lrb,
	void *user_data
#if SUNDIALS_LIB_VERSION < 300
	,
	N_Vector tmpb
#endif
	)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal2(session, cb);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb, Val_unit);
    args[1] = make_spils_solve_arg(rvecb, gammab, deltab, lrb);
    args[2] = NVEC_BACKLINK(zvecb);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODES_BSPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bprecsolvefn_withsens(
	sunrealtype t,
	N_Vector y,
	N_Vector *yS,
	N_Vector yb,
	N_Vector fyb,
	N_Vector rvecb,
	N_Vector zvecb,
	sunrealtype gammab,
	sunrealtype deltab,
	int lrb,
	void *user_data
#if SUNDIALS_LIB_VERSION < 300
	,
	N_Vector tmpb
#endif
	)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb, Val_unit);
    args[1] = make_spils_solve_arg(rvecb, gammab, deltab, lrb);

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[2] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[2], yS);

    args[3] = NVEC_BACKLINK(zvecb);

    cb = CVODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODES_BSPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int bprecsetupfn(
    sunrealtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    sunbooleantype jokb,
    sunbooleantype *jcurPtrB,
    sunrealtype gammab,
    void *user_data
#if SUNDIALS_LIB_VERSION < 300
    ,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b
#endif
    )
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb, Val_unit);
    args[1] = Val_bool(jokb);
    args[2] = caml_copy_double(gammab);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODES_BSPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    /* Update jcurPtr; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcurPtrB = Bool_val (r);
	CAMLreturnT (int, 0);
    }

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bprecsetupfn_withsens(
    sunrealtype t,
    N_Vector y,
    N_Vector *yS,
    N_Vector yb,
    N_Vector fyb,
    sunbooleantype jokb,
    sunbooleantype *jcurPtrB,
    sunrealtype gammab,
    void *user_data
#if SUNDIALS_LIB_VERSION < 300
    ,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b
#endif
    )
{
    CAMLparam0();
    CAMLlocal3(session, bsensext, cb);
    CAMLlocalN(args, 4);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb, Val_unit);

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[1], yS);

    args[2] = Val_bool(jokb);
    args[3] = caml_copy_double(gammab);

    cb = CVODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODES_BSPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    /* Update jcurPtr; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcurPtrB = Bool_val (r);
	CAMLreturnT (int, 0);
    }

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

#if SUNDIALS_LIB_VERSION >= 300
static int bjacsetupfn(
    sunrealtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = sunml_cvodes_make_jac_arg(t, y, yb, fyb, Val_unit);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, arg);

    /* NB: jac_times_vec doesn't accept RecoverableFailure. */
    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

static int bjacsetupfn_withsens(
    sunrealtype t,
    N_Vector y,
    N_Vector *yS,
    N_Vector yb,
    N_Vector fyb,
    void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, bsensext, cb);
    CAMLlocalN(args, 2);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb, Val_unit);

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[1], yS);

    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 2, args);

    /* NB: jac_times_vec doesn't accept RecoverableFailure. */
    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}
#endif

static int bjactimesfn(
    N_Vector vb,
    N_Vector Jvb,
    sunrealtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    void *user_data,
    N_Vector tmpb)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb, NVEC_BACKLINK(tmpb));
    args[1] = NVEC_BACKLINK(vb);
    args[2] = NVEC_BACKLINK(Jvb);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    /* NB: jac_times_vec doesn't accept RecoverableFailure. */
    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bjactimesfn_withsens(
    N_Vector vb,
    N_Vector Jvb,
    sunrealtype t,
    N_Vector y,
    N_Vector *yS,
    N_Vector yb,
    N_Vector fyb,
    void *user_data,
    N_Vector tmpb)
{
    CAMLparam0();
    CAMLlocal3(session, bsensext, cb);
    CAMLlocalN(args, 4);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb, NVEC_BACKLINK(tmpb));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[1], yS);

    args[2] = NVEC_BACKLINK(vb);
    args[3] = NVEC_BACKLINK(Jvb);

    cb = CVODE_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    /* NB: jac_times_vec doesn't accept RecoverableFailure. */
    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}
#endif

#if 300 <= SUNDIALS_LIB_VERSION

static int bjacfn_nosens(
    sunrealtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    SUNMatrix jacb,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));
    args[1] = MAT_BACKLINK(jacb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bjacfn_withsens(
    sunrealtype t,
    N_Vector y,
    N_Vector *ys,
    N_Vector yb,
    N_Vector fyb,
    SUNMatrix jacb,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[1], ys);

    args[2] = MAT_BACKLINK(jacb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#else

static int bjacfn_nosens(
    long int neqb,
    sunrealtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    DlsMat jacb,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(jacb));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bjacfn_withsens(
    long int neqb,
    sunrealtype t,
    N_Vector y,
    N_Vector *ys,
    N_Vector yb,
    N_Vector fyb,
    DlsMat jacb,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal4(session, bsensext, cb, dmat);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(jacb));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[1], ys);

    args[2] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

static int bbandjacfn_nosens(
	long int nb,
	long int mupperb,
	long int mlowerb,
	sunrealtype t,
	N_Vector y,
	N_Vector yb,
	N_Vector fyb,
	DlsMat jacb,
	void *user_data, 	 
	N_Vector tmp1b,
	N_Vector tmp2b,
	N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, bmat);
    value *backref = user_data;
    WEAK_DEREF (session, *backref);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(jacb));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));
    args[1] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bbandjacfn_withsens(
	long int nb,
	long int mupperb,
	long int mlowerb,
	sunrealtype t,
	N_Vector y,
	N_Vector *ys,
	N_Vector yb,
	N_Vector fyb,
	DlsMat jacb,
	void *user_data, 	 
	N_Vector tmp1b,
	N_Vector tmp2b,
	N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal4(session, bsensext, cb, bmat);
    value *backref = user_data;
    int ns;

    WEAK_DEREF (session, *backref);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(jacb));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[2], ys);

    args[2] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

#endif

#if 500 <= SUNDIALS_LIB_VERSION
static int blinsysfn_nosens(
    sunrealtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    SUNMatrix jacb,
    sunbooleantype jokb,
    sunbooleantype *jcurb,
    sunrealtype gammab,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1); // second in pair

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));
    args[1] = MAT_BACKLINK(jacb);
    args[2] = Val_bool(jokb);
    args[3] = caml_copy_double(gammab);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), // linsys_fn -> no/with_sens
				  4, args);

    /* Update jcur; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcurb = Bool_val (r);
	CAMLreturnT(int, 0);
    }
    r = Extract_exception (r);

    CAMLreturnT(int, sunml_cvode_translate_exception (session, r, RECOVERABLE));
}

static int blinsysfn_withsens(
    sunrealtype t,
    N_Vector y,
    N_Vector *ys,
    N_Vector yb,
    N_Vector fyb,
    SUNMatrix jacb,
    sunbooleantype jokb,
    sunbooleantype *jcurb,
    sunrealtype gammab,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 1); // second in pair

    args[0] = sunml_cvodes_make_jac_arg(t, y, yb, fyb,
		sunml_cvode_make_triple_tmp(tmp1b, tmp2b, tmp3b));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[1], ys);

    args[2] = MAT_BACKLINK(jacb);
    args[3] = Val_bool(jokb);
    args[4] = caml_copy_double(gammab);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), // linsys_fn -> no/with_sens
				  5, args);

    /* Update jcur; leave it unchanged if an error occurred.  */
    if (!Is_exception_result (r)) {
	*jcurb = Bool_val (r);
	CAMLreturnT(int, 0);
    }
    r = Extract_exception (r);

    CAMLreturnT(int, sunml_cvode_translate_exception (session, r, RECOVERABLE));
}
#endif

#if 530 <= SUNDIALS_LIB_VERSION
static int bjactimesrhsfn(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
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

/* quadrature interface */

CAMLprim value sunml_cvodes_quad_init(value vdata, value vq0)
{
    CAMLparam2(vdata, vq0);
    CAMLlocal1(r);
    int flag;
    N_Vector q0 = NVEC_VAL(vq0);
    
    flag = CVodeQuadInit(CVODE_MEM_FROM_ML(vdata), quadrhsfn, q0);
    SCHECK_FLAG("CVodeQuadInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quad_reinit(value vdata, value vq0)
{
    CAMLparam2(vdata, vq0);
    CAMLlocal1(r);
    int flag;
    N_Vector q0 = NVEC_VAL(vq0);
    
    flag = CVodeQuadReInit(CVODE_MEM_FROM_ML(vdata), q0);
    SCHECK_FLAG("CVodeQuadReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quad_sv_tolerances(value vdata, value reltol,
					   value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = CVodeQuadSVtolerances(CVODE_MEM_FROM_ML(vdata),
	    Double_val(reltol), atol_nv);
    SCHECK_FLAG("CVodeQuadSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quad_get(value vdata, value vyq)
{
    CAMLparam2(vdata, vyq);
    N_Vector yq = NVEC_VAL(vyq);
    sunrealtype tret;

    int flag = CVodeGetQuad(CVODE_MEM_FROM_ML(vdata), &tret, yq);
    SCHECK_FLAG("CVodeGetQuad", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_cvodes_quad_get_dky(value vdata, value vt, value vk,
				     value vdkyq)
{
    CAMLparam4(vdata, vt, vk, vdkyq);
    N_Vector dkyq = NVEC_VAL(vdkyq);

    int flag = CVodeGetQuadDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkyq);
	    
    SCHECK_FLAG("CVodeGetQuadDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quad_get_err_weights(value vdata, value veqweight)
{
    CAMLparam2(vdata, veqweight);
    N_Vector eqweight = NVEC_VAL(veqweight);

    int flag = CVodeGetQuadErrWeights(CVODE_MEM_FROM_ML(vdata), eqweight);
    SCHECK_FLAG("CVodeGetQuadErrWeights", flag);

    CAMLreturn (Val_unit);
}

/* sensitivity interface */

#if 400 <= SUNDIALS_LIB_VERSION
// hack to work around lack of CVodeGetUserData
typedef struct {
#if 600 <= SUNDIALS_LIB_VERSION
  SUNContext cv_sunctx;
#endif
  sunrealtype cv_uround;
  CVRhsFn cv_f;
  void *cv_user_data;
  //...
} *StartOf_CVodeMem;

static value sunml_cvodes_session_to_value(void *cvode_mem)
{
    value session;
    // void *user_data = CVodeGetUserData(cvode_mem);
    void *user_data = ((StartOf_CVodeMem)cvode_mem)->cv_user_data;

    WEAK_DEREF (session, *(value*)user_data);
    return session;
}

static value sunml_cvodes_bsession_to_value(void *cvode_mem)
{
    CAMLparam0();
    CAMLlocal3(vparent, vsession, vchildptr);

    // for the "B" case, the callback functions are passed the
    // session for the backward integrator directly.
    //
    // The user_data field of Adjoint sessions contains a pointer to
    // the parent session structure (see CVodeCreateB). The actual
    // user data is stored in the associated CVodeBMemRec structure.
    //
    // A CVodeGetUserDataB function would not help, since we do not know
    // what the `which` value is.

    // void *user_data = CVodeGetUserData(cvode_mem);
    void *parent = ((StartOf_CVodeMem)cvode_mem)->cv_user_data;
    void *parent_mem = ((StartOf_CVodeMem)parent)->cv_user_data;
    WEAK_DEREF (vparent, *(value*)parent_mem);
    vchildptr = sunml_wrap_session_pointer(cvode_mem);

    vsession = caml_callback2(
		*caml_named_value("Cvodes.revlookup_bsession"),
		vparent, vchildptr);

    assert (vsession != Val_none);
    vsession = Some_val(vsession);

    return vsession;
}

static void* sunml_cvodes_session_from_value(value vcvode_mem)
{
    return (CVODE_MEM_FROM_ML(vcvode_mem));
}
#endif

CAMLprim value sunml_cvodes_set_nonlinear_solver_sim(value vcvode_mem,
						     value vnlsolv)
{
    CAMLparam2(vcvode_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_cvodes_session_to_value,
				   sunml_cvodes_session_from_value);

    int flag = CVodeSetNonlinearSolverSensSim(cvode_mem, nlsolv);
    CHECK_FLAG ("CVodeSetNonlinearSolverSensSim", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_set_nonlinear_solver_stg(value vcvode_mem,
						     value vnlsolv)
{
    CAMLparam2(vcvode_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_cvodes_session_to_value,
				   sunml_cvodes_session_from_value);

    int flag = CVodeSetNonlinearSolverSensStg(cvode_mem, nlsolv);
    CHECK_FLAG ("CVodeSetNonlinearSolverSensStg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_set_nonlinear_solver_stg1(value vcvode_mem,
						      value vnlsolv)
{
    CAMLparam2(vcvode_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_cvodes_session_to_value,
				   sunml_cvodes_session_from_value);

    int flag = CVodeSetNonlinearSolverSensStg1(cvode_mem, nlsolv);
    CHECK_FLAG ("CVodeSetNonlinearSolverSensStg1", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_sv_tolerances(value vdata, value reltol,
					   value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = sunml_nvector_array_alloc(abstol);

    int flag = CVodeSensSVtolerances(CVODE_MEM_FROM_ML(vdata),
	    Double_val(reltol), atol_nv);
    sunml_nvector_array_free(atol_nv); 
    SCHECK_FLAG("CVodeSensSVtolerances", flag);

    CAMLreturn (Val_unit);
}

static int decode_sens_method(value vmethod)
{
    switch (Tag_val(vmethod)) {
    case VARIANT_CVODES_SENS_METHOD_SIMULTANEOUS:
	return CV_SIMULTANEOUS;

    case VARIANT_CVODES_SENS_METHOD_STAGGERED:
	return CV_STAGGERED;

    case VARIANT_CVODES_SENS_METHOD_STAGGERED1:
	return CV_STAGGERED1;

    default:
	caml_failwith("Illegal sens method.");
    }
}

CAMLprim value sunml_cvodes_sens_free(value vdata)
{
    CAMLparam1(vdata);
    CVodeSensFree(CVODE_MEM_FROM_ML(vdata));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_init(value vdata, value vmethod, value vrhsfn,
				  value vys0)
{
    CAMLparam4(vdata, vmethod, vrhsfn, vys0);
    int ns = Wosize_val (vys0); /* vys0 : nvector array */
    N_Vector *ys0 = sunml_nvector_array_alloc(vys0);
    void *cvode_mem = CVODE_MEM_FROM_ML(vdata);

    CVodeSensFree(cvode_mem);
    int flag = CVodeSensInit(cvode_mem, ns,
			     decode_sens_method(vmethod),
			     ((Bool_val(vrhsfn)) ? sensrhsfn : NULL),
			     ys0);
    sunml_nvector_array_free(ys0); 
    SCHECK_FLAG("CVodeSensInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_init_1(value vdata, value vmethod, value vrhsfn,
				    value vys0)
{
    CAMLparam4(vdata, vmethod, vrhsfn, vys0);
    int ns = Wosize_val (vys0); /* vys0 : nvector array */
    N_Vector *ys0 = sunml_nvector_array_alloc(vys0);
    void *cvode_mem = CVODE_MEM_FROM_ML(vdata);

    CVodeSensFree(cvode_mem);
    int flag = CVodeSensInit1(cvode_mem, ns,
			      decode_sens_method(vmethod),
			      ((Bool_val(vrhsfn)) ? sensrhsfn1 : NULL),
			      ys0);
    sunml_nvector_array_free(ys0); 
    SCHECK_FLAG("CVodeSensInit1", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_reinit(value vdata, value vmethod, value vs0)
{
    CAMLparam3(vdata, vmethod, vs0);
    CAMLlocal1(r);
    int flag;
    N_Vector *s0 = sunml_nvector_array_alloc(vs0);

    flag = CVodeSensReInit(CVODE_MEM_FROM_ML(vdata),
			   decode_sens_method(vmethod),
			   s0);
    sunml_nvector_array_free(s0);
    SCHECK_FLAG("CVodeSensReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_get(value vdata, value vys)
{
    CAMLparam2(vdata, vys);
    N_Vector *ys = sunml_nvector_array_alloc(vys);
    sunrealtype tret;

    int flag = CVodeGetSens(CVODE_MEM_FROM_ML(vdata), &tret, ys);
    sunml_nvector_array_free(ys);
    SCHECK_FLAG("CVodeGetSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_cvodes_sens_get_dky(value vdata, value vt, value vk, value vdkys)
{
    CAMLparam4(vdata, vt, vk, vdkys);
    N_Vector *dkys = sunml_nvector_array_alloc(vdkys);

    int flag = CVodeGetSensDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkys);
    sunml_nvector_array_free(dkys);
    SCHECK_FLAG("CVodeGetSensDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_get1(value vdata, value vis, value vys)
{
    CAMLparam3(vdata, vis, vys);
    N_Vector ys = NVEC_VAL(vys);
    sunrealtype tret;

    int flag = CVodeGetSens1(CVODE_MEM_FROM_ML(vdata), &tret, Int_val(vis), ys);
    SCHECK_FLAG("CVodeGetSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_cvodes_sens_get_dky1(value vdata, value vt, value vk,
				      value vis, value vdkys)
{
    CAMLparam5(vdata, vt, vk, vis, vdkys);
    N_Vector dkys = NVEC_VAL(vdkys);

    int flag = CVodeGetSensDky1(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), Int_val(vis), dkys);

    SCHECK_FLAG("CVodeGetSensDky1", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_get_err_weights(value vdata, value vesweight)
{
    CAMLparam2(vdata, vesweight);
    N_Vector *esweight = sunml_nvector_array_alloc(vesweight);

    int flag = CVodeGetSensErrWeights(CVODE_MEM_FROM_ML(vdata), esweight);
    sunml_nvector_array_free(esweight);
    SCHECK_FLAG("CVodeGetSensErrWeights", flag);

    CAMLreturn (Val_unit);
}

/* sensitivity/quadrature interface */

CAMLprim value sunml_cvodes_quadsens_init(value vdata, value vrhsfn, value vyqs0)
{
    CAMLparam3(vdata, vrhsfn, vyqs0);
    N_Vector *yqs0 = sunml_nvector_array_alloc(vyqs0);

    int flag = CVodeQuadSensInit(CVODE_MEM_FROM_ML(vdata),
				 Bool_val (vrhsfn) ? quadsensrhsfn : NULL,
				 yqs0);
    sunml_nvector_array_free(yqs0); 
    SCHECK_FLAG("CVodeQuadSensInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_reinit(value vdata, value vyqs0)
{
    CAMLparam2(vdata, vyqs0);
    N_Vector *yqs0 = sunml_nvector_array_alloc(vyqs0);

    int flag = CVodeQuadSensReInit(CVODE_MEM_FROM_ML(vdata), yqs0);
    sunml_nvector_array_free(yqs0); 
    SCHECK_FLAG("CVodeQuadSensReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_sv_tolerances(value vdata, value reltol,
					       value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = sunml_nvector_array_alloc(abstol);

    int flag = CVodeQuadSensSVtolerances(CVODE_MEM_FROM_ML(vdata),
	    Double_val(reltol), atol_nv);
    sunml_nvector_array_free(atol_nv); 
    SCHECK_FLAG("CVodeQuadSensSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_get(value vdata, value vyqs)
{
    CAMLparam2(vdata, vyqs);
    N_Vector *yqs = sunml_nvector_array_alloc(vyqs);
    sunrealtype tret;

    int flag = CVodeGetQuadSens(CVODE_MEM_FROM_ML(vdata), &tret, yqs);
    sunml_nvector_array_free(yqs); 
    SCHECK_FLAG("CVodeGetQuadSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_cvodes_quadsens_get1(value vdata, value vis, value vyqs)
{
    CAMLparam3(vdata, vis, vyqs);
    N_Vector yqs = NVEC_VAL(vyqs);
    sunrealtype tret;

    int flag = CVodeGetQuadSens1(CVODE_MEM_FROM_ML(vdata), &tret,
			         Int_val(vis), yqs);
    SCHECK_FLAG("CVodeGetQuadSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_cvodes_quadsens_get_dky(value vdata, value vt, value vk,
					 value vdkyqs)
{
    CAMLparam4(vdata, vt, vk, vdkyqs);
    N_Vector *dkyqs = sunml_nvector_array_alloc(vdkyqs);

    int flag = CVodeGetQuadSensDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
				   Int_val(vk), dkyqs);
    sunml_nvector_array_free(dkyqs); 
    SCHECK_FLAG("CVodeGetQuadSensDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_get_dky1(value vdata, value vt, value vk,
					  value vis, value vdkyqs)
{
    CAMLparam5(vdata, vt, vk, vis, vdkyqs);
    N_Vector dkyqs = NVEC_VAL(vdkyqs);

    int flag = CVodeGetQuadSensDky1(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
				    Int_val(vk), Int_val(vis), dkyqs);
    SCHECK_FLAG("CVodeGetQuadSensDky1", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_get_err_weights(value vdata, value veqweights)
{
    CAMLparam2(vdata, veqweights);
    N_Vector *eqweights = sunml_nvector_array_alloc(veqweights);

    int flag = CVodeGetQuadSensErrWeights(CVODE_MEM_FROM_ML(vdata), eqweights);
    sunml_nvector_array_free(eqweights); 
    SCHECK_FLAG("CVodeGetQuadSensErrWeights", flag);

    CAMLreturn (Val_unit);
}

/* adjoint interface */

static value forward_solver(value vdata, value vtout, value vyret, int onestep)
{
    CAMLparam3(vdata, vtout, vyret);
    CAMLlocal1(ret);
    N_Vector yret = NVEC_VAL(vyret);
    sunrealtype tret;
    int ncheck;
    enum cvode_solver_result_tag solver_result = -1;

    int flag = CVodeF(CVODE_MEM_FROM_ML(vdata), Double_val(vtout), yret,
		      &tret, onestep ? CV_ONE_STEP : CV_NORMAL, &ncheck);
    switch (flag) {
    case CV_SUCCESS:
	solver_result = VARIANT_CVODE_SOLVER_RESULT_SUCCESS;
	break;

    case CV_TSTOP_RETURN:
	solver_result = VARIANT_CVODE_SOLVER_RESULT_STOPTIMEREACHED;
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
	sunml_cvodes_check_flag("CVodeF", flag, CVODE_MEM_FROM_ML(vdata));
    }

    assert (Field (vdata, RECORD_CVODE_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple(3);
    Store_field(ret, 0, caml_copy_double(tret));
    Store_field(ret, 1, Val_int(ncheck));
    Store_field(ret, 2, Val_int(solver_result));

    CAMLreturn(ret);
}

CAMLprim value sunml_cvodes_adj_forward_normal(value vdata, value vtout, value vyret)
{
    CAMLparam3(vdata, vtout, vyret);
    CAMLreturn(forward_solver(vdata, vtout, vyret, 0));
}

CAMLprim value sunml_cvodes_adj_forward_one_step(value vdata, value vtout,
					    value vyret)
{
    CAMLparam3(vdata, vtout, vyret);
    CAMLreturn(forward_solver(vdata, vtout, vyret, 1));
}

CAMLprim value sunml_cvodes_adj_sv_tolerances(value vparent, value vwhich,
					  value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = CVodeSVtolerancesB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				  Double_val(vreltol), atol_nv);
    SCHECK_FLAG("CVodeSStolerancesB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_set_preconditioner(value vparent,
						     value vwhich,
						     value vset_precsetup,
						     value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_precsetup, vusesens);
    void *mem = CVODE_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    int flag;

#if 400 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vusesens)) {
	flag = CVodeSetPreconditionerBS(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn_withsens : NULL,
		    bprecsolvefn_withsens);
	SCHECK_FLAG ("CVodeSetPreconditionerBS", flag);
    } else {
	flag = CVodeSetPreconditionerB(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn : NULL,
		    bprecsolvefn);
	SCHECK_FLAG ("CVodeSetPreconditionerB", flag);
    }
#else
    if (Bool_val(vusesens)) {
#if 260 <= SUNDIALS_LIB_VERSION
	flag = CVSpilsSetPreconditionerBS(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn_withsens : NULL,
		    bprecsolvefn_withsens);
	SCHECK_FLAG ("CVSpilsSetPreconditionerBS", flag);
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    } else {
	flag = CVSpilsSetPreconditionerB(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn : NULL,
		    bprecsolvefn);
	SCHECK_FLAG ("CVSpilsSetPreconditionerB", flag);
    }
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_set_banded_preconditioner(value vparent,
							    value vwhich,
							    value vneqs,
							    value vmupper,
							    value vmlower)
{
    CAMLparam5(vparent, vwhich, vneqs, vmupper, vmlower);
    void *mem = CVODE_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    int flag = CVBandPrecInitB(mem, which, Index_val(vneqs),
			       Index_val(vmupper), Index_val(vmlower));
    SCHECK_FLAG ("CVBandPrecInitB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_jac_times_rhsfn(value vparent,
						    value vwhich,
						    value vhas_rhsfn)
{
    CAMLparam3(vparent, vwhich, vhas_rhsfn);
#if 530 <= SUNDIALS_LIB_VERSION
    void *mem = CVODE_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    CVRhsFn rhsfn = Bool_val (vhas_rhsfn) ? bjactimesrhsfn : NULL;

    int flag = CVodeSetJacTimesRhsFnB(mem, which, rhsfn);
    CHECK_LS_FLAG("CVodeSetJacTimesRhsFnB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION <= 301
/* Work around a bug in Sundials */
int CVSpilsSetJacTimesSetupFnBS(void *cvode_mem, int which,
                                CVSpilsJacTimesSetupFnBS jtsetupBS,
                                CVSpilsJacTimesVecFnBS jtimesBS);
#endif

CAMLprim value sunml_cvodes_adj_spils_set_jac_times(value vparent,
						value vwhich,
						value vhas_setup,
						value vhas_times,
						value vusesens)
{
    CAMLparam5(vparent, vwhich, vhas_setup, vhas_times, vusesens);
    void *mem = CVODE_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    int flag;

#if 400 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vusesens)) {
	flag = CVodeSetJacTimesBS(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn_withsens : NULL,
			Bool_val(vhas_times) ? bjactimesfn_withsens : NULL);
	SCHECK_FLAG ("CVodeSetJacTimesBS", flag);
    } else {
	flag = CVodeSetJacTimesB(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn : NULL,
			Bool_val(vhas_times) ? bjactimesfn : NULL);
	SCHECK_FLAG ("CVodeSetJacTimesB", flag);
    }
#elif 300 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vusesens)) {
#if 301 < SUNDIALS_LIB_VERSION
	flag = CVSpilsSetJacTimesBS(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn_withsens : NULL,
			Bool_val(vhas_times) ? bjactimesfn_withsens : NULL);
	SCHECK_FLAG ("CVSpilsSetJacTimesBS", flag);
#else
	/* Work around a bug in Sundials */
	flag = CVSpilsSetJacTimesSetupFnBS(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn_withsens : NULL,
			Bool_val(vhas_times) ? bjactimesfn_withsens : NULL);
	SCHECK_FLAG ("CVSpilsSetJacTimesSetupFnBS", flag);
#endif
    } else {
	flag = CVSpilsSetJacTimesB(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn : NULL,
			Bool_val(vhas_times) ? bjactimesfn : NULL);
	SCHECK_FLAG ("CVSpilsSetJacTimesB", flag);
    }
#else
    if (Bool_val(vusesens)) {
#if 260 <= SUNDIALS_LIB_VERSION
	flag = CVSpilsSetJacTimesVecFnBS(mem, which,
			Bool_val(vhas_times) ? bjactimesfn_withsens : NULL);
	SCHECK_FLAG ("CVSpilsSetJacTimesVecFnBS", flag);
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    } else {
	flag = CVSpilsSetJacTimesVecFnB(mem, which,
			Bool_val(vhas_times) ? bjactimesfn : NULL);
	SCHECK_FLAG ("CVSpilsSetJacTimesVecFnB", flag);
    }
#endif

    CAMLreturn (Val_unit);
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_cvodes_adj_dls_dense(value vparent, value vwhich,
				      value vnb, value vset_jac, value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_jac, vusesens);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    long nbeqs = Long_val(vnb);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVDenseB (cvode_mem, which, nbeqs);
    SCHECK_FLAG ("CVDenseB", flag);

    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	    flag = CVDlsSetDenseJacFnBS(cvode_mem, which, bjacfn_withsens);
	    SCHECK_FLAG("CVDlsSetDenseJacFnBS", flag);
#else
	    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	} else {
	    flag = CVDlsSetDenseJacFnB(cvode_mem, which, bjacfn_nosens);
	    SCHECK_FLAG("CVDlsSetDenseJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_dls_lapack_dense(value vparent, value vwhich,
					     value vnb, value vset_jac,
					     value vusesens)
{
    CAMLparam3 (vparent, vwhich, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    long nbeqs = Long_val(vnb);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVLapackDenseB (cvode_mem, which, nbeqs);
    SCHECK_FLAG ("CVLapackDenseB", flag);

    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	    flag = CVDlsSetDenseJacFnBS(cvode_mem, which, bjacfn_withsens);
	    SCHECK_FLAG("CVDlsSetDenseJacFnBS", flag);
#else
	    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	} else {
	    flag = CVDlsSetDenseJacFnB(cvode_mem, which, bjacfn_nosens);
	    SCHECK_FLAG("CVDlsSetDenseJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_dls_band (value vparent_which, value vsizes,
				      value vset_jac, value vusesens)
{
    CAMLparam4(vparent_which, vsizes, vset_jac, vusesens);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(Field(vsizes, 0));
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVBandB (cvode_mem, which, nbeqs, Long_val (Field(vsizes, 1)),
		    Long_val (Field(vsizes, 2)));
    SCHECK_FLAG ("CVBandB", flag);

    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	    flag = CVDlsSetBandJacFnBS(cvode_mem, which, bbandjacfn_withsens);
	    SCHECK_FLAG("CVDlsSetBandJacFnBS", flag);
#else
	    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	} else {
	    flag = CVDlsSetBandJacFnB(cvode_mem, which, bbandjacfn_nosens);
	    SCHECK_FLAG("CVDlsSetBandJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_dls_lapack_band (value vparent_which, value vsizes,
					     value vset_jac, value vusesens)
{
    CAMLparam4(vparent_which, vsizes, vset_jac, vusesens);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(Field(vsizes, 0));
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVLapackBandB (cvode_mem, which, nbeqs,
			  Long_val (Field(vsizes, 1)),
			  Long_val (Field(vsizes, 2)));
    SCHECK_FLAG ("CVLapackBandB", flag);

    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	    flag = CVDlsSetBandJacFnBS(cvode_mem, which, bbandjacfn_withsens);
	    SCHECK_FLAG("CVDlsSetBandJacFnBS", flag);
#else
	    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	} else {
	    flag = CVDlsSetBandJacFnB(cvode_mem, which, bbandjacfn_nosens);
	    SCHECK_FLAG("CVDlsSetBandJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_linear_solver (value vparent_which,
		    value vlsolv, value vojmat, value vjac, value vlsf)
{
    CAMLparam5(vparent_which, vlsolv, vojmat, vjac, vlsf);
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal4(vhasjac, vjacsens, vhaslsf, vlsfsens);

    vhasjac = Field(vjac, 0);
    vjacsens = Field(vjac, 1);

    vhaslsf = Field(vlsf, 0);
    vlsfsens = Field(vlsf, 1);

    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    int which = Int_val(Field(vparent_which, 1));
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = CVodeSetLinearSolverB(cvode_mem, which, lsolv, jmat);
    CHECK_LS_FLAG ("CVodeSetLinearSolverB", flag);

    if (Bool_val (vjacsens)) {
	flag = CVodeSetJacFnBS(cvode_mem, which,
			       Bool_val (vhasjac) ? bjacfn_withsens : NULL);
	CHECK_LS_FLAG("CVodeSetJacFnBS", flag);
    } else {
	flag = CVodeSetJacFnB(cvode_mem, which,
			       Bool_val (vhasjac) ? bjacfn_nosens : NULL);
	CHECK_LS_FLAG("CVodeSetJacFnB", flag);
    }

#if 500 <= SUNDIALS_LIB_VERSION
    if (Bool_val (vhaslsf)) {
	if (Bool_val (vlsfsens)) {
	    CVodeSetLinSysFnBS(cvode_mem, which, blinsysfn_withsens);
	    CHECK_LS_FLAG("CVodeSetLinSysFnBS", flag);
	} else {
	    CVodeSetLinSysFnB(cvode_mem, which, blinsysfn_nosens);
	    CHECK_LS_FLAG("CVodeSetLinSysFnB", flag);
	}
    }
#endif

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_dls_set_linear_solver (value vparent_which,
						   value vlsolv, value vjmat,
						   value vhasjac, value vusesens)
{
    CAMLparam5(vparent_which, vlsolv, vjmat, vhasjac, vusesens);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    int which = Int_val(Field(vparent_which, 1));
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVDlsSetLinearSolverB(cvode_mem, which, lsolv, jmat);
    CHECK_FLAG ("CVDlsSetLinearSolverB", flag);

    if (Bool_val (vusesens)) {
	flag = CVDlsSetJacFnBS(cvode_mem, which,
			       Bool_val (vhasjac) ? bjacfn_withsens : NULL);
	SCHECK_FLAG("CVDlsSetJacFnBS", flag);
    } else {
	flag = CVDlsSetJacFnB(cvode_mem, which,
			      Bool_val (vhasjac) ? bjacfn_nosens : NULL);
	SCHECK_FLAG("CVDlsSetJacFnB", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_set_linear_solver (value vparent, value vwhich,
						     value vlsolv)
{
    CAMLparam3(vparent, vwhich, vlsolv);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpilsSetLinearSolverB(cvode_mem, which, lsolv);
    CHECK_FLAG ("CVSpilsSetLinearSolverB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_init_backward(value vparent, value weakref,
					      value vargs, value vwithsens)
{
    CAMLparam4(vparent, weakref, vargs, vwithsens);
    CAMLlocal2(r, vcvode_mem);
    CAMLlocal2(vlmm, viter);
    int flag, lmm_c, which;
    void *parent = CVODE_MEM_FROM_ML(vparent);

    int lmm = Field(vargs, 0);
#if SUNDIALS_LIB_VERSION < 400
    int iter = Field(vargs, 1);
#endif
    sunrealtype tb0 = Double_val(Field(vargs, 2));
    N_Vector initial_nv = NVEC_VAL(Field(vargs, 3));

    switch (Int_val(lmm)) {
    case VARIANT_CVODE_LMM_ADAMS:
	lmm_c = CV_ADAMS;
	break;

    case VARIANT_CVODE_LMM_BDF:
	lmm_c = CV_BDF;
	break;

    default:
	caml_failwith("Illegal lmm value.");
    }

#if 400 <= SUNDIALS_LIB_VERSION
    flag = CVodeCreateB(parent, lmm_c, &which);
#else
    flag = CVodeCreateB(parent, lmm_c,
			Bool_val(iter) ? CV_NEWTON : CV_FUNCTIONAL,
			&which);
#endif
    if (flag != CV_SUCCESS) {
	SCHECK_FLAG("CVodeCreateB", flag);
    }

    if (Bool_val(vwithsens)) {
	flag = CVodeInitBS(parent, which, brhsfn_sens, tb0, initial_nv);
	if (flag != CV_SUCCESS) {
	    SCHECK_FLAG("CVodeInitBS", flag);
	}
    } else {
	flag = CVodeInitB(parent, which, brhsfn, tb0, initial_nv);
	if (flag != CV_SUCCESS) {
	    SCHECK_FLAG("CVodeInitB", flag);
	}
    }

    vcvode_mem = sunml_wrap_session_pointer(CVodeGetAdjCVodeBmem(parent, which));

    value *backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	caml_raise_out_of_memory();
    }
    CVodeSetUserDataB (parent, which, backref);

    r = caml_alloc_tuple (3);
    Store_field (r, 0, vcvode_mem);
    Store_field (r, 1, Val_int(which));
    Store_field (r, 2, (value)backref);

    CAMLreturn(r);
}

CAMLprim value sunml_cvodes_adj_reinit(value vparent, value vwhich,
				   value vtb0, value vyb0)
{
    CAMLparam4(vparent, vwhich, vtb0, vyb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yb0 = NVEC_VAL(vyb0);

    flag = CVodeReInitB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			Double_val(vtb0), yb0);
    SCHECK_FLAG("CVodeReInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_get(value vparent, value vwhich, value vyb)
{
    CAMLparam3(vparent, vwhich, vyb);
    N_Vector yq = NVEC_VAL(vyb);
    sunrealtype tret;

    int flag = CVodeGetB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			 &tret, yq);
    SCHECK_FLAG("CVodeGetB", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_cvodes_adj_get_y(value vdata, value vt, value vy)
{
    CAMLparam3(vdata, vt, vy);

    N_Vector y_nv = NVEC_VAL(vy);

    int flag = CVodeGetAdjY(CVODE_MEM_FROM_ML(vdata), Double_val(vt), y_nv);
    SCHECK_FLAG("CVodeGetAdjY", flag);

    CAMLreturn (Val_unit);
}

/* adjoint/quadrature interface */

CAMLprim value sunml_cvodes_adjquad_initb(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);

    flag = CVodeQuadInitB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			  bquadrhsfn, yqb0);
    SCHECK_FLAG("CVodeQuadInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adjquad_initbs(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = CVodeQuadInitBS(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			   bquadrhsfn_sens, yqb0);
    SCHECK_FLAG("CVodeQuadInitBS", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adjquad_reinit(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = CVodeQuadReInitB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich), yqb0);
    SCHECK_FLAG("CVodeQuadReInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adjquad_get(value vparent, value vwhich, value vyqb)
{
    CAMLparam3(vparent, vwhich, vyqb);
    CAMLlocal1(r);
    N_Vector yqb = NVEC_VAL(vyqb);
    sunrealtype tret;

    int flag = CVodeGetQuadB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			     &tret, yqb);
    SCHECK_FLAG("CVodeGetQuadB", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_cvodes_adjquad_sv_tolerances(value vparent, value vwhich,
					      value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = CVodeQuadSVtolerancesB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				      Double_val(vreltol), atol_nv);
    SCHECK_FLAG("CVodeQuadSVtolerancesB", flag);

    CAMLreturn (Val_unit);
}

/* Interface without nvectors */

#define MAX_ERRMSG_LEN 256

void sunml_cvodes_check_flag(const char *call, int flag, void *cvode_mem)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	    || flag == CV_ROOT_RETURN
	    || flag == CV_TSTOP_RETURN
	    || flag == CV_WARNING) return;

    switch (flag) {
	case CV_TOO_MUCH_WORK:
	    caml_raise_constant(CVODE_EXN(TooMuchWork));

	case CV_TOO_MUCH_ACC:
	    caml_raise_constant(CVODE_EXN(TooMuchAccuracy));

	case CV_ERR_FAILURE:
	    caml_raise_constant(CVODE_EXN(ErrFailure));

	case CV_CONV_FAILURE:
	    caml_raise_constant(CVODE_EXN(ConvergenceFailure));

	/* * */

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

	/* * */

	case CV_ILL_INPUT:
	    caml_raise_constant(CVODE_EXN(IllInput));
	
	case CV_BAD_K:
	    caml_raise_constant(CVODE_EXN(BadK));

	case CV_BAD_T:
	    caml_raise_constant(CVODE_EXN(BadT));

	case CV_TOO_CLOSE:
	    caml_raise_constant(CVODE_EXN(TooClose));

	/* Quadrature */

        case CV_NO_QUAD:
	    caml_raise_constant(CVODES_EXN(QuadNotInitialized));

        case CV_QRHSFUNC_FAIL:
	    caml_raise_constant(CVODES_EXN(QuadRhsFuncFailure));

        case CV_FIRST_QRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(FirstQuadRhsFuncFailure));

        case CV_REPTD_QRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(RepeatedQuadRhsFuncFailure));

        case CV_UNREC_QRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(UnrecoverableQuadRhsFuncFailure));

	/* Sensitivity */

        case CV_NO_SENS:
	    caml_raise_constant(CVODES_EXN(SensNotInitialized));

        case CV_SRHSFUNC_FAIL:
	    caml_raise_constant(CVODES_EXN(SensRhsFuncFailure));

        case CV_FIRST_SRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(FirstSensRhsFuncFailure));

        case CV_REPTD_SRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(RepeatedSensRhsFuncFailure));

        case CV_UNREC_SRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(UnrecoverableSensRhsFuncFailure));

        case CV_BAD_IS:
	    caml_raise_constant(CVODES_EXN(BadSensIdentifier));

	/* Sensitivity > Quadrature */

        case CV_NO_QUADSENS:
	    caml_raise_constant(CVODES_EXN(QuadSensNotInitialized));

        case CV_QSRHSFUNC_FAIL:
	    caml_raise_constant(CVODES_EXN(QuadSensRhsFuncFailure));

        case CV_FIRST_QSRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(FirstQuadSensRhsFuncFailure));

        case CV_REPTD_QSRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(RepeatedQuadSensRhsFuncFailure));

        case CV_UNREC_QSRHSFUNC_ERR:
	    caml_raise_constant(CVODES_EXN(UnrecoverableQuadSensRhsFuncFailure));

	/* Adjoint */

        case CV_NO_ADJ:
	    caml_raise_constant(CVODES_EXN(AdjointNotInitialized));

        case CV_NO_FWD:
	    caml_raise_constant(CVODES_EXN(NoForwardCall));

        case CV_NO_BCK:
	    caml_raise_constant(CVODES_EXN(NoBackwardProblem));

        case CV_BAD_TB0:
	    caml_raise_constant(CVODES_EXN(BadFinalTime));

        case CV_REIFWD_FAIL:
	    caml_raise_constant(CVODES_EXN(ForwardReinitFailure));

        case CV_FWD_FAIL:
	    caml_raise_constant(CVODES_EXN(ForwardFailure));

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

/* basic quadrature interface */

CAMLprim value sunml_cvodes_quad_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = CVodeSetQuadErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quad_ss_tolerances(value vdata,
					   value reltol,
					   value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeQuadSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    SCHECK_FLAG("CVodeQuadSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quad_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_quad_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_quad_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfqevals;
    long int nqetfails;

    flag = CVodeGetQuadStats(CVODE_MEM_FROM_ML(vdata), &nfqevals,
	    &nqetfails);
    SCHECK_FLAG("CVodeGetQuadStats", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(nfqevals));
    Store_field(r, 1, Val_long(nqetfails));

    CAMLreturn(r);
}

/* sensitivity interface */

CAMLprim value sunml_cvodes_sens_set_err_con(value vdata, value verrcons)
{
    CAMLparam2(vdata, verrcons);
    int flag;
    
    flag = CVodeSetSensErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrcons));
    SCHECK_FLAG("CVodeSetSensErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_ss_tolerances(value vdata,
					   value reltol,
					   value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeSensSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), REAL_ARRAY(abstol));
    SCHECK_FLAG("CVodeSensSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSensEEtolerances(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSensEEtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_set_params(value vdata, value vparams)
{
    CAMLparam2(vdata, vparams);
    CAMLlocal3(vp, vpbar, vplist);

    sunrealtype *p = NULL;
    sunrealtype *pbar = NULL;
    int *plist = NULL;
    int i, ns;

    vp = Field(vparams, RECORD_CVODES_SENS_PARAMS_PVALS);
    vpbar = Field(vparams, RECORD_CVODES_SENS_PARAMS_PBAR);
    vplist = Field(vparams, RECORD_CVODES_SENS_PARAMS_PLIST);

    if (vp != Val_none) p = REAL_ARRAY(Some_val(vp));
    if (vpbar != Val_none) pbar = REAL_ARRAY(Some_val(vpbar));

    if (vplist != Val_none) {
	vplist = Some_val(vplist);
	ns = Wosize_val (vplist); /* vplist : int array */
	plist = calloc(ns, sizeof(*plist));

	for (i=0; i < ns; ++i) {
	    plist[i] = Int_val(Field(vplist, i));
	}
    }

    int flag = CVodeSetSensParams(CVODE_MEM_FROM_ML(vdata), p, pbar, plist);
    if (plist != NULL) free(plist);
    SCHECK_FLAG("CVodeSetSensParams", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_toggle_off(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSensToggleOff(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSensToggleOff", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_set_dq_method(value vdata, value vdqtype,
					   value vdqrhomax)
{
    CAMLparam3(vdata, vdqtype, vdqrhomax);
    int dqtype;

    switch (Int_val(vdqtype)) {
    case VARIANT_CVODES_SENS_DQ_METHOD_CENTERED:
	dqtype = CV_CENTERED;
	break;

    case VARIANT_CVODES_SENS_DQ_METHOD_FORWARD:
	dqtype = CV_FORWARD;
	break;

    default:
	caml_failwith("Illegal dq method.");
    }

    int flag = CVodeSetSensDQMethod(CVODE_MEM_FROM_ML(vdata), dqtype,
				    Double_val(vdqrhomax));
    SCHECK_FLAG("CVodeSetSensMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_set_max_nonlin_iters(value vdata, value vmaxcors)
{
    CAMLparam2(vdata, vmaxcors);

    int flag = CVodeSetSensMaxNonlinIters(CVODE_MEM_FROM_ML(vdata),
	    Int_val(vmaxcors));
    SCHECK_FLAG("CVodeSetSensMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_sens_get_num_rhs_evals_sens(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetNumRhsEvalsSens(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetNumRhsEvalsSens", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_sens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_sens_get_num_step_solve_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    long int v;

#if 620 <= SUNDIALS_LIB_VERSION
    int flag;
    flag = CVodeGetNumStepSensSolveFails(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumStepSensSolveFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim void sunml_cvodes_sens_get_num_step_stgr_solve_fails(value vcvode_mem,
							      value vnfails)
{
    CAMLparam2(vcvode_mem, vnfails);

#if 620 <= SUNDIALS_LIB_VERSION
    int flag;
    flag = CVodeGetNumStepStgrSensSolveFails(CVODE_MEM_FROM_ML(vcvode_mem),
					     LONG_ARRAY(vnfails));
    CHECK_FLAG("CVodeGetNumStepStgrSensSolveFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn0;
}

CAMLprim value sunml_cvodes_sens_get_num_lin_solv_setups(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumLinSolvSetups(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_sens_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfsevals;
    long int nfevalss;
    long int nsetfails;
    long int nlinsetupss;

    flag = CVodeGetSensStats(CVODE_MEM_FROM_ML(vdata), &nfsevals,
	    &nfevalss, &nsetfails, &nlinsetupss);
    SCHECK_FLAG("CVodeGetSensStats", flag);

    r = caml_alloc_tuple(RECORD_CVODES_SENS_STATS_SIZE);
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_SENS_EVALS, Val_long(nfsevals));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_RHS_EVALS, Val_long(nfevalss));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_ERR_TEST_FAILS,
							   Val_long(nsetfails));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_LIN_SOLV_SETUPS,
							 Val_long(nlinsetupss));

    CAMLreturn(r);
}

CAMLprim value sunml_cvodes_sens_get_num_nonlin_solv_iters(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumNonlinSolvIters(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_sens_get_num_nonlin_solv_conv_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_sens_get_nonlin_solv_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nsniters;
    long int nsncfails;

    flag = CVodeGetSensNonlinSolvStats(CVODE_MEM_FROM_ML(vdata), &nsniters,
	    &nsncfails);
    SCHECK_FLAG("CVodeGetSensNonlinSolvStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nsniters));
    Store_field(r, 1, Val_long(nsncfails));

    CAMLreturn(r);
}

CAMLprim value sunml_cvodes_sens_get_num_stgr_nonlin_solv_iters(value vdata,
							    value vr)
{
    CAMLparam2(vdata, vr);

    int flag;

    flag = CVodeGetStgrSensNumNonlinSolvIters(CVODE_MEM_FROM_ML(vdata),
					      LONG_ARRAY(vr));
    SCHECK_FLAG("CVodeGetStgrSensNumNonlinSolvIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_sens_get_num_stgr_nonlin_solv_conv_fails(value vdata,
								 value vr)
{
    CAMLparam2(vdata, vr);

    int flag = CVodeGetStgrSensNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vdata),
						      LONG_ARRAY(vr));
    SCHECK_FLAG("CVodeGetStgrSensNumNonlinSolvConvFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_get_current_state_sens(value vdata, value vns)
{
    CAMLparam2(vdata, vns);
    CAMLlocal1(r);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector *yS;

    int flag = CVodeGetCurrentStateSens(CVODE_MEM_FROM_ML(vdata), &yS);
    SCHECK_FLAG("CVodeGetCurrentStateSens", flag);

    r = sunml_wrap_to_nvector_table(Int_val(vns), yS);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (r);
}

CAMLprim value sunml_cvodes_get_current_sens_solve_index(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

#if 500 <= SUNDIALS_LIB_VERSION
    int index;

    int flag = CVodeGetCurrentSensSolveIndex(CVODE_MEM_FROM_ML(vdata), &index);
    SCHECK_FLAG("CVodeGetCurrentSensSolveIndex", flag);

    r = Val_int(index);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (r);
}

CAMLprim value sunml_cvodes_get_nonlin_system_data_sens(value vcvode_mem,
							value vns)
{
    CAMLparam2(vcvode_mem, vns);
    CAMLlocal1(vnv);
#if 540 <= SUNDIALS_LIB_VERSION
    int ns = Int_val(vns);
    sunrealtype tn, gamma, rlS1;
    N_Vector *ySpred, *ySn, *znS1;
    void *user_data;

    int flag = CVodeGetNonlinearSystemDataSens(CVODE_MEM_FROM_ML(vcvode_mem),
		    &tn, &ySpred, &ySn, &gamma, &rlS1, &znS1, &user_data);
    CHECK_FLAG("CVodeGetNonlinearSystemDataSens", flag);

    vnv = caml_alloc_tuple(RECORD_CVODES_NONLIN_SYSTEM_DATA_SIZE);
    Store_field(vnv, RECORD_CVODES_NONLIN_SYSTEM_DATA_TN,
	    caml_copy_double(tn));
    Store_field(vnv, RECORD_CVODES_NONLIN_SYSTEM_DATA_YSPRED,
	    sunml_wrap_to_nvector_table(ns, ySpred));
    Store_field(vnv, RECORD_CVODES_NONLIN_SYSTEM_DATA_YSN,
	    sunml_wrap_to_nvector_table(ns, ySn));
    Store_field(vnv, RECORD_CVODES_NONLIN_SYSTEM_DATA_GAMMA,
	    caml_copy_double(gamma));
    Store_field(vnv, RECORD_CVODES_NONLIN_SYSTEM_DATA_RLS1,
	    caml_copy_double(rlS1));
    Store_field(vnv, RECORD_CVODES_NONLIN_SYSTEM_DATA_ZNS1,
	    sunml_wrap_to_nvector_table(ns, znS1));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim void sunml_cvodes_compute_state_sens(value vcvode_mem,
					       value vyscor, value vysz)
{
    CAMLparam3(vcvode_mem, vyscor, vysz);
#if 540 <= SUNDIALS_LIB_VERSION
    N_Vector *yscor, *ysz;
    int flag;

    yscor = sunml_nvector_array_alloc(vyscor);
    if (yscor == NULL) caml_raise_out_of_memory();
    ysz = sunml_nvector_array_alloc(vysz);
    if (ysz == NULL) {
	sunml_nvector_array_free(yscor);
	caml_raise_out_of_memory();
    }

    flag = CVodeComputeStateSens(CVODE_MEM_FROM_ML(vcvode_mem), yscor, ysz);
    sunml_nvector_array_free(yscor);
    sunml_nvector_array_free(ysz);
    CHECK_FLAG("CVodeComputeStateSens", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_cvodes_compute_state_sens1(value vcvode_mem, value vidx,
						value vycor, value vyz)
{
    CAMLparam4(vcvode_mem, vidx, vycor, vyz);
#if 540 <= SUNDIALS_LIB_VERSION
    int flag = CVodeComputeStateSens1(CVODE_MEM_FROM_ML(vcvode_mem),
				      Int_val(vidx),
				      NVEC_VAL(vycor), NVEC_VAL(vyz));
    CHECK_FLAG("CVodeComputeStateSens1", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

/* sensitivity/quadrature interface */

CAMLprim value sunml_cvodes_quadsens_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = CVodeSetQuadSensErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadSensErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_ss_tolerances(value vdata,
					       value reltol,
					       value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeQuadSensSStolerances(CVODE_MEM_FROM_ML(vdata),
					 Double_val(reltol),
					 REAL_ARRAY(abstol));
    SCHECK_FLAG("CVodeQuadSensSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeQuadSensEEtolerances(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeQuadSensEEtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_quadsens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadSensNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_quadsens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadSensNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_cvodes_quadsens_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfqsevals;
    long int nqsetfails;

    flag = CVodeGetQuadSensStats(CVODE_MEM_FROM_ML(vdata), &nfqsevals,
				 &nqsetfails);
    SCHECK_FLAG("CVodeGetQuadSensStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nfqsevals));
    Store_field(r, 1, Val_long(nqsetfails));

    CAMLreturn(r);
}

/* adjoint interface */

CAMLprim value sunml_cvodes_adj_set_nonlinear_solver(value vparent,
						     value vwhich, value vnlsolv)
{
    CAMLparam3(vparent, vwhich, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);
    int flag;

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_cvodes_bsession_to_value,
				   sunml_cvodes_session_from_value);

    flag = CVodeSetNonlinearSolverB(cvode_mem, Int_val(vwhich), nlsolv);
    CHECK_FLAG ("CVodeSetNonlinearSolverB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_init(value vdata, value vnd, value vinterptype)
{
    CAMLparam3(vdata, vnd, vinterptype);
    int interptype;

    switch(Int_val(vinterptype)) {
    case VARIANT_CVODES_ADJ_INTERPOLATION_POLYNOMIAL:
	interptype = CV_POLYNOMIAL;
	break;

    case VARIANT_CVODES_ADJ_INTERPOLATION_HERMITE:
	interptype = CV_HERMITE;
	break;

    default:
	caml_failwith("Illegal interpolation value.");
    }

    int flag = CVodeAdjInit(CVODE_MEM_FROM_ML(vdata), Long_val(vnd),
			    interptype);
    SCHECK_FLAG("CVodeAdjInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_ss_tolerances(value vparent, value vwhich,
					  value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);

    int flag = CVodeSStolerancesB(CVODE_MEM_FROM_ML(vparent),
		 Int_val(vwhich), Double_val(vreltol), Double_val(vabstol));
    SCHECK_FLAG("CVodeSStolerancesB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_diag(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);

    int flag = CVDiagB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich));
    SCHECK_FLAG("CVDiagB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_spgmr (value vparent, value vwhich,
					 value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpgmrB (cvode_mem, which, sunml_lsolver_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpgmrB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_spbcgs (value vparent, value vwhich,
					  value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpbcgB (cvode_mem, which, sunml_lsolver_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpbcgB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_sptfqmr (value vparent, value vwhich,
					   value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSptfqmrB (cvode_mem, which, sunml_lsolver_precond_type (vtype),
		       Int_val (vmaxl));
    SCHECK_FLAG ("CVSptfqmrB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_functional (value vparent, value vwhich)
{
    CAMLparam2 (vparent, vwhich);
#if SUNDIALS_LIB_VERSION < 400
    int flag = CVodeSetIterTypeB (CVODE_MEM_FROM_ML (vparent), Int_val(vwhich),
				  CV_FUNCTIONAL);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_newton (value vparent, value vwhich)
{
    CAMLparam2 (vparent, vwhich);
#if SUNDIALS_LIB_VERSION < 400
    int flag = CVodeSetIterTypeB (CVODE_MEM_FROM_ML (vparent), Int_val(vwhich),
				  CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_bsession_finalize(value vdata)
{
    if (CVODE_MEM_FROM_ML(vdata) != NULL) {
	value *backref = CVODE_BACKREF_FROM_ML(vdata);
	// NB: CVodeFree() is *not* called: parents free-up backward problems
	sunml_sundials_free_value(backref);
    }
    return Val_unit;
}

CAMLprim value sunml_cvodes_adj_backward_normal(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = CVodeB(CVODE_MEM_FROM_ML(vdata), Double_val(vtbout), CV_NORMAL);
    SCHECK_FLAG("CVodeB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_backward_one_step(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = CVodeB(CVODE_MEM_FROM_ML(vdata), Double_val(vtbout),
		      CV_ONE_STEP);
    SCHECK_FLAG("CVodeB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_no_sensitivity(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSetAdjNoSensi(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSetAdjNoSensi", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_max_ord(value vparent, value vwhich,
					value vmaxord)
{
    CAMLparam3(vparent, vwhich, vmaxord);

    int flag = CVodeSetMaxOrdB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(vmaxord));
    SCHECK_FLAG("CVodeSetMaxOrdB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_max_num_steps(value vparent, value vwhich,
					      value vmxsteps)
{
    CAMLparam3(vparent, vwhich, vmxsteps);

    int flag = CVodeSetMaxNumStepsB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				    Long_val(vmxsteps));
    SCHECK_FLAG("CVodeSetMaxNumStepsB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_init_step(value vparent, value vwhich, value vhin)
{
    CAMLparam3(vparent, vwhich, vhin);

    int flag = CVodeSetInitStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			         Double_val(vhin));
    SCHECK_FLAG("CVodeSetInitStepB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_min_step(value vparent, value vwhich, value vhmin)
{
    CAMLparam3(vparent, vwhich, vhmin);

    int flag = CVodeSetMinStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			        Double_val(vhmin));
    SCHECK_FLAG("CVodeSetMinStepB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_max_step(value vparent, value vwhich, value vhmax)
{
    CAMLparam3(vparent, vwhich, vhmax);

    int flag = CVodeSetMaxStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			        Double_val(vhmax));
    SCHECK_FLAG("CVodeSetMaxStepB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_stab_lim_det(value vparent, value vwhich,
					     value vstldet)
{
    CAMLparam3(vparent, vwhich, vstldet);

    int flag = CVodeSetStabLimDetB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			           Bool_val(vstldet));
    SCHECK_FLAG("CVodeSetStabLimDetB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_constraints(value vparent, value vwhich,
					        value vconstraints)
{
    CAMLparam3(vparent, vwhich, vconstraints);
#if SUNDIALS_LIB_VERSION >= 320
    N_Vector constraints = NVEC_VAL (vconstraints);

    int flag = CVodeSetConstraintsB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			            constraints);
    SCHECK_FLAG("CVodeSetConstraintsB", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_clear_constraints(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);
#if SUNDIALS_LIB_VERSION >= 320
    int flag = CVodeSetConstraintsB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			            NULL);
    SCHECK_FLAG("CVodeSetConstraintsB", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_set_prec_type(value vparent, value vwhich,
						value vptype)
{
    CAMLparam3(vparent, vwhich, vptype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = CVSpilsSetPrecTypeB(CVODE_MEM_FROM_ML(vparent),
				   Int_val(vwhich), sunml_lsolver_precond_type(vptype));
    SCHECK_FLAG("CVSpilsSetPrecTypeB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_set_gs_type(value vparent, value vwhich,
					      value vgstype)
{
    CAMLparam3(vparent, vwhich, vgstype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = CVSpilsSetGSTypeB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				 sunml_lsolver_gs_type(vgstype));
    SCHECK_FLAG("CVSpilsSetGSTypeB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

#if SUNDIALS_LIB_VERSION <= 250
SUNDIALS_EXPORT int CVSpilsSetEpsLinB(void *cvode_mem, int which,
				      sunrealtype eplifacB);
#endif

CAMLprim value sunml_cvodes_adj_set_eps_lin(value vparent, value vwhich,
					    value eplifac)
{
    CAMLparam3(vparent, vwhich, eplifac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetEpsLinB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Double_val(eplifac));
    SCHECK_FLAG("CVodeSetEpsLinB", flag);
#else
    int flag = CVSpilsSetEpsLinB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				 Double_val(eplifac));
    SCHECK_FLAG("CVSpilsSetEpsLinB", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_ls_norm_factor(value vparent, value vwhich,
					           value vfac)
{
    CAMLparam3(vparent, vwhich, vfac);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetLSNormFactorB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			             Double_val(vfac));
    SCHECK_FLAG("CVodeSetLSNormFactorB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_set_linear_solution_scaling(value vparent,
							    value vwhich,
							    value vonoff)
{
    CAMLparam3(vparent, vwhich, vonoff);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = CVodeSetLinearSolutionScalingB(CVODE_MEM_FROM_ML(vparent),
					      Int_val(vwhich),
					      Bool_val(vonoff));
    SCHECK_FLAG("CVodeSetLinearSolutionScalingB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adj_spils_set_maxl(value vparent, value vwhich,
					   value maxl)
{
    CAMLparam3(vparent, vwhich, maxl);
#if SUNDIALS_LIB_VERSION < 300
    int flag = CVSpilsSetMaxlB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(maxl));
    SCHECK_FLAG("CVSpilsSetMaxlB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* adjoint/quadrature interface */

CAMLprim value sunml_cvodes_adjquad_set_err_con(value vparent, value vwhich,
					    value verrconq)
{
    CAMLparam3(vparent, vwhich, verrconq);
    int flag;
    
    flag = CVodeSetQuadErrConB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadErrConB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_adjquad_ss_tolerances(value vparent, value vwhich,
					      value reltol, value abstol)
{
    CAMLparam4(vparent, vwhich, reltol, abstol);

    int flag = CVodeQuadSStolerancesB(CVODE_MEM_FROM_ML(vparent),
				      Int_val(vwhich), Double_val(reltol),
				      Double_val(abstol));
    SCHECK_FLAG("CVodeQuadSStolerancesB", flag);

    CAMLreturn (Val_unit);
}

