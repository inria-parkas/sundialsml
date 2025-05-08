/***********************************************************************
 *                                                                     *
 *               OCaml interface to (serial) Sundials                  *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

/* Interface functions that do not involve NVectors. */

#include "../config.h"
#include <idas/idas.h>
#include <sundials/sundials_band.h>

#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

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
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#endif

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <idas/idas_lapack.h>
#endif

#include "../ida/ida_ml.h"
#include "idas_ml.h"
#include "../sundials/sundials_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../lsolvers/sundials_linearsolver_ml.h"
#include "../lsolvers/sundials_nonlinearsolver_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#define MAX_ERRMSG_LEN 256

CAMLprim value sunml_idas_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (IDAS, exns);
    CAMLreturn (Val_unit);
}

/* Interface with nvectors */

CAMLprim value sunml_idas_alloc_nvector_array(value vn)
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

void sunml_idas_check_flag(const char *call, int flag, void *ida_mem)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == IDA_SUCCESS
	    || flag == IDA_ROOT_RETURN
	    || flag == IDA_TSTOP_RETURN
	    || flag == IDA_WARNING) return;

    switch (flag) {
    case IDA_TOO_MUCH_WORK:
	caml_raise_constant(IDA_EXN(TooMuchWork));

    case IDA_TOO_MUCH_ACC:
	caml_raise_constant(IDA_EXN(TooMuchAccuracy));

    case IDA_ERR_FAIL:
	caml_raise_constant(IDA_EXN(ErrFailure));

    case IDA_CONV_FAIL:
	caml_raise_constant(IDA_EXN(ConvergenceFailure));

    case IDA_LINIT_FAIL:
	caml_raise_constant(IDA_EXN(LinearInitFailure));

    case IDA_LSETUP_FAIL:
	caml_raise_with_arg(IDA_EXN(LinearSetupFailure),
			    sunml_ida_last_lin_exception(ida_mem));

    case IDA_LSOLVE_FAIL:
	caml_raise_with_arg(IDA_EXN(LinearSolveFailure),
			    sunml_ida_last_lin_exception(ida_mem));

    case IDA_RES_FAIL:
	caml_raise_constant(IDA_EXN(ResFuncFailure));

    case IDA_REP_RES_ERR:
	caml_raise_constant(IDA_EXN(RepeatedResFuncFailure));

    case IDA_RTFUNC_FAIL:
	caml_raise_constant(IDA_EXN(RootFuncFailure));

    case IDA_CONSTR_FAIL:
	caml_raise_constant(IDA_EXN(ConstraintFailure));

    case IDA_FIRST_RES_FAIL:
	caml_raise_constant(IDA_EXN(FirstResFuncFailure));

    case IDA_LINESEARCH_FAIL:
	caml_raise_constant(IDA_EXN(LinesearchFailure));

    case IDA_NO_RECOVERY:
	caml_raise_constant(IDA_EXN(NoRecovery));

    case IDA_ILL_INPUT:
	caml_raise_constant(IDA_EXN(IllInput));

    case IDA_BAD_EWT:
	caml_raise_constant(IDA_EXN(BadEwt));

    case IDA_BAD_K:
	caml_raise_constant(IDA_EXN(BadK));

    case IDA_BAD_T:
	caml_raise_constant(IDA_EXN(BadT));

    case IDA_NO_QUAD:
	caml_raise_constant(IDAS_EXN(QuadNotInitialized));

    case IDA_QRHS_FAIL:
	caml_raise_constant(IDAS_EXN(QuadRhsFuncFailure));

    case IDA_FIRST_QRHS_ERR:
	caml_raise_constant(IDAS_EXN(FirstQuadRhsFuncFailure));

    case IDA_REP_QRHS_ERR:
	caml_raise_constant(IDAS_EXN(RepeatedQuadRhsFuncFailure));

    case IDA_NO_SENS:
	caml_raise_constant(IDAS_EXN(SensNotInitialized));

    case IDA_SRES_FAIL:
	caml_raise_constant(IDAS_EXN(SensResFuncFailure));

    case IDA_REP_SRES_ERR:
	caml_raise_constant(IDAS_EXN(RepeatedSensResFuncFailure));

    case IDA_BAD_IS:
	caml_raise_constant(IDAS_EXN(BadSensIdentifier));

    case IDA_NO_QUADSENS:
	caml_raise_constant(IDAS_EXN(QuadSensNotInitialized));

    case IDA_QSRHS_FAIL:
	caml_raise_constant(IDAS_EXN(QuadSensRhsFuncFailure));

    case IDA_FIRST_QSRHS_ERR:
	caml_raise_constant(IDAS_EXN(FirstQuadSensRhsFuncFailure));

    case IDA_REP_QSRHS_ERR:
	caml_raise_constant(IDAS_EXN(RepeatedQuadSensRhsFuncFailure));

    case IDA_NO_ADJ:
	caml_raise_constant(IDAS_EXN(AdjointNotInitialized));

    case IDA_NO_FWD:
	caml_raise_constant(IDAS_EXN(NoForwardCall));

    case IDA_NO_BCK:
	caml_raise_constant(IDAS_EXN(NoBackwardProblem));

    case IDA_REIFWD_FAIL:
	caml_raise_constant(IDAS_EXN(ForwardReinitFailure));

    case IDA_BAD_TB0:
	caml_raise_constant(IDAS_EXN(BadFinalTime));

    case IDA_FWD_FAIL:
	caml_raise_constant(IDAS_EXN(ForwardFailure));

    case IDA_GETY_BADT:
	caml_raise_constant(IDAS_EXN(BadOutputTime));

    default:
	/* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 IDAGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

/* Callbacks */

value sunml_idas_make_jac_arg(sunrealtype t, N_Vector y, N_Vector yp,
			N_Vector yb, N_Vector ypb, N_Vector resb,
			sunrealtype coef, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_IDAS_ADJ_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_YP, NVEC_BACKLINK(yp));
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_YB, NVEC_BACKLINK(yb));
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_YPB, NVEC_BACKLINK(ypb));
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_RESB, NVEC_BACKLINK(resb));
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_COEF, caml_copy_double(coef));
    Store_field(r, RECORD_IDAS_ADJ_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

static int quadrhsfn(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rhsQ,
		     void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal3(session, cb, sensext);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(yp);
    args[3] = NVEC_BACKLINK(rhsQ);

    WEAK_DEREF (session, *(value*)user_data);
    sensext = IDA_SENSEXT_FROM_ML (session);

    // the data payloads inside args[2] and args[3] are only valid during
    // this call, afterward that memory goes back to ida. These bigarrays
    // must not be retained by closure_quadrhsfn! If it wants a permanent
    // copy, then it has to make it manually.

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (IDAS_QUADRHSFN_FROM_EXT (sensext), 4, args);

    CAMLreturnT (int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}


static int sensresfn(int Ns, sunrealtype t,
		     N_Vector y, N_Vector yp, N_Vector resval,
		     N_Vector *yS, N_Vector *ypS, N_Vector *resvalS,
		     void *user_data,
		     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal4(args, session, sensext, vt);
    value *backref = user_data;

    WEAK_DEREF (session, *backref);
    sensext = IDA_SENSEXT_FROM_ML(session);

    args = caml_alloc_tuple (RECORD_IDAS_SENSRESFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_T,
		 caml_copy_double (t));
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_Y,
		 NVEC_BACKLINK(y));
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_YP,
		 NVEC_BACKLINK(yp));
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_RES,
		 NVEC_BACKLINK(resval));
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_SENS,
	         IDAS_SENSARRAY1_FROM_EXT(sensext));
    sunml_nvectors_into_array (Ns, IDAS_SENSARRAY1_FROM_EXT(sensext), yS);
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_SENSP,
	         IDAS_SENSARRAY2_FROM_EXT(sensext));
    sunml_nvectors_into_array (Ns, IDAS_SENSARRAY2_FROM_EXT(sensext), ypS);
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_TMP,
		 sunml_ida_make_triple_tmp (tmp1, tmp2, tmp3));

    sunml_nvectors_into_array (Ns, IDAS_SENSARRAY3_FROM_EXT(sensext), resvalS);

    value r = caml_callback2_exn (IDAS_SENSRESFN_FROM_EXT(sensext), args,
				  IDAS_SENSARRAY3_FROM_EXT(sensext));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int quadsensrhsfn(int ns, sunrealtype t, N_Vector yy, N_Vector yp,
			 N_Vector *yyS, N_Vector *ypS,
			 N_Vector rrQ, N_Vector *rhsvalQS,
			 void *user_data,
		         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal3(args, session, sensext);

    WEAK_DEREF (session, *(value*)user_data);
    sensext = IDA_SENSEXT_FROM_ML(session);

    args = caml_alloc_tuple (RECORD_IDAS_QUADSENSRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_T,
		 caml_copy_double (t));
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_Y,
		 NVEC_BACKLINK (yy));
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_YP,
		 NVEC_BACKLINK (yp));
    sunml_nvectors_into_array (ns, IDAS_SENSARRAY1_FROM_EXT(sensext), yyS);
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_SENS,
		 IDAS_SENSARRAY1_FROM_EXT(sensext));
    sunml_nvectors_into_array (ns, IDAS_SENSARRAY2_FROM_EXT(sensext), ypS);
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_SENSP,
		 IDAS_SENSARRAY2_FROM_EXT(sensext));
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_TMP,
		 sunml_ida_make_triple_tmp (tmp1, tmp2, tmp3));

    sunml_nvectors_into_array (ns, IDAS_SENSARRAY3_FROM_EXT(sensext), rhsvalQS);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(IDAS_QUADSENSRHSFN_FROM_EXT(sensext),
				 args, IDAS_SENSARRAY3_FROM_EXT(sensext));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}


static int bresfn(sunrealtype t, N_Vector y, N_Vector yp,
		  N_Vector yB, N_Vector ypB,
		  N_Vector resvalB, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, bsensext);

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    args = caml_alloc_tuple (RECORD_IDAS_ADJ_BRESFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YP, NVEC_BACKLINK (yp));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YB, NVEC_BACKLINK (yB));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YBP, NVEC_BACKLINK (ypB));

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (IDAS_BRESFN_FROM_EXT(bsensext),
				  args, NVEC_BACKLINK (resvalB));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bresfn_sens(sunrealtype t, N_Vector y, N_Vector yp,
		       N_Vector *yS, N_Vector *ypS,
		       N_Vector yB, N_Vector ypB,
		       N_Vector resvalB, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, bsensext);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));

    args = caml_alloc_tuple (RECORD_IDAS_ADJ_BRESFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YP, NVEC_BACKLINK (yp));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YB, NVEC_BACKLINK (yB));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YBP, NVEC_BACKLINK (ypB));

    sunml_nvectors_into_array (ns, IDAS_BSENSARRAY1_FROM_EXT (bsensext), yS);
    sunml_nvectors_into_array (ns, IDAS_BSENSARRAY2_FROM_EXT (bsensext), ypS);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (IDAS_BRESFN_SENS_FROM_EXT(bsensext), args);
    if (! Is_exception_result (r))
	r = caml_callback3_exn (r, IDAS_BSENSARRAY1_FROM_EXT (bsensext),
				IDAS_BSENSARRAY2_FROM_EXT (bsensext),
				NVEC_BACKLINK (resvalB));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bprecsetupfn(sunrealtype t,
			N_Vector yy,
			N_Vector yp,
			N_Vector yB,
			N_Vector ypB,
			N_Vector resvalB,
			sunrealtype cjB,
			void *user_data
#if SUNDIALS_LIB_VERSION < 300
			,
			N_Vector tmp1B,
			N_Vector tmp2B,
			N_Vector tmp3B
#endif
			)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = sunml_idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB, Val_unit);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDAS_BSPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, arg);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bprecsetupfn_sens(sunrealtype t,
			     N_Vector yy,
			     N_Vector yp,
			     N_Vector *yyS,
			     N_Vector *ypS,
			     N_Vector yB,
			     N_Vector ypB,
			     N_Vector resvalB,
			     sunrealtype cjB,
			     void *user_data
#if SUNDIALS_LIB_VERSION < 300
			     ,
			     N_Vector tmp1B,
			     N_Vector tmp2B,
			     N_Vector tmp3B
#endif
			     )
{
    CAMLparam0();
    CAMLlocal3(session, cb, bsensext);
    CAMLlocalN(args, 3);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB, Val_unit);
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    sunml_nvectors_into_array (ns, args[1], yyS);
    sunml_nvectors_into_array (ns, args[2], ypS);

    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDAS_BSPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int bprecsolvefn(sunrealtype t,
		        N_Vector yy,
			N_Vector yp,
			N_Vector yB,
			N_Vector ypB,
			N_Vector resvalB,
			N_Vector rvecB,
			N_Vector zvecB,
			sunrealtype cjB,
			sunrealtype deltaB,
			void *user_data
#if SUNDIALS_LIB_VERSION < 300
			,
			N_Vector tmpB
#endif
			)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal2(session, cb);

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB, Val_unit);
    args[1] = NVEC_BACKLINK (rvecB);
    args[2] = NVEC_BACKLINK (zvecB);
    args[3] = caml_copy_double (deltaB);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDAS_BSPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bprecsolvefn_sens(sunrealtype t,
			     N_Vector yy,
			     N_Vector yp,
			     N_Vector *yyS,
			     N_Vector *ypS,
			     N_Vector yB,
			     N_Vector ypB,
			     N_Vector resvalB,
			     N_Vector rvecB,
			     N_Vector zvecB,
			     sunrealtype cjB,
			     sunrealtype deltaB,
			     void *user_data
#if SUNDIALS_LIB_VERSION < 300
			     ,
			     N_Vector tmpB
#endif
			     )
{
    CAMLparam0();
    CAMLlocalN(args, 6);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB, Val_unit);

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    sunml_nvectors_into_array (ns, args[1], yyS);
    sunml_nvectors_into_array (ns, args[2], ypS);

    args[3] = NVEC_BACKLINK (rvecB);
    args[4] = NVEC_BACKLINK (zvecB);
    args[5] = caml_copy_double (deltaB);

    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDAS_BSPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 6, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

#if SUNDIALS_LIB_VERSION >= 300
static int bjacsetupfn(sunrealtype t,
		       N_Vector yy,
		       N_Vector yp,
		       N_Vector yyB,
		       N_Vector ypB,
		       N_Vector resvalB,
		       sunrealtype cjB,
		       void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB, Val_unit);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, arg);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

static int bjacsetupfn_sens(sunrealtype t,
			    N_Vector yy,
			    N_Vector yp,
			    N_Vector *yyS,
			    N_Vector *ypS,
			    N_Vector yyB,
			    N_Vector ypB,
			    N_Vector resvalB,
			    sunrealtype cjB,
			    void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB, Val_unit);

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    sunml_nvectors_into_array (ns, args[1], yyS);
    sunml_nvectors_into_array (ns, args[2], ypS);

    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}
#endif

static int bjactimesfn(sunrealtype t,
		       N_Vector yy,
		       N_Vector yp,
		       N_Vector yyB,
		       N_Vector ypB,
		       N_Vector resvalB,
		       N_Vector vB,
		       N_Vector JvB,
		       sunrealtype cjB,
		       void *user_data,
		       N_Vector tmp1B,
		       N_Vector tmp2B)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal2(session, cb);

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_double_tmp (tmp1B, tmp2B));
    args[1] = NVEC_BACKLINK(vB);
    args[2] = NVEC_BACKLINK(JvB);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bjactimesfn_sens(sunrealtype t,
			    N_Vector yy,
			    N_Vector yp,
			    N_Vector *yyS,
			    N_Vector *ypS,
			    N_Vector yyB,
			    N_Vector ypB,
			    N_Vector resvalB,
			    N_Vector vB,
			    N_Vector JvB,
			    sunrealtype cjB,
			    void *user_data,
			    N_Vector tmp1B,
			    N_Vector tmp2B)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_double_tmp (tmp1B, tmp2B));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    sunml_nvectors_into_array (ns, args[1], yyS);
    sunml_nvectors_into_array (ns, args[2], ypS);

    args[3] = NVEC_BACKLINK(vB);
    args[4] = NVEC_BACKLINK(JvB);

    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 5, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}
#endif

#if 300 <= SUNDIALS_LIB_VERSION

static int bjacfn_nosens(sunrealtype t,
			 sunrealtype cjB,
			 N_Vector yy,
			 N_Vector yp,
			 N_Vector yyB,
			 N_Vector ypB,
			 N_Vector resvalB,
			 SUNMatrix JacB,
			 void *user_data,
			 N_Vector tmp1B,
			 N_Vector tmp2B,
			 N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));
    args[1] = MAT_BACKLINK(JacB);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bjacfn_withsens(sunrealtype t,
			   sunrealtype cjB,
			   N_Vector yy,
			   N_Vector yp,
			   N_Vector *yS,
			   N_Vector *ypS,
			   N_Vector yyB,
			   N_Vector ypB,
			   N_Vector resvalB,
			   SUNMatrix JacB,
			   void *user_data,
			   N_Vector tmp1B,
			   N_Vector tmp2B,
			   N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    sunml_nvectors_into_array (ns, args[1], yS);
    sunml_nvectors_into_array (ns, args[2], ypS);

    args[3] = MAT_BACKLINK(JacB);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#else
static int bjacfn_nosens(long int NeqB,
			 sunrealtype t,
			 sunrealtype cjB,
			 N_Vector yy,
			 N_Vector yp,
			 N_Vector yyB,
			 N_Vector ypB,
			 N_Vector resvalB,
			 DlsMat JacB,
			 void *user_data,
			 N_Vector tmp1B,
			 N_Vector tmp2B,
			 N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(JacB));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bjacfn_withsens(long int NeqB,
			   sunrealtype t,
			   sunrealtype cjB,
			   N_Vector yy,
			   N_Vector yp,
			   N_Vector *yS,
			   N_Vector *ypS,
			   N_Vector yyB,
			   N_Vector ypB,
			   N_Vector resvalB,
			   DlsMat JacB,
			   void *user_data,
			   N_Vector tmp1B,
			   N_Vector tmp2B,
			   N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal4(session, bsensext, cb, dmat);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(JacB));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    sunml_nvectors_into_array (ns, args[1], yS);
    sunml_nvectors_into_array (ns, args[2], ypS);

    args[3] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

static int bbandjacfn_nosens(long int NeqB,
			     long int mupperb,
			     long int mlowerb,
			     sunrealtype t,
			     sunrealtype cjB,
			     N_Vector yy,
			     N_Vector yp,
			     N_Vector yyB,
			     N_Vector ypB,
			     N_Vector resvalB,
			     DlsMat JacB,
			     void *user_data,
			     N_Vector tmp1B,
			     N_Vector tmp2B,
			     N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field(cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(JacB));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_triple_tmp(tmp1B, tmp2B, tmp3B));
    args[1] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bbandjacfn_withsens(long int NeqB,
			       long int mupperb,
			       long int mlowerb,
			       sunrealtype t,
			       sunrealtype cjB,
			       N_Vector yy,
			       N_Vector yp,
			       N_Vector *yS,
			       N_Vector *ypS,
			       N_Vector yyB,
			       N_Vector ypB,
			       N_Vector resvalB,
			       DlsMat JacB,
			       void *user_data,
			       N_Vector tmp1B,
			       N_Vector tmp2B,
			       N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal4(session, bsensext, cb, bmat);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field(cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(JacB));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        sunml_ida_make_triple_tmp(tmp1B, tmp2B, tmp3B));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    sunml_nvectors_into_array (ns, args[2], yS);
    sunml_nvectors_into_array (ns, args[3], ypS);

    args[3] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

#endif

static int bquadrhsfn(sunrealtype t, N_Vector y, N_Vector yp,
		      N_Vector yB, N_Vector ypB,
		      N_Vector rhsvalBQ, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, sensext);

    args = caml_alloc_tuple (RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_Y, NVEC_BACKLINK(y));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_YP, NVEC_BACKLINK(yp));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_YB, NVEC_BACKLINK(yB));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_YBP, NVEC_BACKLINK(ypB));

    WEAK_DEREF (session, *(value*)user_data);
    sensext = IDA_SENSEXT_FROM_ML (session);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (IDAS_BQUADRHSFN_FROM_EXT (sensext),
				  args, NVEC_BACKLINK (rhsvalBQ));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bquadrhsfn_sens(sunrealtype t, N_Vector y, N_Vector yp,
			   N_Vector *yS, N_Vector *ypS,
			   N_Vector yB, N_Vector ypB,
			   N_Vector rhsvalBQS, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, sensext);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    sensext = IDA_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(sensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));

    args = caml_alloc_tuple (RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_YP, NVEC_BACKLINK (yp));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_YB, NVEC_BACKLINK (yB));
    Store_field (args, RECORD_IDAS_ADJ_BQUADRHSFN_ARGS_YBP,
		 NVEC_BACKLINK (ypB));

    sunml_nvectors_into_array (ns, IDAS_BSENSARRAY1_FROM_EXT (sensext), yS);
    sunml_nvectors_into_array (ns, IDAS_BSENSARRAY2_FROM_EXT (sensext), ypS);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (IDAS_BQUADRHSFN_SENS_FROM_EXT(sensext), args);
    if (!Is_exception_result (r))
	r = caml_callback3_exn (r, IDAS_BSENSARRAY1_FROM_EXT (sensext),
				IDAS_BSENSARRAY2_FROM_EXT (sensext),
				NVEC_BACKLINK (rhsvalBQS));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#if 530 <= SUNDIALS_LIB_VERSION
static int bjactimesresfn (sunrealtype t, N_Vector y, N_Vector yp,
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
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    CAMLreturnT (int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

/* quadrature interface */

CAMLprim value sunml_idas_quad_init (value vsession, value vyQ0)
{

    CAMLparam2 (vsession, vyQ0);
    N_Vector yQ0 = NVEC_VAL (vyQ0);
    int flag;

    flag = IDAQuadInit (IDA_MEM_FROM_ML (vsession), quadrhsfn, yQ0);
    SCHECK_FLAG ("IDAQuadInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quad_reinit (value vsession, value vyQ0)
{

    CAMLparam2 (vsession, vyQ0);
    N_Vector yQ0 = NVEC_VAL (vyQ0);
    int flag;

    flag = IDAQuadReInit (IDA_MEM_FROM_ML (vsession), yQ0);
    SCHECK_FLAG ("IDAQuadReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quad_sv_tolerances(value vdata, value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = IDAQuadSVtolerances(IDA_MEM_FROM_ML(vdata),
				   Double_val(reltol), atol_nv);
    SCHECK_FLAG("IDAQuadSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quad_get(value vdata, value vyq)
{
    CAMLparam2(vdata, vyq);
    N_Vector yq = NVEC_VAL(vyq);
    sunrealtype tret;

    int flag = IDAGetQuad(IDA_MEM_FROM_ML(vdata), &tret, yq);
    SCHECK_FLAG("IDAGetQuad", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_idas_quad_get_dky(value vdata, value vt, value vk,
				   value vdkyq)
{
    CAMLparam4(vdata, vt, vk, vdkyq);
    N_Vector dkyq = NVEC_VAL(vdkyq);

    int flag = IDAGetQuadDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkyq);
	    
    SCHECK_FLAG("IDAGetQuadDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quad_get_err_weights(value vdata, value veqweight)
{
    CAMLparam2(vdata, veqweight);
    N_Vector eqweight = NVEC_VAL(veqweight);

    int flag = IDAGetQuadErrWeights(IDA_MEM_FROM_ML(vdata), eqweight);
    SCHECK_FLAG("IDAGetQuadErrWeights", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value sunml_idas_quad_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = IDASetQuadErrCon(IDA_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("IDASetQuadErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quad_ss_tolerances(value vdata,
					 value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDAQuadSStolerances(IDA_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    SCHECK_FLAG("IDAQuadSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quad_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadNumRhsEvals(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_quad_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadNumErrTestFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_quad_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfqevals;
    long int nqetfails;

    flag = IDAGetQuadStats(IDA_MEM_FROM_ML(vdata), &nfqevals,
	    &nqetfails);
    SCHECK_FLAG("IDAGetQuadStats", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(nfqevals));
    Store_field(r, 1, Val_long(nqetfails));

    CAMLreturn(r);
}

/* sensitivity interface */

#if 400 <= SUNDIALS_LIB_VERSION
#if SUNDIALS_LIB_VERSION < 630
// hack to work around lack of CVodeGetUserData
typedef struct {
#if 600 <= SUNDIALS_LIB_VERSION
  SUNContext ida_sunctx;
#endif
  sunrealtype ida_uround;
  IDAResFn ida_res;
  void     *ida_user_data;
  //...
} *StartOf_IDAMem;
#endif

static value sunml_idas_session_to_value(void *ida_mem)
{
    value session;
#if 630 <= SUNDIALS_LIB_VERSION
    void *user_data = NULL;
    IDAGetUserData(ida_mem, &user_data);
#else
    void *user_data = ((StartOf_IDAMem)ida_mem)->ida_user_data;
#endif

    WEAK_DEREF (session, *(value*)user_data);
    return session;
}

static value sunml_idas_bsession_to_value(void *ida_mem)
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

#if 630 <= SUNDIALS_LIB_VERSION
    void *parent = NULL;
    void *parent_mem = NULL;
    IDAGetUserData(ida_mem, &parent);
    IDAGetUserData(parent, &parent_mem);
#else
    void *parent = ((StartOf_IDAMem)ida_mem)->ida_user_data;
    void *parent_mem = ((StartOf_IDAMem)parent)->ida_user_data;
#endif
    WEAK_DEREF (vparent, *(value*)parent_mem);
    vchildptr = sunml_wrap_session_pointer(ida_mem);

    vsession = caml_callback2(
		*caml_named_value("Idas.revlookup_bsession"),
		vparent, vchildptr);

    assert (vsession != Val_none);
    vsession = Some_val(vsession);

    CAMLreturn (vsession);
}

static void* sunml_idas_session_from_value(value vida_mem)
{
    return (IDA_MEM_FROM_ML(vida_mem));
}
#endif


CAMLprim value sunml_idas_set_nonlinear_solver_sim(value vida_mem,
						   value vnlsolv)
{
    CAMLparam2(vida_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_idas_session_to_value,
				   sunml_idas_session_from_value);

    int flag = IDASetNonlinearSolverSensSim(ida_mem, nlsolv);
    SCHECK_FLAG ("IDASetNonlinearSolverSensSim", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_set_nonlinear_solver_stg(value vida_mem,
						   value vnlsolv)
{
    CAMLparam2(vida_mem, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_idas_session_to_value,
				   sunml_idas_session_from_value);

    int flag = IDASetNonlinearSolverSensStg(ida_mem, nlsolv);
    SCHECK_FLAG ("IDASetNonlinearSolverSensStg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

static int decode_sens_method(value vmethod)
{
    switch (Int_val(vmethod)) {
    case VARIANT_IDAS_SENS_METHOD_SIMULTANEOUS:
	return IDA_SIMULTANEOUS;

    case VARIANT_IDAS_SENS_METHOD_STAGGERED:
	return IDA_STAGGERED;

    default:
	caml_failwith("Illegal sens method.");
    }
}


CAMLprim value sunml_idas_sens_init(value vdata, value vmethod, value vrhsfn,
				value vyS0, value vypS0)
{
    CAMLparam5(vdata, vmethod, vrhsfn, vyS0, vypS0);
    int ns = Wosize_val (vyS0);	/* vyS0 : nvector array */
    N_Vector *yS0 = sunml_nvector_array_alloc(vyS0);
    N_Vector *ypS0 = sunml_nvector_array_alloc(vypS0);

    int flag = IDASensInit(IDA_MEM_FROM_ML(vdata), ns,
			   decode_sens_method(vmethod),
			   ((Bool_val(vrhsfn)) ? sensresfn : NULL),
			   yS0, ypS0);
    sunml_nvector_array_free(yS0); 
    sunml_nvector_array_free(ypS0); 
    SCHECK_FLAG("IDASensInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_reinit(value vdata, value vmethod, value vyS0,
				  value vypS0)
{
    CAMLparam4(vdata, vmethod, vyS0, vypS0);
    CAMLlocal1(r);
    int flag;
    N_Vector *yS0 = sunml_nvector_array_alloc(vyS0);
    N_Vector *ypS0 = sunml_nvector_array_alloc(vypS0);

    flag = IDASensReInit(IDA_MEM_FROM_ML(vdata),
			 decode_sens_method(vmethod),
			 yS0, ypS0);
    sunml_nvector_array_free(yS0);
    sunml_nvector_array_free(ypS0);
    SCHECK_FLAG("IDASensReInit", flag);

    CAMLreturn (Val_unit);
}

static void sens_calc_ic (void *ida_mem, value session, int icopt, sunrealtype tout1,
			  value vy, value vyp, value vys, value vyps)
{
    CAMLparam5 (session, vy, vyp, vys, vyps);
    CAMLlocal1 (exn);
    int flag;
    N_Vector y, yp;
    N_Vector *ys;
    N_Vector *yps;

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
    y   = (Is_block (vy))  ? NVEC_VAL (Field (vy, 0))  : NULL;
    yp  = (Is_block (vyp)) ? NVEC_VAL (Field (vyp, 0)) : NULL;
    if (y != NULL || yp != NULL) {
	flag = IDAGetConsistentIC (ida_mem, y, yp);
	CHECK_FLAG ("IDAGetConsistentIC", flag);
    }
    ys  = (Is_block (vys))  ? sunml_nvector_array_alloc (Field (vys, 0))  : NULL;
    yps = (Is_block (vyps)) ? sunml_nvector_array_alloc (Field (vyps, 0)) : NULL;
    if (ys != NULL || yps != NULL) {
	flag = IDAGetSensConsistentIC (ida_mem, ys, yps);
	CHECK_FLAG ("IDAGetConsistentIC", flag);
	if (ys)  sunml_nvector_array_free(ys); 
	if (yps) sunml_nvector_array_free(yps); 
    }
    CAMLreturn0;
}

CAMLprim value sunml_ida_sens_calc_ic_y(value vida_mem, value vy, value vys, value tout1)
{
    CAMLparam4 (vida_mem, vy, vys, tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    sens_calc_ic (ida_mem, vida_mem, IDA_Y_INIT, Double_val (tout1),
		  vy, Val_none, vys, Val_none);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_sens_calc_ic_ya_ydp(value vida_mem, value y, value yp,
					 value ys, value yps, value tout1)
{
    CAMLparam5 (vida_mem, y, yp, ys, yps);
    CAMLxparam1 (tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    sens_calc_ic (ida_mem, vida_mem, IDA_YA_YDP_INIT, Double_val (tout1),
		  y, yp, ys, yps);

    CAMLreturn (Val_unit);
}

BYTE_STUB6(sunml_ida_sens_calc_ic_ya_ydp)

CAMLprim value sunml_idas_sens_get(value vdata, value vys)
{
    CAMLparam2(vdata, vys);
    N_Vector *ys = sunml_nvector_array_alloc(vys);
    sunrealtype tret;

    int flag = IDAGetSens(IDA_MEM_FROM_ML(vdata), &tret, ys);
    sunml_nvector_array_free(ys);
    SCHECK_FLAG("IDAGetSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_idas_sens_get_dky(value vdata, value vt, value vk, value vdkys)
{
    CAMLparam4(vdata, vt, vk, vdkys);
    N_Vector *dkys = sunml_nvector_array_alloc(vdkys);

    int flag = IDAGetSensDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkys);
    sunml_nvector_array_free(dkys);
    SCHECK_FLAG("IDAGetSensDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_get1(value vdata, value vis, value vys)
{
    CAMLparam3(vdata, vis, vys);
    N_Vector ys = NVEC_VAL(vys);
    sunrealtype tret;

    int flag = IDAGetSens1(IDA_MEM_FROM_ML(vdata), &tret, Int_val(vis), ys);
    SCHECK_FLAG("IDAGetSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_idas_sens_get_dky1(value vdata, value vt, value vk,
				    value vis, value vdkys)
{
    CAMLparam5(vdata, vt, vk, vis, vdkys);
    N_Vector dkys = NVEC_VAL(vdkys);

    int flag = IDAGetSensDky1(IDA_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), Int_val(vis), dkys);
	    
    SCHECK_FLAG("IDAGetSensDky1", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_get_err_weights(value vdata, value vesweight)
{
    CAMLparam2(vdata, vesweight);
    N_Vector *esweight = sunml_nvector_array_alloc(vesweight);

    int flag = IDAGetSensErrWeights(IDA_MEM_FROM_ML(vdata), esweight);
    sunml_nvector_array_free(esweight);
    SCHECK_FLAG("IDAGetSensErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_set_err_con(value vdata, value verrcons)
{
    CAMLparam2(vdata, verrcons);
    int flag;
    
    flag = IDASetSensErrCon(IDA_MEM_FROM_ML(vdata), Bool_val(verrcons));
    SCHECK_FLAG("IDASetSensErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_ss_tolerances(value vdata,
					 value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDASensSStolerances(IDA_MEM_FROM_ML(vdata),
		 Double_val(reltol), REAL_ARRAY(abstol));
    SCHECK_FLAG("IDASensSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_sv_tolerances(value vdata, value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = sunml_nvector_array_alloc(abstol);

    int flag = IDASensSVtolerances(IDA_MEM_FROM_ML(vdata),
				   Double_val(reltol), atol_nv);
    sunml_nvector_array_free(atol_nv);
    SCHECK_FLAG("IDASensSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASensEEtolerances(IDA_MEM_FROM_ML(vdata));
    SCHECK_FLAG("IDASensEEtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_set_params(value vdata, value vparams)
{
    CAMLparam2(vdata, vparams);
    CAMLlocal3(vp, vpbar, vplist);

    sunrealtype *p = NULL;
    sunrealtype *pbar = NULL;
    int *plist = NULL;
    int i, ns;

    vp = Field(vparams, RECORD_IDAS_SENS_PARAMS_PVALS);
    vpbar = Field(vparams, RECORD_IDAS_SENS_PARAMS_PBAR);
    vplist = Field(vparams, RECORD_IDAS_SENS_PARAMS_PLIST);

    if (vp != Val_none) p = REAL_ARRAY(Some_val(vp));
    if (vpbar != Val_none) pbar = REAL_ARRAY(Some_val(vpbar));

    if (vplist != Val_none) {
	vplist = Some_val(vplist);
	ns = Wosize_val (vplist); /* vplist : int array */
	plist = calloc(ns, sizeof(int));

	for (i=0; i < ns; ++i) {
	    plist[i] = Int_val(Field(vplist, i));
	}
    }

    int flag = IDASetSensParams(IDA_MEM_FROM_ML(vdata), p, pbar, plist);
    if (plist != NULL) free(plist);
    SCHECK_FLAG("IDASetSensParams", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_toggle_off(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASensToggleOff(IDA_MEM_FROM_ML(vdata));
    SCHECK_FLAG("IDASensToggleOff", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_set_dq_method(value vdata, value vdqtype,
					 value vdqrhomax)
{
    CAMLparam3(vdata, vdqtype, vdqrhomax);
    int dqtype;

    switch (Int_val(vdqtype)) {
    case VARIANT_IDAS_SENS_DQ_METHOD_CENTERED:
	dqtype = IDA_CENTERED;
	break;

    case VARIANT_IDAS_SENS_DQ_METHOD_FORWARD:
	dqtype = IDA_FORWARD;
	break;

    default:
	caml_failwith("Illegal dq method.");
    }

    int flag = IDASetSensDQMethod(IDA_MEM_FROM_ML(vdata), dqtype,
				    Double_val(vdqrhomax));
    SCHECK_FLAG("IDASetSensMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_set_max_nonlin_iters(value vdata, value vmaxcors)
{
    CAMLparam2(vdata, vmaxcors);

    int flag = IDASetSensMaxNonlinIters(IDA_MEM_FROM_ML(vdata),
	    Int_val(vmaxcors));
    SCHECK_FLAG("IDASetSensMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_get_num_res_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumResEvals(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumResEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_sens_get_num_res_evals_sens(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetNumResEvalsSens(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetNumResEvalsSens", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_sens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumErrTestFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_sens_get_num_step_solve_fails(value vdata)
{
    CAMLparam1(vdata);
    long int v;

#if 620 <= SUNDIALS_LIB_VERSION
    int flag;
    flag = IDAGetNumStepSensSolveFails(IDA_MEM_FROM_ML(vdata), &v);
    CHECK_FLAG("IDAGetNumStepSensSolveFails", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_sens_get_num_lin_solv_setups(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumLinSolvSetups(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_sens_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfsevals;
    long int nfevalss;
    long int nsetfails;
    long int nlinsetupss;

    flag = IDAGetSensStats(IDA_MEM_FROM_ML(vdata), &nfsevals,
	    &nfevalss, &nsetfails, &nlinsetupss);
    SCHECK_FLAG("IDAGetSensStats", flag);

    r = caml_alloc_tuple(RECORD_IDAS_SENS_STATS_SIZE);
    Store_field(r, RECORD_IDAS_SENS_STATS_NUM_SENS_EVALS, Val_long(nfsevals));
    Store_field(r, RECORD_IDAS_SENS_STATS_NUM_RES_EVALS, Val_long(nfevalss));
    Store_field(r, RECORD_IDAS_SENS_STATS_NUM_ERR_TEST_FAILS,
							   Val_long(nsetfails));
    Store_field(r, RECORD_IDAS_SENS_STATS_NUM_LIN_SOLV_SETUPS,
							 Val_long(nlinsetupss));

    CAMLreturn(r);
}

CAMLprim value sunml_idas_sens_get_num_nonlin_solv_iters(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumNonlinSolvIters(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_sens_get_num_nonlin_solv_conv_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumNonlinSolvConvFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_sens_get_nonlin_solv_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nsniters;
    long int nsncfails;

    flag = IDAGetSensNonlinSolvStats(IDA_MEM_FROM_ML(vdata), &nsniters,
	    &nsncfails);
    SCHECK_FLAG("IDAGetSensNonlinSolvStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nsniters));
    Store_field(r, 1, Val_long(nsncfails));

    CAMLreturn(r);
}

CAMLprim value sunml_idas_sens_get_current_y_sens(value vdata, value vns)
{
    CAMLparam2(vdata, vns);
    CAMLlocal1(r);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector *yyS;

    int flag = IDAGetCurrentYSens(IDA_MEM_FROM_ML(vdata), &yyS);
    SCHECK_FLAG("IDAGetCurrentYSens", flag);

    r = sunml_wrap_to_nvector_table(Int_val(vns), yyS);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (r);
}

CAMLprim value sunml_idas_sens_get_current_yp_sens(value vdata, value vns)
{
    CAMLparam2(vdata, vns);
    CAMLlocal1(r);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector *ypS;

    int flag = IDAGetCurrentYpSens(IDA_MEM_FROM_ML(vdata), &ypS);
    SCHECK_FLAG("IDAGetCurrentYpSens", flag);

    r = sunml_wrap_to_nvector_table(Int_val(vns), ypS);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (r);
}

CAMLprim value sunml_idas_get_nonlin_system_data_sens(value vida_mem,
							value vns)
{
    CAMLparam2(vida_mem, vns);
    CAMLlocal1(vnv);
#if 540 <= SUNDIALS_LIB_VERSION
    int ns = Int_val(vns);
    sunrealtype tn, cj;
    N_Vector *yySpred, *ypSpred, *yySn, *ypSn;
    void *user_data;

    int flag = IDAGetNonlinearSystemDataSens(IDA_MEM_FROM_ML(vida_mem),
		    &tn, &yySpred, &ypSpred, &yySn, &ypSn, &cj, &user_data);
    CHECK_FLAG("IDAGetNonlinearSystemDataSens", flag);

    vnv = caml_alloc_tuple(RECORD_IDAS_NONLIN_SYSTEM_DATA_SIZE);
    Store_field(vnv, RECORD_IDAS_NONLIN_SYSTEM_DATA_TN,
	    caml_copy_double(tn));
    Store_field(vnv, RECORD_IDAS_NONLIN_SYSTEM_DATA_YYSPRED,
	    sunml_wrap_to_nvector_table(ns, yySpred));
    Store_field(vnv, RECORD_IDAS_NONLIN_SYSTEM_DATA_YPSPRED,
	    sunml_wrap_to_nvector_table(ns, ypSpred));
    Store_field(vnv, RECORD_IDAS_NONLIN_SYSTEM_DATA_YYSN,
	    sunml_wrap_to_nvector_table(ns, yySn));
    Store_field(vnv, RECORD_IDAS_NONLIN_SYSTEM_DATA_YPSN,
	    sunml_wrap_to_nvector_table(ns, ypSn));
    Store_field(vnv, RECORD_IDAS_NONLIN_SYSTEM_DATA_CJ,
	    caml_copy_double(cj));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(vnv);
}

CAMLprim value sunml_idas_sens_compute_y_sens (value vida_mem,
					       value vycors, value vys)
{
    CAMLparam3(vida_mem, vycors, vys);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector *ycors = sunml_nvector_array_alloc(vycors);
    N_Vector *ys = sunml_nvector_array_alloc(vys);

    int flag = IDAComputeYSens (IDA_MEM_FROM_ML(vida_mem), ycors, ys);
    sunml_nvector_array_free(ycors); 
    sunml_nvector_array_free(ys); 
    SCHECK_FLAG("IDAComputeYSens", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_sens_compute_yp_sens (value vida_mem,
						value vycors, value vyps)
{
    CAMLparam3(vida_mem, vycors, vyps);

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector *ycors = sunml_nvector_array_alloc(vycors);
    N_Vector *yps = sunml_nvector_array_alloc(vyps);

    int flag = IDAComputeYpSens (IDA_MEM_FROM_ML(vida_mem), ycors, yps);
    sunml_nvector_array_free(ycors); 
    sunml_nvector_array_free(yps); 
    SCHECK_FLAG("IDAComputeYpSens", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

/* sensitivity/quadrature interface */

CAMLprim value sunml_idas_quadsens_init(value vdata, value vrhsfn, value vyqs0)
{
    CAMLparam3(vdata, vrhsfn, vyqs0);
    N_Vector *yqs0 = sunml_nvector_array_alloc(vyqs0);

    int flag = IDAQuadSensInit(IDA_MEM_FROM_ML(vdata),
			       Bool_val (vrhsfn) ? quadsensrhsfn : NULL,
			       yqs0);
    sunml_nvector_array_free(yqs0); 
    SCHECK_FLAG("IDAQuadSensInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_reinit(value vdata, value vyqs0)
{
    CAMLparam2(vdata, vyqs0);
    N_Vector *yqs0 = sunml_nvector_array_alloc(vyqs0);

    int flag = IDAQuadSensReInit(IDA_MEM_FROM_ML(vdata), yqs0);
    sunml_nvector_array_free(yqs0); 
    SCHECK_FLAG("IDAQuadSensReInit", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value sunml_idas_quadsens_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = IDASetQuadSensErrCon(IDA_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("IDASetQuadSensErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_ss_tolerances(value vdata,
					     value reltol,
					     value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDAQuadSensSStolerances(IDA_MEM_FROM_ML(vdata),
					 Double_val(reltol),
					 REAL_ARRAY(abstol));
    SCHECK_FLAG("IDAQuadSensSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_sv_tolerances(value vdata, value reltol,
					     value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = sunml_nvector_array_alloc(abstol);

    int flag = IDAQuadSensSVtolerances(IDA_MEM_FROM_ML(vdata),
				       Double_val(reltol), atol_nv);
    sunml_nvector_array_free(atol_nv);
    SCHECK_FLAG("IDAQuadSensSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDAQuadSensEEtolerances(IDA_MEM_FROM_ML(vdata));
    SCHECK_FLAG("IDAQuadSensEEtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_get(value vdata, value vyqs)
{
    CAMLparam2(vdata, vyqs);
    N_Vector *yqs = sunml_nvector_array_alloc(vyqs);
    sunrealtype tret;

    int flag = IDAGetQuadSens(IDA_MEM_FROM_ML(vdata), &tret, yqs);
    sunml_nvector_array_free(yqs);
    SCHECK_FLAG("IDAGetQuadSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_idas_quadsens_get1(value vdata, value vis, value vyqs)
{
    CAMLparam3(vdata, vis, vyqs);
    N_Vector yqs = NVEC_VAL(vyqs);
    sunrealtype tret;

    int flag = IDAGetQuadSens1(IDA_MEM_FROM_ML(vdata), &tret,
			       Int_val(vis), yqs);
    SCHECK_FLAG("IDAGetQuadSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_idas_quadsens_get_dky(value vdata, value vt, value vk,
				       value vdkyqs)
{
    CAMLparam4(vdata, vt, vk, vdkyqs);
    N_Vector *dkyqs = sunml_nvector_array_alloc(vdkyqs);

    int flag = IDAGetQuadSensDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
				 Int_val(vk), dkyqs);
    sunml_nvector_array_free(dkyqs);
    SCHECK_FLAG("IDAGetQuadSensDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_get_dky1(value vdata, value vt, value vk,
					value vis, value vdkyqs)
{
    CAMLparam5(vdata, vt, vk, vis, vdkyqs);
    N_Vector dkyqs = NVEC_VAL(vdkyqs);

    int flag = IDAGetQuadSensDky1(IDA_MEM_FROM_ML(vdata), Double_val(vt),
				    Int_val(vk), Int_val(vis), dkyqs);
    SCHECK_FLAG("IDAGetQuadSensDky1", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_get_err_weights(value vdata, value veqweights)
{
    CAMLparam2(vdata, veqweights);
    N_Vector *eqweights = sunml_nvector_array_alloc(veqweights);

    int flag = IDAGetQuadSensErrWeights(IDA_MEM_FROM_ML(vdata), eqweights);
    sunml_nvector_array_free(eqweights);
    SCHECK_FLAG("IDAGetQuadSensErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_quadsens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadSensNumRhsEvals(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_quadsens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadSensNumErrTestFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_idas_quadsens_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfqsevals;
    long int nqsetfails;

    flag = IDAGetQuadSensStats(IDA_MEM_FROM_ML(vdata), &nfqsevals,
				 &nqsetfails);
    SCHECK_FLAG("IDAGetQuadSensStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nfqsevals));
    Store_field(r, 1, Val_long(nqsetfails));

    CAMLreturn(r);
}


/* adjoint interface */

CAMLprim value sunml_idas_adj_set_nonlinear_solver(value vparent,
						   value vwhich, value vnlsolv)
{
    CAMLparam3(vparent, vwhich, vnlsolv);
#if 400 <= SUNDIALS_LIB_VERSION
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    SUNNonlinearSolver nlsolv = NLSOLVER_VAL(vnlsolv);
    int flag;

    sunml_nlsolver_set_to_from_mem(nlsolv,
				   sunml_idas_bsession_to_value,
				   sunml_idas_session_from_value);

    flag = IDASetNonlinearSolverB(ida_mem, Int_val(vwhich), nlsolv);
    CHECK_FLAG ("IDASetNonlinearSolverB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}


CAMLprim value sunml_idas_adj_init(value vdata, value vnd, value vinterptype)
{
    CAMLparam3(vdata, vnd, vinterptype);
    int interptype;

    switch(Int_val(vinterptype)) {
    case VARIANT_IDAS_ADJ_INTERPOLATION_POLYNOMIAL:
	interptype = IDA_POLYNOMIAL;
	break;

    case VARIANT_IDAS_ADJ_INTERPOLATION_HERMITE:
	interptype = IDA_HERMITE;
	break;

    default:
	caml_failwith("Illegal interpolation value.");
    }

    int flag = IDAAdjInit(IDA_MEM_FROM_ML(vdata), Long_val(vnd),
			  interptype);
    SCHECK_FLAG("IDAAdjInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_ss_tolerances(value vparent, value vwhich,
					value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);

    int flag = IDASStolerancesB(IDA_MEM_FROM_ML(vparent),
				Int_val(vwhich),
				Double_val(vreltol),
				Double_val(vabstol));
    SCHECK_FLAG("IDASStolerancesB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_spgmr (value vparent, value vwhich,
				       value vmaxl)
{
    CAMLparam3 (vparent, vwhich, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = IDASpgmrB (ida_mem, which, Int_val (vmaxl));
    SCHECK_FLAG ("IDASpgmrB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_spbcgs (value vparent, value vwhich,
				        value vmaxl)
{
    CAMLparam3 (vparent, vwhich, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = IDASpbcgB (ida_mem, which, Int_val (vmaxl));
    SCHECK_FLAG ("IDASpbcgB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_sptfqmr (value vparent, value vwhich,
					 value vmaxl)
{
    CAMLparam3 (vparent, vwhich, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = IDASptfqmrB (ida_mem, which, Int_val (vmaxl));
    SCHECK_FLAG ("IDASptfqmrB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_backward_normal(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = IDASolveB(IDA_MEM_FROM_ML(vdata), Double_val(vtbout),
			 IDA_NORMAL);
    SCHECK_FLAG("IDASolveB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_backward_one_step(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = IDASolveB(IDA_MEM_FROM_ML(vdata), Double_val(vtbout),
			 IDA_ONE_STEP);
    SCHECK_FLAG("IDASolveB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_max_ord(value vparent, value vwhich,
				      value vmaxord)
{
    CAMLparam3(vparent, vwhich, vmaxord);

    int flag = IDASetMaxOrdB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(vmaxord));
    SCHECK_FLAG("IDASetMaxOrdB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_max_num_steps(value vparent, value vwhich,
					    value vmxsteps)
{
    CAMLparam3(vparent, vwhich, vmxsteps);

    int flag = IDASetMaxNumStepsB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				    Long_val(vmxsteps));
    SCHECK_FLAG("IDASetMaxNumStepsB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_init_step(value vparent, value vwhich, value vhin)
{
    CAMLparam3(vparent, vwhich, vhin);

    int flag = IDASetInitStepB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Double_val(vhin));
    SCHECK_FLAG("IDASetInitStepB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_max_step(value vparent, value vwhich, value vhmax)
{
    CAMLparam3(vparent, vwhich, vhmax);

    int flag = IDASetMaxStepB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			      Double_val(vhmax));
    SCHECK_FLAG("IDASetMaxStepB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_constraints(value vparent, value vwhich,
					      value vconstraints)
{
    CAMLparam3(vparent, vwhich, vconstraints);
    N_Vector constraints = NVEC_VAL (vconstraints);

    int flag = IDASetConstraintsB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			          constraints);
    SCHECK_FLAG("IDASetConstraintsB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_clear_constraints(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);

    int flag = IDASetConstraintsB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			          NULL);
    SCHECK_FLAG("IDASetConstraintsB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_set_gs_type(value vparent, value vwhich,
					    value vgstype)
{
    CAMLparam3(vparent, vwhich, vgstype);
#if SUNDIALS_LIB_VERSION < 300
    int flag = IDASpilsSetGSTypeB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				  sunml_lsolver_gs_type(vgstype));
    SCHECK_FLAG("IDASpilsSetGSTypeB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_set_max_restarts(value vparent, value vwhich,
					         value vmaxr)
{
    CAMLparam3(vparent, vwhich, vmaxr);
#if SUNDIALS_LIB_VERSION < 300
    int flag = IDASpilsSetMaxRestartsB(IDA_MEM_FROM_ML(vparent),
				       Int_val(vwhich), Int_val(vmaxr));
    SCHECK_FLAG("IDASpilsSetMaxRestartsB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_eps_lin(value vparent, value vwhich,
					  value eplifac)
{
    CAMLparam3(vparent, vwhich, eplifac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDASetEpsLinB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			     Double_val(eplifac));
    SCHECK_FLAG("IDASetEpsLinB", flag);
#else
    int flag = IDASpilsSetEpsLinB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				  Double_val(eplifac));
    SCHECK_FLAG("IDASpilsSetEpsLinB", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_ls_norm_factor(value vparent, value vwhich,
					         value vfac)
{
    CAMLparam3(vparent, vwhich, vfac);

#if 540 <= SUNDIALS_LIB_VERSION
    int flag = IDASetLSNormFactorB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			           Double_val(vfac));
    SCHECK_FLAG("IDASetLSNormFactorB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_linear_solution_scaling(value vparent,
							  value vwhich,
					                  value vonoff)
{
    CAMLparam3(vparent, vwhich, vonoff);

#if 520 <= SUNDIALS_LIB_VERSION
    int flag = IDASetLinearSolutionScalingB(IDA_MEM_FROM_ML(vparent),
					    Int_val(vwhich),
					    Bool_val(vonoff));
    SCHECK_FLAG("IDASetLinearSolutionScalingB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_increment_factor(value vparent, value vwhich,
						   value dqincfac)
{
    CAMLparam3(vparent, vwhich, dqincfac);

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = IDASetIncrementFactorB(IDA_MEM_FROM_ML(vparent),
	    Int_val(vwhich), Double_val(dqincfac));
    SCHECK_FLAG("IDASetIncrementFactorB", flag);
#else
    int flag = IDASpilsSetIncrementFactorB(IDA_MEM_FROM_ML(vparent),
	    Int_val(vwhich), Double_val(dqincfac));
    SCHECK_FLAG("IDASpilsSetIncrementFactorB", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_set_maxl(value vparent, value vwhich,
					 value maxl)
{
    CAMLparam3(vparent, vwhich, maxl);
#if SUNDIALS_LIB_VERSION < 300
    int flag = IDASpilsSetMaxlB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(maxl));
    SCHECK_FLAG("IDASpilsSetMaxlB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_bsession_finalize(value vdata)
{
    if (IDA_MEM_FROM_ML(vdata) != NULL) {
	value *backref = IDA_BACKREF_FROM_ML(vdata);
	// NB: IDAFree() is *not* called: parents free-up backward problems
	sunml_sundials_free_value(backref);
    }
    return Val_unit;
}

static value forward_solve(value vdata, value vtout, value vy,
			   value vyp, int onestep)
{
    CAMLparam4(vdata, vtout, vy, vyp);
    CAMLlocal1(ret);
    N_Vector y  = NVEC_VAL(vy);
    N_Vector yp = NVEC_VAL(vyp);
    sunrealtype tret;
    int ncheck;
    enum ida_solver_result_tag solver_result = -1;

    int flag = IDASolveF(IDA_MEM_FROM_ML(vdata), Double_val(vtout), &tret,
			 y, yp, onestep ? IDA_ONE_STEP : IDA_NORMAL, &ncheck);
    switch (flag) {
    case IDA_SUCCESS:
	solver_result = VARIANT_IDA_SOLVER_RESULT_SUCCESS;
	break;

    case IDA_TSTOP_RETURN:
	solver_result = VARIANT_IDA_SOLVER_RESULT_STOPTIMEREACHED;
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
	sunml_idas_check_flag("IDASolveF", flag, IDA_MEM_FROM_ML(vdata));
    }

    assert (Field (vdata, RECORD_IDA_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple(3);
    Store_field(ret, 0, caml_copy_double(tret));
    Store_field(ret, 1, Val_int(ncheck));
    Store_field(ret, 2, Val_int(solver_result));

    CAMLreturn(ret);
}

CAMLprim value sunml_idas_adj_forward_normal(value vdata, value vtout,
					 value vy, value vyp)
{
    CAMLparam4(vdata, vtout, vy, vyp);
    CAMLreturn(forward_solve(vdata, vtout, vy, vyp, 0));
}

CAMLprim value sunml_idas_adj_forward_one_step(value vdata, value vtout,
					   value vy, value vyp)
{
    CAMLparam4(vdata, vtout, vy, vyp);
    CAMLreturn(forward_solve(vdata, vtout, vy, vyp, 1));
}

CAMLprim value sunml_idas_adj_sv_tolerances(value vparent, value vwhich,
					value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = IDASVtolerancesB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				Double_val(vreltol), atol_nv);
    SCHECK_FLAG("IDASStolerancesB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_set_preconditioner(value vparent,
						   value vwhich,
						   value vset_precsetup,
						   value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_precsetup, vusesens);
    void *mem = IDA_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    int flag;

#if 400 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vusesens)) {
	flag = IDASetPreconditionerBS(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn_sens : NULL,
		    bprecsolvefn_sens);
	SCHECK_FLAG ("IDASetPreconditionerBS", flag);
    } else {
	flag = IDASetPreconditionerB(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn : NULL,
		    bprecsolvefn);
	SCHECK_FLAG ("IDASetPreconditionerB", flag);
    }
#else
    if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	flag = IDASpilsSetPreconditionerBS(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn_sens : NULL,
		    bprecsolvefn_sens);
	SCHECK_FLAG ("IDASpilsSetPreconditionerBS", flag);
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    } else {
	flag = IDASpilsSetPreconditionerB(
		    mem,
		    which,
		    Bool_val(vset_precsetup) ? bprecsetupfn : NULL,
		    bprecsolvefn);
	SCHECK_FLAG ("IDASpilsSetPreconditionerB", flag);
    }
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_set_jac_times(value vparent,
					      value vwhich,
					      value vhas_setup,
					      value vhas_times,
					      value vusesens)
{
    CAMLparam5(vparent, vwhich, vhas_setup, vhas_times, vusesens);
    void *mem = IDA_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    int flag;

#if   400 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vusesens)) {
	flag = IDASetJacTimesBS(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn_sens : NULL,
			Bool_val(vhas_times) ? bjactimesfn_sens : NULL);
	SCHECK_FLAG ("IDASetJacTimesBS", flag);
    } else {
	flag = IDASetJacTimesB(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn : NULL,
			Bool_val(vhas_times) ? bjactimesfn : NULL);
	SCHECK_FLAG ("IDASetJacTimesB", flag);
    }

#elif 300 <= SUNDIALS_LIB_VERSION
    if (Bool_val(vusesens)) {
	flag = IDASpilsSetJacTimesBS(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn_sens : NULL,
			Bool_val(vhas_times) ? bjactimesfn_sens : NULL);
	SCHECK_FLAG ("IDASpilsSetJacTimesBS", flag);
    } else {
	flag = IDASpilsSetJacTimesB(mem, which,
			Bool_val(vhas_setup) ? bjacsetupfn : NULL,
			Bool_val(vhas_times) ? bjactimesfn : NULL);
	SCHECK_FLAG ("IDASpilsSetJacTimesB", flag);
    }
#else
    if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	flag = IDASpilsSetJacTimesVecFnBS(mem, which,
			Bool_val (vhas_times) ? bjactimesfn_sens : NULL);
	SCHECK_FLAG ("IDASpilsSetJacTimesVecFnBS", flag);
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    } else {
	flag = IDASpilsSetJacTimesVecFnB(mem, which,
			Bool_val (vhas_times) ? bjactimesfn : NULL);
	SCHECK_FLAG ("IDASpilsSetJacTimesVecFnB", flag);
    }
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_jac_times_resfn(value vparent,
						  value vwhich,
						  value vhas_resfn)
{
    CAMLparam3(vparent, vwhich, vhas_resfn);
#if 530 <= SUNDIALS_LIB_VERSION
    void *mem = IDA_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    IDAResFn resfn = Bool_val (vhas_resfn) ? bjactimesresfn : NULL;

    int flag = IDASetJacTimesResFnB(mem, which, resfn);
    CHECK_LS_FLAG("IDASetJacTimesResFnB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_idas_adj_dls_dense(value vparent, value vwhich,
				    value vnb, value vset_jac, value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_jac, vusesens);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    long nbeqs = Long_val(vnb);
    int which = Int_val(vwhich);
    int flag;

    flag = IDADenseB (ida_mem, which, nbeqs);
    SCHECK_FLAG ("IDADenseB", flag);
    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	    flag = IDADlsSetDenseJacFnBS(ida_mem, which, bjacfn_withsens);
	    SCHECK_FLAG("IDADlsSetDenseJacFnBS", flag);
#else
	    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	} else {
	    flag = IDADlsSetDenseJacFnB(ida_mem, which, bjacfn_nosens);
	    SCHECK_FLAG("IDADlsSetDenseJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_dls_lapack_dense(value vparent, value vwhich,
					   value vnb, value vset_jac,
					   value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_jac, vusesens);
#if defined(SUNDIALS_ML_LAPACK) && (SUNDIALS_LIB_VERSION >= 262) && (SUNDIALS_LIB_VERSION < 300)
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    long nbeqs = Long_val(vnb);
    int which = Int_val(vwhich);
    int flag;

    flag = IDALapackDenseB (ida_mem, which, nbeqs);
    SCHECK_FLAG ("IDALapackDenseB", flag);
    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
	    flag = IDADlsSetDenseJacFnBS(ida_mem, which, bjacfn_withsens);
	    SCHECK_FLAG("IDADlsSetDenseJacFnBS", flag);
	} else {
	    flag = IDADlsSetDenseJacFnB(ida_mem, which, bjacfn_nosens);
	    SCHECK_FLAG("IDADlsSetDenseJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_dls_band (value vparent_which, value vsizes,
				    value vset_jac, value vusesens)
{
    CAMLparam4(vparent_which, vsizes, vset_jac, vusesens);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(Field(vsizes, 0));
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = IDABandB (ida_mem, which, nbeqs,
		     Long_val (Field(vsizes, 1)),
		     Long_val (Field(vsizes, 2)));
    SCHECK_FLAG ("IDABandB", flag);
    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	    flag = IDADlsSetBandJacFnBS(ida_mem, which, bbandjacfn_withsens);
	    SCHECK_FLAG("IDADlsSetBandJacFnBS", flag);
#else
	    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	} else {
	    flag = IDADlsSetBandJacFnB(ida_mem, which, bbandjacfn_nosens);
	    SCHECK_FLAG("IDADlsSetBandJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_dls_lapack_band (value vparent_which, value vsizes,
					   value vset_jac, value vusesens)
{
    CAMLparam4(vparent_which, vsizes, vset_jac, vusesens);
#if defined(SUNDIALS_ML_LAPACK) && (SUNDIALS_LIB_VERSION >= 262) && (SUNDIALS_LIB_VERSION < 300)
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(Field(vsizes, 0));
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = IDALapackBandB(ida_mem, which, nbeqs,
			  Long_val (Field(vsizes, 1)),
			  Long_val (Field(vsizes, 2)));
    SCHECK_FLAG ("IDALapackBandB", flag);

    if (Bool_val (vset_jac)) {
	if (Bool_val(vusesens)) {
	    flag = IDADlsSetBandJacFnBS(ida_mem, which, bbandjacfn_withsens);
	    SCHECK_FLAG("IDADlsSetBandJacFnBS", flag);
	} else {
	    flag = IDADlsSetBandJacFnB(ida_mem, which, bbandjacfn_nosens);
	    SCHECK_FLAG("IDADlsSetBandJacFnB", flag);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_linear_solver (value vparent_which,
		    value vlsolv, value vojmat, value vhasjac, value vusesens)
{
    CAMLparam5(vparent_which, vlsolv, vojmat, vhasjac, vusesens);
#if 400 <= SUNDIALS_LIB_VERSION
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparent_which, 0));
    int which = Int_val(Field(vparent_which, 1));
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = IDASetLinearSolverB(ida_mem, which, lsolv, jmat);
    SCHECK_FLAG ("IDASetLinearSolverB", flag);

    if (Bool_val (vusesens)) {
	flag = IDASetJacFnBS(ida_mem, which,
			     Bool_val(vhasjac) ? bjacfn_withsens : NULL);
	SCHECK_FLAG("IDASetJacFnBS", flag);
    } else {
	flag = IDASetJacFnB(ida_mem, which,
			     Bool_val(vhasjac) ? bjacfn_nosens : NULL);
	SCHECK_FLAG("IDASetJacFnB", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_dls_set_linear_solver (value vparent_which,
						 value vlsolv, value vjmat,
						 value vhasjac, value vusesens)
{
    CAMLparam5(vparent_which, vlsolv, vjmat, vhasjac, vusesens);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparent_which, 0));
    int which = Int_val(Field(vparent_which, 1));
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = IDADlsSetLinearSolverB(ida_mem, which, lsolv, jmat);
    SCHECK_FLAG ("IDADlsSetLinearSolverB", flag);

    if (Bool_val (vusesens)) {
	flag = IDADlsSetJacFnBS(ida_mem, which,
			        Bool_val(vhasjac) ? bjacfn_withsens : NULL);
	SCHECK_FLAG("IDADlsSetJacFnBS", flag);
    } else {
	flag = IDADlsSetJacFnB(ida_mem, which,
			       Bool_val(vhasjac) ? bjacfn_nosens : NULL);
	SCHECK_FLAG("IDADlsSetJacFnB", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_spils_set_linear_solver (value vparent, value vwhich,
						   value vlsolv)
{
    CAMLparam3(vparent, vwhich, vlsolv);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    int flag;

    flag = IDASpilsSetLinearSolverB(ida_mem, which, lsolv);
    SCHECK_FLAG ("IDASpilsSetLinearSolverB", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_init_backward(value vparent, value weakref,
					    value vtb0, value vy0, value vyp0,
					    value vwithsens)
{
    CAMLparam5(vparent, weakref, vtb0, vy0, vyp0);
    CAMLxparam1(vwithsens);
    CAMLlocal2(r, vida_mem);
    CAMLlocal2(vlmm, viter);
    int flag, which;
    void *parent = IDA_MEM_FROM_ML(vparent);

    sunrealtype tb0 = Double_val(vtb0);
    N_Vector y0 = NVEC_VAL(vy0);
    N_Vector yp0 = NVEC_VAL(vyp0);

    flag = IDACreateB(parent, &which);
    if (flag != IDA_SUCCESS) {
	SCHECK_FLAG("IDACreateB", flag);
    }

    if (Bool_val(vwithsens)) {
	flag = IDAInitBS(parent, which, bresfn_sens, tb0, y0, yp0);
	if (flag != IDA_SUCCESS) {
	    SCHECK_FLAG("IDAInitBS", flag);
	}
    } else {
	flag = IDAInitB(parent, which, bresfn, tb0, y0, yp0);
	if (flag != IDA_SUCCESS) {
	    SCHECK_FLAG("IDAInitB", flag);
	}
    }

    vida_mem = sunml_wrap_session_pointer(IDAGetAdjIDABmem(parent, which));

    value *backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	caml_raise_out_of_memory();
    }
    IDASetUserDataB (parent, which, backref);

    r = caml_alloc_tuple (3);
    Store_field (r, 0, vida_mem);
    Store_field (r, 1, Val_int(which));
    Store_field (r, 2, (value)backref);

    CAMLreturn(r);
}

BYTE_STUB6(sunml_idas_adj_init_backward)

CAMLprim value sunml_idas_adj_set_id (value vparent, value vwhich, value vid)
{
    CAMLparam3 (vparent, vwhich, vid);
    int flag = IDASetIdB (IDA_MEM_FROM_ML (vparent), Int_val (vwhich),
			  NVEC_VAL (vid));
    SCHECK_FLAG ("IDASetIdB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_set_suppress_alg (value vparent, value vwhich,
					    value vsuppress)
{
    CAMLparam3 (vparent, vwhich, vsuppress);
    int flag = IDASetSuppressAlgB (IDA_MEM_FROM_ML (vparent), Int_val (vwhich),
				   Bool_val (vsuppress));
    SCHECK_FLAG ("IDASetSuppressAlgB", flag);
    CAMLreturn (Val_unit);
}

#if SUNDIALS_LIB_VERSION < 260
// Work around a bug in Sundials 2.5.0
SUNDIALS_EXPORT int IDAAdjSetNoSensi (void *ida_mem);
#endif

CAMLprim value sunml_idas_adj_set_no_sensi (value vsession)
{
    CAMLparam1 (vsession);
    int flag = IDAAdjSetNoSensi (IDA_MEM_FROM_ML (vsession));
    SCHECK_FLAG ("IDAAdjSetNoSensi", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_calc_ic (value vparent, value vwhich,
				   value vtB, value vyB0, value vypB0)
{
    CAMLparam5 (vparent, vwhich, vtB, vyB0, vypB0);
    int flag = IDACalcICB (IDA_MEM_FROM_ML (vparent), Int_val (vwhich),
			   Double_val (vtB), NVEC_VAL (vyB0), NVEC_VAL (vypB0));
    SCHECK_FLAG ("IDACalcICB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_calc_ic_sens (value vparent, value vwhich,
					value vtB, value vyB0, value vypB0,
					value vyS0, value vypS0)
{
    CAMLparam5 (vparent, vwhich, vtB, vyB0, vypB0);
    CAMLxparam2 (vyS0, vypS0);
    N_Vector *yS0 = sunml_nvector_array_alloc (vyS0);
    N_Vector *ypS0 = sunml_nvector_array_alloc (vypS0);
    int flag;

    flag = IDACalcICBS (IDA_MEM_FROM_ML (vparent), Int_val (vwhich),
			Double_val (vtB), NVEC_VAL (vyB0), NVEC_VAL (vypB0),
			yS0, ypS0);
    sunml_nvector_array_free (yS0);
    sunml_nvector_array_free (ypS0);
    SCHECK_FLAG ("IDACalcICBS", flag);

    CAMLreturn (Val_unit);
}

BYTE_STUB7 (sunml_idas_adj_calc_ic_sens)

CAMLprim value sunml_idas_adj_get_consistent_ic (value vparent, value vwhich,
					     value vyB, value vypB)
{
    CAMLparam4 (vparent, vwhich, vyB, vypB);
    N_Vector yB, ypB;

    yB  = Is_block (vyB) ? NVEC_VAL (Field (vyB, 0)) : NULL;
    ypB = Is_block (vypB) ? NVEC_VAL (Field (vypB, 0)) : NULL;

    if (yB || ypB) {
	int flag = IDAGetConsistentICB (IDA_MEM_FROM_ML (vparent),
					Int_val (vwhich),
					yB, ypB);
	SCHECK_FLAG ("IDAGetConsistentICB", flag);
    }

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_reinit(value vparent, value vwhich,
				 value vtB0, value vyB0, value vypB0)
{
    CAMLparam5(vparent, vwhich, vtB0, vyB0, vypB0);
    CAMLlocal1(r);
    int flag;
    N_Vector yB0  = NVEC_VAL(vyB0);
    N_Vector ypB0 = NVEC_VAL(vypB0);

    flag = IDAReInitB (IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
		       Double_val(vtB0), yB0, ypB0);
    SCHECK_FLAG("IDAReInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adj_get(value vparent, value vwhich,
			      value vyB, value vypB)
{
    CAMLparam4(vparent, vwhich, vyB, vypB);
    N_Vector yB = NVEC_VAL(vyB);
    N_Vector ypB = NVEC_VAL(vypB);
    sunrealtype tret;

    int flag = IDAGetB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
		       &tret, yB, ypB);
    SCHECK_FLAG("IDAGetB", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_idas_adj_get_y(value vdata, value vt, value vy, value vyp)
{
    CAMLparam4(vdata, vt, vy, vyp);

    N_Vector y_nv = NVEC_VAL(vy);
    N_Vector yp_nv = NVEC_VAL(vyp);

    int flag = IDAGetAdjY(IDA_MEM_FROM_ML(vdata), Double_val(vt), y_nv, yp_nv);
    SCHECK_FLAG("IdaGetAdjY", flag);

    CAMLreturn (Val_unit);
}


/* adjoint/quadrature interface */

CAMLprim value sunml_idas_adjquad_initb(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = IDAQuadInitB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			bquadrhsfn, yqb0);
    SCHECK_FLAG("IDAQuadInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adjquad_initbs(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = IDAQuadInitBS(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			 bquadrhsfn_sens, yqb0);
    SCHECK_FLAG("IDAQuadInitBS", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adjquad_reinit(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = IDAQuadReInitB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich), yqb0);
    SCHECK_FLAG("IDAQuadReInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adjquad_get(value vparent, value vwhich, value vyqb)
{
    CAMLparam3(vparent, vwhich, vyqb);
    CAMLlocal1(r);
    N_Vector yqb = NVEC_VAL(vyqb);
    sunrealtype tret;

    int flag = IDAGetQuadB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			     &tret, yqb);
    SCHECK_FLAG("IDAGetQuadB", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value sunml_idas_adjquad_sv_tolerances(value vparent, value vwhich,
					    value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = IDAQuadSVtolerancesB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				      Double_val(vreltol), atol_nv);
    SCHECK_FLAG("IDAQuadSVtolerancesB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adjquad_set_err_con(value vparent, value vwhich,
					  value verrconq)
{
    CAMLparam3(vparent, vwhich, verrconq);
    int flag;
    
    flag = IDASetQuadErrConB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			     Bool_val(verrconq));
    SCHECK_FLAG("IDASetQuadErrConB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_idas_adjquad_ss_tolerances(value vparent, value vwhich,
					    value reltol, value abstol)
{
    CAMLparam4(vparent, vwhich, reltol, abstol);

    int flag = IDAQuadSStolerancesB(IDA_MEM_FROM_ML(vparent),
				    Int_val(vwhich), Double_val(reltol),
				    Double_val(abstol));
    SCHECK_FLAG("IDAQuadSStolerancesB", flag);

    CAMLreturn (Val_unit);
}

