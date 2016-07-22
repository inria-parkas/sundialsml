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
#include <idas/idas_dense.h>
#include <idas/idas_direct.h>
#include <idas/idas_band.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_spils.h>

#ifdef SUNDIALS_ML_LAPACK
#include <idas/idas_lapack.h>
#endif

#include "../lsolvers/spils_ml.h"
#include "../ida/ida_ml.h"
#include "idas_ml.h"
#include "../sundials/sundials_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../lsolvers/dls_ml.h"

#define MAX_ERRMSG_LEN 256

CAMLprim value c_idas_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (IDAS, exns);
    CAMLreturn (Val_unit);
}

/* Interface with nvectors */

CAMLprim value c_idas_alloc_nvector_array(value vn)
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

// NB: Normally, we should worry about relinquishing the elements of vy
// after we are finished using them (so as not to block the GC), but we
// instead make the assumption that these elements come from 'within'
// Sundials and thus that they would anyway not be GC-ed.
void idas_wrap_to_nvector_table(int n, value vy, N_Vector *y)
{
    int i;
    for (i = 0; i < n; ++i) {
	Store_field(vy, i, NVEC_BACKLINK(y[i]));
    }
}

static N_Vector *nvector_table_to_array(value vtable)
{
    int ns = Wosize_val (vtable); /* vtable : nvector array */
    N_Vector *r = calloc(ns + 1, sizeof(N_Vector));
    int i;

    for (i=0; i < ns; ++i) {
	r[i] = NVEC_VAL(Field(vtable, i));
    }
    r[ns] = NULL;

    return r;
}

static void free_nvector_array(N_Vector *nvarr)
{
    free(nvarr);
}


void idas_ml_check_flag(const char *call, int flag)
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
	caml_raise_constant(IDA_EXN(LinearSetupFailure));

    case IDA_LSOLVE_FAIL:
	caml_raise_constant(IDA_EXN(LinearSolveFailure));

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

value idas_make_jac_arg(realtype t, N_Vector y, N_Vector yp,
			N_Vector yb, N_Vector ypb, N_Vector resb,
			realtype coef, value tmp)
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

static int quadrhsfn(realtype t, N_Vector y, N_Vector yp, N_Vector rhsQ,
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
    // this call, afterward that memory goes back to cvode. These bigarrays
    // must not be retained by closure_quadrhsfn! If it wants a permanent
    // copy, then it has to make it manually.

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (IDAS_QUADRHSFN_FROM_EXT (sensext), 4, args);

    CAMLreturnT (int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}


static int sensresfn(int Ns, realtype t,
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
    idas_wrap_to_nvector_table (Ns, IDAS_SENSARRAY1_FROM_EXT(sensext), yS);
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_SENSP,
	         IDAS_SENSARRAY2_FROM_EXT(sensext));
    idas_wrap_to_nvector_table (Ns, IDAS_SENSARRAY2_FROM_EXT(sensext), ypS);
    Store_field (args, RECORD_IDAS_SENSRESFN_ARGS_TMP,
		 ida_make_triple_tmp (tmp1, tmp2, tmp3));

    idas_wrap_to_nvector_table (Ns, IDAS_SENSARRAY3_FROM_EXT(sensext), resvalS);

    value r = caml_callback2_exn (IDAS_SENSRESFN_FROM_EXT(sensext), args,
				  IDAS_SENSARRAY3_FROM_EXT(sensext));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int quadsensrhsfn(int ns, realtype t, N_Vector yy, N_Vector yp,
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
    idas_wrap_to_nvector_table (ns, IDAS_SENSARRAY1_FROM_EXT(sensext), yyS);
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_SENS,
		 IDAS_SENSARRAY1_FROM_EXT(sensext));
    idas_wrap_to_nvector_table (ns, IDAS_SENSARRAY2_FROM_EXT(sensext), ypS);
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_SENSP,
		 IDAS_SENSARRAY2_FROM_EXT(sensext));
    Store_field (args, RECORD_IDAS_QUADSENSRHSFN_ARGS_TMP,
		 ida_make_triple_tmp (tmp1, tmp2, tmp3));

    idas_wrap_to_nvector_table (ns,
	    IDAS_SENSARRAY3_FROM_EXT(sensext), rhsvalQS);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(IDAS_QUADSENSRHSFN_FROM_EXT(sensext),
				 args, IDAS_SENSARRAY3_FROM_EXT(sensext));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}


static int bresfn(realtype t, N_Vector y, N_Vector yp,
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

static int bresfn_sens(realtype t, N_Vector y, N_Vector yp,
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

    idas_wrap_to_nvector_table (ns, IDAS_BSENSARRAY1_FROM_EXT (bsensext), yS);
    idas_wrap_to_nvector_table (ns, IDAS_BSENSARRAY2_FROM_EXT (bsensext), ypS);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (IDAS_BRESFN_SENS_FROM_EXT(bsensext), args);
    if (! Is_exception_result (r))
	r = caml_callback3_exn (r, IDAS_BSENSARRAY1_FROM_EXT (bsensext),
				IDAS_BSENSARRAY2_FROM_EXT (bsensext),
				NVEC_BACKLINK (resvalB));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bprecsetupfn(
    realtype t,
    N_Vector yy, N_Vector yp,
    N_Vector yB, N_Vector ypB, N_Vector resvalB,
    realtype cjB, void *user_data,
    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocal3(session, cb, arg);

    arg = idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB,
			    ida_make_triple_tmp(tmp1B, tmp2B, tmp3B));

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
static int bprecsetupfn_sens(
    realtype t,
    N_Vector yy, N_Vector yp,
    N_Vector *yyS, N_Vector *ypS,
    N_Vector yB, N_Vector ypB, N_Vector resvalB,
    realtype cjB, void *user_data,
    N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocal3(session, cb, bsensext);
    CAMLlocalN(args, 3);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));

    args[0] = idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB,
			        ida_make_triple_tmp(tmp1B, tmp2B, tmp3B));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    idas_wrap_to_nvector_table (ns, args[1], yyS);
    idas_wrap_to_nvector_table (ns, args[2], ypS);

    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDAS_BSPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (cb, args[0], args[1], args[2]);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}
#endif

static int bprecsolvefn(
    realtype t,
    N_Vector yy, N_Vector yp,
    N_Vector yB, N_Vector ypB, N_Vector resvalB,
    N_Vector rvecB, N_Vector zvecB,
    realtype cjB, realtype deltaB,
    void *user_data, N_Vector tmpB)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal2(session, cb);

    args[0] = idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB,
			        NVEC_BACKLINK(tmpB));
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
static int bprecsolvefn_sens(
    realtype t,
    N_Vector yy, N_Vector yp,
    N_Vector *yyS, N_Vector *ypS,
    N_Vector yB, N_Vector ypB, N_Vector resvalB,
    N_Vector rvecB, N_Vector zvecB,
    realtype cjB, realtype deltaB,
    void *user_data, N_Vector tmpB)
{
    CAMLparam0();
    CAMLlocalN(args, 6);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    args[0] = idas_make_jac_arg(t, yy, yp, yB, ypB, resvalB, cjB,
			        NVEC_BACKLINK(tmpB));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    idas_wrap_to_nvector_table (ns, args[1], yyS);
    idas_wrap_to_nvector_table (ns, args[2], ypS);

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

static int bjactimesfn(
    realtype t,
    N_Vector yy, N_Vector yp,
    N_Vector yyB, N_Vector ypB, N_Vector resvalB,
    N_Vector vB, N_Vector JvB,
    realtype cjB, void *user_data,
    N_Vector tmp1B, N_Vector tmp2B)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal2(session, cb);

    args[0] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_double_tmp (tmp1B, tmp2B));
    args[1] = NVEC_BACKLINK(vB);
    args[2] = NVEC_BACKLINK(JvB);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (cb, args[0], args[1], args[2]);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bjactimesfn_sens(
    realtype t,
    N_Vector yy, N_Vector yp,
    N_Vector *yyS, N_Vector *ypS,
    N_Vector yyB, N_Vector ypB, N_Vector resvalB,
    N_Vector vB, N_Vector JvB,
    realtype cjB, void *user_data,
    N_Vector tmp1B, N_Vector tmp2B)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    CAMLlocal3(session, bsensext, cb);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    args[0] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_double_tmp (tmp1B, tmp2B));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    idas_wrap_to_nvector_table (ns, args[1], yyS);
    idas_wrap_to_nvector_table (ns, args[2], ypS);

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

static int bjacfn_nosens(long int NeqB, realtype t,
			 realtype cjB, N_Vector yy, N_Vector yp,
			 N_Vector yyB, N_Vector ypB,
			 N_Vector resvalB,
			 DlsMat JacB, void *user_data,
			 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, c_dls_dense_wrap(JacB, 0));
	Store_field(cb, 1, dmat);
    }

    args[0] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (Field(cb, 0), args[0], args[1]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bjacfn_withsens(long int NeqB, realtype t,
			   realtype cjB, N_Vector yy, N_Vector yp,
			   N_Vector *yS, N_Vector *ypS,
			   N_Vector yyB, N_Vector ypB,
			   N_Vector resvalB,
			   DlsMat JacB, void *user_data,
			   N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
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
	Store_some(dmat, c_dls_dense_wrap(JacB, 0));
	Store_field(cb, 1, dmat);
    }

    args[0] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    idas_wrap_to_nvector_table (ns, args[1], yS);
    idas_wrap_to_nvector_table (ns, args[2], ypS);

    args[3] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

static int bbandjacfn_nosens(long int NeqB, long int mupperb, long int mlowerb,
			     realtype t, realtype cjB,
			     N_Vector yy, N_Vector yp,
			     N_Vector yyB, N_Vector ypB,
			     N_Vector resvalB, DlsMat JacB,
			     void *user_data,
			     N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field(cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, c_dls_band_wrap(JacB, 0));
	Store_field(cb, 1, bmat);
    }

    args[0] = caml_alloc_tuple(RECORD_IDAS_ADJ_BANDRANGE_SIZE);
    Store_field(args[0], RECORD_IDAS_ADJ_BANDRANGE_MUPPER, Val_long(mupperb));
    Store_field(args[0], RECORD_IDAS_ADJ_BANDRANGE_MLOWER, Val_long(mlowerb));
    args[1] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_triple_tmp(tmp1B, tmp2B, tmp3B));
    args[2] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (Field(cb, 0), args[0], args[1], args[2]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

#if SUNDIALS_LIB_VERSION >= 260
static int bbandjacfn_withsens(long int NeqB,
			       long int mupperb, long int mlowerb,
			       realtype t, realtype cjB,
			       N_Vector yy, N_Vector yp,
			       N_Vector *yS, N_Vector *ypS,
			       N_Vector yyB, N_Vector ypB,
			       N_Vector resvalB, DlsMat JacB,
			       void *user_data,
			       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    CAMLlocal4(session, bsensext, cb, bmat);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field(cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, c_dls_band_wrap(JacB, 0));
	Store_field(cb, 1, bmat);
    }

    args[0] = caml_alloc_tuple(RECORD_IDAS_ADJ_BANDRANGE_SIZE);
    Store_field(args[0], RECORD_IDAS_ADJ_BANDRANGE_MUPPER, Val_long(mupperb));
    Store_field(args[0], RECORD_IDAS_ADJ_BANDRANGE_MLOWER, Val_long(mlowerb));

    args[1] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_triple_tmp(tmp1B, tmp2B, tmp3B));

    ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[2] = IDAS_BSENSARRAY1_FROM_EXT (bsensext);
    args[3] = IDAS_BSENSARRAY2_FROM_EXT (bsensext);
    idas_wrap_to_nvector_table (ns, args[2], yS);
    idas_wrap_to_nvector_table (ns, args[3], ypS);

    args[4] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 5, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

static int bquadrhsfn(realtype t, N_Vector y, N_Vector yp,
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

static int bquadrhsfn_sens(realtype t, N_Vector y, N_Vector yp,
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

    idas_wrap_to_nvector_table (ns, IDAS_BSENSARRAY1_FROM_EXT (sensext), yS);
    idas_wrap_to_nvector_table (ns, IDAS_BSENSARRAY2_FROM_EXT (sensext), ypS);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (IDAS_BQUADRHSFN_SENS_FROM_EXT(sensext), args);
    if (!Is_exception_result (r))
	r = caml_callback3_exn (r, IDAS_BSENSARRAY1_FROM_EXT (sensext),
				IDAS_BSENSARRAY2_FROM_EXT (sensext),
				NVEC_BACKLINK (rhsvalBQS));

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}


/* quadrature interface */

CAMLprim value c_idas_quad_init (value vsession, value vyQ0)
{

    CAMLparam2 (vsession, vyQ0);
    N_Vector yQ0 = NVEC_VAL (vyQ0);
    int flag;

    flag = IDAQuadInit (IDA_MEM_FROM_ML (vsession), quadrhsfn, yQ0);
    SCHECK_FLAG ("IDAQuadInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quad_reinit (value vsession, value vyQ0)
{

    CAMLparam2 (vsession, vyQ0);
    N_Vector yQ0 = NVEC_VAL (vyQ0);
    int flag;

    flag = IDAQuadReInit (IDA_MEM_FROM_ML (vsession), yQ0);
    SCHECK_FLAG ("IDAQuadReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quad_sv_tolerances(value vdata, value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = IDAQuadSVtolerances(IDA_MEM_FROM_ML(vdata),
				   Double_val(reltol), atol_nv);
    SCHECK_FLAG("IDAQuadSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quad_get(value vdata, value vyq)
{
    CAMLparam2(vdata, vyq);
    N_Vector yq = NVEC_VAL(vyq);
    realtype tret;

    int flag = IDAGetQuad(IDA_MEM_FROM_ML(vdata), &tret, yq);
    SCHECK_FLAG("IDAGetQuad", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_idas_quad_get_dky(value vdata, value vt, value vk,
				   value vdkyq)
{
    CAMLparam4(vdata, vt, vk, vdkyq);
    N_Vector dkyq = NVEC_VAL(vdkyq);

    int flag = IDAGetQuadDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkyq);
	    
    SCHECK_FLAG("IDAGetQuadDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quad_get_err_weights(value vdata, value veqweight)
{
    CAMLparam2(vdata, veqweight);
    N_Vector eqweight = NVEC_VAL(veqweight);

    int flag = IDAGetQuadErrWeights(IDA_MEM_FROM_ML(vdata), eqweight);
    SCHECK_FLAG("IDAGetQuadErrWeights", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value c_idas_quad_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = IDASetQuadErrCon(IDA_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("IDASetQuadErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quad_ss_tolerances(value vdata,
					 value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDAQuadSStolerances(IDA_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    SCHECK_FLAG("IDAQuadSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quad_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadNumRhsEvals(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_quad_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadNumErrTestFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_quad_get_stats(value vdata)
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


CAMLprim value c_idas_sens_init(value vdata, value vmethod, value vrhsfn,
				value vyS0, value vypS0)
{
    CAMLparam5(vdata, vmethod, vrhsfn, vyS0, vypS0);
    int ns = Wosize_val (vyS0);	/* vyS0 : nvector array */
    N_Vector *yS0 = nvector_table_to_array(vyS0);
    N_Vector *ypS0 = nvector_table_to_array(vypS0);

    int flag = IDASensInit(IDA_MEM_FROM_ML(vdata), ns,
			   decode_sens_method(vmethod),
			   ((Bool_val(vrhsfn)) ? sensresfn : NULL),
			   yS0, ypS0);
    free_nvector_array(yS0); 
    free_nvector_array(ypS0); 
    SCHECK_FLAG("IDASensInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_reinit(value vdata, value vmethod, value vyS0,
				  value vypS0)
{
    CAMLparam4(vdata, vmethod, vyS0, vypS0);
    CAMLlocal1(r);
    int flag;
    N_Vector *yS0 = nvector_table_to_array(vyS0);
    N_Vector *ypS0 = nvector_table_to_array(vypS0);

    flag = IDASensReInit(IDA_MEM_FROM_ML(vdata),
			 decode_sens_method(vmethod),
			 yS0, ypS0);
    free_nvector_array(yS0);
    free_nvector_array(ypS0);
    SCHECK_FLAG("IDASensReInit", flag);

    CAMLreturn (Val_unit);
}

static void sens_calc_ic (void *ida_mem, value session, int icopt, realtype tout1,
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
    ys  = (Is_block (vys))  ? nvector_table_to_array (Field (vys, 0))  : NULL;
    yps = (Is_block (vyps)) ? nvector_table_to_array (Field (vyps, 0)) : NULL;
    if (ys != NULL || yps != NULL) {
	flag = IDAGetSensConsistentIC (ida_mem, ys, yps);
	CHECK_FLAG ("IDAGetConsistentIC", flag);
	if (ys)  free_nvector_array(ys); 
	if (yps) free_nvector_array(yps); 
    }
    CAMLreturn0;
}

CAMLprim value c_ida_sens_calc_ic_y(value vida_mem, value vy, value vys, value tout1)
{
    CAMLparam4 (vida_mem, vy, vys, tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    sens_calc_ic (ida_mem, vida_mem, IDA_Y_INIT, Double_val (tout1),
		  vy, Val_none, vys, Val_none);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_sens_calc_ic_ya_ydp(value vida_mem, value y, value yp,
					 value ys, value yps, value tout1)
{
    CAMLparam5 (vida_mem, y, yp, ys, yps);
    CAMLxparam1 (tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    sens_calc_ic (ida_mem, vida_mem, IDA_YA_YDP_INIT, Double_val (tout1),
		  y, yp, ys, yps);

    CAMLreturn (Val_unit);
}

BYTE_STUB6(c_ida_sens_calc_ic_ya_ydp)

CAMLprim value c_idas_sens_get(value vdata, value vys)
{
    CAMLparam2(vdata, vys);
    N_Vector *ys = nvector_table_to_array(vys);
    realtype tret;

    int flag = IDAGetSens(IDA_MEM_FROM_ML(vdata), &tret, ys);
    free_nvector_array(ys);
    SCHECK_FLAG("IDAGetSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_idas_sens_get_dky(value vdata, value vt, value vk, value vdkys)
{
    CAMLparam4(vdata, vt, vk, vdkys);
    N_Vector *dkys = nvector_table_to_array(vdkys);

    int flag = IDAGetSensDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkys);
    free_nvector_array(dkys);
    SCHECK_FLAG("IDAGetSensDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_get1(value vdata, value vis, value vys)
{
    CAMLparam3(vdata, vis, vys);
    N_Vector ys = NVEC_VAL(vys);
    realtype tret;

    int flag = IDAGetSens1(IDA_MEM_FROM_ML(vdata), &tret, Int_val(vis), ys);
    SCHECK_FLAG("IDAGetSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_idas_sens_get_dky1(value vdata, value vt, value vk,
				    value vis, value vdkys)
{
    CAMLparam5(vdata, vt, vk, vis, vdkys);
    N_Vector dkys = NVEC_VAL(vdkys);

    int flag = IDAGetSensDky1(IDA_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), Int_val(vis), dkys);
	    
    SCHECK_FLAG("IDAGetSensDky1", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_get_err_weights(value vdata, value vesweight)
{
    CAMLparam2(vdata, vesweight);
    N_Vector *esweight = nvector_table_to_array(vesweight);

    int flag = IDAGetSensErrWeights(IDA_MEM_FROM_ML(vdata), esweight);
    free_nvector_array(esweight);
    SCHECK_FLAG("IDAGetSensErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_set_err_con(value vdata, value verrcons)
{
    CAMLparam2(vdata, verrcons);
    int flag;
    
    flag = IDASetSensErrCon(IDA_MEM_FROM_ML(vdata), Bool_val(verrcons));
    SCHECK_FLAG("IDASetSensErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_ss_tolerances(value vdata,
					 value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDASensSStolerances(IDA_MEM_FROM_ML(vdata),
		 Double_val(reltol), REAL_ARRAY(abstol));
    SCHECK_FLAG("IDASensSStolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_sv_tolerances(value vdata, value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = nvector_table_to_array(abstol);

    int flag = IDASensSVtolerances(IDA_MEM_FROM_ML(vdata),
				   Double_val(reltol), atol_nv);
    free_nvector_array(atol_nv);
    SCHECK_FLAG("IDASensSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASensEEtolerances(IDA_MEM_FROM_ML(vdata));
    SCHECK_FLAG("IDASensEEtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_set_params(value vdata, value vparams)
{
    CAMLparam2(vdata, vparams);
    CAMLlocal3(vp, vpbar, vplist);

    realtype *p = NULL;
    realtype *pbar = NULL;
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

CAMLprim value c_idas_sens_toggle_off(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASensToggleOff(IDA_MEM_FROM_ML(vdata));
    SCHECK_FLAG("IDASensToggleOff", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_set_dq_method(value vdata, value vdqtype,
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

CAMLprim value c_idas_sens_set_max_nonlin_iters(value vdata, value vmaxcors)
{
    CAMLparam2(vdata, vmaxcors);

    int flag = IDASetSensMaxNonlinIters(IDA_MEM_FROM_ML(vdata),
	    Int_val(vmaxcors));
    SCHECK_FLAG("IDASetSensMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_sens_get_num_res_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumResEvals(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumResEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_sens_get_num_res_evals_sens(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetNumResEvalsSens(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetNumResEvalsSens", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_sens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumErrTestFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_sens_get_num_lin_solv_setups(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumLinSolvSetups(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_sens_get_stats(value vdata)
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

CAMLprim value c_idas_sens_get_num_nonlin_solv_iters(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumNonlinSolvIters(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_sens_get_num_nonlin_solv_conv_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetSensNumNonlinSolvConvFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetSensNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_sens_get_nonlin_solv_stats(value vdata)
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


/* sensitivity/quadrature interface */

CAMLprim value c_idas_quadsens_init(value vdata, value vrhsfn, value vyqs0)
{
    CAMLparam3(vdata, vrhsfn, vyqs0);
    N_Vector *yqs0 = nvector_table_to_array(vyqs0);

    int flag = IDAQuadSensInit(IDA_MEM_FROM_ML(vdata),
			       Bool_val (vrhsfn) ? quadsensrhsfn : NULL,
			       yqs0);
    free_nvector_array(yqs0); 
    SCHECK_FLAG("IDAQuadSensInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quadsens_reinit(value vdata, value vyqs0)
{
    CAMLparam2(vdata, vyqs0);
    N_Vector *yqs0 = nvector_table_to_array(vyqs0);

    int flag = IDAQuadSensReInit(IDA_MEM_FROM_ML(vdata), yqs0);
    free_nvector_array(yqs0); 
    SCHECK_FLAG("IDAQuadSensReInit", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value c_idas_quadsens_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = IDASetQuadSensErrCon(IDA_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("IDASetQuadSensErrCon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quadsens_ss_tolerances(value vdata,
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

CAMLprim value c_idas_quadsens_sv_tolerances(value vdata, value reltol,
					     value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = nvector_table_to_array(abstol);

    int flag = IDAQuadSensSVtolerances(IDA_MEM_FROM_ML(vdata),
				       Double_val(reltol), atol_nv);
    free_nvector_array(atol_nv);
    SCHECK_FLAG("IDAQuadSensSVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quadsens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDAQuadSensEEtolerances(IDA_MEM_FROM_ML(vdata));
    SCHECK_FLAG("IDAQuadSensEEtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quadsens_get(value vdata, value vyqs)
{
    CAMLparam2(vdata, vyqs);
    N_Vector *yqs = nvector_table_to_array(vyqs);
    realtype tret;

    int flag = IDAGetQuadSens(IDA_MEM_FROM_ML(vdata), &tret, yqs);
    free_nvector_array(yqs);
    SCHECK_FLAG("IDAGetQuadSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_idas_quadsens_get1(value vdata, value vis, value vyqs)
{
    CAMLparam3(vdata, vis, vyqs);
    N_Vector yqs = NVEC_VAL(vyqs);
    realtype tret;

    int flag = IDAGetQuadSens1(IDA_MEM_FROM_ML(vdata), &tret,
			       Int_val(vis), yqs);
    SCHECK_FLAG("IDAGetQuadSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_idas_quadsens_get_dky(value vdata, value vt, value vk,
				       value vdkyqs)
{
    CAMLparam4(vdata, vt, vk, vdkyqs);
    N_Vector *dkyqs = nvector_table_to_array(vdkyqs);

    int flag = IDAGetQuadSensDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
				 Int_val(vk), dkyqs);
    free_nvector_array(dkyqs);
    SCHECK_FLAG("IDAGetQuadSensDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quadsens_get_dky1(value vdata, value vt, value vk,
					value vis, value vdkyqs)
{
    CAMLparam5(vdata, vt, vk, vis, vdkyqs);
    N_Vector dkyqs = NVEC_VAL(vdkyqs);

    int flag = IDAGetQuadSensDky1(IDA_MEM_FROM_ML(vdata), Double_val(vt),
				    Int_val(vk), Int_val(vis), dkyqs);
    SCHECK_FLAG("IDAGetQuadSensDky1", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quadsens_get_err_weights(value vdata, value veqweights)
{
    CAMLparam2(vdata, veqweights);
    N_Vector *eqweights = nvector_table_to_array(veqweights);

    int flag = IDAGetQuadSensErrWeights(IDA_MEM_FROM_ML(vdata), eqweights);
    free_nvector_array(eqweights);
    SCHECK_FLAG("IDAGetQuadSensErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_quadsens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadSensNumRhsEvals(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_quadsens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = IDAGetQuadSensNumErrTestFails(IDA_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("IDAGetQuadSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_idas_quadsens_get_stats(value vdata)
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

CAMLprim value c_idas_adj_init(value vdata, value vnd, value vinterptype)
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

CAMLprim value c_idas_adj_ss_tolerances(value vparent, value vwhich,
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

CAMLprim value c_idas_adj_spils_spgmr (value vparent, value vwhich,
				       value vmaxl)
{
    CAMLparam3 (vparent, vwhich, vmaxl);
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = IDASpgmrB (ida_mem, which, Int_val (vmaxl));
    SCHECK_FLAG ("IDASpgmrB", flag);
    CAMLreturn (Val_unit);
}

/* For SPGMR only. */
CAMLprim value c_idas_adj_spils_set_max_restarts (value vida_mem, value vwhich,
						  value vmaxr)
{
    CAMLparam3 (vida_mem, vwhich, vmaxr);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASpilsSetMaxRestartsB (ida_mem, Int_val (vwhich), Int_val (vmaxr));
    CHECK_FLAG ("IDASpilsSetMaxRestartsB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_spils_spbcg (value vparent, value vwhich,
				       value vmaxl)
{
    CAMLparam3 (vparent, vwhich, vmaxl);
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = IDASpbcgB (ida_mem, which, Int_val (vmaxl));
    SCHECK_FLAG ("IDASpbcgB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_spils_sptfqmr (value vparent, value vwhich,
					 value vmaxl)
{
    CAMLparam3 (vparent, vwhich, vmaxl);
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = IDASptfqmrB (ida_mem, which, Int_val (vmaxl));
    SCHECK_FLAG ("IDASptfqmrB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_backward_normal(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = IDASolveB(IDA_MEM_FROM_ML(vdata), Double_val(vtbout),
			 IDA_NORMAL);
    SCHECK_FLAG("IDASolveB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_backward_one_step(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = IDASolveB(IDA_MEM_FROM_ML(vdata), Double_val(vtbout),
			 IDA_ONE_STEP);
    SCHECK_FLAG("IDASolveB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_set_max_ord(value vparent, value vwhich,
				      value vmaxord)
{
    CAMLparam3(vparent, vwhich, vmaxord);

    int flag = IDASetMaxOrdB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(vmaxord));
    SCHECK_FLAG("IDASetMaxOrdB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_set_max_num_steps(value vparent, value vwhich,
					    value vmxsteps)
{
    CAMLparam3(vparent, vwhich, vmxsteps);

    int flag = IDASetMaxNumStepsB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				    Long_val(vmxsteps));
    SCHECK_FLAG("IDASetMaxNumStepsB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_set_init_step(value vparent, value vwhich, value vhin)
{
    CAMLparam3(vparent, vwhich, vhin);

    int flag = IDASetInitStepB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Double_val(vhin));
    SCHECK_FLAG("IDASetInitStepB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_set_max_step(value vparent, value vwhich, value vhmax)
{
    CAMLparam3(vparent, vwhich, vhmax);

    int flag = IDASetMaxStepB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			      Double_val(vhmax));
    SCHECK_FLAG("IDASetMaxStepB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_spils_set_gs_type(value vparent, value vwhich,
					    value vgstype)
{
    CAMLparam3(vparent, vwhich, vgstype);

    int flag = IDASpilsSetGSTypeB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				 spils_gs_type(vgstype));
    SCHECK_FLAG("IDASpilsSetGSTypeB", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value c_idas_adj_spils_set_eps_lin(value vparent, value vwhich,
					    value eplifac)
{
    CAMLparam3(vparent, vwhich, eplifac);

    int flag = IDASpilsSetEpsLinB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				  Double_val(eplifac));
    SCHECK_FLAG("IDASpilsSetEpsLinB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_spils_set_maxl(value vparent, value vwhich,
					 value maxl)
{
    CAMLparam3(vparent, vwhich, maxl);

    int flag = IDASpilsSetMaxlB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(maxl));
    SCHECK_FLAG("IDASpilsSetMaxlB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_bsession_finalize(value vdata)
{
    if (IDA_MEM_FROM_ML(vdata) != NULL) {
	value *backref = IDA_BACKREF_FROM_ML(vdata);
	// NB: IDAFree() is *not* called: parents free-up backward problems
	c_sundials_free_value(backref);
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
    realtype tret;
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
	SCHECK_FLAG ("IDASolveF", flag);
    }

    assert (Field (vdata, RECORD_IDA_SESSION_EXN_TEMP) == Val_none);

    ret = caml_alloc_tuple(3);
    Store_field(ret, 0, caml_copy_double(tret));
    Store_field(ret, 1, Val_int(ncheck));
    Store_field(ret, 2, Val_int(solver_result));

    CAMLreturn(ret);
}

CAMLprim value c_idas_adj_forward_normal(value vdata, value vtout,
					 value vy, value vyp)
{
    CAMLparam4(vdata, vtout, vy, vyp);
    CAMLreturn(forward_solve(vdata, vtout, vy, vyp, 0));
}

CAMLprim value c_idas_adj_forward_one_step(value vdata, value vtout,
					   value vy, value vyp)
{
    CAMLparam4(vdata, vtout, vy, vyp);
    CAMLreturn(forward_solve(vdata, vtout, vy, vyp, 1));
}

CAMLprim value c_idas_adj_sv_tolerances(value vparent, value vwhich,
					value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = IDASVtolerancesB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				Double_val(vreltol), atol_nv);
    SCHECK_FLAG("IDASStolerancesB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_spils_set_preconditioner(value vparent,
						   value vwhich,
						   value vset_precsetup,
						   value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_precsetup, vusesens);
    void *mem = IDA_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    int flag;

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
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_spils_set_jac_times_vec_fn(value vparent,
						     value vwhich,
						     value vset_jac,
						     value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_jac, vusesens);
    void *mem = IDA_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    int flag;

    if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	flag = IDASpilsSetJacTimesVecFnBS(mem, which,
			Bool_val (vset_jac) ? bjactimesfn_sens : NULL);
	SCHECK_FLAG ("IDASpilsSetJacTimesVecFnBS", flag);
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    } else {
	flag = IDASpilsSetJacTimesVecFnB(mem, which,
			Bool_val (vset_jac) ? bjactimesfn : NULL);
	SCHECK_FLAG ("IDASpilsSetJacTimesVecFnB", flag);
    }

    CAMLreturn (Val_unit);
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value c_idas_adj_dls_dense(value vparent, value vwhich,
				    value vnb, value vset_jac, value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_jac, vusesens);
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
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_dls_lapack_dense(value vparent, value vwhich,
					   value vnb, value vset_jac,
					   value vusesens)
{
    CAMLparam4(vparent, vwhich, vset_jac, vusesens);
#if defined(SUNDIALS_ML_LAPACK) && (SUNDIALS_LIB_VERSION >= 262)
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    long nbeqs = Long_val(vnb);
    int which = Int_val(vwhich);
    int flag;

    flag = IDALapackDenseB (ida_mem, which, nbeqs);
    SCHECK_FLAG ("IDALapackDenseB", flag);
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

CAMLprim value c_idas_adj_dls_set_dense_jac_fn(value vparent, value vwhich,
					       value vusesens)
{
    CAMLparam3(vparent, vwhich, vusesens);
    int flag;

    if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	flag = IDADlsSetDenseJacFnBS(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				     bjacfn_withsens);
	SCHECK_FLAG("IDADlsSetDenseJacFnBS", flag);
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    } else {
	flag = IDADlsSetDenseJacFnB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				    bjacfn_nosens);
	SCHECK_FLAG("IDADlsSetDenseJacFnB", flag);
    }
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_dls_clear_dense_jac_fn(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);
    int flag = IDADlsSetDenseJacFnB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				   NULL);
    SCHECK_FLAG("IDADlsSetDenseJacFnB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_dls_band (value vparent_which, value vnb,
				    value vbandrange,
				    value vset_jac, value vusesens)
{
    CAMLparam5(vparent_which, vnb, vbandrange, vset_jac, vusesens);
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(vnb);
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = IDABandB (ida_mem, which, nbeqs,
		     Long_val (Field(vbandrange,
				     RECORD_IDAS_ADJ_BANDRANGE_MUPPER)),
		     Long_val (Field(vbandrange,
				     RECORD_IDAS_ADJ_BANDRANGE_MLOWER)));
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
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_dls_lapack_band (value vparent_which, value vnb,
					   value vbandrange,
					   value vset_jac, value vusesens)
{
    CAMLparam5(vparent_which, vnb, vbandrange, vset_jac, vusesens);
#if defined(SUNDIALS_ML_LAPACK) && (SUNDIALS_LIB_VERSION >= 262)
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(vnb);
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = IDALapackBandB(ida_mem, which, nbeqs,
			  Long_val (Field(vbandrange,
					  RECORD_IDAS_ADJ_BANDRANGE_MUPPER)),
			  Long_val (Field(vbandrange,
					  RECORD_IDAS_ADJ_BANDRANGE_MLOWER)));
    SCHECK_FLAG ("IDALapackBandB", flag);

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

CAMLprim value c_idas_adj_dls_set_band_jac_fn(value vparent, value vwhich,
					      value vusesens)
{
    CAMLparam3(vparent, vwhich, vusesens);
    int flag;

    if (Bool_val(vusesens)) {
#if SUNDIALS_LIB_VERSION >= 260
	flag = IDADlsSetBandJacFnBS(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				    bbandjacfn_withsens);
	SCHECK_FLAG("IDADlsSetBandJacFnBS", flag);
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    } else {
	flag = IDADlsSetBandJacFnB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				   bbandjacfn_nosens);
	SCHECK_FLAG("IDADlsSetBandJacFnB", flag);
    }
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_dls_clear_band_jac_fn(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);
    int flag = IDADlsSetBandJacFnB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				  NULL);
    SCHECK_FLAG("IDADlsSetBandJacFnB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_init_backward(value vparent, value weakref,
					value vtb0, value vy0, value vyp0,
					value vwithsens)
{
    CAMLparam5(vparent, weakref, vtb0, vy0, vyp0);
    CAMLxparam1(vwithsens);
    CAMLlocal2(r, vida_mem);
    CAMLlocal2(vlmm, viter);
    int flag, which;
    void *parent = IDA_MEM_FROM_ML(vparent);

    realtype tb0 = Double_val(vtb0);
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

    vida_mem = caml_alloc_final(1, NULL, sizeof(void *), sizeof(void *) * 15);
    IDA_MEM(vida_mem) = IDAGetAdjIDABmem(parent, which);

    value *backref = c_sundials_malloc_value(weakref);
    if (backref == NULL) {
	caml_raise_out_of_memory();
    }
    IDASetUserDataB (parent, which, backref);

    r = caml_alloc_tuple (4);
    Store_field (r, 0, vida_mem);
    Store_field (r, 1, Val_int(which));
    Store_field (r, 2, (value)backref);
    Store_field (r, 3, 0); // no err_file = NULL; note OCaml doesn't
			   // (seem to) support architectures where
			   // 0 != (value)(void*)NULL.

    CAMLreturn(r);
}

BYTE_STUB6(c_idas_adj_init_backward)

CAMLprim value c_idas_adj_set_id (value vparent, value vwhich, value vid)
{
    CAMLparam3 (vparent, vwhich, vid);
    int flag = IDASetIdB (IDA_MEM_FROM_ML (vparent), Int_val (vwhich),
			  NVEC_VAL (vid));
    SCHECK_FLAG ("IDASetIdB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_set_suppress_alg (value vparent, value vwhich,
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

CAMLprim value c_idas_adj_set_no_sensi (value vsession)
{
    CAMLparam1 (vsession);
    int flag = IDAAdjSetNoSensi (IDA_MEM_FROM_ML (vsession));
    SCHECK_FLAG ("IDAAdjSetNoSensi", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_calc_ic (value vparent, value vwhich,
				   value vtB, value vyB0, value vypB0)
{
    CAMLparam5 (vparent, vwhich, vtB, vyB0, vypB0);
    int flag = IDACalcICB (IDA_MEM_FROM_ML (vparent), Int_val (vwhich),
			   Double_val (vtB), NVEC_VAL (vyB0), NVEC_VAL (vypB0));
    SCHECK_FLAG ("IDACalcICB", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adj_calc_ic_sens (value vparent, value vwhich,
					value vtB, value vyB0, value vypB0,
					value vyS0, value vypS0)
{
    CAMLparam5 (vparent, vwhich, vtB, vyB0, vypB0);
    CAMLxparam2 (vyS0, vypS0);
    N_Vector *yS0 = nvector_table_to_array (vyS0);
    N_Vector *ypS0 = nvector_table_to_array (vypS0);
    int flag;

    flag = IDACalcICBS (IDA_MEM_FROM_ML (vparent), Int_val (vwhich),
			Double_val (vtB), NVEC_VAL (vyB0), NVEC_VAL (vypB0),
			yS0, ypS0);
    free_nvector_array (yS0);
    free_nvector_array (ypS0);
    SCHECK_FLAG ("IDACalcICBS", flag);

    CAMLreturn (Val_unit);
}

BYTE_STUB7 (c_idas_adj_calc_ic_sens)

CAMLprim value c_idas_adj_get_consistent_ic (value vparent, value vwhich,
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

CAMLprim value c_idas_adj_reinit(value vparent, value vwhich,
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

CAMLprim value c_idas_adj_get(value vparent, value vwhich,
			      value vyB, value vypB)
{
    CAMLparam4(vparent, vwhich, vyB, vypB);
    N_Vector yB = NVEC_VAL(vyB);
    N_Vector ypB = NVEC_VAL(vypB);
    realtype tret;

    int flag = IDAGetB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
		       &tret, yB, ypB);
    SCHECK_FLAG("IDAGetB", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_idas_adj_get_y(value vdata, value vt, value vy, value vyp)
{
    CAMLparam4(vdata, vt, vy, vyp);

    N_Vector y_nv = NVEC_VAL(vy);
    N_Vector yp_nv = NVEC_VAL(vyp);

    int flag = IDAGetAdjY(IDA_MEM_FROM_ML(vdata), Double_val(vt), y_nv, yp_nv);
    SCHECK_FLAG("IdaGetAdjY", flag);

    CAMLreturn (Val_unit);
}


/* adjoint/quadrature interface */

CAMLprim value c_idas_adjquad_initb(value vparent, value vwhich, value vyqb0)
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

CAMLprim value c_idas_adjquad_initbs(value vparent, value vwhich, value vyqb0)
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

CAMLprim value c_idas_adjquad_reinit(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = IDAQuadReInitB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich), yqb0);
    SCHECK_FLAG("IDAQuadReInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adjquad_get(value vparent, value vwhich, value vyqb)
{
    CAMLparam3(vparent, vwhich, vyqb);
    CAMLlocal1(r);
    N_Vector yqb = NVEC_VAL(vyqb);
    realtype tret;

    int flag = IDAGetQuadB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			     &tret, yqb);
    SCHECK_FLAG("IDAGetQuadB", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_idas_adjquad_sv_tolerances(value vparent, value vwhich,
					    value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = IDAQuadSVtolerancesB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
				      Double_val(vreltol), atol_nv);
    SCHECK_FLAG("IDAQuadSVtolerancesB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adjquad_set_err_con(value vparent, value vwhich,
					  value verrconq)
{
    CAMLparam3(vparent, vwhich, verrconq);
    int flag;
    
    flag = IDASetQuadErrConB(IDA_MEM_FROM_ML(vparent), Int_val(vwhich),
			     Bool_val(verrconq));
    SCHECK_FLAG("IDASetQuadErrConB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_adjquad_ss_tolerances(value vparent, value vwhich,
					    value reltol, value abstol)
{
    CAMLparam4(vparent, vwhich, reltol, abstol);

    int flag = IDAQuadSStolerancesB(IDA_MEM_FROM_ML(vparent),
				    Int_val(vwhich), Double_val(reltol),
				    Double_val(abstol));
    SCHECK_FLAG("IDAQuadSStolerancesB", flag);

    CAMLreturn (Val_unit);
}

