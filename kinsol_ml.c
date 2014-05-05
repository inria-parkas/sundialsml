/***********************************************************************
 *                                                                     *
 *               OCaml interface to (serial) Sundials                  *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

/* Sundials interface functions that do not involve NVectors. */

#include <kinsol/kinsol.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

/* linear solvers */
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_bbdpre.h>
#include <kinsol/kinsol_spils.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <kinsol/kinsol_lapack.h>
#endif

#include "spils_ml.h"
#include "kinsol_ml.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

void kinsol_ml_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == KIN_SUCCESS
	|| flag == KIN_INITIAL_GUESS_OK
	|| flag == KIN_STEP_LT_STPTOL
	|| flag == KIN_WARNING) return;

    switch (flag) {
    case KIN_MEM_FAIL:
	caml_raise_out_of_memory();

    case KIN_ILL_INPUT:
        caml_raise_constant(*caml_named_value("kinsol_IllInput"));

    case KIN_LINESEARCH_NONCONV:
        caml_raise_constant(*caml_named_value(
		    "kinsol_LineSearchNonConvergence"));

    case KIN_MAXITER_REACHED:
        caml_raise_constant(*caml_named_value("kinsol_MaxIterationsReached"));

    case KIN_MXNEWT_5X_EXCEEDED:
        caml_raise_constant(*caml_named_value("kinsol_MaxNewtonStepExceeded"));

    case KIN_LINESEARCH_BCFAIL:
        caml_raise_constant(*caml_named_value(
		    "kinsol_LineSearchBetaConditionFailure"));

    case KIN_LINSOLV_NO_RECOVERY:
        caml_raise_constant(*caml_named_value("kinsol_LinearSolverNoRecovery"));

    case KIN_LINIT_FAIL:
        caml_raise_constant(*caml_named_value(
		    "kinsol_LinearSolverInitFailure"));

    case KIN_LSETUP_FAIL:
        caml_raise_constant(*caml_named_value("kinsol_LinearSetupFailure"));

    case KIN_LSOLVE_FAIL:
        caml_raise_constant(*caml_named_value("kinsol_LinearSolverFailure"));

    case KIN_SYSFUNC_FAIL:
        caml_raise_constant(*caml_named_value("kinsol_SystemFunctionFailure"));

    case KIN_FIRST_SYSFUNC_ERR:
        caml_raise_constant(*caml_named_value(
		    "kinsol_FirstSystemFunctionFailure"));

    case KIN_REPTD_SYSFUNC_ERR:
        caml_raise_constant(*caml_named_value(
		    "kinsol_RepeatedSystemFunctionFailure"));

    default:
	/* KIN_MEM_NULL, KIN_NO_MALLOC */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: unexpected error code", call);
	caml_failwith(exmsg);
    }
}

/* basic interface */

CAMLprim void c_kinsol_session_finalize(value vdata)
{
    if (KINSOL_MEM_FROM_ML(vdata) != NULL) {
	void *kin_mem = KINSOL_MEM_FROM_ML(vdata);
	value *backref = KINSOL_BACKREF_FROM_ML(vdata);
	KINFree(&kin_mem);
	caml_remove_generational_global_root (backref);
	free (backref);
    }

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_KINSOL_SESSION_ERRFILE));
    if (err_file != NULL) {
	fclose(err_file);
    }
}

/* boiler plate */

CAMLprim void c_kinsol_spils_spgmr(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSpgmr (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSpgmr", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_spbcg(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSpbcg (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSpbcg", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_sptfqmr(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSptfqmr (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSptfqmr", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_set_max_restarts(value vkin_mem, value vmaxrs)
{
    CAMLparam2(vkin_mem, vmaxrs);

    int flag = KINSpilsSetMaxRestarts(KINSOL_MEM_FROM_ML(vkin_mem), Int_val(vmaxrs));
    CHECK_FLAG("KINSetMaxRestarts", flag);

    CAMLreturn0;
}

CAMLprim value c_kinsol_dls_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = KINDlsGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrwLS, &leniwLS);
    CHECK_FLAG("KINDlsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));

    CAMLreturn(r);
}

CAMLprim value c_kinsol_dls_get_num_jac_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINDlsGetNumJacEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINDlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_dls_get_num_func_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINDlsGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINDlsGetNumFuncEvals", flag);

    CAMLreturn(Val_long(r));
}


CAMLprim value c_kinsol_spils_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = KINSpilsGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrw, &leniw);
    CHECK_FLAG("KINSpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
}

CAMLprim value c_kinsol_spils_get_num_lin_iters(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumLinIters(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_conv_fails(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumConvFails(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_prec_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumPrecEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_prec_solves(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumPrecSolves(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_jtimes_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumJtimesEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_func_evals (value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumFuncEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim void c_kinsol_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_KINSOL_SESSION_ERRFILE));

    if (err_file != NULL) {
	fclose(err_file);
	Store_field(vdata, RECORD_KINSOL_SESSION_ERRFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    err_file = fopen(String_val(vpath), mode);
    if (err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = KINSetErrFile(KINSOL_MEM_FROM_ML(vdata), err_file);
    CHECK_FLAG("KINSetErrFile", flag);

    Store_field(vdata, RECORD_KINSOL_SESSION_ERRFILE, Val_long(err_file));

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_info_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    FILE* info_file =
      (FILE *)Long_val(Field(vdata, RECORD_KINSOL_SESSION_INFOFILE));

    if (info_file != NULL) {
	fclose(info_file);
	Store_field(vdata, RECORD_KINSOL_SESSION_INFOFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    info_file = fopen(String_val(vpath), mode);
    if (info_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = KINSetInfoFile(KINSOL_MEM_FROM_ML(vdata), info_file);
    CHECK_FLAG("KINSetInfoFile", flag);

    Store_field(vdata, RECORD_KINSOL_SESSION_INFOFILE, Val_long(info_file));

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_print_level(value vkin_mem, value vplvl)
{
    CAMLparam2(vkin_mem, vplvl);

    int flag = KINSetPrintLevel(KINSOL_MEM_FROM_ML(vkin_mem), Int_val(vplvl));
    CHECK_FLAG("KINSetPrintLevel", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_num_max_iters(value vkin_mem, value vmxiter)
{
    CAMLparam2(vkin_mem, vmxiter);

    int flag = KINSetNumMaxIters(KINSOL_MEM_FROM_ML(vkin_mem), Long_val(vmxiter));
    CHECK_FLAG("KINSetNumMaxIters", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_no_init_setup(value vkin_mem, value vnoinitsetup)
{
    CAMLparam2(vkin_mem, vnoinitsetup);

    int flag = KINSetNoInitSetup(KINSOL_MEM_FROM_ML(vkin_mem), Bool_val(vnoinitsetup));
    CHECK_FLAG("KINSetNoInitSetup", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_max_setup_calls(value vkin_mem, value vmsbset)
{
    CAMLparam2(vkin_mem, vmsbset);

    int flag = KINSetMaxSetupCalls(KINSOL_MEM_FROM_ML(vkin_mem), Long_val(vmsbset));
    CHECK_FLAG("KINSetMaxSetupCalls", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_eta_form(value vkin_mem, value vetachoice)
{
    CAMLparam2(vkin_mem, vetachoice);

    int etachoice;
    if (Is_long(vetachoice)) {
	switch (Int_val(vetachoice)) {
	case VARIANT_KINSOL_ETA_CHOICE1:
	    etachoice = KIN_ETACHOICE1;
	    break;
	}
    } else {
	switch (Tag_val(vetachoice)) {
	case VARIANT_KINSOL_ETA_CHOICE2:
	    etachoice = KIN_ETACHOICE2;
	    break;

	case VARIANT_KINSOL_ETA_CONSTANT:
	    etachoice = KIN_ETACONSTANT;
	    break;
	}
    }

    int flag = KINSetEtaForm(KINSOL_MEM_FROM_ML(vkin_mem), etachoice);
    CHECK_FLAG("KINSetEtaForm", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_eta_const_value(value vkin_mem, value veta)
{
    CAMLparam2(vkin_mem, veta);

    int flag = KINSetEtaConstValue(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(veta));
    CHECK_FLAG("KINSetEtaConstValue", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_eta_params(value vkin_mem, value vegamma, value vealpha)
{
    CAMLparam3(vkin_mem, vegamma, vealpha);

    int flag = KINSetEtaParams(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vegamma), Double_val(vealpha));
    CHECK_FLAG("KINSetEtaParams", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_res_mon_const_value(value vkin_mem, value vomegaconst)
{
    CAMLparam2(vkin_mem, vomegaconst);

    int flag = KINSetResMonConstValue(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vomegaconst));
    CHECK_FLAG("KINSetResMonConstValue", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_res_mon_params(value vkin_mem, value vomegamin, value vomegamax)
{
    CAMLparam3(vkin_mem, vomegamin, vomegamax);

    int flag = KINSetResMonParams(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vomegamin), Double_val(vomegamax));
    CHECK_FLAG("KINSetResMonParams", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_no_min_eps(value vkin_mem, value vnomineps)
{
    CAMLparam2(vkin_mem, vnomineps);

    int flag = KINSetNoMinEps(KINSOL_MEM_FROM_ML(vkin_mem), Bool_val(vnomineps));
    CHECK_FLAG("KINSetNoMinEps", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_max_newton_step(value vkin_mem, value vmxnewtstep)
{
    CAMLparam2(vkin_mem, vmxnewtstep);

    int flag = KINSetMaxNewtonStep(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vmxnewtstep));
    CHECK_FLAG("KINSetMaxNewtonStep", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_max_beta_fails(value vkin_mem, value vmxnbcf)
{
    CAMLparam2(vkin_mem, vmxnbcf);

    int flag = KINSetMaxBetaFails(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vmxnbcf));
    CHECK_FLAG("KINSetMaxBetaFails", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_rel_err_func(value vkin_mem, value vrelfunc)
{
    CAMLparam2(vkin_mem, vrelfunc);

    int flag = KINSetRelErrFunc(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vrelfunc));
    CHECK_FLAG("KINSetRelErrFunc", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_func_norm_tol(value vkin_mem, value vfnormtol)
{
    CAMLparam2(vkin_mem, vfnormtol);

    int flag = KINSetFuncNormTol(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vfnormtol));
    CHECK_FLAG("KINSetFuncNormTol", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_scaled_step_tol(value vkin_mem, value vscsteptol)
{
    CAMLparam2(vkin_mem, vscsteptol);

    int flag = KINSetScaledStepTol(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vscsteptol));
    CHECK_FLAG("KINSetScaledStepTol", flag);

    CAMLreturn0;
}

CAMLprim value c_kinsol_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = KINGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrw, &leniw);
    CHECK_FLAG("KINGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_kinsol_get_num_func_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumFuncEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_num_nonlin_solv_iters(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumNonlinSolvIters(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_num_beta_cond_fails(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumBetaCondFails(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumBetaCondFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_num_backtrack_ops(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumBacktrackOps(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumBacktrackOps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_func_norm(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    realtype v;

    flag = KINGetFuncNorm(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetFuncNorm", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_kinsol_get_step_length(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    realtype v;

    flag = KINGetStepLength(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetStepLength", flag);

    CAMLreturn(caml_copy_double(v));
}

