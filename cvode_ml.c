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

/* Interface functions that do not involve NVectors. */

#include <cvode/cvode.h>
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
#include <cvode/cvode_dense.h>
#include <cvode/cvode_direct.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#include <cvode/cvode_bandpre.h>
#include <cvode/cvode_spils.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <cvode/cvode_lapack.h>
#endif

#include "spils_ml.h"
#include "cvode_ml.h"
#include "sundials_ml.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

void cvode_ml_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	    || flag == CV_ROOT_RETURN
	    || flag == CV_TSTOP_RETURN) return;

    switch (flag) {
	case CV_ILL_INPUT:
	    caml_raise_constant(*caml_named_value("cvode_IllInput"));

	case CV_TOO_CLOSE:
	    caml_raise_constant(*caml_named_value("cvode_TooClose"));

	case CV_TOO_MUCH_WORK:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchWork"));

	case CV_TOO_MUCH_ACC:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchAccuracy"));

	case CV_ERR_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ErrFailure"));

	case CV_CONV_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ConvergenceFailure"));

	case CV_LINIT_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_LinearInitFailure"));

	case CV_LSETUP_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_LinearSetupFailure"));

	case CV_LSOLVE_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_LinearSolveFailure"));

	case CV_RHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_RhsFuncFailure"));

	case CV_FIRST_RHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvode_FirstRhsFuncErr"));

	case CV_REPTD_RHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvode_RepeatedRhsFuncErr"));

	case CV_UNREC_RHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvode_UnrecoverableRhsFuncErr"));

	case CV_RTFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_RootFuncFailure"));

	case CV_BAD_K:
	    caml_raise_constant(*caml_named_value("cvode_BadK"));

	case CV_BAD_T:
	    caml_raise_constant(*caml_named_value("cvode_BadT"));

	case CV_BAD_DKY:
	    caml_raise_constant(*caml_named_value("cvode_BadDky"));

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

/* basic interface */

CAMLprim void c_cvode_session_finalize(value vdata)
{
    if (CVODE_MEM_FROM_ML(vdata) != NULL) {
	void *cvode_mem = CVODE_MEM_FROM_ML(vdata);
	value *backref = CVODE_BACKREF_FROM_ML(vdata);
	CVodeFree(&cvode_mem);
	caml_remove_generational_global_root (backref);
	free (backref);
    }

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_CVODE_SESSION_ERRFILE));
    if (err_file != NULL) {
	fclose(err_file);
    }
}

CAMLprim void c_cvode_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("CVodeSStolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < CVODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = CVodeGetRootInfo(CVODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("CVodeGetRootInfo", flag);

    CAMLreturn0;
}

CAMLprim value c_cvode_get_integrator_stats(value vdata)
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

CAMLprim void c_cvode_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_CVODE_SESSION_ERRFILE));

    if (err_file != NULL) {
	fclose(err_file);
	Store_field(vdata, RECORD_CVODE_SESSION_ERRFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    err_file = fopen(String_val(vpath), mode);
    if (err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = CVodeSetErrFile(CVODE_MEM_FROM_ML(vdata), err_file);
    CHECK_FLAG("CVodeSetErrFile", flag);

    Store_field(vdata, RECORD_CVODE_SESSION_ERRFILE, Val_long(err_file));

    CAMLreturn0;
}

CAMLprim void c_cvode_set_functional (value vdata)
{
    CAMLparam1 (vdata);
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vdata), CV_FUNCTIONAL);
    CHECK_FLAG ("CVodeSetIterType", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < CVODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = CVodeSetRootDirection(CVODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("CVodeSetRootDirection", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_prec_type(value vcvode_mem, value vptype)
{
    CAMLparam2(vcvode_mem, vptype);

    int flag = CVSpilsSetPrecType(CVODE_MEM_FROM_ML(vcvode_mem),
				  spils_precond_type(vptype));
    CHECK_FLAG("CVSpilsSetPrecType", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_diag (value vcvode_mem)
{
    CAMLparam1 (vcvode_mem);
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vcvode_mem), CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVDiag(CVODE_MEM_FROM_ML (vcvode_mem));
    CHECK_FLAG("CVDiag", flag);
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

CAMLprim value c_cvode_get_work_space(value vcvode_mem)
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

CAMLprim value c_cvode_get_num_steps(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumSteps(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_num_lin_solv_setups(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumLinSolvSetups(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_num_err_test_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumErrTestFails(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_last_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    int v;

    flag = CVodeGetLastOrder(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetLastOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_cvode_get_current_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    int v;

    flag = CVodeGetCurrentOrder(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_cvode_get_actual_init_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetActualInitStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_cvode_get_last_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetLastStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetLastStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_cvode_get_current_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetCurrentStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_cvode_get_current_time(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetCurrentTime(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentTime", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim void c_cvode_set_max_ord(value vcvode_mem, value maxord)
{
    CAMLparam2(vcvode_mem, maxord);

    int flag = CVodeSetMaxOrd(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxord));
    CHECK_FLAG("CVodeSetMaxOrd", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_num_steps(value vcvode_mem, value mxsteps)
{
    CAMLparam2(vcvode_mem, mxsteps);


    int flag = CVodeSetMaxNumSteps(CVODE_MEM_FROM_ML(vcvode_mem), Long_val(mxsteps));
    CHECK_FLAG("CVodeSetMaxNumSteps", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_hnil_warns(value vcvode_mem, value mxhnil)
{
    CAMLparam2(vcvode_mem, mxhnil);


    int flag = CVodeSetMaxHnilWarns(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(mxhnil));
    CHECK_FLAG("CVodeSetMaxHnilWarns", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_stab_lim_det(value vcvode_mem, value stldet)
{
    CAMLparam2(vcvode_mem, stldet);


    int flag = CVodeSetStabLimDet(CVODE_MEM_FROM_ML(vcvode_mem), Bool_val(stldet));
    CHECK_FLAG("CVodeSetStabLimDet", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_init_step(value vcvode_mem, value hin)
{
    CAMLparam2(vcvode_mem, hin);


    int flag = CVodeSetInitStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hin));
    CHECK_FLAG("CVodeSetInitStep", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_min_step(value vcvode_mem, value hmin)
{
    CAMLparam2(vcvode_mem, hmin);


    int flag = CVodeSetMinStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hmin));
    CHECK_FLAG("CVodeSetMinStep", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_step(value vcvode_mem, value hmax)
{
    CAMLparam2(vcvode_mem, hmax);


    int flag = CVodeSetMaxStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hmax));
    CHECK_FLAG("CVodeSetMaxStep", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_stop_time(value vcvode_mem, value tstop)
{
    CAMLparam2(vcvode_mem, tstop);


    int flag = CVodeSetStopTime(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(tstop));
    CHECK_FLAG("CVodeSetStopTime", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_err_test_fails(value vcvode_mem, value maxnef)
{
    CAMLparam2(vcvode_mem, maxnef);


    int flag = CVodeSetMaxErrTestFails(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxnef));
    CHECK_FLAG("CVodeSetMaxErrTestFails", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_nonlin_iters(value vcvode_mem, value maxcor)
{
    CAMLparam2(vcvode_mem, maxcor);


    int flag = CVodeSetMaxNonlinIters(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxcor));
    CHECK_FLAG("CVodeSetMaxNonlinIters", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_conv_fails(value vcvode_mem, value maxncf)
{
    CAMLparam2(vcvode_mem, maxncf);


    int flag = CVodeSetMaxConvFails(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxncf));
    CHECK_FLAG("CVodeSetMaxConvFails", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_nonlin_conv_coef(value vcvode_mem, value nlscoef)
{
    CAMLparam2(vcvode_mem, nlscoef);


    int flag = CVodeSetNonlinConvCoef(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(nlscoef));
    CHECK_FLAG("CVodeSetNonlinConvCoef", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_no_inactive_root_warn(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag = CVodeSetNoInactiveRootWarn(CVODE_MEM_FROM_ML(vcvode_mem));
    CHECK_FLAG("CVodeSetNoInactiveRootWarn", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_gs_type(value vcvode_mem, value vgstype)
{
    CAMLparam2(vcvode_mem, vgstype);

    int flag = CVSpilsSetGSType(CVODE_MEM_FROM_ML(vcvode_mem),
				spils_gs_type(vgstype));
    CHECK_FLAG("CVSpilsSetGSType", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_eps_lin(value vcvode_mem, value eplifac)
{
    CAMLparam2(vcvode_mem, eplifac);

    int flag = CVSpilsSetEpsLin(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(eplifac));
    CHECK_FLAG("CVSpilsSetEpsLin", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_maxl(value vcvode_mem, value maxl)
{
    CAMLparam2(vcvode_mem, maxl);

    int flag = CVSpilsSetMaxl(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxl));
    CHECK_FLAG("CVSpilsSetMaxl", flag);

    CAMLreturn0;
}

/* statistic accessor functions */

CAMLprim value c_cvode_get_num_stab_lim_order_reds(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumStabLimOrderReds(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumStabLimOrderReds", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_get_tol_scale_factor(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    realtype r;
    int flag = CVodeGetTolScaleFactor(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetTolScaleFactor", flag);

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value c_cvode_get_num_nonlin_solv_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_get_num_nonlin_solv_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_get_num_g_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumGEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumGEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_dls_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = CVDlsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrwLS, &leniwLS);
    CHECK_FLAG("CVDlsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value c_cvode_dls_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDlsGetNumJacEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVDlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_dls_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDlsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVDlsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_diag_get_work_space(value vcvode_mem)
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

CAMLprim value c_cvode_diag_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDiagGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVDiagGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_bandprec_get_work_space(value vcvode_mem)
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

CAMLprim value c_cvode_bandprec_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVBandPrecGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVBandPrecGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

/* spils functions */

CAMLprim void c_cvode_spils_spgmr (value vcvode_mem, value vmaxl, value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpgmr (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpgmr", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_spils_spbcg (value vcvode_mem, value vmaxl, value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpbcg (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpbcg", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_sptfqmr (value vcvode_mem, value vmaxl,
				     value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSptfqmr (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSptfqmr", flag);

    CAMLreturn0;
}

CAMLprim value c_cvode_spils_get_num_lin_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumLinIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = CVSpilsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_FLAG("CVSpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_cvode_spils_get_num_prec_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_prec_solves(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecSolves(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_jtimes_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumJtimesEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_rhs_evals (value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

