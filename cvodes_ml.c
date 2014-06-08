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

// TODO: remove unneeded header files, add needed ones

#include <cvodes/cvodes.h>
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
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_direct.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_spils.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <cvodes/cvodes_lapack.h>
#endif

#include "spils_ml.h"
#include "cvode_ml.h"
#include "cvodes_ml.h"
#include "sundials_ml.h"

#define MAX_ERRMSG_LEN 256

// TODO: must ensure that these exceptions are registered...
void cvodes_ml_SCHECK_FLAG(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	    || flag == CV_ROOT_RETURN
	    || flag == CV_TSTOP_RETURN) return;

    switch (flag) {
	case CV_TOO_MUCH_WORK:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchWork"));

	case CV_TOO_MUCH_ACC:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchAccuracy"));

	case CV_ERR_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ErrFailure"));

	case CV_CONV_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ConvergenceFailure"));

	/* * */

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

	/* * */

	case CV_ILL_INPUT:
	    caml_raise_constant(*caml_named_value("cvode_IllInput"));
	
	case CV_BAD_K:
	    caml_raise_constant(*caml_named_value("cvode_BadK"));

	case CV_BAD_T:
	    caml_raise_constant(*caml_named_value("cvode_BadT"));

	case CV_BAD_DKY:
	    caml_raise_constant(*caml_named_value("cvode_BadDky"));

	case CV_TOO_CLOSE:
	    caml_raise_constant(*caml_named_value("cvode_TooClose"));

	/* Quadrature */

        case CV_NO_QUAD:
	    caml_raise_constant(*caml_named_value("cvodes_QuadNotInitialized"));

        case CV_QRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_QuadRhsFuncFailure"));

        case CV_FIRST_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstQuadRhsFuncErr"));

        case CV_REPTD_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedQuadRhsFuncErr"));

        case CV_UNREC_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableQuadRhsFuncErr"));

	/* Sensitivity */

        case CV_NO_SENS:
	    caml_raise_constant(*caml_named_value("cvodes_SensNotInitialized"));

        case CV_SRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_SensRhsFuncFailure"));

        case CV_FIRST_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstSensRhsFuncErr"));

        case CV_REPTD_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedSensRhsFuncErr"));

        case CV_UNREC_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableSensRhsFuncErr"));

        case CV_BAD_IS:
	    caml_raise_constant(*caml_named_value("cvodes_BadIS"));

	/* Sensitivity > Quadrature */

        case CV_NO_QUADSENS:
	    caml_raise_constant(*caml_named_value("cvodes_QuadSensNotInitialized"));

        case CV_QSRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_QuadSensRhsFuncFailure"));

        case CV_FIRST_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstQuadSensRhsFuncErr"));

        case CV_REPTD_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedQuadSensRhsFuncErr"));

        case CV_UNREC_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableQuadSensRhsFuncErr"));

	/* Adjoint */

        case CV_NO_ADJ:
	    caml_raise_constant(*caml_named_value("cvodes_AdjointNotInitialized"));

        case CV_NO_FWD:
	    caml_raise_constant(*caml_named_value("cvodes_NoForwardCall"));

        case CV_NO_BCK:
	    caml_raise_constant(*caml_named_value("cvodes_NoBackwardProblem"));

        case CV_BAD_TB0:
	    caml_raise_constant(*caml_named_value("cvodes_BadFinalTime"));

        case CV_REIFWD_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_ForwardReinitializationFailed"));

        case CV_FWD_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_ForwardFailed"));

        case CV_GETY_BADT:
	    caml_raise_constant(*caml_named_value("cvodes_BadOutputTime"));

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

/* basic quadrature interface */

CAMLprim void c_cvodes_quad_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = CVodeSetQuadErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadErrCon", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quad_ss_tolerances(value vdata,
					  value reltol,
					  value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeQuadSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    SCHECK_FLAG("CVodeQuadSStolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_quad_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quad_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quad_get_stats(value vdata)
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

CAMLprim void c_cvodes_sens_set_err_con(value vdata, value verrcons)
{
    CAMLparam2(vdata, verrcons);
    int flag;
    
    flag = CVodeSetSensErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrcons));
    SCHECK_FLAG("CVodeSetSensErrCon", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_ss_tolerances(value vdata,
					  value reltol,
					  value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeSensSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), REAL_ARRAY(abstol));
    SCHECK_FLAG("CVodeSensSStolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSensEEtolerances(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSensEEtolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_set_params(value vdata, value vparams)
{
    CAMLparam2(vdata, vparams);
    CAMLlocal3(vp, vpbar, vplist);

    realtype *p = NULL;
    realtype *pbar = NULL;
    int *plist = NULL;
    int i, ns;

    vp = Field(vparams, RECORD_CVODES_SENS_PARAMS_PVALS);
    vpbar = Field(vparams, RECORD_CVODES_SENS_PARAMS_PBAR);
    vplist = Field(vparams, RECORD_CVODES_SENS_PARAMS_PLIST);

    if (vp != Val_none) p = REAL_ARRAY(Some_val(vp));
    if (vpbar != Val_none) pbar = REAL_ARRAY(Some_val(vpbar));

    if (vplist != Val_none) {
	vplist = Some_val(vplist);
	ns = (int)caml_array_length(vplist);
	plist = calloc(ns, sizeof(int));

	for (i=0; i < ns; ++i) {
	    plist[i] = Int_val(Field(vplist, i));
	}
    }

    int flag = CVodeSetSensParams(CVODE_MEM_FROM_ML(vdata), p, pbar, plist);
    if (plist != NULL) free(plist);
    SCHECK_FLAG("CVodeSetSensParams", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_toggle_off(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSensToggleOff(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSensToggleOff", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_set_dq_method(value vdata, value vdqtype,
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

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_set_max_nonlin_iters(value vdata, value vmaxcors)
{
    CAMLparam2(vdata, vmaxcors);

    int flag = CVodeSetSensMaxNonlinIters(CVODE_MEM_FROM_ML(vdata),
	    Int_val(vmaxcors));
    SCHECK_FLAG("CVodeSetSensMaxNonlinIters", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_sens_get_num_sens_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetNumRhsEvalsSens(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetNumRhsEvalsSens", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_lin_solv_setups(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumLinSolvSetups(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_stats(value vdata)
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
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_RHS_EVALS, Val_long(nfsevals));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_SENS_EVALS, Val_long(nfevalss));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_ERR_TEST_FAILS,
							   Val_long(nsetfails));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_LIN_SOLV_SETUPS,
							 Val_long(nlinsetupss));

    CAMLreturn(r);
}

CAMLprim value c_cvodes_sens_get_num_nonlin_solv_iters(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumNonlinSolvIters(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_nonlin_solv_conv_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_nonlin_solv_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nsniters;
    long int nsncfails;

    flag = CVodeGetSensNonlinSolvStats(CVODE_MEM_FROM_ML(vdata), &nsniters,
	    &nsncfails);
    SCHECK_FLAG("CVodeGetSensNonlinSolvStats", flag);

    r = caml_alloc_tuple(RECORD_CVODES_SENS_NONLIN_STATS_SIZE);
    Store_field(r, RECORD_CVODES_SENS_NONLIN_STATS_NUM_SOLV_ITERS,
							    Val_long(nsniters));
    Store_field(r, RECORD_CVODES_SENS_NONLIN_STATS_NUM_CONV_FAILS,
							   Val_long(nsncfails));

    CAMLreturn(r);
}

CAMLprim void c_cvodes_sens_get_stgr_nonlin_solv_iters(value vdata, value vr)
{
    CAMLparam2(vdata, vr);

    int flag;

    flag = CVodeGetStgrSensNumNonlinSolvIters(CVODE_MEM_FROM_ML(vdata),
					      LONG_ARRAY(vr));
    SCHECK_FLAG("CVodeGetStgrSensNumNonlinSolvIters", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_get_num_stgr_nonlin_solv_conv_fails(value vdata,
							        value vr)
{
    CAMLparam2(vdata, vr);

    int flag = CVodeGetStgrSensNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vdata),
						      LONG_ARRAY(vr));
    SCHECK_FLAG("CVodeGetStgrSensNumNonlinSolvConvFails", flag);

    CAMLreturn0;
}

/* sensitivity/quadrature interface */

CAMLprim void c_cvodes_quadsens_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = CVodeSetQuadSensErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadSensErrCon", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_ss_tolerances(value vdata,
					      value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeQuadSensSStolerances(CVODE_MEM_FROM_ML(vdata),
					 Double_val(reltol),
					 REAL_ARRAY(abstol));
    SCHECK_FLAG("CVodeQuadSensSStolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeQuadSensEEtolerances(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeQuadSensEEtolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_quadsens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadSensNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quadsens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadSensNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quadsens_get_stats(value vdata)
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

CAMLprim void c_cvodes_adj_init(value vdata, value vnd, value vinterptype)
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

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_ss_tolerances(value vparent, value vwhich,
					 value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);

    int flag = CVodeSStolerancesB(CVODE_MEM_FROM_ML(vparent),
		 Int_val(vwhich), Double_val(vreltol), Double_val(vabstol));
    SCHECK_FLAG("CVodeSStolerancesB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_diag(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);

    int flag = CVDiagB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich));
    SCHECK_FLAG("CVDiagB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_spgmr (value vparent, value vwhich,
				        value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpgmrB (cvode_mem, which, spils_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpgmrB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_spbcg (value vparent, value vwhich,
				        value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpbcgB (cvode_mem, which, spils_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpbcgB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_sptfqmr (value vparent, value vwhich,
				          value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSptfqmrB (cvode_mem, which, spils_precond_type (vtype),
		       Int_val (vmaxl));
    SCHECK_FLAG ("CVSptfqmrB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_set_functional (value vparent, value vwhich)
{
    CAMLparam2 (vparent, vwhich);
    int flag = CVodeSetIterTypeB (CVODE_MEM_FROM_ML (vparent), Int_val(vwhich),
				  CV_FUNCTIONAL);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_bsession_finalize(value vdata)
{
    if (CVODE_MEM_FROM_ML(vdata) != NULL) {
	value *backref = CVODE_BACKREF_FROM_ML(vdata);
	// NB: CVodeFree() is *not* called: parents free-up backward problems
	caml_remove_generational_global_root (backref);
	free (backref);
    }
}

CAMLprim void c_cvodes_adj_backward_normal(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = CVodeB(CVODE_MEM_FROM_ML(vdata), Double_val(vtbout), CV_NORMAL);
    SCHECK_FLAG("CVodeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_backward_one_step(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = CVodeB(CVODE_MEM_FROM_ML(vdata), Double_val(vtbout),
		      CV_ONE_STEP);
    SCHECK_FLAG("CVodeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_no_sensitivity(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSetAdjNoSensi(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSetAdjNoSensi", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_max_ord(value vparent, value vwhich,
				       value vmaxord)
{
    CAMLparam3(vparent, vwhich, vmaxord);

    int flag = CVodeSetMaxOrdB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(vmaxord));
    SCHECK_FLAG("CVodeSetMaxOrdB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_max_num_steps(value vparent, value vwhich,
					     value vmxsteps)
{
    CAMLparam3(vparent, vwhich, vmxsteps);

    int flag = CVodeSetMaxNumStepsB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				    Long_val(vmxsteps));
    SCHECK_FLAG("CVodeSetMaxNumStepsB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_init_step(value vparent, value vwhich, value vhin)
{
    CAMLparam3(vparent, vwhich, vhin);

    int flag = CVodeSetInitStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			         Double_val(vhin));
    SCHECK_FLAG("CVodeSetInitStepB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_min_step(value vparent, value vwhich, value vhmin)
{
    CAMLparam3(vparent, vwhich, vhmin);

    int flag = CVodeSetMinStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			        Double_val(vhmin));
    SCHECK_FLAG("CVodeSetMinStepB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_max_step(value vparent, value vwhich, value vhmax)
{
    CAMLparam3(vparent, vwhich, vhmax);

    int flag = CVodeSetMaxStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			        Double_val(vhmax));
    SCHECK_FLAG("CVodeSetMaxStepB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_stab_lim_det(value vparent, value vwhich,
					    value vstldet)
{
    CAMLparam3(vparent, vwhich, vstldet);

    int flag = CVodeSetStabLimDetB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			           Bool_val(vstldet));
    SCHECK_FLAG("CVodeSetStabLimDetB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_prec_type(value vparent, value vwhich,
					 value vptype)
{
    CAMLparam3(vparent, vwhich, vptype);

    int flag = CVSpilsSetPrecTypeB(CVODE_MEM_FROM_ML(vparent),
				   Int_val(vwhich), spils_precond_type(vptype));
    SCHECK_FLAG("CVSpilsSetPrecTypeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_gs_type(value vparent, value vwhich,
				       value vgstype)
{
    CAMLparam3(vparent, vwhich, vgstype);

    int flag = CVSpilsSetGSTypeB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				 spils_gs_type(vgstype));
    SCHECK_FLAG("CVSpilsSetGSTypeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_eps_lin(value vparent, value vwhich,
				       value eplifac)
{
    CAMLparam3(vparent, vwhich, eplifac);

    int flag = CVSpilsSetEpslinB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				 Double_val(eplifac));
    SCHECK_FLAG("CVSpilsSetEpslinB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_maxl(value vparent, value vwhich, value maxl)
{
    CAMLparam3(vparent, vwhich, maxl);

    int flag = CVSpilsSetMaxlB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(maxl));
    SCHECK_FLAG("CVSpilsSetMaxlB", flag);

    CAMLreturn0;
}

/* adjoint/quadrature interface */

CAMLprim void c_cvodes_adjquad_set_err_con(value vparent, value vwhich,
					   value verrconq)
{
    CAMLparam3(vparent, vwhich, verrconq);
    int flag;
    
    flag = CVodeSetQuadErrConB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadErrConB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adjquad_ss_tolerances(value vparent, value vwhich,
					     value reltol, value abstol)
{
    CAMLparam4(vparent, vwhich, reltol, abstol);

    int flag = CVodeQuadSStolerancesB(CVODE_MEM_FROM_ML(vparent),
				      Int_val(vwhich), Double_val(reltol),
				      Double_val(abstol));
    SCHECK_FLAG("CVodeQuadSStolerancesB", flag);

    CAMLreturn0;
}

