/* Aug 2010, Timothy Bourke (INRIA)
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

#include "cvode_serial.h"

CAMLprim value c_last_step_size(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    realtype hlast;

    flag = CVodeGetLastStep(cvode_mem, &hlast);
    CHECK_FLAG("CVodeGetLastStep", flag);

    CAMLreturn(caml_copy_double(hlast));
}

CAMLprim value c_next_step_size(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    realtype hcur;

    flag = CVodeGetCurrentStep(cvode_mem, &hcur);
    CHECK_FLAG("CVodeGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(hcur));
}

CAMLprim value c_get_num_steps(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    long int v;

    flag = CVodeGetNumSteps(cvode_mem, &v);
    CHECK_FLAG("CVodeGetNumSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    long int v;

    flag = CVodeGetNumRhsEvals(cvode_mem, &v);
    CHECK_FLAG("CVodeGetNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_get_num_lin_solv_setups(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    long int v;

    flag = CVodeGetNumLinSolvSetups(cvode_mem, &v);
    CHECK_FLAG("CVodeGetNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_get_num_err_test_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    long int v;

    flag = CVodeGetNumErrTestFails(cvode_mem, &v);
    CHECK_FLAG("CVodeGetNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_get_last_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    int v;

    flag = CVodeGetLastOrder(cvode_mem, &v);
    CHECK_FLAG("CVodeGetLastOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_get_current_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    int v;

    flag = CVodeGetCurrentOrder(cvode_mem, &v);
    CHECK_FLAG("CVodeGetCurrentOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_get_actual_init_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    realtype v;

    flag = CVodeGetActualInitStep(cvode_mem, &v);
    CHECK_FLAG("CVodeGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_get_last_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    realtype v;

    flag = CVodeGetLastStep(cvode_mem, &v);
    CHECK_FLAG("CVodeGetLastStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_get_current_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    realtype v;

    flag = CVodeGetCurrentStep(cvode_mem, &v);
    CHECK_FLAG("CVodeGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_get_current_time(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    realtype v;

    flag = CVodeGetCurrentTime(cvode_mem, &v);
    CHECK_FLAG("CVodeGetCurrentTime", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_set_max_ord(value vcvode_mem, value maxord)
{
    CAMLparam2(vcvode_mem, maxord);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMaxOrd(cvode_mem, Int_val(maxord));
    CHECK_FLAG("CVodeSetMaxOrd", flag);

    CAMLreturn0;
}

CAMLprim value c_set_max_num_steps(value vcvode_mem, value mxsteps)
{
    CAMLparam2(vcvode_mem, mxsteps);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMaxNumSteps(cvode_mem, Long_val(mxsteps));
    CHECK_FLAG("CVodeSetMaxNumSteps", flag);

    CAMLreturn0;
}

CAMLprim value c_set_max_hnil_warns(value vcvode_mem, value mxhnil)
{
    CAMLparam2(vcvode_mem, mxhnil);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMaxHnilWarns(cvode_mem, Int_val(mxhnil));
    CHECK_FLAG("CVodeSetMaxHnilWarns", flag);

    CAMLreturn0;
}

CAMLprim value c_set_stab_lim_det(value vcvode_mem, value stldet)
{
    CAMLparam2(vcvode_mem, stldet);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetStabLimDet(cvode_mem, Bool_val(stldet));
    CHECK_FLAG("CVodeSetStabLimDet", flag);

    CAMLreturn0;
}

CAMLprim value c_set_init_step(value vcvode_mem, value hin)
{
    CAMLparam2(vcvode_mem, hin);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetInitStep(cvode_mem, Double_val(hin));
    CHECK_FLAG("CVodeSetInitStep", flag);

    CAMLreturn0;
}

CAMLprim value c_set_min_step(value vcvode_mem, value hmin)
{
    CAMLparam2(vcvode_mem, hmin);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMinStep(cvode_mem, Double_val(hmin));
    CHECK_FLAG("CVodeSetMinStep", flag);

    CAMLreturn0;
}

CAMLprim value c_set_max_step(value vcvode_mem, value hmax)
{
    CAMLparam2(vcvode_mem, hmax);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMaxStep(cvode_mem, Double_val(hmax));
    CHECK_FLAG("CVodeSetMaxStep", flag);

    CAMLreturn0;
}

CAMLprim value c_set_stop_time(value vcvode_mem, value tstop)
{
    CAMLparam2(vcvode_mem, tstop);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetStopTime(cvode_mem, Double_val(tstop));
    CHECK_FLAG("CVodeSetStopTime", flag);

    CAMLreturn0;
}

CAMLprim value c_set_max_err_test_fails(value vcvode_mem, value maxnef)
{
    CAMLparam2(vcvode_mem, maxnef);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMaxErrTestFails(cvode_mem, Int_val(maxnef));
    CHECK_FLAG("CVodeSetMaxErrTestFails", flag);

    CAMLreturn0;
}

CAMLprim value c_set_max_nonlin_iters(value vcvode_mem, value maxcor)
{
    CAMLparam2(vcvode_mem, maxcor);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMaxNonlinIters(cvode_mem, Int_val(maxcor));
    CHECK_FLAG("CVodeSetMaxNonlinIters", flag);

    CAMLreturn0;
}

CAMLprim value c_set_max_conv_fails(value vcvode_mem, value maxncf)
{
    CAMLparam2(vcvode_mem, maxncf);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetMaxConvFails(cvode_mem, Int_val(maxncf));
    CHECK_FLAG("CVodeSetMaxConvFails", flag);

    CAMLreturn0;
}

CAMLprim value c_set_nonlin_conv_coef(value vcvode_mem, value nlscoef)
{
    CAMLparam2(vcvode_mem, nlscoef);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetNonlinConvCoef(cvode_mem, Double_val(nlscoef));
    CHECK_FLAG("CVodeSetNonlinConvCoef", flag);

    CAMLreturn0;
}

CAMLprim value c_set_no_inactive_root_warn(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVodeSetNoInactiveRootWarn(cvode_mem);
    CHECK_FLAG("CVodeSetNoInactiveRootWarn", flag);

    CAMLreturn0;
}

CAMLprim value c_set_gs_type(value vcvode_mem, value vgstype)
{
    CAMLparam2(vcvode_mem, vgstype);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int gstype;
    switch (Int_val(vgstype)) {
    case VARIANT_GRAMSCHMIDT_TYPE_MODIFIEDGS:
	gstype = MODIFIED_GS;
	break;

    case VARIANT_GRAMSCHMIDT_TYPE_CLASSICALGS:
	gstype = CLASSICAL_GS;
	break;
    }

    int flag = CVSpilsSetGSType(cvode_mem, gstype);
    CHECK_FLAG("CVSpilsSetGSType", flag);

    CAMLreturn0;
}

CAMLprim value c_set_eps_lin(value vcvode_mem, value eplifac)
{
    CAMLparam2(vcvode_mem, eplifac);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVSpilsSetEpsLin(cvode_mem, Double_val(eplifac));
    CHECK_FLAG("CVSpilsSetEpsLin", flag);

    CAMLreturn0;
}

CAMLprim value c_set_maxl(value vcvode_mem, value maxl)
{
    CAMLparam2(vcvode_mem, maxl);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVSpilsSetMaxl(cvode_mem, Int_val(maxl));
    CHECK_FLAG("CVSpilsSetMaxl", flag);

    CAMLreturn0;
}

CAMLprim value c_densematrix_get(value vmatrix, value vij)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Densematrix.get: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Densematrix.get: invalid j");
#endif

    realtype v = DENSE_ELEM(m, i, j);
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_densematrix_set(value vmatrix, value vij, value v)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Densematrix.set: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Densematrix.set: invalid j");
#endif

    DENSE_ELEM(m, i, j) = Double_val(v);
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_bandmatrix_get(value vmatrix, value vij)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Bandmatrix.get: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Bandmatrix.get: invalid j");
#endif

    realtype v = BAND_ELEM(m, i, j);
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_bandmatrix_set(value vmatrix, value vij, value v)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Bandmatrix.set: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Bandmatrix.set: invalid j");
#endif

    BAND_ELEM(m, i, j) = Double_val(v);
    CAMLreturn(caml_copy_double(v));
}


CAMLprim value c_get_num_stab_lim_order_reds(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVodeGetNumStabLimOrderReds(cvode_mem, &r);
    CHECK_FLAG("CVodeGetNumStabLimOrderReds", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_get_tol_scale_factor(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    realtype r;
    int flag = CVodeGetTolScaleFactor(cvode_mem, &r);
    CHECK_FLAG("CVodeGetTolScaleFactor", flag);

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value c_get_num_nonlin_solv_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvIters(cvode_mem, &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_get_num_nonlin_solv_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_get_num_g_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVodeGetNumGEvals(cvode_mem, &r);
    CHECK_FLAG("CVodeGetNumGEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_dls_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVDlsGetNumJacEvals(cvode_mem, &r);
    CHECK_FLAG("CVDlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_dls_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVDlsGetNumRhsEvals(cvode_mem, &r);
    CHECK_FLAG("CVDlsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_diag_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVDiagGetNumRhsEvals(cvode_mem, &r);
    CHECK_FLAG("CVDiagGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_bandprec_get_num_rhs_evals (value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVBandPrecGetNumRhsEvals(cvode_mem, &r);
    CHECK_FLAG("CVBandPrecGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}


CAMLprim value c_spils_get_num_lin_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumLinIters(cvode_mem, &r);
    CHECK_FLAG("CVSpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_spils_get_num_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumConvFails(cvode_mem, &r);
    CHECK_FLAG("CVSpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_spils_get_num_prec_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecEvals(cvode_mem, &r);
    CHECK_FLAG("CVSpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_spils_get_num_prec_solves(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecSolves(cvode_mem, &r);
    CHECK_FLAG("CVSpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_spils_get_num_jtimes_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumJtimesEvals(cvode_mem, &r);
    CHECK_FLAG("CVSpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_spils_get_num_rhs_evals (value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumRhsEvals(cvode_mem, &r);
    CHECK_FLAG("CVSpilsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

