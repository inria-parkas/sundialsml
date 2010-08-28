/* Aug 2010, Timothy Bourke (INRIA)
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

#include "cvode_serial.h"

CAMLprim value c_last_step_size(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);
    int flag;
    realtype hlast;

    flag = CVodeGetLastStep(cvode_mem, &hlast);
    ml_cvode_check_flag("CVodeGetLastStep", flag, NULL);

    CAMLreturn(caml_copy_double(hlast));
}

CAMLprim value c_next_step_size(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);
    int flag;
    realtype hcur;

    flag = CVodeGetCurrentStep(cvode_mem, &hcur);
    ml_cvode_check_flag("CVodeGetCurrentStep", flag, NULL);

    CAMLreturn(caml_copy_double(hcur));
}

CAMLprim value c_set_max_ord(value vcvode_mem, value maxord)
{
    CAMLparam2(vcvode_mem, maxord);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMaxOrd(cvode_mem, Int_val(maxord));
    ml_cvode_check_flag("CVodeSetMaxOrd", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_num_steps(value vcvode_mem, value mxsteps)
{
    CAMLparam2(vcvode_mem, mxsteps);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMaxNumSteps(cvode_mem, Long_val(mxsteps));
    ml_cvode_check_flag("CVodeSetMaxNumSteps", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_hnil_warns(value vcvode_mem, value mxhnil)
{
    CAMLparam2(vcvode_mem, mxhnil);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMaxHnilWarns(cvode_mem, Int_val(mxhnil));
    ml_cvode_check_flag("CVodeSetMaxHnilWarns", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_stability_limit_detection(value vcvode_mem, value stldet)
{
    CAMLparam2(vcvode_mem, stldet);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetStabLimDet(cvode_mem, Bool_val(stldet));
    ml_cvode_check_flag("CVodeSetStabLimDet", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_initial_step_size(value vcvode_mem, value hin)
{
    CAMLparam2(vcvode_mem, hin);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetInitStep(cvode_mem, Double_val(hin));
    ml_cvode_check_flag("CVodeSetInitStep", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_min_abs_step_size(value vcvode_mem, value hmin)
{
    CAMLparam2(vcvode_mem, hmin);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMinStep(cvode_mem, Double_val(hmin));
    ml_cvode_check_flag("CVodeSetMinStep", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_abs_step_size(value vcvode_mem, value hmax)
{
    CAMLparam2(vcvode_mem, hmax);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMaxStep(cvode_mem, Double_val(hmax));
    ml_cvode_check_flag("CVodeSetMaxStep", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_stop_time(value vcvode_mem, value tstop)
{
    CAMLparam2(vcvode_mem, tstop);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetStopTime(cvode_mem, Double_val(tstop));
    ml_cvode_check_flag("CVodeSetStopTime", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_error_test_failures(value vcvode_mem, value maxnef)
{
    CAMLparam2(vcvode_mem, maxnef);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMaxErrTestFails(cvode_mem, Int_val(maxnef));
    ml_cvode_check_flag("CVodeSetMaxErrTestFails", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_nonlinear_iterations(value vcvode_mem, value maxcor)
{
    CAMLparam2(vcvode_mem, maxcor);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMaxNonlinIters(cvode_mem, Int_val(maxcor));
    ml_cvode_check_flag("CVodeSetMaxNonlinIters", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_convergence_failures(value vcvode_mem, value maxncf)
{
    CAMLparam2(vcvode_mem, maxncf);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetMaxConvFails(cvode_mem, Int_val(maxncf));
    ml_cvode_check_flag("CVodeSetMaxConvFails", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_nonlinear_convergence_coeffficient(value vcvode_mem, value nlscoef)
{
    CAMLparam2(vcvode_mem, nlscoef);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetNonlinConvCoef(cvode_mem, Double_val(nlscoef));
    ml_cvode_check_flag("CVodeSetNonlinConvCoef", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_disable_inactive_root_warnings(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVodeSetNoInactiveRootWarn(cvode_mem);
    ml_cvode_check_flag("CVodeSetNoInactiveRootWarn", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_preconditioning_type(value vcvode_mem, value vptype)
{
    CAMLparam2(vcvode_mem, vptype);
    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int ptype;
    switch (Int_val(vptype)) {
    case VARIANT_PRECONDITIONING_TYPE_PRECNONE:
	ptype = PREC_NONE;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECLEFT:
	ptype = PREC_LEFT;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECRIGHT:
	ptype = PREC_RIGHT;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECBOTH:
	ptype = PREC_BOTH;
	break;
    }

    int flag = CVSpilsSetPrecType(cvode_mem, ptype);
    ml_cvode_check_flag("CVSpilsSetPrecType", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_gramschmidt_orthogonalization(value vcvode_mem, value vgstype)
{
    CAMLparam2(vcvode_mem, vgstype);
    void *cvode_mem = ml_cvode_mem(vcvode_mem);

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
    ml_cvode_check_flag("CVSpilsSetGSType", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_eps_linear_convergence_factor(value vcvode_mem, value eplifac)
{
    CAMLparam2(vcvode_mem, eplifac);
    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVSpilsSetEpsLin(cvode_mem, Double_val(eplifac));
    ml_cvode_check_flag("CVSpilsSetEpsLin", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_set_max_subspace_dimension(value vcvode_mem, value maxl)
{
    CAMLparam2(vcvode_mem, maxl);
    void *cvode_mem = ml_cvode_mem(vcvode_mem);

    int flag = CVSpilsSetMaxl(cvode_mem, Int_val(maxl));
    ml_cvode_check_flag("CVSpilsSetMaxl", flag, NULL);

    CAMLreturn0;
}

CAMLprim value c_densematrix_get(value vmatrix, value vij)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);
    realtype v = DENSE_ELEM(m, Int_val(Field(vij, 0)), Int_val(Field(vij, 1)));
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_densematrix_set(value vmatrix, value vij, value v)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);
    DENSE_ELEM(m, Int_val(Field(vij, 0)), Int_val(Field(vij, 1)))
	= Double_val(v);
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_bandmatrix_get(value vmatrix, value vij)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);
    realtype v = BAND_ELEM(m, Int_val(Field(vij, 0)), Int_val(Field(vij, 1)));
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_bandmatrix_set(value vmatrix, value vij, value v)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = (DlsMat)Field(vmatrix, 0);
    BAND_ELEM(m, Int_val(Field(vij, 0)), Int_val(Field(vij, 1)))
	= Double_val(v);
    CAMLreturn(caml_copy_double(v));
}

