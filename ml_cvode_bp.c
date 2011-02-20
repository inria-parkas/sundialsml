/* Aug 2010, Timothy Bourke (INRIA)
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

#include <string.h> /* memcpy */

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
#include <cvode/cvode_band.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#include <cvode/cvode_bandpre.h>
#include <cvode/cvode_spils.h>

#include "ml_cvode.h"

#define INT_ARRAY(v) ((int *)Caml_ba_data_val(v))
#define REAL_ARRAY(v) ((realtype *)Caml_ba_data_val(v))
#define REAL_ARRAY2(v) ((realtype **)Caml_ba_data_val(v))

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

CAMLprim value c_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    long int lenrw;
    long int leniw;

    flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
    CHECK_FLAG("CVodeGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
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

/* Dense matrix functions */

#define DLSMAT(v) (DlsMat)Data_custom_val(v)

static void finalize_dlsmat(value va)
{
    DestroyMat(DLSMAT(va));
}

CAMLprim value c_densematrix_new_dense_mat(value vmn)
{
    CAMLparam1(vmn);
    CAMLlocal1(r);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    DlsMat a = NewDenseMat(m, n);
    mlsize_t approx_size = m * n * sizeof(realtype);

    r = caml_alloc_final(sizeof(DlsMat), &finalize_dlsmat,
			 approx_size, approx_size * 20);

    memcpy(DLSMAT(r), a, sizeof(a));
    CAMLreturn(r);
}

CAMLprim value c_densematrix_print_mat(value va)
{
    CAMLparam1(va);
    PrintMat(DLSMAT(va));
    CAMLreturn0;
}

CAMLprim value c_densematrix_set_to_zero(value va)
{
    CAMLparam1(va);
    SetToZero(DLSMAT(va));
    CAMLreturn0;
}

CAMLprim value c_densematrix_add_identity(value va)
{
    CAMLparam1(va);
    AddIdentity(DLSMAT(va));
    CAMLreturn0;
}

CAMLprim value c_densematrix_copy(value va, value vb)
{
    CAMLparam2(va, vb);
    DenseCopy(DLSMAT(va), DLSMAT(vb));
    CAMLreturn0;
}

CAMLprim value c_densematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
    DenseScale(Double_val(vc), DLSMAT(va));
    CAMLreturn0;
}

CAMLprim value c_densematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);
    int r = DenseGETRF(DLSMAT(va), INT_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(*caml_named_value("cvode_ZeroDiagonalElement"),
			    Val_int(r));
    }
    CAMLreturn0;
}

CAMLprim value c_densematrix_getrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);
    DenseGETRS(DLSMAT(va), INT_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim value c_densematrix_potrf(value va)
{
    CAMLparam1(va);
    DensePOTRF(DLSMAT(va));
    CAMLreturn0;
}

CAMLprim value c_densematrix_potrs(value va, value vb)
{
    CAMLparam2(va, vb);
    DensePOTRS(DLSMAT(va), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim value c_densematrix_geqrf(value va, value vbeta, value vwork)
{
    CAMLparam3(va, vbeta, vwork);
    DenseGEQRF(DLSMAT(va), REAL_ARRAY(vbeta), REAL_ARRAY(vwork));
    CAMLreturn0;
}

CAMLprim value c_densematrix_ormqr(value va, value vormqr)
{
    CAMLparam2(va, vormqr);

    realtype *beta = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_BETA));
    realtype *vn   = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_VN));
    realtype *vm   = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_VM));
    realtype *work = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_WORK));

    DenseORMQR(DLSMAT(va), beta, vn, vm, work);
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

/* Direct dense matrix functions */

#define DDENSEMAT(v) (*((realtype ***)(Data_custom_val(v))))

static void finalize_direct_densemat(value va)
{
    destroyMat(DDENSEMAT(va));
}

CAMLprim value c_densematrix_direct_new_dense_mat(value vmn)
{
    CAMLparam1(vmn);
    CAMLlocal1(vr);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    realtype **a = newDenseMat(m, n);
    mlsize_t approx_size = m * n * sizeof(realtype);

    vr = caml_alloc_final(sizeof(realtype **),
			  &finalize_direct_densemat,
			  approx_size, approx_size * 20);
    realtype ***r = (realtype ***)Data_custom_val(vr);
    *r = a;

    CAMLreturn(vr);
}

CAMLprim value c_densematrix_direct_get(value va, value vij)
{
    CAMLparam2(va, vij);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

    CAMLreturn(caml_copy_double(DDENSEMAT(va)[i][j]));
}

CAMLprim value c_densematrix_direct_set(value va, value vij, value vv)
{
    CAMLparam3(va, vij, vv);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

    DDENSEMAT(va)[i][j] = Double_val(vv);

    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_copy(value va, value vb, value vmn)
{
    CAMLparam3(va, vb, vmn);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    denseCopy(DDENSEMAT(va), DDENSEMAT(vb), m, n);
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_scale(value vc, value va, value vmn)
{
    CAMLparam3(vc, va, vmn);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    denseScale(Double_val(vc), DDENSEMAT(va), m, n);
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_add_identity(value va, value vn)
{
    CAMLparam2(va, vn);
    denseAddIdentity(DDENSEMAT(va), Int_val(vn));
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_getrf(value va, value vmn, value vp)
{
    CAMLparam3(va, vmn, vp);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    int r = denseGETRF(DDENSEMAT(va), m, n, INT_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(*caml_named_value("cvode_ZeroDiagonalElement"),
			    Val_int(r));
    }
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_getrs(value va, value vn,
	value vp, value vb)
{
    CAMLparam4(va, vn, vp, vb);
    denseGETRS(DDENSEMAT(va), Int_val(vn), INT_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_potrf(value va, value vm)
{
    CAMLparam2(va, vm);
    densePOTRF(DDENSEMAT(va), Int_val(vm));
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_potrs(value va, value vm, value vb)
{
    CAMLparam3(va, vm, vb);
    densePOTRS(DDENSEMAT(va), Int_val(vm), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_geqrf(value va, value vmn,
	value vbeta, value vv)
{
    CAMLparam4(va, vmn, vbeta, vv);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    denseGEQRF(DDENSEMAT(va), m, n, REAL_ARRAY(vbeta), REAL_ARRAY(vv));
    CAMLreturn0;
}

CAMLprim value c_densematrix_direct_ormqr(value va, value vmn, value vormqr)
{
    CAMLparam3(va, vmn, vormqr);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    realtype *beta = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_BETA));
    realtype *vn   = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_VN));
    realtype *vm   = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_VM));
    realtype *work = REAL_ARRAY(Field(vormqr, RECORD_DENSEMATRIX_ORMQR_WORK));

    denseORMQR(DDENSEMAT(va), m, n, beta, vn, vm, work);
    CAMLreturn0;
}

/* Band matrix functions */

CAMLprim value c_bandmatrix_new_band_mat(value vsizes)
{
    CAMLparam1(vsizes);
    CAMLlocal1(r);

    int n   = Int_val(Field(vsizes, 0));
    int mu  = Int_val(Field(vsizes, 1));
    int ml  = Int_val(Field(vsizes, 2));
    int smu = Int_val(Field(vsizes, 3));

    DlsMat va = NewBandMat(n, mu, ml, smu);
    mlsize_t approx_size = n * (smu + ml + 2) * sizeof(realtype);

    r = caml_alloc_final(sizeof(DlsMat), &finalize_dlsmat,
			 approx_size, approx_size * 20);

    memcpy(DLSMAT(r), va, sizeof(va));
    CAMLreturn(r);
}

CAMLprim value c_bandmatrix_copy(value va, value vb,
	value vcopymu, value vcopyml)
{
    CAMLparam4(va, vb, vcopymu, vcopyml);
    BandCopy(DLSMAT(va), DLSMAT(vb), Int_val(vcopymu), Int_val(vcopyml));
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
    BandScale(Double_val(vc), DLSMAT(va));
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_gbtrf(value va, value vp)
{
    CAMLparam2(va, vp);
    BandGBTRF(DLSMAT(va), INT_ARRAY(vp));
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_gbtrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);
    BandGBTRS(DLSMAT(va), INT_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
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

CAMLprim value c_bandmatrix_col_get_col(value vmatrix, value vj)
{
    CAMLparam2(vmatrix, vj);
    CAMLlocal1(r);

    DlsMat m = (DlsMat)Field(vmatrix, 0);

    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    if (j < 0 || j >= m->N)
	caml_invalid_argument("Bandmatrix.Col.get_col: invalid j");

    r = caml_alloc(3, Abstract_tag);
    Store_field(r, 1, Val_int(m->mu));
    Store_field(r, 2, Val_int(m->ml));
#else
    r = caml_alloc(1, Abstract_tag);
#endif

    Store_field(r, 0, (value)BAND_COL(m, j));
    CAMLreturn(r);
}

CAMLprim value c_bandmatrix_col_get(value vbcol, value vi, value vj)
{
    CAMLparam3(vbcol, vi, vj);

    realtype *bcol = (realtype *)Field(vbcol, 0);

    int i = Int_val(vi);
    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    int mu = Int_val(Field(vbcol, 1));
    int ml = Int_val(Field(vbcol, 2));

    if (i < (j - mu) || i > (j + ml))
	caml_invalid_argument("Bandmatrix.Col.get: invalid i");
#endif

    CAMLreturn(caml_copy_double(BAND_COL_ELEM(bcol, i, j)));
}

CAMLprim value c_bandmatrix_col_set(value vbcol, value vi, value vj, value ve)
{
    CAMLparam4(vbcol, vi, vj, ve);

    realtype *bcol = (realtype *)Field(vbcol, 0);

    int i = Int_val(vi);
    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    int mu = Int_val(Field(vbcol, 1));
    int ml = Int_val(Field(vbcol, 2));

    if (i < (j - mu) || i > (j + ml))
	caml_invalid_argument("Bandmatrix.Col.set: invalid i");
#endif

    BAND_COL_ELEM(bcol, i, j) = Double_val(ve);

    CAMLreturn(Val_unit);
}

/* Band matrix direct functions */

#define DBANDMAT(v) (*((realtype ***)(Data_custom_val(v))))

static void finalize_direct_bandmat(value va)
{
    destroyMat(DBANDMAT(va));
}

CAMLprim value c_bandmatrix_direct_new_band_mat(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal1(vr);

    int n   = Int_val(Field(vargs, 0));
    int smu = Int_val(Field(vargs, 1));
    int ml  = Int_val(Field(vargs, 2));

    realtype **a = newBandMat(n, smu, ml);
    mlsize_t approx_size = n * (smu + ml + 2) * sizeof(realtype);

    vr = caml_alloc_final(sizeof(realtype **),
			  &finalize_direct_bandmat,
			  approx_size, approx_size * 20);
    realtype ***r = (realtype ***)Data_custom_val(vr);
    *r = a;

    CAMLreturn(vr);
}

CAMLprim value c_bandmatrix_direct_copy(value va, value vb, value vsizes)
{
    CAMLparam3(va, vb, vsizes);

    int n      = Int_val(Field(vsizes, 0));
    int a_smu  = Int_val(Field(vsizes, 1));
    int b_smu  = Int_val(Field(vsizes, 2));
    int copymu = Int_val(Field(vsizes, 3));
    int copyml = Int_val(Field(vsizes, 4));

    bandCopy(DBANDMAT(va), DBANDMAT(vb), n, a_smu, b_smu, copymu, copyml);
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_direct_scale(value vc, value va, value vsizes)
{
    CAMLparam3(vc, va, vsizes);

    int n   = Int_val(Field(vsizes, 0));
    int mu  = Int_val(Field(vsizes, 1));
    int ml  = Int_val(Field(vsizes, 2));
    int smu = Int_val(Field(vsizes, 3));

    bandScale(Double_val(vc), DBANDMAT(va), n, mu, ml, smu);
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_direct_add_identity(value va, value vn, value vsmu)
{
    CAMLparam3(va, vn, vsmu);

    bandAddIdentity(DBANDMAT(va), Int_val(vn), Int_val(vsmu));
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_direct_gbtrf(value va, value vsizes, value vp)
{
    CAMLparam3(va, vsizes, vp);

    int n   = Int_val(Field(vsizes, 0));
    int mu  = Int_val(Field(vsizes, 1));
    int ml  = Int_val(Field(vsizes, 2));
    int smu = Int_val(Field(vsizes, 3));

    bandGBTRF(DBANDMAT(va), n, mu, ml, smu, INT_ARRAY(vp));
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_direct_gbtrs(value va, value vsizes, value vp, value vb)
{
    CAMLparam4(va, vsizes, vp, vb);

    int n   = Int_val(Field(vsizes, 0));
    int smu = Int_val(Field(vsizes, 1));
    int ml  = Int_val(Field(vsizes, 2));

    bandGBTRS(DBANDMAT(va), n, smu, ml, INT_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

/* statistic accessor functions */

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

CAMLprim value c_dls_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int lenrwLS;
    long int leniwLS;

    int flag = CVDlsGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
    CHECK_FLAG("CVDlsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));

    CAMLreturn(r);
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

CAMLprim value c_diag_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int lenrwLS;
    long int leniwLS;

    int flag = CVDiagGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
    CHECK_FLAG("CVDiagGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));

    CAMLreturn(r);
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

CAMLprim value c_bandprec_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int lenrwBP;
    long int leniwBP;

    int flag = CVBandPrecGetWorkSpace(cvode_mem, &lenrwBP, &leniwBP);
    CHECK_FLAG("CVBandPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwBP));
    Store_field(r, 1, Val_int(leniwBP));

    CAMLreturn(r);
}

CAMLprim value c_bandprec_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    long int r;
    int flag = CVBandPrecGetNumRhsEvals(cvode_mem, &r);
    CHECK_FLAG("CVBandPrecGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

/* spils functions */

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

CAMLprim value c_spils_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    int flag;
    long int lenrw;
    long int leniw;

    flag = CVSpilsGetWorkSpace(cvode_mem, &lenrw, &leniw);
    CHECK_FLAG("CVSpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
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

