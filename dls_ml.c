/***********************************************************************
 *                                                                     *
 *     OCaml interface to Sundials (serial) CVODE and IDA solvers      *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) * 
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <sundials/sundials_dense.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

#include <stdio.h>

#define INT_ARRAY(v) ((int *)Caml_ba_data_val(v))
#define LONG_ARRAY(v) ((long int *)Caml_ba_data_val(v))
#define REAL_ARRAY(v) ((realtype *)Caml_ba_data_val(v))
#define REAL_ARRAY2(v) ((realtype **)Caml_ba_data_val(v))

/* Dense matrix functions */

#define DLSMAT(v) (*(DlsMat *)Data_custom_val(v))

static void finalize_dlsmat(value va)
{
    DestroyMat(DLSMAT(va));
}

CAMLprim value c_densematrix_new_dense_mat(value vm, value vn)
{
    CAMLparam2(vm, vn);
    CAMLlocal1(vr);

    int m = Int_val(vm);
    int n = Int_val(vn);

    DlsMat a = NewDenseMat(m, n);
    if (a == NULL)
	caml_failwith("Could not create Dense Matrix.");
    mlsize_t approx_size = m * n * sizeof(realtype);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vr = caml_alloc_final(2, &finalize_dlsmat, approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim void c_densematrix_print_mat(value va)
{
    CAMLparam1(va);
    PrintMat(DLSMAT(va));
    fflush(stdout);
    CAMLreturn0;
}

CAMLprim void c_densematrix_set_to_zero(value va)
{
    CAMLparam1(va);
    SetToZero(DLSMAT(va));
    CAMLreturn0;
}

CAMLprim void c_densematrix_add_identity(value va)
{
    CAMLparam1(va);
    AddIdentity(DLSMAT(va));
    CAMLreturn0;
}

CAMLprim void c_densematrix_copy(value va, value vb)
{
    CAMLparam2(va, vb);
    DenseCopy(DLSMAT(va), DLSMAT(vb));
    CAMLreturn0;
}

CAMLprim void c_densematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
    DenseScale(Double_val(vc), DLSMAT(va));
    CAMLreturn0;
}

CAMLprim void c_densematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);
    int r = DenseGETRF(DLSMAT(va), LONG_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(*caml_named_value("cvode_ZeroDiagonalElement"),
			    Val_int(r));
    }
    CAMLreturn0;
}

CAMLprim void c_densematrix_getrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);
    DenseGETRS(DLSMAT(va), LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_densematrix_potrf(value va)
{
    CAMLparam1(va);
    DensePOTRF(DLSMAT(va));
    CAMLreturn0;
}

CAMLprim void c_densematrix_potrs(value va, value vb)
{
    CAMLparam2(va, vb);
    DensePOTRS(DLSMAT(va), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_densematrix_geqrf(value va, value vbeta, value vwork)
{
    CAMLparam3(va, vbeta, vwork);
    DenseGEQRF(DLSMAT(va), REAL_ARRAY(vbeta), REAL_ARRAY(vwork));
    CAMLreturn0;
}

CAMLprim void c_densematrix_ormqr(value va, value vormqr)
{
    CAMLparam2(va, vormqr);

    realtype *beta = REAL_ARRAY(Field(vormqr, 0));
    realtype *vv   = REAL_ARRAY(Field(vormqr, 1));
    realtype *vw   = REAL_ARRAY(Field(vormqr, 2));
    realtype *work = REAL_ARRAY(Field(vormqr, 3));

    DenseORMQR(DLSMAT(va), beta, vv, vw, work);
    CAMLreturn0;
}
 
CAMLprim value c_densematrix_get(value vmatrix, value vi, value vj)
{
    CAMLparam3(vmatrix, vi, vj);
    DlsMat m = DLSMAT(vmatrix);

    int i = Int_val(vi);
    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Densematrix.get: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Densematrix.get: invalid j");
#endif

    realtype v = DENSE_ELEM(m, i, j);
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_densematrix_set(value vmatrix, value vi, value vj, value v)
{
    CAMLparam4(vmatrix, vi, vj, v);
    DlsMat m = DLSMAT(vmatrix);

    int i = Int_val(vi);
    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Densematrix.set: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Densematrix.set: invalid j");
#endif

    DENSE_ELEM(m, i, j) = Double_val(v);
    CAMLreturn(caml_copy_double(v));
}

/* Array dense matrix functions */

#define DDENSEMAT(v) (*(realtype ***)(Data_custom_val(v)))

static void finalize_arraydensematrix(value va)
{
    destroyMat(DDENSEMAT(va));
}

CAMLprim value c_arraydensematrix_new_dense_mat(value vm, value vn)
{
    CAMLparam2(vm, vn);
    CAMLlocal1(vr);

    int m = Int_val(vm);
    int n = Int_val(vn);

    realtype **a = newDenseMat(m, n);
    if (a == NULL)
	caml_failwith("Could not create Array Dense Matrix.");
    mlsize_t approx_size = m * n * sizeof(realtype);

    vr = caml_alloc_final(2, &finalize_arraydensematrix,
			  approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim value c_arraydensematrix_get(value va, value vi, value vj)
{
    CAMLparam3(va, vi, vj);

    int i = Int_val(vi);
    int j = Int_val(vj);

    CAMLreturn(caml_copy_double(DDENSEMAT(va)[j][i]));
}

CAMLprim void c_arraydensematrix_set(value va, value vi, value vj, value vv)
{
    CAMLparam4(va, vi, vj, vv);

    int i = Int_val(vi);
    int j = Int_val(vj);

    DDENSEMAT(va)[j][i] = Double_val(vv);

    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_copy(value va, value vb, value vm, value vn)
{
    CAMLparam4(va, vb, vm, vn);

    int m = Int_val(vm);
    int n = Int_val(vn);

    denseCopy(DDENSEMAT(va), DDENSEMAT(vb), m, n);
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_scale(value vc, value va, value vm, value vn)
{
    CAMLparam4(vc, va, vm, vn);

    int m = Int_val(vm);
    int n = Int_val(vn);

    denseScale(Double_val(vc), DDENSEMAT(va), m, n);
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_add_identity(value va, value vn)
{
    CAMLparam2(va, vn);
    denseAddIdentity(DDENSEMAT(va), Int_val(vn));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_getrf(value va, value vm, value vn, value vp)
{
    CAMLparam4(va, vm, vn, vp);

    int m = Int_val(vn);
    int n = Int_val(vm);

    int r = denseGETRF(DDENSEMAT(va), m, n, LONG_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(*caml_named_value("cvode_ZeroDiagonalElement"),
			    Val_int(r));
    }
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_getrs(value va, value vn,
	value vp, value vb)
{
    CAMLparam4(va, vn, vp, vb);
    denseGETRS(DDENSEMAT(va), Int_val(vn), LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_potrf(value va, value vm)
{
    CAMLparam2(va, vm);
    densePOTRF(DDENSEMAT(va), Int_val(vm));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_potrs(value va, value vm, value vb)
{
    CAMLparam3(va, vm, vb);
    densePOTRS(DDENSEMAT(va), Int_val(vm), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_geqrf(value va, value vmn,
	value vbeta, value vv)
{
    CAMLparam4(va, vmn, vbeta, vv);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    denseGEQRF(DDENSEMAT(va), m, n, REAL_ARRAY(vbeta), REAL_ARRAY(vv));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_ormqr(value va, value vmn, value vormqr)
{
    CAMLparam3(va, vmn, vormqr);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    realtype *beta = REAL_ARRAY(Field(vormqr, 0));
    realtype *vv   = REAL_ARRAY(Field(vormqr, 1));
    realtype *vw   = REAL_ARRAY(Field(vormqr, 2));
    realtype *work = REAL_ARRAY(Field(vormqr, 3));

    denseORMQR(DDENSEMAT(va), m, n, beta, vv, vw, work);
    CAMLreturn0;
}

/* Band matrix functions */

CAMLprim value c_bandmatrix_new_band_mat(value vn, value vmu,
	value vml, value vsmu)
{
    CAMLparam4(vn, vmu, vml, vsmu);
    CAMLlocal1(vr);

    int n   = Int_val(vn);
    int mu  = Int_val(vmu);
    int ml  = Int_val(vml);
    int smu = Int_val(vsmu);

    DlsMat a = NewBandMat(n, mu, ml, smu);
    if (a == NULL)
	caml_failwith("Could not create Band Matrix.");
    mlsize_t approx_size = n * (smu + ml + 2) * sizeof(realtype);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vr = caml_alloc_final(2, &finalize_dlsmat, approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim void c_bandmatrix_copy(value va, value vb,
	value vcopymu, value vcopyml)
{
    CAMLparam4(va, vb, vcopymu, vcopyml);
    BandCopy(DLSMAT(va), DLSMAT(vb), Int_val(vcopymu), Int_val(vcopyml));
    CAMLreturn0;
}

CAMLprim void c_bandmatrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
    BandScale(Double_val(vc), DLSMAT(va));
    CAMLreturn0;
}

CAMLprim void c_bandmatrix_gbtrf(value va, value vp)
{
    CAMLparam2(va, vp);
    BandGBTRF(DLSMAT(va), LONG_ARRAY(vp));
    CAMLreturn0;
}

CAMLprim void c_bandmatrix_gbtrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);
    BandGBTRS(DLSMAT(va), LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim value c_bandmatrix_get(value vmatrix, value vi, value vj)
{
    CAMLparam3(vmatrix, vi, vj);
    DlsMat m = DLSMAT(vmatrix);

    int i = Int_val(vi);
    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Bandmatrix.get: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Bandmatrix.get: invalid j");
#endif

    realtype v = BAND_ELEM(m, i, j);
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_bandmatrix_set(value vmatrix, value vi, value vj, value v)
{
    CAMLparam4(vmatrix, vi, vj, v);
    DlsMat m = DLSMAT(vmatrix);

    int i = Int_val(vi);
    int j = Int_val(vj);

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

    DlsMat m = DLSMAT(vmatrix);

    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    if (j < 0 || j >= m->N)
	caml_invalid_argument("Bandmatrix.Col.get_col: invalid j");

    r = caml_alloc(4, Abstract_tag);
    Store_field(r, 2, Val_int(m->mu));
    Store_field(r, 3, Val_int(m->ml));
#else
    r = caml_alloc(2, Abstract_tag);
#endif

    Store_field(r, 0, (value)BAND_COL(m, j));
    Store_field(r, 1, vmatrix); /* avoid gc of underlying matrix! */
    CAMLreturn(r);
}

CAMLprim value c_bandmatrix_col_get(value vbcol, value vi, value vj)
{
    CAMLparam3(vbcol, vi, vj);

    realtype *bcol = (realtype *)Field(vbcol, 0);

    int i = Int_val(vi);
    int j = Int_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    int mu = Int_val(Field(vbcol, 2));
    int ml = Int_val(Field(vbcol, 3));

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
    int mu = Int_val(Field(vbcol, 2));
    int ml = Int_val(Field(vbcol, 3));

    if (i < (j - mu) || i > (j + ml))
	caml_invalid_argument("Bandmatrix.Col.set: invalid i");
#endif

    BAND_COL_ELEM(bcol, i, j) = Double_val(ve);

    CAMLreturn(Val_unit);
}

/* Array Band matrix functions */

#define DBANDMAT(v) (*(realtype ***)(Data_custom_val(v)))

static void finalize_arraybandmatrix(value va)
{
    destroyMat(DBANDMAT(va));
}

CAMLprim value c_arraybandmatrix_new_band_mat(value vn, value vsmu, value vml)
{
    CAMLparam3(vn, vsmu, vml);
    CAMLlocal1(vr);

    int n   = Int_val(vn);
    int smu = Int_val(vsmu);
    int ml  = Int_val(vml);

    realtype **a = newBandMat(n, smu, ml);
    if (a == NULL)
	caml_failwith("Could not create Array Band Matrix.");
    mlsize_t approx_size = n * (smu + ml + 2) * sizeof(realtype);

    vr = caml_alloc_final(2, &finalize_arraybandmatrix,
			  approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim void c_arraybandmatrix_copy(value va, value vb, value vsizes)
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

CAMLprim void c_arraybandmatrix_scale(value vc, value va, value vsizes)
{
    CAMLparam3(vc, va, vsizes);

    long int n   = Long_val(Field(vsizes, 0));
    long int mu  = Long_val(Field(vsizes, 1));
    long int ml  = Long_val(Field(vsizes, 2));
    long int smu = Long_val(Field(vsizes, 3));

    bandScale(Double_val(vc), DBANDMAT(va), n, mu, ml, smu);
    CAMLreturn0;
}

CAMLprim void c_arraybandmatrix_add_identity(value va, value vn, value vsmu)
{
    CAMLparam3(va, vn, vsmu);

    bandAddIdentity(DBANDMAT(va), Long_val(vn), Long_val(vsmu));
    CAMLreturn0;
}

CAMLprim void c_arraybandmatrix_gbtrf(value va, value vsizes, value vp)
{
    CAMLparam3(va, vsizes, vp);

    long int n   = Long_val(Field(vsizes, 0));
    long int mu  = Long_val(Field(vsizes, 1));
    long int ml  = Long_val(Field(vsizes, 2));
    long int smu = Long_val(Field(vsizes, 3));

    bandGBTRF(DBANDMAT(va), n, mu, ml, smu, LONG_ARRAY(vp));
    CAMLreturn0;
}

CAMLprim void c_arraybandmatrix_gbtrs(value va, value vsizes, value vp, value vb)
{
    CAMLparam4(va, vsizes, vp, vb);

    long int n   = Long_val(Field(vsizes, 0));
    long int smu = Long_val(Field(vsizes, 1));
    long int ml  = Long_val(Field(vsizes, 2));

    bandGBTRS(DBANDMAT(va), n, smu, ml, LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

