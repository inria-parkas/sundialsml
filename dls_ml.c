/***********************************************************************
 *                                                                     *
 *              Ocaml interface to Sundials CVODE solver               *
 *                                                                     *
 *           Timothy Bourke (INRIA) and Marc Pouzet (LIENS)            *
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under the terms of the GNU Library General Public License, with    *
 *  the special exception on linking described in file LICENSE.        *
 *                                                                     *
 ***********************************************************************/

/* Sundials interface functions that do not involve NVectors. */

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
    DlsMat mat = DLSMAT(va);
    DestroyMat(DLSMAT(va));
}

CAMLprim value c_densematrix_new_dense_mat(value vmn)
{
    CAMLparam1(vmn);
    CAMLlocal1(vr);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

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
 
CAMLprim value c_densematrix_get(value vmatrix, value vij)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = DLSMAT(vmatrix);

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
    CAMLparam3(vmatrix, vij, v);
    DlsMat m = DLSMAT(vmatrix);

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

#define DDENSEMAT(v) (*(realtype ***)(Data_custom_val(v)))

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
    if (a == NULL)
	caml_failwith("Could not create Direct Dense Matrix.");
    mlsize_t approx_size = m * n * sizeof(realtype);

    vr = caml_alloc_final(2, &finalize_direct_densemat,
			  approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim value c_densematrix_direct_get(value va, value vij)
{
    CAMLparam2(va, vij);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

    CAMLreturn(caml_copy_double(DDENSEMAT(va)[j][i]));
}

CAMLprim void c_densematrix_direct_set(value va, value vij, value vv)
{
    CAMLparam3(va, vij, vv);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

    DDENSEMAT(va)[j][i] = Double_val(vv);

    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_copy(value va, value vb, value vmn)
{
    CAMLparam3(va, vb, vmn);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    denseCopy(DDENSEMAT(va), DDENSEMAT(vb), m, n);
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_scale(value vc, value va, value vmn)
{
    CAMLparam3(vc, va, vmn);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    denseScale(Double_val(vc), DDENSEMAT(va), m, n);
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_add_identity(value va, value vn)
{
    CAMLparam2(va, vn);
    denseAddIdentity(DDENSEMAT(va), Int_val(vn));
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_getrf(value va, value vmn, value vp)
{
    CAMLparam3(va, vmn, vp);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    int r = denseGETRF(DDENSEMAT(va), m, n, LONG_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(*caml_named_value("cvode_ZeroDiagonalElement"),
			    Val_int(r));
    }
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_getrs(value va, value vn,
	value vp, value vb)
{
    CAMLparam4(va, vn, vp, vb);
    denseGETRS(DDENSEMAT(va), Int_val(vn), LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_potrf(value va, value vm)
{
    CAMLparam2(va, vm);
    densePOTRF(DDENSEMAT(va), Int_val(vm));
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_potrs(value va, value vm, value vb)
{
    CAMLparam3(va, vm, vb);
    densePOTRS(DDENSEMAT(va), Int_val(vm), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_geqrf(value va, value vmn,
	value vbeta, value vv)
{
    CAMLparam4(va, vmn, vbeta, vv);

    int m = Int_val(Field(vmn, 0));
    int n = Int_val(Field(vmn, 1));

    denseGEQRF(DDENSEMAT(va), m, n, REAL_ARRAY(vbeta), REAL_ARRAY(vv));
    CAMLreturn0;
}

CAMLprim void c_densematrix_direct_ormqr(value va, value vmn, value vormqr)
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

CAMLprim value c_bandmatrix_new_band_mat(value vsizes)
{
    CAMLparam1(vsizes);
    CAMLlocal1(vr);

    int n   = Int_val(Field(vsizes, 0));
    int mu  = Int_val(Field(vsizes, 1));
    int ml  = Int_val(Field(vsizes, 2));
    int smu = Int_val(Field(vsizes, 3));

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

CAMLprim value c_bandmatrix_get(value vmatrix, value vij)
{
    CAMLparam2(vmatrix, vij);
    DlsMat m = DLSMAT(vmatrix);

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
    CAMLparam3(vmatrix, vij, v);
    DlsMat m = DLSMAT(vmatrix);

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

CAMLprim value c_bandmatrix_col_get(value vbcol, value vij)
{
    CAMLparam2(vbcol, vij);

    realtype *bcol = (realtype *)Field(vbcol, 0);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

#if CHECK_MATRIX_ACCESS == 1
    int mu = Int_val(Field(vbcol, 2));
    int ml = Int_val(Field(vbcol, 3));

    if (i < (j - mu) || i > (j + ml))
	caml_invalid_argument("Bandmatrix.Col.get: invalid i");
#endif

    CAMLreturn(caml_copy_double(BAND_COL_ELEM(bcol, i, j)));
}

CAMLprim value c_bandmatrix_col_set(value vbcol, value vij, value ve)
{
    CAMLparam3(vbcol, vij, ve);

    realtype *bcol = (realtype *)Field(vbcol, 0);

    int i = Int_val(Field(vij, 0));
    int j = Int_val(Field(vij, 1));

#if CHECK_MATRIX_ACCESS == 1
    int mu = Int_val(Field(vbcol, 2));
    int ml = Int_val(Field(vbcol, 3));

    if (i < (j - mu) || i > (j + ml))
	caml_invalid_argument("Bandmatrix.Col.set: invalid i");
#endif

    BAND_COL_ELEM(bcol, i, j) = Double_val(ve);

    CAMLreturn(Val_unit);
}

/* Band matrix direct functions */

#define DBANDMAT(v) (*(realtype ***)(Data_custom_val(v)))

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
    if (a == NULL)
	caml_failwith("Could not create Direct Band Matrix.");
    mlsize_t approx_size = n * (smu + ml + 2) * sizeof(realtype);

    vr = caml_alloc_final(2, &finalize_direct_bandmat,
			  approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim void c_bandmatrix_direct_copy(value va, value vb, value vsizes)
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

CAMLprim void c_bandmatrix_direct_scale(value vc, value va, value vsizes)
{
    CAMLparam3(vc, va, vsizes);

    long int n   = Long_val(Field(vsizes, 0));
    long int mu  = Long_val(Field(vsizes, 1));
    long int ml  = Long_val(Field(vsizes, 2));
    long int smu = Long_val(Field(vsizes, 3));

    bandScale(Double_val(vc), DBANDMAT(va), n, mu, ml, smu);
    CAMLreturn0;
}

CAMLprim void c_bandmatrix_direct_add_identity(value va, value vn, value vsmu)
{
    CAMLparam3(va, vn, vsmu);

    bandAddIdentity(DBANDMAT(va), Long_val(vn), Long_val(vsmu));
    CAMLreturn0;
}

CAMLprim void c_bandmatrix_direct_gbtrf(value va, value vsizes, value vp)
{
    CAMLparam3(va, vsizes, vp);

    long int n   = Long_val(Field(vsizes, 0));
    long int mu  = Long_val(Field(vsizes, 1));
    long int ml  = Long_val(Field(vsizes, 2));
    long int smu = Long_val(Field(vsizes, 3));

    bandGBTRF(DBANDMAT(va), n, mu, ml, smu, LONG_ARRAY(vp));
    CAMLreturn0;
}

CAMLprim void c_bandmatrix_direct_gbtrs(value va, value vsizes, value vp, value vb)
{
    CAMLparam4(va, vsizes, vp, vb);

    long int n   = Long_val(Field(vsizes, 0));
    long int smu = Long_val(Field(vsizes, 1));
    long int ml  = Long_val(Field(vsizes, 2));

    bandGBTRS(DBANDMAT(va), n, smu, ml, LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}
