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

#include "sundials_ml.h"

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

    int m = Long_val(vm);
    int n = Long_val(vn);

    DlsMat a = NewDenseMat(m, n);
    if (a == NULL)
	caml_failwith("Could not create Dense Matrix.");
    mlsize_t approx_size = m * n * sizeof(realtype);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vr = caml_alloc_final(2, &finalize_dlsmat, approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim value c_densematrix_size(value va)
{
    CAMLparam1(va);
    CAMLlocal1(vr);

    DlsMat ma = DLSMAT(va);
    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_long(ma->M));
    Store_field(vr, 1, Val_long(ma->N));

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
	caml_raise_with_arg(*caml_named_value("dls_ZeroDiagonalElement"),
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

    int i = Long_val(vi);
    int j = Long_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("DenseMatrix.get: invalid i.");
    if (j < 0 || j >= m->N) caml_invalid_argument("DenseMatrix.get: invalid j.");
#endif

    realtype v = DENSE_ELEM(m, i, j);
    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_densematrix_set(value vmatrix, value vi, value vj, value v)
{
    CAMLparam4(vmatrix, vi, vj, v);
    DlsMat m = DLSMAT(vmatrix);

    int i = Long_val(vi);
    int j = Long_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("DenseMatrix.set: invalid i.");
    if (j < 0 || j >= m->N) caml_invalid_argument("DenseMatrix.set: invalid j.");
#endif

    DENSE_ELEM(m, i, j) = Double_val(v);
    CAMLreturn(caml_copy_double(v));
}

/* Array dense matrix functions */

CAMLprim void c_arraydensematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    denseScale(Double_val(vc), ARRAY2_ACOLS(va), m, n);
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_add_identity(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.add_identity: matrix not square.");
#endif

    denseAddIdentity(ARRAY2_ACOLS(va), m);
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if CHECK_MATRIX_ACCESS == 1
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("ArrayDenseMatrix.getrf: p is too small.");
#endif

    int r = denseGETRF(ARRAY2_ACOLS(va), m, n, LONG_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(*caml_named_value("dls_ZeroDiagonalElement"),
			    Val_int(r));
    }
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_getrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: matrix not square.");
    if (ARRAY1_LEN(vb) < n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: b is too small.");
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: p is too small.");
#endif

    denseGETRS(ARRAY2_ACOLS(va), m, LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_potrf(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.potrf: matrix not square");
#endif

    densePOTRF(ARRAY2_ACOLS(va), m);
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_potrs(value va, value vb)
{
    CAMLparam2(va, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.potrs: matrix not square.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("ArrayDenseMatrix.potrs: b is too small.");
#endif

    densePOTRS(ARRAY2_ACOLS(va), m, REAL_ARRAY(vb));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_geqrf(value va, value vbeta, value vv)
{
    CAMLparam3(va, vbeta, vv);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if CHECK_MATRIX_ACCESS == 1
    if (m < n)
	caml_invalid_argument("ArrayDenseMatrix.geqrf: fewer rows than columns.");
    if (ARRAY1_LEN(vbeta) < n)
	caml_invalid_argument("ArrayDenseMatrix.geqrf: beta is too small.");
    if (ARRAY1_LEN(vv) < m)
	caml_invalid_argument("ArrayDenseMatrix.geqrf: work is too small.");
#endif

    denseGEQRF(ARRAY2_ACOLS(va), m, n, REAL_ARRAY(vbeta), REAL_ARRAY(vv));
    CAMLreturn0;
}

CAMLprim void c_arraydensematrix_ormqr(value va, value vormqr)
{
    CAMLparam2(va, vormqr);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    realtype *beta = REAL_ARRAY(Field(vormqr, 0));
    realtype *vv   = REAL_ARRAY(Field(vormqr, 1));
    realtype *vw   = REAL_ARRAY(Field(vormqr, 2));
    realtype *work = REAL_ARRAY(Field(vormqr, 3));

#if CHECK_MATRIX_ACCESS == 1
    if (m < n)
	caml_invalid_argument("ArrayDenseMatrix.ormqr: fewer rows than columns.");
    if (ARRAY1_LEN(Field(vormqr, 0)) < n)
	caml_invalid_argument("ArrayDenseMatrix.ormqr: beta is too small.");
    if (ARRAY1_LEN(Field(vormqr, 1)) < n)
	caml_invalid_argument("ArrayDenseMatrix.ormqr: v is too small.");
    if (ARRAY1_LEN(Field(vormqr, 2)) < m)
	caml_invalid_argument("ArrayDenseMatrix.ormqr: w is too small.");
    if (ARRAY1_LEN(Field(vormqr, 3)) < m)
	caml_invalid_argument("ArrayDenseMatrix.ormqr: work is too small.");
#endif

    denseORMQR(ARRAY2_ACOLS(va), m, n, beta, vv, vw, work);
    CAMLreturn0;
}

/* Band matrix functions */

CAMLprim value c_bandmatrix_new_band_mat(value vn, value vmu,
	value vml, value vsmu)
{
    CAMLparam4(vn, vmu, vml, vsmu);
    CAMLlocal1(vr);

    int n   = Long_val(vn);
    int mu  = Long_val(vmu);
    int ml  = Long_val(vml);
    int smu = Long_val(vsmu);

    DlsMat a = NewBandMat(n, mu, ml, smu);
    if (a == NULL)
	caml_failwith("Could not create Band Matrix.");
    mlsize_t approx_size = n * (smu + ml + 2) * sizeof(realtype);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vr = caml_alloc_final(2, &finalize_dlsmat, approx_size, approx_size * 20);
    Store_field(vr, 1, (value)a);

    CAMLreturn(vr);
}

CAMLprim value c_bandmatrix_size(value va)
{
    CAMLparam1(va);
    CAMLlocal1(vr);

    DlsMat ma = DLSMAT(va);
    vr = caml_alloc_tuple(3);
    Store_field(vr, 0, Val_long(ma->mu));
    Store_field(vr, 1, Val_long(ma->s_mu));
    Store_field(vr, 2, Val_long(ma->ml));

    CAMLreturn(vr);
}

CAMLprim void c_bandmatrix_copy(value va, value vb,
	value vcopymu, value vcopyml)
{
    CAMLparam4(va, vb, vcopymu, vcopyml);
    BandCopy(DLSMAT(va), DLSMAT(vb), Long_val(vcopymu), Long_val(vcopyml));
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

    int i = Long_val(vi);
    int j = Long_val(vj);

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

    int i = Long_val(vi);
    int j = Long_val(vj);

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

    int j = Long_val(vj);

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

    int i = Long_val(vi);
    int j = Long_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    int mu = Long_val(Field(vbcol, 2));
    int ml = Long_val(Field(vbcol, 3));

    if (i < (j - mu) || i > (j + ml))
	caml_invalid_argument("Bandmatrix.Col.get: invalid i");
#endif

    CAMLreturn(caml_copy_double(BAND_COL_ELEM(bcol, i, j)));
}

CAMLprim value c_bandmatrix_col_set(value vbcol, value vi, value vj, value ve)
{
    CAMLparam4(vbcol, vi, vj, ve);

    realtype *bcol = (realtype *)Field(vbcol, 0);

    int i = Long_val(vi);
    int j = Long_val(vj);

#if CHECK_MATRIX_ACCESS == 1
    int mu = Long_val(Field(vbcol, 2));
    int ml = Long_val(Field(vbcol, 3));

    if (i < (j - mu) || i > (j + ml))
	caml_invalid_argument("Bandmatrix.Col.set: invalid i");
#endif

    BAND_COL_ELEM(bcol, i, j) = Double_val(ve);

    CAMLreturn(Val_unit);
}

/* Array Band matrix functions */

CAMLprim void c_arraybandmatrix_copy(value va, value vb, value vsizes)
{
    CAMLparam3(va, vb, vsizes);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat am = ba->dim[1];
    intnat an = ba->dim[0];

    struct caml_ba_array *bb = ARRAY2_DATA(vb);
    intnat bm = bb->dim[1];
    intnat bn = bb->dim[0];

    int a_smu  = Long_val(Field(vsizes, 0));
    int b_smu  = Long_val(Field(vsizes, 1));
    int copymu = Long_val(Field(vsizes, 2));
    int copyml = Long_val(Field(vsizes, 3));

#if CHECK_MATRIX_ACCESS == 1
    if (am != an)
	caml_invalid_argument("ArrayBandMatrix.copy: source matrix not square.");
    if ((am != bm) || (bm != bn))
	caml_invalid_argument("ArrayBandMatrix.copy: matrix sizes differ.");
#endif

    bandCopy(ARRAY2_ACOLS(va), ARRAY2_ACOLS(vb), am, a_smu, b_smu,
	     copymu, copyml);
    CAMLreturn0;
}

CAMLprim void c_arraybandmatrix_scale(value vc, value va, value vsizes)
{
    CAMLparam3(vc, va, vsizes);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayBandMatrix.scale: matrix not square.");
#endif

    bandScale(Double_val(vc), ARRAY2_ACOLS(va), m, mu, ml, smu);
    CAMLreturn0;
}

CAMLprim void c_arraybandmatrix_add_identity(value va, value vsmu)
{
    CAMLparam2(va, vsmu);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayBandMatrix.add_identity: matrix not square.");
#endif

    bandAddIdentity(ARRAY2_ACOLS(va), n, Long_val(vsmu));
    CAMLreturn0;
}

CAMLprim void c_arraybandmatrix_gbtrf(value va, value vsizes, value vp)
{
    CAMLparam3(va, vsizes, vp);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: matrix not square.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: p is too small.");
#endif

    bandGBTRF(ARRAY2_ACOLS(va), m, mu, ml, smu, LONG_ARRAY(vp));
    CAMLreturn0;
}

CAMLprim void c_arraybandmatrix_gbtrs(value va, value vsizes, value vp, value vb)
{
    CAMLparam4(va, vsizes, vp, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    long int smu = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));

#if CHECK_MATRIX_ACCESS == 1
    if (m != n)
	caml_invalid_argument("ArrayBandMatrix.gbtrs: matrix not square.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: p is too small.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: b is too small.");
#endif

    bandGBTRS(ARRAY2_ACOLS(va), m, smu, ml, LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn0;
}

