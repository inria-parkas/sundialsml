/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) * 
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
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
#include "dls_ml.h"

CAMLprim value c_dls_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (DLS, exns);
    CAMLreturn (Val_unit);
}

/* Shared matrix functions */

static void finalize_dlsmat(value va)
{
    DestroyMat(DLSMAT(va));
}

/* Dense matrix functions */

CAMLprim value c_dls_dense_wrap(DlsMat a, int finalize)
{
    CAMLparam0();
    CAMLlocal3(vv, va, vr);
    mlsize_t approx_size = a->ldim * a->N * sizeof(realtype) + 1;

    va = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, a->data, a->N, a->ldim);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vv = caml_alloc_final(2, finalize ? &finalize_dlsmat : NULL,
			  approx_size, approx_size * 20);
    DLSMAT(vv) = a;

    vr = caml_alloc_tuple(3);
    Store_field(vr, RECORD_DLS_DENSEMATRIX_PAYLOAD, va);
    Store_field(vr, RECORD_DLS_DENSEMATRIX_DLSMAT, vv);
    Store_field(vr, RECORD_DLS_DENSEMATRIX_VALID, Val_bool(1));

    CAMLreturn(vr);
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

    CAMLreturn(c_dls_dense_wrap(a, 1));
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

CAMLprim value c_densematrix_print_mat(value va)
{
    CAMLparam1(va);
    PrintMat(DLSMAT(va));
    fflush(stdout);
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_set_to_zero(value va)
{
    CAMLparam1(va);
    SetToZero(DLSMAT(va));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_add_identity(value va)
{
    CAMLparam1(va);
    AddIdentity(DLSMAT(va));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_copy(value va, value vb)
{
    CAMLparam2(va, vb);
    DenseCopy(DLSMAT(va), DLSMAT(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
    DenseScale(Double_val(vc), DLSMAT(va));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);
    int r = DenseGETRF(DLSMAT(va), LONG_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(DLS_EXN(ZeroDiagonalElement),
			    Val_int(r));
    }
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_getrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);
    DenseGETRS(DLSMAT(va), LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_potrf(value va)
{
    CAMLparam1(va);
    DensePOTRF(DLSMAT(va));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_potrs(value va, value vb)
{
    CAMLparam2(va, vb);
    DensePOTRS(DLSMAT(va), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_geqrf(value va, value vbeta, value vwork)
{
    CAMLparam3(va, vbeta, vwork);
    DenseGEQRF(DLSMAT(va), REAL_ARRAY(vbeta), REAL_ARRAY(vwork));
    CAMLreturn (Val_unit);
}

CAMLprim value c_densematrix_ormqr(value va, value vormqr)
{
    CAMLparam2(va, vormqr);

    realtype *beta = REAL_ARRAY(Field(vormqr, 0));
    realtype *vv   = REAL_ARRAY(Field(vormqr, 1));
    realtype *vw   = REAL_ARRAY(Field(vormqr, 2));
    realtype *work = REAL_ARRAY(Field(vormqr, 3));

    DenseORMQR(DLSMAT(va), beta, vv, vw, work);
    CAMLreturn (Val_unit);
}
 
CAMLprim value c_densematrix_get(value vmatrix, value vi, value vj)
{
    CAMLparam3(vmatrix, vi, vj);
    DlsMat m = DLSMAT(vmatrix);

    int i = Long_val(vi);
    int j = Long_val(vj);

#if SUNDIALS_ML_SAFE == 1
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

#if SUNDIALS_ML_SAFE == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("DenseMatrix.set: invalid i.");
    if (j < 0 || j >= m->N) caml_invalid_argument("DenseMatrix.set: invalid j.");
#endif

    DENSE_ELEM(m, i, j) = Double_val(v);
    CAMLreturn(caml_copy_double(v));
}

/* Array dense matrix functions */

CAMLprim value c_arraydensematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    denseScale(Double_val(vc), ARRAY2_ACOLS(va), m, n);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_add_identity(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];

    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.add_identity: matrix not square.");
#endif

    denseAddIdentity(ARRAY2_ACOLS(va), m);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("ArrayDenseMatrix.getrf: p is too small.");
#endif

    int r = denseGETRF(ARRAY2_ACOLS(va), m, n, LONG_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(DLS_EXN(ZeroDiagonalElement),
			    Val_int(r));
    }
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_getrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: matrix not square.");
    if (ARRAY1_LEN(vb) < n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: b is too small.");
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: p is too small.");
#endif

    denseGETRS(ARRAY2_ACOLS(va), m, LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_getrs_off(value va, value vp,
					    value vb, value vboff)
{
    CAMLparam4(va, vp, vb, vboff);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat boff = Int_val(vboff);

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: matrix not square.");
    if (ARRAY1_LEN(vb) - boff < n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: b is too small.");
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("ArrayDenseMatrix.getrs: p is too small.");
#endif

    denseGETRS(ARRAY2_ACOLS(va), m, LONG_ARRAY(vp), REAL_ARRAY(vb) + boff);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_potrf(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.potrf: matrix not square");
#endif

    densePOTRF(ARRAY2_ACOLS(va), m);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_potrs(value va, value vb)
{
    CAMLparam2(va, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("ArrayDenseMatrix.potrs: matrix not square.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("ArrayDenseMatrix.potrs: b is too small.");
#endif

    densePOTRS(ARRAY2_ACOLS(va), m, REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_geqrf(value va, value vbeta, value vv)
{
    CAMLparam3(va, vbeta, vv);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (m < n)
	caml_invalid_argument("ArrayDenseMatrix.geqrf: fewer rows than columns.");
    if (ARRAY1_LEN(vbeta) < n)
	caml_invalid_argument("ArrayDenseMatrix.geqrf: beta is too small.");
    if (ARRAY1_LEN(vv) < m)
	caml_invalid_argument("ArrayDenseMatrix.geqrf: work is too small.");
#endif

    denseGEQRF(ARRAY2_ACOLS(va), m, n, REAL_ARRAY(vbeta), REAL_ARRAY(vv));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_ormqr(value va, value vormqr)
{
    CAMLparam2(va, vormqr);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    realtype *beta = REAL_ARRAY(Field(vormqr, 0));
    realtype *vv   = REAL_ARRAY(Field(vormqr, 1));
    realtype *vw   = REAL_ARRAY(Field(vormqr, 2));
    realtype *work = REAL_ARRAY(Field(vormqr, 3));

#if SUNDIALS_ML_SAFE == 1
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
    CAMLreturn (Val_unit);
}

/* Band matrix functions */

CAMLprim value c_dls_band_wrap(DlsMat a, int finalize)
{
    CAMLparam0();
    CAMLlocal3(vv, va, vr);
    mlsize_t approx_size = a->ldim * a->N * sizeof(realtype) + 2;

    va = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, a->data, a->N, a->ldim);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vv = caml_alloc_final(2, finalize ? &finalize_dlsmat : NULL,
			  approx_size, approx_size * 20);
    DLSMAT(vv) = a;

    vr = caml_alloc_tuple(4);
    Store_field(vr, RECORD_DLS_BANDMATRIX_PAYLOAD, va);
    Store_field(vr, RECORD_DLS_BANDMATRIX_DLSMAT, vv);
    Store_field(vr, RECORD_DLS_BANDMATRIX_SMU, Val_int(a->s_mu));
    Store_field(vr, RECORD_DLS_BANDMATRIX_VALID, Val_bool(1));

    CAMLreturn(vr);
}

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

    CAMLreturn(c_dls_band_wrap(a, 1));
}

CAMLprim value c_bandmatrix_size(value va)
{
    CAMLparam1(va);
    CAMLlocal1(vr);

    DlsMat ma = DLSMAT(va);
    vr = caml_alloc_tuple(3);
    Store_field(vr, 0, Val_long(ma->N));
    Store_field(vr, 1, Val_long(ma->mu));
    Store_field(vr, 2, Val_long(ma->ml));
    Store_field(vr, 3, Val_long(ma->s_mu));

    CAMLreturn(vr);
}

CAMLprim value c_bandmatrix_copy(value va, value vb,
				 value vcopymu, value vcopyml)
{
    CAMLparam4(va, vb, vcopymu, vcopyml);
    BandCopy(DLSMAT(va), DLSMAT(vb), Long_val(vcopymu), Long_val(vcopyml));
    CAMLreturn (Val_unit);
}

CAMLprim value c_bandmatrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
    BandScale(Double_val(vc), DLSMAT(va));
    CAMLreturn (Val_unit);
}

CAMLprim value c_bandmatrix_gbtrf(value va, value vp)
{
    CAMLparam2(va, vp);
    BandGBTRF(DLSMAT(va), LONG_ARRAY(vp));
    CAMLreturn (Val_unit);
}

CAMLprim value c_bandmatrix_gbtrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);
    BandGBTRS(DLSMAT(va), LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_bandmatrix_get(value vmatrix, value vi, value vj)
{
    CAMLparam3(vmatrix, vi, vj);
    DlsMat m = DLSMAT(vmatrix);

    int i = Long_val(vi);
    int j = Long_val(vj);

#if SUNDIALS_ML_SAFE == 1
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

#if SUNDIALS_ML_SAFE == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("Bandmatrix.set: invalid i");
    if (j < 0 || j >= m->N) caml_invalid_argument("Bandmatrix.set: invalid j");
#endif

    BAND_ELEM(m, i, j) = Double_val(v);
    CAMLreturn(caml_copy_double(v));
}

/* Array Band matrix functions */

CAMLprim value c_arraybandmatrix_copy(value va, value vb, value vsizes)
{
    CAMLparam3(va, vb, vsizes);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat am = ba->dim[0];

    int a_smu  = Long_val(Field(vsizes, 0));
    int b_smu  = Long_val(Field(vsizes, 1));
    int copymu = Long_val(Field(vsizes, 2));
    int copyml = Long_val(Field(vsizes, 3));

#if SUNDIALS_ML_SAFE == 1
    intnat an = ba->dim[1];
    struct caml_ba_array *bb = ARRAY2_DATA(vb);

    intnat bm = bb->dim[0];
    intnat bn = bb->dim[1];

    if (an < copymu + copyml + 1)
	caml_invalid_argument("ArrayBandMatrix.copy: source matrix too small.");
    if (bn < copymu + copyml + 1)
	caml_invalid_argument("ArrayBandMatrix.copy: destination matrix too small.");
    if ((am != bm) || (bm != bn))
	caml_invalid_argument("ArrayBandMatrix.copy: matrix sizes differ.");
#endif

    bandCopy(ARRAY2_ACOLS(va), ARRAY2_ACOLS(vb), am, a_smu, b_smu,
	     copymu, copyml);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraybandmatrix_scale(value vc, value va, value vsizes)
{
    CAMLparam3(vc, va, vsizes);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("ArrayBandMatrix.scale: matrix badly sized.");
#endif

    bandScale(Double_val(vc), ARRAY2_ACOLS(va), m, mu, ml, smu);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraybandmatrix_add_identity(value va, value vsmu)
{
    CAMLparam2(va, vsmu);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];
    intnat smu = Long_val(vsmu);

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n <= smu)
	caml_invalid_argument("ArrayBandMatrix.add_identity: matrix badly sized.");
#endif

    bandAddIdentity(ARRAY2_ACOLS(va), m, smu);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraybandmatrix_gbtrf(value va, value vsizes, value vp)
{
    CAMLparam3(va, vsizes, vp);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: p is too small.");
#endif

    bandGBTRF(ARRAY2_ACOLS(va), m, mu, ml, smu, LONG_ARRAY(vp));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraybandmatrix_gbtrs(value va, value vsizes, value vp, value vb)
{
    CAMLparam4(va, vsizes, vp, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int smu = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < smu + ml + 1)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: p is too small.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("ArrayBandMatrix.gbtrf: b is too small.");
#endif

    bandGBTRS(ARRAY2_ACOLS(va), m, smu, ml, LONG_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

