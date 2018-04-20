/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2018 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include "../config.h"

#include <stdio.h>
#include <stdbool.h>

#if SUNDIALS_LIB_VERSION >= 300
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_sparse.h>

#else // SUNDIALS_LIB_VERSION < 300
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>
#endif

#include <sundials/sundials_dense.h>
#include <sundials/sundials_band.h>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>

#include "../sundials/sundials_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../lsolvers/matrix_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

extern CAMLprim value caml_ba_blit(value vsrc, value vdst);

CAMLprim void ml_mat_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (MATRIX, exns);
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Matrix.Dense
 */

static void finalize_mat_content_dense(value va)
{
    MAT_CONTENT_DENSE_TYPE content = MAT_CONTENT_DENSE(va);

    if (content->cols != NULL)
	free(content->cols);
    free(content);
    // A->data is destroyed by the associated bigarray finalizer
}

CAMLprim value ml_matrix_dense_create(value vm, value vn)
{
    CAMLparam2(vm, vn);
    CAMLlocal3(vdata, vcptr, vr);
    sundials_ml_index j;
    MAT_CONTENT_DENSE_TYPE content = NULL;
    sundials_ml_index m = Long_val(vm);
    sundials_ml_index n = Long_val(vn);

    // columns first
    vdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, NULL, n, m);

    // Setup the C-side content
    content = (MAT_CONTENT_DENSE_TYPE) malloc(sizeof *content);
    if (content == NULL) caml_raise_out_of_memory();

    content->M = m;
    content->N = n;
    content->ldata = m*n;
    content->data = REAL_ARRAY(vdata);
    content->cols = NULL;
    content->cols = (realtype **) malloc(n * sizeof(realtype *));
    if (content->cols == NULL) {
	free(content);
	caml_raise_out_of_memory();
    }
    for (j=0; j<n; j++)
	content->cols[j] = content->data + j * m;

#if SUNDIALS_LIB_VERSION < 300
    content->ldim = m;
    content->type = SUNDIALS_DENSE;
#endif

    // Setup the OCaml-side
    vcptr = caml_alloc_final(1, &finalize_mat_content_dense, 1, 20);
    MAT_CONTENT_DENSE(vcptr) = content;

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vdata);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR, vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID, Val_bool(1));

    CAMLreturn(vr);
}

#if SUNDIALS_LIB_VERSION < 300
CAMLprim value c_matrix_dense_wrap(DlsMat a)
{
    CAMLparam0();
    CAMLlocal3(vcptr, vcontent, vr);

    vcontent = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, a->data, a->N, a->ldim);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vcptr = caml_alloc_final(1, NULL, 1, 20);
    DLSMAT(vcptr) = a;

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vcontent);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR,  vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID,   Val_bool(1));

    CAMLreturn(vr);
}
#endif

/*
CAMLprim value ml_matrix_dense_get(value vcptr, value vi, value vj)
{
    CAMLparam3(vcptr, vi, vj);
    sundials_ml_index i = Long_val(vi);
    sundials_ml_index j = Long_val(vj);
    MAT_CONTENT_DENSE_TYPE m = MAT_CONTENT_DENSE(vcptr);

#if SUNDIALS_ML_SAFE == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("i");
    if (j < 0 || j >= m->N) caml_invalid_argument("j");
#endif

    CAMLreturn(caml_copy_double(m->cols[j][i]));
}
*/

/*
CAMLprim void ml_matrix_dense_set(value vcptr, value vi, value vj, value vv)
{
    CAMLparam4(vcptr, vi, vj, vv);
    sundials_ml_index i = Long_val(vi);
    sundials_ml_index j = Long_val(vj);
    MAT_CONTENT_DENSE_TYPE m = MAT_CONTENT_DENSE(vcptr);

#if SUNDIALS_ML_SAFE == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("i");
    if (j < 0 || j >= m->N) caml_invalid_argument("j");
#endif

    m->cols[j][i] = Double_val(vv);
    CAMLreturn0;
}
*/

CAMLprim void ml_matrix_dense_scale_add(value vc, value vcptra, value vcptrb)
{
    CAMLparam3(vc, vcptra, vcptrb);
    realtype c = Double_val(vc);
    int i, j;
    MAT_CONTENT_DENSE_TYPE contenta = MAT_CONTENT_DENSE(vcptra);
    MAT_CONTENT_DENSE_TYPE contentb = MAT_CONTENT_DENSE(vcptrb);

#if SUNDIALS_ML_SAFE == 1
    if ((contenta->M != contentb->M) || (contenta->N != contentb->N))
	caml_raise_constant(MATRIX_EXN(IncompatibleArguments));
#endif

    // adapted from SUNMatScaleAdd_Dense
    for (j=0; j < contenta->N; j++)
	for (i=0; i < contenta->M; i++)
	    contenta->cols[j][i] =
		c * contenta->cols[j][i] + contentb->cols[j][i];

    CAMLreturn0;
}

CAMLprim void ml_matrix_dense_scale_addi(value vc, value vcptra)
{
    CAMLparam2(vc, vcptra);
    realtype c = Double_val(vc);
    int i, j;
    MAT_CONTENT_DENSE_TYPE a = MAT_CONTENT_DENSE(vcptra);

    // adapted from SUNMatScaleAddI_Dense
    for (j=0; j < a->N; j++)
	for (i=0; i < a->M; i++)
	{
	    a->cols[j][i] *= c;
	    if (i == j)
		a->cols[j][i] += 1.0;
	}

    CAMLreturn0;
}

CAMLprim void ml_matrix_dense_matvec(value vcptra, value vx, value vy)
{
    CAMLparam3(vcptra, vx, vy);
    sundials_ml_index i, j;
    realtype *col_j, *xd, *yd;
    MAT_CONTENT_DENSE_TYPE a = MAT_CONTENT_DENSE(vcptra);

    xd = REAL_ARRAY(vx);
    yd = REAL_ARRAY(vy);

    // Adapted from SUNMatMatvec_Dense
    for (i=0; i < a->M; i++)
	yd[i] = 0.0;
    for(j=0; j < a->N; j++) {
	col_j = a->cols[j];
	for (i=0; i < a->M; i++)
	    yd[i] += col_j[i] * xd[j];
    }

    CAMLreturn0;
}

CAMLprim value ml_matrix_dense_space(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_DENSE_TYPE content = MAT_CONTENT_DENSE(vcptr);

    // Directly adapted from SUNMatSpace_Dense
    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_long(content->ldata));
    Store_field(vr, 1, Val_long(3 + content->N));

    CAMLreturn(vr);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Matrix.Band
 */

static void finalize_mat_content_band(value vcptra)
{
    MAT_CONTENT_BAND_TYPE content = MAT_CONTENT_BAND(vcptra);

    if (content->cols != NULL)
	free(content->cols);
    free(content);
    // content->data is destroyed by the associated bigarray finalizer
}

static int matrix_band_create_vcptr(sundials_ml_index n,
	sundials_ml_index mu, sundials_ml_index ml, sundials_ml_index smu,
	value *pvdata, value *pvcptr)
{
    MAT_CONTENT_BAND_TYPE content = NULL;
    sundials_ml_index colSize = smu + ml + 1;
    sundials_ml_index j;

    // columns first
    *pvdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, NULL, n, colSize);

    // Setup the C-side content
    content = (MAT_CONTENT_BAND_TYPE) malloc(sizeof *content);
    if (content == NULL) return 0;

    content->M = n;
    content->N = n;
    content->mu = mu;
    content->ml = ml;
    content->s_mu = smu;
    content->ldim = colSize;
    content->ldata = n * colSize;
    content->data = REAL_ARRAY(*pvdata);
    content->ldim = colSize;

    content->cols = NULL;
    content->cols = (realtype **) malloc(n * sizeof(realtype *));
    if (content->cols == NULL) {
	free(content);
	return 0;
    }
    for (j=0; j < n; j++) content->cols[j] = content->data + j * colSize;

#if SUNDIALS_LIB_VERSION < 300
    content->type = SUNDIALS_BAND;
#endif

    // Setup the OCaml-side
    *pvcptr = caml_alloc_final(1, &finalize_mat_content_band, 1, 20);
    MAT_CONTENT_BAND(*pvcptr) = content;

    return 1;
}

static void matrix_band_realloc(sundials_ml_index n, sundials_ml_index mu,
			        sundials_ml_index ml, sundials_ml_index smu,
			        value va, bool free_cols)
{
    CAMLparam1(va);
    CAMLlocal5(vpayload, vcptr, vnewdata, vnewdims, vnewpayload);
    MAT_CONTENT_BAND_TYPE content = NULL;
    sundials_ml_index colSize = smu + ml + 1;
    sundials_ml_index j;
    realtype **newcols;

    vcptr    = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    vpayload = Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD);

    content = MAT_CONTENT_BAND(vcptr);

    // columns first
    vnewdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, NULL, n, colSize);

    newcols = (realtype **) malloc(n * sizeof(realtype *));
    if (newcols == NULL) caml_raise_out_of_memory();

    // Reconfigure the C-side content
    content->M = n;
    content->N = n;
    content->mu = mu;
    content->ml = ml;
    content->s_mu = smu;
    content->ldim = colSize;
    content->ldata = n * colSize;
    content->data = REAL_ARRAY(vnewdata); // old array is gc'ed by bigarray
    content->ldim = colSize;

    if (free_cols) free(content->cols);
    content->cols = newcols;
    for (j=0; j < n; j++) content->cols[j] = content->data + j * colSize;

    // Reconfigure the OCaml-side
    vnewdims = caml_alloc_tuple(RECORD_MAT_BANDDIMENSIONS_SIZE);
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_N,   Val_long(n));
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_MU,  Val_long(mu));
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_SMU, Val_long(smu));
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_ML,  Val_long(ml));

    vnewpayload = caml_alloc_tuple(RECORD_MAT_BANDDATA_SIZE);
    Store_field(vnewpayload, RECORD_MAT_BANDDATA_DATA, vnewdata);
    Store_field(vnewpayload, RECORD_MAT_BANDDATA_DIMS, vnewdims);

    Store_field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vnewpayload);

    CAMLreturn0;
}

static value ml_matrix_band_create_mat(sundials_ml_index n,
	sundials_ml_index mu, sundials_ml_index ml, sundials_ml_index smu)
{
    CAMLparam0();
    CAMLlocal5(vdata, vcptr, vr, vdims, vpayload);

    if (!matrix_band_create_vcptr(n, mu, ml, smu, &vdata, &vcptr))
	caml_raise_out_of_memory();

    vdims = caml_alloc_tuple(RECORD_MAT_BANDDIMENSIONS_SIZE);
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_N,   Val_long(n));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_MU,  Val_long(mu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_SMU, Val_long(smu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_ML,  Val_long(ml));

    vpayload = caml_alloc_tuple(RECORD_MAT_BANDDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DATA, vdata);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DIMS, vdims);

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR, vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID, Val_bool(1));

    CAMLreturn(vr);
}

CAMLprim value ml_matrix_band_create(value vdims)
{
    CAMLparam1(vdims);
    CAMLlocal4(vdata, vcptr, vr, vpayload);

    sundials_ml_index n   =
	Long_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_N));
    sundials_ml_index mu  =
	Long_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_MU));
    sundials_ml_index ml  =
	Long_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_ML));
    sundials_ml_index smu =
	Long_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_SMU));

    if (!matrix_band_create_vcptr(n, mu, ml, smu, &vdata, &vcptr))
	caml_raise_out_of_memory();

    vpayload = caml_alloc_tuple(RECORD_MAT_BANDDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DATA, vdata);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DIMS, vdims);

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR,  vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID,   Val_bool(1));

    CAMLreturn(vr);
}

/*
CAMLprim value ml_matrix_band_get(value vcptr, value vi, value vj)
{
    CAMLparam3(vcptr, vi, vj);
    sundials_ml_index i = Long_val(vi);
    sundials_ml_index j = Long_val(vj);
    MAT_CONTENT_BAND_TYPE m = MAT_CONTENT_BAND(vcptr);

    sundials_ml_index ii = i - j + m->s_mu;

#if SUNDIALS_ML_SAFE == 1
    if (ii < 0 || ii >= m->ldim) caml_invalid_argument("i");
    if (j < 0 || j >= m->N) caml_invalid_argument("j");
#endif

    CAMLreturn(caml_copy_double(m->cols[j][ii]));
}
*/

/*
CAMLprim void ml_matrix_band_set(value vcptr, value vi, value vj, value v)
{
    CAMLparam4(vcptr, vi, vj, v);
    sundials_ml_index i = Long_val(vi);
    sundials_ml_index j = Long_val(vj);
    MAT_CONTENT_BAND_TYPE m = MAT_CONTENT_BAND(vcptr);
    sundials_ml_index ii = i - j + m->s_mu;

#if SUNDIALS_ML_SAFE == 1
    if (ii < 0 || ii >= m->ldim) caml_invalid_argument("i");
    if (j < 0 || j >= m->N) caml_invalid_argument("j");
#endif

    m->cols[j][ii] = Double_val(v);
    CAMLreturn0;
}
*/

// Adapted directly from SUNMatCopy_Band
CAMLprim void ml_matrix_band_copy(value vcptra, value vb)
{
    CAMLparam2(vcptra, vb);
    CAMLlocal3(vcptrb, vpayloadb, vdatab);
    sundials_ml_index i, j;
    realtype *A_colj, *B_colj;
    MAT_CONTENT_BAND_TYPE A, B;

    vcptrb = Field(vb, RECORD_MAT_MATRIXCONTENT_RAWPTR);

    A = MAT_CONTENT_BAND(vcptra);
    B = MAT_CONTENT_BAND(vcptrb);

#if SUNDIALS_ML_SAFE == 1
    if (A->M != B->M || A->N != B->N)
	caml_raise_constant(MATRIX_EXN(IncompatibleArguments));
#endif

    /* Grow B if A's bandwidth is larger */
    if ( (A->mu > B->mu) || (A->ml > B->ml) )
	matrix_band_realloc(B->N, SUNMAX(B->mu, A->mu), SUNMAX(B->ml, A->ml),
		SUNMAX(B->s_mu, A->s_mu), vb, 1);

    /* Perform operation */
    for (i=0; i < B->ldata; i++)
	B->data[i] = 0.0;

    for (j=0; j < B->N; j++) {
	B_colj = B->cols[j] + B->s_mu;
	A_colj = A->cols[j] + A->s_mu;
	for (i=-A->mu; i <= A->ml; i++)
	    B_colj[i] = A_colj[i];
    }

    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
static int csmat_band_copy(SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal3(vcontenta, vcontentb, r);
    static value * closure_f = NULL;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = MAT_BACKLINK(B);

    if (closure_f == NULL) closure_f = caml_named_value("ml_matrix_band_copy");

    r = caml_callback2_exn(*closure_f,
			   Field(vcontenta, RECORD_MAT_MATRIXCONTENT_RAWPTR),
			   vcontentb);
    if (Is_exception_result(r))
	sundials_ml_warn_discarded_exn(Extract_exception (r), "Band.blit");

    CAMLreturnT(int, Is_exception_result(r));
}
#endif

// Adapted directly from SMScaleAddNew_Band
static void matrix_band_scale_add_new(value vc, value va, value vcptrb)
{
    CAMLparam3(vc, va, vcptrb);
    CAMLlocal2(vcptra, voldpayload);
    realtype c = Double_val(vc);
    sundials_ml_index i, j, A_N, A_ml, A_mu, A_s_mu, new_mu, new_ml;
    realtype **A_cols;
    realtype *A_colj, *B_colj, *C_colj;
    MAT_CONTENT_BAND_TYPE B, C;

    vcptra      = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    // ensure that the old data array is not gc-ed until we are done.
    voldpayload = Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD);

    B = MAT_CONTENT_BAND(vcptrb);
    C = MAT_CONTENT_BAND(vcptra);

    // take original values from A (before reallocation)
    A_cols = C->cols;
    A_N    = C->N;
    A_mu   = C->mu;
    A_s_mu = C->s_mu;
    A_ml   = C->ml;

    /* create new matrix large enough to hold both A and B */
    new_mu = SUNMAX(B->mu, A_mu);
    new_ml = SUNMAX(B->ml, A_ml);
    matrix_band_realloc(B->N, new_mu, new_ml,
			SUNMIN(A_N - 1, new_mu + new_ml), va, 0);

    /* scale/add c*A into new matrix */
    for (j=0; j < A_N; j++) {
	A_colj = A_cols[j] + A_s_mu;
	C_colj = C->cols[j] + C->s_mu;

	for (i= -A_mu; i <= A_ml; i++)
	    C_colj[i] = c * A_colj[i];
    }
    free(A_cols);

    /* add B into new matrix */
    for (j=0; j < B->N; j++) {
	B_colj = B->cols[j] + B->s_mu;
	C_colj = C->cols[j] + C->s_mu;

	for (i= - B->mu; i <= B->ml; i++)
	    C_colj[i] += B_colj[i];
    }

    CAMLreturn0;
}

CAMLprim void ml_matrix_band_scale_add(value vc, value va, value vcptrb)
{
    CAMLparam3(vc, va, vcptrb);
    CAMLlocal1(vcptra);
    realtype c;
    realtype *A_colj, *B_colj;
    MAT_CONTENT_BAND_TYPE A, B;
    sundials_ml_index i, j;

    vcptra = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    A = MAT_CONTENT_BAND(vcptra);
    B = MAT_CONTENT_BAND(vcptrb);

    if ( (B->mu > A->mu) || (B->ml > A->ml) ) {
	matrix_band_scale_add_new(vc, va, vcptrb);

    } else {
	c = Double_val(vc);

	for (j=0; j < A->N; j++) {
	    A_colj = A->cols[j] + A->s_mu;
	    B_colj = B->cols[j] + B->s_mu;
	    for (i = -B->mu; i <= B->ml; i++)
		A_colj[i] = c * A_colj[i] + B_colj[i];
	}
    }

    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
static int csmat_band_scale_add(realtype c, SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal3(vcontenta, vcontentb, r);
    static value * closure_f = NULL;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = MAT_BACKLINK(B);

    if (closure_f == NULL)
	closure_f = caml_named_value("ml_matrix_band_scale_add");

    r = caml_callback3_exn(*closure_f,
			   caml_copy_double(c),
			   vcontenta,
			   Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR));
    if (Is_exception_result(r))
	sundials_ml_warn_discarded_exn(Extract_exception (r), "Band.scale_add");

    CAMLreturnT(int, Is_exception_result(r));
}
#endif

CAMLprim void ml_matrix_band_scale_addi(value vc, value va)
{
    CAMLparam2(vc, va);
    int i, j;
    realtype *A_colj;
    realtype c = Double_val(vc);
    MAT_CONTENT_BAND_TYPE a = MAT_CONTENT_BAND(va);

    // adapted from SUNMatScaleAddI_Band
    for (j=0; j < a->N; j++) {
	A_colj = a->cols[j] + a->s_mu;
	for (i = -a->mu; i <= a->ml; i++)
	    A_colj[i] *= c;
	a->cols[j][a->s_mu] += 1.0;
    }

    CAMLreturn0;
}

CAMLprim void ml_matrix_band_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);
    sundials_ml_index i, j, is, ie;
    realtype *col_j, *xd, *yd;
    MAT_CONTENT_BAND_TYPE a = MAT_CONTENT_BAND(va);

    xd = REAL_ARRAY(vx);
    yd = REAL_ARRAY(vy);

    // Adapted directly from SUNMatMatvec_Band
    for (i=0; i < a->M; i++)
	yd[i] = 0.0;
    for(j=0; j < a->N; j++) {
	col_j = a->cols[j] + a->s_mu;
	is = SUNMAX(0, j - a->mu);
	ie = SUNMIN(a->M - 1, j + a->ml);
	for (i=is; i <= ie; i++)
	    yd[i] += col_j[i-j]*xd[j];
    }

    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION < 300
CAMLprim value c_matrix_band_wrap(DlsMat a)
{
    CAMLparam0();
    CAMLlocal3(vcptr, vcontent, vr);

    vcontent = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, a->data, a->N, a->ldim);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vcptr = caml_alloc_final(1, NULL, 1, 20);
    DLSMAT(vcptr) = a;

    vdims = caml_alloc_tuple(RECORD_MAT_BANDDIMENSIONS_SIZE);
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_N,   Val_long(a->N));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_MU,  Val_long(a->mu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_SMU, Val_long(a->s_mu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_ML,  Val_long(a->ml));

    vpayload = caml_alloc_tuple(RECORD_MAT_BANDDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DATA, vcontent);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DIMS, vdims);


    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR,  vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID,   Val_bool(1));

    CAMLreturn(vr);
}
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Matrix.Sparse
 */

/* In Sundials < 3.0.0, synchronization between the C values and the OCaml
 * payload is not guaranteed since the SparseAddMat and SparseAddIdentity
 * functions may reallocate the underlying storage (with no way of notifying
 * this library). The ml_matrix_sparse_rewrap() function can be called to
 * resynchronize the OCaml values cached in the payload with the underlying
 * C data structures.
 *
 * In Sundials >= 3.0.0, we override the SUNSparseFromDenseMatrix,
 * SUNSparseFromBandMatrix, SUNMatClone_Sparse, SUNMatScaleAddI_Sparse, and
 * SUNMatScaleAdd_Sparse functions so that they update the OCaml payload
 * when the C data structures are reallocated.
 */

static void finalize_mat_content_sparse(value vcptra)
{
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptra);

#if SUNDIALS_LIB_VERSION >= 300
    /* indexvals, indexptrs, and data are freed when the corresponding
       bigarrays are finalized */
    free(content);

#elif SUNDIALS_LIB_VERSION >= 270
    SparseDestroyMat(SLSMAT(va));

#else // SUNDIALS_LIB_VERSION < 270
    DestroySparseMat(SLSMAT(va));

#endif
}

static void zero_sparse(MAT_CONTENT_SPARSE_TYPE A)
{
    sundials_ml_index i, *indexvals, *indexptrs;
    realtype *data;

    data = A->data;
#if SUNDIALS_LIB_VERSION >= 270
    indexptrs = A->indexptrs;
    indexvals = A->indexvals;
#else
    indexptrs = A->rowvals;
    indexvals = A->colptrs;
#endif

    /* Perform operation */
    for (i=0; i < A->NNZ; i++) {
	data[i] = 0.0;
	indexvals[i] = 0;
    }
    for (i=0; i < A->NP; i++)
	indexptrs[i] = 0;
    indexptrs[A->NP] = 0;
}

static void matrix_sparse_create_vcptr(sundials_ml_index m, sundials_ml_index n,
				       sundials_ml_index nnz, int sformat,
				       value *pvdata, value *pvidxvals,
				       value *pvidxptrs, value *pvcptr)
{
    CAMLparam0();

#if SUNDIALS_LIB_VERSION >= 300
    SUNMatrixContent_Sparse content = NULL;
    sundials_ml_index np = (sformat == CSC_MAT) ? n : m;

    // columns first
    *pvdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, NULL, nnz);
    *pvidxvals = caml_ba_alloc_dims(BIGARRAY_INT, 1, NULL, nnz);
    *pvidxptrs = caml_ba_alloc_dims(BIGARRAY_INT, 1, NULL, np + 1);

    // Setup the C-side content
    content = (SUNMatrixContent_Sparse)
		malloc(sizeof *content);
    if (content == NULL) caml_raise_out_of_memory();

    content->sparsetype = sformat;
    content->M = m;
    content->N = n;
    content->NNZ = nnz;
    content->NP = np;
    switch(sformat){
    case CSC_MAT:
	content->rowvals = &(content->indexvals);
	content->colptrs = &(content->indexptrs);
	/* CSR indices */
	content->colvals = NULL;
	content->rowptrs = NULL;
	break;
    case CSR_MAT:
	content->colvals = &(content->indexvals);
	content->rowptrs = &(content->indexptrs);
	/* CSC indices */
	content->rowvals = NULL;
	content->colptrs = NULL;
	break;
    }
    content->data = REAL_ARRAY(*pvdata);
    content->indexvals = INDEX_ARRAY(*pvidxvals);
    content->indexptrs = INDEX_ARRAY(*pvidxptrs);
    zero_sparse(content); // reproduce effect of callocs in Sundials code

    // Setup the OCaml-side
    *pvcptr = caml_alloc_final(1, &finalize_mat_content_sparse, 1, 20);
    MAT_CONTENT_SPARSE(*pvcptr) = content;

#else // SUNDIALS_LIB_VERSION < 300 (As per c_sparsematrix_new_sparse_mat)

#if SUNDIALS_LIB_VERSION >= 270
    SlsMat a = SparseNewMat(m, n, nnz, sformat);
#else
    SlsMat a = NewSparseMat(m, n, nnz);
#endif
    if (a == NULL) caml_raise_out_of_memory();

#if SUNDIALS_LIB_VERSION >= 270
    *pvidxvals = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->indexvals, a->NNZ);
    *pvidxptrs = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->indexptrs, a->NP + 1);
#else
    *pvidxvals = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->rowvals, a->NNZ);
    *pvidxptrs = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->colptrs, a->N + 1);
#endif
    *pvdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, a->data, a->NNZ);

    *pvcptr = caml_alloc_final(1, finalize_mat_content_sparse, 1, 20);
    DLSMAT(*pvcptr) = a;

#endif

    CAMLreturn0;
}

static value matrix_sparse_create_mat(sundials_ml_index m, sundials_ml_index n,
	sundials_ml_index nnz, int sformat)
{
    CAMLparam0();
    CAMLlocal5(vdata, vidxvals, vidxptrs, vcptr, vpayload);
    CAMLlocal1(vr);

    matrix_sparse_create_vcptr(m, n, nnz, sformat,
			       &vdata, &vidxvals, &vidxptrs, &vcptr);

    vpayload = caml_alloc_tuple(RECORD_MAT_SPARSEDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXVALS, vidxvals);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXPTRS, vidxptrs);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_DATA,    vdata);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_SFORMAT,
	    MAT_TO_SFORMAT(sformat));

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR,  vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID,   Val_bool(1));

    CAMLreturn(vr);
}

CAMLprim value ml_matrix_sparse_create(value vm, value vn, value vnnz,
				      value vsformat)
{
    CAMLparam4(vm, vn, vnnz, vsformat);
    CAMLlocal1(vr);

    vr = matrix_sparse_create_mat(Long_val(vm), Long_val(vn), Long_val(vnnz),
				  MAT_FROM_SFORMAT(vsformat));

    CAMLreturn(vr);
}

// Adapted directly from SUNSparseFromDenseMatrix
CAMLprim value ml_matrix_sparse_from_dense(value vsformat, value vcptrad,
					  value vdroptol)
{
    CAMLparam3(vsformat, vcptrad, vdroptol);
    CAMLlocal2(vr, vcptr);
    int sformat = MAT_FROM_SFORMAT(vsformat);
    realtype droptol = Double_val(vdroptol);
    realtype *As_data;

    sundials_ml_index i, j, nnz, *As_indexvals, *As_indexptrs;
    sundials_ml_index M, N;
    MAT_CONTENT_DENSE_TYPE Ad = MAT_CONTENT_DENSE(vcptrad);
    MAT_CONTENT_SPARSE_TYPE As;

    /* set size of new matrix */
    M = Ad->M;
    N = Ad->N;

    /* determine total number of nonzeros */
    nnz = 0;
    for (j=0; j < N; j++)
	for (i=0; i < M; i++)
	    nnz += SUNRabs(Ad->cols[j][i]) > droptol;

    /* allocate sparse matrix */
    vr = matrix_sparse_create_mat(M, N, nnz, sformat);
    vcptr = Field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    As = MAT_CONTENT_SPARSE(vcptr);

    As_data = As->data;
#if SUNDIALS_LIB_VERSION >= 270
    As_indexptrs = As->indexptrs;
    As_indexvals = As->indexvals;
#else
    As_indexptrs = As->rowvals;
    As_indexvals = As->colptrs;
#endif

    /* copy nonzeros from Ad into As, based on CSR/CSC type */
    nnz = 0;
    if (sformat == CSC_MAT) {
	for (j=0; j < N; j++) {
	    As_indexptrs[j] = nnz;
	    for (i=0; i<M; i++) {
		if ( SUNRabs(Ad->cols[j][i]) > droptol ) {
		    As_indexvals[nnz] = i;
		    As_data[nnz++] = Ad->cols[j][i];
		}
	    }
	}
	As_indexptrs[N] = nnz;
    } else {       /* CSR_MAT */
	for (i=0; i < M; i++) {
	    As_indexptrs[i] = nnz;
	    for (j=0; j < N; j++) {
		if ( SUNRabs(Ad->cols[j][i]) > droptol ) {
		    As_indexvals[nnz] = j;
		    As_data[nnz++] = Ad->cols[j][i];
		}
	    }
	}
	As_indexptrs[M] = nnz;
    }

    CAMLreturn(vr);
}

// Adapted directly from SUNSparseFromBandMatrix
CAMLprim value ml_matrix_sparse_from_band(value vsformat, value vcptrab,
					 value vdroptol)
{
    CAMLparam3(vsformat, vcptrab, vdroptol);
    CAMLlocal2(vr, vcptr);
    int sformat = MAT_FROM_SFORMAT(vsformat);
    realtype droptol = Double_val(vdroptol);
    realtype *As_data;

    sundials_ml_index i, j, nnz, *As_indexptrs, *As_indexvals;
    sundials_ml_index M, N;
    MAT_CONTENT_BAND_TYPE Ab = MAT_CONTENT_BAND(vcptrab);
    MAT_CONTENT_SPARSE_TYPE As;

    /* set size of new matrix */
    M = Ab->M;
    N = Ab->N;

    /* determine total number of nonzeros */
    nnz = 0;
    for (j=0; j < N; j++)
	for (i=SUNMAX(0, j - Ab->mu); i <= SUNMIN(M-1, j + Ab->ml); i++)
	    nnz += (SUNRabs(Ab->cols[j][i - j + Ab->s_mu]) > droptol);

    /* allocate sparse matrix */
    vr = matrix_sparse_create_mat(M, N, nnz, sformat);
    vcptr = Field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    As = MAT_CONTENT_SPARSE(vcptr);

    As_data = As->data;
#if SUNDIALS_LIB_VERSION >= 270
    As_indexptrs = As->indexptrs;
    As_indexvals = As->indexvals;
#else
    As_indexptrs = As->rowvals;
    As_indexvals = As->colptrs;
#endif

    /* copy nonzeros from Ab into As, based on CSR/CSC type */
    nnz = 0;
    if (sformat == CSC_MAT) {
	for (j=0; j < N; j++) {
	    As_indexptrs[j] = nnz;
	    for (i=SUNMAX(0, j - Ab->mu); i <= SUNMIN(M-1, j + Ab->ml); i++) {
		if ( SUNRabs(Ab->cols[j][i - j + Ab->s_mu]) > droptol ) {
		    As_indexvals[nnz] = i;
		    As_data[nnz++] = Ab->cols[j][i - j + Ab->s_mu];
		}
	    }
	}
	As_indexptrs[N] = nnz;
    } else {       /* CSR_MAT */
	for (i=0; i < M; i++) {
	    As_indexptrs[i] = nnz;
	    for (j=SUNMAX(0, i - Ab->ml); j <= SUNMIN(N-1, i + Ab->s_mu); j++) {
		if ( SUNRabs(Ab->cols[j][i - j + Ab->s_mu] ) > droptol ) {
		    As_indexvals[nnz] = j;
		    As_data[nnz++] = Ab->cols[j][i - j + Ab->s_mu];
		}
	    }
	}
	As_indexptrs[M] = nnz;
    }

    CAMLreturn(vr);
}

CAMLprim value ml_matrix_sparse_size(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_SPARSE_TYPE a = MAT_CONTENT_SPARSE(vcptr);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_long(a->M));
    Store_field(vr, 1, Val_long(a->N));

    CAMLreturn(vr);
}

CAMLprim value ml_matrix_sparse_dims(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_long(content->NNZ));
    Store_field(vr, 1, Val_long(content->NP));

    CAMLreturn(vr);
}

#if SUNDIALS_LIB_VERSION < 300
static value sparse_wrap_payload(SlsMat a)
{
    CAMLparam0;
    CAMLlocal4(vpayload, vidxvals, vidxptrs, vdata);
    int format;

#if SUNDIALS_LIB_VERSION >= 270
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->indexvals, a->NNZ);
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->indexptrs, a->NP + 1);
    format = a->sparsetype;
#else
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->rowvals, a->NNZ);
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INT, 1, a->colptrs, a->N + 1);
    format = CSC_MAT;
#endif
    vdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, content->data, content->NNZ);

    vpayload = caml_alloc_tuple(RECORD_MAT_SPARSEDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXVALS, vidxvals);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXPTRS, vidxptrs);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_DATA,    vdata);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_SFORMAT, MAT_TO_SFORMAT(format));

    CAMLreturn(vpayload);
}
#endif

static void matrix_sparse_resize(value va, sundials_ml_index nnz,
	int copy, int free)
{
    CAMLparam1(va);
    CAMLlocal5(vcptr, vdata, vidxvals, vpayload, vnewpayload);
    CAMLlocal1(vidxptrs);
    sundials_ml_index old_nnz, i;
    sundials_ml_index *old_indexvals, *new_indexvals;
    sundials_ml_index *old_indexptrs, *new_indexptrs;
    realtype *old_data, *new_data;
    MAT_CONTENT_SPARSE_TYPE A;

    vcptr    = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    // holding vpayload ensures that the underlying arrays are not gc-ed.
    vpayload = Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD);
    A = MAT_CONTENT_SPARSE(vcptr);

    if (nnz <= 0)
	nnz = A->NNZ;

#if SUNDIALS_LIB_VERSION >= 270
    old_indexptrs = A->indexptrs;
    old_indexvals = A->indexvals;
#else
    old_indexptrs = A->rowvals;
    old_indexvals = A->colptrs;
#endif
    old_data = A->data;
    old_nnz = A->NNZ;

#if SUNDIALS_LIB_VERSION >= 300
    vdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, NULL, nnz);
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INT, 1, NULL, nnz);
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INT, 1, NULL, A->NP + 1);

    /* Strictly speaking, it is not necessary to reallocate and recopy the
       indexptrs. We do it for two reasons:
       (1) it is conceptually simpler to deep-copy the sparse payload
       (2) it simplifies functions like ml_matrix_sparse_scale_add,
	   since they can "pretend" that resize duplicates the
	   matrix argument if they cache the old underlying arrays
	   and call with free=0. */

    vnewpayload = caml_alloc_tuple(RECORD_MAT_SPARSEDATA_SIZE);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_IDXVALS, vidxvals);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_IDXPTRS, vidxptrs);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_DATA,    vdata);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_SFORMAT,
		MAT_TO_SFORMAT(A->sparsetype));

    Store_field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vnewpayload);

    new_data = REAL_ARRAY(vdata);
    new_indexptrs = INDEX_ARRAY(vidxptrs);
    new_indexvals = INDEX_ARRAY(vidxvals);

    A->data = new_data;
    A->indexptrs = new_indexptrs;
    A->indexvals = new_indexvals;
#else
    new_indexvals = (int *) malloc(nnz * sizeof(int));
    if (new_indexvals == NULL) caml_raise_out_of_memory();

    new_indexptrs = (int *) malloc((A->NP + 1) * sizeof(int));
    if (new_indexptrs == NULL) {
	free(new_indexvals);
	caml_raise_out_of_memory();
    }

    new_data = (realtype *) malloc(nnz * sizeof(realtype));
    if (new_data == NULL) {
	free(new_indexptrs);
	free(new_indexvals);
	caml_raise_out_of_memory();
    }

#if SUNDIALS_LIB_VERSION >= 270
    A->indexptrs = new_indexptrs;
    A->indexvals = new_indexvals;
#else
    A->rowvals = new_indexptrs;
    A->colptrs = new_indexvals;
#endif
    A->data = new_data;

    vnewpayload = sparse_wrap_payload(A);
    Store_field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vnewpayload);
#endif
    A->NNZ = nnz;

    if (copy) {
	for (i=0; i < SUNMIN(old_nnz, nnz); i++) {
	    new_data[i] = old_data[i];
	    new_indexvals[i] = old_indexvals[i];
	}
	for (; i < nnz; i++) {
	    new_data[i] = 0.0;
	    new_indexvals[i] = 0;
	}
	for (i=0; i < A->NP + 1; i++) {
	    new_indexptrs[i] = old_indexptrs[i];
	}
    }

#if SUNDIALS_LIB_VERSION < 300
    if (free) {
	free(old_data);
	free(old_indexptrs);
	free(old_indexvals);
	// For >= 300, the OCaml bigarrays ensure gc of these arrays
    }
#endif

    CAMLreturn0;
}


// Adapted directly from SUNMatScaleAdd_Sparse
CAMLprim void ml_matrix_sparse_scale_add(value vc, value va, value vcptrb)
{
    CAMLparam3(vc, va, vcptrb);
    CAMLlocal5(vcptra, vdatac, vidxvalsc, vidxptrsc, vcptrc);
    CAMLlocal1(vpayload);
    realtype c = Double_val(vc);
    sundials_ml_index j, i, p, nz;
    bool newmat;
    sundials_ml_index *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
    realtype *x, *Ax, *Bx, *Cx;
    sundials_ml_index M, N;
    sundials_ml_index *A_indexptrs, *A_indexvals, *B_indexptrs, *B_indexvals;
    MAT_CONTENT_SPARSE_TYPE A, B;

    vcptra = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    A = MAT_CONTENT_SPARSE(vcptra);
    B = MAT_CONTENT_SPARSE(vcptrb);

    /* Perform operation */

    /* if A is CSR matrix, transpose M and N */
    if (A->sparsetype == CSC_MAT) {
	M = A->M;
	N = A->N;
    } else {
	M = A->N;
	N = A->M;
    }

    /* create work arrays for row indices and nonzero column values */
    w = (sundials_ml_index *) malloc(M * sizeof(sundials_ml_index));
    x = (realtype *) malloc(M * sizeof(realtype));

#if SUNDIALS_LIB_VERSION >= 270
    A_indexptrs = A->indexptrs;
    A_indexvals = A->indexvals;
    B_indexptrs = B->indexptrs;
    B_indexvals = B->indexvals;
#else
    A_indexptrs = A->rowvals;
    A_indexvals = A->colptrs;
    B_indexptrs = B->rowvals;
    B_indexvals = B->colptrs;
#endif

    /* determine if A already contains the sparsity pattern of B,
       and calculate the number of non-zeroes required if we create
       a new matrix (instead of reallocating). */
    newmat = 0;
    nz = 0;
    for (j=0; j < N; j++) {

	/* clear work array */
	for (i=0; i < M; i++)
	    w[i] = 0;

	/* scan column of A, incrementing w by one */
	for (i=A_indexptrs[j]; i < A_indexptrs[j+1]; i++) {
	    w[A_indexvals[i]] += 1;
	    nz += 1;
	}

	/* scan column of B, decrementing w by one */
	for (i=B_indexptrs[j]; i < B_indexptrs[j+1]; i++) {
	    if (w[B_indexvals[i]] == 0) nz += 1;
	    w[B_indexvals[i]] -= 1;
	}

	/* if any entry of w is negative, A doesn't contain B's sparsity */
	if (!newmat) {
	    for (i=0; i < M; i++)
		if (w[i] < 0) {
		    newmat = 1;
		    break;
		}
	}
    }

    /* perform operation */

    /*   case 1: A already contains sparsity pattern of B */
    if (!newmat) {

	/* iterate through columns, adding matrices */
	for (j=0; j < N; j++) {

	    /* clear work array */
	    for (i=0; i < M; i++)
		x[i] = 0.0;

	    /* scan column of B, updating work array */
	    for (i = B_indexptrs[j]; i < B_indexptrs[j+1]; i++)
		x[B_indexvals[i]] = B->data[i];

	    /* scan column of A, updating entries appropriately array */
	    for (i = A_indexptrs[j]; i < A_indexptrs[j+1]; i++)
		A->data[i] = c*A->data[i] + x[A_indexvals[i]];
	}

    /*   case 2: A does not already contain B's sparsity */
    } else {

	/* access data from (pre resize) CSR structures (return if failure) */
	Ap = A_indexptrs;
	Ai = A_indexvals;
	Ax = A->data;
	Bp = B_indexptrs;
	Bi = B_indexvals;
	Bx = B->data;

	/* reallocate memory within A */
	matrix_sparse_resize(va, nz, 0, 0); // no-copy, no-free
	zero_sparse(A);

	/* access data from (post resize) CSR structures (return if failure) */
	Cp = A_indexptrs;
	Ci = A_indexvals;
	Cx = A->data;

	/* initialize total nonzero count */
	nz = 0;

	/* iterate through columns */
	for (j=0; j < N; j++) {

	    /* set current column pointer to current # nonzeros */
	    Cp[j] = nz;

	    /* clear out temporary arrays for this column */
	    for (i=0; i<M; i++) {
		w[i] = 0;
		x[i] = 0.0;
	    }

	    /* iterate down column of A, collecting nonzeros */
	    for (p=Ap[j]; p < Ap[j+1]; p++) {
		w[Ai[p]] += 1;         /* indicate that row is filled */
		x[Ai[p]] = c*Ax[p];    /* collect/scale value */
	    }

	    /* iterate down column of B, collecting nonzeros */
	    for (p=Bp[j]; p<Bp[j+1]; p++) {
		w[Bi[p]] += 1;       /* indicate that row is filled */
		x[Bi[p]] += Bx[p];   /* collect value */
	    }

	    /* fill entries of C with this column's data */
	    for (i=0; i < M; i++) {
		if ( w[i] > 0 ) {
		    Ci[nz] = i;
		    Cx[nz++] = x[i];
		}
	    }
	}

	/* indicate end of data */
	Cp[N] = nz;

#if SUNDIALS_LIB_VERSION < 300
	// free the old underlying arrays (for >= 300, the gc does it)
	free(Ap);
	free(Ai);
	free(Ax);
#endif
    }

    /* clean up */
    free(w);
    free(x);

    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
static int csmat_sparse_scale_add(realtype c, SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal3(vcontenta, vcontentb, r);
    static value * closure_f = NULL;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = MAT_BACKLINK(B);

    if (closure_f == NULL)
	closure_f = caml_named_value("ml_matrix_sparse_scale_add");

    r = caml_callback3_exn(*closure_f,
			   caml_copy_double(c),
			   vcontenta,
			   Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR));
    if (Is_exception_result(r))
	sundials_ml_warn_discarded_exn(Extract_exception (r), "Sparse.scale_add");

    CAMLreturnT(int, Is_exception_result(r));
}
#endif

// Adapted directly from SUNMatScaleAddI_Sparse
CAMLprim void ml_matrix_sparse_scale_addi(value vc, value va)
{
    CAMLparam2(vc, va);
    CAMLlocal2(vcptr, vpayload);
    CAMLlocal4(vdatac, vidxvalsc, vidxptrsc, vcptrc);
    realtype c = Double_val(vc);
    sundials_ml_index j, i, p, nz;
    bool newmat, found;
    sundials_ml_index *w, *Ap, *Ai, *Cp, *Ci;
    realtype *x, *Ax, *Cx;
    sundials_ml_index M, N, *A_indexptrs, *A_indexvals;
    MAT_CONTENT_SPARSE_TYPE A;

    vcptr = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    A = MAT_CONTENT_SPARSE(vcptr);

#if SUNDIALS_LIB_VERSION >= 270
    A_indexptrs = A->indexptrs;
    A_indexvals = A->indexvals;
#else
    A_indexptrs = A->rowvals;
    A_indexvals = A->colptrs;
#endif

    /* Perform operation */

    if (A->sparsetype == CSC_MAT) {
	M = A->M;
	N = A->N;
    } else {
	M = A->N;
	N = A->M;
    }

    /* determine if A already contains values on the diagonal (hence
       no memory allocation necessary), and calculate the number of non-zeroes
       required if we create a new matrix (instead of reallocating). */
    newmat = 0;
    nz = SUNMIN(A->N, A->M);
    for (j=0; j < N; j++) {
	/* scan column (row if CSR) of A, searching for diagonal value */
	found = 0;
	for (i=A_indexptrs[j]; i < A_indexptrs[j+1]; i++) {
	    if (A_indexvals[i] == j) {
		found = 1;
	    } else {
		nz += 1;
	    }
	}
	/* if no diagonal found, signal new matrix */
	if (!found) newmat = 1;
    }

    /* perform operation */

    /*   case 1: A already contains a diagonal */
    if (!newmat) {
	/* iterate through columns, adding 1.0 to diagonal */
	for (j=0; j < SUNMIN(A->N, A->M); j++)
	    for (i=A_indexptrs[j]; i < A_indexptrs[j+1]; i++)
		if (A_indexvals[i] == j) {
		    A->data[i] = 1.0 + c * A->data[i];
		} else {
		    A->data[i] = c * A->data[i];
		}

    /*   case 2: A does not already contain a diagonal */
    } else {
	/* create work arrays for row indices and nonzero column values */
	w = (sundials_ml_index *) malloc(A->M * sizeof(sundials_ml_index));
	x = (realtype *) malloc(A->M * sizeof(realtype));

	/* access data from (pre resize) CSR structures (return if failure) */
	Ap = A_indexptrs;
	Ai = A_indexvals;
	Ax = A->data;

	/* reallocate memory within A */
	matrix_sparse_resize(va, nz, 0, 0); // no-copy, no-free
	zero_sparse(A);

	/* access data from (post resize) CSR structures (return if failure) */
	Cp = A_indexptrs;
	Ci = A_indexvals;
	Cx = A->data;

	/* initialize total nonzero count */
	nz = 0;

	/* iterate through columns (rows for CSR) */
	for (j=0; j < N; j++) {

	    /* set current column (row) pointer to current # nonzeros */
	    Cp[j] = nz;

	    /* clear out temporary arrays for this column (row) */
	    for (i=0; i < M; i++) {
		w[i] = 0;
		x[i] = 0.0;
	    }

	    /* iterate down column (along row) of A, collecting nonzeros */
	    for (p=Ap[j]; p < Ap[j+1]; p++) {
		w[Ai[p]] += 1;         /* indicate that row is filled */
		x[Ai[p]] = c*Ax[p];    /* collect/scale value */
	    }

	    /* add identity to this column (row) */
	    if (j < M) {
		w[j] += 1;     /* indicate that row is filled */
		x[j] += 1.0;  /* update value */
	    }

	    /* fill entries of C with this column's (row's) data */
	    for (i=0; i < M; i++) {
		if ( w[i] > 0 ) {
		    Ci[nz] = i;
		    Cx[nz++] = x[i];
		}
	    }
	}

	/* indicate end of data */
	Cp[N] = nz;

#if SUNDIALS_LIB_VERSION < 300
	// free the old underlying arrays (for >= 300, the gc does it)
	free(Ap);
	free(Ai);
	free(Ax);
#endif

	/* clean up */
	free(w);
	free(x);
    }

    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
static int csmat_sparse_scale_addi(realtype c, SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(vcontenta, r);
    static value * closure_f = NULL;

    vcontenta = MAT_BACKLINK(A);

    if (closure_f == NULL)
	closure_f = caml_named_value("ml_matrix_sparse_scale_addi");

    r = caml_callback2_exn(*closure_f, caml_copy_double(c), vcontenta);
    if (Is_exception_result(r))
	sundials_ml_warn_discarded_exn(Extract_exception (r),
				       "Sparse.scale_addi");

    CAMLreturnT(int, Is_exception_result(r));
}
#endif

// Adapted directly from SUNMatMatvec_Sparse, Matvec_SparseCSC, and
// Matvec_SparseCSR
CAMLprim void ml_matrix_sparse_matvec(value vcptra, value vx, value vy)
{
    CAMLparam3(vcptra, vx, vy);
    MAT_CONTENT_SPARSE_TYPE A = MAT_CONTENT_SPARSE(vcptra);
    sundials_ml_index i, j;
    sundials_ml_index *Ap, *Ai;
    realtype *Ax, *xd, *yd;

    xd = REAL_ARRAY(vx);
    yd = REAL_ARRAY(vy);

#if SUNDIALS_LIB_VERSION >= 270
    Ap = A->indexptrs;
    Ai = A->indexvals;
#else
    Ap = A->rowvals;
    Ai = A->colptrs;
#endif
    Ax = A->data;

    /* initialize result */
    for (i=0; i < A->M; i++)
	yd[i] = 0.0;

    if (A->sparsetype == CSC_MAT) {
	/* iterate through matrix columns */
	for (j=0; j < A->N; j++) {
	    /* iterate down column of A, performing product */
	    for (i=Ap[j]; i < Ap[j+1]; i++)
		yd[Ai[i]] += Ax[i]*xd[j];
	}

    } else {
	/* iterate through matrix rows */
	for (i=0; i < A->M; i++) {
	    /* iterate along row of A, performing product */
	    for (j=Ap[i]; j < Ap[i+1]; j++)
		yd[i] += Ax[j]*xd[Ai[j]];
	}
    }

    CAMLreturn0;
}

CAMLprim void ml_matrix_sparse_resize(value va, value vnnz, value vcopy)
{
    CAMLparam3(va, vnnz, vcopy);
    matrix_sparse_resize(va, Long_val(vnnz), Bool_val(vcopy), 1);
    CAMLreturn0;
}

// Adapted directly from SUNMatCopy_Sparse
CAMLprim void ml_matrix_sparse_copy(value vcptra, value vb)
{
    CAMLparam2(vcptra, vb);
    CAMLlocal1(vcptrb);
    MAT_CONTENT_SPARSE_TYPE A = MAT_CONTENT_SPARSE(vcptra);
    MAT_CONTENT_SPARSE_TYPE B;
    sundials_ml_index i, A_nz, *A_indexptrs, *A_indexvals;
    sundials_ml_index *B_indexvals, *B_indexptrs;

    vcptrb = Field(vb, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    B = MAT_CONTENT_SPARSE(vcptrb);

#if SUNDIALS_LIB_VERSION >= 270
    A_indexptrs = A->indexptrs;
    A_indexvals = A->indexvals;
    B_indexptrs = B->indexptrs;
    B_indexvals = B->indexvals;
#else
    A_indexptrs = A->rowvals;
    A_indexvals = A->colptrs;
    B_indexptrs = B->rowvals;
    B_indexvals = B->colptrs;
#endif

    /* Perform operation */
    A_nz = A_indexptrs[A->NP];

    /* ensure that B is allocated with at least as
    much memory as we have nonzeros in A */
    if (B->NNZ < A_nz) matrix_sparse_resize(vb, A_nz, 0, 1); // no-copy, free

    /* zero out B so that copy works correctly */
    zero_sparse(B);

    /* copy the data and row indices over */
    for (i=0; i < A_nz; i++){
	B->data[i] = A->data[i];
	B_indexvals[i] = A_indexvals[i];
    }

    /* copy the column pointers over */
    for (i=0; i < A->NP; i++) {
	B_indexptrs[i] = A_indexptrs[i];
    }
    B_indexptrs[A->NP] = A_nz;

    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
static int csmat_sparse_copy(SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal3(vcontenta, vcontentb, r);
    static value * closure_f = NULL;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = MAT_BACKLINK(B);

    if (closure_f == NULL)
	closure_f = caml_named_value("ml_matrix_sparse_copy");

    r = caml_callback2_exn(*closure_f,
			   Field(vcontenta, RECORD_MAT_MATRIXCONTENT_RAWPTR),
			   vcontentb);
    if (Is_exception_result(r))
	sundials_ml_warn_discarded_exn(Extract_exception (r), "Sparse.blit");

    CAMLreturnT(int, Is_exception_result(r));
}
#endif

CAMLprim value ml_matrix_sparse_space(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);

    // Directly adapted from SUNMatSpace_Dense
    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_long(content->NNZ));
    Store_field(vr, 1, Val_long(10 + content->NP + content->NNZ));

    CAMLreturn(vr);
}

// Sundials < 3.0.0
CAMLprim void ml_matrix_sparse_set_idx(value vcptr, value vj, value vidx)
{
    CAMLparam3(vcptr, vj, vidx);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_index j = Long_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NP + 1) caml_invalid_argument("j");
#endif

#if SUNDIALS_LIB_VERSION >= 270
    content->indexptrs[j] = Long_val(vidx);
#else
    content->rowvals[j] = Long_val(vidx);
#endif

    CAMLreturn0;
}

// Sundials < 3.0.0
CAMLprim value ml_matrix_sparse_get_idx(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_index j = Long_val(vj);
    realtype r;

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NP + 1) caml_invalid_argument("j");
#endif

#if SUNDIALS_LIB_VERSION >= 270
    r = content->indexptrs[j];
#else
    r = content->rowvals[j];
#endif

    CAMLreturn(Val_long(r));
}

// Sundials < 3.0.0
CAMLprim void ml_matrix_sparse_set_data(value vcptr, value vj, value vv)
{
    CAMLparam3(vcptr, vj, vv);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_index j = Long_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    content->data[j] = Double_val(vv);

    CAMLreturn0;
}

// Sundials < 3.0.0
CAMLprim value ml_matrix_sparse_get_val(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_index j = Long_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    CAMLreturn(Val_long(content->indexvals[j]));
}

// Sundials < 3.0.0
CAMLprim void ml_matrix_sparse_set_val(value vcptr, value vj, value vv)
{
    CAMLparam3(vcptr, vj, vv);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_index j = Long_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    content->indexvals[j] = Long_val(vv);

    CAMLreturn0;
}

// Sundials < 3.0.0
CAMLprim value ml_matrix_sparse_get_data(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_index j = Long_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    CAMLreturn(caml_copy_double(content->data[j]));
}

// Sundials < 3.0.0
// Reconnects the OCaml payload to the underlying C data if necessary
CAMLprim value ml_matrix_sparse_rewrap(value vm)
{
    CAMLparam1(vm);
    CAMLlocal2(vcptr, vpayload);

#if SUNDIALS_LIB_VERSION >= 300
    caml_failwith("ml_matrix_sparse_rewrap should not be called!");
#else
    MAT_CONTENT_SPARSE_TYPE content;
    void *ba_data;
    sundials_ml_index *idxptrs, *idxvals;
    int format;

    vcptr = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    content = MAT_CONTENT_SPARSE(vcptr);

    vpayload = Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD);
    ba_data = Caml_ba_data_val(Field(vpayload, RECORD_MAT_SPARSEDATA_DATA));

    // Has the underlying data array changed?
    if (A->data <> ba_data)
	vpayload = sparse_wrap_payload(content);
#endif

    CAMLreturn(vpayload);
}

#if SUNDIALS_LIB_VERSION < 300
CAMLprim value c_matrix_sparse_wrap(SlsMat a)
{
    CAMLparam0();
    CAMLlocal3(vpayload, vcptr, vr);

    vpayload = sparse_wrap_payload(a);

    vcptr = caml_alloc_final(1, NULL, 1, 20);
    SLSMAT(vcptr) = a;

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR,  vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID,   Val_bool(1));

    CAMLreturn(vr);
}
#endif


// Sundials < 3.0.0
CAMLprim void ml_matrix_sparse_set_to_zero(value vcptr)
{
    CAMLparam1(vcptr);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);

    /* Perform operation */
    zero_sparse(content);

    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Matrix
 */

#if SUNDIALS_LIB_VERSION >= 300
/*
 *  Type	    content		    backlink
 *  --------+-------------------------------------------------
 *  dense   | MAT_CONTENT_DENSE_TYPE	    Dense.t  (matrix_content)
 *  band    | MAT_CONTENT_BAND_TYPE	    Band.t   (matrix_content)
 *  sparse  | MAT_CONTENT_SPARSE_TYPE	    Sparse.t (matrix_content)
 *  custom  | matrix_ops record		    'm (custom data)
 */
static SUNMatrix alloc_smat(void *content, value backlink,
			    bool content_is_value)
{
    SUNMatrix smat;

    /* Alloc memory in C heap */
    smat = (SUNMatrix)malloc(sizeof(struct csmat));
    if (smat == NULL) return NULL;

    smat->ops = (SUNMatrix_Ops)malloc(sizeof(struct _generic_SUNMatrix_Ops));
    if (smat->ops == NULL) { free(smat); return(NULL); }

    smat->content = content;
    if (content_is_value) // NB: free_smat does not unregister the content
	caml_register_generational_global_root((value *)&(smat->content));

    MAT_BACKLINK(smat) = backlink;
    caml_register_generational_global_root(&MAT_BACKLINK(smat));

    return smat;
}

#define MAT_OP_TABLE(smat)  ((smat)->content)
#define GET_OP(smat, x) (Field((value)MAT_OP_TABLE(smat), x))

static void free_smat(SUNMatrix smat)
{
    caml_remove_generational_global_root(&MAT_BACKLINK(smat));
    // smat->content is not owned by the sunmatrix
    free(smat->ops);
    free(smat);
}

static void free_custom_smat(SUNMatrix smat)
{
    caml_remove_generational_global_root((value *)&MAT_OP_TABLE(smat));
    free_smat(smat);
}

static void finalize_caml_smat(value vsmat)
{
    free_smat(MAT_CVAL(vsmat));
}

static void finalize_caml_custom_smat(value vsmat)
{
    free_custom_smat(MAT_CVAL(vsmat));
}

static void csmat_clone_ops(SUNMatrix dst, SUNMatrix src)
{
    SUNMatrix_Ops ops = (SUNMatrix_Ops)dst->ops;

    ops->getid     = src->ops->getid;
    ops->clone     = src->ops->clone;
    ops->destroy   = src->ops->destroy;
    ops->zero      = src->ops->zero;
    ops->copy      = src->ops->copy;
    ops->scaleadd  = src->ops->scaleadd;
    ops->scaleaddi = src->ops->scaleaddi;
    ops->matvec    = src->ops->matvec;
    ops->space     = src->ops->space;
}

static SUNMatrix csmat_dense_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(vcontenta, vcontentb);
    SUNMatrix B;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = ml_matrix_dense_create(Val_long(SM_ROWS_D(A)),
				       Val_long(SM_COLUMNS_D(A)));

    caml_ba_blit(Field(vcontenta, RECORD_MAT_MATRIXCONTENT_PAYLOAD),
		 Field(vcontentb, RECORD_MAT_MATRIXCONTENT_PAYLOAD));

    B = alloc_smat(
	    MAT_CONTENT(Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR)),
	    vcontentb, false);
    csmat_clone_ops(B, A);

    CAMLreturnT(SUNMatrix, B);
}

static SUNMatrix csmat_band_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal4(vcontenta, vpayloada, vcontentb, vpayloadb);
    SUNMatrix B;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = ml_matrix_band_create_mat(SM_COLUMNS_B(A), SM_UBAND_B(A),
				          SM_LBAND_B(A), SM_SUBAND_B(A));

    vpayloada = Field(vcontenta, RECORD_MAT_MATRIXCONTENT_PAYLOAD);
    vpayloadb = Field(vcontentb, RECORD_MAT_MATRIXCONTENT_PAYLOAD);

    caml_ba_blit(Field(vpayloada, RECORD_MAT_BANDDATA_DATA),
		 Field(vpayloadb, RECORD_MAT_BANDDATA_DATA));

    B = alloc_smat(
	    MAT_CONTENT(Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR)),
	    vcontentb, false);
    csmat_clone_ops(B, A);

    CAMLreturnT(SUNMatrix, B);
}

static SUNMatrix csmat_sparse_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal4(vcontenta, vpayloada, vcontentb, vpayloadb);
    SUNMatrix B;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = matrix_sparse_create_mat(SM_ROWS_S(A), SM_COLUMNS_S(A),
					 SM_NNZ_S(A), SM_SPARSETYPE_S(A));

    vpayloada = Field(vcontenta, RECORD_MAT_MATRIXCONTENT_PAYLOAD);
    vpayloadb = Field(vcontentb, RECORD_MAT_MATRIXCONTENT_PAYLOAD);

    caml_ba_blit(Field(vpayloada, RECORD_MAT_SPARSEDATA_IDXVALS),
		 Field(vpayloadb, RECORD_MAT_SPARSEDATA_IDXVALS));
    caml_ba_blit(Field(vpayloada, RECORD_MAT_SPARSEDATA_IDXPTRS),
		 Field(vpayloadb, RECORD_MAT_SPARSEDATA_IDXPTRS));
    caml_ba_blit(Field(vpayloada, RECORD_MAT_SPARSEDATA_DATA),
		 Field(vpayloadb, RECORD_MAT_SPARSEDATA_DATA));

    B = alloc_smat(
	    MAT_CONTENT(Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR)),
	    vcontentb, false);
    csmat_clone_ops(B, A);

    CAMLreturnT(SUNMatrix, B);
}

static SUNMatrix csmat_custom_clone(SUNMatrix A);
static SUNMatrix_ID csmat_custom_getid(SUNMatrix A);
static int csmat_custom_zero(SUNMatrix A);
static int csmat_custom_copy(SUNMatrix A, SUNMatrix B);
static int csmat_custom_scale_add(realtype c, SUNMatrix A, SUNMatrix B);
static int csmat_custom_scale_addi(realtype c, SUNMatrix A);
static int csmat_custom_matvec(SUNMatrix A, N_Vector x, N_Vector y);
static int csmat_custom_space(SUNMatrix A, long int *lenrw, long int *leniw);

#endif

CAMLprim value ml_matrix_wrap(value vid, value vcontent, value vpayload)
{
    CAMLparam3(vid, vcontent, vpayload);
    CAMLlocal1(vr);

#if SUNDIALS_LIB_VERSION >= 300
    int mat_id = Int_val(vid);
    SUNMatrix smat;

    /* Create matrix */
    if (mat_id == MATRIX_ID_CUSTOM) {
	smat = alloc_smat((void *)vcontent, vpayload, true);
    } else {
	smat = alloc_smat(MAT_CONTENT(vcontent), vpayload, false);
    }
    if (smat == NULL) caml_raise_out_of_memory();

    /* Attach operations */
    switch (mat_id) {
    case MATRIX_ID_DENSE:
	smat->ops->clone       = csmat_dense_clone;       // ours
	smat->ops->destroy     = free_smat;  // ours (only called for c clones)
	smat->ops->getid       = SUNMatGetID_Dense;
	smat->ops->zero        = SUNMatZero_Dense;
	smat->ops->copy        = SUNMatCopy_Dense;
	smat->ops->scaleadd    = SUNMatScaleAdd_Dense;
	smat->ops->scaleaddi   = SUNMatScaleAddI_Dense;
	smat->ops->matvec      = SUNMatMatvec_Dense;
	smat->ops->space       = SUNMatSpace_Dense;
	break;

    case MATRIX_ID_BAND:
	smat->ops->clone       = csmat_band_clone;        // ours
	smat->ops->destroy     = free_smat;  // ours (only called for c clones)
	smat->ops->getid       = SUNMatGetID_Band;
	smat->ops->zero        = SUNMatZero_Band;
	smat->ops->copy        = csmat_band_copy;         // ours
	smat->ops->scaleadd    = csmat_band_scale_add;    // ours
	smat->ops->scaleaddi   = SUNMatScaleAddI_Band;
	smat->ops->matvec      = SUNMatMatvec_Band;
	smat->ops->space       = SUNMatSpace_Band;
	break;

    case MATRIX_ID_SPARSE:
	smat->ops->clone       = csmat_sparse_clone;      // ours
	smat->ops->destroy     = free_smat;  // ours (only called for c clones)
	smat->ops->getid       = SUNMatGetID_Sparse;
	smat->ops->zero        = SUNMatZero_Sparse;
	smat->ops->copy        = csmat_sparse_copy;	 // ours
	smat->ops->scaleadd    = csmat_sparse_scale_add;  // ours
	smat->ops->scaleaddi   = csmat_sparse_scale_addi; // ours
	smat->ops->matvec      = SUNMatMatvec_Sparse;
	smat->ops->space       = SUNMatSpace_Sparse;
	break;

    case MATRIX_ID_CUSTOM:
	smat->ops->clone       = csmat_custom_clone;
	smat->ops->destroy     = free_custom_smat;  // ours (only called for
						   //	    c clones)
	smat->ops->getid       = csmat_custom_getid;
	smat->ops->zero        = csmat_custom_zero;
	smat->ops->copy        = csmat_custom_copy;
	smat->ops->scaleadd    = csmat_custom_scale_add;
	smat->ops->scaleaddi   = csmat_custom_scale_addi;
	smat->ops->matvec      = csmat_custom_matvec;
	smat->ops->space       = csmat_custom_space;
	break;
    }

    // Setup the OCaml-side
    vr = caml_alloc_final(1,
	    (mat_id == MATRIX_ID_CUSTOM)
		? &finalize_caml_custom_smat
		: &finalize_caml_smat,
	    1, 20);
    MAT_CVAL(vr) = smat;

#else // SUNDIALS_LIB_VERSION < 300
    vr = Val_unit;
#endif
    CAMLreturn(vr);
}

#if SUNDIALS_LIB_VERSION >= 300

static SUNMatrix csmat_custom_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(mlop, vcontentb);
    SUNMatrix B;

    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_CLONE);

    vcontentb = caml_callback_exn(mlop, MAT_BACKLINK(A));
    if (Is_exception_result (vcontentb)) {
	sundials_ml_warn_discarded_exn (Extract_exception (vcontentb),
				    "user-defined matrix operation m_clone");
	fputs ("Sundials/ML has no sensible value to return to Sundials, "
	       "and incorrect values risk memory corruption.  Abort.", stderr);
	fflush (stderr);
	abort ();
    }

    B = alloc_smat(MAT_OP_TABLE(A), vcontentb, true);
    csmat_clone_ops(B, A);

    CAMLreturnT(SUNMatrix, B);
}

static SUNMatrix_ID csmat_custom_getid(SUNMatrix A)
{
    return SUNMATRIX_CUSTOM;
}

static int csmat_custom_zero(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_ZERO);

    r = caml_callback_exn(mlop, MAT_BACKLINK(A));
    if (Is_exception_result (r))
	sundials_ml_warn_discarded_exn (Extract_exception (r),
				    "user-defined matrix operation m_zero");

    CAMLreturnT(int, Is_exception_result (r));
}

static int csmat_custom_copy(SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_COPY);

    r = caml_callback2_exn(mlop, MAT_BACKLINK(A), MAT_BACKLINK(B));
    if (Is_exception_result (r))
	sundials_ml_warn_discarded_exn (Extract_exception (r),
				    "user-defined matrix operation m_copy");

    CAMLreturnT(int, Is_exception_result (r));
}

static int csmat_custom_scale_add(realtype c, SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_SCALE_ADD);

    r = caml_callback3_exn(mlop, caml_copy_double(c), MAT_BACKLINK(A),
			   MAT_BACKLINK(B));
    if (Is_exception_result (r))
	sundials_ml_warn_discarded_exn (Extract_exception (r),
				"user-defined matrix operation m_scale_add");

    CAMLreturnT(int, Is_exception_result (r));
}

static int csmat_custom_scale_addi(realtype c, SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_SCALE_ADDI);

    r = caml_callback2_exn(mlop, caml_copy_double(c), MAT_BACKLINK(A));
    if (Is_exception_result (r))
	sundials_ml_warn_discarded_exn (Extract_exception (r),
				"user-defined matrix operation m_scale_addi");

    CAMLreturnT(int, Is_exception_result (r));
}

static int csmat_custom_matvec(SUNMatrix A, N_Vector x, N_Vector y)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_MATVEC);

    r = caml_callback3_exn(mlop, MAT_BACKLINK(A), NVEC_BACKLINK(x),
			   NVEC_BACKLINK(y));
    if (Is_exception_result (r))
	sundials_ml_warn_discarded_exn (Extract_exception (r),
				"user-defined matrix operation m_matvec");

    CAMLreturnT(int, Is_exception_result (r));
}

static int csmat_custom_space(SUNMatrix A, long int *lenrw, long int *leniw)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_SPACE);

    r = caml_callback_exn(mlop, MAT_BACKLINK(A));
    if (Is_exception_result (r))
	sundials_ml_warn_discarded_exn (Extract_exception (r),
				"user-defined matrix operation m_space");

    CAMLreturnT(int, Is_exception_result (r));
}
#endif

CAMLprim void ml_matrix_scale_add(value vc, value va, value vb)
{
    CAMLparam3(vc, va, vb);
    if (SUNMatScaleAdd(Double_val(vc), MAT_VAL(va), MAT_VAL(vb)))
	caml_failwith("SUNMatScaleAdd");
    CAMLreturn0;
}

CAMLprim void ml_matrix_scale_addi(value vc, value va)
{
    CAMLparam2(vc, va);
    if (SUNMatScaleAddI(Double_val(vc), MAT_VAL(va)))
	caml_failwith("SUNMatScaleAddI");
    CAMLreturn0;
}

CAMLprim void ml_matrix_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);
    if (SUNMatMatvec(MAT_VAL(va), NVEC_VAL(vx), NVEC_VAL(vy)))
	caml_failwith("SUNMatMatvec");
    CAMLreturn0;
}

CAMLprim void ml_matrix_zero(value va)
{
    CAMLparam1(va);
    if (SUNMatZero(MAT_VAL(va)))
	caml_failwith("SUNMatZero");
    CAMLreturn0;
}

CAMLprim void ml_matrix_copy(value va, value vb)
{
    CAMLparam2(va, vb);
    if (SUNMatCopy(MAT_VAL(va), MAT_VAL(vb)))
	caml_failwith("SUNMatCopy");
    CAMLreturn0;
}

CAMLprim value ml_matrix_space(value va)
{
    CAMLparam1(va);
    CAMLlocal1(vr);
    long int lenrw, leniw;

    if (SUNMatSpace(MAT_VAL(va), &lenrw, &leniw))
	caml_failwith("SUNMatSpace");

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_long(lenrw));
    Store_field(vr, 1, Val_long(leniw));

    CAMLreturn(vr);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Array matrices
 */

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
	caml_invalid_argument("matrix not square.");
#endif

    denseAddIdentity(ARRAY2_ACOLS(va), m);
    CAMLreturn (Val_unit);
}

#if SUNDIALS_LIB_VERSION >= 260
CAMLprim value c_arraydensematrix_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);
    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vx) < n)
	caml_invalid_argument("x array too small.");
    if (ARRAY1_LEN(vy) < m)
	caml_invalid_argument("y array too small.");
#endif
    denseMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx), REAL_ARRAY(vy), m, n);
    CAMLreturn (Val_unit);
}
#else
CAMLprim value c_arraydensematrix_matvec(value va, value vx, value vy)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

CAMLprim value c_arraydensematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("pivot array too small.");
#endif

    int r = denseGETRF(ARRAY2_ACOLS(va), m, n, INDEX_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
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
	caml_invalid_argument("matrix not square.");
    if (ARRAY1_LEN(vb) < n)
	caml_invalid_argument("solution vector too small.");
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("pivot array too small.");
#endif

    denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp), REAL_ARRAY(vb));
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
	caml_invalid_argument("matrix not square.");
    if (ARRAY1_LEN(vb) - boff < n)
	caml_invalid_argument("b is too small.");
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("p is too small.");
#endif

    denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp), REAL_ARRAY(vb) + boff);
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
	caml_invalid_argument("matrix not square");
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
	caml_invalid_argument("matrix not square.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("b is too small.");
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
	caml_invalid_argument("fewer rows than columns.");
    if (ARRAY1_LEN(vbeta) < n)
	caml_invalid_argument("beta is too small.");
    if (ARRAY1_LEN(vv) < m)
	caml_invalid_argument("work is too small.");
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
	caml_invalid_argument("fewer rows than columns.");
    if (ARRAY1_LEN(Field(vormqr, 0)) < n)
	caml_invalid_argument("beta is too small.");
    if (ARRAY1_LEN(Field(vormqr, 1)) < n)
	caml_invalid_argument("v is too small.");
    if (ARRAY1_LEN(Field(vormqr, 2)) < m)
	caml_invalid_argument("w is too small.");
    if (ARRAY1_LEN(Field(vormqr, 3)) < m)
	caml_invalid_argument("work is too small.");
#endif

    denseORMQR(ARRAY2_ACOLS(va), m, n, beta, vv, vw, work);
    CAMLreturn (Val_unit);
}

/* Array Band matrix functions */

CAMLprim value c_arraybandmatrix_copy(value va, value vb, value vsizes)
{
    CAMLparam3(va, vb, vsizes);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat am = ba->dim[0];

    long int a_smu  = Long_val(Field(vsizes, 0));
    long int b_smu  = Long_val(Field(vsizes, 1));
    long int copymu = Long_val(Field(vsizes, 2));
    long int copyml = Long_val(Field(vsizes, 3));

#if SUNDIALS_ML_SAFE == 1
    intnat an = ba->dim[1];
    struct caml_ba_array *bb = ARRAY2_DATA(vb);

    intnat bm = bb->dim[0];
    intnat bn = bb->dim[1];

    if (an < copymu + copyml + 1)
	caml_invalid_argument("source matrix too small.");
    if (bn < copymu + copyml + 1)
	caml_invalid_argument("destination matrix too small.");
    if ((am != bm) || (bm != bn))
	caml_invalid_argument("matrix sizes differ.");
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
	caml_invalid_argument("matrix badly sized.");
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
	caml_invalid_argument("matrix badly sized.");
#endif

    bandAddIdentity(ARRAY2_ACOLS(va), m, smu);
    CAMLreturn (Val_unit);
}

#if SUNDIALS_LIB_VERSION >= 260
CAMLprim value c_arraybandmatrix_matvec(value va, value vsizes,
					value vx, value vy)
{
    CAMLparam4(va, vsizes, vx, vy);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("matrix badly sized.");

    if (ARRAY1_LEN(vx) < n)
	caml_invalid_argument("x array too small.");
    if (ARRAY1_LEN(vy) < m)
	caml_invalid_argument("y array too small.");
#endif
    bandMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx), REAL_ARRAY(vy),
	       m, mu, ml, smu);
    CAMLreturn (Val_unit);
}
#else
CAMLprim value c_arraybandmatrix_matvec(value va, value vx, value vy,
					value vsizes)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

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
	caml_invalid_argument("matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("p is too small.");
#endif

    bandGBTRF(ARRAY2_ACOLS(va), m, mu, ml, smu, INDEX_ARRAY(vp));
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
	caml_invalid_argument("matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("p is too small.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("b is too small.");
#endif

    bandGBTRS(ARRAY2_ACOLS(va), m, smu, ml, INDEX_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

