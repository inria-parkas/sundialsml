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

#elif SUNDIALS_LIB_VERSION >= 260 // 260 <= SUNDIALS_LIB_VERSION < 300
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>
#else // SUNDIALS_LIB_VERSION < 260
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#endif

#include <sundials/sundials_dense.h>
#include <sundials/sundials_band.h>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>

#include "../sundials/sundials_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

extern CAMLprim value caml_ba_blit(value vsrc, value vdst);
extern CAMLprim value caml_ba_fill(value vb, value vinit);

enum ocaml_values_index {
    IX_dense_ops = 0,
    IX_band_ops,
    IX_sparse_ops,
    NUM_OCAML_VALUES
};

static value ocaml_values[NUM_OCAML_VALUES];

CAMLprim void sunml_mat_init_module (value exns, value matrix_ops)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (MATRIX, exns);
    REGISTER_OCAML_VALUES (matrix_ops);
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

CAMLprim value sunml_matrix_dense_create(value vm, value vn)
{
    CAMLparam2(vm, vn);
    CAMLlocal3(vdata, vcptr, vr);
    sundials_ml_index j;
    MAT_CONTENT_DENSE_TYPE content = NULL;
    sundials_ml_index m = Index_val(vm);
    sundials_ml_index n = Index_val(vn);

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
    content->cols = (sunrealtype **) malloc(n * sizeof(sunrealtype *));
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
CAMLprim value sunml_matrix_dense_wrap(DlsMat a)
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
    sundials_ml_index i = Index_val(vi);
    sundials_ml_index j = Index_val(vj);
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
    sundials_ml_index i = Index_val(vi);
    sundials_ml_index j = Index_val(vj);
    MAT_CONTENT_DENSE_TYPE m = MAT_CONTENT_DENSE(vcptr);

#if SUNDIALS_ML_SAFE == 1
    if (i < 0 || i >= m->M) caml_invalid_argument("i");
    if (j < 0 || j >= m->N) caml_invalid_argument("j");
#endif

    m->cols[j][i] = Double_val(vv);
    CAMLreturn0;
}
*/

CAMLprim void sunml_matrix_dense_scale_add(value vc, value vcptra, value vcptrb)
{
    CAMLparam3(vc, vcptra, vcptrb);
    sunrealtype c = Double_val(vc);
    sundials_ml_index i, j;
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

CAMLprim void sunml_matrix_dense_scale_addi(value vc, value vcptra)
{
    CAMLparam2(vc, vcptra);
    sunrealtype c = Double_val(vc);
    sundials_ml_index i, j;
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

CAMLprim void sunml_matrix_dense_matvec(value vcptra, value vx, value vy)
{
    CAMLparam3(vcptra, vx, vy);
    sundials_ml_index i, j;
    sunrealtype *col_j, *xd, *yd;
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

CAMLprim value sunml_matrix_dense_space(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_DENSE_TYPE content = MAT_CONTENT_DENSE(vcptr);

    // Directly adapted from SUNMatSpace_Dense
    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_index(content->ldata));
    Store_field(vr, 1, Val_index(3 + content->N));

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

static bool matrix_band_create_vcptr(sundials_ml_index n,
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
    content->cols = (sunrealtype **) malloc(n * sizeof(sunrealtype *));
    if (content->cols == NULL) {
	free(content);
	return false;
    }
    for (j=0; j < n; j++) content->cols[j] = content->data + j * colSize;

#if SUNDIALS_LIB_VERSION < 300
    content->type = SUNDIALS_BAND;
#endif

    // Setup the OCaml-side
    *pvcptr = caml_alloc_final(1, &finalize_mat_content_band, 1, 20);
    MAT_CONTENT_BAND(*pvcptr) = content;

    return true;
}

// reallocate the storage underlying a matrix_content (Band.t)
static bool matrix_band_realloc(sundials_ml_index n, sundials_ml_index mu,
			        sundials_ml_index ml, sundials_ml_index smu,
			        value va, bool free_cols)
{
    CAMLparam1(va);
    CAMLlocal5(vpayload, vcptr, vnewdata, vnewdims, vnewpayload);
    MAT_CONTENT_BAND_TYPE content = NULL;
    sundials_ml_index colSize = smu + ml + 1;
    sundials_ml_index j;
    sunrealtype **newcols;

    vcptr    = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    vpayload = Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD);

    content = MAT_CONTENT_BAND(vcptr);

    // columns first
    vnewdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, NULL, n, colSize);
    caml_ba_fill(vnewdata, caml_copy_double(0.));

    newcols = (sunrealtype **) malloc(n * sizeof(sunrealtype *));
    if (newcols == NULL) CAMLreturnT(bool, false);

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
    // very slight chance of raising an exception and messing things up
    // no caml_alloc_shr_no_track_noexc in OCaml 4.03.0
    // https://caml.inria.fr/mantis/view.php?id=7572
    vnewdims = caml_alloc_tuple(RECORD_MAT_BANDDIMENSIONS_SIZE);
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_N,   Val_index(n));
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_MU,  Val_index(mu));
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_SMU, Val_index(smu));
    Store_field(vnewdims, RECORD_MAT_BANDDIMENSIONS_ML,  Val_index(ml));

    // ditto
    vnewpayload = caml_alloc_tuple(RECORD_MAT_BANDDATA_SIZE);
    Store_field(vnewpayload, RECORD_MAT_BANDDATA_DATA, vnewdata);
    Store_field(vnewpayload, RECORD_MAT_BANDDATA_DIMS, vnewdims);

    Store_field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vnewpayload);

    CAMLreturnT(bool, true);
}

#if 300 <= SUNDIALS_LIB_VERSION
static value sunml_matrix_band_create_mat(sundials_ml_index n,
	sundials_ml_index mu, sundials_ml_index ml, sundials_ml_index smu)
{
    CAMLparam0();
    CAMLlocal5(vdata, vcptr, vr, vdims, vpayload);

    if (!matrix_band_create_vcptr(n, mu, ml, smu, &vdata, &vcptr))
	caml_raise_out_of_memory();

    vdims = caml_alloc_tuple(RECORD_MAT_BANDDIMENSIONS_SIZE);
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_N,   Val_index(n));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_MU,  Val_index(mu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_SMU, Val_index(smu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_ML,  Val_index(ml));

    vpayload = caml_alloc_tuple(RECORD_MAT_BANDDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DATA, vdata);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DIMS, vdims);

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_RAWPTR, vcptr);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_VALID, Val_bool(1));

    CAMLreturn(vr);
}
#endif

CAMLprim value sunml_matrix_band_create(value vdims)
{
    CAMLparam1(vdims);
    CAMLlocal4(vdata, vcptr, vr, vpayload);

    sundials_ml_index n   =
	Index_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_N));
    sundials_ml_index mu  =
	Index_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_MU));
    sundials_ml_index ml  =
	Index_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_ML));
    sundials_ml_index smu =
	Index_val(Field(vdims, RECORD_MAT_BANDDIMENSIONS_SMU));

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
    sundials_ml_index i = Index_val(vi);
    sundials_ml_index j = Index_val(vj);
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
    sundials_ml_index i = Index_val(vi);
    sundials_ml_index j = Index_val(vj);
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
CAMLprim void sunml_matrix_band_copy(value vcptra, value vb)
{
    CAMLparam2(vcptra, vb);
    CAMLlocal3(vcptrb, vpayloadb, vdatab);
    sundials_ml_index i, j;
    sunrealtype *A_colj, *B_colj;
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
	if (! matrix_band_realloc(B->N, SUNMAX(B->mu, A->mu),
				        SUNMAX(B->ml, A->ml),
					SUNMAX(B->s_mu, A->s_mu), vb, 1) )
	    caml_raise_out_of_memory();

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
// Adapted directly from SUNMatCopy_Band
static int csmat_band_copy(SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal1(vcontentb);
    sundials_ml_index i, j;
    sunrealtype *A_colj, *B_colj;

    vcontentb = MAT_BACKLINK(B);

    if ((SM_ROWS_B(A) != SM_ROWS_B(B)) || (SM_COLUMNS_B(A) != SM_COLUMNS_B(B)))
	CAMLreturnT(int, 1);

    /* Grow B if A's bandwidth is larger */
    if ((SM_UBAND_B(A) > SM_UBAND_B(B)) || (SM_LBAND_B(A) > SM_LBAND_B(B)))
	if (! matrix_band_realloc(SM_COLUMNS_B(B),
				  SUNMAX(SM_UBAND_B(A),  SM_UBAND_B(B)),
				  SUNMAX(SM_LBAND_B(A),  SM_LBAND_B(B)),
				  SUNMAX(SM_SUBAND_B(A), SM_SUBAND_B(B)),
				  vcontentb, 1) )
	    CAMLreturnT(int, 1);

    /* Perform operation */
    if (SUNMatZero_Band(B) != 0)
	CAMLreturnT(int, 1);
    for (j=0; j<SM_COLUMNS_B(B); j++) {
	B_colj = SM_COLUMN_B(B,j);
	A_colj = SM_COLUMN_B(A,j);
	for (i=-SM_UBAND_B(A); i<=SM_LBAND_B(A); i++)
	    B_colj[i] = A_colj[i];
    }

    CAMLreturnT(int, 0);
}
#endif

// Adapted directly from SMScaleAddNew_Band
static int matrix_band_scale_add_new(sunrealtype c, value va, value vcptrb)
{
    CAMLparam2(va, vcptrb);
    CAMLlocal2(vcptra, voldpayload);
    sundials_ml_index i, j, A_N, A_ml, A_mu, A_s_mu, new_mu, new_ml;
    sunrealtype **A_cols;
    sunrealtype *A_colj, *B_colj, *C_colj;
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
    if (! matrix_band_realloc(B->N, new_mu, new_ml,
			      SUNMIN(A_N - 1, new_mu + new_ml), va, 0) )
	CAMLreturnT(int, 1); // failure

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

    CAMLreturnT(int, 0); // success
}

CAMLprim void sunml_matrix_band_scale_add(value vc, value va, value vcptrb)
{
    CAMLparam3(vc, va, vcptrb);
    CAMLlocal1(vcptra);
    sunrealtype c = Double_val(vc);
    sunrealtype *A_colj, *B_colj;
    MAT_CONTENT_BAND_TYPE A, B;
    sundials_ml_index i, j;

    vcptra = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    A = MAT_CONTENT_BAND(vcptra);
    B = MAT_CONTENT_BAND(vcptrb);

    if ( (B->mu > A->mu) || (B->ml > A->ml) ) {
	if (matrix_band_scale_add_new(c, va, vcptrb) != 0)
	    caml_raise_out_of_memory();

    } else {
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
// Adapted directly from SMScaleAddNew_Band
static int csmat_band_scale_add(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    sundials_ml_index i, j;
    sunrealtype *A_colj, *B_colj;

    /* Verify that A and B are compatible */
    if ((SM_ROWS_B(A) != SM_ROWS_B(B)) || (SM_COLUMNS_B(A) != SM_COLUMNS_B(B)))
	CAMLreturnT(int, 1);

    /* Call separate routine in B has larger bandwidth(s) than A */
    if ((SM_UBAND_B(B) > SM_UBAND_B(A)) || (SM_LBAND_B(B) > SM_LBAND_B(A)))
	CAMLreturnT(int,
	    matrix_band_scale_add_new(c, MAT_BACKLINK(A),
	      Field(MAT_BACKLINK(B), RECORD_MAT_MATRIXCONTENT_RAWPTR)));

    /* Otherwise, perform operation in-place */
    for (j=0; j<SM_COLUMNS_B(A); j++) {
	A_colj = SM_COLUMN_B(A,j);
	B_colj = SM_COLUMN_B(B,j);
	for (i=-SM_UBAND_B(B); i<=SM_LBAND_B(B); i++)
	    A_colj[i] = c*A_colj[i] + B_colj[i];
    }

    CAMLreturnT(int, 0);
}
#endif

CAMLprim void sunml_matrix_band_scale_addi(value vc, value va)
{
    CAMLparam2(vc, va);
    sundials_ml_index i, j;
    sunrealtype *A_colj;
    sunrealtype c = Double_val(vc);
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

CAMLprim void sunml_matrix_band_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);
    sundials_ml_index i, j, is, ie;
    sunrealtype *col_j, *xd, *yd;
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
CAMLprim value sunml_matrix_band_wrap(DlsMat a)
{
    CAMLparam0();
    CAMLlocal3(vcptr, vcontent, vr);
    CAMLlocal2(vdims, vpayload);

    vcontent = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, a->data, a->N, a->ldim);

    /* a DlsMat is a pointer to a struct _DlsMat */
    vcptr = caml_alloc_final(1, NULL, 1, 20);
    DLSMAT(vcptr) = a;

    vdims = caml_alloc_tuple(RECORD_MAT_BANDDIMENSIONS_SIZE);
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_N,   Val_index(a->N));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_MU,  Val_index(a->mu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_SMU, Val_index(a->s_mu));
    Store_field(vdims, RECORD_MAT_BANDDIMENSIONS_ML,  Val_index(a->ml));

    vpayload = caml_alloc_tuple(RECORD_MAT_BANDDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DATA, vcontent);
    Store_field(vpayload, RECORD_MAT_BANDDATA_DIMS, vdims);

    vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
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
 * this library). The sunml_matrix_sparse_rewrap() function can be called to
 * resynchronize the OCaml values cached in the payload with the underlying
 * C data structures.
 *
 * In Sundials >= 3.0.0, we override the SUNSparseFromDenseMatrix,
 * SUNSparseFromBandMatrix, SUNMatClone_Sparse, SUNMatScaleAddI_Sparse, and
 * SUNMatScaleAdd_Sparse functions so that they update the OCaml payload
 * when the C data structures are reallocated.
 */
#if SUNDIALS_LIB_VERSION < 260
CAMLprim value sunml_matrix_sparse_create(value vm, value vn, value vnnz,
				          value vsformat)
{
    CAMLparam4(vm, vn, vnnz, vsformat);
    CAMLreturn(Val_unit);
}

// Adapted directly from SUNSparseFromDenseMatrix
CAMLprim value sunml_matrix_sparse_from_dense(value vsformat, value vcptrad,
					   value vdroptol)
{
    CAMLparam3(vsformat, vcptrad, vdroptol);
    CAMLreturn(Val_unit);
}

// Adapted directly from SUNSparseFromBandMatrix
CAMLprim value sunml_matrix_sparse_from_band(value vsformat, value vcptrab,
					  value vdroptol)
{
    CAMLparam3(vsformat, vcptrab, vdroptol);
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_matrix_sparse_size(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_matrix_sparse_dims(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLreturn(Val_unit);
}

CAMLprim void sunml_matrix_sparse_scale_add(value vc, value va, value vcptrb)
{
    CAMLparam3(vc, va, vcptrb);
    CAMLreturn0;
}

CAMLprim void sunml_matrix_sparse_scale_addi(value vc, value va)
{
    CAMLparam2(vc, va);
    CAMLreturn0;
}

CAMLprim void sunml_matrix_sparse_matvec(value vcptra, value vx, value vy)
{
    CAMLparam3(vcptra, vx, vy);
    CAMLreturn0;
}

CAMLprim void sunml_matrix_sparse_resize(value va, value vnnz, value vcopy)
{
    CAMLparam3(va, vnnz, vcopy);
    CAMLreturn0;
}

CAMLprim void sunml_matrix_sparse_copy(value vcptra, value vb)
{
    CAMLparam2(vcptra, vb);
    CAMLreturn0;
}

CAMLprim value sunml_matrix_sparse_space(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLreturn(Val_unit);
}

CAMLprim void sunml_matrix_sparse_set_idx(value vcptr, value vj, value vidx)
{
    CAMLparam3(vcptr, vj, vidx);
    CAMLreturn0;
}

CAMLprim value sunml_matrix_sparse_get_idx(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    CAMLreturn(Val_unit);
}

CAMLprim void sunml_matrix_sparse_set_data(value vcptr, value vj, value vv)
{
    CAMLparam3(vcptr, vj, vv);
    CAMLreturn0;
}

CAMLprim value sunml_matrix_sparse_get_val(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    CAMLreturn(Val_unit);
}

CAMLprim void sunml_matrix_sparse_set_val(value vcptr, value vj, value vv)
{
    CAMLparam3(vcptr, vj, vv);
    CAMLreturn0;
}

CAMLprim value sunml_matrix_sparse_get_data(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_matrix_sparse_rewrap(value vm)
{
    CAMLparam1(vm);
    CAMLreturn(Val_unit);
}

CAMLprim void sunml_matrix_sparse_set_to_zero(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLreturn0;
}
#else

static void finalize_mat_content_sparse(value vcptra)
{
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptra);

#if SUNDIALS_LIB_VERSION >= 300
    /* indexvals, indexptrs, and data are freed when the corresponding
       bigarrays are finalized */
    free(content);

#elif SUNDIALS_LIB_VERSION >= 270
    SparseDestroyMat(content);

#else // SUNDIALS_LIB_VERSION < 270
    DestroySparseMat(content);

#endif
}

static void zero_sparse(MAT_CONTENT_SPARSE_TYPE A)
{
    sundials_ml_smat_index i, *indexvals, *indexptrs;
    sunrealtype *data;

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
	indexvals[i] = 0L;
    }
    for (i=0; i < A->NP; i++)
	indexptrs[i] = 0L;
    indexptrs[A->NP] = 0L;
}

static bool matrix_sparse_create_vcptr(sundials_ml_smat_index m,
				       sundials_ml_smat_index n,
				       sundials_ml_smat_index nnz,
				       int sformat,
				       value *pvdata, value *pvidxvals,
				       value *pvidxptrs, value *pvcptr)
{
    CAMLparam0();

#if SUNDIALS_LIB_VERSION >= 300
    sundials_ml_smat_index *indexvals, *indexptrs;
    sunrealtype *data;
    sundials_ml_smat_index np = (sformat == CSC_MAT) ? n : m;
    SUNMatrixContent_Sparse content = NULL;

    data = (sunrealtype *) malloc(nnz * sizeof(sunrealtype));
    if (data == NULL) CAMLreturnT(bool, false);

    indexvals = (sundials_ml_smat_index *)
		    malloc(nnz * sizeof(sundials_ml_smat_index));
    if (indexvals == NULL) {
	free(data);
	CAMLreturnT(bool, false);
    }

    indexptrs = (sundials_ml_smat_index *)
	malloc((np + 1) * sizeof(sundials_ml_smat_index));
    if (indexptrs == NULL) {
	free(data);
	free(indexvals);
	CAMLreturnT(bool, false);
    }

    // allocate the array memory manually (above) to minimize the chances of an
    // out of memory exception (a very slight chance remains...)
    *pvdata = caml_ba_alloc_dims(BIGARRAY_FLOAT | CAML_BA_MANAGED,
				 1, data, nnz);
    *pvidxvals = caml_ba_alloc_dims(BIGARRAY_INDEX | CAML_BA_MANAGED,
				 1, indexvals, nnz);
    *pvidxptrs = caml_ba_alloc_dims(BIGARRAY_INDEX | CAML_BA_MANAGED,
				 1, indexptrs, np + 1);

    // Setup the C-side content
    content = (SUNMatrixContent_Sparse) malloc(sizeof *content);
    if (content == NULL) CAMLreturnT(bool, false);

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
    if (a == NULL) CAMLreturnT(int, false);

#if SUNDIALS_LIB_VERSION >= 270
    *pvidxvals = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->indexvals, a->NNZ);
    *pvidxptrs = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->indexptrs, a->NP + 1);
#else
    *pvidxvals = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->rowvals, a->NNZ);
    *pvidxptrs = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->colptrs, a->N + 1);
#endif
    *pvdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, a->data, a->NNZ);

    *pvcptr = caml_alloc_final(1, finalize_mat_content_sparse, 1, 20);
    SLSMAT(*pvcptr) = a;

#endif

    CAMLreturnT(bool, true);
}

static bool matrix_sparse_create_mat(sundials_ml_smat_index m,
				     sundials_ml_smat_index n,
				     sundials_ml_smat_index nnz,
				     int sformat, value *vr)
{
    CAMLparam0();
    CAMLlocal5(vdata, vidxvals, vidxptrs, vcptr, vpayload);

    if (! matrix_sparse_create_vcptr(m, n, nnz, sformat,
	       &vdata, &vidxvals, &vidxptrs, &vcptr) )
	CAMLreturnT(bool, false);

    vpayload = caml_alloc_tuple(RECORD_MAT_SPARSEDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXVALS, vidxvals);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXPTRS, vidxptrs);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_DATA,    vdata);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_SFORMAT,
	    MAT_TO_SFORMAT(sformat));

    *vr = caml_alloc_tuple(RECORD_MAT_MATRIXCONTENT_SIZE);
    Store_field(*vr, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vpayload);
    Store_field(*vr, RECORD_MAT_MATRIXCONTENT_RAWPTR,  vcptr);
    Store_field(*vr, RECORD_MAT_MATRIXCONTENT_VALID,   Val_bool(1));

    CAMLreturnT(bool, true);
}

CAMLprim value sunml_matrix_sparse_create(value vm, value vn, value vnnz,
				       value vsformat)
{
    CAMLparam4(vm, vn, vnnz, vsformat);
    CAMLlocal1(vr);

    if (! matrix_sparse_create_mat(SmatIndex_val(vm), SmatIndex_val(vn), SmatIndex_val(vnnz),
				   MAT_FROM_SFORMAT(vsformat), &vr) )
	caml_raise_out_of_memory();

    CAMLreturn(vr);
}

// Adapted directly from SUNSparseFromDenseMatrix
CAMLprim value sunml_matrix_sparse_from_dense(value vsformat, value vcptrad,
					   value vdroptol)
{
    CAMLparam3(vsformat, vcptrad, vdroptol);
    CAMLlocal2(vr, vcptr);
    int sformat = MAT_FROM_SFORMAT(vsformat);
    sunrealtype droptol = Double_val(vdroptol);
    sunrealtype *As_data;

    sundials_ml_smat_index i, j, nnz, *As_indexvals, *As_indexptrs;
    sundials_ml_smat_index M, N;
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
    if (! matrix_sparse_create_mat(M, N, nnz, sformat, &vr) )
	caml_raise_out_of_memory();
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
CAMLprim value sunml_matrix_sparse_from_band(value vsformat, value vcptrab,
					  value vdroptol)
{
    CAMLparam3(vsformat, vcptrab, vdroptol);
    CAMLlocal2(vr, vcptr);
    int sformat = MAT_FROM_SFORMAT(vsformat);
    sunrealtype droptol = Double_val(vdroptol);
    sunrealtype *As_data;

    sundials_ml_smat_index i, j, nnz, *As_indexptrs, *As_indexvals;
    sundials_ml_smat_index M, N;
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
    if (! matrix_sparse_create_mat(M, N, nnz, sformat, &vr) )
	caml_raise_out_of_memory();
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

CAMLprim value sunml_matrix_sparse_size(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_SPARSE_TYPE a = MAT_CONTENT_SPARSE(vcptr);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_index(a->M));
    Store_field(vr, 1, Val_index(a->N));

    CAMLreturn(vr);
}

CAMLprim value sunml_matrix_sparse_dims(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_index(content->NNZ));
    Store_field(vr, 1, Val_index(content->NP));

    CAMLreturn(vr);
}

#if SUNDIALS_LIB_VERSION < 300
static value sparse_wrap_payload(SlsMat a)
{
    CAMLparam0();
    CAMLlocal4(vpayload, vidxvals, vidxptrs, vdata);
    int format;

#if SUNDIALS_LIB_VERSION >= 270
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->indexvals, a->NNZ);
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->indexptrs, a->NP + 1);
    format = a->sparsetype;
#else
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->rowvals, a->NNZ);
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INDEX, 1, a->colptrs, a->N + 1);
    format = CSC_MAT;
#endif
    vdata = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, a->data, a->NNZ);

    vpayload = caml_alloc_tuple(RECORD_MAT_SPARSEDATA_SIZE);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXVALS, vidxvals);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_IDXPTRS, vidxptrs);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_DATA,    vdata);
    Store_field(vpayload, RECORD_MAT_SPARSEDATA_SFORMAT, MAT_TO_SFORMAT(format));

    CAMLreturn(vpayload);
}
#endif

static bool matrix_sparse_resize(value va, sundials_ml_smat_index nnz,
	bool docopy, bool dofree)
{
    CAMLparam1(va);
    CAMLlocal5(vcptr, vdata, vidxvals, vpayload, vnewpayload);
    CAMLlocal1(vidxptrs);
    sundials_ml_smat_index old_nnz, i, np;
    sundials_ml_smat_index *old_indexvals, *new_indexvals;
    sundials_ml_smat_index *old_indexptrs, *new_indexptrs;
    sunrealtype *old_data, *new_data;
    MAT_CONTENT_SPARSE_TYPE A;

    vcptr    = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    // holding vpayload ensures that the underlying arrays are not gc-ed.
    vpayload = Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD);
    A = MAT_CONTENT_SPARSE(vcptr);

#if SUNDIALS_LIB_VERSION >= 270
    old_indexptrs = A->indexptrs;
    old_indexvals = A->indexvals;
    np = A->NP;
#else
    old_indexptrs = A->rowvals;
    old_indexvals = A->colptrs;
    np = A->N;
#endif
    old_data = A->data;
    old_nnz = A->NNZ;

    if (nnz <= 0)
	nnz = old_indexptrs[np];

    new_data = (sunrealtype *) malloc(nnz * sizeof(sunrealtype));
    if (new_data == NULL) CAMLreturnT(bool, false);

    new_indexvals = (sundials_ml_smat_index *)
			malloc(nnz * sizeof(sundials_ml_smat_index));
    if (new_indexvals == NULL) {
	free(new_data);
	CAMLreturnT(bool, false);
    }

    new_indexptrs = (sundials_ml_smat_index *)
			malloc((np + 1) * sizeof(sundials_ml_smat_index));
    if (new_indexptrs == NULL) {
	free(new_data);
	free(new_indexvals);
	CAMLreturnT(bool, false);
    }

#if SUNDIALS_LIB_VERSION >= 300
    // allocate the array memory manually (above) to minimize the chances of an
    // out of memory exception (a very slight chance remains...)
    vdata = caml_ba_alloc_dims(BIGARRAY_FLOAT | CAML_BA_MANAGED, 1,
		new_data, nnz);
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INDEX | CAML_BA_MANAGED, 1,
		new_indexvals, nnz);
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INDEX | CAML_BA_MANAGED, 1,
		new_indexptrs, np + 1);

    /* Strictly speaking, it is not necessary to reallocate and recopy the
       indexptrs. We do it for two reasons:
       (1) it is conceptually simpler to deep-copy the sparse payload
       (2) it simplifies functions like sunml_matrix_sparse_scale_add,
	   since they can "pretend" that resize duplicates the
	   matrix argument if they cache the old underlying arrays
	   and call with free=0. */

    // very slight chance of raising an out of memory exception and messing
    // things up badly
    // no caml_alloc_shr_no_track_noexc in OCaml 4.03.0
    // https://caml.inria.fr/mantis/view.php?id=7572
    vnewpayload = caml_alloc_tuple(RECORD_MAT_SPARSEDATA_SIZE);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_IDXVALS, vidxvals);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_IDXPTRS, vidxptrs);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_DATA,    vdata);
    Store_field(vnewpayload, RECORD_MAT_SPARSEDATA_SFORMAT,
		MAT_TO_SFORMAT(A->sparsetype));

    Store_field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD, vnewpayload);

    A->data = new_data;
    A->indexptrs = new_indexptrs;
    A->indexvals = new_indexvals;
#else

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

    if (docopy) {
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
    if (dofree) {
	free(old_data);
	free(old_indexptrs);
	free(old_indexvals);
	// For >= 300, the OCaml bigarrays ensure gc of these arrays
    }
#endif

    CAMLreturnT(bool, true);
}

// Adapted directly from SUNMatScaleAdd_Sparse
static bool matrix_sparse_scale_add(sunrealtype c, value va, value vcptrb)
{
    CAMLparam2(va, vcptrb);
    CAMLlocal5(vcptra, vdatac, vidxvalsc, vidxptrsc, vcptrc);
    CAMLlocal1(vpayload);
    sundials_ml_smat_index j, i, p, newvals, cend, nz;
    bool newmat;
    sundials_ml_smat_index *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
    sunrealtype *x, *Ax, *Bx, *Cx;
    sundials_ml_smat_index M, N;
    sundials_ml_smat_index *A_indexptrs, *A_indexvals, *B_indexptrs, *B_indexvals;
    MAT_CONTENT_SPARSE_TYPE A, B;

    vcptra = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    A = MAT_CONTENT_SPARSE(vcptra);
    B = MAT_CONTENT_SPARSE(vcptrb);

    /* Perform operation */

#if SUNDIALS_LIB_VERSION >= 270
    /* if A is CSR matrix, transpose M and N */
    if (A->sparsetype == CSC_MAT) {
	M = A->M;
	N = A->N;
    } else {
	M = A->N;
	N = A->M;
    }
#else
    M = A->M;
    N = A->N;
#endif

    /* create work arrays for row indices and nonzero column values */
    w = (sundials_ml_smat_index *) malloc(M * sizeof(sundials_ml_smat_index));
    if (w == NULL) CAMLreturnT(int, 0);

    x = (sunrealtype *) malloc(M * sizeof(sunrealtype));
    if (x == NULL) {
	free(w);
	CAMLreturnT(bool, false);
    }

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

    /* determine if A already contains the sparsity pattern of B */
    newvals = 0;
    for (j=0; j < N; j++) {

	/* clear work array */
	for (i=0; i < M; i++)
	    w[i] = 0;

	/* scan column of A, incrementing w by one */
	for (i=A_indexptrs[j]; i < A_indexptrs[j+1]; i++)
	    w[A_indexvals[i]] += 1;

	/* scan column of B, decrementing w by one */
	for (i=B_indexptrs[j]; i < B_indexptrs[j+1]; i++)
	    w[B_indexvals[i]] -= 1;

	/* if any entry of w is negative, A doesn't contain B's sparsity,
	   so increment necessary storage counter */
	for (i = 0; i < M; i++) {
	  if (w[i] < 0)  newvals += 1;
	}
    }

    /* If extra nonzeros required, check whether A has sufficient storage space
       for new nonzero entries (so B can be inserted into existing storage) */
    newmat = (newvals > (A->NNZ - A_indexptrs[N]));

    /* perform operation based on existing/necessary structure */

    /* case 1: A already contains sparsity pattern of B */
    if (newvals == 0) {
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

    /*   case 2: A has sufficient storage, but does not already contain B's sparsity */
    } else if (!newmat) {
        /* determine storage location where last column (row) should end */
        nz = A_indexptrs[N] + newvals;

        /* store pointer past last column (row) from original A,
           and store updated value in revised A */
        cend = A_indexptrs[N];
        A_indexptrs[N] = nz;

        /* iterate through columns (rows) backwards */
        for (j = N-1; j >= 0; j--) {

	    /* clear out temporary arrays for this column (row) */
	    for (i = 0; i < M; i++) {
		w[i] = 0;
		x[i] = 0.0;
	    }

	    /* iterate down column (row) of A, collecting nonzeros */
	    for (p = A_indexptrs[j]; p < cend; p++) {
		w[A_indexvals[p]] += 1; /* indicate that row (column) is filled */
		x[A_indexvals[p]] = c*A->data[p];    /* collect/scale value */
	    }

	    /* iterate down column of B, collecting nonzeros */
	    for (p = B_indexptrs[j]; p < B_indexptrs[j+1]; p++) {
		w[B_indexvals[p]] += 1;          /* indicate that row is filled */
		x[B_indexvals[p]] += B->data[p]; /* collect value */
	    }

	    /* fill entries of A with this column's (row's) data */
	    for (i = M-1; i >= 0; i--) {
	        if (w[i] > 0) {
		    A_indexvals[--nz] = i;
		    A->data[nz] = x[i];
	        }
	    }

	    /* store ptr past this col (row) from orig A,
	       update value for new A */
	    cend = A_indexptrs[j];
	    A_indexptrs[j] = nz;
        }

    /*   case 3: A must be reallocated with sufficient storage */
    } else {
	/* access data from (pre resize) CSR structures (return if failure) */
	Ap = A_indexptrs;
	Ai = A_indexvals;
	Ax = A->data;
	Bp = B_indexptrs;
	Bi = B_indexvals;
	Bx = B->data;

	/* reallocate memory within A - no-copy, no-free */
	if (! matrix_sparse_resize(va, Ap[N] + newvals, 0, 0) ) {
	    free(w);
	    free(x);
	    CAMLreturnT(bool, false);
	}
	zero_sparse(A);

	/* access data from (post resize) CSR structures (return if failure) */
#if SUNDIALS_LIB_VERSION >= 270
	Cp = A->indexptrs;
	Ci = A->indexvals;
#else
	Cp = A->colptrs;
	Ci = A->rowvals;
#endif
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

    CAMLreturnT(bool, true);
}

#if false
CAMLprim void ml_debug_sparse(value va)
{
    CAMLparam1(va);
    CAMLlocal5(vpayload, vidxvals, vidxptrs, vdata, vcptr);
    MAT_CONTENT_SPARSE_TYPE A;
    struct caml_ba_array *ba;

    vpayload = Field(va, RECORD_MAT_MATRIXCONTENT_PAYLOAD);
    vcptr = Field(va, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    A = MAT_CONTENT_SPARSE(vcptr);

    vidxvals = Field(vpayload, RECORD_MAT_SPARSEDATA_IDXVALS);
    vidxptrs = Field(vpayload, RECORD_MAT_SPARSEDATA_IDXPTRS);
    vdata = Field(vpayload, RECORD_MAT_SPARSEDATA_DATA);

    printf("**ml_debug_sparse\n");
    printf("--vmat=%p vpayload=%p mat=%p mat->NP=%lld\n",
	    (void *)va, (void *)vpayload, (void *)A, A->NP);
    printf("--mat->indexptrs[mat->NP]=%lld\n", A->indexptrs[A->NP]);

    ba = Caml_ba_array_val(vidxvals);
    printf("--vidxvals=%p data=%p (len=%6ld) (A->indexvals=%p)\n",
	    (void *)vidxvals, (void *)ba->data, ba->dim[0],
	    (void *)A->indexvals);
    ba = Caml_ba_array_val(vidxptrs);
    printf("--vidxptrs=%p data=%p (len=%6ld) (A->indexptrs=%p)\n",
	    (void *)vidxptrs, (void *)ba->data, ba->dim[0],
	    (void *)A->indexptrs);
    ba = Caml_ba_array_val(vdata);
    printf("--   vdata=%p data=%p (len=%6ld)      (A->data=%p)\n",
	    (void *)vdata, (void *)ba->data, ba->dim[0], (void *)A->data);
    fflush(stdout);

    CAMLreturn0;
}
#endif

CAMLprim void sunml_matrix_sparse_scale_add(value vc, value va, value vcptrb)
{
    CAMLparam3(vc, va, vcptrb);
    if (! matrix_sparse_scale_add(Double_val(vc), va, vcptrb) )
	caml_raise_out_of_memory();
    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
static int csmat_sparse_scale_add(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
    return (matrix_sparse_scale_add(c, MAT_BACKLINK(A),
	      Field(MAT_BACKLINK(B), RECORD_MAT_MATRIXCONTENT_RAWPTR))
	    ? 0 : 1);
}
#endif

// Adapted directly from SUNMatScaleAddI_Sparse
static bool matrix_sparse_scale_addi(sunrealtype c, value va)
{
    CAMLparam1(va);
    CAMLlocal2(vcptr, vpayload);
    CAMLlocal4(vdatac, vidxvalsc, vidxptrsc, vcptrc);
    sundials_ml_smat_index j, i, p, nz, cend, newvals;
    bool newmat, found;
    sundials_ml_smat_index *w, *Ap, *Ai, *Cp, *Ci;
    sunrealtype *x, *Ax, *Cx;
    sundials_ml_smat_index M, N, *A_indexptrs, *A_indexvals;
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

#if SUNDIALS_LIB_VERSION >= 270
    if (A->sparsetype == CSC_MAT) {
	M = A->M;
	N = A->N;
    } else {
	M = A->N;
	N = A->M;
    }
#else
    M = A->M;
    N = A->N;
#endif

    /* determine if A already contains values on the diagonal (hence
       no memory allocation necessary), and calculate the number of non-zeroes
       required if we create a new matrix (instead of reallocating). */
    newvals = 0;
    for (j=0; j < SUNMIN(M, N); j++) {
	/* scan column (row if CSR) of A, searching for diagonal value */
	found = 0;
	for (i=A_indexptrs[j]; i < A_indexptrs[j+1]; i++) {
	    if (A_indexvals[i] == j) {
		found = 1;
		break;
	    }
	}
	/* if no diagonal found, signal new matrix */
	if (!found) newvals += 1;
    }

    /* If extra nonzeros required, check whether matrix has sufficient
       storage space for new nonzero entries  (so I can be inserted into
       existing storage) */
    newmat = (newvals > (A->NNZ - A_indexptrs[N]));

    /* perform operation based on existing/necessary structure */

    /*   case 1: A already contains a diagonal */
    if (newvals == 0) {
	/* iterate through columns, adding 1.0 to diagonal */
	for (j=0; j < SUNMIN(M, N); j++)
	    for (i=A_indexptrs[j]; i < A_indexptrs[j+1]; i++)
		if (A_indexvals[i] == j) {
		    A->data[i] = 1.0 + c * A->data[i];
		} else {
		    A->data[i] = c * A->data[i];
		}

    /*   case 2: A has sufficient storage,
                 but does not already contain a diagonal */
    } else if (!newmat) {

      /* create work arrays for nonzero indices and values in a
         single column (row) */
      w = (sundials_ml_smat_index *) malloc(M * sizeof(sundials_ml_smat_index));
      if (w == NULL) CAMLreturnT(bool, false);
      x = (sunrealtype *) malloc(M * sizeof(sunrealtype));
      if (x == NULL) {
	  free(w);
	  CAMLreturnT(bool, false);
      }

      /* determine storage location where last column (row) should end */
      nz = A_indexptrs[N] + newvals;

      /* store pointer past last column (row) from original A,
         and store updated value in revised A */
      cend = A_indexptrs[N];
      A_indexptrs[N] = nz;

      /* iterate through columns (rows) backwards */
      for (j = N-1; j >= 0; j--) {

        /* clear out temporary arrays for this column (row) */
        for (i = 0; i < M; i++) {
          w[i] = 0;
          x[i] = 0.0;
        }

        /* iterate down column (row) of A, collecting nonzeros */
        for (p = A_indexptrs[j]; p < cend; p++) {
          w[A_indexvals[p]] += 1;    /* indicate that row (column) is filled */
          x[A_indexvals[p]] = c * A->data[p];  /* collect/scale value */
        }

        /* add identity to this column (row) */
        if (j < M) {
          w[j] += 1;     /* indicate that row (column) is filled */
          x[j] += 1.0;   /* update value */
        }

        /* fill entries of A with this column's (row's) data */
        for (i = M-1; i >= 0; i--) {
          if (w[i] > 0) {
            A_indexvals[--nz] = i;
            A->data[nz] = x[i];
          }
        }

        /* store ptr past this col (row) from orig A, update value for new A */
        cend = A_indexptrs[j];
        A_indexptrs[j] = nz;
      }

      /* clean up */
      free(w);
      free(x);

    /*   case 3: A must be reallocated with sufficient storage */
    } else {
	/* create work arrays for row indices and nonzero column values */
	w = (sundials_ml_smat_index *) malloc(M * sizeof(sundials_ml_smat_index));
	if (w == NULL) CAMLreturnT(bool, false);

	x = (sunrealtype *) malloc(M * sizeof(sunrealtype));
	if (x == NULL) {
	    free(w);
	    CAMLreturnT(bool, false);
	}

	/* access data from (pre resize) CSR structures (return if failure) */
	Ap = A_indexptrs;
	Ai = A_indexvals;
	Ax = A->data;

	/* reallocate memory within A -- no-copy, no-free */
	if (! matrix_sparse_resize(va, Ap[N] + newvals, 0, 0)) {
	    free(w);
	    free(x);
	    CAMLreturnT(bool, false);
	}
	zero_sparse(A);

	/* access data from (post resize) CSR structures (return if failure) */
#if SUNDIALS_LIB_VERSION >= 270
	Cp = A->indexptrs;
	Ci = A->indexvals;
#else
	Cp = A->colptrs;
	Ci = A->rowvals;
#endif
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
		if (w[i] > 0) {
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

    CAMLreturnT(bool, true);
}

CAMLprim void sunml_matrix_sparse_scale_addi(value vc, value va)
{
    CAMLparam2(vc, va);
    if (! matrix_sparse_scale_addi(Double_val(vc), va))
	caml_raise_out_of_memory();
    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
static int csmat_sparse_scale_addi(sunrealtype c, SUNMatrix A)
{
    return (matrix_sparse_scale_addi(c, MAT_BACKLINK(A)) ? 0 : 1);
}
#endif

// Adapted directly from SUNMatMatvec_Sparse, Matvec_SparseCSC, and
// Matvec_SparseCSR
CAMLprim void sunml_matrix_sparse_matvec(value vcptra, value vx, value vy)
{
    CAMLparam3(vcptra, vx, vy);
    MAT_CONTENT_SPARSE_TYPE A = MAT_CONTENT_SPARSE(vcptra);
    sundials_ml_smat_index i, j;
    sundials_ml_smat_index *Ap, *Ai;
    sunrealtype *Ax, *xd, *yd;

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

CAMLprim void sunml_matrix_sparse_resize(value va, value vnnz, value vcopy)
{
    CAMLparam3(va, vnnz, vcopy);
    if (! matrix_sparse_resize(va, SmatIndex_val(vnnz), Bool_val(vcopy), 1))
	caml_raise_out_of_memory();
    CAMLreturn0;
}

// Adapted directly from SUNMatCopy_Sparse
static bool matrix_sparse_copy(value vcptra, value vb)
{
    CAMLparam2(vcptra, vb);
    CAMLlocal1(vcptrb);
    MAT_CONTENT_SPARSE_TYPE A = MAT_CONTENT_SPARSE(vcptra);
    MAT_CONTENT_SPARSE_TYPE B;
    sundials_ml_smat_index i, A_nz, *A_indexptrs, *A_indexvals;
    sundials_ml_smat_index *B_indexvals, *B_indexptrs;

    vcptrb = Field(vb, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    B = MAT_CONTENT_SPARSE(vcptrb);

#if SUNDIALS_LIB_VERSION >= 270
    A_indexptrs = A->indexptrs;
    A_indexvals = A->indexvals;
#else
    A_indexptrs = A->rowvals;
    A_indexvals = A->colptrs;
#endif

    /* Perform operation */
    A_nz = A_indexptrs[A->NP];

    /* ensure that B is allocated with at least as
    much memory as we have nonzeros in A */
    if (B->NNZ < A_nz)
	if (! matrix_sparse_resize(vb, A_nz, 0, 1)) // no-copy, free
	    CAMLreturnT(bool, false);

#if SUNDIALS_LIB_VERSION >= 270
    B_indexptrs = B->indexptrs;
    B_indexvals = B->indexvals;
#else
    B_indexptrs = B->rowvals;
    B_indexvals = B->colptrs;
#endif

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

    CAMLreturnT(bool, true);
}

CAMLprim void sunml_matrix_sparse_copy(value vcptra, value vb)
{
    CAMLparam2(vcptra, vb);
    if (! matrix_sparse_copy(vcptra, vb))
	caml_raise_out_of_memory();
    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 300
// Adapted directly from SUNMatCopy_Sparse
static int csmat_sparse_copy(SUNMatrix A, SUNMatrix B)
{
    return (matrix_sparse_copy(
		Field(MAT_BACKLINK(A), RECORD_MAT_MATRIXCONTENT_RAWPTR),
		MAT_BACKLINK(B)) ? 0 : 1);
}
#endif

CAMLprim value sunml_matrix_sparse_space(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(vr);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);

    // Directly adapted from SUNMatSpace_Dense
    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_index(content->NNZ));
    Store_field(vr, 1, Val_index(10 + content->NP + content->NNZ));

    CAMLreturn(vr);
}

#if 520 <= SUNDIALS_LIB_VERSION
// Adapted directly from sunmatrix_sparse.c:
// - since not exposed by Sundials
// - apply directly to sparse content (not SUNMatrix)
static int csmat_sparse_format_convert(const MAT_CONTENT_SPARSE_TYPE A,
				       MAT_CONTENT_SPARSE_TYPE B)
{
    sunrealtype *Ax, *Bx;
    sundials_ml_index *Ap, *Aj;
    sundials_ml_index *Bp, *Bi;
    sundials_ml_index n_row, n_col, nnz;
    sundials_ml_index n, col, csum, row, last;

    Ap = A->indexptrs;
    Aj = A->indexvals;
    Ax = A->data;
    
    n_row = (A->sparsetype == CSR_MAT) ? A->M : A->N;
    n_col = (A->sparsetype == CSR_MAT) ? A->N : A->M;

    Bp = B->indexptrs;
    Bi = B->indexvals;
    Bx = B->data;

    nnz = Ap[n_row];

    zero_sparse(B);

    /* compute number of non-zero entries per column (if CSR) or per row (if CSC) of A */
    for (n = 0; n < nnz; n++)
    {
        Bp[Aj[n]]++;
    }

    /* cumualtive sum the nnz per column to get Bp[] */
    for (col = 0, csum = 0; col < n_col; col++)
    {
        sundials_ml_index temp  = Bp[col];
        Bp[col] = csum;
        csum += temp;
    }
    Bp[n_col] = nnz;

    for (row = 0; row < n_row; row++)
    {
        sundials_ml_index jj;
        for (jj = Ap[row]; jj < Ap[row+1]; jj++)
        {
            sundials_ml_index col  = Aj[jj];
            sundials_ml_index dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }

    for (col = 0, last = 0; col <= n_col; col++)
    {
        sundials_ml_index temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }

    return 0;
}
#endif

CAMLprim value sunml_matrix_sparse_tocsr(value vcptra)
{
    CAMLparam1(vcptra);
    CAMLlocal2(vcptrb, vb);
#if 520 <= SUNDIALS_LIB_VERSION
    MAT_CONTENT_SPARSE_TYPE a = MAT_CONTENT_SPARSE(vcptra);
    MAT_CONTENT_SPARSE_TYPE b;

    if (! matrix_sparse_create_mat(a->M, a->N, a->NNZ, CSR_MAT, &vb) )
	caml_raise_out_of_memory();

    vcptrb = Field(vb, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    b = MAT_CONTENT_SPARSE(vcptrb);
    csmat_sparse_format_convert(a, b);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vb);
}

CAMLprim value sunml_matrix_sparse_tocsc(value vcptra)
{
    CAMLparam1(vcptra);
    CAMLlocal2(vcptrb, vb);
#if 520 <= SUNDIALS_LIB_VERSION
    MAT_CONTENT_SPARSE_TYPE a = MAT_CONTENT_SPARSE(vcptra);
    MAT_CONTENT_SPARSE_TYPE b;

    if (! matrix_sparse_create_mat(a->M, a->N, a->NNZ, CSC_MAT, &vb) )
	caml_raise_out_of_memory();

    vcptrb = Field(vb, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    b = MAT_CONTENT_SPARSE(vcptrb);
    csmat_sparse_format_convert(a, b);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vb);
}

// Sundials < 3.0.0
CAMLprim void sunml_matrix_sparse_set_idx(value vcptr, value vj, value vidx)
{
    CAMLparam3(vcptr, vj, vidx);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_smat_index j = SmatIndex_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NP + 1) caml_invalid_argument("j");
#endif

#if SUNDIALS_LIB_VERSION >= 270
    content->indexptrs[j] = SmatIndex_val(vidx);
#else
    content->rowvals[j] = SmatIndex_val(vidx);
#endif

    CAMLreturn0;
}

// Sundials < 3.0.0
CAMLprim value sunml_matrix_sparse_get_idx(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_smat_index j = SmatIndex_val(vj);
    sundials_ml_index r;

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NP + 1) caml_invalid_argument("j");
#endif

#if SUNDIALS_LIB_VERSION >= 270
    r = content->indexptrs[j];
#else
    r = content->rowvals[j];
#endif

    CAMLreturn(Val_index(r));
}

// Sundials < 3.0.0
CAMLprim void sunml_matrix_sparse_set_data(value vcptr, value vj, value vv)
{
    CAMLparam3(vcptr, vj, vv);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_smat_index j = SmatIndex_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    content->data[j] = Double_val(vv);

    CAMLreturn0;
}

// Sundials < 3.0.0
CAMLprim value sunml_matrix_sparse_get_val(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_smat_index j = SmatIndex_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    CAMLreturn(Val_SmatIndex(content->indexvals[j]));
}

// Sundials < 3.0.0
CAMLprim void sunml_matrix_sparse_set_val(value vcptr, value vj, value vv)
{
    CAMLparam3(vcptr, vj, vv);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_smat_index j = SmatIndex_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    content->indexvals[j] = SmatIndex_val(vv);

    CAMLreturn0;
}

// Sundials < 3.0.0
CAMLprim value sunml_matrix_sparse_get_data(value vcptr, value vj)
{
    CAMLparam2(vcptr, vj);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);
    sundials_ml_smat_index j = SmatIndex_val(vj);

#if SUNDIALS_ML_SAFE == 1
    if (j < 0 || j >= content->NNZ) caml_invalid_argument("j");
#endif

    CAMLreturn(caml_copy_double(content->data[j]));
}

// Sundials < 3.0.0
// Reconnects the OCaml payload to the underlying C data if necessary
CAMLprim value sunml_matrix_sparse_rewrap(value vm)
{
    CAMLparam1(vm);
    CAMLlocal2(vcptr, vpayload);

#if SUNDIALS_LIB_VERSION >= 300
    caml_failwith("sunml_matrix_sparse_rewrap should not be called!");
#else
    MAT_CONTENT_SPARSE_TYPE content;
    void *ba_data;

    vcptr = Field(vm, RECORD_MAT_MATRIXCONTENT_RAWPTR);
    content = MAT_CONTENT_SPARSE(vcptr);

    vpayload = Field(vm, RECORD_MAT_MATRIXCONTENT_PAYLOAD);
    ba_data = Caml_ba_data_val(Field(vpayload, RECORD_MAT_SPARSEDATA_DATA));

    // Has the underlying data array changed?
    if (content->data != ba_data)
	vpayload = sparse_wrap_payload(content);
#endif

    CAMLreturn(vpayload);
}

#if SUNDIALS_LIB_VERSION < 300
CAMLprim value sunml_matrix_sparse_wrap(SlsMat a)
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
CAMLprim void sunml_matrix_sparse_set_to_zero(value vcptr)
{
    CAMLparam1(vcptr);
    MAT_CONTENT_SPARSE_TYPE content = MAT_CONTENT_SPARSE(vcptr);

    /* Perform operation */
    zero_sparse(content);

    CAMLreturn0;
}

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Matrix
 */

#if 300 <= SUNDIALS_LIB_VERSION
/*
 *  Type	    content		    backlink
 *  --------+-------------------------------------------------
 *  dense   | MAT_CONTENT_DENSE_TYPE	    Dense.t  (matrix_content)
 *  band    | MAT_CONTENT_BAND_TYPE	    Band.t   (matrix_content)
 *  sparse  | MAT_CONTENT_SPARSE_TYPE	    Sparse.t (matrix_content)
 *  custom  | matrix_ops record		    'm (custom data)
 */
static SUNMatrix alloc_smat(void *content, value backlink, value context,
			    bool content_is_value)
{
    SUNMatrix smat;

    /* Alloc memory in C heap */
    smat = (SUNMatrix)calloc(1, sizeof(struct csmat));
    if (smat == NULL) return NULL;

    smat->ops = (SUNMatrix_Ops)calloc(1, sizeof(struct _generic_SUNMatrix_Ops));
    if (smat->ops == NULL) { free(smat); return(NULL); }

    smat->content = content;
    if (content_is_value) // NB: free_smat does not unregister the content
	caml_register_generational_global_root((value *)&(smat->content));

    MAT_BACKLINK(smat) = backlink;
    caml_register_generational_global_root(&MAT_BACKLINK(smat));

    MAT_CONTEXT(smat) = context;
    caml_register_generational_global_root(&MAT_CONTEXT(smat));

    return smat;
}

/*
 * For custom matrices (Custom.t, ArrayDense.t, ArrayBand.t) only, the
 * content field stores an OCaml pair
 *	((ops : ('m, 'd) matrix_ops), (id : ('k, 'm, 'nd, 'nk) id))
 *
 * MAT_CUSTOM_CONTENT accesses the pair.
 * MAT_OP_TABLE accesses the first element of the pair.
 * MAT_OCAML_ID accesses the second element of the pair.
 *
 * GET_OP returns specific fields within the first element of the pair.
 */
#define MAT_CUSTOM_CONTENT(smat)  ((smat)->content)
#define MAT_OP_TABLE(smat) (Field((value)MAT_CUSTOM_CONTENT(smat), 0))
#define MAT_OCAML_ID(smat) (Field((value)MAT_CUSTOM_CONTENT(smat), 1))
#define GET_OP(smat, x) (Field((value)MAT_OP_TABLE(smat), x))

static void free_smat(SUNMatrix smat)
{
    caml_remove_generational_global_root(&MAT_BACKLINK(smat));
    caml_remove_generational_global_root(&MAT_CONTEXT(smat));
    // smat->content is not owned by the sunmatrix
    free(smat->ops);
    free(smat);
}

static void free_custom_smat(SUNMatrix smat)
{
    caml_remove_generational_global_root((value *)&MAT_CUSTOM_CONTENT(smat));
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
#if 500 <= SUNDIALS_LIB_VERSION
    ops->matvecsetup = src->ops->matvecsetup;
#endif
    ops->matvec    = src->ops->matvec;
    ops->space     = src->ops->space;
}

static SUNMatrix csmat_dense_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(vcontenta, vcontentb);
    SUNMatrix B;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = sunml_matrix_dense_create(Val_index(SM_ROWS_D(A)),
				       Val_index(SM_COLUMNS_D(A)));

    B = alloc_smat(
	    MAT_CONTENT(Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR)),
	    vcontentb, MAT_CONTEXT(A), false);
    csmat_clone_ops(B, A);
#if 600 <= SUNDIALS_LIB_VERSION
    B->sunctx = A->sunctx;
#endif

    CAMLreturnT(SUNMatrix, B);
}

static SUNMatrix csmat_band_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal4(vcontenta, vpayloada, vcontentb, vpayloadb);
    SUNMatrix B;

    vcontenta = MAT_BACKLINK(A);
    vcontentb = sunml_matrix_band_create_mat(SM_COLUMNS_B(A), SM_UBAND_B(A),
				          SM_LBAND_B(A), SM_SUBAND_B(A));

    B = alloc_smat(
	    MAT_CONTENT(Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR)),
	    vcontentb, MAT_CONTEXT(A), false);
    csmat_clone_ops(B, A);
#if 600 <= SUNDIALS_LIB_VERSION
    B->sunctx = A->sunctx;
#endif

    CAMLreturnT(SUNMatrix, B);
}

static SUNMatrix csmat_sparse_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal4(vcontenta, vpayloada, vcontentb, vpayloadb);
    SUNMatrix B;

    vcontenta = MAT_BACKLINK(A);
    if (! matrix_sparse_create_mat(SM_ROWS_S(A), SM_COLUMNS_S(A),
				   SM_NNZ_S(A), SM_SPARSETYPE_S(A),
				   &vcontentb) )
	CAMLreturnT(SUNMatrix, NULL);

    B = alloc_smat(
	    MAT_CONTENT(Field(vcontentb, RECORD_MAT_MATRIXCONTENT_RAWPTR)),
	    vcontentb, MAT_CONTEXT(A), false);
    csmat_clone_ops(B, A);
#if 600 <= SUNDIALS_LIB_VERSION
    B->sunctx = A->sunctx;
#endif

    CAMLreturnT(SUNMatrix, B);
}

static SUNMatrix csmat_custom_clone(SUNMatrix A);
static SUNMatrix_ID csmat_custom_getid(SUNMatrix A);
static int csmat_custom_zero(SUNMatrix A);
static int csmat_custom_copy(SUNMatrix A, SUNMatrix B);
static int csmat_custom_scale_add(sunrealtype c, SUNMatrix A, SUNMatrix B);
static int csmat_custom_scale_addi(sunrealtype c, SUNMatrix A);
#if 500 <= SUNDIALS_LIB_VERSION
static int csmat_custom_matvecsetup(SUNMatrix A);
#endif
static int csmat_custom_matvec(SUNMatrix A, N_Vector x, N_Vector y);
static int csmat_custom_space(SUNMatrix A, long int *lenrw, long int *leniw);

#endif

CAMLprim value sunml_matrix_wrap(value vid, value vcontent, value vpayload,
				 value vhasmatvecsetup, value vcontext)
{
    CAMLparam5(vid, vcontent, vpayload, vhasmatvecsetup, vcontext);
    CAMLlocal1(vr);

#if 300 <= SUNDIALS_LIB_VERSION
    int mat_id = Int_val(vid);
    SUNMatrix smat;
    int content_is_value =
	(mat_id == MATRIX_ID_CUSTOM)
	|| (mat_id == MATRIX_ID_ARRAYDENSE)
	|| (mat_id == MATRIX_ID_ARRAYBAND);

    /* Create matrix */
    if (content_is_value) {
	smat = alloc_smat((void *)vcontent, vpayload, vcontext, true);
    } else {
	smat = alloc_smat(MAT_CONTENT(vcontent), vpayload, vcontext, false);
    }
    if (smat == NULL) caml_raise_out_of_memory();
#if 600 <= SUNDIALS_LIB_VERSION
    smat->sunctx = ML_CONTEXT(vcontext);
#endif

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
#if 500 <= SUNDIALS_LIB_VERSION
	smat->ops->matvecsetup = NULL;
#endif
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
#if 500 <= SUNDIALS_LIB_VERSION
	smat->ops->matvecsetup = NULL;
#endif
	smat->ops->matvec      = SUNMatMatvec_Band;
	smat->ops->space       = SUNMatSpace_Band;
	break;

    case MATRIX_ID_SPARSE_CSC:
    case MATRIX_ID_SPARSE_CSR:
	smat->ops->clone       = csmat_sparse_clone;      // ours
	smat->ops->destroy     = free_smat;  // ours (only called for c clones)
	smat->ops->getid       = SUNMatGetID_Sparse;
	smat->ops->zero        = SUNMatZero_Sparse;
	smat->ops->copy        = csmat_sparse_copy;	 // ours
	smat->ops->scaleadd    = csmat_sparse_scale_add;  // ours
	smat->ops->scaleaddi   = csmat_sparse_scale_addi; // ours
#if 500 <= SUNDIALS_LIB_VERSION
	smat->ops->matvecsetup = NULL;
#endif
	smat->ops->matvec      = SUNMatMatvec_Sparse;
	smat->ops->space       = SUNMatSpace_Sparse;
	break;

    case MATRIX_ID_CUSTOM:
    case MATRIX_ID_ARRAYDENSE:
    case MATRIX_ID_ARRAYBAND:
	smat->ops->clone       = csmat_custom_clone;
	smat->ops->destroy     = free_custom_smat;  // ours (only called for
						   //	    c clones)
	smat->ops->getid       = csmat_custom_getid;
	smat->ops->zero        = csmat_custom_zero;
	smat->ops->copy        = csmat_custom_copy;
	smat->ops->scaleadd    = csmat_custom_scale_add;
	smat->ops->scaleaddi   = csmat_custom_scale_addi;
	smat->ops->matvec      = csmat_custom_matvec;
#if 500 <= SUNDIALS_LIB_VERSION
	smat->ops->matvecsetup =
	    Bool_val(vhasmatvecsetup) ? csmat_custom_matvecsetup
				      : NULL;
#endif
	smat->ops->space       = csmat_custom_space;
	break;
    }

    // Setup the OCaml-side
    vr = caml_alloc_final(1,
	    content_is_value
		? &finalize_caml_custom_smat
		: &finalize_caml_smat,
	    1, 20);
    MAT_CVAL(vr) = smat;

#else // SUNDIALS_LIB_VERSION < 300
    vr = Val_unit;
#endif
    CAMLreturn(vr);
}

#if 650 <= SUNDIALS_LIB_VERSION
CAMLprim value sunml_matrix_wrap_any(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal4(vr, vmat, vid, vrawptr);

    // Create a new SUNMatrix value that shares its data with the
    // original one.

    SUNMatrix B = NULL;
    vmat = caml_alloc_tuple(RECORD_MAT_MATRIX_SIZE);
    Store_field(vmat, RECORD_MAT_MATRIX_PAYLOAD, MAT_BACKLINK(A));

    switch (SUNMatGetID(A)) {
	case SUNMATRIX_DENSE:
	case SUNMATRIX_BAND:
	case SUNMATRIX_SPARSE:
	    vrawptr = caml_alloc_final(1, finalize_caml_smat, 1, 20);
	    B = alloc_smat(A->content, MAT_BACKLINK(A),
			   MAT_CONTEXT(A), false);
	    break;

	case SUNMATRIX_CUSTOM:
	    vrawptr = caml_alloc_final(1, finalize_caml_custom_smat, 1, 20);
	    B = alloc_smat(MAT_CUSTOM_CONTENT(A), MAT_BACKLINK(A),
			   MAT_CONTEXT(A), true);
	    break;

	case SUNMATRIX_MAGMADENSE:
	case SUNMATRIX_ONEMKLDENSE:
	case SUNMATRIX_SLUNRLOC:
	case SUNMATRIX_CUSPARSE:
	case SUNMATRIX_GINKGO:
	case SUNMATRIX_KOKKOSDENSE:
	    caml_failwith("sunml_matrix_wrap_any: unexpected id");
    }
    MAT_CVAL(vrawptr) = B;
    Store_field(vmat, RECORD_MAT_MATRIX_RAWPTR, vrawptr);
    csmat_clone_ops(B, A);

    // Decode the fields to create an "any" matrix.
    enum mat_matrix_id_tag matid;
    switch (SUNMatGetID(B)) {
	case SUNMATRIX_DENSE:
	    Store_field(vmat, RECORD_MAT_MATRIX_ID, MATRIX_ID_DENSE);
	    Store_field(vmat, RECORD_MAT_MATRIX_MATOPS, OCAML_VALUE(dense_ops));
	    matid = MATRIX_ID_DENSE;
	    break;

	case SUNMATRIX_BAND:
	    Store_field(vmat, RECORD_MAT_MATRIX_ID, MATRIX_ID_BAND);
	    Store_field(vmat, RECORD_MAT_MATRIX_MATOPS, OCAML_VALUE(band_ops));
	    matid = MATRIX_ID_BAND;
	    break;

	case SUNMATRIX_SPARSE:
	    switch (SM_SPARSETYPE_S(B)) {
	    case CSC_MAT:
		matid = MATRIX_ID_SPARSE_CSC;
		break;
	    case CSR_MAT:
		matid = MATRIX_ID_SPARSE_CSR;
		break;
	    default:
		caml_failwith("sunml_matrix_wrap_any: unexpected sparse type");
	    }
	    Store_field(vmat, RECORD_MAT_MATRIX_ID, matid);
	    Store_field(vmat, RECORD_MAT_MATRIX_MATOPS, OCAML_VALUE(sparse_ops));
	    break;

	case SUNMATRIX_CUSTOM:
	    vid = MAT_OCAML_ID(B);
	    Store_field(vmat, RECORD_MAT_MATRIX_ID, vid);
	    Store_field(vmat, RECORD_MAT_MATRIX_MATOPS, MAT_OP_TABLE(B));
	    matid = Int_val(vid);
	    break;

	case SUNMATRIX_MAGMADENSE:
	case SUNMATRIX_ONEMKLDENSE:
	case SUNMATRIX_SLUNRLOC:
	case SUNMATRIX_CUSPARSE:
	case SUNMATRIX_GINKGO:
	case SUNMATRIX_KOKKOSDENSE:
	    caml_failwith("sunml_matrix_wrap_any: unexpected id");
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    vr = caml_alloc(1, matid);
#pragma GCC diagnostic pop
    Store_field(vr, 0, vmat);

    CAMLreturn(vr);
}
#endif

#if 300 <= SUNDIALS_LIB_VERSION

static SUNMatrix csmat_custom_clone(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(mlop, vcontentb);
    SUNMatrix B;

    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_CLONE);

    vcontentb = caml_callback_exn(mlop, MAT_BACKLINK(A));
    if (Is_exception_result (vcontentb))
	CAMLreturnT(SUNMatrix, NULL);

    B = alloc_smat(MAT_CUSTOM_CONTENT(A), vcontentb, MAT_CONTEXT(A), true);
    csmat_clone_ops(B, A);
#if 600 <= SUNDIALS_LIB_VERSION
    B->sunctx = A->sunctx;
#endif

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
    if (Is_exception_result (r)) {
	r = Extract_exception(r);
	CAMLreturnT(int, 1);
    }

    CAMLreturnT(int, 0);
}

static int csmat_custom_copy(SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_COPY);

    r = caml_callback2_exn(mlop, MAT_BACKLINK(A), MAT_BACKLINK(B));
    if (Is_exception_result (r)) {
	r = Extract_exception(r);
	CAMLreturnT(int, 1);
    }

    CAMLreturnT(int, 0);
}

static int csmat_custom_scale_add(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_SCALE_ADD);

    r = caml_callback3_exn(mlop, caml_copy_double(c), MAT_BACKLINK(A),
			   MAT_BACKLINK(B));
    if (Is_exception_result (r)) {
	r = Extract_exception(r);
	CAMLreturnT(int, 1);
    }

    CAMLreturnT(int, 0);
}

static int csmat_custom_scale_addi(sunrealtype c, SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_SCALE_ADDI);

    r = caml_callback2_exn(mlop, caml_copy_double(c), MAT_BACKLINK(A));
    if (Is_exception_result (r)) {
	r = Extract_exception(r);
	CAMLreturnT(int, 1);
    }

    CAMLreturnT(int, 0);
}

#if 500 <= SUNDIALS_LIB_VERSION
static int csmat_custom_matvecsetup(SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = Some_val (GET_OP(A, RECORD_MAT_MATRIXOPS_MATVEC_SETUP));

    r = caml_callback_exn(mlop, MAT_BACKLINK(A));
    if (Is_exception_result (r)) {
	r = Extract_exception(r);
	CAMLreturnT(int, 1);
    }

    CAMLreturnT(int, 0);
}
#endif

static int csmat_custom_matvec(SUNMatrix A, N_Vector x, N_Vector y)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_MATVEC);

    r = caml_callback3_exn(mlop, MAT_BACKLINK(A), NVEC_BACKLINK(x),
			   NVEC_BACKLINK(y));
    if (Is_exception_result (r)) {
	r = Extract_exception(r);
	CAMLreturnT(int, 1);
    }

    CAMLreturnT(int, 0);
}

static int csmat_custom_space(SUNMatrix A, long int *lenrw, long int *leniw)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(A, RECORD_MAT_MATRIXOPS_SPACE);

    r = caml_callback_exn(mlop, MAT_BACKLINK(A));
    if (Is_exception_result (r)) {
	r = Extract_exception(r);
	CAMLreturnT(int, 1);
    }

    *lenrw = Long_val(Field(r, 0));
    *leniw = Long_val(Field(r, 1));

    CAMLreturnT(int, 0);
}
#endif

CAMLprim void sunml_matrix_scale_add(value vc, value va, value vb)
{
    CAMLparam3(vc, va, vb);
#if 300 <= SUNDIALS_LIB_VERSION
    if (SUNMatScaleAdd(Double_val(vc), MAT_VAL(va), MAT_VAL(vb)))
	caml_failwith("SUNMatScaleAdd");
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_matrix_scale_addi(value vc, value va)
{
    CAMLparam2(vc, va);
#if 300 <= SUNDIALS_LIB_VERSION
    if (SUNMatScaleAddI(Double_val(vc), MAT_VAL(va)))
	caml_failwith("SUNMatScaleAddI");
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_matrix_matvecsetup(value va)
{
    CAMLparam1(va);
#if 500 <= SUNDIALS_LIB_VERSION
    SUNMatrix A = MAT_VAL(va);

    if (A->ops->matvecsetup) {
	if (SUNMatMatvecSetup(A))
	    caml_failwith("SUNMatMatvecsetup");
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_matrix_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);
#if 300 <= SUNDIALS_LIB_VERSION
    if (SUNMatMatvec(MAT_VAL(va), NVEC_VAL(vx), NVEC_VAL(vy)))
	caml_failwith("SUNMatMatvec");
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_matrix_zero(value va)
{
    CAMLparam1(va);
#if 300 <= SUNDIALS_LIB_VERSION
    if (SUNMatZero(MAT_VAL(va)))
	caml_failwith("SUNMatZero");
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_matrix_copy(value va, value vb)
{
    CAMLparam2(va, vb);
#if 300 <= SUNDIALS_LIB_VERSION
    if (SUNMatCopy(MAT_VAL(va), MAT_VAL(vb)))
	caml_failwith("SUNMatCopy");
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim value sunml_matrix_space(value va)
{
    CAMLparam1(va);
    CAMLlocal1(vr);
#if 300 <= SUNDIALS_LIB_VERSION
    long int lenrw, leniw;

    if (SUNMatSpace(MAT_VAL(va), &lenrw, &leniw))
	caml_failwith("SUNMatSpace");

    vr = caml_alloc_tuple(2);
    Store_field(vr, 0, Val_long(lenrw));
    Store_field(vr, 1, Val_long(leniw));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim void sunml_matrix_print_dense(value vm, value vfile)
{
    CAMLparam2(vm, vfile);
#if 300 <= SUNDIALS_LIB_VERSION
    SUNDenseMatrix_Print(MAT_VAL(vm), ML_CFILE(vfile));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_matrix_print_band(value vm, value vfile)
{
    CAMLparam2(vm, vfile);
#if 300 <= SUNDIALS_LIB_VERSION
    SUNBandMatrix_Print(MAT_VAL(vm), ML_CFILE(vfile));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_matrix_print_sparse(value vm, value vfile)
{
    CAMLparam2(vm, vfile);
#if 300 <= SUNDIALS_LIB_VERSION
    SUNSparseMatrix_Print(MAT_VAL(vm), ML_CFILE(vfile));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Array matrices
 */

CAMLprim value sunml_arraydensematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseScale(Double_val(vc), ARRAY2_ACOLS(va), m, n);
#else
    denseScale(Double_val(vc), ARRAY2_ACOLS(va), m, n);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraydensematrix_add_identity(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];

    if (m != n)
	caml_invalid_argument("matrix not square.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseAddIdentity(ARRAY2_ACOLS(va), m);
#else
    denseAddIdentity(ARRAY2_ACOLS(va), m);
#endif
    CAMLreturn (Val_unit);
}

#if 260 <= SUNDIALS_LIB_VERSION
CAMLprim value sunml_arraydensematrix_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);
    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vx) < n)
	caml_invalid_argument("x array too small.");
    if (ARRAY1_LEN(vy) < m)
	caml_invalid_argument("y array too small.");
#endif
#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx),
			  REAL_ARRAY(vy), m, n);
#else
    denseMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx), REAL_ARRAY(vy), m, n);
#endif
    CAMLreturn (Val_unit);
}
#else
CAMLprim value sunml_arraydensematrix_matvec(value va, value vx, value vy)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

CAMLprim value sunml_arraydensematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("pivot array too small.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    sundials_ml_index r = SUNDlsMat_denseGETRF(ARRAY2_ACOLS(va), m, n,
					       INDEX_ARRAY(vp));
#else
    sundials_ml_index r = denseGETRF(ARRAY2_ACOLS(va), m, n, INDEX_ARRAY(vp));
#endif

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_index(r));
    }
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraydensematrix_getrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);

    struct caml_ba_array *ba = ARRAY2_BA(va);
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

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp),
			 REAL_ARRAY(vb));
#else
    denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp), REAL_ARRAY(vb));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraydensematrix_getrs_off(value va, value vp,
					    value vb, value vboff)
{
    CAMLparam4(va, vp, vb, vboff);

    struct caml_ba_array *ba = ARRAY2_BA(va);
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

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp),
			 REAL_ARRAY(vb) + boff);
#else
    denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp), REAL_ARRAY(vb) + boff);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraydensematrix_potrf(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("matrix not square");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_densePOTRF(ARRAY2_ACOLS(va), m);
#else
    densePOTRF(ARRAY2_ACOLS(va), m);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraydensematrix_potrs(value va, value vb)
{
    CAMLparam2(va, vb);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("matrix not square.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("b is too small.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_densePOTRS(ARRAY2_ACOLS(va), m, REAL_ARRAY(vb));
#else
    densePOTRS(ARRAY2_ACOLS(va), m, REAL_ARRAY(vb));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraydensematrix_geqrf(value va, value vbeta, value vv)
{
    CAMLparam3(va, vbeta, vv);

    struct caml_ba_array *ba = ARRAY2_BA(va);
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

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseGEQRF(ARRAY2_ACOLS(va), m, n, REAL_ARRAY(vbeta),
			 REAL_ARRAY(vv));
#else
    denseGEQRF(ARRAY2_ACOLS(va), m, n, REAL_ARRAY(vbeta), REAL_ARRAY(vv));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraydensematrix_ormqr(value va, value vormqr)
{
    CAMLparam2(va, vormqr);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    sunrealtype *beta = REAL_ARRAY(Field(vormqr, 0));
    sunrealtype *vv   = REAL_ARRAY(Field(vormqr, 1));
    sunrealtype *vw   = REAL_ARRAY(Field(vormqr, 2));
    sunrealtype *work = REAL_ARRAY(Field(vormqr, 3));

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

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseORMQR(ARRAY2_ACOLS(va), m, n, beta, vv, vw, work);
#else
    denseORMQR(ARRAY2_ACOLS(va), m, n, beta, vv, vw, work);
#endif
    CAMLreturn (Val_unit);
}

/* Array Band matrix functions */

CAMLprim value sunml_arraybandmatrix_copy(value va, value vb, value vsizes)
{
    CAMLparam3(va, vb, vsizes);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat am = ba->dim[0];

    sundials_ml_index a_smu  = Index_val(Field(vsizes, 0));
    sundials_ml_index b_smu  = Index_val(Field(vsizes, 1));
    sundials_ml_index copymu = Index_val(Field(vsizes, 2));
    sundials_ml_index copyml = Index_val(Field(vsizes, 3));

#if SUNDIALS_ML_SAFE == 1
    intnat an = ba->dim[1];
    struct caml_ba_array *bb = ARRAY2_BA(vb);

    intnat bm = bb->dim[0];
    intnat bn = bb->dim[1];

    if (an < copymu + copyml + 1)
	caml_invalid_argument("source matrix too small.");
    if (bn < copymu + copyml + 1)
	caml_invalid_argument("destination matrix too small.");
    if (am != bm)
	caml_invalid_argument("matrix sizes differ.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_bandCopy(ARRAY2_ACOLS(va), ARRAY2_ACOLS(vb), am, a_smu, b_smu,
		       copymu, copyml);
#else
    bandCopy(ARRAY2_ACOLS(va), ARRAY2_ACOLS(vb), am, a_smu, b_smu,
	     copymu, copyml);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraybandmatrix_scale(value vc, value va, value vsizes)
{
    CAMLparam3(vc, va, vsizes);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[0];

    sundials_ml_index smu = Index_val(Field(vsizes, 0));
    sundials_ml_index mu  = Index_val(Field(vsizes, 1));
    sundials_ml_index ml  = Index_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("matrix badly sized.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_bandScale(Double_val(vc), ARRAY2_ACOLS(va), m, mu, ml, smu);
#else
    bandScale(Double_val(vc), ARRAY2_ACOLS(va), m, mu, ml, smu);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraybandmatrix_add_identity(value vsmu, value va)
{
    CAMLparam2(vsmu, va);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[0];
    sundials_ml_index smu = Index_val(vsmu);

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n <= smu)
	caml_invalid_argument("matrix badly sized.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_bandAddIdentity(ARRAY2_ACOLS(va), m, smu);
#else
    bandAddIdentity(ARRAY2_ACOLS(va), m, smu);
#endif
    CAMLreturn (Val_unit);
}

#if 260 <= SUNDIALS_LIB_VERSION
CAMLprim value sunml_arraybandmatrix_matvec(value va, value vsizes,
					value vx, value vy)
{
    CAMLparam4(va, vsizes, vx, vy);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[0];

    sundials_ml_index smu = Index_val(Field(vsizes, 0));
    sundials_ml_index mu  = Index_val(Field(vsizes, 1));
    sundials_ml_index ml  = Index_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("matrix badly sized.");

    if (ARRAY1_LEN(vx) < m)
	caml_invalid_argument("x array too small.");
    if (ARRAY1_LEN(vy) < m)
	caml_invalid_argument("y array too small.");
#endif
#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_bandMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx), REAL_ARRAY(vy),
			 m, mu, ml, smu);
#else
    bandMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx), REAL_ARRAY(vy),
	       m, mu, ml, smu);
#endif
    CAMLreturn (Val_unit);
}
#else
CAMLprim value sunml_arraybandmatrix_matvec(value va, value vx, value vy,
					value vsizes)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

CAMLprim value sunml_arraybandmatrix_gbtrf(value va, value vsizes, value vp)
{
    CAMLparam3(va, vsizes, vp);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[0];

    sundials_ml_index smu = Index_val(Field(vsizes, 0));
    sundials_ml_index mu  = Index_val(Field(vsizes, 1));
    sundials_ml_index ml  = Index_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("p is too small.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_bandGBTRF(ARRAY2_ACOLS(va), m, mu, ml, smu, INDEX_ARRAY(vp));
#else
    bandGBTRF(ARRAY2_ACOLS(va), m, mu, ml, smu, INDEX_ARRAY(vp));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arraybandmatrix_gbtrs(value va, value vsizes, value vp, value vb)
{
    CAMLparam4(va, vsizes, vp, vb);

    struct caml_ba_array *ba = ARRAY2_BA(va);
    intnat m = ba->dim[0];

    sundials_ml_index smu = Index_val(Field(vsizes, 0));
    sundials_ml_index ml  = Index_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < smu + ml + 1)
	caml_invalid_argument("matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("p is too small.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("b is too small.");
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_bandGBTRS(ARRAY2_ACOLS(va), m, smu, ml, INDEX_ARRAY(vp),
			REAL_ARRAY(vb));
#else
    bandGBTRS(ARRAY2_ACOLS(va), m, smu, ml, INDEX_ARRAY(vp), REAL_ARRAY(vb));
#endif
    CAMLreturn (Val_unit);
}

