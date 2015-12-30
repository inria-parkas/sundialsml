/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include "config.h"

#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

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
#include "sls_ml.h"
#include "dls_ml.h"

static void finalize_slsmat(value va)
{
    DestroySparseMat(SLSMAT(va));
}

CAMLprim value c_sls_sparse_wrap(SlsMat a, int finalize)
{
    CAMLparam0();
    CAMLlocal4(vcolptrs, vrowvals, vdata, vv);
    CAMLlocal1(vr);

    mlsize_t approx_size = (a->M + a->N) * sizeof(long int)
			    + a->NNZ * sizeof(realtype) + 5;

    vcolptrs = caml_ba_alloc_dims(BIGARRAY_INT,   1, a->colptrs, a->N + 1);
    vrowvals = caml_ba_alloc_dims(BIGARRAY_INT,   1, a->rowvals, a->NNZ);
    vdata    = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, a->data,    a->NNZ);

    /* a SlsMat is a pointer to a struct _SlsMat */
    vv = caml_alloc_final(2, finalize ? &finalize_slsmat : NULL,
			  approx_size, approx_size * 20);
    SLSMAT(vv) = a;

    vr = caml_alloc_tuple(RECORD_SLS_SPARSEMATRIX_SIZE);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_COLPTRS, vcolptrs);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_ROWVALS, vrowvals);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_DATA,    vdata);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_SLSMAT,  vv);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_VALID,   Val_bool(1));

    CAMLreturn(vr);
}

CAMLprim value c_sparsematrix_realloc(value vma, value vnnz)
{
    CAMLparam1(vma);
    CAMLlocal2(vrowvals, vdata);

    SlsMat a = SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT));

    struct caml_ba_array *ba_colptrs =
	Caml_ba_array_val(Field(vma, RECORD_SLS_SPARSEMATRIX_COLPTRS));
    struct caml_ba_array *ba_rowvals =
	Caml_ba_array_val(Field(vma, RECORD_SLS_SPARSEMATRIX_ROWVALS));
    struct caml_ba_array *ba_data =
	Caml_ba_array_val(Field(vma, RECORD_SLS_SPARSEMATRIX_DATA));

    int nnz = Int_val(vnnz);
    if (nnz > 0) {
	a->rowvals = realloc(a->rowvals, nnz*sizeof(int));
	a->data = realloc(a->data, nnz*sizeof(realtype));
	if (a->rowvals == NULL || a->data == NULL)
	    caml_raise_out_of_memory();
	a->NNZ = nnz;
    } else {
	ba_colptrs->data = a->colptrs;
    }

    ba_rowvals->data = a->rowvals;
    ba_rowvals->dim[0] = a->NNZ;

    ba_data->data = a->data;
    ba_data->dim[0] = a->NNZ;

    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_new_sparse_mat(value vm, value vn, value vnnz)
{
    CAMLparam3(vm, vn, vnnz);

    int m = Int_val(vm);
    int n = Int_val(vn);
    int nnz = Int_val(vn);

    SlsMat a = NewSparseMat(m, n, nnz);
    if (a == NULL)
	caml_raise_out_of_memory();

    CAMLreturn(c_sls_sparse_wrap(a, 1));
}

CAMLprim value c_sparsematrix_size(value va)
{
    CAMLparam1(va);
    CAMLlocal1(vr);

    SlsMat ma = SLSMAT(va);
    vr = caml_alloc_tuple(3);
    Store_field(vr, 0, Val_int(ma->M));
    Store_field(vr, 1, Val_int(ma->N));
    Store_field(vr, 2, Val_int(ma->NNZ));

    CAMLreturn(vr);
}

CAMLprim value c_sparsematrix_print_mat(value va)
{
    CAMLparam1(va);
    PrintSparseMat(SLSMAT(va));
    fflush(stdout);
    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_set_to_zero(value va)
{
    CAMLparam1(va);
    SlsSetToZero(SLSMAT(va));
    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_convert_dls(value va)
{
    CAMLparam1(va);
    CAMLlocal1(vr);

    DlsMat da = DLSMAT(va);
    SlsMat ma = SlsConvertDls(da);

    CAMLreturn(c_sls_sparse_wrap(ma, 1));
}

CAMLprim value c_sparsematrix_add_identity(value vma)
{
    CAMLparam1(vma);

    AddIdentitySparseMat(SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT)));
    c_sparsematrix_realloc(vma, Val_int(0));

    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_copy(value va, value vmb)
{
    CAMLparam2(va, vmb);

    CopySparseMat(SLSMAT(va),
		  SLSMAT(Field(vmb, RECORD_SLS_SPARSEMATRIX_SLSMAT)));
    c_sparsematrix_realloc(vmb, Val_int(0));

    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
    ScaleSparseMat(Double_val(vc), SLSMAT(va));
    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_add(value vma, value vb)
{
    CAMLparam2(vma, vb);

    SlsAddMat(SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT)), SLSMAT(vb));
    c_sparsematrix_realloc(vma, Val_int(0));

    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);

    SlsMat a = SLSMAT(va);

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vx) != a->N)
	caml_invalid_argument("SparseMatrix.matvec: x has wrong size.");
    if (ARRAY1_LEN(vy) != a->M)
	caml_invalid_argument("SparseMatrix.matvec: y has wrong size.");
#endif

    SlsMatvec(a, REAL_ARRAY(vx), REAL_ARRAY(vy));
    /* assert (r == 0) */

    CAMLreturn (Val_unit);
}

