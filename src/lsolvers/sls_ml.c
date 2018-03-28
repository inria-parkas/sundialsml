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

#include "../config.h"

#include <stdio.h>

#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include "../sundials/sundials_ml.h"
#include "../lsolvers/sls_ml.h"
#include "../lsolvers/dls_ml.h"

static void finalize_slsmat(value va)
{
#if SUNDIALS_LIB_VERSION >= 270
    SparseDestroyMat(SLSMAT(va));
#else
    DestroySparseMat(SLSMAT(va));
#endif
}

CAMLprim value c_sls_sparse_wrap(SlsMat a, int finalize, value vformat)
{
    CAMLparam1(vformat);
    CAMLlocal4(vidxptrs, vidxvals, vdata, vv);
    CAMLlocal1(vr);

#if SUNDIALS_LIB_VERSION >= 270
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INT,   1, a->indexptrs, a->N + 1);
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INT,   1, a->indexvals, a->NNZ);
#else
    vidxptrs = caml_ba_alloc_dims(BIGARRAY_INT,   1, a->colptrs, a->N + 1);
    vidxvals = caml_ba_alloc_dims(BIGARRAY_INT,   1, a->rowvals, a->NNZ);
#endif
    vdata    = caml_ba_alloc_dims(BIGARRAY_FLOAT, 1, a->data,    a->NNZ);

    /* a SlsMat is a pointer to a struct _SlsMat */
    vv = caml_alloc_final(2, finalize ? &finalize_slsmat : NULL, 1, 20);
    SLSMAT(vv) = a;

    vr = caml_alloc_tuple(RECORD_SLS_SPARSEMATRIX_SIZE);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_IDXPTRS, vidxptrs);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_IDXVALS, vidxvals);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_DATA,    vdata);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_SLSMAT,  vv);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_SFORMAT, vformat);
    Store_field(vr, RECORD_SLS_SPARSEMATRIX_VALID,   Val_bool(1));

    CAMLreturn(vr);
}

CAMLprim value c_sparsematrix_realloc(value vma, value vnnz)
{
    CAMLparam1(vma);
    CAMLlocal2(vidxvals, vdata);

    SlsMat a = SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT));

    struct caml_ba_array *ba_idxptrs =
	Caml_ba_array_val(Field(vma, RECORD_SLS_SPARSEMATRIX_IDXPTRS));
    struct caml_ba_array *ba_idxvals =
	Caml_ba_array_val(Field(vma, RECORD_SLS_SPARSEMATRIX_IDXVALS));
    struct caml_ba_array *ba_data =
	Caml_ba_array_val(Field(vma, RECORD_SLS_SPARSEMATRIX_DATA));

    int nnz = Int_val(vnnz);
    if (nnz > 0) {
#if SUNDIALS_LIB_VERSION >= 270
	a->indexvals = realloc(a->indexvals, nnz*sizeof(int));
	if (a->indexvals == NULL)
	    caml_raise_out_of_memory();
#else
	a->rowvals = realloc(a->rowvals, nnz*sizeof(int));
	if (a->rowvals == NULL)
	    caml_raise_out_of_memory();
#endif
	a->data = realloc(a->data, nnz*sizeof(realtype));
	if (a->data == NULL)
	    caml_raise_out_of_memory();
	a->NNZ = nnz;
    } else {
#if SUNDIALS_LIB_VERSION >= 270
	ba_idxptrs->data = a->indexptrs;
#else
	ba_idxptrs->data = a->colptrs;
#endif
    }

#if SUNDIALS_LIB_VERSION >= 270
    ba_idxvals->data = a->indexvals;
#else
    ba_idxvals->data = a->rowvals;
#endif
    ba_idxvals->dim[0] = a->NNZ;

    ba_data->data = a->data;
    ba_data->dim[0] = a->NNZ;

    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_new_sparse_mat(value vm, value vn, value vnnz,
					     value vformat)
{
    CAMLparam4(vm, vn, vnnz, vformat);

    int m = Int_val(vm);
    int n = Int_val(vn);
    int nnz = Int_val(vnnz);

#if SUNDIALS_LIB_VERSION >= 270
    SlsMat a = SparseNewMat(m, n, nnz, Int_val(vformat));
#else
    SlsMat a = NewSparseMat(m, n, nnz);
#endif
    if (a == NULL)
	caml_raise_out_of_memory();

    CAMLreturn(c_sls_sparse_wrap(a, 1, vformat));
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

CAMLprim value c_sparsematrix_print_mat(value vlogfile, value va)
{
    CAMLparam2(vlogfile, va);
#if SUNDIALS_LIB_VERSION >= 270
    SparsePrintMat(SLSMAT(va), ML_CFILE(vlogfile));
#else
    if (ML_CFILE(vlogfile) != stdout)
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    PrintSparseMat(SLSMAT(va));
#endif
    fflush(ML_CFILE(vlogfile));
    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_set_to_zero(value va)
{
    CAMLparam1(va);
#if SUNDIALS_LIB_VERSION >= 270
    SparseSetMatToZero(SLSMAT(va));
#else
    SlsSetToZero(SLSMAT(va));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_convert_dls(value vformat, value va)
{
    CAMLparam2(vformat, va);
    CAMLlocal1(vr);

    DlsMat da = DLSMAT(va);
#if SUNDIALS_LIB_VERSION >= 270
    SlsMat ma = SparseFromDenseMat(da, Int_val(vformat));
#else
    SlsMat ma = SlsConvertDls(da);
#endif

    CAMLreturn(c_sls_sparse_wrap(ma, 1, vformat));
}

CAMLprim value c_sparsematrix_add_identity(value vma)
{
    CAMLparam1(vma);

#if SUNDIALS_LIB_VERSION >= 270
    SparseAddIdentityMat(SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT)));
#else
    AddIdentitySparseMat(SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT)));
#endif
    c_sparsematrix_realloc(vma, Val_int(0));

    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_copy(value va, value vmb)
{
    CAMLparam2(va, vmb);

#if SUNDIALS_LIB_VERSION >= 270
    SparseCopyMat(SLSMAT(va),
		  SLSMAT(Field(vmb, RECORD_SLS_SPARSEMATRIX_SLSMAT)));
#else
    CopySparseMat(SLSMAT(va),
		  SLSMAT(Field(vmb, RECORD_SLS_SPARSEMATRIX_SLSMAT)));
#endif
    c_sparsematrix_realloc(vmb, Val_int(0));

    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);
#if SUNDIALS_LIB_VERSION >= 270
    SparseScaleMat(Double_val(vc), SLSMAT(va));
#else
    ScaleSparseMat(Double_val(vc), SLSMAT(va));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_sparsematrix_add(value vma, value vb)
{
    CAMLparam2(vma, vb);

#if SUNDIALS_LIB_VERSION >= 270
    SparseAddMat(SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT)), SLSMAT(vb));
#else
    SlsAddMat(SLSMAT(Field(vma, RECORD_SLS_SPARSEMATRIX_SLSMAT)), SLSMAT(vb));
#endif
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

#if SUNDIALS_LIB_VERSION >= 270
    SparseMatvec(a, REAL_ARRAY(vx), REAL_ARRAY(vy));
#else
    SlsMatvec(a, REAL_ARRAY(vx), REAL_ARRAY(vy));
#endif
    /* assert (r == 0) */

    CAMLreturn (Val_unit);
}

