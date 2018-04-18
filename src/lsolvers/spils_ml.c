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

#include "../sundials/sundials_ml.h"

#include <sundials/sundials_config.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_nvector.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include "../lsolvers/matrix_ml.h"
#include "../nvectors/nvector_ml.h"

CAMLprim value c_spils_modified_gs(value vv, value vh, value vk, value vp)
{
    CAMLparam4(vv, vh, vk, vp);

    int p = Int_val(vp);
    int k = Int_val(vk);
    int i;
    realtype new_vk_norm;
    N_Vector* v;

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < k) || (hn < k))
	caml_invalid_argument("Spils.modified_gs: h is too small.");
    if (Wosize_val (vv) < k)
	caml_invalid_argument("Spils.modified_gs: v is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_out_of_memory();
    for (i = 0; i <= p; ++i) {
	v[i] = NVEC_VAL(Field(vv, k - p + i));
    }

    ModifiedGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm);

    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value c_spils_classical_gs(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal3(vv, vh, vs);

    int k = Int_val(Field(vargs, 2));
    int p = Int_val(Field(vargs, 3));
    N_Vector temp;

    vv = Field(vargs, 0);
    vh = Field(vargs, 1);
    vs = Field(vargs, 5);

    int i;
    realtype new_vk_norm;
    N_Vector* v;

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < k) || (hn < k))
	caml_invalid_argument("Spils.classical_gs: h is too small.");
    if (Wosize_val (vv) < k)
	caml_invalid_argument("Spils.classical_gs: v is too small.");
    if (ARRAY1_LEN(vs) < k)
	caml_invalid_argument("Spils.classical_gs: s is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_out_of_memory();
    for (i = 0; i <= p; ++i) {
	v[i] = NVEC_VAL(Field(vv, k - p + i));
    }

    temp = NVEC_VAL(Field(vargs, 4));
    ClassicalGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm,
	        temp, REAL_ARRAY(vs));

    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value c_spils_qr_fact(value vh, value vq, value vnewjob)
{
    CAMLparam3(vh, vq, vnewjob);
    int r;
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hn = bh->dim[0];

#if SUNDIALS_ML_SAFE == 1
    intnat hm = bh->dim[1];

    if ((hm < hn + 1))
	caml_invalid_argument("Spils.qr_fact: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * hn)
	caml_invalid_argument("Spils.qr_fact: q is too small.");
#endif

    r = QRfact(hn, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), Bool_val(vnewjob));

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

CAMLprim value c_spils_qr_sol(value vh, value vq, value vb)
{
    CAMLparam3(vh, vq, vb);
    int r;
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hn = bh->dim[0];

#if SUNDIALS_ML_SAFE == 1
    intnat hm = bh->dim[1];

    if (hm < hn + 1)
	caml_invalid_argument("Spils.qr_sol: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * hn)
	caml_invalid_argument("Spils.qr_sol: q is too small.");
    if (ARRAY1_LEN(vb) < hn + 1)
	caml_invalid_argument("Spils.qr_sol: b is too small.");
#endif

    r = QRsol(hn, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), REAL_ARRAY(vb));

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

