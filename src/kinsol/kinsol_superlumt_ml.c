/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2015 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>

#include "../sundials/sundials_ml.h"

#ifndef SUNDIALS_ML_SUPERLUMT
CAMLprim value c_kinsol_superlumt_init (value vkin_mem, value vneqs,
				        value vnnz, value vnthreads)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_kinsol_superlumt_set_ordering (value vkin_mem, value vorder)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_kinsol_superlumt_get_num_jac_evals(value vkin_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "kinsol_ml.h"
#include "../lsolvers/sls_ml.h"

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_sparse.h>
#include <kinsol/kinsol_superlumt.h>

enum kinsol_superlumt_ordering_tag {
  VARIANT_KINSOL_SUPERLUMT_NATURAL    = 0,
  VARIANT_KINSOL_SUPERLUMT_MINDEGPROD = 1,
  VARIANT_KINSOL_SUPERLUMT_MINDEGSUM  = 2,
  VARIANT_KINSOL_SUPERLUMT_COLAMD     = 3,
};

static int jacfn(
	N_Vector u,
	N_Vector fu,	     
	SlsMat Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, smat);

    WEAK_DEREF (session, *(value*)user_data);
    args[0] = kinsol_make_jac_arg(u, fu, kinsol_make_double_tmp(tmp1, tmp2));

    cb = KINSOL_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    smat = Field(cb, 1);
    if (smat == Val_none) {
#if SUNDIALS_LIB_VERSION >= 270
	Store_some(smat, c_sls_sparse_wrap(Jac, 0, Val_int(Jac->sparsetype)));
#else
	Store_some(smat, c_sls_sparse_wrap(Jac, 0, Val_int(0)));
#endif
	Store_field(cb, 1, smat);

	args[1] = Some_val(smat);
    } else {
	args[1] = Some_val(smat);
	c_sparsematrix_realloc(args[1], 0);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, UNRECOVERABLE));
}

CAMLprim value c_kinsol_superlumt_init (value vkin_mem, value vneqs,
				        value vnnz, value vnthreads)
{
    CAMLparam4(vkin_mem, vneqs, vnnz, vnthreads);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSuperLUMT (kin_mem, Int_val(vnthreads), Int_val(vneqs),
			 Int_val(vnnz));
    CHECK_FLAG ("KINSuperLUMT", flag);
    flag = KINSlsSetSparseJacFn(kin_mem, jacfn);
    CHECK_FLAG("KINSlsSetSparseJacFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_kinsol_superlumt_set_ordering (value vkin_mem, value vorder)
{
    CAMLparam2(vkin_mem, vorder);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);

    int flag = KINSuperLUMTSetOrdering (kin_mem, Int_val(vorder));
    CHECK_FLAG ("KINSuperLUMTSetOrdering", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_kinsol_superlumt_get_num_jac_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);

    long int r;
    int flag = KINSlsGetNumJacEvals(kin_mem, &r);
    CHECK_FLAG("KINSlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

#endif
