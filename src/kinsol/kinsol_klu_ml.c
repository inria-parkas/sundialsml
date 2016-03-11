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
#include "kinsol_ml.h"
#include "../lsolvers/sls_ml.h"

#ifndef SUNDIALS_ML_KLU
CAMLprim value c_kinsol_klu_init (value vkin_mem, value vneqs, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_kinsol_klu_set_ordering (value vkin_mem, value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_kinsol_klu_reinit (value vkin_mem, value vn, value vnnz,
				   value vrealloc)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_kinsol_klu_get_num_jac_evals(value vkin_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_sparse.h>
#include <kinsol/kinsol_klu.h>

enum kinsol_klu_ordering_tag {
  VARIANT_CVODE_KLU_AMD     = 0,
  VARIANT_CVODE_KLU_COLAMD  = 1,
  VARIANT_CVODE_KLU_NATURAL = 2,
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
	Store_some(smat, c_sls_sparse_wrap(Jac, 0));
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

CAMLprim value c_kinsol_klu_init (value vkin_mem, value vneqs, value vnnz)
{
    CAMLparam3(vkin_mem, vneqs, vnnz);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINKLU (kin_mem, Int_val(vneqs), Int_val(vnnz));
    CHECK_FLAG ("KINKLU", flag);
    flag = KINSlsSetSparseJacFn(kin_mem, jacfn);
    CHECK_FLAG("KINSlsSetSparseJacFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_kinsol_klu_set_ordering (value vkin_mem, value vordering)
{
    CAMLparam2(vkin_mem, vordering);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);

    int flag = KINKLUSetOrdering (kin_mem, Int_val(vordering));
    CHECK_FLAG ("KINKLUSetOrdering", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_kinsol_klu_reinit (value vkin_mem, value vn, value vnnz,
				   value vrealloc)
{
    CAMLparam4(vkin_mem, vn, vnnz, vrealloc);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);

    int reinit_type;
    if (Bool_val(vrealloc)) {
	reinit_type = 1;
    } else {
	reinit_type = 2;
    }

    int flag = KINKLUReInit (kin_mem, Int_val(vn), Int_val(vnnz), reinit_type);
    CHECK_FLAG ("KINKLUReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_kinsol_klu_get_num_jac_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);

    long int r;
    int flag = KINSlsGetNumJacEvals(kin_mem, &r);
    CHECK_FLAG("KINSlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

#endif
