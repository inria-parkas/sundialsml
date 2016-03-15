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
#include "ida_ml.h"
#include "../lsolvers/sls_ml.h"

#ifndef SUNDIALS_ML_KLU
CAMLprim value c_ida_klu_init (value vida_mem, value vneqs, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_ida_klu_set_ordering (value vida_mem, value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_ida_klu_reinit (value vida_mem, value vn, value vnnz,
				   value vrealloc)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_ida_klu_get_num_jac_evals(value vida_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef SUNDIALSML_WITHSENS
/* IDAS (with sensitivity) */

#include <idas/idas.h>
#include <idas/idas_sparse.h>
#include <idas/idas_klu.h>

#else
/* IDA (without sensitivity) */

#include <ida/ida.h>
#include <ida/ida_sparse.h>
#include <ida/ida_klu.h>

#endif

enum ida_klu_ordering_tag {
  VARIANT_IDA_KLU_AMD     = 0,
  VARIANT_IDA_KLU_COLAMD  = 1,
  VARIANT_IDA_KLU_NATURAL = 2,
};

static int jacfn (realtype t, realtype coef,
		  N_Vector y, N_Vector yp, N_Vector res,
		  SlsMat jac, void *user_data,
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, smat);

    WEAK_DEREF (session, *(value*)user_data);
    args[0] = ida_make_jac_arg (t, coef, y, yp, res,
				ida_make_triple_tmp (tmp1, tmp2, tmp3));

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    smat = Field(cb, 1);

    if (smat == Val_none) {
	Store_some(smat, c_sls_sparse_wrap(jac, 0));
	Store_field(cb, 1, smat);

	args[1] = Some_val(smat);
    } else {
	args[1] = Some_val(smat);
	c_sparsematrix_realloc(args[1], 0);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

CAMLprim value c_ida_klu_init (value vida_mem, value vneqs, value vnnz)
{
    CAMLparam3(vida_mem, vneqs, vnnz);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDAKLU (ida_mem, Int_val(vneqs), Int_val(vnnz));
    CHECK_FLAG ("IDAKLU", flag);
    flag = IDASlsSetSparseJacFn(ida_mem, jacfn);
    CHECK_FLAG("IDASlsSetSparseJacFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_klu_set_ordering (value vida_mem, value vordering)
{
    CAMLparam2(vida_mem, vordering);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    int flag = IDAKLUSetOrdering (ida_mem, Int_val(vordering));
    CHECK_FLAG ("IDAKLUSetOrdering", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_klu_reinit (value vida_mem, value vn, value vnnz,
				   value vrealloc)
{
    CAMLparam4(vida_mem, vn, vnnz, vrealloc);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    int reinit_type;
    if (Bool_val(vrealloc)) {
	reinit_type = 1;
    } else {
	reinit_type = 2;
    }

    int flag = IDAKLUReInit (ida_mem, Int_val(vn), Int_val(vnnz), reinit_type);
    CHECK_FLAG ("IDAKLUReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_klu_get_num_jac_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDASlsGetNumJacEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

#endif
