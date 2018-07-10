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

#ifndef SUNDIALS_ML_KLU
CAMLprim value sunml_ida_klu_init (value vida_mem, value vformat,
			       value vneqs, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_ida_klu_set_ordering (value vida_mem, value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_ida_klu_reinit (value vida_mem, value vn, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_ida_klu_get_num_jac_evals(value vida_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "ida_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#ifdef SUNDIALSML_WITHSENS
/* IDAS (with sensitivity) */

#include <idas/idas.h>

#if SUNDIALS_LIB_VERSION < 300
#include <idas/idas_sparse.h>
#include <idas/idas_klu.h>
#endif

#else
/* IDA (without sensitivity) */

#include <ida/ida.h>

#if SUNDIALS_LIB_VERSION < 300
#include <ida/ida_sparse.h>
#include <ida/ida_klu.h>
#endif

#endif

#if SUNDIALS_LIB_VERSION < 300
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
	Store_some(smat, c_matrix_sparse_wrap(jac));
	Store_field(cb, 1, smat);

	args[1] = Some_val(smat);
    } else {
	args[1] = Some_val(smat);
	sunml_matrix_sparse_rewrap(args[1]);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

CAMLprim value sunml_ida_klu_init (value vida_mem, value vformat,
			       value vneqs, value vnnz)
{
    CAMLparam4(vida_mem, vformat, vneqs, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

#if SUNDIALS_LIB_VERSION >= 270
    flag = IDAKLU (ida_mem, Int_val(vneqs), Int_val(vnnz),
		   MAT_FROM_SFORMAT(vformat));
#else
    flag = IDAKLU (ida_mem, Int_val(vneqs), Int_val(vnnz));
#endif
    CHECK_FLAG ("IDAKLU", flag);
    flag = IDASlsSetSparseJacFn(ida_mem, jacfn);
    CHECK_FLAG("IDASlsSetSparseJacFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_klu_set_ordering (value vida_mem, value vordering)
{
    CAMLparam2(vida_mem, vordering);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    int flag = IDAKLUSetOrdering (ida_mem, Int_val(vordering));
    CHECK_FLAG ("IDAKLUSetOrdering", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_klu_reinit (value vida_mem, value vn, value vnnz)
{
    CAMLparam3(vida_mem, vn, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int nnz = Int_val(vnnz);

    int flag = IDAKLUReInit (ida_mem, Int_val(vn), nnz, (nnz > 0) ? 1 : 2);
    CHECK_FLAG ("IDAKLUReInit", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_ida_klu_get_num_jac_evals(value vida_mem)
{
    CAMLparam1(vida_mem);
    long int r = 0;
#if SUNDIALS_LIB_VERSION < 300
    int flag = IDASlsGetNumJacEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASlsGetNumJacEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

#endif
