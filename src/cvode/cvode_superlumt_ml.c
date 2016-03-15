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
CAMLprim value c_cvode_superlumt_init (value vcvode_mem, value vneqs,
				       value vnnz, value vnthreads)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_cvode_superlumt_set_ordering (value vcvode_mem, value vorder)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_cvode_superlumt_get_num_jac_evals(value vcvode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "cvode_ml.h"
#include "../lsolvers/sls_ml.h"

#ifdef SUNDIALSML_WITHSENS
/* CVODES (with sensitivity) */

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_superlumt.h>

#else
/* CVODE (without sensitivity) */

#include <cvode/cvode.h>
#include <cvode/cvode_sparse.h>
#include <cvode/cvode_superlumt.h>

#endif

enum cvode_superlumt_ordering_tag {
  VARIANT_CVODE_SUPERLUMT_NATURAL    = 0,
  VARIANT_CVODE_SUPERLUMT_MINDEGPROD = 1,
  VARIANT_CVODE_SUPERLUMT_MINDEGSUM  = 2,
  VARIANT_CVODE_SUPERLUMT_COLAMD     = 3,
};

static int jacfn(
	realtype t,
	N_Vector y,
	N_Vector fy,
	SlsMat Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, smat);

    WEAK_DEREF (session, *(value*)user_data);
    args[0] = cvode_make_jac_arg (t, y, fy,
				  cvode_make_triple_tmp (tmp1, tmp2, tmp3));

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
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

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

CAMLprim value c_cvode_superlumt_init (value vcvode_mem, value vneqs,
				       value vnnz, value vnthreads)
{
    CAMLparam4(vcvode_mem, vneqs, vnnz, vnthreads);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVSuperLUMT (cvode_mem, Int_val(vnthreads), Int_val(vneqs),
			Int_val(vnnz));
    CHECK_FLAG ("CVSuperLUMT", flag);
    flag = CVSlsSetSparseJacFn(cvode_mem, jacfn);
    CHECK_FLAG("CVSlsSetSparseJacFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_cvode_superlumt_set_ordering (value vcvode_mem, value vorder)
{
    CAMLparam2(vcvode_mem, vorder);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    int flag = CVSuperLUMTSetOrdering (cvode_mem, Int_val(vorder));
    CHECK_FLAG ("CVSuperLUMTSetOrdering", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_cvode_superlumt_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    long int r;
    int flag = CVSlsGetNumJacEvals(cvode_mem, &r);
    CHECK_FLAG("CVSlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

#endif
