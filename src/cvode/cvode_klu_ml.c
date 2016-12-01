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
CAMLprim value c_cvode_klu_init (value vcvode_mem, value vformat,
			 	 value vneqs, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_cvode_klu_set_ordering (value vcvode_mem, value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_cvode_klu_reinit (value vcvode_mem, value vn, value vnnz,
				   value vrealloc)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_cvode_klu_get_num_jac_evals(value vcvode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "cvode_ml.h"
#include "../lsolvers/sls_ml.h"

#ifdef SUNDIALSML_WITHSENS
/* CVODES (with sensitivity) */

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_klu.h>

#else
/* CVODE (without sensitivity) */

#include <cvode/cvode.h>
#include <cvode/cvode_sparse.h>
#include <cvode/cvode_klu.h>

#endif

enum cvode_klu_ordering_tag {
  VARIANT_CVODE_KLU_AMD     = 0,
  VARIANT_CVODE_KLU_COLAMD  = 1,
  VARIANT_CVODE_KLU_NATURAL = 2,
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

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

CAMLprim value c_cvode_klu_init (value vcvode_mem, value vformat,
				 value vneqs, value vnnz)
{
    CAMLparam4(vcvode_mem, vformat, vneqs, vnnz);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

#if SUNDIALS_LIB_VERSION >= 270
    flag = CVKLU (cvode_mem, Int_val(vneqs), Int_val(vnnz), Int_val(vformat));
#else
    flag = CVKLU (cvode_mem, Int_val(vneqs), Int_val(vnnz));
#endif
    CHECK_FLAG ("CVKLU", flag);
    flag = CVSlsSetSparseJacFn(cvode_mem, jacfn);
    CHECK_FLAG("CVSlsSetSparseJacFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_cvode_klu_set_ordering (value vcvode_mem, value vordering)
{
    CAMLparam2(vcvode_mem, vordering);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    int flag = CVKLUSetOrdering (cvode_mem, Int_val(vordering));
    CHECK_FLAG ("CVKLUSetOrdering", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_cvode_klu_reinit (value vcvode_mem, value vn, value vnnz,
				   value vrealloc)
{
    CAMLparam4(vcvode_mem, vn, vnnz, vrealloc);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    int reinit_type;
    if (Bool_val(vrealloc)) {
	reinit_type = 1;
    } else {
	reinit_type = 2;
    }

    int flag = CVKLUReInit (cvode_mem, Int_val(vn), Int_val(vnnz), reinit_type);
    CHECK_FLAG ("CVKLUReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_cvode_klu_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    long int r;
    int flag = CVSlsGetNumJacEvals(cvode_mem, &r);
    CHECK_FLAG("CVSlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

#endif
