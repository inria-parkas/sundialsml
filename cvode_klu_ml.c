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

#include <cvode/cvode.h>
#include <cvode/cvode_sparse.h>
#include <cvode/cvode_klu.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>

#include "sundials_ml.h"
#include "cvode_ml.h"
#include "sls_ml.h"

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
    args[0] = c_cvode_make_jac_arg (t, y, fy,
				    c_cvode_make_triple_tmp (tmp1, tmp2, tmp3));

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
    value r = caml_callback2_exn (Field(cb, 0), args[0], args[1]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

CAMLprim value c_cvode_klu_init (value vcvode_mem, value vneqs, value vnnz)
{
    CAMLparam3(vcvode_mem, vneqs, vnnz);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVKLU (cvode_mem, Int_val(vneqs), Int_val(vnnz));
    CHECK_FLAG ("CVKLU", flag);
    flag = CVSlsSetSparseJacFn(CVODE_MEM_FROM_ML(vcvode_mem), jacfn);
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

