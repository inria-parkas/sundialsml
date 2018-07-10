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
CAMLprim value sunml_cvode_klu_init (value vcvode_mem, value vformat,
			 	 value vneqs, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_cvode_klu_set_ordering (value vcvode_mem, value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_cvode_klu_reinit (value vcvode_mem, value vn, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_cvode_klu_get_num_jac_evals(value vcvode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "cvode_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#ifdef SUNDIALSML_WITHSENS
/* CVODES (with sensitivity) */

#include <cvodes/cvodes.h>

#if SUNDIALS_LIB_VERSION < 300
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_klu.h>
#endif

#else
/* CVODE (without sensitivity) */

#include <cvode/cvode.h>

#if SUNDIALS_LIB_VERSION < 300
#include <cvode/cvode_sparse.h>
#include <cvode/cvode_klu.h>
#endif

#endif

#if SUNDIALS_LIB_VERSION < 300
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

    // always rewrap without caching (simplified backwards compatibility)
    args[1] = c_matrix_sparse_wrap(Jac);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

CAMLprim value sunml_cvode_klu_init (value vcvode_mem, value vformat,
				 value vneqs, value vnnz)
{
    CAMLparam4(vcvode_mem, vformat, vneqs, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

#if SUNDIALS_LIB_VERSION >= 270
    flag = CVKLU (cvode_mem, Int_val(vneqs), Int_val(vnnz),
		  MAT_FROM_SFORMAT(vformat));
#else
    flag = CVKLU (cvode_mem, Int_val(vneqs), Int_val(vnnz));
#endif
    CHECK_FLAG ("CVKLU", flag);
    flag = CVSlsSetSparseJacFn(cvode_mem, jacfn);
    CHECK_FLAG("CVSlsSetSparseJacFn", flag);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_klu_set_ordering (value vcvode_mem, value vordering)
{
    CAMLparam2(vcvode_mem, vordering);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    int flag = CVKLUSetOrdering (cvode_mem, Int_val(vordering));
    CHECK_FLAG ("CVKLUSetOrdering", flag);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_klu_reinit (value vcvode_mem, value vn, value vnnz)
{
    CAMLparam3(vcvode_mem, vn, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int nnz = Int_val(vnnz);

    int flag = CVKLUReInit (cvode_mem, Int_val(vn), nnz, (nnz > 0) ? 1 : 2);
    CHECK_FLAG ("CVKLUReInit", flag);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvode_klu_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    long int r = 0;
#if SUNDIALS_LIB_VERSION < 300
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);

    int flag = CVSlsGetNumJacEvals(cvode_mem, &r);
    CHECK_FLAG("CVSlsGetNumJacEvals", flag);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

#endif
