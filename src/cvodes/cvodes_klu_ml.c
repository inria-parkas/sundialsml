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
#include "../cvode/cvode_ml.h"
#include "cvodes_ml.h"
#include "../lsolvers/sls_ml.h"

#ifndef SUNDIALS_ML_KLU
CAMLprim value c_cvodes_klub_init (value vparent, value vwhich,
				   value vneqs, value vnnz, value vusesens)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_klu.h>

static int jacfn_nosens( /* CVSlsSparseJacFnB */
    realtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    SlsMat jacb,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, smat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = cvodes_make_jac_arg (t, y, yb, fyb,
			    cvode_make_triple_tmp (tmp1b, tmp2b, tmp3b));

    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, c_sls_sparse_wrap(jacb, 0));
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

static int jacfn_withsens( /* CVSlsSparseJacFnBS */
    realtype t,
    N_Vector y,
    N_Vector *ys,
    N_Vector yb,
    N_Vector fyb,
    SlsMat jacb,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal4(session, bsensext, cb, smat);
    int ns;

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = CVODE_SENSEXT_FROM_ML(session);

    cb = CVODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = cvodes_make_jac_arg (t, y, yb, fyb,
			    cvode_make_triple_tmp (tmp1b, tmp2b, tmp3b));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    cvodes_wrap_to_nvector_table(ns, args[1], ys);

    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, c_sls_sparse_wrap(jacb, 0));
	Store_field(cb, 1, smat);

	args[2] = Some_val(smat);
    } else {
	args[2] = Some_val(smat);
	c_sparsematrix_realloc(args[2], 0);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}


#if SUNDIALS_LIB_VERSION < 262
// Work around a bug in Sundials 2.6.0 and 2.6.1
SUNDIALS_EXPORT int CVKLUB(void *, int, int, int);
#endif

CAMLprim value c_cvodes_klub_init (value vparent, value vwhich,
				   value vneqs, value vnnz, value vusesens)
{
    CAMLparam5(vparent, vwhich, vneqs, vnnz, vusesens);
    void *mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVKLUB (mem, which, Int_val(vneqs), Int_val(vnnz));
    CHECK_FLAG ("CVKLUB", flag);
    if (Bool_val(vusesens)) {
	flag = CVSlsSetSparseJacFnBS(mem, which, jacfn_withsens);
	CHECK_FLAG("CVSlsSetSparseJacFnBS", flag);
    } else {
	flag = CVSlsSetSparseJacFnB (mem, which, jacfn_nosens);
	CHECK_FLAG("CVSlsSetSparseJacFnB", flag);
    }

    CAMLreturn (Val_unit);
}

#endif
