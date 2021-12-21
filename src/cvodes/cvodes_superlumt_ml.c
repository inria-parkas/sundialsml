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
#include "../nvectors/nvector_ml.h"

#ifndef SUNDIALS_ML_SUPERLUMT
CAMLprim value sunml_cvodes_superlumtb_init (value vparent_which,
		 			 value vneqs, value vnnz,
					 value vnthreads, value vusesens)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "../cvode/cvode_ml.h"
#include "cvodes_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#include <cvodes/cvodes.h>

#if SUNDIALS_LIB_VERSION < 300
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_superlumt.h>
#endif

#if SUNDIALS_LIB_VERSION < 300
static int jacfn_nosens( /* CVSlsSparseJacFnB */
    sunrealtype t,
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

    args[0] = sunml_cvodes_make_jac_arg (t, y, yb, fyb,
			    sunml_cvode_make_triple_tmp (tmp1b, tmp2b, tmp3b));

    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, sunml_matrix_sparse_wrap(jacb));
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

static int jacfn_withsens( /* CVSlsSparseJacFnBS */
    sunrealtype t,
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

    args[0] = sunml_cvodes_make_jac_arg (t, y, yb, fyb,
			    sunml_cvode_make_triple_tmp (tmp1b, tmp2b, tmp3b));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    sunml_nvectors_into_array(ns, args[1], ys);

    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, sunml_matrix_sparse_wrap(jacb));
	Store_field(cb, 1, smat);

	args[2] = Some_val(smat);
    } else {
	args[2] = Some_val(smat);
	sunml_matrix_sparse_rewrap(args[2]);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

CAMLprim value sunml_cvodes_superlumtb_init (value vparent_which,
					 value vneqs, value vnnz,
					 value vnthreads, value vusesens)
{
    CAMLparam5(vparent_which, vneqs, vnnz, vnthreads, vusesens);
#if SUNDIALS_LIB_VERSION < 300
    void *mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = CVSuperLUMTB (mem, which, Int_val(vnthreads), Int_val(vneqs),
			 Int_val(vnnz));
    CHECK_FLAG ("CVSuperLUMTB", flag);
    if (Bool_val(vusesens)) {
	flag = CVSlsSetSparseJacFnBS(mem, which, jacfn_withsens);
	CHECK_FLAG("CVSlsSetSparseJacFnBS", flag);
    } else {
	flag = CVSlsSetSparseJacFnB (mem, which, jacfn_nosens);
	CHECK_FLAG("CVSlsSetSparseJacFnB", flag);
    }

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

#endif
