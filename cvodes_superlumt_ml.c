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

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_superlumt.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>

#include "sundials_ml.h"
#include "cvode_ml.h"
#include "cvodes_ml.h"
#include "sls_ml.h"

static int bjacfn_nosens( /* CVSlsSparseJacFnB */
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

    args[0] = c_cvodes_make_jac_arg (t, y, yb, fyb,
			    c_cvode_make_triple_tmp (tmp1b, tmp2b, tmp3b));

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
    value r = caml_callback2_exn (Field(cb, 0), args[0], args[1]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bjacfn_withsens( /* CVSlsSparseJacFnBS */
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

    args[0] = c_cvodes_make_jac_arg (t, y, yb, fyb,
			    c_cvode_make_triple_tmp (tmp1b, tmp2b, tmp3b));

    ns = Int_val(Field(bsensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = CVODES_BSENSARRAY_FROM_EXT(bsensext);
    c_cvodes_wrap_to_nvector_table(ns, args[1], ys);

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
    value r = caml_callback3_exn (Field(cb, 0), args[0], args[1], args[2]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

CAMLprim value c_cvode_superlumtb_init (value vparent_which,
					value vneqs, value vnnz,
					value vnthreads, value vusesens)
{
    CAMLparam5(vparent_which, vneqs, vnnz, vnthreads, vusesens);
    void *mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = CVSuperLUMTB (mem, which, Int_val(vnthreads), Int_val(vneqs),
			 Int_val(vnnz));
    CHECK_FLAG ("CVSuperLUMTB", flag);
    if (Bool_val(vusesens)) {
	flag = CVSlsSetSparseJacFnBS(mem, which, bjacfn_withsens);
	CHECK_FLAG("CVSlsSetSparseJacFnBS", flag);
    } else {
	flag = CVSlsSetSparseJacFnB (mem, which, bjacfn_nosens);
	CHECK_FLAG("CVSlsSetSparseJacFnB", flag);
    }


    CAMLreturn (Val_unit);
}

