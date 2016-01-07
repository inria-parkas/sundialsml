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

#include <idas/idas.h>
#include <idas/idas_sparse.h>
#include <idas/idas_klu.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>

#include "../sundials/sundials_ml.h"
#include "../ida/ida_ml.h"
#include "idas_ml.h"
#include "../lsolvers/sls_ml.h"

static int jacfn_nosens( /* IDASlsSparseJacFnB */
	realtype t,
	realtype cjB,
	N_Vector yy,
	N_Vector yp,
	N_Vector yyB,
	N_Vector ypB,
	N_Vector resvalB,
	SlsMat jacB,
	void *user_data,
	N_Vector tmp1B,
	N_Vector tmp2B,
	N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, smat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));

    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, c_sls_sparse_wrap(jacB, 0));
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

static int jacfn_withsens( /* IDASlsSparseJacFnB */
	realtype t,
	realtype cjB,
	N_Vector yy,
	N_Vector yp,
	N_Vector *ys,
	N_Vector *yps,
	N_Vector yyB,
	N_Vector ypB,
	N_Vector resvalB,
	SlsMat jacB,
	void *user_data,
	N_Vector tmp1B,
	N_Vector tmp2B,
	N_Vector tmp3B)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal4(session, bsensext, cb, smat);

    WEAK_DEREF (session, *(value*)user_data);
    bsensext = IDA_SENSEXT_FROM_ML(session);

    cb = IDA_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = idas_make_jac_arg(t, yy, yp, yyB, ypB, resvalB, cjB,
			        ida_make_triple_tmp (tmp1B, tmp2B, tmp3B));

    int ns = Int_val(Field(bsensext, RECORD_IDAS_BWD_SESSION_NUMSENSITIVITIES));
    args[1] = IDAS_BSENSARRAY1_FROM_EXT(bsensext);
    idas_wrap_to_nvector_table(ns, args[1], ys);
    args[2] = IDAS_BSENSARRAY2_FROM_EXT(bsensext);
    idas_wrap_to_nvector_table(ns, args[2], yps);

    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, c_sls_sparse_wrap(jacB, 0));
	Store_field(cb, 1, smat);

	args[3] = Some_val(smat);
    } else {
	args[3] = Some_val(smat);
	c_sparsematrix_realloc(args[3], 0);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

CAMLprim value c_idas_klub_init (value vparent, value vwhich,
				 value vneqs, value vnnz, value vusesens)
{
    CAMLparam5(vparent, vwhich, vneqs, vnnz, vusesens);
    void *mem = IDA_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = IDAKLUB (mem, which, Int_val(vneqs), Int_val(vnnz));
    CHECK_FLAG ("IDAKLUB", flag);
    if (Bool_val(vusesens)) {
	flag = IDASlsSetSparseJacFnBS(mem, which, jacfn_withsens);
	CHECK_FLAG("IDASlsSetSparseJacFnBS", flag);
    } else {
	flag = IDASlsSetSparseJacFnB (mem, which, jacfn_nosens);
	CHECK_FLAG("IDASlsSetSparseJacFnB", flag);
    }

    CAMLreturn (Val_unit);
}

