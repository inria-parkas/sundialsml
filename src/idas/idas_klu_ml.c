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
CAMLprim value sunml_idas_klub_init (value vparent_which, value vformat,
				 value vneqs, value vnnz, value vusesens)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "../ida/ida_ml.h"
#include "idas_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#include <idas/idas.h>

#if SUNDIALS_LIB_VERSION < 300
#include <idas/idas_sparse.h>
#include <idas/idas_klu.h>
#endif

#if SUNDIALS_LIB_VERSION < 300
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
	Store_some(smat, c_matrix_sparse_wrap(jacB));
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
	Store_some(smat, c_matrix_sparse_wrap(jacB));
	Store_field(cb, 1, smat);

	args[3] = Some_val(smat);
    } else {
	args[3] = Some_val(smat);
	sunml_matrix_sparse_rewrap(args[3]);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

CAMLprim value sunml_idas_klub_init (value vparent_which, value vformat,
				 value vneqs, value vnnz, value vusesens)
{
    CAMLparam5(vparent_which, vformat, vneqs, vnnz, vusesens);
#if SUNDIALS_LIB_VERSION < 300
    void *mem = IDA_MEM_FROM_ML (Field(vparent_which, 0));
    int which = Int_val(Field(vparent_which, 1));
    int flag;

#if SUNDIALS_LIB_VERSION >= 270
    flag = IDAKLUB (mem, which, Int_val(vneqs), Int_val(vnnz),
		    Int_val(vformat));
#else
    flag = IDAKLUB (mem, which, Int_val(vneqs), Int_val(vnnz));
#endif
    CHECK_FLAG ("IDAKLUB", flag);
    if (Bool_val(vusesens)) {
	flag = IDASlsSetSparseJacFnBS(mem, which, jacfn_withsens);
	CHECK_FLAG("IDASlsSetSparseJacFnBS", flag);
    } else {
	flag = IDASlsSetSparseJacFnB (mem, which, jacfn_nosens);
	CHECK_FLAG("IDASlsSetSparseJacFnB", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

#endif
