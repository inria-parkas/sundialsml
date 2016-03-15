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
CAMLprim value c_arkode_superlumt_init (value varkode_mem, value vneqs,
					value vnnz, value vnthreads)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_superlumt_set_ordering (value varkode_mem,
						value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_superlumt_get_num_jac_evals(value varkode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_mass_superlumt_init (value varkode_mem, value vneqs,
					     value vnnz, value vnthreads)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_mass_superlumt_set_ordering (value varkode_mem,
					             value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_superlumt_get_num_mass_evals(value varkode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "arkode_ml.h"
#include "../lsolvers/sls_ml.h"

#include <arkode/arkode.h>
#include <arkode/arkode_sparse.h>
#include <arkode/arkode_superlumt.h>

enum arkode_superlumt_ordering_tag {
  VARIANT_ARKODE_SUPERLUMT_NATURAL    = 0,
  VARIANT_ARKODE_SUPERLUMT_MINDEGPROD = 1,
  VARIANT_ARKODE_SUPERLUMT_MINDEGSUM  = 2,
  VARIANT_ARKODE_SUPERLUMT_COLAMD     = 3,
};

static int jacfn(realtype t,
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
    args[0] = arkode_make_jac_arg (t, y, fy,
				  arkode_make_triple_tmp (tmp1, tmp2, tmp3));

    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
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

CAMLprim value c_arkode_superlumt_init (value varkode_mem, value vneqs,
					value vnnz, value vnthreads)
{
    CAMLparam4(varkode_mem, vneqs, vnnz, vnthreads);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKSuperLUMT (arkode_mem,
			 Int_val(vnthreads),
			 Int_val(vneqs),
			 Int_val(vnnz));
    CHECK_FLAG ("ARKSuperLUMT", flag);
    flag = ARKSlsSetSparseJacFn(arkode_mem, jacfn);
    CHECK_FLAG("ARKSlsSetSparseJacFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_superlumt_set_ordering (value varkode_mem,
						value vordering)
{
    CAMLparam2(varkode_mem, vordering);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKSuperLUMTSetOrdering (arkode_mem, Int_val(vordering));
    CHECK_FLAG ("ARKSuperLUMTSetOrdering", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_superlumt_get_num_jac_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    long int r;
    int flag = ARKSlsGetNumJacEvals(arkode_mem, &r);
    CHECK_FLAG("ARKSlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

/* Mass matrix rouintes */

static int massfn(realtype t,
		  SlsMat M,
		  void *user_data,
		  N_Vector tmp1,
		  N_Vector tmp2,
		  N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 3);
    CAMLlocal3(session, cb, smat);

    WEAK_DEREF (session, *(value*)user_data);
    args[0] = caml_copy_double(t);
    args[1] = arkode_make_triple_tmp (tmp1, tmp2, tmp3);

    cb = ARKODE_MASS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, c_sls_sparse_wrap(M, 0));
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

CAMLprim value c_arkode_mass_superlumt_init (value varkode_mem, value vneqs,
					     value vnnz, value vnthreads)
{
    CAMLparam4(varkode_mem, vneqs, vnnz, vnthreads);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSuperLUMT (arkode_mem,
			     Int_val(vnthreads),
			     Int_val(vneqs),
			     Int_val(vnnz),
			     massfn);
    CHECK_FLAG ("ARKMassSuperLUMT", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_mass_superlumt_set_ordering (value varkode_mem,
					             value vordering)
{
    CAMLparam2(varkode_mem, vordering);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKMassSuperLUMTSetOrdering (arkode_mem, Int_val(vordering));
    CHECK_FLAG ("ARKMassSuperLUMTSetOrdering", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_superlumt_get_num_mass_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    long int r;
    int flag = ARKSlsGetNumMassEvals(arkode_mem, &r);
    CHECK_FLAG("ARKSlsGetNumMassEvals", flag);

    CAMLreturn(Val_long(r));
}

#endif
