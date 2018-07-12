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
CAMLprim value sunml_arkode_superlumt_init (value varkode_mem, value vneqs,
					value vnnz, value vnthreads)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_arkode_superlumt_set_ordering (value varkode_mem,
						value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_arkode_superlumt_get_num_jac_evals(value varkode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_arkode_mass_superlumt_init (value varkode_mem, value vneqs,
					     value vnnz, value vnthreads)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_arkode_mass_superlumt_set_ordering (value varkode_mem,
					             value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value sunml_arkode_superlumt_get_num_mass_evals(value varkode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "arkode_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#include <arkode/arkode.h>

#if SUNDIALS_LIB_VERSION < 300
#include <arkode/arkode_sparse.h>
#include <arkode/arkode_superlumt.h>
#endif

#if SUNDIALS_LIB_VERSION < 300
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
    args[0] = sunml_arkode_make_jac_arg (t, y, fy,
				  sunml_arkode_make_triple_tmp (tmp1, tmp2, tmp3));

    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, sunml_matrix_sparse_wrap(Jac));
	Store_field(cb, 1, smat);
	args[1] = Some_val(smat);
    } else {
	args[1] = Some_val(smat);
	sunml_matrix_sparse_rewrap(args[1]);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (Field(cb, 0), args[0], args[1]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

CAMLprim value sunml_arkode_superlumt_init (value varkode_mem, value vneqs,
					value vnnz, value vnthreads)
{
    CAMLparam4(varkode_mem, vneqs, vnnz, vnthreads);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKSuperLUMT (arkode_mem,
			 Int_val(vnthreads),
			 Int_val(vneqs),
			 Int_val(vnnz));
    CHECK_FLAG ("ARKSuperLUMT", flag);
    flag = ARKSlsSetSparseJacFn(arkode_mem, jacfn);
    CHECK_FLAG("ARKSlsSetSparseJacFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_superlumt_set_ordering (value varkode_mem,
						value vordering)
{
    CAMLparam2(varkode_mem, vordering);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKSuperLUMTSetOrdering (arkode_mem, Int_val(vordering));
    CHECK_FLAG ("ARKSuperLUMTSetOrdering", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_superlumt_get_num_jac_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r = 0;
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKSlsGetNumJacEvals(arkode_mem, &r);
    CHECK_FLAG("ARKSlsGetNumJacEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

/* Mass matrix routines */

#if SUNDIALS_LIB_VERSION < 300
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
    args[1] = sunml_arkode_make_triple_tmp (tmp1, tmp2, tmp3);

    cb = ARKODE_MASS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    smat = Field(cb, 1);
    if (smat == Val_none) {
	Store_some(smat, sunml_matrix_sparse_wrap(M));
	Store_field(cb, 1, smat);
	args[2] = Some_val(smat);
    } else {
	args[2] = Some_val(smat);
	sunml_matrix_sparse_rewrap(args[2]);
    }

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (Field(cb, 0), args[0], args[1], args[2]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}
#endif

CAMLprim value sunml_arkode_mass_superlumt_init (value varkode_mem, value vneqs,
					     value vnnz, value vnthreads)
{
    CAMLparam4(varkode_mem, vneqs, vnnz, vnthreads);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKMassSuperLUMT (arkode_mem,
			     Int_val(vnthreads),
			     Int_val(vneqs),
			     Int_val(vnnz),
			     massfn);
    CHECK_FLAG ("ARKMassSuperLUMT", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_mass_superlumt_set_ordering (value varkode_mem,
					             value vordering)
{
    CAMLparam2(varkode_mem, vordering);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKMassSuperLUMTSetOrdering (arkode_mem, Int_val(vordering));
    CHECK_FLAG ("ARKMassSuperLUMTSetOrdering", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_arkode_superlumt_get_num_mass_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    long int r = 0;
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKSlsGetNumMassEvals(arkode_mem, &r);
    CHECK_FLAG("ARKSlsGetNumMassEvals", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_long(r));
}

#endif
