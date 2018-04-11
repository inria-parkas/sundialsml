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
CAMLprim value c_arkode_klu_init (value varkode_mem, value vformat,
				  value vneqs, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_klu_set_ordering (value varkode_mem, value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_klu_reinit (value varkode_mem, value vn, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_klu_get_num_jac_evals(value varkode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_mass_klu_init (value varkode_mem, value vformat,
				       value vneqs, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_mass_klu_set_ordering (value varkode_mem,
					       value vordering)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_mass_klu_reinit (value varkode_mem,
					 value vn, value vnnz)
{ CAMLparam0(); CAMLreturn (Val_unit); }

CAMLprim value c_arkode_klu_get_num_mass_evals(value varkode_mem)
{ CAMLparam0(); CAMLreturn (Val_unit); }
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "arkode_ml.h"
#include "../lsolvers/matrix_ml.h"

#include <arkode/arkode.h>

#if SUNDIALS_LIB_VERSION < 300
#include <arkode/arkode_sparse.h>
#include <arkode/arkode_klu.h>
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
    args[0] = arkode_make_jac_arg (t, y, fy,
				  arkode_make_triple_tmp (tmp1, tmp2, tmp3));

    cb = ARKODE_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    smat = Field(cb, 1);
    if (smat == Val_none) {
#if SUNDIALS_LIB_VERSION >= 270
	Store_some(smat, c_matrix_sparse_wrap(Jac, 0, Val_int(Jac->sparsetype)));
#else
	Store_some(smat, c_matrix_sparse_wrap(Jac, 0, Val_int(0)));
#endif
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
#endif

CAMLprim value c_arkode_klu_init (value varkode_mem, value vformat,
				  value vneqs, value vnnz)
{
    CAMLparam4(varkode_mem, vformat, vneqs, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

#if SUNDIALS_LIB_VERSION >= 270
    flag = ARKKLU (arkode_mem, Int_val(vneqs), Int_val(vnnz), Int_val(vformat));
#else
    flag = ARKKLU (arkode_mem, Int_val(vneqs), Int_val(vnnz));
#endif
    CHECK_FLAG ("ARKKLU", flag);
    flag = ARKSlsSetSparseJacFn(arkode_mem, jacfn);
    CHECK_FLAG("ARKSlsSetSparseJacFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_klu_set_ordering (value varkode_mem, value vordering)
{
    CAMLparam2(varkode_mem, vordering);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKKLUSetOrdering (arkode_mem, Int_val(vordering));
    CHECK_FLAG ("ARKKLUSetOrdering", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_klu_reinit (value varkode_mem, value vn, value vnnz)
{
    CAMLparam3(varkode_mem, vn, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int nnz = Int_val(vnnz);

    int flag = ARKKLUReInit (arkode_mem, Int_val(vn), nnz, (nnz > 0) ? 1 : 2);
    CHECK_FLAG ("ARKKLUReInit", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_klu_get_num_jac_evals(value varkode_mem)
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

/* Mass matrix rouintes */

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
    args[1] = arkode_make_triple_tmp (tmp1, tmp2, tmp3);

    cb = ARKODE_MASS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);
    smat = Field(cb, 1);
    if (smat == Val_none) {
#if SUNDIALS_LIB_VERSION >= 270
	Store_some(smat, c_matrix_sparse_wrap(M, 0, Val_int(M->sparsetype)));
#else
	Store_some(smat, c_matrix_sparse_wrap(M, 0, Val_int(0)));
#endif
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
#endif

CAMLprim value c_arkode_mass_klu_init (value varkode_mem, value vformat,
				       value vneqs, value vnnz)
{
    CAMLparam4(varkode_mem, vformat, vneqs, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

#if SUNDIALS_LIB_VERSION >= 270
    flag = ARKMassKLU (arkode_mem, Int_val(vneqs), Int_val(vnnz),
		       MAT_FROM_SFORMAT(vformat), massfn);
#else
    flag = ARKMassKLU (arkode_mem, Int_val(vneqs), Int_val(vnnz), massfn);
#endif
    CHECK_FLAG ("ARKMassKLU", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_mass_klu_set_ordering (value varkode_mem,
					       value vordering)
{
    CAMLparam2(varkode_mem, vordering);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);

    int flag = ARKMassKLUSetOrdering (arkode_mem, Int_val(vordering));
    CHECK_FLAG ("ARKMassKLUSetOrdering", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_mass_klu_reinit (value varkode_mem, value vn,
					 value vnnz)
{
    CAMLparam3(varkode_mem, vn, vnnz);
#if SUNDIALS_LIB_VERSION < 300
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int nnz = Int_val(vnnz);

    int flag = ARKMassKLUReInit (arkode_mem, Int_val(vn), nnz,
				 (nnz > 0) ? 1 : 2);
    CHECK_FLAG ("ARKMassKLUReInit", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_klu_get_num_mass_evals(value varkode_mem)
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
