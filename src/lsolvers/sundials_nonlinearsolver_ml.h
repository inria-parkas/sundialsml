/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2020 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#ifndef _NLSOLVER_ML_H__
#define _NLSOLVER_ML_H__

#include <caml/mlvalues.h>

#if SUNDIALS_LIB_VERSION >= 400
#include <sundials/sundials_nonlinearsolver.h>

// turn a cptr into SUNNonlinearSolver
#define NLSOLVER_VAL(v) (*(SUNNonlinearSolver *)Data_custom_val(v))

/*
  Extend the nonlinear solver (nls) struct with extra fields required by
  OCaml.

  callbacks - provides access to OCaml implementations of the four
	      callbacks: sysfn, lsetupfn, lsolvefn, convtestfn.

  When integrators invoke setup or solve on an nls, they pass the active
  session as the "void *mem" argument. To permit callbacks into OCaml, we
  override the setup and solve functions to transform the mem argument
  into an nls_callback_value that is passed to the original setup and
  solve functions.

  C callbacks access the first field of an nls_callback_value and ignore the
  other two fields. OCaml callbacks ignore the first field and use the two
  other fields.

  The to_value function is set by the calling integrator and provides a
  means of transforming the C session to an OCaml session value.

  The from_value function is also set by the calling integrator to transform
  and OCaml session value into a C session. This permits the OCaml
  implementations of the setup and solve functions to invoke the
  integrator's callbacks with the required value.

  The to_value function is only used on callbacks from integrators and thus
  may be set to NULL for an nls used directly from OCaml. The from_value
  function is only required for integrator callbacks; if one is not
  specified, the C mem arguments is simply set to NULL since it is unused.
 */

struct sunml_nls_callbacks {
    SUNNonlinSolSysFn sysfn;
    SUNNonlinSolLSetupFn lsetupfn;
    SUNNonlinSolLSolveFn lsolvefn;
    SUNNonlinSolConvTestFn ctestfn;
};
typedef struct sunml_nls_callbacks* p_sunml_nls_callbacks;

struct sunml_nls {
    struct _generic_SUNNonlinearSolver nls;
    value callbacks;

    int (*orig_setup)(SUNNonlinearSolver, N_Vector, void*);
    int (*orig_solve)(SUNNonlinearSolver, N_Vector, N_Vector, N_Vector, realtype,
		      booleantype, void*);
    value (*to_value)(void *);
    void* (*from_value)(value);

    int (*orig_setsysfn)(SUNNonlinearSolver, SUNNonlinSolSysFn);
    int (*orig_setlsetupfn)(SUNNonlinearSolver, SUNNonlinSolLSetupFn);
    int (*orig_setlsolvefn)(SUNNonlinearSolver, SUNNonlinSolLSolveFn);
#if 500 <= SUNDIALS_LIB_VERSION
    int (*orig_setctestfn)(SUNNonlinearSolver, SUNNonlinSolConvTestFn, void*);
#else
    int (*orig_setctestfn)(SUNNonlinearSolver, SUNNonlinSolConvTestFn);
#endif

    struct sunml_nls_callbacks c_callbacks;
};
typedef struct sunml_nls* p_sunml_nls;

void sunml_nlsolver_set_to_from_mem(SUNNonlinearSolver nls,
				    value (*to_value)(void *),
				    void* (*from_value)(value));

#define NLS_EXTENDED(nls) ((struct sunml_nls *)nls)
#define NLS_CALLBACKS(nls) ((NLS_EXTENDED(nls))->callbacks)

void sunml_nlsolver_check_flag(const char *call, int flag);
#define NLS_CHECK_FLAG(call, flag) if (flag != SUN_NLS_SUCCESS) \
				 sunml_nlsolver_check_flag(call, flag)

#endif

enum nlsolver_callbacks_index {
  RECORD_NLSOLVER_CALLBACKS_SYSFN	   = 0,
  RECORD_NLSOLVER_CALLBACKS_LSETUPFN,
  RECORD_NLSOLVER_CALLBACKS_LSOLVEFN,
  RECORD_NLSOLVER_CALLBACKS_CONVTESTFN,
  RECORD_NLSOLVER_CALLBACKS_SIZE
};

enum nlsolver_ops_index {
    RECORD_NLSOLVER_OPS_TYPE = 0,
    RECORD_NLSOLVER_OPS_INIT,
    RECORD_NLSOLVER_OPS_SETUP,
    RECORD_NLSOLVER_OPS_SOLVE,
    RECORD_NLSOLVER_OPS_SET_SYS_FN,
    RECORD_NLSOLVER_OPS_SET_LSETUP_FN,
    RECORD_NLSOLVER_OPS_SET_LSOLVE_FN,
    RECORD_NLSOLVER_OPS_SET_CONVTEST_FN,
    RECORD_NLSOLVER_OPS_SET_MAX_ITERS,
    RECORD_NLSOLVER_OPS_SET_INFO_FILE,
    RECORD_NLSOLVER_OPS_SET_PRINT_LEVEL,
    RECORD_NLSOLVER_OPS_GET_NUM_ITERS,
    RECORD_NLSOLVER_OPS_GET_CUR_ITER,
    RECORD_NLSOLVER_OPS_GET_NUM_CONV_FAILS,
};

enum nlsolver_index {
  RECORD_NLSOLVER_RAWPTR	   = 0,
  RECORD_NLSOLVER_SOLVER,
  RECORD_NLSOLVER_ATTACHED,
  RECORD_NLSOLVER_SIZE
};

enum nlsolver_convtest_tag {
  VARIANT_NLSOLVER_CONVTEST_SUCCESS        = 0,
  VARIANT_NLSOLVER_CONVTEST_CONTINUE,
  VARIANT_NLSOLVER_CONVTEST_RECOVER,
};

enum nlsolver_type {
  VARIANT_NLSOLVER_TYPE_ROOT_FIND        = 0,
  VARIANT_NLSOLVER_TYPE_FIXED_POINT
};

/* This enum must list exceptions in the same order as the call to
 * c_init_module in sundials_NonlinearSolver.ml.  */
enum nlsolver_exn_index {
    NLSOLVER_EXN_VectorOpError,
    NLSOLVER_EXN_IncorrectUse,
    NLSOLVER_EXN_ExtFail,
    NLSOLVER_EXN_SET_SIZE
};

#define NLSOLVER_EXN(name)     REGISTERED_EXN(NLSOLVER, name)
#define NLSOLVER_EXN_TAG(name) REGISTERED_EXN_TAG(NLSOLVER, name)

#endif
