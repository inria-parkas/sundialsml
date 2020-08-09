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

#define NLMEM_VAL(v) (*(void **)Data_custom_val(v))

struct sunml_nls {
    struct _generic_SUNNonlinearSolver nls;
    value callbacks;
};

#define NLS_CALLBACKS(nls) (((struct sunml_nls *)nls)->callbacks)

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

enum nlsolver_solver_tag {
  VARIANT_NLSOLVER_SOLVER_NEWTONSOLVER            = 0, // constant constrs
  VARIANT_NLSOLVER_SOLVER_FIXEDPOINTSOLVER        = 0, // non-constant constrs
  VARIANT_NLSOLVER_SOLVER_CUSTOMSOLVER            = 1,
  VARIANT_NLSOLVER_SOLVER_CUSTOMSENSSOLVER        = 2
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
    RECORD_NLSOLVER_OPS_GET_NUM_ITERS,
    RECORD_NLSOLVER_OPS_GET_CUR_ITER,
    RECORD_NLSOLVER_OPS_GET_NUM_CONV_FAILS,
};

enum nlsolver_index {
  RECORD_NLSOLVER_RAWPTR	   = 0,
  RECORD_NLSOLVER_SOLVER,
  RECORD_NLSOLVER_CALLBACKS,
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
    NLSOLVER_EXN_SET_SIZE
};

#define NLSOLVER_EXN(name)     REGISTERED_EXN(NLSOLVER, name)
#define NLSOLVER_EXN_TAG(name) REGISTERED_EXN_TAG(NLSOLVER, name)

#endif
