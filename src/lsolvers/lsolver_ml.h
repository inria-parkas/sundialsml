/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2018 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#ifndef _LSOLVER_ML_H__
#define _LSOLVER_ML_H__

#include <caml/mlvalues.h>

#if SUNDIALS_LIB_VERSION >= 300

// turn a cptr into SUNLinearSolver
#define LSOLVER_VAL(v) (*(SUNLinearSolver *)Data_custom_val(v))

#endif

enum preconditioning_type_tag {
  VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_NONE  = 0,
  VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_LEFT,
  VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_RIGHT,
  VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_BOTH,
};

int lsolver_precond_type(value);
int lsolver_gs_type(value);

enum lsolver_iterative_solver_tag {
    VARIANT_LSOLVER_ITERATIVE_SOLVER_SPBCGS = 0,
    VARIANT_LSOLVER_ITERATIVE_SOLVER_SPFGMR,
    VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR,
    VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR,
    VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG,
};

// values must match convention for SUNKLUSetOrdering
enum lsolver_klu_ordering_tag {
    VARIANT_LSOLVER_KLU_ORDERING_AMD     = 0,
    VARIANT_LSOLVER_KLU_ORDERING_COLAMD  = 1,
    VARIANT_LSOLVER_KLU_ORDERING_NATURAL = 2,
};

// values must match convention for SUNSuperLUMTSetOrdering
enum lsolver_superlumt_ordering_tag {
    VARIANT_LSOLVER_SUPERLUMT_ORDERING_NATURAL       = 0,
    VARIANT_LSOLVER_SUPERLUMT_ORDERING_MINDEGREEPROD = 1,
    VARIANT_LSOLVER_SUPERLUMT_ORDERING_MINDEGREESUM  = 2,
    VARIANT_LSOLVER_SUPERLUMT_ORDERING_COLAMD	     = 3,
};

enum lsolver_gramschmidt_type_tag {
    VARIANT_LSOLVER_GRAMSCHMIDT_TYPE_MODIFIEDGS = 0,
    VARIANT_LSOLVER_GRAMSCHMIDT_TYPE_CLASSICALGS,
};

// values must match Lsolver_impl.Custom.ops type
enum lsolver_ops_index {
    RECORD_LSOLVER_OPS_INIT = 0,
    RECORD_LSOLVER_OPS_SETUP,
    RECORD_LSOLVER_OPS_SOLVE,
    RECORD_LSOLVER_OPS_SET_ATIMES,
    RECORD_LSOLVER_OPS_SET_PRECONDITIONER,
    RECORD_LSOLVER_OPS_SET_SCALING_VECTORS,
    RECORD_LSOLVER_OPS_GET_NUM_ITERS,
    RECORD_LSOLVER_OPS_GET_RES_NORM,
    RECORD_LSOLVER_OPS_GET_RES_ID,
    RECORD_LSOLVER_OPS_GET_WORK_SPACE,
};

// values must match Lsolver_impl.Custom.has_ops type
enum lsolver_hasops_index {
    RECORD_LSOLVER_HASOPS_SET_ATIMES = 0,
    RECORD_LSOLVER_HASOPS_SET_PRECONDITIONER,
    RECORD_LSOLVER_HASOPS_SET_SCALING_VECTORS,
    RECORD_LSOLVER_HASOPS_GET_NUM_ITERS,
    RECORD_LSOLVER_HASOPS_GET_RES_NORM,
    RECORD_LSOLVER_HASOPS_GET_RES_ID,
    RECORD_LSOLVER_HASOPS_GET_WORK_SPACE,
};

enum lsolver_preconditioning_type_tag {
  VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_NONE  = 0,
  VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_LEFT,
  VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_RIGHT,
  VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_BOTH,
};

/* This enum must list exceptions in the same order as the call to
 * c_init_module in lsolver.ml.  */
enum lsolver_exn_index {
    LSOLVER_EXN_UnrecoverableFailure = 0,
    LSOLVER_EXN_MatrixNotSquare,
    LSOLVER_EXN_InvalidArgument,
    LSOLVER_EXN_ATimesFailure,
    LSOLVER_EXN_PSetFailure,
    LSOLVER_EXN_PSolveFailure,
    LSOLVER_EXN_GSFailure,
    LSOLVER_EXN_QRSolFailure,
    LSOLVER_EXN_ResReduced,
    LSOLVER_EXN_ConvFailure,
    LSOLVER_EXN_QRfactFailure,
    LSOLVER_EXN_LUfactFailure,
    LSOLVER_EXN_PackageFailure,
    LSOLVER_EXN_IllegalPrecType,
    LSOLVER_EXN_InternalFailure,
    LSOLVER_EXN_SET_SIZE
};


#define LSOLVER_EXN(name)     REGISTERED_EXN(LSOLVER, name)
#define LSOLVER_EXN_TAG(name) REGISTERED_EXN_TAG(LSOLVER, name)

#endif

