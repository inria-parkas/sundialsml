/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

#ifndef _IDA_ML_H__
#define _IDA_ML_H__

#include "sundials_ml.h"

void ida_ml_check_flag(const char *call, int flag);
void ida_ml_set_linear_solver(void *ida_mem, value ls, int n);

#define CHECK_FLAG(call, flag) if (flag != IDA_SUCCESS) \
				 ida_ml_check_flag(call, flag)


/* Indices into the Ida_*.session type.  This enum must be in the same order as
 * the session type's member declaration.  The session data structure is shared
 * in four parts across the OCaml and C heaps.  See cvode_ml.h for details.
 */
enum ida_index {
    RECORD_IDA_SESSION_MEM = 0,
    RECORD_IDA_SESSION_BACKREF,
    RECORD_IDA_SESSION_NROOTS,
    RECORD_IDA_SESSION_ERRFILE,
    RECORD_IDA_SESSION_EXN_TEMP,
    RECORD_IDA_SESSION_RESFN,
    RECORD_IDA_SESSION_ROOTSFN,
    RECORD_IDA_SESSION_ERRH,
    RECORD_IDA_SESSION_ERRW,
    RECORD_IDA_SESSION_LS_CALLBACKS,
    RECORD_IDA_SESSION_SAFETY_FLAGS,
    RECORD_IDA_SESSION_SIZE	/* This has to come last. */
};

#define IDA_MEM_FROM_ML(v) ((void *)Field((v), RECORD_IDA_SESSION_MEM))
#define IDA_BACKREF_FROM_ML(v) ((value*)Field((v), RECORD_IDA_SESSION_BACKREF))
#define IDA_NEQS_FROM_ML(v)    Long_val(Field((v), RECORD_IDA_SESSION_NEQS))
#define IDA_NROOTS_FROM_ML(v)  Long_val(Field((v), RECORD_IDA_SESSION_NROOTS))
#define IDA_ERRFILE_FROM_ML(v) Long_val(Field((v), RECORD_IDA_SESSION_ERRFILE))
#define IDA_RESFN_FROM_ML(v)   (Field((v), RECORD_IDA_SESSION_RESFN))
#define IDA_ROOTSFN_FROM_ML(v) (Field((v), RECORD_IDA_SESSION_ROOTSFN))
#define IDA_ERRH_FROM_ML(v)    (Field((v), RECORD_IDA_SESSION_ERRH))
#define IDA_ERRW_FROM_ML(v)    (Field((v), RECORD_IDA_SESSION_ERRW))
#define IDA_LS_CALLBACKS_FROM_ML(v)     (Field((v), RECORD_IDA_SESSION_LS_CALLBACKS))

#if SAFETY_CHECKS

/* Get the values of all flags in one go.  */
#define IDA_SAFETY_FLAGS(v)  \
    Int_val(Field((v), RECORD_IDA_SESSION_SAFETY_FLAGS))

/* Get the value of a particular flag or flags.  */
#define IDA_TEST_SAFETY_FLAG(v, f) (IDA_SAFETY_FLAGS(v) & (f))

/* Specify the values of all flags in one go.  */
#define IDA_WRITE_SAFETY_FLAGS(v, f) \
    Store_field ((v), RECORD_IDA_SESSION_SAFETY_FLAGS, Val_int ((f)))

/* This is like &=.  Set all flags to 0 except those set in m.  */
#define IDA_MASK_SAFETY_FLAGS(v, m) \
    IDA_WRITE_SAFETY_FLAGS(v, IDA_SAFETY_FLAGS (v) & (m))

/* This is like |=.  Set a particular flag or flags to 1.  */
#define IDA_SET_SAFETY_FLAG(v, m)				\
    IDA_WRITE_SAFETY_FLAGS(v, IDA_SAFETY_FLAGS(v) | (m))

/* This is like &= with ~.  Set a particular flag or flags to 0.  */
#define IDA_UNSET_SAFETY_FLAG(v, f)				\
    IDA_WRITE_SAFETY_FLAGS(v, IDA_SAFETY_FLAGS(v) & ~(m))

#else
/* All uses of the safety flags field must be guarded by #if SAFETY_CHECKS.
   If you see BUG__SAFETY_CHECKS_IS_FALSE, it means you forgot that guard.  */
#define IDA_SAFETY_FLAGS BUG__SAFETY_CHECKS_IS_FALSE
#define IDA_TEST_SAFETY_FLAGS BUG__SAFETY_CHECKS_IS_FALSE
#define IDA_WRITE_SAFETY_FLAGS BUG__SAFETY_CHECKS_IS_FALSE
#define IDA_MASK_SAFETY_FLAGS BUG__SAFETY_CHECKS_IS_FALSE
#define IDA_SET_SAFETY_FLAG BUG__SAFETY_CHECKS_IS_FALSE
#define IDA_UNSET_SAFETY_FLAG BUG__SAFETY_CHECKS_IS_FALSE
#endif

/* Flags kept in IDA_SAFETY_FLAGS.  There must be no more than 31 of them.
   All flags must be 0 for a new session.  */
enum ida_safety_flags {
    IDA_SAFETY_FLAG_SOLVING = 1,	/* IDASolve() has been called.  */
    IDA_SAFETY_FLAG_ID_SET = 2,		/* IDASetId() has been called. */

    /* Set of flags that IDAReInit() doesn't change.  */
    IDA_SAFETY_FLAG_REINIT_KEEPS =
      IDA_SAFETY_FLAG_ID_SET
};

enum ida_linear_solver_tag {
    /* all constructors have arguments */
#if IDA_ML_BIGARRAYS
    VARIANT_IDA_LINEAR_SOLVER_DENSE = 0,
    VARIANT_IDA_LINEAR_SOLVER_LAPACKDENSE,
    VARIANT_IDA_LINEAR_SOLVER_BAND,
    VARIANT_IDA_LINEAR_SOLVER_LAPACKBAND,
#endif
    VARIANT_IDA_LINEAR_SOLVER_SPGMR
#if !IDA_ML_BIGARRAYS
    = 0
#endif
    ,
    VARIANT_IDA_LINEAR_SOLVER_SPBCG,
    VARIANT_IDA_LINEAR_SOLVER_SPTFQMR
};

enum ida_integrator_stats_index {
    RECORD_IDA_INTEGRATOR_STATS_STEPS = 0,
    RECORD_IDA_INTEGRATOR_STATS_RES_EVALS,
    RECORD_IDA_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS,
    RECORD_IDA_INTEGRATOR_STATS_ERROR_TEST_FAILURES,
    RECORD_IDA_INTEGRATOR_STATS_LAST_INTERNAL_ORDER,
    RECORD_IDA_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER,
    RECORD_IDA_INTEGRATOR_STATS_INITIAL_STEP_SIZE,
    RECORD_IDA_INTEGRATOR_STATS_LAST_STEP_SIZE,
    RECORD_IDA_INTEGRATOR_STATS_NEXT_STEP_SIZE,
    RECORD_IDA_INTEGRATOR_STATS_INTERNAL_TIME,
    RECORD_IDA_INTEGRATOR_STATS_SIZE /* This has to come last. */
};

enum ida_jacobian_arg_index {
    RECORD_IDA_JACOBIAN_ARG_JAC_T = 0,
    RECORD_IDA_JACOBIAN_ARG_JAC_Y,
    RECORD_IDA_JACOBIAN_ARG_JAC_YP,
    RECORD_IDA_JACOBIAN_ARG_JAC_RES,
    RECORD_IDA_JACOBIAN_ARG_JAC_COEF,
    RECORD_IDA_JACOBIAN_ARG_JAC_TMP,
    RECORD_IDA_JACOBIAN_ARG_SIZE /* This has to come last. */
};

enum ida_solver_result_tag {
    VARIANT_IDA_SOLVER_RESULT_CONTINUE = 0,
    VARIANT_IDA_SOLVER_RESULT_ROOTSFOUND,
    VARIANT_IDA_SOLVER_RESULT_STOPTIMEREACHED,
};

enum ida_bandrange_index {
  RECORD_IDA_BANDRANGE_MUPPER = 0,
  RECORD_IDA_BANDRANGE_MLOWER,
  RECORD_IDA_BANDRANGE_SIZE
};

#endif /* _IDA_ML_H__ */
