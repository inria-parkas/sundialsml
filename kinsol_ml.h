/***********************************************************************
 *                                                                     *
 *               OCaml interface to (serial) Sundials                  *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

/*
 * This module defines all constants and functions which do not depend on
 * the representation of continuous state vectors, i.e., those that are
 * shared between the bigarray and nvector versions of ida_ml_nvec.
 *
 */

#ifndef _KINSOL_ML_H__
#define _KINSOL_ML_H__

#include "sundials_ml.h"

void kinsol_ml_check_flag(const char *call, int flag);

#define CHECK_FLAG(call, flag) if (flag != KIN_SUCCESS) \
				 kinsol_ml_check_flag(call, flag)

/* Indices into the Kinsol_*.session type.  This enum must be in the same order as
 * the session type's member declaration.  The session data structure is shared
 * in four parts across the OCaml and C heaps.  See cvode_ml.h for details.
 */
enum kinsol_index {
    RECORD_KINSOL_SESSION_MEM = 0,
    RECORD_KINSOL_SESSION_BACKREF,
    RECORD_KINSOL_SESSION_NEQS,
    RECORD_KINSOL_SESSION_ERRFILE,
    RECORD_KINSOL_SESSION_INFOFILE,
    RECORD_KINSOL_SESSION_EXN_TEMP,
    RECORD_KINSOL_SESSION_SYSFN,
    RECORD_KINSOL_SESSION_ERRH,
    RECORD_KINSOL_SESSION_INFOH,
#if KINSOL_ML_BIGARRAY
    RECORD_KINSOL_SESSION_JACFN,
    RECORD_KINSOL_SESSION_BANDJACFN,
#endif
    RECORD_KINSOL_SESSION_PRESETUPFN,
    RECORD_KINSOL_SESSION_PRESOLVEFN,
    RECORD_KINSOL_SESSION_JACTIMESFN,
    RECORD_KINSOL_SESSION_SIZE	/* This has to come last. */
};

#define KINSOL_MEM_FROM_ML(v)   ((void *)Field((v), RECORD_KINSOL_SESSION_MEM))
#define KINSOL_BACKREF_FROM_ML(v) \
    ((value *)(Field((v), RECORD_KINSOL_SESSION_BACKREF)))
#define KINSOL_NEQS_FROM_ML(v)  Long_val(Field((v), RECORD_KINSOL_SESSION_NEQS))
#define KINSOL_PRESETUPFN_FROM_ML(v) Field((v), RECORD_KINSOL_SESSION_PRESETUPFN)

enum kinsol_spils_prec_solve_arg_index {
  RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_USCALE = 0,
  RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_FSCALE,
  RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_SIZE
};

enum kinsol_spils_callbacks_index {
  RECORD_KINSOL_SPILS_CALLBACKS_PREC_SOLVE_FN = 0,
  RECORD_KINSOL_SPILS_CALLBACKS_PREC_SETUP_FN,
  RECORD_KINSOL_SPILS_CALLBACKS_JAC_TIMES_VEC_FN,
  RECORD_KINSOL_SPILS_CALLBACKS_SIZE
};

enum kinsol_linear_solver_tag {
#if KINSOL_ML_BIGARRAY
  VARIANT_KINSOL_LINEAR_SOLVER_DENSE = 0,
  VARIANT_KINSOL_LINEAR_SOLVER_LAPACKDENSE,
  VARIANT_KINSOL_LINEAR_SOLVER_BAND,
  VARIANT_KINSOL_LINEAR_SOLVER_LAPACKBAND,
#endif
  VARIANT_KINSOL_LINEAR_SOLVER_SPGMR
#if !KINSOL_ML_BIGARRAY
  = 0
#endif
  ,
  VARIANT_KINSOL_LINEAR_SOLVER_SPBCG,
  VARIANT_KINSOL_LINEAR_SOLVER_SPTFQMR,
};

enum kinsol_jacobian_arg_index {
  RECORD_KINSOL_JACOBIAN_ARG_JAC_U   = 0,
  RECORD_KINSOL_JACOBIAN_ARG_JAC_FU,
  RECORD_KINSOL_JACOBIAN_ARG_JAC_TMP,
  RECORD_KINSOL_JACOBIAN_ARG_SIZE
};

enum kinsol_print_level_tag {
  VARIANT_KINSOL_PRINT_LEVEL_NO_INFORMATION = 0,
  VARIANT_KINSOL_PRINT_LEVEL_SHOW_SCALED_NORMS,
  VARIANT_KINSOL_PRINT_LEVEL_SHOW_SCALED_DFNORM,
  VARIANT_KINSOL_PRINT_LEVEL_SHOW_GLOBAL_VALUES,
};

enum kinsol_eta_params_index {
  RECORD_KINSOL_ETA_PARAMS_EGAMMA = 0,
  RECORD_KINSOL_ETA_PARAMS_EALPHA,
  RECORD_KINSOL_ETA_PARAMS_SIZE
};

enum kinsol_eta_choice_tag {
  /* untagged: */
  VARIANT_KINSOL_ETA_CHOICE1 = 0,
  /* tagged: */
  VARIANT_KINSOL_ETA_CHOICE2 = 0,
  VARIANT_KINSOL_ETA_CONSTANT
};

enum kinsol_result_tag {
  VARIANT_KINSOL_RESULT_SUCCESS = 0,
  VARIANT_KINSOL_RESULT_INITIAL_GUESS_OK,
  VARIANT_KINSOL_RESULT_STOPPED_ON_STEP_TOL
};

#endif
