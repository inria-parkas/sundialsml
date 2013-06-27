/***********************************************************************
 *                                                                     *
 *              Ocaml interface to Sundials CVODE solver               *
 *                                                                     *
 *           Timothy Bourke (INRIA) and Marc Pouzet (LIENS)            *
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under the terms of the GNU Library General Public License, with    *
 *  the special exception on linking described in file LICENSE.        *
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
    RECORD_IDA_SESSION_NEQS,
    RECORD_IDA_SESSION_NROOTS,
    RECORD_IDA_SESSION_ERRFILE,
    RECORD_IDA_SESSION_EXN_TEMP,
    RECORD_IDA_SESSION_RESFN,
    RECORD_IDA_SESSION_ROOTSFN,
    RECORD_IDA_SESSION_ERRH,
    RECORD_IDA_SESSION_ERRW,
    RECORD_IDA_SESSION_JACFN,
    RECORD_IDA_SESSION_BANDJACFN,
    RECORD_IDA_SESSION_PRESETUPFN,
    RECORD_IDA_SESSION_PRESOLVEFN,
    RECORD_IDA_SESSION_JACTIMESFN,
    RECORD_IDA_SESSION_SIZE	/* This has to come last. */
};

#define IDA_MEM_FROM_ML(v) ((void *)Field((v), RECORD_IDA_SESSION_MEM))
#define IDA_BACKREF_FROM_ML(v) ((value*)Field((v), RECORD_IDA_SESSION_NEQS))
#define IDA_NEQS_FROM_ML(v)    Long_val(Field((v), RECORD_IDA_SESSION_NEQS))
#define IDA_NROOTS_FROM_ML(v)  Long_val(Field((v), RECORD_IDA_SESSION_NROOTS))
#define IDA_ERRFILE_FROM_ML(v) Long_val(Field((v), RECORD_IDA_SESSION_ERRFILE))
#define IDA_RESFN_FROM_ML(v)   (Field((v), RECORD_IDA_SESSION_RESFN))
#define IDA_ROOTSFN_FROM_ML(v) (Field((v), RECORD_IDA_SESSION_ROOTSFN))
#define IDA_ERRH_FROM_ML(v)    (Field((v), RECORD_IDA_SESSION_ERRH))
#define IDA_ERRW_FROM_ML(v)    (Field((v), RECORD_IDA_SESSION_ERRW))
#define IDA_JACFN_FROM_ML(v)     (Field((v), RECORD_IDA_SESSION_JACFN))
#define IDA_BANDJACFN_FROM_ML(v) (Field((v), RECORD_IDA_SESSION_BANDJACFN))

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

enum ida_error_details_index {
    RECORD_IDA_ERROR_DETAILS_ERROR_CODE = 0,
    RECORD_IDA_ERROR_DETAILS_MODULE_NAME,
    RECORD_IDA_ERROR_DETAILS_FUNCTION_NAME,
    RECORD_IDA_ERROR_DETAILS_ERROR_MESSAGE,
    RECORD_IDA_ERROR_DETAILS_SIZE /* This has to come last. */
};

enum ida_linear_solver_tag {
    /* untagged: */
    VARIANT_IDA_LINEAR_SOLVER_DENSE = 0,
    VARIANT_IDA_LINEAR_SOLVER_LAPACKDENSE,
    /* tagged: */
    VARIANT_IDA_LINEAR_SOLVER_BAND = 0,
    VARIANT_IDA_LINEAR_SOLVER_LAPACKBAND,
    VARIANT_IDA_LINEAR_SOLVER_SPGMR,
    VARIANT_IDA_LINEAR_SOLVER_SPBCG,
    VARIANT_IDA_LINEAR_SOLVER_SPTFQMR,
};

enum ida_solver_result_tag {
    VARIANT_IDA_SOLVER_RESULT_CONTINUE = 0,
    VARIANT_IDA_SOLVER_RESULT_ROOTSFOUND,
    VARIANT_IDA_SOLVER_RESULT_STOPTIMEREACHED,
};

enum ida_gramschmidt_type {
    VARIANT_IDA_GRAMSCHMIDT_TYPE_MODIFIEDGS  = 0,
    VARIANT_IDA_GRAMSCHMIDT_TYPE_CLASSICALGS,
};

#endif /* _IDA_ML_H__ */
