/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#ifndef _IDA_ML_H__
#define _IDA_ML_H__

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <caml/mlvalues.h>

#include "../sundials/sundials_ml.h"

void sunml_ida_check_flag(const char *call, int flag, void *ida_mem);
#if 400 <= SUNDIALS_LIB_VERSION
void sunml_ida_check_ls_flag(const char *call, int flag);
#else
void sunml_ida_check_dls_flag(const char *call, int flag);
void sunml_ida_check_spils_flag(const char *call, int flag);
#endif

value sunml_ida_make_jac_arg(realtype t, realtype coef, N_Vector y, N_Vector yp,
		       N_Vector res, value tmp);
value sunml_ida_make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
value sunml_ida_make_double_tmp(N_Vector tmp1, N_Vector tmp2);

value sunml_ida_last_lin_exception(void *ida_mem);

#define CHECK_FLAG(call, flag) if (flag != IDA_SUCCESS) \
				 sunml_ida_check_flag(call, flag, NULL)
#if 400 <= SUNDIALS_LIB_VERSION
#define CHECK_LS_FLAG(call, flag) if (flag != IDALS_SUCCESS) \
				 sunml_ida_check_ls_flag(call, flag)
#else
#define CHECK_SPILS_FLAG(call, flag) if (flag != IDASPILS_SUCCESS) \
				 sunml_ida_check_spils_flag(call, flag)
#define CHECK_DLS_FLAG(call, flag) if (flag != IDADLS_SUCCESS) \
				 sunml_ida_check_dls_flag(call, flag)
#endif

typedef enum {
    UNRECOVERABLE = 0,
    RECOVERABLE = 1
} recoverability;

int sunml_ida_translate_exception (value session, value exn_result,
			     recoverability recoverable);

/* Check return value of a callback.  The common (no-error) case is
 * handled without setting up a new caml frame with CAMLparam.
 * Evaluates to:
 *   0 if result is not an exception
 *   1 if result is RecoverableFailure and recoverable == RECOVERABLE
 *  -1 otherwise, and records the exception in result in the session */
#define CHECK_EXCEPTION(session, result, recoverable)		\
    (Is_exception_result (result)				\
     ? sunml_ida_translate_exception (session, result, recoverable)	\
     : 0)

/* Indices into the Ida_*.session type.  This enum must be in the same order as
 * the session type's member declaration.  The session data structure is shared
 * in four parts across the OCaml and C heaps.  See cvode_ml.h for details.
 */
enum ida_index {
    RECORD_IDA_SESSION_MEM = 0,
    RECORD_IDA_SESSION_BACKREF,
    RECORD_IDA_SESSION_NROOTS,
    RECORD_IDA_SESSION_CHECKVEC,
    RECORD_IDA_SESSION_EXN_TEMP,
    RECORD_IDA_SESSION_ID_SET,
    RECORD_IDA_SESSION_RESFN,
    RECORD_IDA_SESSION_ROOTSFN,
    RECORD_IDA_SESSION_ERRH,
    RECORD_IDA_SESSION_ERRW,
    RECORD_IDA_SESSION_LS_SOLVER,
    RECORD_IDA_SESSION_LS_CALLBACKS,
    RECORD_IDA_SESSION_LS_PRECFNS,
    RECORD_IDA_SESSION_NLS_SOLVER,
    RECORD_IDA_SESSION_SENSEXT,
    RECORD_IDA_SESSION_SIZE	/* This has to come last. */
};

#define IDA_MEM(v) (*(void **)Data_custom_val(v))
#define IDA_MEM_FROM_ML(v) (IDA_MEM(Field((v), RECORD_IDA_SESSION_MEM)))
#define IDA_BACKREF_FROM_ML(v) ((value*)Field((v), RECORD_IDA_SESSION_BACKREF))
#define IDA_NEQS_FROM_ML(v)    Long_val(Field((v), RECORD_IDA_SESSION_NEQS))
#define IDA_NROOTS_FROM_ML(v)  Long_val(Field((v), RECORD_IDA_SESSION_NROOTS))
#define IDA_ERRFILE_FROM_ML(v) Long_val(Field((v), RECORD_IDA_SESSION_ERRFILE))
#define IDA_RESFN_FROM_ML(v)   (Field((v), RECORD_IDA_SESSION_RESFN))
#define IDA_ROOTSFN_FROM_ML(v) (Field((v), RECORD_IDA_SESSION_ROOTSFN))
#define IDA_ERRH_FROM_ML(v)    (Field((v), RECORD_IDA_SESSION_ERRH))
#define IDA_ERRW_FROM_ML(v)    (Field((v), RECORD_IDA_SESSION_ERRW))
#define IDA_LS_CALLBACKS_FROM_ML(v)     (Field((v), RECORD_IDA_SESSION_LS_CALLBACKS))
#define IDA_LS_PRECFNS_FROM_ML(v)     (Field((v), RECORD_IDA_SESSION_LS_PRECFNS))
#define IDA_SENSEXT_FROM_ML(v)     (Field(Field((v), RECORD_IDA_SESSION_SENSEXT), 0))

enum ida_spils_precfns_index {
    RECORD_IDA_SPILS_PRECFNS_PREC_SOLVE_FN = 0,
    RECORD_IDA_SPILS_PRECFNS_PREC_SETUP_FN,
    RECORD_IDA_SPILS_PRECFNS_SIZE
};

enum ida_bbd_precfns_index {
  RECORD_IDA_BBD_PRECFNS_LOCAL_FN = 0,
  RECORD_IDA_BBD_PRECFNS_COMM_FN,
  RECORD_IDA_BBD_PRECFNS_SIZE
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
    VARIANT_IDA_SOLVER_RESULT_SUCCESS = 0,
    VARIANT_IDA_SOLVER_RESULT_ROOTSFOUND,
    VARIANT_IDA_SOLVER_RESULT_STOPTIMEREACHED,
};

enum ida_bandrange_index {
  RECORD_IDA_BANDRANGE_MUPPER = 0,
  RECORD_IDA_BANDRANGE_MLOWER,
  RECORD_IDA_BANDRANGE_SIZE
};

enum ida_bandwidths_index {
  RECORD_IDA_BANDWIDTHS_MUDQ = 0,
  RECORD_IDA_BANDWIDTHS_MLDQ,
  RECORD_IDA_BANDWIDTHS_MUKEEP,
  RECORD_IDA_BANDWIDTHS_MLKEEP,
  RECORD_IDA_BANDWIDTHS_SIZE
};

enum ida_nonlin_system_data_index {
  RECORD_IDA_NONLIN_SYSTEM_DATA_TN = 0,
  RECORD_IDA_NONLIN_SYSTEM_DATA_YYPRED,
  RECORD_IDA_NONLIN_SYSTEM_DATA_YPPRED,
  RECORD_IDA_NONLIN_SYSTEM_DATA_YYN,
  RECORD_IDA_NONLIN_SYSTEM_DATA_YPN,
  RECORD_IDA_NONLIN_SYSTEM_DATA_RES,
  RECORD_IDA_NONLIN_SYSTEM_DATA_CJ,
  RECORD_IDA_NONLIN_SYSTEM_DATA_SIZE
};

/* This enum must list exceptions in the same order as the call to
 * c_register_exns in ida.ml.  */
enum ida_exn_index {
    IDA_EXN_IllInput = 0,
    IDA_EXN_TooMuchWork,
    IDA_EXN_TooMuchAccuracy,
    IDA_EXN_ErrFailure,
    IDA_EXN_ConvergenceFailure,
    IDA_EXN_LinearInitFailure,
    IDA_EXN_LinearSetupFailure,
    IDA_EXN_LinearSolveFailure,
    IDA_EXN_NonlinearFailure,
    IDA_EXN_NonlinearInitFailure,
    IDA_EXN_NonlinearSetupFailure,
    IDA_EXN_NonlinearSetupRecoverable,
    IDA_EXN_ResFuncFailure,
    IDA_EXN_FirstResFuncFailure,
    IDA_EXN_RepeatedResFuncFailure,
    IDA_EXN_RootFuncFailure,
    IDA_EXN_ConstraintFailure,

    IDA_EXN_LinesearchFailure,
    IDA_EXN_NoRecovery,
    IDA_EXN_BadEwt,

    IDA_EXN_BadK,
    IDA_EXN_BadT,
    IDA_EXN_VectorOpErr,
    IDA_EXN_SET_SIZE,
};

#define IDA_EXN(name)     REGISTERED_EXN(IDA, name)
#define IDA_EXN_TAG(name) REGISTERED_EXN_TAG(IDA, name)


#endif /* _IDA_ML_H__ */
