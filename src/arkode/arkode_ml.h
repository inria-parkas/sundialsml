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

#ifndef _ARKODE_ML_H__
#define _ARKODE_ML_H__

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <caml/mlvalues.h>

/*
 * The session data structure is shared in four parts across the OCaml and C
 * heaps:
 *
 *           C HEAP                 .             OCAML HEAP
 * ---------------------------------.-----------------------------------
 *                                  .       (Program using Sundials/ML)
 *                                  .                            |
 *              +---------------+   .   +-------------------+    |
 *              | generational  +------>| weak ref : Weak.t |    |
 *              | global root   |   .   +-------------------+    |
 *              | (type: value) |   .   | ~ ~ ~ ~ ~ ~ ~ ~ ~ |    |
 *              +---------------+   .   +----------------+--+    |
 *                      /|\  /|\    .                    |       |
 *                       |    |     .                    |       |
 *                       |    |     .                   \|/     \|/
 *                       |    |     .                 +----------------+
 *                       |    |     .                 |  session       |
 *   +-------------+     |    |     .                 +----------------+
 *   | arkode_mem  |<---------------------------------+ arkode         |
 *   +-------------+     |    +-----------------------+ backref        |
 *   |    ...      |     |          .                 | nroots         |
 *   |ark_user_data+-----+          .                 | ls_callbacks   |
 *   |    ...      |                .                 | ...            |
 *   +-------------+                .                 +----------------+
 *
 *  * An arkode_mem structure is allocated by ARKodeInit for each session. It
 *    is the "C side" of the session data structure.
 *
 *  * The "OCaml side" of the session data structure is a record which contains
 *    a pointer to arkode_mem, several data fields, and the callback closures.
 *    It is returned to users of the library and used like any other abstract
 *    data type in OCaml.
 *
 *  * arkode_mem holds an indirect reference to the session record as user data
 *    (set by ARKodeSetUserData).  It cannot directly point to the record
 *    because the GC can change the record's address.  Instead, user data
 *    points to a global root which the GC updates whenever it relocates the
 *    session.  We cannot simply point ark_user_data to the weak reference and
 *    make it a global root, because we have no direct access to the members of
 *    arkode_mem.
 *
 *  * The global root points to a weak reference (a Weak.t of size 1) which
 *    points to the session record.  The root is destroyed when the session
 *    record is GC'ed -- note that if the root referenced the session record
 *    via a non-weak pointer the session would never be GC'ed, hence the root
 *    would never be destroyed either.
 *
 * 1. c_arkode_init() on the C side creates arkode_mem and the global root, and
 *    the OCaml side wraps that in a session record.  The OCaml side associates
 *    that record with a finalizer that unregisters the global root and frees
 *    all the C-side memory.
 *
 * 2. Callback functions (the right-hand side function, root function, etc)
 *    access the session through the user data.  This is the only way they can
 *    access the session.  The weak pointer is guaranteed to be alive during
 *    callback because the session record is alive.  The session record is
 *    captured in a C stack variable of type value when control enters the C
 *    stub that initiated the callback.
 *
 * 3. Other functions, like those that query integrator statistics, access the
 *    session record directly.
 *
 * 4. Eventually, when the user program abandons all references to the session
 *    record, the GC can reclaim the record because the only remaining direct
 *    reference to it is the weak pointer.
 *
 * NB: ark_user_data can't point directly to the session, unlike how in nvectors
 * (see nvector_ml.h) the backlink points directly to the payload.  This is
 * because the session contains the closure ls_callback, which may close over
 * the session in reasonable use cases.  In nvectors, by contrast, the payload
 * should be just an array of float's.
 */

void sunml_arkode_check_flag(const char *call, int flag, void *arkode_mem);
#if 400 <= SUNDIALS_LIB_VERSION
void sunml_arkode_check_ls_flag(const char *call, int flag);
#else
void sunml_arkode_check_dls_flag(const char *call, int flag);
void sunml_arkode_check_spils_flag(const char *call, int flag);
#endif

value sunml_arkode_make_jac_arg(sunrealtype t, N_Vector y, N_Vector fy, value tmp);
value sunml_arkode_make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

#define CHECK_FLAG(call, flag) if (flag != ARK_SUCCESS) \
				 sunml_arkode_check_flag(call, flag, NULL)
#if 400 <= SUNDIALS_LIB_VERSION
#define CHECK_LS_FLAG(call, flag) if (flag != ARKLS_SUCCESS) \
				 sunml_arkode_check_ls_flag(call, flag)
#else
#define CHECK_SPILS_FLAG(call, flag) if (flag != ARKLS_SUCCESS) \
				 sunml_arkode_check_ls_flag(call, flag)
#define CHECK_DLS_FLAG(call, flag) if (flag != ARKLS_SUCCESS) \
				 sunml_arkode_check_ls_flag(call, flag)
#endif

typedef enum {
    UNRECOVERABLE = 0,
    RECOVERABLE = 1
} recoverability;

int sunml_arkode_translate_exception (value session, value exn_result,
			        recoverability recoverable);

/* Check return value of a callback.  The common (no-error) case is
 * handled without setting up a new caml frame with CAMLparam.
 * Evaluates to:
 *   0 if result is not an exception
 *   1 if result is RecoverableFailure and recoverable == RECOVERABLE
 *  -1 otherwise, and records the exception in result in the session */
#define CHECK_EXCEPTION(session, result, recoverable)			\
    (Is_exception_result (result)					\
     ? sunml_arkode_translate_exception (session,				\
				   result = Extract_exception (result),	\
				   recoverable)				\
     : 0)

/* Interface with OCaml types */

/* Indices into the Arkode_*.session type.  This enum must be in the same order
 * as the session type's member declaration.  */
enum arkode_session_index {
    RECORD_ARKODE_SESSION_ARKODE = 0,
    RECORD_ARKODE_SESSION_BACKREF,
    RECORD_ARKODE_SESSION_NROOTS,
    RECORD_ARKODE_SESSION_CHECKVEC,
    RECORD_ARKODE_SESSION_USES_RESV,
    RECORD_ARKODE_SESSION_CONTEXT,
    RECORD_ARKODE_SESSION_EXN_TEMP,
    RECORD_ARKODE_SESSION_PROBLEM,
    RECORD_ARKODE_SESSION_RHSFN1,
    RECORD_ARKODE_SESSION_RHSFN2,
    RECORD_ARKODE_SESSION_ROOTSFN,
    RECORD_ARKODE_SESSION_ERRH,
    RECORD_ARKODE_SESSION_ERRW,
    RECORD_ARKODE_SESSION_RESW,
    RECORD_ARKODE_SESSION_ERROR_FILE,
    RECORD_ARKODE_SESSION_DIAG_FILE,
    RECORD_ARKODE_SESSION_ADAPTFN,
    RECORD_ARKODE_SESSION_STABFN,
    RECORD_ARKODE_SESSION_RESIZEFN,
    RECORD_ARKODE_SESSION_POSTSTEPFN,
    RECORD_ARKODE_SESSION_POSTSTAGEFN,
    RECORD_ARKODE_SESSION_STAGEPREDICTFN,
    RECORD_ARKODE_SESSION_PREINNERFN,
    RECORD_ARKODE_SESSION_POSTINNERFN,
    RECORD_ARKODE_SESSION_PREINNERARRAY,
    RECORD_ARKODE_SESSION_LINSOLVER,
    RECORD_ARKODE_SESSION_LS_SOLVER,
    RECORD_ARKODE_SESSION_LS_CALLBACKS,
    RECORD_ARKODE_SESSION_LS_PRECFNS,
    RECORD_ARKODE_SESSION_MASS_SOLVER,
    RECORD_ARKODE_SESSION_MASS_CALLBACKS,
    RECORD_ARKODE_SESSION_MASS_PRECFNS,
    RECORD_ARKODE_SESSION_NLS_SOLVER,
    RECORD_ARKODE_SESSION_NLS_RHSFN,
    RECORD_ARKODE_SESSION_INNER_SESSION,
    RECORD_ARKODE_SESSION_SIZE,
};

#define ARKODE_MEM(v) (SUNML_MEM(v))
#define ARKODE_MEM_FROM_ML(v) \
    (ARKODE_MEM(Field((v), RECORD_ARKODE_SESSION_ARKODE)))
#define ARKODE_BACKREF_FROM_ML(v) \
    ((value *)(Field((v), RECORD_ARKODE_SESSION_BACKREF)))
#define ARKODE_NROOTS_FROM_ML(v) \
    Long_val(Field((v), RECORD_ARKODE_SESSION_NROOTS))
#define ARKODE_PROBLEM(v) \
    Int_val(Field((v), RECORD_ARKODE_SESSION_PROBLEM))
#define ARKODE_ROOTSFN_FROM_ML(v) \
    Field((v), RECORD_ARKODE_SESSION_ROOTSFN)
#define ARKODE_LS_CALLBACKS_FROM_ML(v) \
    Field((v), RECORD_ARKODE_SESSION_LS_CALLBACKS)
#define ARKODE_LS_PRECFNS_FROM_ML(v) \
    Field((v), RECORD_ARKODE_SESSION_LS_PRECFNS)
#define ARKODE_MASS_CALLBACKS_FROM_ML(v) \
    Field((v), RECORD_ARKODE_SESSION_MASS_CALLBACKS)
#define ARKODE_MASS_PRECFNS_FROM_ML(v) \
    Field((v), RECORD_ARKODE_SESSION_MASS_PRECFNS)

enum arkode_problem_type_tag {
  VARIANT_ARKODE_PROBLEM_TYPE_IMPLICIT_ONLY = 0,
  VARIANT_ARKODE_PROBLEM_TYPE_EXPLICIT_ONLY,
  VARIANT_ARKODE_PROBLEM_TYPE_IMPLICIT_AND_EXPLICIT,
};

enum arkode_spils_precfns_index {
  RECORD_ARKODE_SPILS_PRECFNS_PREC_SOLVE_FN = 0,
  RECORD_ARKODE_SPILS_PRECFNS_PREC_SETUP_FN,
  RECORD_ARKODE_SPILS_PRECFNS_SIZE
};

enum arkode_spils_mass_precfns_index {
  RECORD_ARKODE_SPILS_MASS_PRECFNS_PREC_SOLVE_FN = 0,
  RECORD_ARKODE_SPILS_MASS_PRECFNS_PREC_SETUP_FN,
  RECORD_ARKODE_SPILS_MASS_PRECFNS_SIZE
};

enum arkode_alternate_callbacks_index {
  RECORD_ARKODE_ALTERNATE_CALLBACKS_LINIT = 0,
  RECORD_ARKODE_ALTERNATE_CALLBACKS_LSETUP,
  RECORD_ARKODE_ALTERNATE_CALLBACKS_LSOLVE,
  RECORD_ARKODE_ALTERNATE_CALLBACKS_SIZE
};

enum arkode_alternate_mass_callbacks_index {
  RECORD_ARKODE_ALTERNATE_MASS_CALLBACKS_MINIT = 0,
  RECORD_ARKODE_ALTERNATE_MASS_CALLBACKS_MSETUP,
  RECORD_ARKODE_ALTERNATE_MASS_CALLBACKS_MSOLVE,
  RECORD_ARKODE_ALTERNATE_MASS_CALLBACKS_SIZE
};

enum arkode_alternate_lsetup_args_index {
    RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_CONV_FAIL = 0,
    RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_Y,
    RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_RHS,
    RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_TMP,
    RECORD_ARKODE_ALTERNATE_LSETUP_ARGS_SIZE
};

enum arkode_alternate_lsolve_args_index {
    RECORD_ARKODE_ALTERNATE_LSOLVE_ARGS_Y = 0,
    RECORD_ARKODE_ALTERNATE_LSOLVE_ARGS_RHS,
    RECORD_ARKODE_ALTERNATE_LSOLVE_ARGS_SIZE
};

enum arkode_bbd_precfns_index {
  RECORD_ARKODE_BBD_PRECFNS_LOCAL_FN = 0,
  RECORD_ARKODE_BBD_PRECFNS_COMM_FN,
  RECORD_ARKODE_BBD_PRECFNS_SIZE
};

enum arkode_bandrange_index {
  RECORD_ARKODE_BANDRANGE_MUPPER = 0,
  RECORD_ARKODE_BANDRANGE_MLOWER,
  RECORD_ARKODE_BANDRANGE_SIZE
};

enum arkode_solver_result_tag {
  VARIANT_ARKODE_SOLVER_RESULT_SUCCESS = 0,
  VARIANT_ARKODE_SOLVER_RESULT_ROOTSFOUND,
  VARIANT_ARKODE_SOLVER_RESULT_STOPTIMEREACHED,
};

enum arkode_ark_timestepper_stats_index {
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_EXP_STEPS = 0,
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_ACC_STEPS,
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_STEP_ATTEMPTS,
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_NUM_NFE_EVALS,
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_NUM_NFI_EVALS,
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_LIN_SOLV_SETUPS,
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_NUM_ERR_TEST_FAILS,
  RECORD_ARKODE_ARK_TIMESTEPPER_STATS_SIZE
};

enum arkode_erk_timestepper_stats_index {
  RECORD_ARKODE_ERK_TIMESTEPPER_STATS_EXP_STEPS = 0,
  RECORD_ARKODE_ERK_TIMESTEPPER_STATS_ACC_STEPS,
  RECORD_ARKODE_ERK_TIMESTEPPER_STATS_STEP_ATTEMPTS,
  RECORD_ARKODE_ERK_TIMESTEPPER_STATS_NUM_NF_EVALS,
  RECORD_ARKODE_ERK_TIMESTEPPER_STATS_NUM_ERR_TEST_FAILS,
  RECORD_ARKODE_ERK_TIMESTEPPER_STATS_SIZE
};

enum arkode_step_stats_index {
  RECORD_ARKODE_STEP_STATS_STEPS = 0,
  RECORD_ARKODE_STEP_STATS_ACTUAL_INIT_STEP,
  RECORD_ARKODE_STEP_STATS_LAST_STEP,
  RECORD_ARKODE_STEP_STATS_CURRENT_STEP,
  RECORD_ARKODE_STEP_STATS_CURRENT_TIME,
  RECORD_ARKODE_STEP_STATS_SIZE
};

enum arkode_jacobian_arg_index {
  RECORD_ARKODE_JACOBIAN_ARG_JAC_T = 0,
  RECORD_ARKODE_JACOBIAN_ARG_JAC_Y,
  RECORD_ARKODE_JACOBIAN_ARG_JAC_FY,
  RECORD_ARKODE_JACOBIAN_ARG_JAC_TMP,
  RECORD_ARKODE_JACOBIAN_ARG_SIZE
};

enum arkode_spils_solve_arg_index {
  RECORD_ARKODE_SPILS_SOLVE_ARG_RHS = 0,
  RECORD_ARKODE_SPILS_SOLVE_ARG_GAMMA,
  RECORD_ARKODE_SPILS_SOLVE_ARG_DELTA,
  RECORD_ARKODE_SPILS_SOLVE_ARG_LEFT,
  RECORD_ARKODE_SPILS_SOLVE_ARG_SIZE
};

enum arkode_spils_mass_solve_arg_index {
  RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_RHS = 0,
  RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_DELTA,
  RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_LEFT,
  RECORD_ARKODE_SPILS_MASS_SOLVE_ARG_SIZE
};

enum arkode_alternate_conv_fail_index {
  VARIANT_ARKODE_ALTERNATE_NO_FAILURES = 0,
  VARIANT_ARKODE_ALTERNATE_FAIL_BAD_J,
  VARIANT_ARKODE_ALTERNATE_FAIL_OTHER
};

enum arkode_bandblock_bandwidths_index {
  RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MUDQ = 0,
  RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MLDQ,
  RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MUKEEP,
  RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MLKEEP,
  RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_SIZE
};

enum arkode_butcher_table_index {
  RECORD_ARKODE_BUTCHER_TABLE_METHOD_ORDER = 0,
  RECORD_ARKODE_BUTCHER_TABLE_STAGES,
  RECORD_ARKODE_BUTCHER_TABLE_STAGE_VALUES,
  RECORD_ARKODE_BUTCHER_TABLE_STAGE_TIMES,
  RECORD_ARKODE_BUTCHER_TABLE_COEFFICIENTS,
  RECORD_ARKODE_BUTCHER_TABLE_EMBEDDING,
  RECORD_ARKODE_BUTCHER_TABLE_SIZE
};

enum arkode_sprk_methodid_tag {
  VARIANT_ARKODE_SPRK_METHODID_EULER_1_1 = 0,
  VARIANT_ARKODE_SPRK_METHODID_LEAPFROG_2_2,
  VARIANT_ARKODE_SPRK_METHODID_PSEUDO_LEAPFROG_2_2,
  VARIANT_ARKODE_SPRK_METHODID_RUTH_3_3,
  VARIANT_ARKODE_SPRK_METHODID_MCLACHLAN_2_2,
  VARIANT_ARKODE_SPRK_METHODID_MCLACHLAN_3_3,
  VARIANT_ARKODE_SPRK_METHODID_CANDY_ROZMUS_4_4,
  VARIANT_ARKODE_SPRK_METHODID_MCLACHLAN_4_4,
  VARIANT_ARKODE_SPRK_METHODID_MCLACHLAN_5_6,
  VARIANT_ARKODE_SPRK_METHODID_YOSHIDA_6_8,
  VARIANT_ARKODE_SPRK_METHODID_SUZUKI_UMENO_8_16,
  VARIANT_ARKODE_SPRK_METHODID_SOFRONIOU_10_36,
};

enum arkode_sprk_table_index {
  RECORD_ARKODE_SPRK_TABLE_METHOD_ORDER = 0,
  RECORD_ARKODE_SPRK_TABLE_STAGES,
  RECORD_ARKODE_SPRK_TABLE_COEFFICIENTS1,
  RECORD_ARKODE_SPRK_TABLE_COEFFICIENTS2,
  RECORD_ARKODE_SPRK_TABLE_SIZE
};

enum arkode_adaptivity_args_index {
  RECORD_ARKODE_ADAPTIVITY_ARGS_H1 = 0,
  RECORD_ARKODE_ADAPTIVITY_ARGS_H2,
  RECORD_ARKODE_ADAPTIVITY_ARGS_H3,
  RECORD_ARKODE_ADAPTIVITY_ARGS_E1,
  RECORD_ARKODE_ADAPTIVITY_ARGS_E2,
  RECORD_ARKODE_ADAPTIVITY_ARGS_E3,
  RECORD_ARKODE_ADAPTIVITY_ARGS_Q,
  RECORD_ARKODE_ADAPTIVITY_ARGS_P,
  RECORD_ARKODE_ADAPTIVITY_ARGS_SIZE
};

enum arkode_adaptivity_params_index {
  RECORD_ARKODE_ADAPTIVITY_PARAMS_KS = 0,
  RECORD_ARKODE_ADAPTIVITY_PARAMS_METHOD_ORDER,
  RECORD_ARKODE_ADAPTIVITY_PARAMS_SIZE
};

enum arkode_adaptivity_method_tag {
  VARIANT_ARKODE_ADAPTIVITY_METHOD_PIDCONTROLLER = 0,
  VARIANT_ARKODE_ADAPTIVITY_METHOD_PICONTROLLER = 1,
  VARIANT_ARKODE_ADAPTIVITY_METHOD_ICONTROLLER = 2,
  VARIANT_ARKODE_ADAPTIVITY_METHOD_EXPLICITGUSTAFSSON = 3,
  VARIANT_ARKODE_ADAPTIVITY_METHOD_IMPLICITGUSTAFSSON = 4,
  VARIANT_ARKODE_ADAPTIVITY_METHOD_IMEXGUSTAFSSON = 5,
  VARIANT_ARKODE_ADAPTIVITY_METHOD_ADAPTIVITYFN,
};

enum arkode_predictor_method_tag {
  VARIANT_ARKODE_PREDICTOR_METHOD_TRIVIALPREDICTOR = 0,
  VARIANT_ARKODE_PREDICTOR_METHOD_MAXIMUMORDERPREDICTOR = 1,
  VARIANT_ARKODE_PREDICTOR_METHOD_VARIABLEORDERPREDICTOR = 2,
  VARIANT_ARKODE_PREDICTOR_METHOD_CUTOFFORDERPREDICTOR = 3,
  VARIANT_ARKODE_PREDICTOR_METHOD_BOOTSTRAPPREDICTOR = 4,
  VARIANT_ARKODE_PREDICTOR_METHOD_MINIMALCORRECTIONPREDICTOR = 5,
};

enum arkode_nonlin_system_data_index {
  RECORD_ARKODE_NONLIN_SYSTEM_DATA_TCUR = 0,
  RECORD_ARKODE_NONLIN_SYSTEM_DATA_ZPRED,
  RECORD_ARKODE_NONLIN_SYSTEM_DATA_ZI,
  RECORD_ARKODE_NONLIN_SYSTEM_DATA_FI,
  RECORD_ARKODE_NONLIN_SYSTEM_DATA_GAMMA,
  RECORD_ARKODE_NONLIN_SYSTEM_DATA_SDATA,
  RECORD_ARKODE_NONLIN_SYSTEM_DATA_SIZE
};

enum arkode_mri_coupling_index {
  RECORD_ARKODE_MRI_COUPLING_CPTR = 0,
  RECORD_ARKODE_MRI_COUPLING_NMAT,
  RECORD_ARKODE_MRI_COUPLING_STAGES,
  RECORD_ARKODE_MRI_COUPLING_METHOD_ORDER,
  RECORD_ARKODE_MRI_COUPLING_EMBEDDING,
  RECORD_ARKODE_MRI_COUPLING_EXPLICIT_MATRICES,
  RECORD_ARKODE_MRI_COUPLING_IMPLICIT_MATRICES,
  RECORD_ARKODE_MRI_COUPLING_ABSCISSAE,
  RECORD_ARKODE_MRI_COUPLING_SIZE
};

enum arkode_mri_coupling_table_tag {
  VARIANT_ARKODE_MRI_COUPLING_KW3 = 0,
  VARIANT_ARKODE_MRI_COUPLING_GARK_ERK33a,
  VARIANT_ARKODE_MRI_COUPLING_GARK_ERK45a,
  VARIANT_ARKODE_MRI_COUPLING_GARK_IRK21a,
  VARIANT_ARKODE_MRI_COUPLING_GARK_ESDIRK34a,
  VARIANT_ARKODE_MRI_COUPLING_GARK_ESDIRK46a,
  VARIANT_ARKODE_MRI_COUPLING_IMEX_GARK3a,
  VARIANT_ARKODE_MRI_COUPLING_IMEX_GARK3b,
  VARIANT_ARKODE_MRI_COUPLING_IMEX_GARK4,
};

enum arkode_mri_istepper_full_rhsfn_mode_tag {
  VARIANT_ARKODE_MRI_ISTEPPER_FULL_RHSFN_MODE_START = 0,
  VARIANT_ARKODE_MRI_ISTEPPER_FULL_RHSFN_MODE_END,
  VARIANT_ARKODE_MRI_ISTEPPER_FULL_RHSFN_MODE_OTHER,
};

enum arkode_mri_istepper_index {
  RECORD_ARKODE_MRI_ISTEPPER_RAWPTR = 0,
  RECORD_ARKODE_MRI_ISTEPPER_VAL,
  RECORD_ARKODE_MRI_ISTEPPER_ICHECKVEC,
  RECORD_ARKODE_MRI_ISTEPPER_SIZE
};

enum arkode_mri_istepper_tag {
  /* constructors with arguments */
  VARIANT_ARKODE_MRI_ISTEPPER_ARKSTEP = 0,
  VARIANT_ARKODE_MRI_ISTEPPER_CUSTOM,
  /* constructors without arguments */
  VARIANT_ARKODE_MRI_ISTEPPER_SUNDIALS = 0,
};

enum arkode_mri_istepper_callbacks_index {
  RECORD_ARKODE_MRI_ISTEPPER_CALLBACKS_EVOLVE_FN = 0,
  RECORD_ARKODE_MRI_ISTEPPER_CALLBACKS_FULL_RHS_FN,
  RECORD_ARKODE_MRI_ISTEPPER_CALLBACKS_RESET_FN,
  RECORD_ARKODE_MRI_ISTEPPER_CALLBACKS_SIZE
};

enum arkode_mri_istepper_forcing_data_index {
  RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_TSHIFT = 0,
  RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_TSCALE,
  RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_FORCING,
  RECORD_ARKODE_MRI_ISTEPPER_FORCING_DATA_SIZE
};

/* This enum must list exceptions in the same order as the call to
 * c_register_exns in arkode.ml.  */
enum arkode_exn_index {
    ARKODE_EXN_IllInput = 0,
    ARKODE_EXN_TooClose,
    ARKODE_EXN_TooMuchWork,
    ARKODE_EXN_TooMuchAccuracy,
    ARKODE_EXN_InnerStepFail,
    ARKODE_EXN_ErrFailure,
    ARKODE_EXN_ConvergenceFailure,
    ARKODE_EXN_LinearInitFailure,
    ARKODE_EXN_LinearSetupFailure,
    ARKODE_EXN_LinearSolveFailure,
    ARKODE_EXN_NonlinearInitFailure,
    ARKODE_EXN_NonlinearSetupFailure,
    ARKODE_EXN_NonlinearSetupRecoverable,
    ARKODE_EXN_NonlinearOperationError,
    ARKODE_EXN_MassInitFailure,
    ARKODE_EXN_MassSetupFailure,
    ARKODE_EXN_MassSolveFailure,
    ARKODE_EXN_MassMultFailure,
    ARKODE_EXN_RhsFuncFailure,
    ARKODE_EXN_FirstRhsFuncFailure,
    ARKODE_EXN_RepeatedRhsFuncFailure,
    ARKODE_EXN_UnrecoverableRhsFuncFailure,
    ARKODE_EXN_RootFuncFailure,
    ARKODE_EXN_PostprocStepFailure,
    ARKODE_EXN_BadK,
    ARKODE_EXN_BadT,
    ARKODE_EXN_VectorOpErr,
    ARKODE_EXN_SET_SIZE
};

#define ARKODE_EXN(name)     REGISTERED_EXN(ARKODE, name)
#define ARKODE_EXN_TAG(name) REGISTERED_EXN_TAG(ARKODE, name)


#endif
