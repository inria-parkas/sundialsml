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

#ifndef _CVODES_ML_H__
#define _CVODES_ML_H__

/* The sessions of CVODES work almost exactly as described in cvode_ml.h.
 *
 * 1. Standard (quadrature, sensitivity, sensitivity/quadrature,
 *    adjoint-forward) sessions are created by the basic cvode code. Extra
 *    details are added and accessed through the sensext field.
 *
 * 2. In CVodes backward sessions are 'piggy-backed' onto parent sessions;
 *    the parent maintains a linked list and each session is identified by a
 *    'which-id' (basically its order in the list). For all set* functions,
 *    one must pass the pair of parent and which-id, but for many get*
 *    functions, one passes the underlying cvode_mem record associated with
 *    the backward session (normally accessed through CVodeGetAdjCVodeBmem).
 *
 *    In the OCaml interface to CVodes we create a session object for each
 *    backward session that works almost exactly like the basic CVode
 *    sessions (as described in cvode_ml.h); i.e., we set the user-data for
 *    the backward session (CVodeSetUserDataB) to a weak pointer to a
 *    session object on the OCaml side. The sensext field is used to track
 *    parent, which-id, and callbacks. The first two are used for set*
 *    interface functions and the session itself is used for get* interface
 *    functions.
 *
 *    There are two small differences. First, on the OCaml side, we type
 *    each such session as a 'bsession' (rather than as a 'session'), even
 *    though the underlying representation is compatible with a standard
 *    session (excepting the sensext field). Second, the finalizer neither
 *    calls CVodeFree, since the parent takes care of this, nor closes an
 *    error file, since one cannot be created (backward sessions inherit
 *    from the parent at the time of their creation). In the parent we
 *    maintain a list of child (backward) sessions and in each child we
 *    maintain a reference to the parent so that all are garbage collected
 *    at once (such loops are not a problem since they are completely on the
 *    OCaml side). If we did not maintain a link from the parent to a child,
 *    a child could be garbage collected and then the callback (triggered by
 *    a CVodeB on the parent) would fail.
 *
 */

void cvodes_ml_check_flag(const char *call, int flag);

// NB: overrides CHECK_FLAG macro in cvode_ml.h
#define SCHECK_FLAG(call, flag) if (flag != CV_SUCCESS) \
				  cvodes_ml_check_flag(call, flag)

enum cvodes_fwd_session_index {
    RECORD_CVODES_FWD_SESSION_QUADRHSFN = 0,
    RECORD_CVODES_FWD_SESSION_NUMSENSITIVITIES,
    RECORD_CVODES_FWD_SESSION_SENSARRAY1,
    RECORD_CVODES_FWD_SESSION_SENSARRAY2,
    RECORD_CVODES_FWD_SESSION_SENSPVALS,
    RECORD_CVODES_FWD_SESSION_SENSRHSFN,
    RECORD_CVODES_FWD_SESSION_SENSRHSFN1,
    RECORD_CVODES_FWD_SESSION_QUADSENSRHSFN,
    RECORD_CVODES_FWD_SESSION_BSESSIONS,
    RECORD_CVODES_FWD_SESSION_SIZE
};

#define CVODES_SENSARRAY1_FROM_EXT(v)    Field((v), RECORD_CVODES_FWD_SESSION_SENSARRAY1)
#define CVODES_SENSARRAY2_FROM_EXT(v)    Field((v), RECORD_CVODES_FWD_SESSION_SENSARRAY2)
#define CVODES_SENSRHSFN_FROM_EXT(v)     Field((v), RECORD_CVODES_FWD_SESSION_SENSRHSFN)
#define CVODES_QUADSENSRHSFN_FROM_EXT(v) Field((v), RECORD_CVODES_FWD_SESSION_QUADSENSRHSFN)

enum cvodes_bwd_session_index {
    RECORD_CVODES_BWD_SESSION_PARENT = 0,
    RECORD_CVODES_BWD_SESSION_WHICH,
    RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES,
    RECORD_CVODES_BWD_SESSION_BSENSARRAY,
    RECORD_CVODES_BWD_SESSION_BRHSFN,
    RECORD_CVODES_BWD_SESSION_BRHSFN_SENS,
    RECORD_CVODES_BWD_SESSION_BQUADRHSFN,
    RECORD_CVODES_BWD_SESSION_BQUADRHSFN_SENS,
    RECORD_CVODES_BWD_SESSION_SIZE
};

#define CVODES_PARENT_FROM_BEXT(v) Field((v), RECORD_CVODES_BWD_SESSION_PARENT)
#define CVODES_WHICH_FROM_BEXT(v)  Field((v), RECORD_CVODES_BWD_SESSION_WHICH)

#define CVODES_BSENSARRAY_FROM_EXT(v)      Field((v), RECORD_CVODES_BWD_SESSION_BSENSARRAY)
#define CVODES_BRHSFN_FROM_EXT(v)          Field((v), RECORD_CVODES_BWD_SESSION_BRHSFN)
#define CVODES_BQUADRHSFN_FROM_EXT(v)      Field((v), RECORD_CVODES_BWD_SESSION_BQUADRHSFN)
#define CVODES_BRHSFN_SENS_FROM_EXT(v)     Field((v), RECORD_CVODES_BWD_SESSION_BRHSFN_SENS)
#define CVODES_BQUADRHSFN_SENS_FROM_EXT(v) Field((v), RECORD_CVODES_BWD_SESSION_BQUADRHSFN_SENS)

enum cvodes_sens_method {
    VARIANT_CVODES_SENS_METHOD_SIMULTANEOUS = 0,
    VARIANT_CVODES_SENS_METHOD_STAGGERED,
    VARIANT_CVODES_SENS_METHOD_STAGGERED1
};

enum cvodes_sens_params_index {
    RECORD_CVODES_SENS_PARAMS_PVALS   = 0,
    RECORD_CVODES_SENS_PARAMS_PBAR,
    RECORD_CVODES_SENS_PARAMS_PLIST,
    RECORD_CVODES_SENS_PARAMS_SIZE
};

enum cvodes_sens_dq_method {
    VARIANT_CVODES_SENS_DQ_METHOD_CENTERED = 0,
    VARIANT_CVODES_SENS_DQ_METHOD_FORWARD
};

enum cvodes_sens_stats_index {
    RECORD_CVODES_SENS_STATS_NUM_RHS_EVALS = 0,
    RECORD_CVODES_SENS_STATS_NUM_SENS_EVALS,
    RECORD_CVODES_SENS_STATS_NUM_ERR_TEST_FAILS,
    RECORD_CVODES_SENS_STATS_NUM_LIN_SOLV_SETUPS,
    RECORD_CVODES_SENS_STATS_SIZE
};

enum cvodes_adj_interpolation {
    VARIANT_CVODES_ADJ_INTERPOLATION_POLYNOMIAL = 0,
    VARIANT_CVODES_ADJ_INTERPOLATION_HERMITE
};

enum cvodes_adj_jacobian_arg_index {
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_T   = 0,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_Y,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_YB,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_FYB,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_TMP,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_SIZE
};

enum cvodes_adj_bandrange_index {
    RECORD_CVODES_ADJ_BANDRANGE_MUPPER = 0,
    RECORD_CVODES_ADJ_BANDRANGE_MLOWER,
    RECORD_CVODES_ADJ_BANDRANGE_SIZE
};

enum cvodes_adj_solve_arg_index {
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_RVEC   = 0,
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_GAMMA,
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_DELTA,
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_LR,
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_SIZE
};

enum cvodes_enum_index {
    CVODES_EXN_QuadNotInitialized,
    CVODES_EXN_QuadRhsFuncFailure,
    CVODES_EXN_FirstQuadRhsFuncFailure,
    CVODES_EXN_RepeatedQuadRhsFuncFailure,
    CVODES_EXN_UnrecoverableQuadRhsFuncFailure,

    CVODES_EXN_SensNotInitialized,
    CVODES_EXN_SensRhsFuncFailure,
    CVODES_EXN_FirstSensRhsFuncFailure,
    CVODES_EXN_RepeatedSensRhsFuncFailure,
    CVODES_EXN_UnrecoverableSensRhsFuncFailure,
    CVODES_EXN_BadSensIdentifier,

    CVODES_EXN_QuadSensNotInitialized,
    CVODES_EXN_QuadSensRhsFuncFailure,
    CVODES_EXN_FirstQuadSensRhsFuncFailure,
    CVODES_EXN_RepeatedQuadSensRhsFuncFailure,
    CVODES_EXN_UnrecoverableQuadSensRhsFuncFailure,

    CVODES_EXN_AdjointNotInitialized,
    CVODES_EXN_NoForwardCall,
    CVODES_EXN_ForwardReinitFailure,
    CVODES_EXN_ForwardFailure,
    CVODES_EXN_NoBackwardProblem,
    CVODES_EXN_BadFinalTime,
    CVODES_EXN_BadOutputTime,
    CVODES_EXN_SET_SIZE
};

#define CVODES_EXN(name) (Field (Field (sundials_ml_exn_table,	\
					CVODES_EXN_SET),	\
				 CVODES_EXN_ ## name))

#endif
