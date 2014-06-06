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

#ifndef _CVODES_ML_H__
#define _CVODES_ML_H__

#include "cvode_ml.h"

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
 *    from the parent at the time of their creation). Note that a parent
 *    will always be garbage collected after all of its children (since the
 *    latter hold references to the former and not vice-versa).
 */

void cvodes_ml_check_flag(const char *call, int flag);

// NB: overrides CHECK_FLAG macro in cvode_ml.h
#define CHECK_FLAG(call, flag) if (flag != CV_SUCCESS) \
				 cvodes_ml_check_flag(call, flag)

enum cvodes_sensext {
    VARIANT_CVODES_SENSEXT_NONE = 0,
    VARIANT_CVODES_SENSEXT_FWD,
    VARIANT_CVODES_SENSEXT_BWD,

enum cvodes_fwd_session_index {
    RECORD_CVODES_FWD_SESSION_QUADRHSFN = 0,
    RECORD_CVODES_FWD_SESSION_SENSARRAY1,
    RECORD_CVODES_FWD_SESSION_SENSARRAY2,
    RECORD_CVODES_FWD_SESSION_SENSPVALS,
    RECORD_CVODES_FWD_SESSION_SENSRHSFN,
    RECORD_CVODES_FWD_SESSION_SENSRHSFN1,
    RECORD_CVODES_FWD_SESSION_QUADSENSRHSFN,
    RECORD_CVODES_FWD_SESSION_SIZE
}

#define CVODES_SENSARRAY1_FROM_EXT(v)    Field((v), RECORD_CVODES_FWD_SESSION_SENSARRAY1)
#define CVODES_SENSARRAY2_FROM_EXT(v)    Field((v), RECORD_CVODES_FWD_SESSION_SENSARRAY2)
#define CVODES_SENSRHSFN_FROM_EXT(v)     Field((v), RECORD_CVODES_FWD_SESSION_SENSRHSFN)
#define CVODES_QUADSENSRHSFN_FROM_EXT(v) Field((v), RECORD_CVODES_FWD_SESSION_QUADSENSRHSFN)

enum cvodes_bwd_session_index {
    RECORD_CVODES_BWD_SESSION_PARENT = 0,
    RECORD_CVODES_BWD_SESSION_WHICH,
    RECORD_CVODES_BWD_SESSION_BSENSARRAY,
    RECORD_CVODES_BWD_SESSION_BRHSFN,
    RECORD_CVODES_BWD_SESSION_BRHSFN1,
    RECORD_CVODES_BWD_SESSION_BQUADRHSFN,
    RECORD_CVODES_BWD_SESSION_BQUADRHSFN1,
    RECORD_CVODES_BWD_SESSION_BPRESETUPFN,
    RECORD_CVODES_BWD_SESSION_BPRESOLVEFN,
    RECORD_CVODES_BWD_SESSION_BJACTIMESFN,
    RECORD_CVODES_BWD_SESSION_SIZE
}

#define CVODES_BPARENT_FROM_EXT(v) Field((v), RECORD_CVODES_BWD_SESSION_PARENT)
#define CVODES_BWHICH_FROM_EXT(v)  Field((v), RECORD_CVODES_BWD_SESSION_WHICH)

#define CVODES_BSENSARRAY_FROM_EXT(v)  Field((v), RECORD_CVODES_BWD_SESSION_BSENSARRAY)
#define CVODES_BRHSFN1_FROM_EXT(v)     Field((v), RECORD_CVODES_BWD_SESSION_BRHSFN1)
#define CVODES_BQUADRHSFN1_FROM_EXT(v) Field((v), RECORD_CVODES_BWD_SESSION_BQUADRHSFN1)
#define CVODES_BPRESETUPFN_FROM_EXT(v) Field((v), RECORD_CVODES_BWD_SESSION_BPRESETUPFN)

enum cvodes_quad_tol_tag {
    VARIANT_CVODES_QUAD_TOL_NONE = 0,
    VARIANT_CVODES_QUAD_TOL_SS,
    VARIANT_CVODES_QUAD_TOL_SV,
}

enum cvodes_sens_tol_tag {
    VARIANT_CVODES_SENS_TOL_SS = 0,
    VARIANT_CVODES_SENS_TOL_SV,
    VARIANT_CVODES_SENS_TOL_EE,
}

enum cvodes_sens_method {
    VARIANT_CVODES_SENS_METHOD_SIMULTANEOUS = 0,
    VARIANT_CVODES_SENS_METHOD_STAGGERED,
    VARIANT_CVODES_SENS_METHOD_STAGGERED1,
}

enum cvodes_sens_rhsfn {
    VARIANT_CVODES_SENS_RHSFN_ALLATONCE = 0,
    VARIANT_CVODES_SENS_RHSFN_ONEBYONE,
}

enum cvodes_sens_params_index {
    RECORD_CVODES_SENS_PARAMS_PVALS   = 0,
    RECORD_CVODES_SENS_PARAMS_PBAR,
    RECORD_CVODES_SENS_PARAMS_PLIST,
    RECORD_CVODES_SENS_PARAMS_SIZE
}

enum cvodes_sens_dq_method {
    VARIANT_CVODES_SENS_DQ_METHOD_CENTERED = 0,
    VARIANT_CVODES_SENS_DQ_METHOD_FORWARD,
}

enum cvodes_sens_stats_index {
    RECORD_CVODES_SENS_STATS_NUM_RHS_EVALS = 0,
    RECORD_CVODES_SENS_STATS_NUM_SENS_EVALS,
    RECORD_CVODES_SENS_STATS_NUM_ERR_TEST_FAILS,
    RECORD_CVODES_SENS_STATS_NUM_LIN_SOLV_SETUPS,
    RECORD_CVODES_SENS_STATS_SIZE
}

enum cvodes_sens_nonlin_stats_index {
    RECORD_CVODES_SENS_NONLIN_STATS_NUM_SOLV_ITERS = 0,
    RECORD_CVODES_SENS_NONLIN_STATS_NUM_CONV_FAILS,
    RECORD_CVODES_SENS_NONLIN_STATS_SIZE
}

enum cvodes_sens_quad_tol_tag {
    VARIANT_CVODES_SENS_QUAD_TOL_NONE = 0,
    VARIANT_CVODES_SENS_QUAD_TOL_SS,
    VARIANT_CVODES_SENS_QUAD_TOL_SV,
    VARIANT_CVODES_SENS_QUAD_TOL_EE,
}

enum cvodes_adj_interpolation {
    VARIANT_CVODES_ADJ_INTERPOLATION_POLYNOMIAL = 0,
    VARIANT_CVODES_ADJ_INTERPOLATION_HERMITE,
}

enum cvodes_adj_brhsfn {
    VARIANT_CVODES_ADJ_BRHSFN_BASIC = 0,
    VARIANT_CVODES_ADJ_BRHSFN_WITHSENS,
}

enum cvodes_adj_jacobian_arg_index {
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_T   = 0,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_Y,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_YB,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_FYB,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_TMP,
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_SIZE
}

enum cvodes_adj_iter_tag {
    VARIANT_CVODES_ADJ_ITER_TAG_NEWTON = 0,
    VARIANT_CVODES_ADJ_ITER_TAG_FUNCTIONAL,
};

enum cvodes_linear_solver_tag {
    /* untagged: */
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_DIAG = 0,
    /* tagged: */
#if CVODE_ML_BIGARRAYS
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_DENSE = 0,
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_LAPACKDENSE,
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_BAND,
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_LAPACKBAND,
#endif
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_SPGMR
#if !CVODE_ML_BIGARRAYS
    = 0
#endif
    ,
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_SPBCG,
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_SPTFQMR,
#if CVODE_ML_BIGARRAYS
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_BANDED_SPGMR,
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_BANDED_SPBCG,
    VARIANT_CVODES_ADJ_LINEAR_SOLVER_BANDED_SPTFQMR,
#endif
};

enum cvodes_adj_bandrange_index {
    RECORD_CVODES_ADJ_BANDRANGE_MUPPER = 0,
    RECORD_CVODES_ADJ_BANDRANGE_MLOWER,
    RECORD_CVODES_ADJ_BANDRANGE_SIZE
};

enum cvodes_adj_spils_params_index {
    RECORD_CVODE_SPILS_PARAMS_MAXL = 0,
    RECORD_CVODE_SPILS_PARAMS_PREC_TYPE,
    RECORD_CVODE_SPILS_PARAMS_SIZE
};

enum cvodes_adj_spils_callbacks_index {
    RECORD_CVODES_ADJ_SPILS_CALLBACKS_PREC_SOLVE_FN = 0,
    RECORD_CVODES_ADJ_SPILS_CALLBACKS_PREC_SETUP_FN,
    RECORD_CVODES_ADJ_SPILS_CALLBACKS_JAC_TIMES_VEC_FN,
    RECORD_CVODES_ADJ_SPILS_CALLBACKS_SIZE
};

enum cvodes_adj_solve_arg_index {
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_RVECB   = 0,
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_GAMMAB,
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_DELTAB,
    RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_SIZE
};

enum cvodes_adj_tol_tag {
    VARIANT_CVODES_ADJ_TOL_SS = 0,
    VARIANT_CVODES_ADJ_TOL_SV,
}

enum cvodes_adj_quad_bquadrhsfn {
    VARIANT_CVODES_ADJ_BRHSFN_BASIC = 0,
    VARIANT_CVODES_ADJ_BRHSFN_WITHSENS,
}

enum cvodes_adj_quad_tol_tag {
    VARIANT_CVODES_ADJ_TOL_NONE = 0,
    VARIANT_CVODES_ADJ_TOL_SS,
    VARIANT_CVODES_ADJ_TOL_SV,
}

#endif
