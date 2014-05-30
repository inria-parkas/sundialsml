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

// TODO: document how it is all supposed to work

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
    RECORD_CVODES_FWD_SESSION_SENSPVALS,
    RECORD_CVODES_FWD_SESSION_SENSRHSFN,
    RECORD_CVODES_FWD_SESSION_SENSRHSFN1,
    RECORD_CVODES_FWD_SESSION_QUADSENSRHSFN,
    RECORD_CVODES_FWD_SESSION_SIZE
}

enum cvodes_bwd_session_index {
    RECORD_CVODES_BWD_SESSION_PARENT = 0,
    RECORD_CVODES_BWD_SESSION_WHICH,
    RECORD_CVODES_BWD_SESSION_BRHSFN,
    RECORD_CVODES_BWD_SESSION_BRHSFN1,
    RECORD_CVODES_BWD_SESSION_BQUADRHSFN,
    RECORD_CVODES_BWD_SESSION_BQUADRHSFN1,
    RECORD_CVODES_BWD_SESSION_BPRESETUPFN,
    RECORD_CVODES_BWD_SESSION_BPRESOLVEFN,
    RECORD_CVODES_BWD_SESSION_BJACTIMESFN,
    RECORD_CVODES_BWD_SESSION_SIZE
}

// TODO: add accessor macros

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
    RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_FY,
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
    RECORD_CVODE_SPILS_PARAMS_PRETYPE,
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
