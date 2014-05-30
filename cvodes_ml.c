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

/* Interface functions that do not involve NVectors. */

// TODO: remove unneeded header files, add needed ones
// TODO: standardize on function names

#include <cvodes/cvodes.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

/* linear solvers */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_direct.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_spils.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <cvodes/cvodes_lapack.h>
#endif

#include "spils_ml.h"
#include "cvodes_ml.h"
#include "sundials_ml.h"

// TODO: must ensure that these exceptions are registered...
void cvodes_ml_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	    || flag == CV_ROOT_RETURN
	    || flag == CV_TSTOP_RETURN) return;

    switch (flag) {
	case CV_TOO_MUCH_WORK:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchWork"));

	case CV_TOO_MUCH_ACC:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchAccuracy"));

	case CV_ERR_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ErrFailure"));

	case CV_CONV_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ConvergenceFailure"));

	/* * */

	case CV_LINIT_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_LinearInitFailure"));

	case CV_LSETUP_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_LinearSetupFailure"));

	case CV_LSOLVE_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_LinearSolveFailure"));

	case CV_RHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_RhsFuncFailure"));

	case CV_FIRST_RHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvode_FirstRhsFuncErr"));

	case CV_REPTD_RHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvode_RepeatedRhsFuncErr"));

	case CV_UNREC_RHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvode_UnrecoverableRhsFuncErr"));

	case CV_RTFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvode_RootFuncFailure"));

	/* * */

	case CV_ILL_INPUT:
	    caml_raise_constant(*caml_named_value("cvode_IllInput"));
	
	case CV_BAD_K:
	    caml_raise_constant(*caml_named_value("cvode_BadK"));

	case CV_BAD_T:
	    caml_raise_constant(*caml_named_value("cvode_BadT"));

	case CV_BAD_DKY:
	    caml_raise_constant(*caml_named_value("cvode_BadDky"));

	case CV_TOO_CLOSE:
	    caml_raise_constant(*caml_named_value("cvode_TooClose"));

	/* Quadrature */

        case CV_NO_QUAD:
	    caml_raise_constant(*caml_named_value("cvodes_QuadNotInitialized"));

        case CV_QRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_QuadRhsFuncFailure"));

        case CV_FIRST_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstQuadRhsFuncErr"));

        case CV_REPTD_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedQuadRhsFuncErr"));

        case CV_UNREC_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableQuadRhsFuncErr"));

	/* Sensitivity */

        case CV_NO_SENS:
	    caml_raise_constant(*caml_named_value("cvodes_SensNotInitialized"));

        case CV_SRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_SensRhsFuncFailure"));

        case CV_FIRST_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstSensRhsFuncErr"));

        case CV_REPTD_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedSensRhsFuncErr"));

        case CV_UNREC_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableSensRhsFuncErr"));

        case CV_BAD_IS:
	    caml_raise_constant(*caml_named_value("cvodes_BadIS"));

	/* Sensitivity > Quadrature */

        case CV_NO_QUADSENS:
	    caml_raise_constant(*caml_named_value("cvodes_QuadSensNotInitialized"));

        case CV_QSRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_QuadSensRhsFuncFailure"));

        case CV_FIRST_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstQuadSensRhsFuncErr"));

        case CV_REPTD_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedQuadSensRhsFuncErr"));

        case CV_UNREC_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableQuadSensRhsFuncErr"));

	/* Adjoint */

        case CV_NO_ADJ:
	    caml_raise_constant(*caml_named_value("cvodes_AdjointNotInitialized"));

        case CV_NO_FWD:
	    caml_raise_constant(*caml_named_value("cvodes_NoForwardCall"));

        case CV_NO_BCK:
	    caml_raise_constant(*caml_named_value("cvodes_NoBackwardProblem"));

        case CV_BAD_TB0:
	    caml_raise_constant(*caml_named_value("cvodes_BadTB0"));

        case CV_REIFWD_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_ForwardReinitializationFailed"));

        case CV_FWD_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_ForwardFailed"));

        case CV_GETY_BADT:
	    caml_raise_constant(*caml_named_value("cvodes_BadT"));

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

/* basic interface */

