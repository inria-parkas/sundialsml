/***********************************************************************
 *                                                                     *
 *              Ocaml interface to Sundials CVODE solver               *
 *                                                                     *
 *       Timothy Bourke (INRIA Rennes) and Marc Pouzet (LIENS)         *
 *                                                                     *
 *  Copyright 2011 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under the terms of the GNU Library General Public License, with    *
 *  the special exception on linking described in file LICENSE.        *
 *                                                                     *
 ***********************************************************************/

#ifndef _CVODE_ML_H__
#define _CVODE_ML_H__

/* Configuration options */
#define CHECK_MATRIX_ACCESS 1

void cvode_ml_check_flag(const char *call, int flag, void *to_free);

#define CHECK_FLAG(call, flag) if (flag != CV_SUCCESS) \
				 cvode_ml_check_flag(call, flag, NULL)

// TODO: Is there any risk that the Ocaml GC will try to free the
//	 closures? Do we have to somehow record that we're using them,
//	 and then release them again in the free routine?
//	 SEE: cvode_ml_data_alloc and finalize
struct cvode_ml_data {
    void *cvode_mem;
    long int neq;
    intnat num_roots;
    value *closure_rhsfn;
    value *closure_rootsfn;
    value *closure_errh;
    value *closure_errw;

    value *closure_jacfn;
    value *closure_bandjacfn;
    value *closure_presetupfn;
    value *closure_presolvefn;
    value *closure_jactimesfn;

    FILE *err_file;
};
typedef struct cvode_ml_data* cvode_ml_data_p;

value cvode_ml_data_alloc(void* cvode_mem);
void cvode_ml_set_linear_solver(void *cvode_mem, value ls, int n);

#define CVODE_DATA(v) ((cvode_ml_data_p)Data_custom_val(v))
#define CVODE_DATA_FROM_ML(name, v) \
    cvode_ml_data_p (name) = CVODE_DATA(v); \
    if ((name)->cvode_mem == NULL) caml_failwith("This session has been freed")

#define CVODE_MEM_FROM_ML(name, v) \
    void *(name) = (CVODE_DATA(v))->cvode_mem; \
    if ((name) == NULL) caml_failwith("This session has been freed")

value cvode_ml_big_real();

/* Interface with Ocaml types */

#define BIGARRAY_FLOAT (CAML_BA_FLOAT64 | CAML_BA_C_LAYOUT)
#define BIGARRAY_INT (CAML_BA_INT32 | CAML_BA_C_LAYOUT)

#define VARIANT_LMM_ADAMS 0
#define VARIANT_LMM_BDF   1

#define RECORD_BANDRANGE_MUPPER 0
#define RECORD_BANDRANGE_MLOWER 1

#define RECORD_SPRANGE_PRETYPE 0
#define RECORD_SPRANGE_MAXL    1

/* untagged: */
#define VARIANT_LINEAR_SOLVER_DENSE	    0
#define VARIANT_LINEAR_SOLVER_LAPACKDENSE   1
#define VARIANT_LINEAR_SOLVER_DIAG	    2
/* tagged: */
#define VARIANT_LINEAR_SOLVER_BAND		    0
#define VARIANT_LINEAR_SOLVER_LAPACKBAND	    1
#define VARIANT_LINEAR_SOLVER_SPGMR		    2
#define VARIANT_LINEAR_SOLVER_SPBCG		    3
#define VARIANT_LINEAR_SOLVER_SPTFQMR		    4
#define VARIANT_LINEAR_SOLVER_BANDED_SPGMR	    5
#define VARIANT_LINEAR_SOLVER_BANDED_SPBCG	    6
#define VARIANT_LINEAR_SOLVER_BANDED_SPTFQMR	    7

#define VARIANT_SOLVER_RESULT_CONTINUE		0
#define VARIANT_SOLVER_RESULT_ROOTSFOUND	1
#define VARIANT_SOLVER_RESULT_STOPTIMEREACHED	2

#define RECORD_INTEGRATOR_STATS_STEPS			0
#define RECORD_INTEGRATOR_STATS_RHS_EVALS		1
#define RECORD_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS	2
#define RECORD_INTEGRATOR_STATS_ERROR_TEST_FAILURES	3
#define RECORD_INTEGRATOR_STATS_LAST_INTERNAL_ORDER	4
#define RECORD_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER	5
#define RECORD_INTEGRATOR_STATS_INITIAL_STEP_SIZE	6
#define RECORD_INTEGRATOR_STATS_LAST_STEP_SIZE		7
#define RECORD_INTEGRATOR_STATS_NEXT_STEP_SIZE		8
#define RECORD_INTEGRATOR_STATS_INTERNAL_TIME		9

#define RECORD_ERROR_DETAILS_ERROR_CODE	    0
#define RECORD_ERROR_DETAILS_MODULE_NAME    1
#define RECORD_ERROR_DETAILS_FUNCTION_NAME  2
#define RECORD_ERROR_DETAILS_ERROR_MESSAGE  3

#define RECORD_JACOBIAN_ARG_JAC_T	0
#define RECORD_JACOBIAN_ARG_JAC_Y	1
#define RECORD_JACOBIAN_ARG_JAC_FY	2
#define RECORD_JACOBIAN_ARG_JAC_TMP	3

#define RECORD_SPILS_SOLVE_ARG_RHS	0
#define RECORD_SPILS_SOLVE_ARG_GAMMA	1
#define RECORD_SPILS_SOLVE_ARG_DELTA	2
#define RECORD_SPILS_SOLVE_ARG_LEFT	3

#define VARIANT_HANDLER_RHSFN		0
#define VARIANT_HANDLER_ROOTSFN		1
#define VARIANT_HANDLER_ERRORHANDLER	2
#define VARIANT_HANDLER_ERRORWEIGHT	3
#define VARIANT_HANDLER_JACFN		4
#define VARIANT_HANDLER_BANDJACFN	5
#define VARIANT_HANDLER_PRESETUPFN	6
#define VARIANT_HANDLER_PRESOLVEFN	7
#define VARIANT_HANDLER_JACTIMESFN	8

#define VARIANT_PRECOND_TYPE_PRECNONE	0
#define VARIANT_PRECOND_TYPE_PRECLEFT	1
#define VARIANT_PRECOND_TYPE_PRECRIGHT	2
#define VARIANT_PRECOND_TYPE_PRECBOTH	3

#define VARIANT_GRAMSCHMIDT_TYPE_MODIFIEDGS  0
#define VARIANT_GRAMSCHMIDT_TYPE_CLASSICALGS 1

#define VARIANT_PRECONDITIONING_TYPE_PRECNONE	0
#define VARIANT_PRECONDITIONING_TYPE_PRECLEFT	1
#define VARIANT_PRECONDITIONING_TYPE_PRECRIGHT	2
#define VARIANT_PRECONDITIONING_TYPE_PRECBOTH	3

#define VARIANT_GRAMSCHMIDT_TYPE_MODIFIEDGS	0
#define VARIANT_GRAMSCHMIDT_TYPE_CLASSICALGS	1

#endif

