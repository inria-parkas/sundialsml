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

/*
 * This module defines all constants and functions which do not depend on
 * the representation of continuous state vectors, i.e., those that are
 * shared between the bigarray and nvector versions of cvode_ml_nvec.
 *
 */

#ifndef _CVODE_ML_H__
#define _CVODE_ML_H__

/* Configuration options */
#define CHECK_MATRIX_ACCESS 1

void cvode_ml_check_flag(const char *call, int flag, void *to_free);

#define CHECK_FLAG(call, flag) if (flag != CV_SUCCESS) \
				 cvode_ml_check_flag(call, flag, NULL)


/*
 * The session data structure is shared in four parts across the OCaml and C
 * heaps:
 *
 *           C HEAP                 .           OCAML HEAP             |
 * ---------------------------------.-----------------------------------
 *                                  .
 *                                  .    (Program using Sundials/ML)
 *                                  .                |
 *                                  .               \|/
 *                                  .          +------------+
 *                                  .          |  session   | (Tag = Abstract)
 *                                  .          +------------+
 *             +-------------------------------+ cvode      |
 *             |                +--------------+ user data  |
 *            \|/               |   .          | err_file   |
 *       +------------+         |   .          +------------+
 *       | cvode_mem  |         |   .
 *       +------------+         |   .
 *       |    ...     |         |   .
 *    ---+cv_user_data|         |   .
 *    |  |    ...     |         |   .
 *    |  +------------+         |   .      _____________
 *    |_______     _____________+   .     /            |
 *            |   |                 .    /            \|/
 *           \|/ \|/                .   /     +----------------+
 *       +-----------+       _______.__/      |   callbacks    | (Tag = 0)
 *       | user_data |      /       .         +----------------+
 *       +-----------+     /        .         | closure_rhsfn  |
 *       | neqs      |    /         .         | closure_rootfn |
 *       | nroots    |   /          .         |     ...        |
 *       | callbacks +--/           .         +----------------+
 *       +-----------+              .
 *                                  .
 *
 *  Why 4 parts?
 *
 *  * cvode_mem.cv_user_data is write only. The Sundials API provides no way
 *    to read this value outside of the callback functions (we thus maintain
 *    our own link from session).
 *
 *  * user_data is registered with Sundials and then passed to the callback
 *    functions. It must thus be placed in the C Heap (so that it is not
 *    moved by the garbage collector).
 *
 *  * user_data.callbacks is linked to a record in the OCaml heap that
 *    records the various user callback functions for the associated
 *    session. This field is registered as a global root, so that call
 *    backs and the associated closures are not reclaimed; and so that
 *    it is updated after generation copies and compaction.
 *
 *    The closures could be stored directly in user_data, but then each
 *    would have to be registered as a global root.
 *
 *  The init command:
 *  1. creates a session in the OCaml heap.
 *  2. calls CVodeInit to create a cvode_mem in the C heap.
 *  3. creates a callback in the OCaml heap.
 *  4. calls caml_stat_alloc to create a user_data in the C heap.
 *  5. registers user_data.callbacks as a global root.
 *
 *  A finalizer is associated with session, on reclamation it:
 *  1. calls CVodeFree to free cvode_mem.
 *  2. unregisters user_data.callbacks as a global root.
 *  3. calls caml_stat_free to free user_data.
 *  4. closes session.err_file if necessary.
 */

struct cvode_ml_user_data {
    long int neq;
    intnat num_roots;
    value callbacks;
};

typedef struct cvode_ml_user_data* cvode_ml_user_data_p;

#define USER_DATA_NUM_ROOTS(ud) ((cvode_ml_user_data_p)(ud))->num_roots
#define USER_DATA_NEQ(ud) ((cvode_ml_user_data_p)(ud))->neq

#define USER_DATA_RHSFN(ud)   \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 0))

#define USER_DATA_ROOTSFN(ud) \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 1))

#define USER_DATA_ERRH(ud)    \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 2))

#define USER_DATA_ERRW(ud)    \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 3))

#define USER_DATA_JACFN(ud)      \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 4))

#define USER_DATA_BANDJACFN(ud)  \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 5))

#define USER_DATA_PRESETUPFN(ud) \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 6))

#define USER_DATA_PRESOLVEFN(ud) \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 7))

#define USER_DATA_JACTIMESFN(ud) \
    (Field(((cvode_ml_user_data_p)(ud))->callbacks, 8))


#define USER_DATA_SET_RHSFN(ud, f)      \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 0, (f))

#define USER_DATA_SET_ROOTSFN(ud, f)    \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 1, (f))

#define USER_DATA_SET_ERRH(ud, f)       \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 2, (f))

#define USER_DATA_SET_ERRW(ud, f)       \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 3, (f))

#define USER_DATA_SET_JACFN(ud, f)      \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 4, (f))

#define USER_DATA_SET_BANDJACFN(ud, f)  \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 5, (f))

#define USER_DATA_SET_PRESETUPFN(ud, f) \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 6, (f))

#define USER_DATA_SET_PRESOLVEFN(ud, f) \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 7, (f))

#define USER_DATA_SET_JACTIMESFN(ud, f) \
    Store_field(((cvode_ml_user_data_p)(ud))->callbacks, 8, (f))

#define USER_DATA_NUMCALLBACKS 9

struct cvode_ml_session {
    void *cvode_mem;
    cvode_ml_user_data_p user_data;
    FILE *err_file;
};
typedef struct cvode_ml_session* cvode_ml_session_p;

#define CVODE_SESSION(v) ((cvode_ml_session_p)Data_custom_val(v))
#define CVODE_SESSION_FROM_ML(name, v) \
    cvode_ml_session_p (name) = CVODE_SESSION(v); \
    if ((name)->cvode_mem == NULL) caml_failwith("This session has been freed")

#define CVODE_MEM_FROM_ML(name, v) \
    void *(name) = (CVODE_SESSION(v))->cvode_mem; \
    if ((name) == NULL) caml_failwith("This session has been freed")

value cvode_ml_session_alloc(void* cvode_mem);
void cvode_ml_set_linear_solver(void *cvode_mem, value ls, int n);

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

