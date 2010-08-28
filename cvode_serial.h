/* Aug 2010, Timothy Bourke (INRIA)
 *
 * Ocaml interface to the Sundials 2.4.0 CVode solver for serial NVectors.
 *
 */

#ifndef __CVODE_SERIAL_H__

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

/* linear solvers */
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>

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
#define VARIANT_LINEAR_SOLVER_BAND	    0
#define VARIANT_LINEAR_SOLVER_LAPACKBAND    1
#define VARIANT_LINEAR_SOLVER_SPGMR	    2
#define VARIANT_LINEAR_SOLVER_SPBCG	    3
#define VARIANT_LINEAR_SOLVER_SPTFQMR	    4

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
#define VARIANT_HANDLER_JACFN		3
#define VARIANT_HANDLER_BANDJACFN	4
#define VARIANT_HANDLER_PRESETUPFN	5
#define VARIANT_HANDLER_PRESOLVEFN	6
#define VARIANT_HANDLER_JACTIMESFN	7

#define VARIANT_HANDLER_JACFN		3
#define VARIANT_HANDLER_BANDJACFN	4
#define VARIANT_HANDLER_PRESETUPFN	5
#define VARIANT_HANDLER_PRESOLVEFN	6
#define VARIANT_HANDLER_JACTIMESFN	7

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

/*
 * main functions
 */
 
CAMLprim value c_init(value lmm, value iter, value initial, value num_roots,
		      value t0);
CAMLprim value c_neqs(value vdata);
CAMLprim value c_nroots(value vdata);
CAMLprim value c_set_tolerances(value vdata, value reltol, value abstol);
CAMLprim value c_reinit(value vdata, value t0, value y0);
CAMLprim value c_get_roots(value vdata, value roots);
CAMLprim value c_free(value vdata);
CAMLprim value c_advance(value vdata, value nextt, value y);
CAMLprim value c_step(value vdata, value nextt, value y);
CAMLprim value c_get_dky(value vdata, value vt, value vk, value vy);
CAMLprim value c_integrator_stats(value vdata);
CAMLprim value c_last_step_size(value vdata);
CAMLprim value c_next_step_size(value vdata);

/*
 * optional input functions
 */

CAMLprim value c_set_error_file(value vdata, value vpath, value vtrunc);
CAMLprim value c_register_handler(value vdata, value handler);
/* call c_register_handler first: */
CAMLprim value c_enable_error_handler(value vdata);
CAMLprim value c_set_max_ord(value vdata, value maxord);
CAMLprim value c_set_max_num_steps(value vdata, value mxsteps);
CAMLprim value c_set_max_hnil_warns(value vdata, value mxhnil);
CAMLprim value c_set_stability_limit_detection(value vdata, value stldet);
CAMLprim value c_set_initial_step_size(value vdata, value hin);
CAMLprim value c_set_min_abs_step_size(value vdata, value hmin);
CAMLprim value c_set_max_abs_step_size(value vdata, value hmax);
CAMLprim value c_set_stop_time(value vdata, value tstop);
CAMLprim value c_set_max_error_test_failures(value vdata, value maxnef);
CAMLprim value c_set_max_nonlinear_iterations(value vdata, value maxcor);
CAMLprim value c_set_max_convergence_failures(value vdata, value maxncf);
CAMLprim value c_set_nonlinear_convergence_coeffficient(value vdata, value nlscoef);
CAMLprim value c_set_nonlinear_iteration_type(value vdata, value iter);
CAMLprim value c_set_root_direction(value vdata, value rootdirs);
CAMLprim value c_disable_inactive_root_warnings(value vdata);

/*
 * direct linear solvers optional input functions
 */

/* call c_register_handler first: */
CAMLprim value c_enable_dense_jacobian_fn(value vdata);
/* call c_register_handler first: */
CAMLprim value c_enable_band_jacobian_fn(value vdata);

/* iterative linear solvers optional input functions */

/* call c_register_handler for both functions first: */
CAMLprim value c_enable_preconditioner_fns(value vdata);
/* call c_register_handler first: */
CAMLprim value c_enable_jacobian_times_vector_fn(value vdata);
CAMLprim value c_set_preconditioning_type(value vdata, value vptype);
CAMLprim value c_set_gramschmidt_orthogonalization(value vdata, value vgstype);
CAMLprim value c_set_eps_linear_convergence_factor(value vdata, value eplifac);
CAMLprim value c_set_max_subspace_dimension(value vdata, value maxl);

/*
 * functions for the abstract types Densematrix.t and Bandmatrix.t
 */

CAMLprim value c_densematrix_get(value vmatrix, value vij);
CAMLprim value c_densematrix_set(value vmatrix, value vij, value v);
CAMLprim value c_bandmatrix_get(value vmatrix, value vij);
CAMLprim value c_bandmatrix_set(value vmatrix, value vij, value v);

/*
 * Internals
 */
void ml_cvode_check_flag(const char *call, int flag, void *to_free);
typedef struct ml_cvode_data* ml_cvode_data_p;
value ml_cvode_data_alloc(mlsize_t approx_size);
void *ml_cvode_mem(value vdata);

#endif

