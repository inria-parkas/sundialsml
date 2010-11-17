/* Aug 2010, Timothy Bourke (INRIA)
 *
 * Trickier definitions for Sundials interface that do not involve NVectors.
 *
 */

/*
 * TODO:
 * - call Gc.full_major () from the f and roots routines to see if we get
 *   any segmentation fault problems.
 * - see notes throughout program related to garbage collection.
 */

#include <cvode/cvode.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>
#include <caml/alloc.h>

/* linear solvers */
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#include <cvode/cvode_bandpre.h>

#include "ml_cvode.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

static const char *callback_ocaml_names[] = {
    "cvode_serial_callback_rhsfn",
    "cvode_serial_callback_rootsfn",
    "cvode_serial_callback_errorhandler",
    "cvode_serial_callback_errorweight",
    "cvode_serial_callback_jacfn",
    "cvode_serial_callback_bandjacfn",
    "cvode_serial_callback_presetupfn",
    "cvode_serial_callback_presolvefn",
    "cvode_serial_callback_jactimesfn",
};

static void finalize_closure(value** closure_field)
{
    if (*closure_field != NULL) {
	caml_remove_generational_global_root(*closure_field);
	*closure_field = NULL;
    }
}

static void finalize(value vdata)
{
    CVODE_DATA_FROM_ML(data, vdata);

    // TODO:
    // The Ocaml Manual (18.9.1) says:
    // ``Note: the finalize, compare, hash, serialize and deserialize
    // functions attached to custom block descriptors must never trigger a
    // garbage collection. Within these functions, do not call any of the
    // Caml allocation functions, and do not perform a callback into Caml
    // code. Do not use CAMLparam to register the parameters to these
    // functions, and do not use CAMLreturn to return the result.''
    //
    // But, obviously, we're calling two caml functions. We need to find out
    // if this is ok.
    finalize_closure(&(data->closure_rhsfn));
    finalize_closure(&(data->closure_rootsfn));
    finalize_closure(&(data->closure_errh));
    finalize_closure(&(data->closure_errw));

    finalize_closure(&(data->closure_jacfn));
    finalize_closure(&(data->closure_bandjacfn));
    finalize_closure(&(data->closure_presetupfn));
    finalize_closure(&(data->closure_presolvefn));
    finalize_closure(&(data->closure_jactimesfn));

    if (data->cvode_mem != NULL) {
	CVodeFree(&(data->cvode_mem));
    }

    if (data->err_file != NULL) {
	fclose(data->err_file);
	data->err_file = NULL;
    }
}

// TODO:
// The Ocaml Manual (18.9.3) says:
// ``The contents of custom blocks are not scanned by the garbage collector,
// and must therefore not contain any pointer inside the Caml heap. In other
// terms, never store a Caml value in a custom block, and do not use Field,
// Store_field nor caml_modify to access the data part of a custom block.
// Conversely, any C data structure (not containing heap pointers) can be
// stored in a custom block.''
//
// But, obviously, we're storing two closure values in the struct. We need
// to find out if and when this is ok.
//
value ml_cvode_data_alloc(void* cvode_mem)
{
    CAMLparam0();
    CAMLlocal1(r);

    mlsize_t approx_size = 0;
    long int lenrw, leniw;
    if (CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw) == CV_SUCCESS) {
    	approx_size = lenrw * sizeof(realtype) + leniw * sizeof(long int);
    }

    r = caml_alloc_final(sizeof(struct ml_cvode_data),
			 &finalize, approx_size, approx_size * 5);

    ml_cvode_data_p d = CVODE_DATA(r);
    d->cvode_mem = cvode_mem;

    d->err_file = NULL;
    d->closure_rhsfn = NULL;
    d->closure_rootsfn = NULL;
    d->closure_errh = NULL;
    d->closure_errw = NULL;
    d->closure_jacfn = NULL;
    d->closure_bandjacfn = NULL;
    d->closure_presetupfn = NULL;
    d->closure_presolvefn = NULL;
    d->closure_jactimesfn = NULL;
    
    CAMLreturn(r);
}

void ml_cvode_check_flag(const char *call, int flag, void *to_free)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	|| flag == CV_ROOT_RETURN
	|| flag == CV_TSTOP_RETURN) return;

    if (to_free != NULL) free(to_free);

    switch (flag) {
    case CV_ILL_INPUT:
	caml_raise_constant(*caml_named_value("cvode_IllInput"));
	break;

    case CV_TOO_CLOSE:
	caml_raise_constant(*caml_named_value("cvode_TooClose"));
	break;

    case CV_TOO_MUCH_WORK:
	caml_raise_constant(*caml_named_value("cvode_TooMuchWork"));
	break;

    case CV_TOO_MUCH_ACC:
	caml_raise_constant(*caml_named_value("cvode_TooMuchAccuracy"));
	break;

    case CV_ERR_FAILURE:
	caml_raise_constant(*caml_named_value("cvode_ErrFailure"));
	break;

    case CV_CONV_FAILURE:
	caml_raise_constant(*caml_named_value("cvode_ConvergenceFailure"));
	break;

    case CV_LINIT_FAIL:
	caml_raise_constant(*caml_named_value("cvode_LinearInitFailure"));
	break;

    case CV_LSETUP_FAIL:
	caml_raise_constant(*caml_named_value("cvode_LinearSetupFailure"));
	break;

    case CV_LSOLVE_FAIL:
	caml_raise_constant(*caml_named_value("cvode_LinearSolveFailure"));
	break;

    case CV_RHSFUNC_FAIL:
	caml_raise_constant(*caml_named_value("cvode_RhsFuncFailure"));
	break;

    case CV_FIRST_RHSFUNC_ERR:
	caml_raise_constant(*caml_named_value("cvode_FirstRhsFuncError"));
	break;

    case CV_REPTD_RHSFUNC_ERR:
	caml_raise_constant(*caml_named_value("cvode_RepeatedRhsFuncError"));
	break;

    case CV_UNREC_RHSFUNC_ERR:
	caml_raise_constant(*caml_named_value("cvode_UnrecoverableRhsFuncError"));
	break;

    case CV_RTFUNC_FAIL:
	caml_raise_constant(*caml_named_value("cvode_RootFuncFailure"));
	break;

    case CV_BAD_K:
	caml_raise_constant(*caml_named_value("cvode_BadK"));
	break;

    case CV_BAD_T:
	caml_raise_constant(*caml_named_value("cvode_BadT"));
	break;

    case CV_BAD_DKY:
	caml_raise_constant(*caml_named_value("cvode_BadDky"));
	break;

    default:
	/* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 CVodeGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

static int precond_type(value vptype)
{
    CAMLparam1(vptype);

    int ptype;
    switch (Int_val(vptype)) {
    case VARIANT_PRECONDITIONING_TYPE_PRECNONE:
	ptype = PREC_NONE;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECLEFT:
	ptype = PREC_LEFT;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECRIGHT:
	ptype = PREC_RIGHT;
	break;

    case VARIANT_PRECONDITIONING_TYPE_PRECBOTH:
	ptype = PREC_BOTH;
	break;
    }

    CAMLreturn(ptype);
}

/* callbacks */

static void errh(
	int error_code,
	const char *module,
	const char *func,
	char *msg,
	void *eh_data)
{
    CAMLparam0();
    CAMLlocal1(a);

    value *closure_errh = ((ml_cvode_data_p)eh_data)->closure_errh;

    a = caml_alloc_tuple(4);
    Store_field(a, RECORD_ERROR_DETAILS_ERROR_CODE, Val_int(error_code));
    Store_field(a, RECORD_ERROR_DETAILS_MODULE_NAME, caml_copy_string(module));
    Store_field(a, RECORD_ERROR_DETAILS_FUNCTION_NAME, caml_copy_string(func));
    Store_field(a, RECORD_ERROR_DETAILS_ERROR_MESSAGE, caml_copy_string(msg));

    caml_callback(*closure_errh, a);

    CAMLreturn0;
}

CAMLprim value c_enable_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
 
    int flag = CVodeSetErrHandlerFn(data->cvode_mem, errh, (void *)data);
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn0;
}

/* basic interface */

void set_linear_solver(void *cvode_mem, value ls, int n)
{
    int flag;

    if (Is_block(ls)) {
	int field0 = Field(Field(ls, 0), 0); /* mupper, pretype */
	int field1 = Field(Field(ls, 0), 1); /* mlower, maxl */
	value sprange, bandrange;

	switch (Tag_val(ls)) {
	case VARIANT_LINEAR_SOLVER_BAND:
	    flag = CVBand(cvode_mem, n, Int_val(field0), Int_val(field1));
	    CHECK_FLAG("CVBand", flag);
	    break;

#if SUNDIALS_BLAS_LAPACK == 1
	case VARIANT_LINEAR_SOLVER_LAPACKBAND:
	    flag = CVLapackBand(cvode_mem, n, Int_val(field0), Int_val(field1));
	    CHECK_FLAG("CVLapackBand", flag);
	    break;
#endif

	case VARIANT_LINEAR_SOLVER_SPGMR:
	    flag = CVSpgmr(cvode_mem, precond_type(field0), Int_val(field1));
	    CHECK_FLAG("CVSpgmr", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_SPBCG:
	    flag = CVSpbcg(cvode_mem, precond_type(field0), Int_val(field1));
	    CHECK_FLAG("CVSpbcg", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_SPTFQMR:
	    flag = CVSptfqmr(cvode_mem, precond_type(field0), Int_val(field1));
	    CHECK_FLAG("CVSPtfqmr", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_BANDED_SPGMR:
	    sprange = Field(ls, 0);
	    bandrange = Field(ls, 1);

	    flag = CVSpgmr(cvode_mem, precond_type(Field(sprange, 0)),
				      Int_val(Field(sprange, 1)));
	    CHECK_FLAG("CVSpgmr", flag);

	    flag = CVBandPrecInit(cvode_mem, n, Int_val(Field(bandrange, 0)),
						Int_val(Field(bandrange, 1)));
	    CHECK_FLAG("CVBandPrecInit", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_BANDED_SPBCG:
	    sprange = Field(ls, 0);
	    bandrange = Field(ls, 1);

	    flag = CVSpbcg(cvode_mem, precond_type(Field(sprange, 0)),
				      Int_val(Field(sprange, 1)));
	    CHECK_FLAG("CVSpbcg", flag);

	    flag = CVBandPrecInit(cvode_mem, n, Int_val(Field(bandrange, 0)),
						Int_val(Field(bandrange, 1)));
	    CHECK_FLAG("CVBandPrecInit", flag);
	    break;

	case VARIANT_LINEAR_SOLVER_BANDED_SPTFQMR:
	    sprange = Field(ls, 0);
	    bandrange = Field(ls, 1);

	    flag = CVSptfqmr(cvode_mem, precond_type(Field(sprange, 0)),
				        Int_val(Field(sprange, 1)));
	    CHECK_FLAG("CVSptfqmr", flag);

	    flag = CVBandPrecInit(cvode_mem, n, Int_val(Field(bandrange, 0)),
						Int_val(Field(bandrange, 1)));
	    CHECK_FLAG("CVBandPrecInit", flag);
	    break;

	default:
	    caml_failwith("Illegal linear solver block value.");
	    break;
	}

    } else {
	switch (Int_val(ls)) {
	case VARIANT_LINEAR_SOLVER_DENSE:
	    flag = CVDense(cvode_mem, n);
	    CHECK_FLAG("CVDense", flag);
	    break;

#if SUNDIALS_BLAS_LAPACK == 1
	case VARIANT_LINEAR_SOLVER_LAPACKDENSE:
	    flag = CVLapackDense(cvode_mem, n);
	    CHECK_FLAG("CVLapackDense", flag);
	    break;
#endif

	case VARIANT_LINEAR_SOLVER_DIAG:
	    flag = CVDiag(cvode_mem);
	    CHECK_FLAG("CVDiag", flag);
	    break;

	default:
	    caml_failwith("Illegal linear solver value.");
	    break;
	}
    }
}
 
CAMLprim value c_neqs(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    CAMLreturn(Val_int(data->neq));
}

CAMLprim value c_nroots(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    CAMLreturn(Val_int(data->num_roots));
}

CAMLprim value c_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    CVODE_DATA_FROM_ML(data, vdata);

    int flag = CVodeSStolerances(data->cvode_mem,
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("CVodeSStolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    CVODE_DATA_FROM_ML(data, vdata);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = Caml_ba_data_val(roots);

    if (roots_l < data->num_roots) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = CVodeGetRootInfo(data->cvode_mem, roots_d);
    CHECK_FLAG("CVodeGetRootInfo", flag);

    CAMLreturn0;
}

CAMLprim value c_free(value vdata)
{
    CAMLparam1(vdata);
    finalize(vdata);
    Store_field(vdata, 1, (value)NULL);
    CAMLreturn0;
}

CAMLprim value c_get_integrator_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    CVODE_DATA_FROM_ML(data, vdata);
    int flag;

    long int nsteps;
    long int nfevals;    
    long int nlinsetups;
    long int netfails;

    int qlast;
    int qcur;	 

    realtype hinused;
    realtype hlast;
    realtype hcur;
    realtype tcur;

    flag = CVodeGetIntegratorStats(data->cvode_mem,
	&nsteps,
	&nfevals,    
	&nlinsetups,
	&netfails,
	&qlast,
	&qcur,	 
	&hinused,
	&hlast,
	&hcur,
	&tcur
    ); 
    CHECK_FLAG("CVodeGetIntegratorStats", flag);

    r = caml_alloc_tuple(10);
    Store_field(r, RECORD_INTEGRATOR_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_INTEGRATOR_STATS_RHS_EVALS, Val_long(nfevals));
    Store_field(r, RECORD_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS, Val_long(nlinsetups));
    Store_field(r, RECORD_INTEGRATOR_STATS_ERROR_TEST_FAILURES, Val_long(netfails));

    Store_field(r, RECORD_INTEGRATOR_STATS_LAST_INTERNAL_ORDER, Val_int(qlast));
    Store_field(r, RECORD_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER, Val_int(qcur));

    Store_field(r, RECORD_INTEGRATOR_STATS_INITIAL_STEP_SIZE, caml_copy_double(hinused));
    Store_field(r, RECORD_INTEGRATOR_STATS_LAST_STEP_SIZE, caml_copy_double(hlast));
    Store_field(r, RECORD_INTEGRATOR_STATS_NEXT_STEP_SIZE, caml_copy_double(hcur));
    Store_field(r, RECORD_INTEGRATOR_STATS_INTERNAL_TIME, caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim value c_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    CVODE_DATA_FROM_ML(data, vdata);

    if (data->err_file != NULL) {
	fclose(data->err_file);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    data->err_file = fopen(String_val(vpath), mode);
    if (data->err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = CVodeSetErrFile(data->cvode_mem, data->err_file);
    CHECK_FLAG("CVodeSetErrFile", flag);

    CAMLreturn0;
}

CAMLprim value c_register_handler(value vdata, value handler)
{
    CAMLparam2(vdata, handler);

    CVODE_DATA_FROM_ML(data, vdata);
    const char* ocaml_name = callback_ocaml_names[Int_val(handler)];
    value **handler_field;

    switch (Int_val(handler)) {
    case VARIANT_HANDLER_RHSFN:
	handler_field = &(data->closure_rhsfn);
	break;

    case VARIANT_HANDLER_ROOTSFN:
	handler_field = &(data->closure_rootsfn);
	break;

    case VARIANT_HANDLER_ERRORHANDLER:
	handler_field = &(data->closure_errh);
	break;

    case VARIANT_HANDLER_ERRORWEIGHT:
	handler_field = &(data->closure_errw);
	break;

    case VARIANT_HANDLER_JACFN:
	handler_field = &(data->closure_jacfn);
	break;

    case VARIANT_HANDLER_BANDJACFN:
	handler_field = &(data->closure_bandjacfn);
	break;

    case VARIANT_HANDLER_PRESETUPFN:
	handler_field = &(data->closure_presetupfn);
	break;

    case VARIANT_HANDLER_PRESOLVEFN:
	handler_field = &(data->closure_presolvefn);
	break;

    case VARIANT_HANDLER_JACTIMESFN:
	handler_field = &(data->closure_jactimesfn);
	break;

    default:
	break;
    }

    if ((*handler_field) != NULL) {
	caml_remove_generational_global_root(*handler_field);
    }
    (*handler_field) = caml_named_value(ocaml_name);
    // TODO: check if this call is necessary and ok:
    caml_register_generational_global_root(*handler_field);

    CAMLreturn0;
}

CAMLprim value c_set_iter_type(value vdata, value iter)
{
    CAMLparam2(vdata, iter);

    CVODE_DATA_FROM_ML(data, vdata);

    int iter_c;
    if (Is_block(iter)) {
	iter_c = CV_NEWTON;
    } else {
	iter_c = CV_FUNCTIONAL;
    }

    int flag = CVodeSetIterType(data->cvode_mem, iter_c);
    CHECK_FLAG("CVodeSetIterType", flag);

    if (iter_c == CV_NEWTON) {
	set_linear_solver(data->cvode_mem, Field(iter, 0), data->neq);
    }

    CAMLreturn0;
}

CAMLprim value c_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    CVODE_DATA_FROM_ML(data, vdata);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = Caml_ba_data_val(rootdirs);

    if (rootdirs_l < data->num_roots) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = CVodeSetRootDirection(data->cvode_mem, rootdirs_d);
    CHECK_FLAG("CVodeSetRootDirection", flag);

    CAMLreturn0;
}

CAMLprim value c_set_prec_type(value vcvode_mem, value vptype)
{
    CAMLparam2(vcvode_mem, vptype);
    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);

    int flag = CVSpilsSetPrecType(cvode_mem, precond_type(vptype));
    CHECK_FLAG("CVSpilsSetPrecType", flag);

    CAMLreturn0;
}

value ml_cvode_big_real()
{
    CAMLparam0();
    CAMLreturn(caml_copy_double(BIG_REAL));
}

