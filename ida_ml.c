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

/* Sundials IDA interface functions that do not involve NVectors. */

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/fail.h>
#include <caml/alloc.h>
#include <caml/callback.h>
#include <caml/bigarray.h>

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <ida/ida_band.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_sptfqmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_lapack.h>
#include <sundials/sundials_config.h>

#include <stdio.h>
#include "ida_ml.h"

#define MAX_ERRMSG_LEN 256

void ida_ml_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == IDA_SUCCESS
	|| flag == IDA_ROOT_RETURN
	|| flag == IDA_TSTOP_RETURN) return;

    switch (flag) {
    case IDA_ILL_INPUT:
	caml_raise_constant(*caml_named_value("ida_IllInput"));
	break;

    case IDA_TOO_MUCH_WORK:
	caml_raise_constant(*caml_named_value("ida_TooMuchWork"));
	break;

    case IDA_TOO_MUCH_ACC:
	caml_raise_constant(*caml_named_value("ida_TooMuchAccuracy"));
	break;

    case IDA_LINIT_FAIL:
	caml_raise_constant(*caml_named_value("ida_LinearInitFailure"));
	break;

    case IDA_LSETUP_FAIL:
	caml_raise_constant(*caml_named_value("ida_LinearSetupFailure"));
	break;

    case IDA_LSOLVE_FAIL:
	caml_raise_constant(*caml_named_value("ida_LinearSolveFailure"));
	break;

    case IDA_RES_FAIL:
	caml_raise_constant(*caml_named_value("ida_ResFuncFailure"));
	break;

    case IDA_REP_RES_ERR:
	caml_raise_constant(*caml_named_value("ida_RepeatedResFuncError"));
	break;

    case IDA_RTFUNC_FAIL:
	caml_raise_constant(*caml_named_value("ida_RootFuncFailure"));
	break;

    case IDA_BAD_K:
	caml_raise_constant(*caml_named_value("ida_BadK"));
	break;

    case IDA_BAD_T:
	caml_raise_constant(*caml_named_value("ida_BadT"));
	break;

    case IDA_BAD_DKY:
	caml_raise_constant(*caml_named_value("ida_BadDky"));
	break;

    default:
	/* e.g. IDA_MEM_NULL, IDA_ILL_INPUT, IDA_MEM_FAIL */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 IDAGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

void ida_ml_set_linear_solver(void *ida_mem, value ls, int n)
{
    int flag;

    if (Is_block(ls)) {
	value arg = Field (ls, 0);

	switch (Tag_val(ls)) {
	case VARIANT_IDA_LINEAR_SOLVER_BAND:
	    flag = IDABand(ida_mem, n,
			   Long_val(Field (arg, 0)),
			   Long_val(Field (arg, 1)));
	    CHECK_FLAG("IDABand", flag);
	    break;

	case VARIANT_IDA_LINEAR_SOLVER_LAPACKBAND:
#if SUNDIALS_BLAS_LAPACK == 1
	    flag = IDALapackBand(ida_mem, n,
				 Long_val(Field (arg, 0)),
				 Long_val(Field (arg, 1)));
	    CHECK_FLAG("IDALapackBand", flag);
#else
	    caml_failwith("Lapack solvers are not available.");
#endif
	    break;

	case VARIANT_IDA_LINEAR_SOLVER_SPGMR:
	    flag = IDASpgmr(ida_mem, Int_val(arg));
	    CHECK_FLAG("IDASpgmr", flag);
	    break;

	case VARIANT_IDA_LINEAR_SOLVER_SPBCG:
	    flag = IDASpbcg(ida_mem, Int_val(arg));
	    CHECK_FLAG("IDASpbcg", flag);
	    break;

	case VARIANT_IDA_LINEAR_SOLVER_SPTFQMR:
	    flag = IDASptfqmr(ida_mem, Int_val(arg));
	    CHECK_FLAG("IDASPtfqmr", flag);
	    break;

	default:
	    caml_failwith("Illegal linear solver block value.");
	    break;
	}

    } else {
	switch (Int_val(ls)) {
	case VARIANT_IDA_LINEAR_SOLVER_DENSE:
	    flag = IDADense(ida_mem, n);
	    CHECK_FLAG("IDADense", flag);
	    break;

	case VARIANT_IDA_LINEAR_SOLVER_LAPACKDENSE:
#if SUNDIALS_BLAS_LAPACK == 1
	    flag = IDALapackDense(ida_mem, n);
	    CHECK_FLAG("IDALapackDense", flag);
#else
	    caml_failwith("Lapack solvers are not available.");
#endif
	    break;

	default:
	    caml_failwith("Illegal linear solver value.");
	    break;
	}
    }
}

CAMLprim void c_ida_session_finalize(value vdata)
{
    if (IDA_MEM_FROM_ML(vdata) != NULL) {
	void *ida_mem = IDA_MEM_FROM_ML(vdata);
	IDAFree(&ida_mem);
    }

    FILE* err_file = (FILE *)Long_val(Field(vdata, RECORD_IDA_SESSION_ERRFILE));
    if (err_file != NULL) {
	fclose(err_file);
    }
}

 
CAMLprim void c_ida_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDASStolerances(IDA_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("IDASStolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_ida_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < IDA_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = IDAGetRootInfo(IDA_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("IDAGetRootInfo", flag);

    CAMLreturn0;
}

CAMLprim value c_ida_get_integrator_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

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

    flag = IDAGetIntegratorStats(IDA_MEM_FROM_ML(vdata),
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
    CHECK_FLAG("IDAGetIntegratorStats", flag);

    r = caml_alloc_tuple(10);
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_RES_EVALS, Val_long(nfevals));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS, Val_long(nlinsetups));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_ERROR_TEST_FAILURES, Val_long(netfails));

    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_LAST_INTERNAL_ORDER, Val_int(qlast));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER, Val_int(qcur));

    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_INITIAL_STEP_SIZE, caml_copy_double(hinused));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_LAST_STEP_SIZE, caml_copy_double(hlast));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_NEXT_STEP_SIZE, caml_copy_double(hcur));
    Store_field(r, RECORD_IDA_INTEGRATOR_STATS_INTERNAL_TIME, caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim value c_ida_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    FILE* err_file = (FILE *)Long_val(Field(vdata, RECORD_IDA_SESSION_ERRFILE));

    if (err_file != NULL) {
	fclose(err_file);
	Store_field(vdata, RECORD_IDA_SESSION_ERRFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    err_file = fopen(String_val(vpath), mode);
    if (err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = IDASetErrFile(IDA_MEM_FROM_ML(vdata), err_file);
    CHECK_FLAG("IDASetErrFile", flag);

    Store_field(vdata, RECORD_IDA_SESSION_ERRFILE, Val_long(err_file));

    CAMLreturn0;
}

CAMLprim value c_ida_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < IDA_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = IDASetRootDirection(IDA_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("IDASetRootDirection", flag);

    CAMLreturn0;
}

/* FIXME: move this to sundials_ml.c? */
value ida_ml_big_real()
{
    CAMLparam0();
    CAMLreturn(caml_copy_double(BIG_REAL));
}

value ida_ml_unit_roundoff()
{
    CAMLparam0();
    CAMLreturn(caml_copy_double(UNIT_ROUNDOFF));
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

CAMLprim value c_ida_get_work_space(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = IDAGetWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrw, &leniw);
    CHECK_FLAG("IDAGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
}


CAMLprim value c_ida_get_num_steps(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumSteps(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_ida_get_num_res_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumResEvals(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumResEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_ida_get_num_lin_solv_setups(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumLinSolvSetups(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_ida_get_num_err_test_fails(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDAGetNumErrTestFails(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_ida_get_last_order(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    int v;

    flag = IDAGetLastOrder(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetLastOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_ida_get_current_order(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    int v;

    flag = IDAGetCurrentOrder(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetCurrentOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_ida_get_actual_init_step(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    realtype v;

    flag = IDAGetActualInitStep(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_ida_get_last_step(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1 (tmp);

    int flag;
    realtype v;

    flag = IDAGetLastStep(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetLastStep", flag);

    tmp = caml_copy_double(v);
    CAMLreturn(tmp);
}

CAMLprim value c_ida_get_current_step(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    realtype v;

    flag = IDAGetCurrentStep(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_ida_get_current_time(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    realtype v;

    flag = IDAGetCurrentTime(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDAGetCurrentTime", flag);

    CAMLreturn(caml_copy_double(v));
}


CAMLprim void c_ida_set_max_ord(value vida_mem, value maxord)
{
    CAMLparam2(vida_mem, maxord);


    int flag = IDASetMaxOrd(IDA_MEM_FROM_ML(vida_mem), Int_val(maxord));
    CHECK_FLAG("IDASetMaxOrd", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_max_num_steps(value vida_mem, value mxsteps)
{
    CAMLparam2(vida_mem, mxsteps);


    int flag = IDASetMaxNumSteps(IDA_MEM_FROM_ML(vida_mem), Long_val(mxsteps));
    CHECK_FLAG("IDASetMaxNumSteps", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_init_step(value vida_mem, value hin)
{
    CAMLparam2(vida_mem, hin);


    int flag = IDASetInitStep(IDA_MEM_FROM_ML(vida_mem), Double_val(hin));
    CHECK_FLAG("IDASetInitStep", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_max_step(value vida_mem, value hmax)
{
    CAMLparam2(vida_mem, hmax);


    int flag = IDASetMaxStep(IDA_MEM_FROM_ML(vida_mem), Double_val(hmax));
    CHECK_FLAG("IDASetMaxStep", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_stop_time(value vida_mem, value tstop)
{
    CAMLparam2(vida_mem, tstop);


    int flag = IDASetStopTime(IDA_MEM_FROM_ML(vida_mem), Double_val(tstop));
    CHECK_FLAG("IDASetStopTime", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_max_err_test_fails(value vida_mem, value maxnef)
{
    CAMLparam2(vida_mem, maxnef);


    int flag = IDASetMaxErrTestFails(IDA_MEM_FROM_ML(vida_mem), Int_val(maxnef));
    CHECK_FLAG("IDASetMaxErrTestFails", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_max_nonlin_iters(value vida_mem, value maxcor)
{
    CAMLparam2(vida_mem, maxcor);


    int flag = IDASetMaxNonlinIters(IDA_MEM_FROM_ML(vida_mem), Int_val(maxcor));
    CHECK_FLAG("IDASetMaxNonlinIters", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_max_conv_fails(value vida_mem, value maxncf)
{
    CAMLparam2(vida_mem, maxncf);


    int flag = IDASetMaxConvFails(IDA_MEM_FROM_ML(vida_mem), Int_val(maxncf));
    CHECK_FLAG("IDASetMaxConvFails", flag);

    CAMLreturn0;
}

CAMLprim void c_ida_set_nonlin_conv_coef(value vida_mem, value nlscoef)
{
    CAMLparam2(vida_mem, nlscoef);


    int flag = IDASetNonlinConvCoef(IDA_MEM_FROM_ML(vida_mem), Double_val(nlscoef));
    CHECK_FLAG("IDASetNonlinConvCoef", flag);

    CAMLreturn0;
}


CAMLprim void c_ida_set_no_inactive_root_warn(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag = IDASetNoInactiveRootWarn(IDA_MEM_FROM_ML(vida_mem));
    CHECK_FLAG("IDASetNoInactiveRootWarn", flag);

    CAMLreturn0;
}

CAMLprim value c_ida_get_tol_scale_factor(value vida_mem)
{
    CAMLparam1(vida_mem);

    realtype r;
    int flag = IDAGetTolScaleFactor(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetTolScaleFactor", flag);

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value c_ida_get_num_nonlin_solv_iters(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDAGetNumNonlinSolvIters(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_get_num_nonlin_solv_conv_fails(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDAGetNumNonlinSolvConvFails(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_get_num_g_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDAGetNumGEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDAGetNumGEvals", flag);

    CAMLreturn(Val_long(r));
}


CAMLprim value c_ida_dls_get_work_space(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = IDADlsGetWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrwLS, &leniwLS);
    CHECK_FLAG("IDADlsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));

    CAMLreturn(r);
}


CAMLprim value c_ida_dls_get_num_jac_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDADlsGetNumJacEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDADlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_dls_get_num_res_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDADlsGetNumResEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDADlsGetNumResEvals", flag);

    CAMLreturn(Val_long(r));
}
