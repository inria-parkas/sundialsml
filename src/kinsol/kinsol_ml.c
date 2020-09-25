/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include "../config.h"

#include <kinsol/kinsol.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

/* linear solvers */
#if   400 <= SUNDIALS_LIB_VERSION
#include <kinsol/kinsol_ls.h>
#else
#include <kinsol/kinsol_direct.h>
#include <kinsol/kinsol_spils.h>
#endif

#if SUNDIALS_LIB_VERSION < 410
#include <kinsol/kinsol_impl.h>
#endif

#if SUNDIALS_LIB_VERSION < 300
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#if 260 <= SUNDIALS_LIB_VERSION
#include <kinsol/kinsol_spfgmr.h>
#endif
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#endif

#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
#include <kinsol/kinsol_lapack.h>
#endif

#include "../lsolvers/sundials_linearsolver_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"
#include "kinsol_ml.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

CAMLprim value sunml_kinsol_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (KINSOL, exns);
    CAMLreturn (Val_unit);
}

int sunml_kinsol_translate_exception(value session, value r,
			       recoverability recoverable)
{
    CAMLparam2(session, r);
    CAMLlocal1(exn);

    if (!Is_exception_result(r)) CAMLreturnT (int, 0);

    r = Extract_exception(r);

    if (recoverable && Field(r, 0) == SUNDIALS_EXN_TAG (RecoverableFailure))
	CAMLreturnT (int, 1);

    /* Unrecoverable error.  Save the exception and return -1.  */
    exn = caml_alloc_small (1,0);
    Field (exn, 0) = r;
    Store_field (session, RECORD_KINSOL_SESSION_EXN_TEMP, exn);
    CAMLreturnT (int, -1);
}


static void errh(
	int error_code,
	const char *module,
	const char *func,
	char *msg,
	void *eh_data)
{
    CAMLparam0();
    CAMLlocal2(session, a);
    value *backref = eh_data;

    a = caml_alloc_tuple(RECORD_SUNDIALS_ERROR_DETAILS_SIZE);
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_CODE,
                Val_int(error_code));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_MODULE_NAME,
                caml_copy_string(module));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_FUNCTION_NAME,
                caml_copy_string(func));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_MESSAGE,
                caml_copy_string(msg));

    WEAK_DEREF (session, *backref);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (Field(session, RECORD_KINSOL_SESSION_ERRH), a);
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined error handler");

    CAMLreturn0;
}

CAMLprim value sunml_kinsol_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = KINSetErrHandlerFn(KINSOL_MEM_FROM_ML(vdata), errh,
				  KINSOL_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("KINSetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = KINSetErrHandlerFn(KINSOL_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("KINSetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

static void infoh(
	const char *module,
	const char *func,
	char *msg,
	void *ih_data)
{
    CAMLparam0();
    CAMLlocal2(session, a);
    value *backref = ih_data;

    a = caml_alloc_tuple(RECORD_SUNDIALS_ERROR_DETAILS_SIZE);
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_CODE, Val_int(0));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_MODULE_NAME,
                caml_copy_string(module));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_FUNCTION_NAME,
                caml_copy_string(func));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_MESSAGE,
                caml_copy_string(msg));

    WEAK_DEREF (session, *backref);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (Field (session, RECORD_KINSOL_SESSION_INFOH),
				 a);
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined info handler");

    CAMLreturn0;
}

CAMLprim value sunml_kinsol_set_info_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = KINSetInfoHandlerFn(KINSOL_MEM_FROM_ML(vdata), infoh,
				   KINSOL_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("KINSetInfoHandlerFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_clear_info_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = KINSetInfoHandlerFn(KINSOL_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("KINSetInfoHandlerFn", flag);

    CAMLreturn (Val_unit);
}

static int sysfn(N_Vector uu, N_Vector val, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(vuu, vval, session);

    vuu = NVEC_BACKLINK(uu);
    vval = NVEC_BACKLINK(val);

    WEAK_DEREF (session, *(value*)user_data);

    // The data payloads inside vuu and vval are only valid during this
    // call, afterward that memory goes back to kinsol. These bigarrays must
    // not be retained by closure_rhsfn! If it wants a permanent copy, then
    // it has to make it manually.

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(KINSOL_SYSFN_FROM_ML (session), vuu, vval);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static value make_prec_solve_arg(N_Vector uscale, N_Vector fscale)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_SIZE);
    Store_field(r, RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_USCALE,
	        NVEC_BACKLINK(uscale));
    Store_field(r, RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_FSCALE,
		NVEC_BACKLINK(fscale));

    CAMLreturn(r);
}

value sunml_kinsol_make_jac_arg(N_Vector u, N_Vector fu, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_KINSOL_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_U, NVEC_BACKLINK(u));
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_FU, NVEC_BACKLINK(fu));
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

value sunml_kinsol_make_double_tmp(N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, NVEC_BACKLINK(tmp1));
    Store_field(r, 1, NVEC_BACKLINK(tmp2));
    CAMLreturn(r);
}

#if 300 <= SUNDIALS_LIB_VERSION
static int jacfn(
	N_Vector u,
	N_Vector fu,	     
	SUNMatrix Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocalN (args, 2);
    CAMLlocal2(session, cb);

    WEAK_DEREF (session, *(value*)user_data);
    cb = KINSOL_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    args[0] = sunml_kinsol_make_jac_arg(u, fu, sunml_kinsol_make_double_tmp(tmp1, tmp2));
    args[1] = MAT_BACKLINK(Jac);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, UNRECOVERABLE));
}

#else
/* Dense and band Jacobians only work with serial NVectors.  */
static int jacfn(
	long int n,
	N_Vector u,
	N_Vector fu,	     
	DlsMat Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocalN (args, 2);
    CAMLlocal3(session, cb, dmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = KINSOL_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    dmat = Field(cb, 1);
    if (dmat == Val_none) {
	Store_some(dmat, sunml_matrix_dense_wrap(Jac));
	Store_field(cb, 1, dmat);
    }

    args[0] = sunml_kinsol_make_jac_arg(u, fu, sunml_kinsol_make_double_tmp(tmp1, tmp2));
    args[1] = Some_val(dmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, UNRECOVERABLE));
}

static int bandjacfn(
	long int N,
	long int mupper,
	long int mlower, 	 
	N_Vector u,
	N_Vector fu, 	 
	DlsMat Jac,
	void *user_data, 	 
	N_Vector tmp1,
	N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal3(session, cb, bmat);

    WEAK_DEREF (session, *(value*)user_data);
    cb = KINSOL_LS_CALLBACKS_FROM_ML(session);
    cb = Field (cb, 0);

    bmat = Field(cb, 1);
    if (bmat == Val_none) {
	Store_some(bmat, sunml_matrix_band_wrap(Jac));
	Store_field(cb, 1, bmat);
    }

    args[0] = sunml_kinsol_make_jac_arg(u, fu, sunml_kinsol_make_double_tmp(tmp1, tmp2));
    args[1] = Some_val(bmat);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (Field(cb, 0), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, UNRECOVERABLE));
}
#endif

static int precsetupfn(
    N_Vector uu,
    N_Vector uscale,
    N_Vector fu,
    N_Vector fscale,
    void *user_data
#if SUNDIALS_LIB_VERSION < 300
    ,
    N_Vector tmp1,
    N_Vector tmp2
#endif
    )
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 2);

    args[0] = sunml_kinsol_make_jac_arg(uu, fu, Val_unit);
    args[1] = make_prec_solve_arg(uscale, fscale);

    WEAK_DEREF (session, *(value*)user_data);
    cb = KINSOL_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_KINSOL_SPILS_PRECFNS_PREC_SETUP_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int precsolvefn(
	N_Vector uu,
	N_Vector uscale,
	N_Vector fu,
	N_Vector fscale,
	N_Vector vv,
	void *user_data
#if SUNDIALS_LIB_VERSION < 300
	,
	N_Vector tmp
#endif
	)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN(args, 3);

    args[0] = sunml_kinsol_make_jac_arg(uu, fu, Val_unit);
    args[1] = make_prec_solve_arg(uscale, fscale);
    args[2] = NVEC_BACKLINK(vv);

    WEAK_DEREF (session, *(value*)user_data);
    cb = KINSOL_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_KINSOL_SPILS_PRECFNS_PREC_SOLVE_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(cb, 3, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int jactimesfn(
	N_Vector v,
	N_Vector Jv,
	N_Vector u,
	booleantype *new_uu,
	void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, cb);
    CAMLlocalN (args, 4);

    args[0] = NVEC_BACKLINK(v);
    args[1] = NVEC_BACKLINK(Jv);
    args[2] = NVEC_BACKLINK(u);
    args[3] = Val_bool(*new_uu);

    WEAK_DEREF (session, *(value*)user_data);
    cb = KINSOL_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    if (!Is_exception_result (r)) {
	*new_uu = Bool_val (r);
	CAMLreturnT (int, 0);
    }

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value sunml_kinsol_dls_dense (value vkin_mem, value vset_jac)
{
    CAMLparam2(vkin_mem, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINDense (kin_mem, neqs);
    CHECK_FLAG ("KINDense", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetDenseJacFn(KINSOL_MEM_FROM_ML(vkin_mem), jacfn);
	CHECK_DLS_FLAG("KINDlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_dls_lapack_dense (value vkin_mem, value vset_jac)
{
    CAMLparam2 (vkin_mem, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINLapackDense (kin_mem, neqs);
    CHECK_FLAG ("KINLapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetDenseJacFn (KINSOL_MEM_FROM_ML (vkin_mem), jacfn);
	CHECK_DLS_FLAG("KINDlsSetDenseJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_dls_band (value vkin_mem,
				      value vmupper,
				      value vmlower,
				      value vset_jac)
{
    CAMLparam4(vkin_mem, vmupper, vmlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINBand (kin_mem, neqs, Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("KINBand", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vkin_mem), bandjacfn);
	CHECK_DLS_FLAG("KINDlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_dls_lapack_band (value vkin_mem, value vmupper,
					     value vmlower, value vset_jac)
{
    CAMLparam4(vkin_mem, vmupper, vmlower, vset_jac);
#if SUNDIALS_LIB_VERSION < 300 && defined SUNDIALS_ML_LAPACK
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINLapackBand (kin_mem, neqs,
			  Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("KINLapackBand", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vkin_mem), bandjacfn);
	CHECK_DLS_FLAG("KINDlsSetBandJacFn", flag);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_linear_solver (value vkin_mem, value vlsolv,
					       value vojmat, value vhasjac)
{
    CAMLparam4(vkin_mem, vlsolv, vojmat, vhasjac);
#if 400 <= SUNDIALS_LIB_VERSION
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = (vojmat == Val_none) ? NULL : MAT_VAL(Some_val(vojmat));
    int flag;

    flag = KINSetLinearSolver(kin_mem, lsolv, jmat);
    CHECK_LS_FLAG ("KINSetLinearSolver", flag);
    flag = KINSetJacFn(kin_mem, Bool_val(vhasjac) ? jacfn : NULL);
    CHECK_LS_FLAG("KINSetJacFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_dls_set_linear_solver (value vkin_mem, value vlsolv,
						   value vjmat, value vhasjac)
{
    CAMLparam4(vkin_mem, vlsolv, vjmat, vhasjac);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    SUNMatrix jmat = MAT_VAL(vjmat);
    int flag;

    flag = KINDlsSetLinearSolver(kin_mem, lsolv, jmat);
    CHECK_DLS_FLAG ("KINDlsSetLinearSolver", flag);
    flag = KINDlsSetJacFn(kin_mem, Bool_val(vhasjac) ? jacfn : NULL);
    CHECK_DLS_FLAG("KINDlsSetJacFn", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_spils_set_linear_solver (value vkin_mem, value vlsolv)
{
    CAMLparam2(vkin_mem, vlsolv);
#if 300 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 400
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    SUNLinearSolver lsolv = LSOLVER_VAL(vlsolv);
    int flag;

    flag = KINSpilsSetLinearSolver(kin_mem, lsolv);
    CHECK_SPILS_FLAG ("KINSpilsSetLinearSolver", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_spils_set_preconditioner (value vsession,
						  value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    void *mem = KINSOL_MEM_FROM_ML (vsession);

#if 400 <= SUNDIALS_LIB_VERSION
    KINLsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = KINSetPreconditioner (mem, setup, precsolvefn);
    CHECK_LS_FLAG ("KINSetPreconditioner", flag);
#else
    KINSpilsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;
    int flag = KINSpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_SPILS_FLAG ("KINSpilsSetPreconditioner", flag);
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_spils_set_jac_times_vec_fn(value vdata, value vset_jac)
{
    CAMLparam2(vdata, vset_jac);
#if 400 <= SUNDIALS_LIB_VERSION
    KINLsJacTimesVecFn jac = Bool_val (vset_jac) ? jactimesfn : NULL;
    int flag = KINSetJacTimesVecFn(KINSOL_MEM_FROM_ML(vdata), jac);
    CHECK_LS_FLAG("KINSetJacTimesVecFn", flag);
#else
    KINSpilsJacTimesVecFn jac = Bool_val (vset_jac) ? jactimesfn : NULL;
    int flag = KINSpilsSetJacTimesVecFn(KINSOL_MEM_FROM_ML(vdata), jac);
    CHECK_SPILS_FLAG("KINSpilsSetJacTimesVecFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_no_res_mon(value vkin_mem, value vnonniresmon)
{
    CAMLparam2(vkin_mem, vnonniresmon);

    int flag = KINSetNoResMon(KINSOL_MEM_FROM_ML(vkin_mem),
			      Bool_val(vnonniresmon));
    CHECK_FLAG("KINSetNoResMon", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_max_sub_setup_calls(value vkin_mem, value vmsbsetsub)
{
    CAMLparam2(vkin_mem, vmsbsetsub);

    int flag = KINSetMaxSubSetupCalls(KINSOL_MEM_FROM_ML(vkin_mem),
			              Long_val(vmsbsetsub));
    CHECK_FLAG("KINSetMaxSubSetupCalls", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_constraints(value vkin_mem, value vconstraints)
{
    CAMLparam2(vkin_mem, vconstraints);
    int flag;
    N_Vector constraints = NVEC_VAL(vconstraints);
    flag = KINSetConstraints(KINSOL_MEM_FROM_ML(vkin_mem), constraints);

    CHECK_FLAG("KINSetConstraints", flag);

    CAMLreturn (Val_unit);
}

/* basic interface */

/* KINCreate() + KINInit().  */
CAMLprim value sunml_kinsol_init(value weakref, value vtemp,
			     value vomaxiters, value vomaa)
{
    CAMLparam4(weakref, vtemp, vomaxiters, vomaa);
    CAMLlocal2(r, vkin_mem);
    int flag;
    value *backref;
    void *kin_mem;
    N_Vector temp;

    kin_mem = KINCreate();
    if (kin_mem == NULL)
	caml_failwith("KINCreate returned NULL");

    vkin_mem = caml_alloc_final(1, NULL, 1, 5);
    KINSOL_MEM(vkin_mem) = kin_mem;

    if (vomaxiters != Val_none) {
	flag = KINSetNumMaxIters(kin_mem, Long_val(Some_val(vomaxiters)));
	CHECK_FLAG("KINSetNumMaxIters", flag);
    }

#if 260 <= SUNDIALS_LIB_VERSION
    if (vomaa != Val_none) {
	flag = KINSetMAA(kin_mem, Long_val(Some_val(vomaa)));
	CHECK_FLAG("KINSetMAA", flag);
    }
#endif

    temp = NVEC_VAL(vtemp);
    flag = KINInit(kin_mem, sysfn, temp);
    if (flag != KIN_SUCCESS) {
	/* As of SUNDIALS 2.5.0, KINInit frees kin_mem upon failure,
	 * but only if the failure is due to allocation of vectors.
	 * IDAInit and KINInit never frees the mem pointer.
	 * Confusing :( */
	if (flag != KIN_MEM_FAIL) KINFree (&kin_mem);
	CHECK_FLAG("KINInit", flag);
    }

    backref = sunml_sundials_malloc_value(weakref);
    if (backref == NULL) {
	KINFree (&kin_mem);
	caml_raise_out_of_memory();
    }
    KINSetUserData (kin_mem, backref);

    r = caml_alloc_tuple (2);
    Store_field (r, 0, vkin_mem);
    Store_field (r, 1, (value)backref);

    CAMLreturn(r);
}

CAMLprim value sunml_kinsol_solve(value vdata, value vu, value vstrategy,
	 		      value vuscale, value vfscale)
{
    CAMLparam5(vdata, vu, vstrategy, vuscale, vfscale);
    CAMLlocal1(ret);
    int flag;
    enum kinsol_result_tag result = -1;

    N_Vector u = NVEC_VAL (vu);
    N_Vector uscale = NVEC_VAL (vuscale);
    N_Vector fscale = NVEC_VAL (vfscale);

    int strategy = KIN_NONE;
    switch (Int_val(vstrategy)) {
    case VARIANT_KINSOL_STRATEGY_NEWTON:
	strategy = KIN_NONE;
	break;

    case VARIANT_KINSOL_STRATEGY_LINESEARCH:
	strategy = KIN_LINESEARCH;
	break;

    case VARIANT_KINSOL_STRATEGY_PICARD:
#if 260 <= SUNDIALS_LIB_VERSION
	strategy = KIN_PICARD;
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	break;

    case VARIANT_KINSOL_STRATEGY_FIXEDPOINT:
#if 260 <= SUNDIALS_LIB_VERSION
	strategy = KIN_FP;
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
	break;
    }

    flag = KINSol(KINSOL_MEM_FROM_ML(vdata), u, strategy, uscale, fscale);
    CHECK_FLAG("KINSol", flag);

    switch (flag) {
    case KIN_SUCCESS:
	result = VARIANT_KINSOL_RESULT_SUCCESS;
	break;

    case KIN_INITIAL_GUESS_OK:
	result = VARIANT_KINSOL_RESULT_INITIAL_GUESS_OK;
	break;

    case KIN_STEP_LT_STPTOL:
	result = VARIANT_KINSOL_RESULT_STOPPED_ON_STEP_TOL;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_KINSOL_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_KINSOL_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	sunml_kinsol_check_flag("KIN", flag, KINSOL_MEM_FROM_ML(vdata));
    }

    assert (Field (vdata, RECORD_KINSOL_SESSION_EXN_TEMP) == Val_none);

    CAMLreturn (Val_int (result));
}

static value sunml_kinsol_last_lin_exception(void *kin_mem)
{
#if 400 <= SUNDIALS_LIB_VERSION
    long int lsflag;

    if (kin_mem != NULL
	  && KINGetLastLinFlag(kin_mem, &lsflag) == KINLS_SUCCESS)
	return sunml_lsolver_exception_from_flag(lsflag);
#endif
    return Val_none;
}

void sunml_kinsol_check_flag(const char *call, int flag, void *kin_mem)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == KIN_SUCCESS
	|| flag == KIN_INITIAL_GUESS_OK
	|| flag == KIN_STEP_LT_STPTOL
	|| flag == KIN_WARNING) return;

    switch (flag) {
    case KIN_MEM_FAIL:
	caml_raise_out_of_memory();

    case KIN_ILL_INPUT:
        caml_raise_constant(KINSOL_EXN(IllInput));

    case KIN_LINESEARCH_NONCONV:
        caml_raise_constant(KINSOL_EXN(LineSearchNonConvergence));

    case KIN_MAXITER_REACHED:
        caml_raise_constant(KINSOL_EXN(MaxIterationsReached));

    case KIN_MXNEWT_5X_EXCEEDED:
        caml_raise_constant(KINSOL_EXN(MaxNewtonStepExceeded));

    case KIN_LINESEARCH_BCFAIL:
        caml_raise_constant(KINSOL_EXN(LineSearchBetaConditionFailure));

    case KIN_LINSOLV_NO_RECOVERY:
        caml_raise_constant(KINSOL_EXN(LinearSolverNoRecovery));

    case KIN_LINIT_FAIL:
        caml_raise_constant(KINSOL_EXN(LinearSolverInitFailure));

    case KIN_LSETUP_FAIL:
        caml_raise_with_arg(KINSOL_EXN(LinearSetupFailure),
			    sunml_kinsol_last_lin_exception(kin_mem));

    case KIN_LSOLVE_FAIL:
        caml_raise_with_arg(KINSOL_EXN(LinearSolveFailure),
			    sunml_kinsol_last_lin_exception(kin_mem));

    case KIN_SYSFUNC_FAIL:
        caml_raise_constant(KINSOL_EXN(SystemFunctionFailure));

    case KIN_FIRST_SYSFUNC_ERR:
        caml_raise_constant(KINSOL_EXN(FirstSystemFunctionFailure));

    case KIN_REPTD_SYSFUNC_ERR:
        caml_raise_constant(KINSOL_EXN(RepeatedSystemFunctionFailure));

#if 400 <= SUNDIALS_LIB_VERSION
    case KIN_VECTOROP_ERR:
	caml_raise_constant(KINSOL_EXN(VectorOpErr));
#endif

    default:
	/* KIN_MEM_NULL, KIN_NO_MALLOC */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 KINGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

#if 400 <= SUNDIALS_LIB_VERSION
void sunml_kinsol_check_ls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == KINLS_SUCCESS) return;

    switch (flag) {
	case KINLS_ILL_INPUT:
	    caml_raise_constant(KINSOL_EXN(IllInput));

	case KINLS_MEM_FAIL:
	    caml_raise_out_of_memory();

	case KINLS_SUNMAT_FAIL:
	case KINLS_JACFUNC_ERR:
	case KINLS_SUNLS_FAIL:
	default:
	    /* e.g. KINLS_MEM_NULL, KINLS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    KINGetLinReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#else
void sunml_kinsol_check_dls_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == KINDLS_SUCCESS) return;

    switch (flag) {
	case KINDLS_ILL_INPUT:
	    caml_raise_constant(KINSOL_EXN(IllInput));

	case KINDLS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if 300 <= SUNDIALS_LIB_VERSION
	case KINDLS_SUNMAT_FAIL:
	case KINDLS_JACFUNC_ERR:
#endif
	default:
	    /* e.g. KINDLS_MEM_NULL, KINDLS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    KINDlsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

void sunml_kinsol_check_spils_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == KINSPILS_SUCCESS) return;

    switch (flag) {
	case KINSPILS_ILL_INPUT:
	    caml_raise_constant(KINSOL_EXN(IllInput));

	case KINSPILS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if 300 <= SUNDIALS_LIB_VERSION
	case KINSPILS_SUNLS_FAIL:
#endif
	default:
	    /* e.g. KINSPILS_MEM_NULL, KINSPILS_PMEM_NULL, KINSPILS_LMEM_NULL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		     KINSpilsGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}
#endif

/* basic interface */

CAMLprim value sunml_kinsol_session_finalize(value vdata)
{
    if (KINSOL_MEM_FROM_ML(vdata) != NULL) {
	void *kin_mem = KINSOL_MEM_FROM_ML(vdata);
	value *backref = KINSOL_BACKREF_FROM_ML(vdata);
	KINFree(&kin_mem);
	sunml_sundials_free_value(backref);
    }

    return Val_unit;
}

/* boiler plate */

CAMLprim value sunml_kinsol_spils_set_max_restarts(value vkin_mem, value vmaxr)
{
    CAMLparam2(vkin_mem, vmaxr);
#if SUNDIALS_LIB_VERSION < 300
    int flag = KINSpilsSetMaxRestarts(KINSOL_MEM_FROM_ML(vkin_mem),
				      Int_val(vmaxr));
    CHECK_SPILS_FLAG("KINSpilsSetMaxRestarts", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_spils_spgmr(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSpgmr (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSpgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_spils_spfgmr(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
#if 260 <= SUNDIALS_LIB_VERSION && SUNDIALS_LIB_VERSION < 300
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSpfgmr (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSpfgmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_spils_spbcgs(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSpbcg (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSpbcg", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_spils_sptfqmr(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
#if SUNDIALS_LIB_VERSION < 300
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSptfqmr (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSptfqmr", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_dls_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetLinWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrwLS, &leniwLS);
    CHECK_LS_FLAG("KINGetLinWorkSpace", flag);
#else
    int flag = KINDlsGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrwLS, &leniwLS);
    CHECK_DLS_FLAG("KINDlsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));

    CAMLreturn(r);
}

CAMLprim value sunml_kinsol_get_num_jac_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumJacEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetNumJacEvals", flag);
#else
    int flag = KINDlsGetNumJacEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_DLS_FLAG("KINDlsGetNumJacEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_kinsol_dls_get_num_func_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumLinFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetLinNumFuncEvals", flag);
#else
    int flag = KINDlsGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_DLS_FLAG("KINDlsGetNumFuncEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}


CAMLprim value sunml_kinsol_spils_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    long int lenrw;
    long int leniw;

#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetLinWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrw, &leniw);
    CHECK_LS_FLAG("KINGetLinWorkSpace", flag);
#else
    int flag = KINSpilsGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrw, &leniw);
    CHECK_SPILS_FLAG("KINSpilsGetWorkSpace", flag);
#endif

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_kinsol_get_num_lin_iters(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumLinIters(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetNumLinIters", flag);
#else
    int flag = KINSpilsGetNumLinIters(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_SPILS_FLAG("KINSpilsGetNumLinIters", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_kinsol_get_num_lin_conv_fails(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumLinConvFails(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetNumLinConvFails", flag);
#else
    int flag = KINSpilsGetNumConvFails(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_SPILS_FLAG("KINSpilsGetNumConvFails", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_kinsol_get_num_prec_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumPrecEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetNumPrecEvals", flag);
#else
    int flag = KINSpilsGetNumPrecEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_SPILS_FLAG("KINSpilsGetNumPrecEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_kinsol_get_num_prec_solves(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumPrecSolves(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetNumPrecSolves", flag);
#else
    int flag = KINSpilsGetNumPrecSolves(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_SPILS_FLAG("KINSpilsGetNumPrecSolves", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_kinsol_get_num_jtimes_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumJtimesEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetNumJtimesEvals", flag);
#else
    int flag = KINSpilsGetNumJtimesEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_SPILS_FLAG("KINSpilsGetNumJtimesEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_kinsol_spils_get_num_func_evals (value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
#if 400 <= SUNDIALS_LIB_VERSION
    int flag = KINGetNumLinFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_LS_FLAG("KINGetNumLinFuncEvals", flag);
#else
    int flag = KINSpilsGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_SPILS_FLAG("KINSpilsGetNumFuncEvals", flag);
#endif

    CAMLreturn(Val_long(r));
}

CAMLprim value sunml_kinsol_set_error_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

    int flag = KINSetErrFile(KINSOL_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("KINSetErrFile", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_info_file(value vdata, value vfile)
{
    CAMLparam2(vdata, vfile);

    int flag = KINSetInfoFile(KINSOL_MEM_FROM_ML(vdata), ML_CFILE(vfile));
    CHECK_FLAG("KINSetInfoFile", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_print_level(value vkin_mem, value vplvl)
{
    CAMLparam2(vkin_mem, vplvl);

    int flag = KINSetPrintLevel(KINSOL_MEM_FROM_ML(vkin_mem), Int_val(vplvl));
    CHECK_FLAG("KINSetPrintLevel", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_damping_aa(value vkin_mem, value vbeta)
{
    CAMLparam2(vkin_mem, vbeta);

    int flag = KINSetDampingAA(KINSOL_MEM_FROM_ML(vkin_mem), Double_val(vbeta));
    CHECK_FLAG("KINSetDampingAA", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_num_max_iters(value vkin_mem, value vmxiter)
{
    CAMLparam2(vkin_mem, vmxiter);

    int flag = KINSetNumMaxIters(KINSOL_MEM_FROM_ML(vkin_mem), Long_val(vmxiter));
    CHECK_FLAG("KINSetNumMaxIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_maa(value vkin_mem, value vmaa)
{
    CAMLparam2(vkin_mem, vmaa);
#if 260 <= SUNDIALS_LIB_VERSION

    int flag = KINSetMAA(KINSOL_MEM_FROM_ML(vkin_mem), Long_val(vmaa));
    CHECK_FLAG("KINSetMAA", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_no_init_setup(value vkin_mem, value vnoinitsetup)
{
    CAMLparam2(vkin_mem, vnoinitsetup);

    int flag = KINSetNoInitSetup(KINSOL_MEM_FROM_ML(vkin_mem), Bool_val(vnoinitsetup));
    CHECK_FLAG("KINSetNoInitSetup", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_max_setup_calls(value vkin_mem, value vmsbset)
{
    CAMLparam2(vkin_mem, vmsbset);

    int flag = KINSetMaxSetupCalls(KINSOL_MEM_FROM_ML(vkin_mem), Long_val(vmsbset));
    CHECK_FLAG("KINSetMaxSetupCalls", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_eta_form(value vkin_mem, value vetachoice)
{
    CAMLparam2(vkin_mem, vetachoice);

    int etachoice = -1;
    if (Is_long(vetachoice)) {
	switch (Int_val(vetachoice)) {
	case VARIANT_KINSOL_ETA_CHOICE1:
	    etachoice = KIN_ETACHOICE1;
	    break;
	}
    } else {
	switch (Tag_val(vetachoice)) {
	case VARIANT_KINSOL_ETA_CHOICE2:
	    etachoice = KIN_ETACHOICE2;
	    break;

	case VARIANT_KINSOL_ETA_CONSTANT:
	    etachoice = KIN_ETACONSTANT;
	    break;
	}
    }

    int flag = KINSetEtaForm(KINSOL_MEM_FROM_ML(vkin_mem), etachoice);
    CHECK_FLAG("KINSetEtaForm", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_eta_const_value(value vkin_mem, value veta)
{
    CAMLparam2(vkin_mem, veta);

    int flag = KINSetEtaConstValue(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(veta));
    CHECK_FLAG("KINSetEtaConstValue", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_eta_params(value vkin_mem, value vegamma, value vealpha)
{
    CAMLparam3(vkin_mem, vegamma, vealpha);

    int flag = KINSetEtaParams(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vegamma), Double_val(vealpha));
    CHECK_FLAG("KINSetEtaParams", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_res_mon_const_value(value vkin_mem, value vomegaconst)
{
    CAMLparam2(vkin_mem, vomegaconst);

    int flag = KINSetResMonConstValue(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vomegaconst));
    CHECK_FLAG("KINSetResMonConstValue", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_res_mon_params(value vkin_mem, value vomegamin, value vomegamax)
{
    CAMLparam3(vkin_mem, vomegamin, vomegamax);

    int flag = KINSetResMonParams(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vomegamin), Double_val(vomegamax));
    CHECK_FLAG("KINSetResMonParams", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_no_min_eps(value vkin_mem, value vnomineps)
{
    CAMLparam2(vkin_mem, vnomineps);

    int flag = KINSetNoMinEps(KINSOL_MEM_FROM_ML(vkin_mem), Bool_val(vnomineps));
    CHECK_FLAG("KINSetNoMinEps", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_max_newton_step(value vkin_mem, value vmxnewtstep)
{
    CAMLparam2(vkin_mem, vmxnewtstep);

    int flag = KINSetMaxNewtonStep(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vmxnewtstep));
    CHECK_FLAG("KINSetMaxNewtonStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_max_beta_fails(value vkin_mem, value vmxnbcf)
{
    CAMLparam2(vkin_mem, vmxnbcf);

    int flag = KINSetMaxBetaFails(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vmxnbcf));
    CHECK_FLAG("KINSetMaxBetaFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_rel_err_func(value vkin_mem, value vrelfunc)
{
    CAMLparam2(vkin_mem, vrelfunc);

    int flag = KINSetRelErrFunc(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vrelfunc));
    CHECK_FLAG("KINSetRelErrFunc", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_func_norm_tol(value vkin_mem, value vfnormtol)
{
    CAMLparam2(vkin_mem, vfnormtol);

    int flag = KINSetFuncNormTol(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vfnormtol));
    CHECK_FLAG("KINSetFuncNormTol", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_set_scaled_step_tol(value vkin_mem, value vscsteptol)
{
    CAMLparam2(vkin_mem, vscsteptol);

    int flag = KINSetScaledStepTol(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vscsteptol));
    CHECK_FLAG("KINSetScaledStepTol", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_kinsol_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = KINGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrw, &leniw);
    CHECK_FLAG("KINGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value sunml_kinsol_get_num_func_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumFuncEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_kinsol_get_num_nonlin_solv_iters(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumNonlinSolvIters(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_kinsol_get_num_beta_cond_fails(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumBetaCondFails(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumBetaCondFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_kinsol_get_num_backtrack_ops(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumBacktrackOps(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumBacktrackOps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value sunml_kinsol_get_func_norm(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    realtype v;

    flag = KINGetFuncNorm(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetFuncNorm", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value sunml_kinsol_get_step_length(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    realtype v;

    flag = KINGetStepLength(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetStepLength", flag);

    CAMLreturn(caml_copy_double(v));
}

