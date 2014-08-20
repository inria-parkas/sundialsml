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

#include <kinsol/kinsol.h>
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
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_spils.h>
#include <kinsol/kinsol_impl.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <kinsol/kinsol_lapack.h>
#endif

#include "dls_ml.h"
#include "spils_ml.h"
#include "kinsol_ml.h"

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

static void errh(
	int error_code,
	const char *module,
	const char *func,
	char *msg,
	void *eh_data)
{
    CAMLparam0();
    CAMLlocal1(a);
    value *backref = eh_data;

    CAML_FN (call_errh);

    a = caml_alloc_tuple(RECORD_SUNDIALS_ERROR_DETAILS_SIZE);
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_CODE,
                Val_int(error_code));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_MODULE_NAME,
                caml_copy_string(module));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_FUNCTION_NAME,
                caml_copy_string(func));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_MESSAGE,
                caml_copy_string(msg));

    caml_callback2(*call_errh, *backref, a);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = KINSetErrHandlerFn(KINSOL_MEM_FROM_ML(vdata), errh,
				  KINSOL_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("KINSetErrHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = KINSetErrHandlerFn(KINSOL_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("KINSetErrHandlerFn", flag);

    CAMLreturn0;
}

static void infoh(
	const char *module,
	const char *func,
	char *msg,
	void *ih_data)
{
    CAMLparam0();
    CAMLlocal1(a);
    value *backref = ih_data;

    CAML_FN (call_infoh);

    a = caml_alloc_tuple(RECORD_SUNDIALS_ERROR_DETAILS_SIZE);
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_CODE, Val_int(0));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_MODULE_NAME,
                caml_copy_string(module));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_FUNCTION_NAME,
                caml_copy_string(func));
    Store_field(a, RECORD_SUNDIALS_ERROR_DETAILS_ERROR_MESSAGE,
                caml_copy_string(msg));

    caml_callback2(*call_infoh, *backref, a);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_info_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = KINSetInfoHandlerFn(KINSOL_MEM_FROM_ML(vdata), infoh,
				   KINSOL_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("KINSetInfoHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_clear_info_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = KINSetInfoHandlerFn(KINSOL_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("KINSetInfoHandlerFn", flag);

    CAMLreturn0;
}

static int sysfn(N_Vector uu, N_Vector val, void *user_data)
{
    CAMLparam0();
    CAMLlocal2(vuu, vval);
    int r;
    value *backref = user_data;
    CAML_FN (call_sysfn);

    vuu = NVEC_BACKLINK(uu);
    vval = NVEC_BACKLINK(val);

    // The data payloads inside vuu and vval are only valid during this
    // call, afterward that memory goes back to kinsol. These bigarrays must
    // not be retained by closure_rhsfn! If it wants a permanent copy, then
    // it has to make it manually.
    r = Int_val (caml_callback3(*call_sysfn, *backref, vuu, vval));

    CAMLreturnT(int, r);
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

static value make_jac_arg(N_Vector u, N_Vector fu, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_KINSOL_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_U, NVEC_BACKLINK(u));
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_FU, NVEC_BACKLINK(fu));
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

static value make_double_tmp(N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, NVEC_BACKLINK(tmp1));
    Store_field(r, 1, NVEC_BACKLINK(tmp2));
    CAMLreturn(r);
}

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
    CAMLlocalN (args, 3);
    int r;
    value *backref = user_data;
    CAML_FN (call_jacfn);

    args[0] = *backref;
    args[1] = make_jac_arg(u, fu, make_double_tmp(tmp1, tmp2));
    args[2] = caml_alloc_final (2, NULL, 0, 1);
    DLSMAT(args[2]) = Jac;

    r = Int_val (caml_callbackN (*call_jacfn,
				 sizeof (args) / sizeof (*args),
				 args));

    CAMLreturnT(int, r);
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
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_bandjacfn);

    args[0] = *backref;
    args[1] = caml_alloc_tuple(RECORD_KINSOL_BANDRANGE_SIZE);
    Store_field(args[1], RECORD_KINSOL_BANDRANGE_MUPPER, Val_long(mupper));
    Store_field(args[1], RECORD_KINSOL_BANDRANGE_MLOWER, Val_long(mlower));
    args[2] = make_jac_arg(u, fu, make_double_tmp(tmp1, tmp2));
    args[3] = caml_alloc_final(2, NULL, 0, 1);
    DLSMAT(args[3]) = Jac;

    r = Int_val (caml_callbackN(*call_bandjacfn,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int precsetupfn(
    N_Vector uu,
    N_Vector uscale,
    N_Vector fu,
    N_Vector fscale,
    void *user_data,
    N_Vector tmp1,
    N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal2(session, r);
    CAMLlocalN(args, 3);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_precsetupfn);

    args[0] = *backref;
    args[1] = make_jac_arg(uu, fu, make_double_tmp(tmp1, tmp2));
    args[2] = make_prec_solve_arg(uscale, fscale);

    retcode = Int_val (caml_callbackN(*call_precsetupfn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    CAMLreturnT(int, retcode);
}

static int precsolvefn(
	N_Vector uu,
	N_Vector uscale,
	N_Vector fu,
	N_Vector fscale,
	N_Vector vv,
	void *user_data,
	N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal1(rv);
    CAMLlocalN(args, 4);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_precsolvefn);

    args[0] = *backref;
    args[1] = make_jac_arg(uu, fu, NVEC_BACKLINK(tmp));
    args[2] = make_prec_solve_arg(uscale, fscale);
    args[3] = NVEC_BACKLINK(vv);

    retcode = Int_val (caml_callbackN(*call_precsolvefn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    CAMLreturnT(int, retcode);
}

static int jactimesfn(
	N_Vector v,
	N_Vector Jv,
	N_Vector u,
	booleantype *new_uu,
	void *user_data)
{
    CAMLparam0();
    CAMLlocal1(vr);
    CAMLlocalN (args, 5);
    value *backref = user_data;
    CAML_FN (call_jactimesfn);
    int r;

    args[0] = *backref;
    args[1] = NVEC_BACKLINK(v);
    args[2] = NVEC_BACKLINK(Jv);
    args[3] = NVEC_BACKLINK(u);
    args[4] = Val_bool(*new_uu);

    vr = caml_callbackN (*call_jactimesfn,
			 sizeof (args) / sizeof (*args),
			 args);
    r = Field(vr, 1);
    if (r == 0) *new_uu = Bool_val(Field(vr, 0));

    CAMLreturnT(int, r);
}

static int linit(KINMem kin_mem)
{
    CAMLparam0();
    int r;
    value *backref = kin_mem->kin_user_data;
    CAML_FN (call_linit);

    r = Int_val (caml_callback(*call_linit, *backref));

    CAMLreturnT(int, r);
}

static int lsetup(KINMem kin_mem)
{
    CAMLparam0();
    int r;
    value *backref = kin_mem->kin_user_data;
    CAML_FN (call_lsetup);

    r = Int_val (caml_callback(*call_lsetup, *backref));

    CAMLreturnT(int, r);
}

static int lsolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal1(vr);
    value *backref = kin_mem->kin_user_data;
    CAML_FN (call_lsolve);

    args[0] = *backref;
    args[1] = NVEC_BACKLINK(x);
    args[2] = NVEC_BACKLINK(b);

    vr = caml_callbackN(*call_lsolve, sizeof (args) / sizeof (*args), args);
    *res_norm = Double_val(Field(vr, 0));

    CAMLreturnT(int, Int_val(Field(vr, 1)));
}

static void lfree(KINMem kin_mem)
{
    CAMLparam0();
    value *backref = kin_mem->kin_user_data;
    CAML_FN (call_lfree);

    caml_callback(*call_lfree, *backref);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_alternate (value vkin_mem, value vhas_init,
				      value vhas_setup, value vhas_free)
{
    CAMLparam4(vkin_mem, vhas_init, vhas_setup, vhas_free);
    KINMem kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);

    kin_mem->kin_linit  = Bool_val(vhas_init)  ? linit : NULL;
    kin_mem->kin_lsetup = Bool_val(vhas_setup) ? lsetup : NULL;
    kin_mem->kin_setupNonNull = Bool_val(vhas_setup);
    kin_mem->kin_lsolve = lsolve;
    kin_mem->kin_lfree  = Bool_val(vhas_free)  ? lfree : NULL;
    kin_mem->kin_lmem   = NULL;

    CAMLreturn0;
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim void c_kinsol_dls_dense (value vkin_mem, value vset_jac)
{
    CAMLparam2(vkin_mem, vset_jac);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINDense (kin_mem, neqs);
    CHECK_FLAG ("KINDense", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetDenseJacFn(KINSOL_MEM_FROM_ML(vkin_mem), jacfn);
	CHECK_FLAG("KINDlsSetDenseJacFn", flag);
    }
    CAMLreturn0;
}

CAMLprim void c_kinsol_dls_lapack_dense (value vkin_mem, value vset_jac)
{
    CAMLparam2 (vkin_mem, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINLapackDense (kin_mem, neqs);
    CHECK_FLAG ("KINLapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetDenseJacFn (KINSOL_MEM_FROM_ML (vkin_mem), jacfn);
	CHECK_FLAG("KINDlsSetDenseJacFn", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn0;
}

CAMLprim void c_kinsol_dls_set_dense_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetDenseJacFn(KINSOL_MEM_FROM_ML(vdata), jacfn);
    CHECK_FLAG("KINDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_kinsol_dls_clear_dense_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetDenseJacFn(KINSOL_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("KINDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_kinsol_dls_band (value vkin_mem,
				 value vmupper,
				 value vmlower,
				 value vset_jac)
{
    CAMLparam4(vkin_mem, vmupper, vmlower, vset_jac);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINBand (kin_mem, neqs, Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("KINBand", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vkin_mem), bandjacfn);
	CHECK_FLAG("KINDlsSetBandJacFn", flag);
    }
    CAMLreturn0;
}

CAMLprim void c_kinsol_dls_lapack_band (value vkin_mem, value vmupper,
					value vmlower, value vset_jac)
{
    CAMLparam4(vkin_mem, vmupper, vmlower, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINLapackBand (kin_mem, neqs,
			  Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("KINLapackBand", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vkin_mem), bandjacfn);
	CHECK_FLAG("KINDlsSetBandJacFn", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn0;
}

CAMLprim void c_kinsol_dls_set_band_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vdata), bandjacfn);
    CHECK_FLAG("KINDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_kinsol_dls_clear_band_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("KINDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_set_preconditioner (value vsession,
						value vset_precsetup)
{
    CAMLparam2 (vsession, vset_precsetup);
    int flag;
    void *mem = KINSOL_MEM_FROM_ML (vsession);
    KINSpilsPrecSetupFn setup = Bool_val (vset_precsetup) ? precsetupfn : NULL;

    flag = KINSpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_FLAG ("KINSpilsSetPreconditioner", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_set_jac_times_vec_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINSpilsSetJacTimesVecFn(KINSOL_MEM_FROM_ML(vdata), jactimesfn);
    CHECK_FLAG("KINSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_clear_jac_times_vec_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINSpilsSetJacTimesVecFn(KINSOL_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("KINSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

CAMLprim void c_kinsol_set_no_res_mon(value vkin_mem, value vnonniresmon)
{
    CAMLparam2(vkin_mem, vnonniresmon);

    int flag = KINSetNoResMon(KINSOL_MEM_FROM_ML(vkin_mem),
			      Bool_val(vnonniresmon));
    CHECK_FLAG("KINSetNoResMon", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_max_sub_setup_calls(value vkin_mem, value vmsbsetsub)
{
    CAMLparam2(vkin_mem, vmsbsetsub);

    int flag = KINSetMaxSubSetupCalls(KINSOL_MEM_FROM_ML(vkin_mem),
			              Long_val(vmsbsetsub));
    CHECK_FLAG("KINSetMaxSubSetupCalls", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_constraints(value vkin_mem, value vconstraints)
{
    CAMLparam2(vkin_mem, vconstraints);
    int flag;
    N_Vector constraints;

    // TODO
#if SAFETY_CHECKS && KINSOL_ML_BIGARRAYS
    int neqs = KINSOL_NEQS_FROM_ML (vkin_mem);

    if (neqs != Caml_ba_array_val(vconstraints)->dim[0])
	caml_invalid_argument (
	    "Kinsol.set_constraints: constraints vector has incorrect length");
#endif

    constraints = NVEC_VAL(vconstraints);
    flag = KINSetConstraints(KINSOL_MEM_FROM_ML(vkin_mem), constraints);

    CHECK_FLAG("KINSetConstraints", flag);

    CAMLreturn0;
}

/* basic interface */

/* KINCreate() + KINInit().  */
CAMLprim value c_kinsol_init(value weakref, value vtemp)
{
    CAMLparam2(weakref, vtemp);
    CAMLlocal1(r);
    int flag;
    value *backref;
    void *kin_mem;
    N_Vector temp;

    if (sizeof(int) != 4) {
	caml_failwith("The library assumes that an int (in C) has 32-bits.");
    }

    kin_mem = KINCreate();
    if (kin_mem == NULL)
	caml_failwith("KINCreate returned NULL");
    temp = NVEC_VAL(vtemp);
    flag = KINInit(kin_mem, sysfn, temp);
    if (flag != KIN_SUCCESS) {
	KINFree (kin_mem);
	CHECK_FLAG("KINInit", flag);
    }

    backref = malloc (sizeof (*backref));
    if (backref == NULL) {
	KINFree (kin_mem);
	caml_raise_out_of_memory();
    }
    *backref = weakref;
    caml_register_generational_global_root (backref);
    KINSetUserData (kin_mem, backref);

    r = caml_alloc_tuple (4);
    Store_field (r, 0, (value)kin_mem);
    Store_field (r, 1, (value)backref);
    Store_field (r, 2, Val_long (0)); /* no err_file = NULL */
    Store_field (r, 3, Val_long (0)); /* no info_file = NULL */

    CAMLreturn(r);
}

CAMLprim value c_kinsol_solve(value vdata, value vu, value vlinesearch,
	 		      value vuscale, value vfscale)
{
    CAMLparam5(vdata, vu, vlinesearch, vuscale, vfscale);
    CAMLlocal1(ret);
    int flag;
    enum kinsol_result_tag result;
    N_Vector u;
    N_Vector uscale;
    N_Vector fscale;

    // TODO
#if SAFETY_CHECKS && KINSOL_ML_BIGARRAYS
    int neqs = KINSOL_NEQS_FROM_ML (vdata);

    /* This can't be checked for generic nvectors.  */
    if (neqs != Caml_ba_array_val(vu)->dim[0])
	caml_invalid_argument ("Kinsol.solve: u vector has incorrect length");
    if (neqs != Caml_ba_array_val(vuscale)->dim[0])
	caml_invalid_argument ("Kinsol.solve: uscale vector has incorrect length");
    if (neqs != Caml_ba_array_val(vfscale)->dim[0])
	caml_invalid_argument ("Kinsol.solve: fscale vector has incorrect length");
#endif

    u = NVEC_VAL (vu);
    uscale = NVEC_VAL (vuscale);
    fscale = NVEC_VAL (vfscale);
    flag = KINSol(KINSOL_MEM_FROM_ML(vdata),
		  u,
		  Bool_val (vlinesearch) ? KIN_LINESEARCH : KIN_NONE,
		  uscale,
		  fscale);

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
	CHECK_FLAG ("KIN", flag);
    }

    /* Hmm...should this go in the production code or not?  */
    if (Is_block (Field (vdata, RECORD_KINSOL_SESSION_EXN_TEMP)))
	abort ();

    CAMLreturn (Val_int (result));
}


void kinsol_ml_check_flag(const char *call, int flag)
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
        caml_raise_constant(*caml_named_value("kinsol_IllInput"));

    case KIN_LINESEARCH_NONCONV:
        caml_raise_constant(*caml_named_value(
		    "kinsol_LineSearchNonConvergence"));

    case KIN_MAXITER_REACHED:
        caml_raise_constant(*caml_named_value("kinsol_MaxIterationsReached"));

    case KIN_MXNEWT_5X_EXCEEDED:
        caml_raise_constant(*caml_named_value("kinsol_MaxNewtonStepExceeded"));

    case KIN_LINESEARCH_BCFAIL:
        caml_raise_constant(*caml_named_value(
		    "kinsol_LineSearchBetaConditionFailure"));

    case KIN_LINSOLV_NO_RECOVERY:
        caml_raise_constant(*caml_named_value("kinsol_LinearSolverNoRecovery"));

    case KIN_LINIT_FAIL:
        caml_raise_constant(*caml_named_value(
		    "kinsol_LinearSolverInitFailure"));

    case KIN_LSETUP_FAIL:
        caml_raise_constant(*caml_named_value("kinsol_LinearSetupFailure"));

    case KIN_LSOLVE_FAIL:
        caml_raise_constant(*caml_named_value("kinsol_LinearSolverFailure"));

    case KIN_SYSFUNC_FAIL:
        caml_raise_constant(*caml_named_value("kinsol_SystemFunctionFailure"));

    case KIN_FIRST_SYSFUNC_ERR:
        caml_raise_constant(*caml_named_value(
		    "kinsol_FirstSystemFunctionFailure"));

    case KIN_REPTD_SYSFUNC_ERR:
        caml_raise_constant(*caml_named_value(
		    "kinsol_RepeatedSystemFunctionFailure"));

    default:
	/* KIN_MEM_NULL, KIN_NO_MALLOC */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: unexpected error code", call);
	caml_failwith(exmsg);
    }
}

/* basic interface */

CAMLprim void c_kinsol_session_finalize(value vdata)
{
    if (KINSOL_MEM_FROM_ML(vdata) != NULL) {
	void *kin_mem = KINSOL_MEM_FROM_ML(vdata);
	value *backref = KINSOL_BACKREF_FROM_ML(vdata);
	KINFree(&kin_mem);
	caml_remove_generational_global_root (backref);
	free (backref);
    }

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_KINSOL_SESSION_ERRFILE));
    if (err_file != NULL) {
	fclose(err_file);
    }

    FILE* info_file =
      (FILE *)Long_val(Field(vdata, RECORD_KINSOL_SESSION_INFOFILE));
    if (info_file != NULL) {
	fclose(info_file);
    }
}

/* boiler plate */

CAMLprim void c_kinsol_spils_spgmr(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSpgmr (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSpgmr", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_spbcg(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSpbcg (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSpbcg", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_sptfqmr(value vkin_mem, value vmaxl)
{
    CAMLparam2(vkin_mem, vmaxl);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINSptfqmr (kin_mem, Int_val (vmaxl));
    CHECK_FLAG ("KINSptfqmr", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_spils_set_max_restarts(value vkin_mem, value vmaxrs)
{
    CAMLparam2(vkin_mem, vmaxrs);

    int flag = KINSpilsSetMaxRestarts(KINSOL_MEM_FROM_ML(vkin_mem), Int_val(vmaxrs));
    CHECK_FLAG("KINSetMaxRestarts", flag);

    CAMLreturn0;
}

CAMLprim value c_kinsol_dls_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = KINDlsGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrwLS, &leniwLS);
    CHECK_FLAG("KINDlsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));

    CAMLreturn(r);
}

CAMLprim value c_kinsol_dls_get_num_jac_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINDlsGetNumJacEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINDlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_dls_get_num_func_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINDlsGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINDlsGetNumFuncEvals", flag);

    CAMLreturn(Val_long(r));
}


CAMLprim value c_kinsol_spils_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = KINSpilsGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrw, &leniw);
    CHECK_FLAG("KINSpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
}

CAMLprim value c_kinsol_spils_get_num_lin_iters(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumLinIters(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_conv_fails(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumConvFails(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_prec_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumPrecEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_prec_solves(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumPrecSolves(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_jtimes_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumJtimesEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_kinsol_spils_get_num_func_evals (value vkin_mem)
{
    CAMLparam1(vkin_mem);

    long int r;
    int flag = KINSpilsGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &r);
    CHECK_FLAG("KINSpilsGetNumFuncEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim void c_kinsol_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_KINSOL_SESSION_ERRFILE));

    if (err_file != NULL) {
	fclose(err_file);
	Store_field(vdata, RECORD_KINSOL_SESSION_ERRFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    err_file = fopen(String_val(vpath), mode);
    if (err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = KINSetErrFile(KINSOL_MEM_FROM_ML(vdata), err_file);
    CHECK_FLAG("KINSetErrFile", flag);

    Store_field(vdata, RECORD_KINSOL_SESSION_ERRFILE, Val_long(err_file));

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_info_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    FILE* info_file =
      (FILE *)Long_val(Field(vdata, RECORD_KINSOL_SESSION_INFOFILE));

    if (info_file != NULL) {
	fclose(info_file);
	Store_field(vdata, RECORD_KINSOL_SESSION_INFOFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    info_file = fopen(String_val(vpath), mode);
    if (info_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = KINSetInfoFile(KINSOL_MEM_FROM_ML(vdata), info_file);
    CHECK_FLAG("KINSetInfoFile", flag);

    Store_field(vdata, RECORD_KINSOL_SESSION_INFOFILE, Val_long(info_file));

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_print_level(value vkin_mem, value vplvl)
{
    CAMLparam2(vkin_mem, vplvl);

    int flag = KINSetPrintLevel(KINSOL_MEM_FROM_ML(vkin_mem), Int_val(vplvl));
    CHECK_FLAG("KINSetPrintLevel", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_num_max_iters(value vkin_mem, value vmxiter)
{
    CAMLparam2(vkin_mem, vmxiter);

    int flag = KINSetNumMaxIters(KINSOL_MEM_FROM_ML(vkin_mem), Long_val(vmxiter));
    CHECK_FLAG("KINSetNumMaxIters", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_no_init_setup(value vkin_mem, value vnoinitsetup)
{
    CAMLparam2(vkin_mem, vnoinitsetup);

    int flag = KINSetNoInitSetup(KINSOL_MEM_FROM_ML(vkin_mem), Bool_val(vnoinitsetup));
    CHECK_FLAG("KINSetNoInitSetup", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_max_setup_calls(value vkin_mem, value vmsbset)
{
    CAMLparam2(vkin_mem, vmsbset);

    int flag = KINSetMaxSetupCalls(KINSOL_MEM_FROM_ML(vkin_mem), Long_val(vmsbset));
    CHECK_FLAG("KINSetMaxSetupCalls", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_eta_form(value vkin_mem, value vetachoice)
{
    CAMLparam2(vkin_mem, vetachoice);

    int etachoice;
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

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_eta_const_value(value vkin_mem, value veta)
{
    CAMLparam2(vkin_mem, veta);

    int flag = KINSetEtaConstValue(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(veta));
    CHECK_FLAG("KINSetEtaConstValue", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_eta_params(value vkin_mem, value vegamma, value vealpha)
{
    CAMLparam3(vkin_mem, vegamma, vealpha);

    int flag = KINSetEtaParams(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vegamma), Double_val(vealpha));
    CHECK_FLAG("KINSetEtaParams", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_res_mon_const_value(value vkin_mem, value vomegaconst)
{
    CAMLparam2(vkin_mem, vomegaconst);

    int flag = KINSetResMonConstValue(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vomegaconst));
    CHECK_FLAG("KINSetResMonConstValue", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_res_mon_params(value vkin_mem, value vomegamin, value vomegamax)
{
    CAMLparam3(vkin_mem, vomegamin, vomegamax);

    int flag = KINSetResMonParams(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vomegamin), Double_val(vomegamax));
    CHECK_FLAG("KINSetResMonParams", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_no_min_eps(value vkin_mem, value vnomineps)
{
    CAMLparam2(vkin_mem, vnomineps);

    int flag = KINSetNoMinEps(KINSOL_MEM_FROM_ML(vkin_mem), Bool_val(vnomineps));
    CHECK_FLAG("KINSetNoMinEps", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_max_newton_step(value vkin_mem, value vmxnewtstep)
{
    CAMLparam2(vkin_mem, vmxnewtstep);

    int flag = KINSetMaxNewtonStep(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vmxnewtstep));
    CHECK_FLAG("KINSetMaxNewtonStep", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_max_beta_fails(value vkin_mem, value vmxnbcf)
{
    CAMLparam2(vkin_mem, vmxnbcf);

    int flag = KINSetMaxBetaFails(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vmxnbcf));
    CHECK_FLAG("KINSetMaxBetaFails", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_rel_err_func(value vkin_mem, value vrelfunc)
{
    CAMLparam2(vkin_mem, vrelfunc);

    int flag = KINSetRelErrFunc(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vrelfunc));
    CHECK_FLAG("KINSetRelErrFunc", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_func_norm_tol(value vkin_mem, value vfnormtol)
{
    CAMLparam2(vkin_mem, vfnormtol);

    int flag = KINSetFuncNormTol(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vfnormtol));
    CHECK_FLAG("KINSetFuncNormTol", flag);

    CAMLreturn0;
}

CAMLprim void c_kinsol_set_scaled_step_tol(value vkin_mem, value vscsteptol)
{
    CAMLparam2(vkin_mem, vscsteptol);

    int flag = KINSetScaledStepTol(KINSOL_MEM_FROM_ML(vkin_mem),
	    Double_val(vscsteptol));
    CHECK_FLAG("KINSetScaledStepTol", flag);

    CAMLreturn0;
}

CAMLprim value c_kinsol_get_work_space(value vkin_mem)
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

CAMLprim value c_kinsol_get_num_func_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumFuncEvals(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumFuncEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_num_nonlin_solv_iters(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumNonlinSolvIters(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_num_beta_cond_fails(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumBetaCondFails(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumBetaCondFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_num_backtrack_ops(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINGetNumBacktrackOps(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetNumBacktrackOps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_kinsol_get_func_norm(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    realtype v;

    flag = KINGetFuncNorm(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetFuncNorm", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_kinsol_get_step_length(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    realtype v;

    flag = KINGetStepLength(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINGetStepLength", flag);

    CAMLreturn(caml_copy_double(v));
}

