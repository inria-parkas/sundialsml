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

#include <cvode/cvode.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
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
#include <cvode/cvode_bandpre.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_nvector.h>

#include "dls_ml.h"
#include "spils_ml.h"
#include "sundials_ml.h"
#include "cvode_ml.h"
#include "nvector_ml.h"

#if SUNDIALS_BLAS_LAPACK == 1
#include <cvode/cvode_lapack.h>
#endif

#include <stdio.h>
#define MAX_ERRMSG_LEN 256

/* callbacks */

#define CVTYPE(fname) c_cvode_ ## fname
#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define CVTYPESTR(fname) QUOTE(CVTYPE(fname))

#define CAML_FN(name)					\
    static value *name;					\
    if (name == NULL)					\
	name = caml_named_value (CVTYPESTR (name));

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

CAMLprim void c_cvode_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = CVodeSetErrHandlerFn(CVODE_MEM_FROM_ML(vdata), errh,
				    CVODE_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSetErrHandlerFn(CVODE_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("CVodeSetErrHandlerFn", flag);

    CAMLreturn0;
}

static int check_exception(value session, value r)
{
    CAMLparam2(session, r);
    CAMLlocal1(exn);

    static value *recoverable_failure = NULL;
    if (recoverable_failure == NULL) {
	recoverable_failure =
	    caml_named_value("cvode_RecoverableFailure");
    }

    if (!Is_exception_result(r)) return 0;

    r = Extract_exception(r);

    if (Field(r, 0) == *recoverable_failure)
	CAMLreturnT (int, 1);

    /* Unrecoverable error.  Save the exception and return -1.  */
    exn = caml_alloc_small (1,0);
    Field (exn, 0) = r;
    Store_field (session, RECORD_CVODE_SESSION_EXN_TEMP, exn);
    CAMLreturnT (int, -1);
}

static int rhsfn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_rhsfn);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(ydot);

    r = Int_val (caml_callbackN(*call_rhsfn,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, r);
    CAMLlocalN(args, 3);

    value *backref = user_data;
    intnat nroots;

    /* The length of gout is only available at the nroots field of the session
     * structure, so a dereference of the backref is unavoidable.  Therefore,
     * we do all of the setup here and directly call the user-supplied OCaml
     * function without going through an OCaml trampoline.  */
    WEAK_DEREF (session, *backref);

    nroots = CVODE_NROOTS_FROM_ML (session);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);
    args[2] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    r = caml_callbackN_exn (CVODE_ROOTSFN_FROM_ML (session),
			    sizeof (args) / sizeof (*args),
			    args);

    CAMLreturnT(int, check_exception(session, r));
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    int r;
    value *backref = user_data;
    CAML_FN (call_errw);

    args[0] = *backref;
    args[1] = NVEC_BACKLINK (y);
    args[2] = NVEC_BACKLINK (ewt);

    r = Int_val (caml_callbackN (*call_errw,
				 sizeof (args) / sizeof (*args),
				 args));

    CAMLreturnT (int, r);
}

static value make_jac_arg(realtype t, N_Vector y, N_Vector fy, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_CVODE_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_FY, NVEC_BACKLINK(fy));
    Store_field(r, RECORD_CVODE_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

static value make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(3);
    Store_field(r, 0, NVEC_BACKLINK(tmp1));
    Store_field(r, 1, NVEC_BACKLINK(tmp2));
    Store_field(r, 2, NVEC_BACKLINK(tmp3));
    CAMLreturn(r);
}

/* Dense and band Jacobians only work with serial NVectors.  */
static int jacfn(
	long int n,
	realtype t,
	N_Vector y,
	N_Vector fy,	     
	DlsMat Jac,
	void *user_data,
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN (args, 3);
    int r;
    value *backref = user_data;
    CAML_FN (call_jacfn);

    args[0] = *backref;
    args[1] = make_jac_arg (t, y, fy, make_triple_tmp (tmp1, tmp2, tmp3));
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
	realtype t,
	N_Vector y,
	N_Vector fy,
	DlsMat Jac,
	void *user_data, 	 
	N_Vector tmp1,
	N_Vector tmp2,
	N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_bandjacfn);

    args[0] = *backref;
    args[1] = caml_alloc_tuple(RECORD_CVODE_BANDRANGE_SIZE);
    Store_field(args[1], RECORD_CVODE_BANDRANGE_MUPPER, Val_long(mupper));
    Store_field(args[1], RECORD_CVODE_BANDRANGE_MLOWER, Val_long(mlower));
    args[2] = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));
    args[3] = caml_alloc_final(2, NULL, 0, 1);
    DLSMAT(args[3]) = Jac;

    r = Int_val (caml_callbackN(*call_bandjacfn,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int presetupfn(
    realtype t,
    N_Vector y,
    N_Vector fy,
    booleantype jok,
    booleantype *jcurPtr,
    realtype gamma,
    void *user_data,
    N_Vector tmp1,
    N_Vector tmp2,
    N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal2(session, r);
    CAMLlocalN(args, 3);
    value *backref = user_data;

    /* The presetup function must return a boolean (in addition to possible
     * exceptions), so, we do all of the setup here and directly call the
     * user-supplied OCaml function without going through an OCaml
     * trampoline.  */
    WEAK_DEREF (session, *backref);

    args[0] = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));
    args[1] = Val_bool(jok);
    args[2] = caml_copy_double(gamma);

    r = caml_callbackN_exn(CVODE_PRESETUPFN_FROM_ML (session),
                           sizeof (args) / sizeof (*args),
                           args);

    if (!Is_exception_result(r)) {
	*jcurPtr = Bool_val(r);
    }

    CAMLreturnT(int, check_exception(session, r));
}

static value make_spils_solve_arg(
	N_Vector r,
	realtype gamma,
	realtype delta,
	int lr)

{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(RECORD_CVODE_SPILS_SOLVE_ARG_SIZE);
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_RHS, NVEC_BACKLINK(r));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_GAMMA,
                caml_copy_double(gamma));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_DELTA,
                caml_copy_double(delta));
    Store_field(v, RECORD_CVODE_SPILS_SOLVE_ARG_LEFT,
                lr == 1 ? Val_true : Val_false);

    CAMLreturn(v);
}

static int presolvefn(
	realtype t,
	N_Vector y,
	N_Vector fy,
	N_Vector r,
	N_Vector z,
	realtype gamma,
	realtype delta,
	int lr,
	void *user_data,
	N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal1(rv);
    CAMLlocalN(args, 4);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_presolvefn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, fy, NVEC_BACKLINK(tmp));
    args[2] = make_spils_solve_arg(r, gamma, delta, lr);
    args[3] = NVEC_BACKLINK(z);

    retcode = Int_val (caml_callbackN(*call_presolvefn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    CAMLreturnT(int, retcode);
}

static int jactimesfn(
    N_Vector v,
    N_Vector Jv,
    realtype t,
    N_Vector y,
    N_Vector fy,
    void *user_data,
    N_Vector tmp)
{
    CAMLparam0();
    CAMLlocal1(r);
    CAMLlocalN(args, 4);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_jactimesfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, fy, NVEC_BACKLINK(tmp));
    args[2] = NVEC_BACKLINK(v);
    args[3] = NVEC_BACKLINK(Jv);

    retcode = Int_val (caml_callbackN(*call_jactimesfn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    CAMLreturnT(int, retcode);
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim void c_cvode_dls_dense (value vcvode_mem, value vneqs, value vset_jac)
{
    CAMLparam3(vcvode_mem, vneqs, vset_jac);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    long neqs = Long_val(vneqs);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVDense (cvode_mem, neqs);
    CHECK_FLAG ("CVDense", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetDenseJacFn(CVODE_MEM_FROM_ML(vcvode_mem), jacfn);
	CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    }
    CAMLreturn0;
}

CAMLprim void c_cvode_dls_lapack_dense (value vcvode_mem, value vneqs,
					value vset_jac)
{
    CAMLparam3 (vcvode_mem, vneqs, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVLapackDense (cvode_mem, neqs);
    CHECK_FLAG ("CVLapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetDenseJacFn (cvode_mem, jacfn);
	CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn0;
}

CAMLprim void c_cvode_dls_set_dense_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetDenseJacFn(CVODE_MEM_FROM_ML(vdata), jacfn);
    CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_dls_clear_dense_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetDenseJacFn(CVODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_dls_band (value vcvode_mem_neqs,
				value vmupper, value vmlower,
				value vset_jac)
{
    CAMLparam4(vcvode_mem_neqs, vmupper, vmlower, vset_jac);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vcvode_mem_neqs, 0));
    long neqs = Long_val (Field(vcvode_mem_neqs, 1));
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVBand (cvode_mem, neqs, Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("CVBand", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetBandJacFn(cvode_mem, bandjacfn);
	CHECK_FLAG("CVDlsSetBandJacFn", flag);
    }
    CAMLreturn0;
}

CAMLprim void c_cvode_dls_lapack_band (value vcvode_mem_neqs, value vmupper,
					value vmlower, value vset_jac)
{
    CAMLparam4(vcvode_mem_neqs, vmupper, vmlower, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vcvode_mem_neqs, 0));
    long neqs = Long_val(Field(vcvode_mem_neqs, 1));
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVLapackBand (cvode_mem, neqs,
			 Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("CVLapackBand", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetBandJacFn(cvode_mem, bandjacfn);
	CHECK_FLAG("CVDlsSetBandJacFn", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn0;
}

CAMLprim void c_cvode_dls_set_band_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetBandJacFn(CVODE_MEM_FROM_ML(vdata), bandjacfn);
    CHECK_FLAG("CVDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_dls_clear_band_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVDlsSetBandJacFn(CVODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("CVDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_spils_banded_spgmr (value vcvode_mem_neqs,
					  value vmupper, value vmlower,
					  value vmaxl, value vtype)
{
    CAMLparam5 (vcvode_mem_neqs, vmupper, vmlower, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vcvode_mem_neqs, 0));
    long neqs = Long_val(Field(vcvode_mem_neqs, 1));
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpgmr (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpgmr", flag);
    flag = CVBandPrecInit (cvode_mem, neqs,
			   Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("CVBandPrecInit", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_spils_banded_spbcg (value vcvode_mem_neqs,
					  value vmupper, value vmlower,
					  value vmaxl, value vtype)
{
    CAMLparam5 (vcvode_mem_neqs, vmupper, vmlower, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vcvode_mem_neqs, 0));
    long neqs = Long_val(Field(vcvode_mem_neqs, 1));
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpbcg (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpbcg", flag);
    flag = CVBandPrecInit (cvode_mem, neqs,
			   Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("CVBandPrecInit", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_spils_banded_sptfqmr (value vcvode_mem_neqs,
					    value vmupper, value vmlower,
					    value vmaxl, value vtype)
{
    CAMLparam5 (vcvode_mem_neqs, vmupper, vmlower, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vcvode_mem_neqs, 0));
    long neqs = Long_val(Field(vcvode_mem_neqs, 1));
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSptfqmr (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSptfqmr", flag);
    flag = CVBandPrecInit (cvode_mem, neqs,
			   Long_val (vmupper), Long_val (vmlower));
    CHECK_FLAG ("CVBandPrecInit", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_preconditioner (value vsession,
						value vset_presetup,
						value vset_jac)
{
    CAMLparam3 (vsession, vset_presetup, vset_jac);
    int flag;
    void *mem = CVODE_MEM_FROM_ML (vsession);
    CVSpilsPrecSetupFn setup = Bool_val (vset_presetup) ? presetupfn : NULL;

    flag = CVSpilsSetPreconditioner (mem, setup, presolvefn);
    CHECK_FLAG ("CVSpilsSetPreconditioner", flag);
    if (Bool_val (vset_jac)) {
	flag = CVSpilsSetJacTimesVecFn (mem, jactimesfn);
	CHECK_FLAG ("CVSpilsSetJacTimesVecFn", flag);
    }

    CAMLreturn0;
}

CAMLprim void c_cvode_wf_tolerances (value vdata)
{
    CAMLparam1(vdata);
 
    int flag = CVodeWFtolerances(CVODE_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("CVodeWFtolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_preconditioner(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVSpilsSetPreconditioner(CVODE_MEM_FROM_ML(vdata),
					presetupfn, presolvefn);
    CHECK_FLAG("CVSpilsSetPreconditioner", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_set_jac_times_vec_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVSpilsSetJacTimesVecFn(CVODE_MEM_FROM_ML(vdata), jactimesfn);
    CHECK_FLAG("CVSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_clear_jac_times_vec_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = CVSpilsSetJacTimesVecFn(CVODE_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("CVSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

/* basic interface */

/* CVodeCreate() + CVodeInit().  */
CAMLprim value c_cvode_init(value weakref, value lmm, value iter, value initial,
			    value t0)
{
    CAMLparam5(weakref, lmm, iter, initial, t0);
    CAMLlocal1(r);

    if (sizeof(int) != 4) {
	caml_failwith("The library assumes that an int (in C) has 32-bits.");
    }

    int flag;

    int lmm_c;
    switch (Int_val(lmm)) {
    case VARIANT_CVODE_LMM_ADAMS:
	lmm_c = CV_ADAMS;
	break;

    case VARIANT_CVODE_LMM_BDF:
	lmm_c = CV_BDF;
	break;

    default:
	caml_failwith("Illegal lmm value.");
	break;
    }

    int iter_c;
    if (Is_block(iter)) {
	iter_c = CV_NEWTON;
    } else {
	iter_c = CV_FUNCTIONAL;
    }

    void *cvode_mem = CVodeCreate(lmm_c, iter_c);
    if (cvode_mem == NULL)
	caml_failwith("CVodeCreate returned NULL");

    N_Vector initial_nv = NVEC_VAL(initial);
    flag = CVodeInit(cvode_mem, rhsfn, Double_val(t0), initial_nv);
    if (flag != CV_SUCCESS) {
	CVodeFree (cvode_mem);
	CHECK_FLAG("CVodeInit", flag);
    }

    value *backref;
    backref = malloc (sizeof (*backref));
    if (backref == NULL) {
	CVodeFree (cvode_mem);
	caml_raise_out_of_memory();
    }
    *backref = weakref;
    caml_register_generational_global_root (backref);
    CVodeSetUserData (cvode_mem, backref);

    r = caml_alloc_tuple (3);
    Store_field (r, 0, (value)cvode_mem);
    Store_field (r, 1, (value)backref);
    Store_field (r, 2, Val_long (0)); /* no err_file = NULL */

    CAMLreturn(r);
}

/* Set the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim void c_cvode_root_init (value vdata, value vnroots)
{
    CAMLparam2 (vdata, vnroots);
    void *cvode_mem = CVODE_MEM_FROM_ML (vdata);
    int nroots = Int_val (vnroots);
    int flag = CVodeRootInit (cvode_mem, nroots, roots);
    CHECK_FLAG ("CVodeRootInit", flag);
    Store_field (vdata, RECORD_CVODE_SESSION_NROOTS, vnroots);
    CAMLreturn0;
}

CAMLprim void c_cvode_sv_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = CVodeSVtolerances(CVODE_MEM_FROM_ML(vdata),
				 Double_val(reltol), atol_nv);
    CHECK_FLAG("CVodeSVtolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_reinit(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    // TODO:
#if SAFETY_CHECKS && CVODE_ML_BIGARRAYS
    /* This can't be checked for generic nvectors.  */
    if (CVODE_NEQS_FROM_ML (vdata) != Caml_ba_array_val(y0)->dim[0])
	caml_invalid_argument ("Cvode.reinit: y vector has incorrect length");
#endif

    N_Vector y0_nv = NVEC_VAL(y0);
    int flag = CVodeReInit(CVODE_MEM_FROM_ML(vdata), Double_val(t0), y0_nv);
    CHECK_FLAG("CVodeReInit", flag);

    CAMLreturn0;
}

static value solver(value vdata, value nextt, value vy, int onestep)
{
    CAMLparam3(vdata, nextt, vy);
    CAMLlocal1(ret);
    realtype tret;
    int flag;
    N_Vector y;
    enum cvode_solver_result_tag result;

    // TODO:
#if SAFETY_CHECKS && CVODE_ML_BIGARRAYS
    /* This can't be checked for generic nvectors.  */
    if (CVODE_NEQS_FROM_ML (vdata) != Caml_ba_array_val(vy)->dim[0])
	caml_invalid_argument ("Cvode.solve: y vector has incorrect length");
#endif

    y = NVEC_VAL (vy);
    // Caml_ba_data_val(y) must not be shifted by the OCaml GC during this
    // function call, which calls Caml through the callback f.  Is this
    // guaranteed?
    flag = CVode (CVODE_MEM_FROM_ML (vdata), Double_val (nextt), y, &tret,
		  onestep ? CV_ONE_STEP : CV_NORMAL);

    switch (flag) {
    case CV_SUCCESS:
	result = VARIANT_CVODE_SOLVER_RESULT_CONTINUE;
	break;

    case CV_ROOT_RETURN:
	result = VARIANT_CVODE_SOLVER_RESULT_ROOTSFOUND;
	break;

    case CV_TSTOP_RETURN:
	result = VARIANT_CVODE_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_CVODE_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_CVODE_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	CHECK_FLAG ("CVode", flag);
    }

    /* Hmm...should this go in the production code or not?  */
    if (Is_block (Field (vdata, RECORD_CVODE_SESSION_EXN_TEMP)))
	abort ();

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}

CAMLprim value c_cvode_solve_normal(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value c_cvode_solve_one_step(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 1));
}

CAMLprim void c_cvode_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVEC_VAL(vy);

    int flag = CVodeGetDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
			   Int_val(vk), y_nv);
    CHECK_FLAG("CVodeGetDky", flag);
    
    CAMLreturn0;
}

CAMLprim void c_cvode_get_err_weights(value vcvode_mem, value verrws)
{
    CAMLparam2(vcvode_mem, verrws);

    N_Vector errws_nv = NVEC_VAL(verrws);

    int flag = CVodeGetErrWeights(CVODE_MEM_FROM_ML(vcvode_mem), errws_nv);
    CHECK_FLAG("CVodeGetErrWeights", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_get_est_local_errors(value vcvode_mem, value vele)
{
    CAMLparam2(vcvode_mem, vele);

    N_Vector ele_nv = NVEC_VAL(vele);

    int flag = CVodeGetEstLocalErrors(CVODE_MEM_FROM_ML(vcvode_mem), ele_nv);
    CHECK_FLAG("CVodeGetEstLocalErrors", flag);

    CAMLreturn0;
}


void cvode_ml_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	    || flag == CV_ROOT_RETURN
	    || flag == CV_TSTOP_RETURN) return;

    switch (flag) {
	case CV_ILL_INPUT:
	    caml_raise_constant(*caml_named_value("cvode_IllInput"));

	case CV_TOO_CLOSE:
	    caml_raise_constant(*caml_named_value("cvode_TooClose"));

	case CV_TOO_MUCH_WORK:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchWork"));

	case CV_TOO_MUCH_ACC:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchAccuracy"));

	case CV_ERR_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ErrFailure"));

	case CV_CONV_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ConvergenceFailure"));

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

	case CV_BAD_K:
	    caml_raise_constant(*caml_named_value("cvode_BadK"));

	case CV_BAD_T:
	    caml_raise_constant(*caml_named_value("cvode_BadT"));

	case CV_BAD_DKY:
	    caml_raise_constant(*caml_named_value("cvode_BadDky"));

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

/* basic interface */

CAMLprim void c_cvode_session_finalize(value vdata)
{
    if (CVODE_MEM_FROM_ML(vdata) != NULL) {
	void *cvode_mem = CVODE_MEM_FROM_ML(vdata);
	value *backref = CVODE_BACKREF_FROM_ML(vdata);
	CVodeFree(&cvode_mem);
	caml_remove_generational_global_root (backref);
	free (backref);
    }

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_CVODE_SESSION_ERRFILE));
    if (err_file != NULL) {
	fclose(err_file);
    }
}

CAMLprim void c_cvode_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("CVodeSStolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_get_root_info(value vdata, value roots)
{
    CAMLparam2(vdata, roots);

    int roots_l = Caml_ba_array_val(roots)->dim[0];
    int *roots_d = INT_ARRAY(roots);

    if (roots_l < CVODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("roots array is too short");
    }

    int flag = CVodeGetRootInfo(CVODE_MEM_FROM_ML(vdata), roots_d);
    CHECK_FLAG("CVodeGetRootInfo", flag);

    CAMLreturn0;
}

CAMLprim value c_cvode_get_integrator_stats(value vdata)
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

    flag = CVodeGetIntegratorStats(CVODE_MEM_FROM_ML(vdata),
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

    r = caml_alloc_tuple(RECORD_CVODE_INTEGRATOR_STATS_SIZE);
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_STEPS, Val_long(nsteps));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_RHS_EVALS, Val_long(nfevals));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_LINEAR_SOLVER_SETUPS, Val_long(nlinsetups));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_ERROR_TEST_FAILURES, Val_long(netfails));

    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_LAST_INTERNAL_ORDER, Val_int(qlast));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_NEXT_INTERNAL_ORDER, Val_int(qcur));

    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_INITIAL_STEP_SIZE, caml_copy_double(hinused));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_LAST_STEP_SIZE, caml_copy_double(hlast));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_NEXT_STEP_SIZE, caml_copy_double(hcur));
    Store_field(r, RECORD_CVODE_INTEGRATOR_STATS_INTERNAL_TIME, caml_copy_double(tcur));

    CAMLreturn(r);
}

CAMLprim void c_cvode_set_error_file(value vdata, value vpath, value vtrunc)
{
    CAMLparam3(vdata, vpath, vtrunc);

    FILE* err_file =
      (FILE *)Long_val(Field(vdata, RECORD_CVODE_SESSION_ERRFILE));

    if (err_file != NULL) {
	fclose(err_file);
	Store_field(vdata, RECORD_CVODE_SESSION_ERRFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    err_file = fopen(String_val(vpath), mode);
    if (err_file == NULL) {
	uerror("fopen", vpath);
    }

    int flag = CVodeSetErrFile(CVODE_MEM_FROM_ML(vdata), err_file);
    CHECK_FLAG("CVodeSetErrFile", flag);

    Store_field(vdata, RECORD_CVODE_SESSION_ERRFILE, Val_long(err_file));

    CAMLreturn0;
}

CAMLprim void c_cvode_set_functional (value vdata)
{
    CAMLparam1 (vdata);
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vdata), CV_FUNCTIONAL);
    CHECK_FLAG ("CVodeSetIterType", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_set_root_direction(value vdata, value rootdirs)
{
    CAMLparam2(vdata, rootdirs);

    int rootdirs_l = Caml_ba_array_val(rootdirs)->dim[0];
    int *rootdirs_d = INT_ARRAY(rootdirs);

    if (rootdirs_l < CVODE_NROOTS_FROM_ML(vdata)) {
	caml_invalid_argument("root directions array is too short");
    }

    int flag = CVodeSetRootDirection(CVODE_MEM_FROM_ML(vdata), rootdirs_d);
    CHECK_FLAG("CVodeSetRootDirection", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_prec_type(value vcvode_mem, value vptype)
{
    CAMLparam2(vcvode_mem, vptype);

    int flag = CVSpilsSetPrecType(CVODE_MEM_FROM_ML(vcvode_mem),
				  spils_precond_type(vptype));
    CHECK_FLAG("CVSpilsSetPrecType", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_diag (value vcvode_mem)
{
    CAMLparam1 (vcvode_mem);
    int flag = CVodeSetIterType (CVODE_MEM_FROM_ML (vcvode_mem), CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVDiag(CVODE_MEM_FROM_ML (vcvode_mem));
    CHECK_FLAG("CVDiag", flag);
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Boiler plate definitions for Sundials interface.
 *
 */

CAMLprim value c_cvode_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = CVodeGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_FLAG("CVodeGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_cvode_get_num_steps(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumSteps(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumSteps", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_num_lin_solv_setups(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumLinSolvSetups(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_num_err_test_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVodeGetNumErrTestFails(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvode_get_last_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    int v;

    flag = CVodeGetLastOrder(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetLastOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_cvode_get_current_order(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    int v;

    flag = CVodeGetCurrentOrder(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentOrder", flag);

    CAMLreturn(Val_int(v));
}

CAMLprim value c_cvode_get_actual_init_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetActualInitStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetActualInitStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_cvode_get_last_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetLastStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetLastStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_cvode_get_current_step(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetCurrentStep(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentStep", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim value c_cvode_get_current_time(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    realtype v;

    flag = CVodeGetCurrentTime(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVodeGetCurrentTime", flag);

    CAMLreturn(caml_copy_double(v));
}

CAMLprim void c_cvode_set_max_ord(value vcvode_mem, value maxord)
{
    CAMLparam2(vcvode_mem, maxord);

    int flag = CVodeSetMaxOrd(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxord));
    CHECK_FLAG("CVodeSetMaxOrd", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_num_steps(value vcvode_mem, value mxsteps)
{
    CAMLparam2(vcvode_mem, mxsteps);


    int flag = CVodeSetMaxNumSteps(CVODE_MEM_FROM_ML(vcvode_mem), Long_val(mxsteps));
    CHECK_FLAG("CVodeSetMaxNumSteps", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_hnil_warns(value vcvode_mem, value mxhnil)
{
    CAMLparam2(vcvode_mem, mxhnil);


    int flag = CVodeSetMaxHnilWarns(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(mxhnil));
    CHECK_FLAG("CVodeSetMaxHnilWarns", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_stab_lim_det(value vcvode_mem, value stldet)
{
    CAMLparam2(vcvode_mem, stldet);


    int flag = CVodeSetStabLimDet(CVODE_MEM_FROM_ML(vcvode_mem), Bool_val(stldet));
    CHECK_FLAG("CVodeSetStabLimDet", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_init_step(value vcvode_mem, value hin)
{
    CAMLparam2(vcvode_mem, hin);


    int flag = CVodeSetInitStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hin));
    CHECK_FLAG("CVodeSetInitStep", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_min_step(value vcvode_mem, value hmin)
{
    CAMLparam2(vcvode_mem, hmin);


    int flag = CVodeSetMinStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hmin));
    CHECK_FLAG("CVodeSetMinStep", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_step(value vcvode_mem, value hmax)
{
    CAMLparam2(vcvode_mem, hmax);


    int flag = CVodeSetMaxStep(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(hmax));
    CHECK_FLAG("CVodeSetMaxStep", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_stop_time(value vcvode_mem, value tstop)
{
    CAMLparam2(vcvode_mem, tstop);


    int flag = CVodeSetStopTime(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(tstop));
    CHECK_FLAG("CVodeSetStopTime", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_err_test_fails(value vcvode_mem, value maxnef)
{
    CAMLparam2(vcvode_mem, maxnef);


    int flag = CVodeSetMaxErrTestFails(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxnef));
    CHECK_FLAG("CVodeSetMaxErrTestFails", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_nonlin_iters(value vcvode_mem, value maxcor)
{
    CAMLparam2(vcvode_mem, maxcor);


    int flag = CVodeSetMaxNonlinIters(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxcor));
    CHECK_FLAG("CVodeSetMaxNonlinIters", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_max_conv_fails(value vcvode_mem, value maxncf)
{
    CAMLparam2(vcvode_mem, maxncf);


    int flag = CVodeSetMaxConvFails(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxncf));
    CHECK_FLAG("CVodeSetMaxConvFails", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_nonlin_conv_coef(value vcvode_mem, value nlscoef)
{
    CAMLparam2(vcvode_mem, nlscoef);


    int flag = CVodeSetNonlinConvCoef(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(nlscoef));
    CHECK_FLAG("CVodeSetNonlinConvCoef", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_set_no_inactive_root_warn(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag = CVodeSetNoInactiveRootWarn(CVODE_MEM_FROM_ML(vcvode_mem));
    CHECK_FLAG("CVodeSetNoInactiveRootWarn", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_gs_type(value vcvode_mem, value vgstype)
{
    CAMLparam2(vcvode_mem, vgstype);

    int flag = CVSpilsSetGSType(CVODE_MEM_FROM_ML(vcvode_mem),
				spils_gs_type(vgstype));
    CHECK_FLAG("CVSpilsSetGSType", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_eps_lin(value vcvode_mem, value eplifac)
{
    CAMLparam2(vcvode_mem, eplifac);

    int flag = CVSpilsSetEpsLin(CVODE_MEM_FROM_ML(vcvode_mem), Double_val(eplifac));
    CHECK_FLAG("CVSpilsSetEpsLin", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_set_maxl(value vcvode_mem, value maxl)
{
    CAMLparam2(vcvode_mem, maxl);

    int flag = CVSpilsSetMaxl(CVODE_MEM_FROM_ML(vcvode_mem), Int_val(maxl));
    CHECK_FLAG("CVSpilsSetMaxl", flag);

    CAMLreturn0;
}

/* statistic accessor functions */

CAMLprim value c_cvode_get_num_stab_lim_order_reds(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumStabLimOrderReds(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumStabLimOrderReds", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_get_tol_scale_factor(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    realtype r;
    int flag = CVodeGetTolScaleFactor(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetTolScaleFactor", flag);

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value c_cvode_get_num_nonlin_solv_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_get_num_nonlin_solv_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_get_nonlin_solv_stats(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int nniters, nncfails;
    int flag = CVodeGetNonlinSolvStats(CVODE_MEM_FROM_ML(vcvode_mem),
				       &nniters, &nncfails);
    CHECK_FLAG("CVodeGetNonlinSolvStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nniters));
    Store_field(r, 1, Val_long(nncfails));

    CAMLreturn(r);
}

CAMLprim value c_cvode_get_num_g_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVodeGetNumGEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVodeGetNumGEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_dls_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = CVDlsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrwLS, &leniwLS);
    CHECK_FLAG("CVDlsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value c_cvode_dls_get_num_jac_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDlsGetNumJacEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVDlsGetNumJacEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_dls_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDlsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVDlsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_diag_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwLS;
    long int leniwLS;

    int flag = CVDiagGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrwLS, &leniwLS);
    CHECK_FLAG("CVDiagGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrwLS));
    Store_field(r, 1, Val_long(leniwLS));

    CAMLreturn(r);
}

CAMLprim value c_cvode_diag_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVDiagGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVDiagGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_bandprec_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    long int lenrwBP;
    long int leniwBP;

    int flag = CVBandPrecGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrwBP, &leniwBP);
    CHECK_FLAG("CVBandPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrwBP));
    Store_field(r, 1, Val_long(leniwBP));

    CAMLreturn(r);
}

CAMLprim value c_cvode_bandprec_get_num_rhs_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVBandPrecGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVBandPrecGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

/* spils functions */

CAMLprim void c_cvode_spils_spgmr (value vcvode_mem, value vmaxl, value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpgmr (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpgmr", flag);
    CAMLreturn0;
}

CAMLprim void c_cvode_spils_spbcg (value vcvode_mem, value vmaxl, value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSpbcg (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSpbcg", flag);

    CAMLreturn0;
}

CAMLprim void c_cvode_spils_sptfqmr (value vcvode_mem, value vmaxl,
				     value vtype)
{
    CAMLparam3 (vcvode_mem, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVodeSetIterType (cvode_mem, CV_NEWTON);
    CHECK_FLAG ("CVodeSetIterType", flag);
    flag = CVSptfqmr (cvode_mem, spils_precond_type (vtype), Int_val (vmaxl));
    CHECK_FLAG ("CVSptfqmr", flag);

    CAMLreturn0;
}

CAMLprim value c_cvode_spils_get_num_lin_iters(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumLinIters(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_conv_fails(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumConvFails(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = CVSpilsGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_FLAG("CVSpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_cvode_spils_get_num_prec_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_prec_solves(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumPrecSolves(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_jtimes_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumJtimesEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_cvode_spils_get_num_rhs_evals (value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    long int r;
    int flag = CVSpilsGetNumRhsEvals(CVODE_MEM_FROM_ML(vcvode_mem), &r);
    CHECK_FLAG("CVSpilsGetNumRhsEvals", flag);

    CAMLreturn(Val_long(r));
}

