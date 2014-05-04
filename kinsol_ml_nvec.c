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

/* The parts of the Sundials interface that distinguish between Serial
   NVectors (handled by Bigarrays) and generic NVectors (handled by a
   wrapper type). */

#include <kinsol/kinsol.h>
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
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <sundials/sundials_config.h>

#include "spils_ml.h"
#include "kinsol_ml.h"
#include "nvector_ml.h"

#if SUNDIALS_BLAS_LAPACK == 1
#include <kinsol/kinsol_lapack.h>
#endif

#ifdef RESTRICT_INTERNAL_PRECISION
#ifdef __GNUC__
#include <fpu_control.h>
#endif
#endif

// Call with KINSOL_ML_BIGARRAYS to compile for the Serial NVector to
// Bigarray interface code.

#ifdef KINSOL_ML_BIGARRAYS

#define CVTYPE(fname) c_ba_kinsol_ ## fname
#include <nvector/nvector_serial.h>

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0

#define NVECTORIZE_VAL(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDVAL(nv) N_VDestroy(nv)

#else

#define CVTYPE(fname) c_nvec_kinsol_ ## fname
#include <sundials/sundials_nvector.h>

#define WRAP_NVECTOR(v) NVEC_DATA(v)
#define RELINQUISH_WRAPPEDNV(v) {}

#define NVECTORIZE_VAL(v) NVEC_VAL(v)
#define RELINQUISH_NVECTORIZEDVAL(nv) {}

#endif

#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define CVTYPESTR(fname) QUOTE(CVTYPE(fname))

/* callbacks */

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

CAMLprim void CVTYPE(set_err_handler_fn)(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = KINSetErrHandlerFn(KINSOL_MEM_FROM_ML(vdata), errh,
				  KINSOL_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("KINSetErrHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(clear_err_handler_fn)(value vdata)
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

CAMLprim void CVTYPE(set_info_handler_fn)(value vdata)
{
    CAMLparam1(vdata);
 
    int flag = KINSetInfoHandlerFn(KINSOL_MEM_FROM_ML(vdata), infoh,
				   KINSOL_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("KINSetInfoHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(clear_info_handler_fn)(value vdata)
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
    CAML_FN (call_rhsfn);

    vuu = WRAP_NVECTOR(uu);
    vval = WRAP_NVECTOR(val);

    // The data payloads inside vuu and vval are only valid during this
    // call, afterward that memory goes back to kinsol. These bigarrays must
    // not be retained by closure_rhsfn! If it wants a permanent copy, then
    // it has to make it manually.
    r = Int_val (caml_callback3(*call_rhsfn, *backref, vuu, vval));

    RELINQUISH_WRAPPEDNV(vuu);
    RELINQUISH_WRAPPEDNV(vval);

    CAMLreturnT(int, r);
}

static value make_prec_solve_arg(N_Vector uscale, N_Vector fscale)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_SIZE);
    Store_field(r, RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_USCALE,
	        WRAP_NVECTOR(uscale));
    Store_field(r, RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_FSCALE,
		WRAP_NVECTOR(fscale));

    CAMLreturn(r);
}

static void relinquish_prec_solve_arg(value arg)
{
    CAMLparam1(arg);

    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_USCALE));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_KINSOL_SPILS_PREC_SOLVE_ARG_FSCALE));

    CAMLreturn0;
}

static value make_jac_arg(N_Vector u, N_Vector fu, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_KINSOL_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_U, WRAP_NVECTOR(u));
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_FU, WRAP_NVECTOR(fu));
    Store_field(r, RECORD_KINSOL_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

static value make_double_tmp(N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, WRAP_NVECTOR(tmp1));
    Store_field(r, 1, WRAP_NVECTOR(tmp2));
    CAMLreturn(r);
}

#define DOUBLE 1
#define SINGLE 0
static void relinquish_jac_arg(value arg, int doublearg)
{
    CAMLparam1(arg);
    CAMLlocal1(tmp);

    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_KINSOL_JACOBIAN_ARG_JAC_U));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_KINSOL_JACOBIAN_ARG_JAC_FU));

    tmp = Field(arg, RECORD_KINSOL_JACOBIAN_ARG_JAC_TMP);

    if (doublearg) {
	RELINQUISH_WRAPPEDNV(Field(tmp, 0));
	RELINQUISH_WRAPPEDNV(Field(tmp, 1));
    } else {
	RELINQUISH_WRAPPEDNV(tmp);
    }

    CAMLreturn0;
}

/* Dense and band Jacobians only work with serial NVectors.  */
#if KINSOL_ML_BIGARRAYS
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
    Store_field (args[2], 1, (value)Jac);

    r = Int_val (caml_callbackN (*call_jacfn,
				 sizeof (args) / sizeof (*args),
				 args));

    relinquish_jac_arg(args[1], DOUBLE);
    // note: matrix is also invalid after the callback

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
    CAMLlocalN(args, 5);
    int r;
    value *backref = user_data;
    CAML_FN (call_bandjacfn);

    args[0] = *backref;
    args[1] = make_jac_arg(u, fu, make_double_tmp(tmp1, tmp2));
    args[2] = Val_long(mupper);
    args[3] = Val_long(mlower);
    args[4] = caml_alloc_final(2, NULL, 0, 1);
    Store_field (args[4], 1, (value)Jac);

    r = Int_val (caml_callbackN(*call_bandjacfn,
                                sizeof (args) / sizeof (*args),
                                args));

    relinquish_jac_arg(args[1], DOUBLE);
    // note: args[4] is also invalid after the callback

    CAMLreturnT(int, r);
}
#endif	/* KINSOL_ML_BIGARRAYS */

static int presetupfn(
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
    CAML_FN (call_presetupfn);

    args[0] = *backref;
    args[1] = make_jac_arg(uu, fu, make_double_tmp(tmp1, tmp2));
    args[2] = make_prec_solve_arg(uscale, fscale);

    retcode = Int_val (caml_callbackN(*call_presetupfn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    relinquish_jac_arg(args[1], DOUBLE);
    relinquish_prec_solve_arg(args[2]);

    CAMLreturnT(int, retcode);
}

static int presolvefn(
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
    CAMLlocalN(args, 3);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_presolvefn);

    args[0] = *backref;
    args[1] = make_jac_arg(uu, fu, WRAP_NVECTOR(tmp));
    args[2] = make_prec_solve_arg(uscale, fscale);

    retcode = Int_val (caml_callbackN(*call_presolvefn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    relinquish_jac_arg(args[1], SINGLE);
    relinquish_prec_solve_arg(args[2]);

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
    args[1] = WRAP_NVECTOR(v);
    args[2] = WRAP_NVECTOR(Jv);
    args[3] = WRAP_NVECTOR(u);
    args[4] = Val_bool(*new_uu);

    vr = caml_callbackN (*call_jactimesfn,
			 sizeof (args) / sizeof (*args),
			 args);
    r = Field(vr, 1);
    if (r == 0) *new_uu = Bool_val(Field(vr, 0));

    RELINQUISH_WRAPPEDNV(args[1]);
    RELINQUISH_WRAPPEDNV(args[2]);
    RELINQUISH_WRAPPEDNV(args[3]);

    CAMLreturnT(int, r);
}

/* Dense and Band can only be used with serial NVectors.  */
#ifdef KINSOL_ML_BIGARRAYS
CAMLprim void CVTYPE(dls_dense) (value vkin_mem, value vset_jac)
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

CAMLprim void CVTYPE(dls_lapack_dense) (value vkin_mem, value vset_jac)
{
    CAMLparam2 (vkin_mem, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    long neqs = KINSOL_NEQS_FROM_ML (vkin_mem);
    int flag;

    flag = KINLapackDense (kin_mem, neqs);
    CHECK_FLAG ("KINCVLapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = KINDlsSetDenseJacFn (KINSOL_MEM_FROM_ML (vkin_mem), jacfn);
	CHECK_FLAG("KINDlsSetDenseJacFn", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_set_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetDenseJacFn(KINSOL_MEM_FROM_ML(vdata), jacfn);
    CHECK_FLAG("KINDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_clear_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetDenseJacFn(KINSOL_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("KINDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_band) (value vkin_mem,
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

CAMLprim void CVTYPE(dls_lapack_band) (value vkin_mem, value vmupper,
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

CAMLprim void CVTYPE(dls_set_band_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vdata), bandjacfn);
    CHECK_FLAG("KINDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(dls_clear_band_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINDlsSetBandJacFn(KINSOL_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("KINDlsSetBandJacFn", flag);
    CAMLreturn0;
}
#endif	/* KINSOL_ML_BIGARRAYS */

CAMLprim void CVTYPE(spils_set_preconditioner) (value vsession,
						value vset_presetup)
{
    CAMLparam2 (vsession, vset_presetup);
    int flag;
    void *mem = KINSOL_MEM_FROM_ML (vsession);
    KINSpilsPrecSetupFn setup = Bool_val (vset_presetup) ? presetupfn : NULL;

    flag = KINSpilsSetPreconditioner (mem, setup, presolvefn);
    CHECK_FLAG ("KINCVSpilsSetPreconditioner", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(spils_set_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINSpilsSetJacTimesVecFn(KINSOL_MEM_FROM_ML(vdata), jactimesfn);
    CHECK_FLAG("KINSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

CAMLprim void CVTYPE(spils_clear_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = KINSpilsSetJacTimesVecFn(KINSOL_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("KINSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}


#ifdef KINSOL_ML_BIGARRAYS
CAMLprim void CVTYPE(set_no_res_mon)(value vkin_mem, value vnonniresmon)
{
    CAMLparam2(vkin_mem, vnonniresmon);

    int flag = KINSetNoResMon(KINSOL_MEM_FROM_ML(vkin_mem),
			      Bool_val(vnonniresmon));
    CHECK_FLAG("KINSetNoResMon", flag);

    CAMLreturn0;
}

CAMLprim void CVTYPE(set_max_sub_setup_calls)(value vkin_mem, value vmsbsetsub)
{
    CAMLparam2(vkin_mem, vmsbsetsub);

    int flag = KINSetMaxSubSetupCalls(KINSOL_MEM_FROM_ML(vkin_mem),
			              Long_val(vmsbsetsub));
    CHECK_FLAG("KINSetMaxSubSetupCalls", flag);

    CAMLreturn0;
}
#endif	/* KINSOL_ML_BIGARRAYS */

CAMLprim void CVTYPE(set_constraints)(value vkin_mem, value vconstraints)
{
    CAMLparam2(vkin_mem, vconstraints);
    int flag;
    N_Vector constraints;

#if SAFETY_CHECKS && KINSOL_ML_BIGARRAYS
    int neqs = KINSOL_NEQS_FROM_ML (vkin_mem);

    if (neqs != Caml_ba_array_val(vconstraints)->dim[0])
	caml_invalid_argument (
	    "Kinsol.set_constraints: constraints vector has incorrect length");
#endif

    constraints = NVECTORIZE_VAL(vconstraints);
    flag = KINSetConstraints(KINSOL_MEM_FROM_ML(vkin_mem), constraints);
    RELINQUISH_NVECTORIZEDVAL(constraints);

    CHECK_FLAG("KINSetConstraints", flag);

    CAMLreturn0;
}

/* basic interface */

/* KINCreate() + KINInit().  */
CAMLprim value CVTYPE(init)(value weakref, value vtemp)
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

#ifdef RESTRICT_INTERNAL_PRECISION
#ifdef __GNUC__
    fpu_control_t fpu_cw;
    _FPU_GETCW(fpu_cw);
    fpu_cw = (fpu_cw & ~_FPU_EXTENDED & ~_FPU_SINGLE) | _FPU_DOUBLE;
    _FPU_SETCW(fpu_cw);
#endif
#endif

    kin_mem = KINCreate();
    if (kin_mem == NULL)
	caml_failwith("KINCreate returned NULL");

    temp = NVECTORIZE_VAL(vtemp);
    flag = KINInit(kin_mem, sysfn, temp);
    RELINQUISH_NVECTORIZEDVAL(temp);
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

CAMLprim value CVTYPE(solve)(value vdata, value vu, value vlinesearch,
	 		     value vuscale, value vfscale)
{
    CAMLparam5(vdata, vu, vlinesearch, vuscale, vfscale);
    CAMLlocal1(ret);
    int flag;
    enum kinsol_result_tag result;
    N_Vector u;
    N_Vector uscale;
    N_Vector fscale;

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

    u = NVECTORIZE_VAL (vu);
    uscale = NVECTORIZE_VAL (vuscale);
    fscale = NVECTORIZE_VAL (vfscale);
    flag = KINSol(KINSOL_MEM_FROM_ML(vdata),
		  u,
		  Bool_val (vlinesearch) ? KIN_LINESEARCH : KIN_NONE,
		  uscale,
		  fscale);
    RELINQUISH_NVECTORIZEDVAL (u);
    RELINQUISH_NVECTORIZEDVAL (uscale);
    RELINQUISH_NVECTORIZEDVAL (fscale);

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

