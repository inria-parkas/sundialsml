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

/* The parts of the Sundials IDA interface that distinguish between Serial
   NVectors (handled by Bigarrays) and generic NVectors (handled by a wrapper
   type). */

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/fail.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/bigarray.h>

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <ida/ida_band.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_lapack.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_sptfqmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_lapack.h>
#include <sundials/sundials_config.h>

#include <nvector/nvector_serial.h>

#include "ida_ml.h"
#include "nvector_ml.h"

// Call with IDA_ML_BIGARRAYS to compile for the Serial NVector to
// Bigarray interface code.

#ifdef IDA_ML_BIGARRAYS

#define IDATYPE(fname) c_ba_ida_ ## fname
#include <nvector/nvector_serial.h>

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0

#define NVECTORIZE_VAL(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDVAL(nv) N_VDestroy(nv)

#else

#define IDATYPE(fname) c_nvec_ida_ ## fname
#include <sundials/sundials_nvector.h>

#define WRAP_NVECTOR(v) NVEC_DATA(v)
#define RELINQUISH_WRAPPEDNV(v) {}

#define NVECTORIZE_VAL(v) NVEC_VAL(v)
#define RELINQUISH_NVECTORIZEDVAL(nv) {}

#endif

#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define IDATYPESTR(fname) QUOTE(IDATYPE(fname))

#define CAML_FN(name)					\
    static value *name;					\
    if (name == NULL)					\
	name = caml_named_value (IDATYPESTR (name));

#define Val_none (Val_int(0))

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

    a = caml_alloc_tuple(RECORD_IDA_ERROR_DETAILS_SIZE);
    Store_field(a, RECORD_IDA_ERROR_DETAILS_ERROR_CODE,
		Val_int(error_code));
    Store_field(a, RECORD_IDA_ERROR_DETAILS_MODULE_NAME,
		caml_copy_string(module));
    Store_field(a, RECORD_IDA_ERROR_DETAILS_FUNCTION_NAME,
		caml_copy_string(func));
    Store_field(a, RECORD_IDA_ERROR_DETAILS_ERROR_MESSAGE,
		caml_copy_string(msg));

    caml_callback2_exn (*call_errh, *backref, a);

    CAMLreturn0;
}

CAMLprim void IDATYPE(set_err_handler_fn)(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), errh, NULL);
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(clear_err_handler_fn)(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn0;
}

static int resfn (realtype t, N_Vector y, N_Vector yp,
		  N_Vector resval, void *user_data)
{
    CAMLparam0 ();
    CAMLlocalN (args, 5);
    int r;
    value *backref = (value *)user_data;
    CAML_FN (call_resfn);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = WRAP_NVECTOR (y);
    args[3] = WRAP_NVECTOR (yp);
    args[4] = WRAP_NVECTOR (resval);

    r = Int_val (caml_callbackN (*call_resfn,
				 sizeof (args) / sizeof (*args),
				 args));

    RELINQUISH_WRAPPEDNV (args[2]);
    RELINQUISH_WRAPPEDNV (args[3]);
    RELINQUISH_WRAPPEDNV (args[4]);

    CAMLreturnT (int, r);
}

static value make_jac_arg(realtype t, realtype coef, N_Vector y, N_Vector yp,
			  N_Vector res, value tmp)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_IDA_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_COEF, caml_copy_double(coef));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_Y, WRAP_NVECTOR(y));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_YP, WRAP_NVECTOR(yp));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_RES, WRAP_NVECTOR(res));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_TMP, tmp);

    CAMLreturn(r);
}

static value make_triple_tmp(N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(3);
    Store_field(r, 0, WRAP_NVECTOR(tmp1));
    Store_field(r, 1, WRAP_NVECTOR(tmp2));
    Store_field(r, 2, WRAP_NVECTOR(tmp3));
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

static void relinquish_jac_arg(value arg, int tmp_size)
{
    CAMLparam0();
    CAMLlocal1(tmp);

    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_Y));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_YP));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_RES));

    tmp = Field(arg, RECORD_IDA_JACOBIAN_ARG_JAC_TMP);

    switch (tmp_size) {
    case 1:
	RELINQUISH_WRAPPEDNV (tmp);
	break;
    case 3:
	RELINQUISH_WRAPPEDNV (Field (tmp, 2));
	/*FALLTHROUGH*/
    case 2:
	RELINQUISH_WRAPPEDNV (Field (tmp, 1));
	RELINQUISH_WRAPPEDNV (Field (tmp, 0));
	break;
    default:
	abort ();
    }

    CAMLreturn0;
}

static int jacfn (long int neq, realtype t, realtype coef,
		  N_Vector y, N_Vector yp, N_Vector res,
		  DlsMat jac, void *user_data,
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 3);
    int r;
    value *backref = (value *)user_data;
    CAML_FN (call_jacfn);

    args[0] = *backref;
    args[1] = make_jac_arg (t, coef, y, yp, res,
			    make_triple_tmp (tmp1, tmp2, tmp3));
    args[2] = caml_alloc_final (2, NULL, 0, 1);
    Store_field (args[2], 1, (value)jac);

    r = Int_val (caml_callbackN (*call_jacfn,
				 sizeof (args) / sizeof (*args),
				 args));

    relinquish_jac_arg (args[1], 3);

    CAMLreturnT (int, r);
}

static int bandjacfn (long int neq, long int mupper, long int mlower,
		      realtype t, realtype coef, N_Vector y, N_Vector yp,
		      N_Vector res, DlsMat jac, void *user_data,
		      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 5);
    int r;
    value *backref = (value *)user_data;
    CAML_FN (call_bandjacfn);

    args[0] = *backref;
    args[1] = make_jac_arg (t, coef, y, yp, res,
			    make_triple_tmp (tmp1, tmp2, tmp3));
    args[2] = Val_int (mupper);
    args[3] = Val_int (mlower);
    args[4] = caml_alloc_final (2, NULL, 0, 1);
    Store_field (args[4], 1, (value)jac);

    r = Int_val (caml_callbackN (*call_bandjacfn,
				 sizeof (args) / sizeof (*args),
				 args));

    relinquish_jac_arg (args[1], 3);

    CAMLreturnT (int, r);
}

CAMLprim void IDATYPE (dump) (value v)
{
    CAMLparam1 (v);
    fprintf (stderr, "#- %p ", (void*)v);
    if (Is_block (v)) {
	int i, n;
	n = Wosize_val (v);
	if (n > RECORD_IDA_SESSION_SIZE) n = RECORD_IDA_SESSION_SIZE;
	fprintf (stderr, " -> [");
	for (i = 0; i < n - 1; ++i)
	    fprintf (stderr, "%p, ", (void*)Field (v, i));
	if (i < n)
	    fprintf (stderr, "%p", (void*)Field (v, i));
	if (++i < Wosize_val (v))
	    fprintf (stderr, ", ...");
	fprintf (stderr, "]");
    }
    fprintf (stderr, "-#\n");
    fflush (stderr);
    CAMLreturn0;
}

CAMLprim value caml_weak_get (value ar, value n);
static int rootsfn (realtype t, N_Vector y, N_Vector yp,
		    realtype *gout, void *user_data)
{
    CAMLparam0 ();
    CAMLlocalN (args, 5);
    value *backref = (value *)user_data;
    intnat nroots = NV_LENGTH_S (y);
    int r;
    CAML_FN (call_rootsfn);

    args[0] = *backref;
    args[1] = caml_copy_double (t);
    args[2] = WRAP_NVECTOR (y);
    args[3] = WRAP_NVECTOR (yp);
    args[4] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    r = Int_val (caml_callbackN (*call_rootsfn,
				 sizeof (args) / sizeof (*args),
				 args));

    RELINQUISH_WRAPPEDNV (args[2]);
    RELINQUISH_WRAPPEDNV (args[3]);

    CAMLreturnT (int, r);
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    int r;
    value *backref = user_data;
    CAML_FN (call_errw);

    args[0] = *backref;
    args[1] = WRAP_NVECTOR(y);
    args[2] = WRAP_NVECTOR(ewt);

    r = Int_val (caml_callbackN (*call_errw,
				 sizeof (args) / sizeof (*args),
				 args));

    RELINQUISH_WRAPPEDNV(args[1]);
    RELINQUISH_WRAPPEDNV(args[2]);

    CAMLreturnT(int, r);
}

static int presetupfn(
    realtype t,
    N_Vector y,
    N_Vector yp,
    N_Vector res,
    realtype cj,
    void *user_data,
    N_Vector tmp1,
    N_Vector tmp2,
    N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    int r;
    value *backref = user_data;
    CAML_FN (call_presetupfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, cj, y, yp, res,
			   make_triple_tmp(tmp1, tmp2, tmp3));
    r = Int_val (caml_callbackN (*call_presetupfn,
				 sizeof (args) / sizeof (*args),
				 args));

    relinquish_jac_arg(args[1], 3);

    CAMLreturnT(int, r);
}

static int presolvefn(
	realtype t,
	N_Vector y,
	N_Vector yp,
	N_Vector res,
	N_Vector r,
	N_Vector z,
	realtype cj,
	realtype delta,
	void *user_data,
	N_Vector tmp)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    value *backref = user_data;
    int rv;
    CAML_FN (call_presolvefn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, cj, y, yp, res, WRAP_NVECTOR(tmp));
    args[2] = WRAP_NVECTOR(r);
    args[3] = WRAP_NVECTOR(z);
    args[4] = caml_copy_double (delta);

    rv = Int_val (caml_callbackN (*call_presolvefn,
				  sizeof (args) / sizeof (*args),
				  args));

    relinquish_jac_arg(args[1], 1);
    RELINQUISH_WRAPPEDNV(args[2]);
    RELINQUISH_WRAPPEDNV(args[3]);

    CAMLreturnT (int, rv);
}

static int jactimesfn(
    realtype t,
    N_Vector y,
    N_Vector yp,
    N_Vector res,
    N_Vector v,
    N_Vector Jv,
    realtype cj,
    void *user_data,
    N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_jactimesfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, cj, y, yp, res, make_double_tmp (tmp1, tmp2));
    args[2] = WRAP_NVECTOR(v);
    args[3] = WRAP_NVECTOR(Jv);

    r = Int_val (caml_callbackN (*call_jactimesfn,
				 sizeof (args) / sizeof (*args),
				 args));

    relinquish_jac_arg(args[1], 2);
    RELINQUISH_WRAPPEDNV(args[2]);
    RELINQUISH_WRAPPEDNV(args[3]);

    CAMLreturnT (int, r);
}



CAMLprim void IDATYPE(wf_tolerances)(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDAWFtolerances(IDA_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("IDAWFtolerances", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_set_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetDenseJacFn(IDA_MEM_FROM_ML(vdata), jacfn);
    CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_clear_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetDenseJacFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_set_band_jac_fn)(value vdata, value fbandjacfn)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vdata), bandjacfn);
    CHECK_FLAG("IDADlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(dls_clear_band_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDADlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(set_preconditioner)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDASpilsSetPreconditioner(IDA_MEM_FROM_ML(vdata),
					 presetupfn, presolvefn);
    CHECK_FLAG("IDASpilsSetPreconditioner", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(set_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDASpilsSetJacTimesVecFn(IDA_MEM_FROM_ML(vdata), jactimesfn);
    CHECK_FLAG("IDASpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

CAMLprim void IDATYPE(clear_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDASpilsSetJacTimesVecFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDASpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

/* Sets the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim void IDATYPE(root_init) (value vida_mem, value vnroots)
{
    CAMLparam2 (vida_mem, vnroots);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int nroots = Int_val (vnroots);
    int flag;
    flag = IDARootInit (ida_mem, nroots, rootsfn);
    CHECK_FLAG ("IDARootInit", flag);
    Store_field (vida_mem, RECORD_IDA_SESSION_NROOTS, vnroots);
    CAMLreturn0;
}

/* IDACreate + IDAInit + IDASetUser.  The vy and vyp vectors must have the same
 * length (not checked).  Returns an IDAMem (= void*) along with a global root
 * wrapping weakref that is attached to the IDAMem.  Before doing the remaining
 * initialization, these pointers have to be wrapped in a session so that if
 * initialization fails the GC frees IDAMem and the weak global root.  Doing
 * that cleanup in C would be ugly, e.g. we wouldn't be able to reuse
 * c_ida_set_linear_solver() so we have to duplicate it or hack some ad-hoc
 * extensions to that function.  */
CAMLprim value IDATYPE (init) (value weakref, value vt0, value vy, value vyp)
{
    CAMLparam4 (weakref, vy, vyp, vt0);
    CAMLlocal1 (r);
    int flag;
    N_Vector y, yp;
    void *ida_mem = NULL;
    value *backref = NULL;

    ida_mem = IDACreate ();
    if (ida_mem == NULL)
	caml_failwith ("IDACreate failed");

    y = NVECTORIZE_VAL (vy);
    yp = NVECTORIZE_VAL (vyp);
    flag = IDAInit (ida_mem, resfn, Double_val (vt0), y, yp);
    RELINQUISH_NVECTORIZEDVAL (y);
    RELINQUISH_NVECTORIZEDVAL (yp);
    if (flag != IDA_SUCCESS) {
	IDAFree (ida_mem);
	CHECK_FLAG ("IDAInit", flag);
    }

    /* NB: the user data must be set before specifying the linear solver, so
     * don't move this code down.  */
    backref = malloc (sizeof (*backref));
    if (backref == NULL) {
	IDAFree (&ida_mem);
	caml_failwith ("Out of memory");
    }
    *backref = weakref;
    caml_register_generational_global_root (backref);
    IDASetUserData (ida_mem, backref);

    r = caml_alloc_tuple(3);
    Store_field(r, 0, (value)ida_mem);
    Store_field(r, 1, (value)backref);
    Store_field(r, 2, (value)NULL); // no err_file = NULL

    CAMLreturn (r);
}

CAMLprim void IDATYPE(sv_tolerances) (value ida_mem, value vrtol, value vavtol)
{
    CAMLparam3 (ida_mem, vrtol, vavtol);
    N_Vector avtol;
    int flag;

    avtol = NVECTORIZE_VAL (vavtol);
    flag = IDASVtolerances (IDA_MEM_FROM_ML (ida_mem),
			    Double_val (vrtol), avtol);
    RELINQUISH_NVECTORIZEDVAL (avtol);
    CHECK_FLAG ("IDASVtolerances", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(reinit)(value vdata, value t0, value y0, value yp0)
{
    CAMLparam4(vdata, t0, y0, yp0);

    N_Vector y0_nv = NVECTORIZE_VAL(y0);
    N_Vector yp0_nv = NVECTORIZE_VAL(yp0);
    int flag = IDAReInit(IDA_MEM_FROM_ML(vdata), Double_val(t0),
			 y0_nv, yp0_nv);
    RELINQUISH_NVECTORIZEDVAL(yp0_nv);
    RELINQUISH_NVECTORIZEDVAL(y0_nv);
    CHECK_FLAG("IDAReInit", flag);

    CAMLreturn0;
}

static value solve (value vdata, value nextt, value vy, value vyp, int onestep)
{
    CAMLparam4 (vdata, nextt, vy, vyp);
    CAMLlocal1 (ret);
    void *ida_mem = IDA_MEM_FROM_ML (vdata);
    realtype tret;
    int flag;
    N_Vector y, yp;
    enum ida_solver_result_tag result;

    y = NVECTORIZE_VAL (vy);
    yp = NVECTORIZE_VAL (vyp);
    flag = IDASolve (ida_mem, Double_val (nextt), &tret, y, yp,
	             onestep ? IDA_ONE_STEP : IDA_NORMAL);
    RELINQUISH_NVECTORIZEDVAL (y);
    RELINQUISH_NVECTORIZEDVAL (yp);

    switch (flag) {
    case IDA_SUCCESS:
	result = VARIANT_IDA_SOLVER_RESULT_CONTINUE;
	break;

    case IDA_ROOT_RETURN:
	result = VARIANT_IDA_SOLVER_RESULT_ROOTSFOUND;
	break;

    case IDA_TSTOP_RETURN:
	result = VARIANT_IDA_SOLVER_RESULT_STOPTIMEREACHED;
	break;

    default:
	/* If an exception was recorded, propagate it.  This accounts for
	 * almost all failures except for repeated recoverable failures in the
	 * residue function.  */
	ret = Field (vdata, RECORD_IDA_SESSION_EXN_TEMP);
	if (Is_block (ret)) {
	    Store_field (vdata, RECORD_IDA_SESSION_EXN_TEMP, Val_none);
	    /* FIXME: In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	CHECK_FLAG ("IDASolve", flag);
    }

    /* Hmm...should this go in the production code or not?  */
    if (Is_block (Field (vdata, RECORD_IDA_SESSION_EXN_TEMP)))
	abort ();

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}


CAMLprim value IDATYPE(normal)(value vdata, value nextt, value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 0));
}

CAMLprim value IDATYPE(one_step)(value vdata, value nextt, value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 1));
}

CAMLprim void IDATYPE(get_dky)(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVECTORIZE_VAL(vy);

    int flag = IDAGetDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
			 Int_val(vk), y_nv);
    CHECK_FLAG("IDAGetDky", flag);
    RELINQUISH_NVECTORIZEDVAL(y_nv);

    CAMLreturn0;
}

CAMLprim void IDATYPE(get_err_weights)(value vida_mem, value verrws)
{
    CAMLparam2(vida_mem, verrws);

    N_Vector errws_nv = NVECTORIZE_VAL(verrws);

    int flag = IDAGetErrWeights(IDA_MEM_FROM_ML(vida_mem), errws_nv);
    RELINQUISH_NVECTORIZEDVAL(errws_nv);
    CHECK_FLAG("IDAGetErrWeights", flag);

    CAMLreturn0;
}

CAMLprim void IDATYPE(get_est_local_errors)(value vida_mem, value vele)
{
    CAMLparam2(vida_mem, vele);

    N_Vector ele_nv = NVECTORIZE_VAL(vele);

    int flag = IDAGetEstLocalErrors(IDA_MEM_FROM_ML(vida_mem), ele_nv);
    RELINQUISH_NVECTORIZEDVAL(ele_nv);
    CHECK_FLAG("IDAGetEstLocalErrors", flag);

    CAMLreturn0;
}


CAMLprim void IDATYPE(set_id) (value vida_mem, value vid)
{
    CAMLparam2(vida_mem, vid);
    N_Vector id;

    id = NVECTORIZE_VAL (vid);
    int flag = IDASetId (IDA_MEM_FROM_ML(vida_mem), id);
    RELINQUISH_NVECTORIZEDVAL (id);
    CHECK_FLAG("IDASetId", flag);

    CAMLreturn0;
}


static void calc_ic (void *ida_mem, value session, int icopt, realtype tout1)
{
    CAMLparam1 (session);
    CAMLlocal1 (exn);
    int flag;

    flag = IDACalcIC (ida_mem, icopt, tout1);

    if (flag < 0) {
	/* If an exception is saved in the session, grab it and re-raise. */
	exn = Field (session, RECORD_IDA_SESSION_EXN_TEMP);
	if (Is_block (exn)) {
	    Store_field (session, RECORD_IDA_SESSION_EXN_TEMP, Val_none);
	    /* FIXME: In bytecode, caml_raise() duplicates some parts of
	     * the stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (exn, 0));
	}
	/* Otherwise, raise a generic exception like Ida.ResFuncFailure */
	CHECK_FLAG ("IDACalcIC", flag);
    }
    CAMLreturn0;
}

CAMLprim void IDATYPE(calc_ic_y)(value vida_mem, value tout1)
{
    CAMLparam2 (vida_mem, tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    calc_ic (ida_mem, vida_mem, IDA_Y_INIT, Double_val (tout1));

    CAMLreturn0;
}

CAMLprim void IDATYPE(calc_ic_ya_ydp)(value vida_mem, value vid, value tout1)
{
    CAMLparam3 (vida_mem, vid, tout1);
    int flag;
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    N_Vector id = NVECTORIZE_VAL (vid);
    flag = IDASetId (ida_mem, id);
    RELINQUISH_NVECTORIZEDVAL (id);
    CHECK_FLAG ("IDASetId", flag);

    calc_ic (ida_mem, vida_mem, IDA_YA_YDP_INIT, Double_val (tout1));
    CAMLreturn0;
}

CAMLprim void IDATYPE(set_constraints) (value vida_mem, value vconstraints)
{
    CAMLparam2(vida_mem, vconstraints);
    int flag;

    N_Vector constraints = NVECTORIZE_VAL (vconstraints);
    flag = IDASetConstraints (IDA_MEM_FROM_ML (vida_mem), constraints);
    RELINQUISH_NVECTORIZEDVAL (constraints);
    CHECK_FLAG ("IDASetConstraints", flag);

    CAMLreturn0;
}
