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

#include <errno.h>
#include <string.h>

/* Sundials IDA interface functions that do not involve NVectors. */

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/fail.h>
#include <caml/alloc.h>
#include <caml/callback.h>
#include <caml/bigarray.h>
#include <caml/unixsupport.h>

#ifdef SUNDIALSML_WITHSENS
/* IDAS (with sensitivity) */

#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <idas/idas_band.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_impl.h>
#include <sundials/sundials_config.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <idas/idas_lapack.h>
#endif

#else  /* IDA (without sensitivity) */

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <ida/ida_band.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_sptfqmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_impl.h>
#include <sundials/sundials_config.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <ida/ida_lapack.h>
#endif

#endif

#include <stdio.h>
#include "spils_ml.h"
#include "ida_ml.h"

#define MAX_ERRMSG_LEN 256

#include "ida_ml.h"
#include "nvector_ml.h"
#include "dls_ml.h"


enum callback_index {
    IX_call_resfn = 0,
    IX_call_errh,
    IX_call_errw,
    IX_call_jacfn,
    IX_call_bandjacfn,
    IX_call_precsetupfn,
    IX_call_precsolvefn,
    IX_call_jactimesfn,
    IX_call_linit,
    IX_call_lsetup,
    IX_call_lsolve,
    NUM_CALLBACKS
};

static value callbacks[NUM_CALLBACKS];

CAMLprim value c_ida_init_module (value cbs, value exns)
{
    CAMLparam2 (cbs, exns);
    REGISTER_EXNS (IDA, exns);
    REGISTER_CALLBACKS (cbs);
    CAMLreturn (Val_unit);
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

    caml_callback2_exn (CAML_FN(call_errh), *backref, a);

    CAMLreturn0;
}

CAMLprim value c_ida_set_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), errh,
				  IDA_BACKREF_FROM_ML(vdata));
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_clear_err_handler_fn(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDASetErrHandlerFn(IDA_MEM_FROM_ML(vdata), NULL, NULL);
    CHECK_FLAG("IDASetErrHandlerFn", flag);

    CAMLreturn (Val_unit);
}
#include <stdio.h>
static int resfn (realtype t, N_Vector y, N_Vector yp,
		  N_Vector resval, void *user_data)
{
    CAMLparam0 ();
    CAMLlocalN (args, 5);
    int r;
    value *backref = user_data;

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK (y);
    args[3] = NVEC_BACKLINK (yp);
    args[4] = NVEC_BACKLINK (resval);

    r = Int_val (caml_callbackN (CAML_FN(call_resfn),
				 sizeof (args) / sizeof (*args),
				 args));

    CAMLreturnT (int, r);
}

static value make_jac_arg(realtype t, realtype coef, N_Vector y, N_Vector yp,
			  N_Vector res, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_IDA_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_COEF, caml_copy_double(coef));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_YP, NVEC_BACKLINK(yp));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_RES, NVEC_BACKLINK(res));
    Store_field(r, RECORD_IDA_JACOBIAN_ARG_JAC_TMP, tmp);

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
static int jacfn (long int neq, realtype t, realtype coef,
		  N_Vector y, N_Vector yp, N_Vector res,
		  DlsMat jac, void *user_data,
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 3);
    int r;
    value *backref = user_data;

    args[0] = *backref;
    args[1] = make_jac_arg (t, coef, y, yp, res,
			    make_triple_tmp (tmp1, tmp2, tmp3));
    args[2] = c_dls_wrap(jac, 0); // TODO: cache for efficiency!

    r = Int_val (caml_callbackN (CAML_FN(call_jacfn),
				 sizeof (args) / sizeof (*args),
				 args));
    c_dls_relinquish(Field(args[2], 1)); // TODO: cache for efficiency!

    CAMLreturnT (int, r);
}

static int bandjacfn (long int neq, long int mupper, long int mlower,
		      realtype t, realtype coef, N_Vector y, N_Vector yp,
		      N_Vector res, DlsMat jac, void *user_data,
		      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0 ();
    CAMLlocalN (args, 4);
    int r;
    value *backref = user_data;

    args[0] = *backref;
    args[1] = caml_alloc_tuple(RECORD_IDA_BANDRANGE_SIZE);
    Store_field(args[1], RECORD_IDA_BANDRANGE_MUPPER, Val_long(mupper));
    Store_field(args[1], RECORD_IDA_BANDRANGE_MLOWER, Val_long(mlower));
    args[2] = make_jac_arg (t, coef, y, yp, res,
			    make_triple_tmp (tmp1, tmp2, tmp3));
    args[3] = c_dls_wrap(jac, 0); // TODO: cache for efficiency!

    r = Int_val (caml_callbackN (CAML_FN(call_bandjacfn),
				 sizeof (args) / sizeof (*args),
				 args));
    c_dls_relinquish(Field(args[3], 1)); // TODO: cache for efficiency!

    CAMLreturnT (int, r);
}

static int check_exception(value session, value r)
{
    CAMLparam2(session, r);
    CAMLlocal1(exn);

    if (!Is_exception_result(r)) return 0;

    r = Extract_exception(r);

    if (Field(r, 0) == SUNDIALS_EXN (RecoverableFailure))
	CAMLreturnT (int, 1);

    /* Unrecoverable error.  Save the exception and return -1.  */
    exn = caml_alloc_small (1,0);
    Field (exn, 0) = r;
    Store_field (session, RECORD_IDA_SESSION_EXN_TEMP, exn);
    CAMLreturnT (int, -1);
}

static int rootsfn (realtype t, N_Vector y, N_Vector yp,
		    realtype *gout, void *user_data)
{
    CAMLparam0 ();
    CAMLlocal2 (session, r);
    CAMLlocalN (args, 4);
    value *backref = (value *)user_data;
    intnat nroots;

    /* The length of gout is only available at the nroots field of the session
     * structure, so a dereference of the backref is unavoidable.  Therefore,
     * we will do all of the setup here and directly call the user-supplied
     * OCaml function without going through an OCaml trampoline.  */
    session = sundials_ml_weak_get (*backref, Val_int (0));
    if (!Is_block (session))
	caml_failwith ("Internal error: weak reference is dead");
    session = Field (session, 0);

    nroots = IDA_NROOTS_FROM_ML (session);

    args[0] = caml_copy_double (t);
    args[1] = NVEC_BACKLINK (y);
    args[2] = NVEC_BACKLINK (yp);
    args[3] = caml_ba_alloc (BIGARRAY_FLOAT, 1, gout, &nroots);

    r = caml_callbackN_exn (IDA_ROOTSFN_FROM_ML (session),
			    sizeof (args) / sizeof (*args),
			    args);

    CAMLreturnT (int, check_exception (session, r));
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    int r;
    value *backref = user_data;

    args[0] = *backref;
    args[1] = NVEC_BACKLINK (y);
    args[2] = NVEC_BACKLINK (ewt);

    r = Int_val (caml_callbackN (CAML_FN(call_errw),
				 sizeof (args) / sizeof (*args),
				 args));

    CAMLreturnT(int, r);
}

static int precsetupfn(
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

    args[0] = *backref;
    args[1] = make_jac_arg(t, cj, y, yp, res,
			   make_triple_tmp(tmp1, tmp2, tmp3));
    r = Int_val (caml_callbackN (CAML_FN(call_precsetupfn),
				 sizeof (args) / sizeof (*args),
				 args));

    CAMLreturnT(int, r);
}

static int precsolvefn(
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

    args[0] = *backref;
    args[1] = make_jac_arg(t, cj, y, yp, res, NVEC_BACKLINK (tmp));
    args[2] = NVEC_BACKLINK (r);
    args[3] = NVEC_BACKLINK (z);
    args[4] = caml_copy_double (delta);

    rv = Int_val (caml_callbackN (CAML_FN(call_precsolvefn),
				  sizeof (args) / sizeof (*args),
				  args));

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

    args[0] = *backref;
    args[1] = make_jac_arg (t, cj, y, yp, res, make_double_tmp (tmp1, tmp2));
    args[2] = NVEC_BACKLINK (v);
    args[3] = NVEC_BACKLINK (Jv);

    r = Int_val (caml_callbackN (CAML_FN(call_jactimesfn),
				 sizeof (args) / sizeof (*args),
				 args));

    CAMLreturnT (int, r);
}

static int linit(IDAMem ida_mem)
{
    CAMLparam0();
    int r;
    value *backref = ida_mem->ida_user_data;

    r = Int_val (caml_callback(CAML_FN(call_linit), *backref));

    CAMLreturnT(int, r);
}

static int lsetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp, N_Vector resp,
		  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    int r;
    value *backref = ida_mem->ida_user_data;

    args[0] = *backref;
    args[1] = NVEC_BACKLINK(yyp);
    args[2] = NVEC_BACKLINK(ypp);
    args[3] = NVEC_BACKLINK(resp);
    args[4] = make_triple_tmp(tmp1, tmp2, tmp3);

    r = Int_val (caml_callbackN(CAML_FN(call_lsetup),
				sizeof (args) / sizeof (*args),
				args));

    CAMLreturnT(int, r);
}

static int lsolve(IDAMem ida_mem, N_Vector b, N_Vector weight, N_Vector ycur,
		  N_Vector ypcur, N_Vector rescur)
{
    CAMLparam0();
    CAMLlocalN(args, 6);
    int r;
    value *backref = ida_mem->ida_user_data;

    args[0] = *backref;
    args[1] = NVEC_BACKLINK(b);
    args[2] = NVEC_BACKLINK(weight);
    args[3] = NVEC_BACKLINK(ycur);
    args[4] = NVEC_BACKLINK(ypcur);
    args[5] = NVEC_BACKLINK(rescur);

    r = Int_val (caml_callbackN(CAML_FN(call_lsolve),
				sizeof (args) / sizeof (*args),
				args));

    CAMLreturnT(int, r);
}

CAMLprim value c_ida_set_alternate (value vida_mem, value vhas_init,
				    value vhas_setup)
{
    CAMLparam3(vida_mem, vhas_init, vhas_setup);
    IDAMem ida_mem = IDA_MEM_FROM_ML (vida_mem);

    ida_mem->ida_linit   = Bool_val(vhas_init)  ? linit : NULL;
    ida_mem->ida_lsetup  = Bool_val(vhas_setup) ? lsetup : NULL;
    ida_mem->ida_setupNonNull = Bool_val(vhas_setup);
    ida_mem->ida_lsolve = lsolve;
    ida_mem->ida_lmem   = NULL;

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_get_cj (value vida_mem)
{
    CAMLparam1 (vida_mem);
    IDAMem ida_mem = IDA_MEM_FROM_ML (vida_mem);
    CAMLreturn (caml_copy_double (ida_mem->ida_cj));
}

CAMLprim value c_ida_get_cjratio (value vida_mem)
{
    CAMLparam1 (vida_mem);
    IDAMem ida_mem = IDA_MEM_FROM_ML (vida_mem);
    CAMLreturn (caml_copy_double (ida_mem->ida_cjratio));
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim value c_ida_dls_dense (value vida_mem, value vneqs, value vset_jac)
{
    CAMLparam3(vida_mem, vneqs, vset_jac);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDADense (ida_mem, neqs);
    CHECK_FLAG ("IDADense", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetDenseJacFn(IDA_MEM_FROM_ML(vida_mem), jacfn);
	CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    }
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_dls_lapack_dense (value vida_mem, value vneqs,
				       value vset_jac)
{
    CAMLparam3 (vida_mem, vneqs, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDALapackDense (ida_mem, neqs);
    CHECK_FLAG ("IDALapackDense", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetDenseJacFn (IDA_MEM_FROM_ML (vida_mem), jacfn);
	CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_dls_set_dense_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetDenseJacFn(IDA_MEM_FROM_ML(vdata), jacfn);
    CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_dls_clear_dense_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetDenseJacFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDADlsSetDenseJacFn", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_dls_band (value vida_mem, value vneqs,
			       value mupper, value mlower, value vset_jac)
{
    CAMLparam5(vida_mem, vneqs, mupper, mlower, vset_jac);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDABand (ida_mem, neqs, Long_val (mupper), Long_val (mlower));
    CHECK_FLAG ("IDABand", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vida_mem), bandjacfn);
	CHECK_FLAG("IDADlsSetBandJacFn", flag);
    }
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_dls_lapack_band (value vida_mem, value vneqs,
				      value mupper, value mlower,
				      value vset_jac)
{
    CAMLparam5(vida_mem, vneqs, mupper, mlower, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    long neqs = Long_val (vneqs);
    int flag;

    flag = IDALapackBand (ida_mem, neqs, Long_val (mupper), Long_val (mlower));
    CHECK_FLAG ("IDALapackBand", flag);
    if (Bool_val (vset_jac)) {
	flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vida_mem), bandjacfn);
	CHECK_FLAG("IDADlsSetBandJacFn", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_dls_set_band_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vdata), bandjacfn);
    CHECK_FLAG("IDADlsSetBandJacFn", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_dls_clear_band_jac_fn(value vdata)
{
    CAMLparam1(vdata);
    int flag = IDADlsSetBandJacFn(IDA_MEM_FROM_ML(vdata), NULL);
    CHECK_FLAG("IDADlsSetBandJacFn", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_set_preconditioner (value vsession,
					       value vset_presetup)
{
    CAMLparam2 (vsession, vset_presetup);
    void *mem = IDA_MEM_FROM_ML (vsession);
    IDASpilsPrecSetupFn setup = Bool_val (vset_presetup) ? precsetupfn : NULL;
    int flag = IDASpilsSetPreconditioner (mem, setup, precsolvefn);
    CHECK_FLAG ("IDASpilsSetPreconditioner", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_set_jac_times_vec_fn(value vdata, value vset_jac)
{
    CAMLparam2(vdata, vset_jac);
    IDASpilsJacTimesVecFn jac = Bool_val (vset_jac) ? jactimesfn : NULL;
    int flag = IDASpilsSetJacTimesVecFn(IDA_MEM_FROM_ML(vdata), jac);
    CHECK_FLAG("IDASpilsSetJacTimesVecFn", flag);
    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_spgmr (value vida_mem, value vmaxl)
{
    CAMLparam2 (vida_mem, vmaxl);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASpgmr (ida_mem, Int_val (vmaxl));
    CHECK_FLAG ("IDASpgmr", flag);

    CAMLreturn (Val_unit);
}

/* For SPGMR only.  */
CAMLprim value c_ida_spils_set_max_restarts (value vida_mem, value vmaxr)
{
    CAMLparam2 (vida_mem, vmaxr);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASpilsSetMaxRestarts (ida_mem, Int_val (vmaxr));
    CHECK_FLAG ("IDASpilsSetMaxRestarts", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_spbcg (value vida_mem, value vmaxl)
{
    CAMLparam2 (vida_mem, vmaxl);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASpbcg (ida_mem, Int_val (vmaxl));
    CHECK_FLAG ("IDASpbcg", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_sptfqmr (value vida_mem, value vmaxl)
{
    CAMLparam2 (vida_mem, vmaxl);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDASptfqmr (ida_mem, Int_val (vmaxl));
    CHECK_FLAG ("IDASptfqmr", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_wf_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = IDAWFtolerances(IDA_MEM_FROM_ML(vdata), errw);
    CHECK_FLAG("IDAWFtolerances", flag);

    CAMLreturn (Val_unit);
}

/* Sets the root function to a generic trampoline and set the number of
 * roots.  */
CAMLprim value c_ida_root_init (value vida_mem, value vnroots)
{
    CAMLparam2 (vida_mem, vnroots);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int nroots = Int_val (vnroots);
    int flag;
    flag = IDARootInit (ida_mem, nroots, rootsfn);
    CHECK_FLAG ("IDARootInit", flag);
    Store_field (vida_mem, RECORD_IDA_SESSION_NROOTS, vnroots);
    CAMLreturn (Val_unit);
}

/* IDACreate + IDAInit + IDASetUser.  The vy and vyp vectors must have the same
 * length (not checked).  Returns an IDAMem (= void*) along with a global root
 * wrapping weakref that is attached to the IDAMem.  Before doing the remaining
 * initialization, these pointers have to be wrapped in a session so that if
 * initialization fails the GC frees IDAMem and the weak global root.  Doing
 * that cleanup in C would be ugly, e.g. we wouldn't be able to reuse
 * c_ida_set_linear_solver() so we have to duplicate it or hack some ad-hoc
 * extensions to that function.  */
CAMLprim value c_ida_init (value weakref, value vt0, value vy, value vyp)
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

    y = NVEC_VAL (vy);
    yp = NVEC_VAL (vyp);
    flag = IDAInit (ida_mem, resfn, Double_val (vt0), y, yp);
    if (flag != IDA_SUCCESS) {
	IDAFree (ida_mem);
	CHECK_FLAG ("IDAInit", flag);
    }

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
    Store_field(r, 2, Val_long (0)); // no err_file = NULL

    CAMLreturn (r);
}

CAMLprim value c_ida_sv_tolerances (value ida_mem, value vrtol, value vavtol)
{
    CAMLparam3 (ida_mem, vrtol, vavtol);
    N_Vector avtol;
    int flag;

    avtol = NVEC_VAL (vavtol);
    flag = IDASVtolerances (IDA_MEM_FROM_ML (ida_mem),
			    Double_val (vrtol), avtol);
    CHECK_FLAG ("IDASVtolerances", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_reinit(value vdata, value t0, value y0, value yp0)
{
    CAMLparam4(vdata, t0, y0, yp0);

    N_Vector y0_nv = NVEC_VAL(y0);
    N_Vector yp0_nv = NVEC_VAL(yp0);
    int flag = IDAReInit(IDA_MEM_FROM_ML(vdata), Double_val(t0),
			 y0_nv, yp0_nv);
    CHECK_FLAG("IDAReInit", flag);

#if SAFETY_CHECKS
    IDA_MASK_SAFETY_FLAGS(vdata, IDA_SAFETY_FLAG_REINIT_KEEPS);
#endif

    CAMLreturn (Val_unit);
}

static value solve (value vdata, value nextt, value vy, value vyp, int onestep)
{
    CAMLparam4 (vdata, nextt, vy, vyp);
    CAMLlocal1 (ret);
    void *ida_mem = IDA_MEM_FROM_ML (vdata);
    realtype tret;
    int flag;
    N_Vector y, yp;
    enum ida_solver_result_tag result = -1;

#if SAFETY_CHECKS && IDA_ML_BIGARRAYS
    /* This can't be checked for generic nvectors.  */
    if (IDA_NEQS_FROM_ML (vdata) != Caml_ba_array_val(vy)->dim[0])
	caml_invalid_argument ("Ida.solve: y vector has incorrect length");
    if (IDA_NEQS_FROM_ML (vdata) != Caml_ba_array_val(vyp)->dim[0])
	caml_invalid_argument ("Ida.solve: y' vector has incorrect length");
#endif

    y = NVEC_VAL (vy);
    yp = NVEC_VAL (vyp);
    flag = IDASolve (ida_mem, Double_val (nextt), &tret, y, yp,
	             onestep ? IDA_ONE_STEP : IDA_NORMAL);

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
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (ret, 0));
	}
	CHECK_FLAG ("IDASolve", flag);
    }

#if SAFETY_CHECKS
    IDA_SET_SAFETY_FLAG (vdata, IDA_SAFETY_FLAG_SOLVING);
#endif

    /* Hmm...should this go in the production code or not?  */
    if (Is_block (Field (vdata, RECORD_IDA_SESSION_EXN_TEMP)))
	abort ();

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, caml_copy_double (tret));
    Store_field (ret, 1, Val_int (result));

    CAMLreturn (ret);
}


CAMLprim value c_ida_solve_normal (value vdata, value nextt,
				   value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 0));
}

CAMLprim value c_ida_solve_one_step (value vdata, value nextt,
				     value y, value yp)
{
    CAMLparam4(vdata, nextt, y, yp);
    CAMLreturn(solve(vdata, nextt, y, yp, 1));
}

CAMLprim value c_ida_get_dky(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    N_Vector y_nv = NVEC_VAL(vy);

    int flag = IDAGetDky(IDA_MEM_FROM_ML(vdata), Double_val(vt),
			 Int_val(vk), y_nv);
    CHECK_FLAG("IDAGetDky", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_get_err_weights(value vida_mem, value verrws)
{
    CAMLparam2(vida_mem, verrws);

    N_Vector errws_nv = NVEC_VAL(verrws);

    int flag = IDAGetErrWeights(IDA_MEM_FROM_ML(vida_mem), errws_nv);
    CHECK_FLAG("IDAGetErrWeights", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_get_est_local_errors(value vida_mem, value vele)
{
    CAMLparam2(vida_mem, vele);

    N_Vector ele_nv = NVEC_VAL(vele);

    int flag = IDAGetEstLocalErrors(IDA_MEM_FROM_ML(vida_mem), ele_nv);
    CHECK_FLAG("IDAGetEstLocalErrors", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value c_ida_set_id (value vida_mem, value vid)
{
    CAMLparam2(vida_mem, vid);
    N_Vector id;

    id = NVEC_VAL (vid);
    int flag = IDASetId (IDA_MEM_FROM_ML(vida_mem), id);
    CHECK_FLAG("IDASetId", flag);

#if SAFETY_CHECKS
    IDA_SET_SAFETY_FLAG (vida_mem, IDA_SAFETY_FLAG_ID_SET);
#endif

    CAMLreturn (Val_unit);
}


static void calc_ic (void *ida_mem, value session, int icopt, realtype tout1,
		     value vy, value vyp)
{
    CAMLparam3 (session, vy, vyp);
    CAMLlocal1 (exn);
    int flag;
    N_Vector y, yp;

#if SAFETY_CHECKS
    if (IDA_SAFETY_FLAGS (session) & IDA_SAFETY_FLAG_SOLVING)
	caml_invalid_argument ("Ida.calc_ic: called after Ida.solve_*");
#endif

    flag = IDACalcIC (ida_mem, icopt, tout1);

    if (flag < 0) {
	/* If an exception is saved in the session, grab it and re-raise. */
	exn = Field (session, RECORD_IDA_SESSION_EXN_TEMP);
	if (Is_block (exn)) {
	    Store_field (session, RECORD_IDA_SESSION_EXN_TEMP, Val_none);
	    /* In bytecode, caml_raise() duplicates some parts of the
	     * stacktrace.  This does not seem to happen in native code
	     * execution.  */
	    caml_raise (Field (exn, 0));
	}
	/* Otherwise, raise a generic exception like Ida.ResFuncFailure */
	CHECK_FLAG ("IDACalcIC", flag);
    }

    /* Retrieve the calculated initial conditions if y,yp are given.  */
    y  = (Is_block (vy))  ? NVEC_VAL (Field (vy, 0))  : NULL;
    yp = (Is_block (vyp)) ? NVEC_VAL (Field (vyp, 0)) : NULL;
    if (y != NULL || yp != NULL) {
	flag = IDAGetConsistentIC (ida_mem, y, yp);
	CHECK_FLAG ("IDAGetConsistentIC", flag);
    }
    CAMLreturn0;
}

CAMLprim value c_ida_calc_ic_y(value vida_mem, value vy, value tout1)
{
    CAMLparam3 (vida_mem, vy, tout1);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    calc_ic (ida_mem, vida_mem, IDA_Y_INIT, Double_val (tout1), vy, Val_none);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_calc_ic_ya_ydp(value vida_mem, value y, value yp,
				    value vid, value tout1)
{
    CAMLparam5 (vida_mem, y, yp, vid, tout1);
    int flag;
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);

    N_Vector id = NVEC_VAL (vid);
    flag = IDASetId (ida_mem, id);
    CHECK_FLAG ("IDASetId", flag);

#if SAFETY_CHECKS
    IDA_SET_SAFETY_FLAG (vida_mem, IDA_SAFETY_FLAG_ID_SET);
#endif

    calc_ic (ida_mem, vida_mem, IDA_YA_YDP_INIT, Double_val (tout1), y, yp);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_constraints (value vida_mem, value vconstraints)
{
    CAMLparam2(vida_mem, vconstraints);
    int flag;

    N_Vector constraints = NVEC_VAL (vconstraints);
    flag = IDASetConstraints (IDA_MEM_FROM_ML (vida_mem), constraints);
    CHECK_FLAG ("IDASetConstraints", flag);

    CAMLreturn (Val_unit);
}

void ida_ml_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == IDA_SUCCESS
	|| flag == IDA_ROOT_RETURN
	|| flag == IDA_TSTOP_RETURN) return;

    switch (flag) {
    case IDA_ILL_INPUT:
	caml_raise_constant(IDA_EXN(IllInput));

    case IDA_CONV_FAIL:
	caml_raise_constant(IDA_EXN(ConvergenceFailure));

    case IDA_TOO_MUCH_WORK:
	caml_raise_constant(IDA_EXN(TooMuchWork));

    case IDA_TOO_MUCH_ACC:
	caml_raise_constant(IDA_EXN(TooMuchAccuracy));

    case IDA_LINIT_FAIL:
	caml_raise_constant(IDA_EXN(LinearInitFailure));

    case IDA_LSETUP_FAIL:
	caml_raise_constant(IDA_EXN(LinearSetupFailure));

    case IDA_LSOLVE_FAIL:
	caml_raise_constant(IDA_EXN(LinearSolveFailure));

    case IDA_BAD_EWT:
	caml_raise_constant(IDA_EXN(BadEwt));

    case IDA_NO_RECOVERY:
	caml_raise_constant(IDA_EXN(NoRecovery));

    case IDA_RES_FAIL:
	caml_raise_constant(IDA_EXN(ResFuncFailure));

    case IDA_FIRST_RES_FAIL:
	caml_raise_constant(IDA_EXN(FirstResFuncFailure));

    case IDA_REP_RES_ERR:
	caml_raise_constant(IDA_EXN(RepeatedResFuncFailure));

    case IDA_RTFUNC_FAIL:
	caml_raise_constant(IDA_EXN(RootFuncFailure));

    case IDA_CONSTR_FAIL:
	caml_raise_constant(IDA_EXN(ConstraintFailure));

    case IDA_BAD_K:
	caml_raise_constant(IDA_EXN(BadK));

    case IDA_BAD_T:
	caml_raise_constant(IDA_EXN(BadT));

    case IDA_BAD_DKY:
	caml_raise_constant(IDA_EXN(BadDky));

    default:
	/* e.g. IDA_MEM_NULL, IDA_MEM_FAIL */
	snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		 IDAGetReturnFlagName(flag));
	caml_failwith(exmsg);
    }
}

CAMLprim value c_ida_session_finalize(value vdata)
{
    if (IDA_MEM_FROM_ML(vdata) != NULL) {
	void *ida_mem = IDA_MEM_FROM_ML(vdata);
	value *backref = IDA_BACKREF_FROM_ML(vdata);
	IDAFree(&ida_mem);
	caml_remove_generational_global_root (backref);
	free (backref);
    }

    FILE* err_file = (FILE *)Long_val(Field(vdata, RECORD_IDA_SESSION_ERRFILE));
    if (err_file != NULL) {
	fclose(err_file);
    }
    return Val_unit;
}

 
CAMLprim value c_ida_ss_tolerances(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = IDASStolerances(IDA_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    CHECK_FLAG("IDASStolerances", flag);

    CAMLreturn (Val_unit);
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

    CAMLreturn (Val_unit);
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

    r = caml_alloc_tuple(RECORD_IDA_INTEGRATOR_STATS_SIZE);
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

    FILE* err_file = (FILE*)Long_val(Field(vdata, RECORD_IDA_SESSION_ERRFILE));

    if (err_file != NULL) {
	fclose(err_file);
	Store_field(vdata, RECORD_IDA_SESSION_ERRFILE, 0);
    }
    char *mode = Bool_val(vtrunc) ? "w" : "a";
    err_file = fopen(String_val(vpath), mode);
    if (err_file == NULL) {
	// uerror("fopen", vpath); /* depends on unix.cma */
	caml_failwith(strerror(errno));
    }

    int flag = IDASetErrFile(IDA_MEM_FROM_ML(vdata), err_file);
    CHECK_FLAG("IDASetErrFile", flag);

    Store_field(vdata, RECORD_IDA_SESSION_ERRFILE, Val_long(err_file));

    CAMLreturn (Val_unit);
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

    CAMLreturn (Val_unit);
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

CAMLprim value c_ida_get_num_backtrack_ops (value vida_mem)
{
    CAMLparam1 (vida_mem);
    int flag;
    long nbo;
    flag = IDAGetNumBacktrackOps (IDA_MEM_FROM_ML (vida_mem), &nbo);
    CHECK_FLAG ("IDAGetNumBcktrackOps", flag);
    CAMLreturn (Val_int (nbo));
}

CAMLprim value c_ida_set_max_ord(value vida_mem, value maxord)
{
    CAMLparam2(vida_mem, maxord);


    int flag = IDASetMaxOrd(IDA_MEM_FROM_ML(vida_mem), Int_val(maxord));
    CHECK_FLAG("IDASetMaxOrd", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_max_num_steps(value vida_mem, value mxsteps)
{
    CAMLparam2(vida_mem, mxsteps);


    int flag = IDASetMaxNumSteps(IDA_MEM_FROM_ML(vida_mem), Long_val(mxsteps));
    CHECK_FLAG("IDASetMaxNumSteps", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_init_step(value vida_mem, value hin)
{
    CAMLparam2(vida_mem, hin);


    int flag = IDASetInitStep(IDA_MEM_FROM_ML(vida_mem), Double_val(hin));
    CHECK_FLAG("IDASetInitStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_max_step(value vida_mem, value hmax)
{
    CAMLparam2(vida_mem, hmax);


    int flag = IDASetMaxStep(IDA_MEM_FROM_ML(vida_mem), Double_val(hmax));
    CHECK_FLAG("IDASetMaxStep", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_stop_time(value vida_mem, value tstop)
{
    CAMLparam2(vida_mem, tstop);


    int flag = IDASetStopTime(IDA_MEM_FROM_ML(vida_mem), Double_val(tstop));
    CHECK_FLAG("IDASetStopTime", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_max_err_test_fails(value vida_mem, value maxnef)
{
    CAMLparam2(vida_mem, maxnef);


    int flag = IDASetMaxErrTestFails(IDA_MEM_FROM_ML(vida_mem), Int_val(maxnef));
    CHECK_FLAG("IDASetMaxErrTestFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_max_nonlin_iters(value vida_mem, value maxcor)
{
    CAMLparam2(vida_mem, maxcor);


    int flag = IDASetMaxNonlinIters(IDA_MEM_FROM_ML(vida_mem), Int_val(maxcor));
    CHECK_FLAG("IDASetMaxNonlinIters", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_max_conv_fails(value vida_mem, value maxncf)
{
    CAMLparam2(vida_mem, maxncf);


    int flag = IDASetMaxConvFails(IDA_MEM_FROM_ML(vida_mem), Int_val(maxncf));
    CHECK_FLAG("IDASetMaxConvFails", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_nonlin_conv_coef(value vida_mem, value nlscoef)
{
    CAMLparam2(vida_mem, nlscoef);


    int flag = IDASetNonlinConvCoef(IDA_MEM_FROM_ML(vida_mem), Double_val(nlscoef));
    CHECK_FLAG("IDASetNonlinConvCoef", flag);

    CAMLreturn (Val_unit);
}


CAMLprim value c_ida_set_no_inactive_root_warn(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag = IDASetNoInactiveRootWarn(IDA_MEM_FROM_ML(vida_mem));
    CHECK_FLAG("IDASetNoInactiveRootWarn", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_set_suppress_alg (value vida_mem, value vb)
{
    CAMLparam2(vida_mem, vb);

#if SAFETY_CHECKS
    if (! IDA_TEST_SAFETY_FLAG (vida_mem, IDA_SAFETY_FLAG_ID_SET))
	caml_invalid_argument ("Ida.set_suppress_alg: var types not set");
#endif

    int flag = IDASetSuppressAlg(IDA_MEM_FROM_ML(vida_mem),
				 Bool_val (vb));
    CHECK_FLAG("IDASetSuppressAlg", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_set_gs_type(value vida_mem, value vgstype)
{
    CAMLparam2(vida_mem, vgstype);

    int flag = IDASpilsSetGSType(IDA_MEM_FROM_ML(vida_mem),
				 spils_gs_type(vgstype));
    CHECK_FLAG("IDASpilsSetGSType", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_set_eps_lin(value vida_mem, value eplifac)
{
    CAMLparam2(vida_mem, eplifac);

    int flag = IDASpilsSetEpsLin(IDA_MEM_FROM_ML(vida_mem),
				 Double_val(eplifac));
    CHECK_FLAG("IDASpilsSetEpsLin", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_spils_set_maxl(value vida_mem, value maxl)
{
    CAMLparam2(vida_mem, maxl);

    int flag = IDASpilsSetMaxl(IDA_MEM_FROM_ML(vida_mem), Int_val(maxl));
    CHECK_FLAG("IDASpilsSetMaxl", flag);

    CAMLreturn (Val_unit);
}


/* statistic accessor functions */

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

CAMLprim value c_ida_get_nonlin_solv_stats(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(ret);

    long int nniters, nncfails;
    int flag = IDAGetNonlinSolvStats(IDA_MEM_FROM_ML(vida_mem),
				     &nniters, &nncfails);
    CHECK_FLAG("IDAGetNonlinSolvStats", flag);

    ret = caml_alloc_tuple (2);
    Store_field (ret, 0, Val_long (nniters));
    Store_field (ret, 1, Val_long (nncfails));

    CAMLreturn(ret);
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

/* spils functions */

CAMLprim value c_ida_spils_get_num_lin_iters(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDASpilsGetNumLinIters(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASpilsGetNumLinIters", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_spils_get_num_conv_fails(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDASpilsGetNumConvFails(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASpilsGetNumConvFails", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_spils_get_work_space(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = IDASpilsGetWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrw, &leniw);
    CHECK_FLAG("IDASpilsGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_int(lenrw));
    Store_field(r, 1, Val_int(leniw));

    CAMLreturn(r);
}

CAMLprim value c_ida_spils_get_num_prec_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDASpilsGetNumPrecEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASpilsGetNumPrecEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_spils_get_num_prec_solves(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDASpilsGetNumPrecSolves(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASpilsGetNumPrecSolves", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_spils_get_num_jtimes_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDASpilsGetNumJtimesEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASpilsGetNumJtimesEvals", flag);

    CAMLreturn(Val_long(r));
}

CAMLprim value c_ida_spils_get_num_res_evals (value vida_mem)
{
    CAMLparam1(vida_mem);

    long int r;
    int flag = IDASpilsGetNumResEvals(IDA_MEM_FROM_ML(vida_mem), &r);
    CHECK_FLAG("IDASpilsGetNumResEvals", flag);

    CAMLreturn(Val_long(r));
}

