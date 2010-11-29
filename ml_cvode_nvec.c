
/* Nov 2010, Timothy Bourke (INRIA)
 *
 * The parts of the Sundials interface that distinguish between Serial
 * NVectors (handled by Bigarrays) and generic NVectors (handled by a
 * wrapper type).
 *
 */

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

#include "ml_cvode.h"
#include "ml_nvector.h"

#ifdef RESTRICT_INTERNAL_PRECISION
#ifdef __GNUC__
#include <fpu_control.h>
#endif
#endif

// Call with _ML_CVODE_BIGARRAYS to compile for the Serial NVector to
// Bigarray interface code.

#ifdef ML_CVODE_BIGARRAYS

#define CVTYPE(fname) c_ba_ ## fname
#include <nvector/nvector_serial.h>

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0

#define NVECTORIZE_VAL(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDVAL(nv) N_VDestroy(nv)

#else

#define CVTYPE(fname) c_nvec_ ## fname
#include <sundials/sundials_nvector.h>

#define WRAP_NVECTOR(v) NVEC_DATA(v)
#define RELINQUISH_WRAPPEDNV(v) {}

#define NVECTORIZE_VAL(v) NVEC_VAL(v)
#define RELINQUISH_NVECTORIZEDVAL(nv) {}

#endif

/* callbacks */

static int check_exception(value r)
{
    CAMLparam0();
    CAMLlocal1(exn);

    static value *recoverable_failure = NULL;

    if (!Is_exception_result(r)) return 0;

    if (recoverable_failure == NULL) {
	recoverable_failure =
	    caml_named_value("cvode_RecoverableFailure");
    }

    exn = Extract_exception(r);
    CAMLreturn((Field(exn, 0) == *recoverable_failure) ? 1 : -1);
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(y_d, ydot_d, r);

    value *closure_rhsfn = ((ml_cvode_data_p)user_data)->closure_rhsfn;

    y_d = WRAP_NVECTOR(y);
    ydot_d = WRAP_NVECTOR(ydot);

    // TODO: the data payloads inside y_d and ydot_d are only valid
    //	     during this call, afterward that memory goes back to cvode.
    //	     These bigarrays must not be retained by closure_rhsfn! If
    //	     it wants a permanent copy, then it has to make it manually.
    //
    //       Eventually y_d and ydot_d will be reclaimed by the ocaml gc,
    //       which should not, however, free the attached payload.
    r = caml_callback3_exn(*closure_rhsfn, caml_copy_double(t), y_d, ydot_d);

    RELINQUISH_WRAPPEDNV(y_d);
    RELINQUISH_WRAPPEDNV(ydot_d);

    CAMLreturn(check_exception(r));
}

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(y_d, gout_d, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    y_d = WRAP_NVECTOR(y);

    gout_d = caml_ba_alloc(BIGARRAY_FLOAT, 1, gout, &(data->num_roots));

    // TODO: see notes for f()
    r = caml_callback3_exn(*(data->closure_rootsfn), caml_copy_double(t),
				 y_d, gout_d);

    RELINQUISH_WRAPPEDNV(y_d);
    Caml_ba_array_val(gout_d)->dim[0] = 0;

    CAMLreturn(check_exception(r));
}

static int errw(N_Vector y, N_Vector ewt, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(y_d, ewt_d, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    y_d = WRAP_NVECTOR(y);
    ewt_d = WRAP_NVECTOR(ewt);

    // TODO: see notes for f()
    r = caml_callback2_exn(*(data->closure_errw), y_d, ewt_d);

    RELINQUISH_WRAPPEDNV(y_d);
    RELINQUISH_WRAPPEDNV(ewt_d);

    CAMLreturn(check_exception(r));
}

static value make_jac_arg(realtype t, N_Vector y, N_Vector fy, value tmp)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_tuple(4);
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_Y, WRAP_NVECTOR(y));
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_FY, WRAP_NVECTOR(fy));
    Store_field(r, RECORD_JACOBIAN_ARG_JAC_TMP, tmp);

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

static void relinquish_jac_arg(value arg, int triple)
{
    CAMLparam0();
    CAMLlocal1(tmp);

    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_JACOBIAN_ARG_JAC_Y));
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_JACOBIAN_ARG_JAC_FY));

    tmp = Field(arg, RECORD_JACOBIAN_ARG_JAC_TMP);

    if (triple) {
	RELINQUISH_WRAPPEDNV(Field(tmp, 0));
	RELINQUISH_WRAPPEDNV(Field(tmp, 1));
	RELINQUISH_WRAPPEDNV(Field(tmp, 2));
    } else {
	RELINQUISH_WRAPPEDNV(tmp);
    }

    CAMLreturn0;
}

static int jacfn(
	int n,
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
    CAMLlocal3(arg, r, matrix);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));

    matrix = caml_alloc(1, Abstract_tag);
    Store_field(matrix, 0, (value)Jac);

    r = caml_callback2_exn(*(data->closure_jacfn), arg, matrix);
    relinquish_jac_arg(arg, 1);
    // note: matrix is also invalid after the callback

    CAMLreturn(check_exception(r));
}

static int bandjacfn(
	int N,
	int mupper,
	int mlower, 	 
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
    CAMLlocal1(r);
    CAMLlocalN(args, 4);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    args[0] = Val_int(mupper);
    args[1] = Val_int(mlower);
    args[2] = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));
    args[3] = caml_alloc(1, Abstract_tag);
    Store_field(args[3], 0, (value)Jac);

    r = caml_callbackN_exn(*(data->closure_bandjacfn), 4, args);

    relinquish_jac_arg(args[2], 1);
    // note: matrix is also invalid after the callback

    CAMLreturn(check_exception(r));
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
    CAMLlocal2(arg, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, make_triple_tmp(tmp1, tmp2, tmp3));

    r = caml_callback3_exn(*(data->closure_presetupfn),
	    arg, Val_bool(jok), caml_copy_double(gamma));

    relinquish_jac_arg(arg, 1);

    if (!Is_exception_result(r)) {
	*jcurPtr = Bool_val(r);
    }

    CAMLreturn(check_exception(r));
}

static value make_spils_solve_arg(
	N_Vector r,
	realtype gamma,
	realtype delta,
	int lr)

{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(4);
    Store_field(v, RECORD_SPILS_SOLVE_ARG_RHS, WRAP_NVECTOR(r));
    Store_field(v, RECORD_SPILS_SOLVE_ARG_GAMMA, caml_copy_double(gamma));
    Store_field(v, RECORD_SPILS_SOLVE_ARG_DELTA, caml_copy_double(delta));
    Store_field(v, RECORD_SPILS_SOLVE_ARG_LEFT, lr == 1 ? Val_true : Val_false);

    CAMLreturn(v);
}

static relinquish_spils_solve_arg(value arg)
{
    CAMLparam0();
    RELINQUISH_WRAPPEDNV(Field(arg, RECORD_SPILS_SOLVE_ARG_RHS));
    CAMLreturn0;
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
    CAMLlocal4(arg, solvearg, zv, rv);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, WRAP_NVECTOR(tmp));
    solvearg = make_spils_solve_arg(r, gamma, delta, lr);
    zv = WRAP_NVECTOR(z);

    rv = caml_callback3_exn(*(data->closure_presolvefn), arg, solvearg, zv);

    relinquish_jac_arg(arg, 0);
    relinquish_spils_solve_arg(solvearg);
    RELINQUISH_WRAPPEDNV(zv);

    CAMLreturn(check_exception(rv));
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
    CAMLlocal4(arg, varg, jvarg, r);

    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    arg = make_jac_arg(t, y, fy, WRAP_NVECTOR(tmp));
    varg = WRAP_NVECTOR(v);
    jvarg = WRAP_NVECTOR(Jv);

    r = caml_callback3_exn(*(data->closure_jactimesfn), arg, varg, jvarg);

    relinquish_jac_arg(arg, 0);
    RELINQUISH_WRAPPEDNV(varg);
    RELINQUISH_WRAPPEDNV(jvarg);

    CAMLreturn(check_exception(r));
}

CAMLprim value CVTYPE(wf_tolerances)(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
 
    int flag = CVodeWFtolerances(data->cvode_mem, errw);
    CHECK_FLAG("CVodeWFtolerances", flag);

    CAMLreturn0;
}

CAMLprim value CVTYPE(dls_enable_dense_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVDlsSetDenseJacFn(data->cvode_mem, jacfn);
    CHECK_FLAG("CVDlsSetDenseJacFn", flag);
    CAMLreturn0;
}

CAMLprim value CVTYPE(dls_enable_band_jac_fn)(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVDlsSetBandJacFn(data->cvode_mem, bandjacfn);
    CHECK_FLAG("CVDlsSetBandJacFn", flag);
    CAMLreturn0;
}

CAMLprim value CVTYPE(enable_preconditioner)(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVSpilsSetPreconditioner(data->cvode_mem,
	    presetupfn, presolvefn);
    CHECK_FLAG("CVSpilsSetPreconditioner", flag);
    CAMLreturn0;
}

CAMLprim value CVTYPE(enable_jac_times_vec_fn)(value vdata)
{
    CAMLparam1(vdata);
    CVODE_DATA_FROM_ML(data, vdata);
    int flag = CVSpilsSetJacTimesVecFn(data->cvode_mem, jactimesfn);
    CHECK_FLAG("CVSpilsSetJacTimesVecFn", flag);
    CAMLreturn0;
}

/* basic interface */

CAMLprim value CVTYPE(init)(value lmm, value iter, value initial,
	value num_roots, value t0)
{
    CAMLparam4(lmm, iter, initial, num_roots);
    CAMLlocal1(vdata);

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

    int flag;

    int lmm_c;
    switch (Int_val(lmm)) {
    case VARIANT_LMM_ADAMS:
	lmm_c = CV_ADAMS;
	break;

    case VARIANT_LMM_BDF:
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

    N_Vector initial_nv = NVECTORIZE_VAL(initial);

    void *cvode_mem = CVodeCreate(lmm_c, iter_c);

    vdata = ml_cvode_data_alloc(cvode_mem);
    CVODE_DATA_FROM_ML(data, vdata);

    data->neq = Caml_ba_array_val(initial)->dim[0];
    data->num_roots = Int_val(num_roots);

    if (data->cvode_mem == NULL) {
	free(data);
	caml_failwith("CVodeCreate returned NULL");
	CAMLreturn0;
    }

    flag = CVodeInit(data->cvode_mem, f, Double_val(t0), initial_nv);
    RELINQUISH_NVECTORIZEDVAL(initial_nv);
    ml_cvode_check_flag("CVodeInit", flag, data);

    if (data->num_roots > 0) {
	flag = CVodeRootInit(data->cvode_mem, data->num_roots, roots);
	ml_cvode_check_flag("CVodeRootInit", flag, data);
    }

    CVodeSetUserData(data->cvode_mem, (void *)data);

    // setup linear solvers (if necessary)
    if (iter_c == CV_NEWTON) {
	set_linear_solver(data->cvode_mem, Field(iter, 0), data->neq);
    }

    // default tolerances
    flag = CVodeSStolerances(data->cvode_mem, RCONST(1.0e-4), RCONST(1.0e-8));
    ml_cvode_check_flag("CVodeSStolerances", flag, data);

    CAMLreturn(vdata);
}

CAMLprim value CVTYPE(sv_tolerances)(value vdata, value reltol, value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    CVODE_DATA_FROM_ML(data, vdata);

    N_Vector atol_nv = NVECTORIZE_VAL(abstol);

    int flag = CVodeSVtolerances(data->cvode_mem, Double_val(reltol), atol_nv);
    RELINQUISH_NVECTORIZEDVAL(atol_nv);
    CHECK_FLAG("CVodeSVtolerances", flag);

    CAMLreturn0;
}

CAMLprim value CVTYPE(reinit)(value vdata, value t0, value y0)
{
    CAMLparam3(vdata, t0, y0);

    CVODE_DATA_FROM_ML(data, vdata);

    N_Vector y0_nv = NVECTORIZE_VAL(y0);
    int flag = CVodeReInit(data->cvode_mem, Double_val(t0), y0_nv);
    RELINQUISH_NVECTORIZEDVAL(y0_nv);
    CHECK_FLAG("CVodeReInit", flag);

    CAMLreturn0;
}

static value solver(value vdata, value nextt, value y, int onestep)
{
    CAMLparam2(vdata, nextt);
    CAMLlocal1(r);

    realtype t = 0.0;
    CVODE_DATA_FROM_ML(data, vdata);

    N_Vector y_nv = NVECTORIZE_VAL(y);

    // TODO:
    // The payload of y (a big array) must not be shifted by the Ocaml GC
    // during this function call, even though Caml will be reentered
    // through the callback f. Is this guaranteed?
    int flag = CVode(data->cvode_mem, Double_val(nextt), y_nv, &t,
		     onestep ? CV_ONE_STEP : CV_NORMAL);
    RELINQUISH_NVECTORIZEDVAL(y_nv);
    CHECK_FLAG("CVode", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, caml_copy_double(t));

    switch (flag) {
    case CV_ROOT_RETURN:
	Store_field(r, 1, Val_int(VARIANT_SOLVER_RESULT_ROOTSFOUND));
	break;

    case CV_TSTOP_RETURN:
	Store_field(r, 1, Val_int(VARIANT_SOLVER_RESULT_STOPTIMEREACHED));
	break;

    default:
	Store_field(r, 1, Val_int(VARIANT_SOLVER_RESULT_CONTINUE));
    }

    CAMLreturn(r);
}

CAMLprim value CVTYPE(normal)(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value CVTYPE(one_step)(value vdata, value nextt, value y)
{
    CAMLparam3(vdata, nextt, y);
    CAMLreturn(solver(vdata, nextt, y, 1));
}

CAMLprim value CVTYPE(get_dky)(value vdata, value vt, value vk, value vy)
{
    CAMLparam4(vdata, vt, vk, vy);

    CVODE_DATA_FROM_ML(data, vdata);
    N_Vector y_nv = NVECTORIZE_VAL(vy);

    int flag = CVodeGetDky(data->cvode_mem, Double_val(vt), Int_val(vk), y_nv);
    CHECK_FLAG("CVodeGetDky", flag);
    RELINQUISH_NVECTORIZEDVAL(y_nv);
    
    CAMLreturn0;
}

CAMLprim value CVTYPE(get_err_weights)(value vcvode_mem, value verrws)
{
    CAMLparam2(vcvode_mem, verrws);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    N_Vector errws_nv = NVECTORIZE_VAL(verrws);

    int flag = CVodeGetErrWeights(cvode_mem, errws_nv);
    RELINQUISH_NVECTORIZEDVAL(errws_nv);
    CHECK_FLAG("CVodeGetErrWeights", flag);

    CAMLreturn0;
}

CAMLprim value CVTYPE(get_est_local_errors)(value vcvode_mem, value vele)
{
    CAMLparam2(vcvode_mem, vele);

    CVODE_MEM_FROM_ML(cvode_mem, vcvode_mem);
    N_Vector ele_nv = NVECTORIZE_VAL(vele);

    int flag = CVodeGetEstLocalErrors(cvode_mem, ele_nv);
    RELINQUISH_NVECTORIZEDVAL(ele_nv);
    CHECK_FLAG("CVodeGetEstLocalErrors", flag);

    CAMLreturn0;
}

CAMLprim value CVTYPE(vmax_norm)(value u)
{
    CAMLparam0();
    CAMLlocal1(r);

    N_Vector u_nv = NVECTORIZE_VAL(u);
    r = caml_copy_double(N_VMaxNorm(u_nv));
    RELINQUISH_NVECTORIZEDVAL(u_nv);

    CAMLreturn(r);
}

