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

#include <cvodes/cvodes.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>
#include <sundials/sundials_nvector.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

/* linear solvers */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_direct.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_spils.h>

#if SUNDIALS_BLAS_LAPACK == 1
#include <cvodes/cvodes_lapack.h>
#endif

#include "dls_ml.h"
#include "spils_ml.h"
#include "sundials_ml.h"
#include "cvode_ml.h"
#include "cvodes_ml.h"
#include "nvector_ml.h"

/* Interface with nvectors */

CAMLprim value c_cvodes_alloc_nvector_array(value vn)
{
    CAMLparam1(vn);
    value r;
    int n = Int_val(vn);

    if (n > 0) {
	r = caml_alloc_tuple(n);
    } else {
	r = Val_int(0);
    }

    CAMLreturn(r);
}

// NB: Normally, we should worry about relinquishing the elements of vy
// after we are finished using them (so as not to block the GC), but we
// instead make the assumption that these elements come from 'within'
// Sundials and thus that they would anyway not be GC-ed.
static void wrap_to_nvector_table(int n, value vy, N_Vector *y)
{
    int i;
    for (i = 0; i < n; ++i) {
	Store_field(vy, i, NVEC_BACKLINK(y[i]));
    }
}

static N_Vector *nvector_table_to_array(value vtable)
{
    int ns = Wosize_val (vtable); /* vtable : nvector array */
    N_Vector *r = calloc(ns + 1, sizeof(N_Vector));
    int i;

    for (i=0; i < ns; ++i) {
	r[i] = NVEC_VAL(Field(vtable, i));
    }
    r[ns] = NULL;

    return r;
}

static void free_nvector_array(N_Vector *nvarr)
{
    free(nvarr);
}

static int check_exception(value session, value r)
{
    CAMLparam2(session, r);
    CAMLlocal1(exn);

    if (!Is_exception_result(r)) return 0;

    r = Extract_exception(r);

    /* Unrecoverable error.  Save the exception and return -1.  */
    exn = caml_alloc_small (1,0);
    Field (exn, 0) = r;
    Store_field (session, RECORD_CVODE_SESSION_EXN_TEMP, exn);
    CAMLreturnT (int, -1);
}

/* Callbacks */

static int quadrhsfn(realtype t, N_Vector y, N_Vector yQdot, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_quadrhsfn);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(yQdot);

    // the data payloads inside args[2] and args[3] are only valid during
    // this call, afterward that memory goes back to cvode. These bigarrays
    // must not be retained by closure_quadrhsfn! If it wants a permanent
    // copy, then it has to make it manually.
    r = Int_val (caml_callbackN(*call_quadrhsfn,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int sensrhsfn(int ns, realtype t, N_Vector y, N_Vector ydot,
		     N_Vector *ys, N_Vector *ysdot, void *user_data,
		     N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal3(session, sensext, r);
    CAMLlocalN(args, 7);
    value *backref = user_data;

    /* We need to dereference on the C side, so that we can get access to
     * the values used to pass arrays of nvectors. */
    WEAK_DEREF (session, *backref);
    sensext = CVODE_SENSEXT_FROM_ML(session);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(ydot);
    args[3] = CVODES_SENSARRAY1_FROM_EXT(sensext);
    args[4] = CVODES_SENSARRAY2_FROM_EXT(sensext);
    args[5] = NVEC_BACKLINK(tmp1);
    args[6] = NVEC_BACKLINK(tmp2);

    wrap_to_nvector_table(ns, args[3], ys);
    wrap_to_nvector_table(ns, args[4], ysdot);

    // The data payloads inside args[1..6] are only valid during this call,
    // afterward that memory goes back to cvode. These bigarrays must not be
    // retained by closure_quadrhsfn! If it wants a permanent copy, then it
    // has to make it manually.
    r = caml_callbackN_exn(CVODES_SENSRHSFN_FROM_EXT(sensext),
			   sizeof (args) / sizeof (*args),
                           args);

    CAMLreturnT(int, check_exception(session, r));
}

static int sensrhsfn1(int ns, realtype t, N_Vector y, N_Vector ydot,
		      int is, N_Vector ys, N_Vector ysdot, void *user_data,
		      N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocalN(args, 9);
    int r;
    value *backref = user_data;
    CAML_FN (call_sensrhsfn1);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(ydot);
    args[4] = Val_int(is);
    args[5] = NVEC_BACKLINK(ys);
    args[6] = NVEC_BACKLINK(ysdot);
    args[7] = NVEC_BACKLINK(tmp1);
    args[8] = NVEC_BACKLINK(tmp2);

    // The data payloads inside args[2..3, 5..7] are only valid during this call,
    // afterward that memory goes back to cvode. These bigarrays must not be
    // retained by closure_quadrhsfn! If it wants a permanent copy, then it
    // has to make it manually.
    r = Int_val (caml_callbackN(*call_sensrhsfn1,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int quadsensrhsfn(int ns, realtype t, N_Vector y, N_Vector *ys,
		         N_Vector yqdot, N_Vector *yqsdot, void *user_data,
		         N_Vector tmp1, N_Vector tmp2)
{
    CAMLparam0();
    CAMLlocal3(session, sensext, r);
    CAMLlocalN(args, 7);
    value *backref = user_data;

    /* We need to dereference on the C side, so that we can get access to
     * the values used to pass arrays of nvectors. */
    WEAK_DEREF (session, *backref);
    sensext = CVODE_SENSEXT_FROM_ML(session);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = CVODES_SENSARRAY1_FROM_EXT(sensext);
    args[3] = NVEC_BACKLINK(yqdot);
    args[4] = CVODES_SENSARRAY2_FROM_EXT(sensext);
    args[5] = NVEC_BACKLINK(tmp1);
    args[6] = NVEC_BACKLINK(tmp2);

    wrap_to_nvector_table(ns, args[2], ys);
    wrap_to_nvector_table(ns, args[4], yqsdot);

    // The data payloads inside args[2..7] are only valid during this call,
    // afterward that memory goes back to cvode. These bigarrays must not be
    // retained by closure_quadrhsfn! If it wants a permanent copy, then it
    // has to make it manually.
    r = caml_callbackN_exn(CVODES_QUADSENSRHSFN_FROM_EXT(sensext),
		           sizeof (args) / sizeof (*args),
                           args);

    CAMLreturnT(int, check_exception(session, r));
}

static int brhsfn(realtype t, N_Vector y, N_Vector yb, N_Vector ybdot,
		  void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    int r;
    value *backref = user_data;
    CAML_FN (call_brhsfn);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(yb);
    args[4] = NVEC_BACKLINK(ybdot);

    // The data payloads inside args[2..4] are only valid during this call,
    // afterward that memory goes back to cvode. These bigarrays must not be
    // retained by closure_quadrhsfn! If it wants a permanent copy, then it
    // has to make it manually.
    r = Int_val (caml_callbackN(*call_brhsfn,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int brhsfn1(realtype t, N_Vector y, N_Vector *ys, N_Vector yb,
		   N_Vector ybdot, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(session, sensext, r);
    CAMLlocalN(args, 5);
    value *backref = user_data;
    int ns;

    /* We need to dereference on the C side, so that we can get access to
     * the values used to pass arrays of nvectors. */
    WEAK_DEREF (session, *backref);
    sensext = CVODE_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(sensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = CVODES_BSENSARRAY_FROM_EXT(sensext);
    args[3] = NVEC_BACKLINK(yb);
    args[4] = NVEC_BACKLINK(ybdot);

    wrap_to_nvector_table(ns, args[2], ys);

    // The data payloads inside args[2..5] are only valid during this call,
    // afterward that memory goes back to cvode. These bigarrays must not be
    // retained by closure_quadrhsfn! If it wants a permanent copy, then it
    // has to make it manually.
    r = caml_callbackN_exn(CVODES_BRHSFN1_FROM_EXT(sensext),
                           sizeof (args) / sizeof (*args),
                           args);

    CAMLreturnT(int, check_exception(session, r));
}

static int bquadrhsfn(realtype t, N_Vector y, N_Vector yb, N_Vector qbdot,
		      void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    int r;
    value *backref = user_data;
    CAML_FN (call_bquadrhsfn);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(yb);
    args[4] = NVEC_BACKLINK(qbdot);

    // The data payloads inside args[2..4] are only valid during this call,
    // afterward that memory goes back to cvode. These bigarrays must not be
    // retained by closure_quadrhsfn! If it wants a permanent copy, then it
    // has to make it manually.
    r = Int_val (caml_callbackN(*call_bquadrhsfn,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int bquadrhsfn1(realtype t, N_Vector y, N_Vector *ys, N_Vector yb,
		       N_Vector qbdot, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    CAMLlocal3(session, sensext, r);
    int ns;
    value *backref = user_data;

    /* We need to dereference on the C side, so that we can get access to
     * the values used to pass arrays of nvectors. */
    WEAK_DEREF (session, *backref);
    sensext = CVODE_SENSEXT_FROM_ML(session);
    ns = Int_val(Field(sensext, RECORD_CVODES_BWD_SESSION_NUMSENSITIVITIES));

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = CVODES_BSENSARRAY_FROM_EXT(sensext);
    args[3] = NVEC_BACKLINK(yb);
    args[4] = NVEC_BACKLINK(qbdot);

    wrap_to_nvector_table(ns, args[2], ys);

    // The data payloads inside args[2..5] are only valid during this call,
    // afterward that memory goes back to cvode. These bigarrays must not be
    // retained by closure_quadrhsfn! If it wants a permanent copy, then it
    // has to make it manually.
    r = caml_callbackN(CVODES_BQUADRHSFN1_FROM_EXT(sensext),
                       sizeof (args) / sizeof (*args),
                       args);

    CAMLreturnT(int, check_exception(session, r));
}

static value make_jac_arg(realtype t, N_Vector y, N_Vector yb,
			   N_Vector fyb, value tmp)
{
    CAMLparam1(tmp);
    CAMLlocal1(r);

    r = caml_alloc_tuple(RECORD_CVODES_ADJ_JACOBIAN_ARG_SIZE);
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_T, caml_copy_double(t));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_Y, NVEC_BACKLINK(y));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_YB, NVEC_BACKLINK(yb));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_FYB, NVEC_BACKLINK(fyb));
    Store_field(r, RECORD_CVODES_ADJ_JACOBIAN_ARG_JAC_TMP, tmp);

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

static value make_spils_solve_arg(
	N_Vector rvecb,
	realtype gammab,
	realtype deltab,
	int lrb)

{
    CAMLparam0();
    CAMLlocal1(v);

    v = caml_alloc_tuple(RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_SIZE);
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_RVEC, NVEC_BACKLINK(rvecb));
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_GAMMA,
                caml_copy_double(gammab));
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_DELTA,
                caml_copy_double(deltab));
    Store_field(v, RECORD_CVODES_ADJ_SPILS_SOLVE_ARG_LR, Val_bool(lrb=1));

    CAMLreturn(v);
}

static int bprecsolvefn(
	realtype t,
	N_Vector y,
	N_Vector yb,
	N_Vector fyb,
	N_Vector rvecb,
	N_Vector zvecb,
	realtype gammab,
	realtype deltab,
	int lrb,
	void *user_data,
	N_Vector tmpb)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_bprecsolvefn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, yb, fyb, NVEC_BACKLINK(tmpb));
    args[2] = make_spils_solve_arg(rvecb, gammab, deltab, lrb);
    args[3] = NVEC_BACKLINK(zvecb);

    retcode = Int_val (caml_callbackN(*call_bprecsolvefn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    CAMLreturnT(int, retcode);
}

static int bprecsetupfn(
    realtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    booleantype jokb,
    booleantype *jcurPtrB,
    realtype gammab,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocal2(session, vr);
    CAMLlocalN(args, 4);
    value *backref = user_data;
    CAML_FN (call_bprecsetupfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, yb, fyb, make_triple_tmp(tmp1b, tmp2b, tmp3b));
    args[2] = Val_bool(jokb);
    args[3] = caml_copy_double(gammab);

    vr = caml_callbackN(*call_bprecsetupfn,
			sizeof (args) / sizeof (*args),
			args);
    *jcurPtrB = Bool_val(Field(vr, 0));

    CAMLreturnT(int, Int_val(Field(vr, 1)));
}

static int bjactimesfn(
    N_Vector vb,
    N_Vector Jvb,
    realtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    void *user_data,
    N_Vector tmpb)
{
    CAMLparam0();
    CAMLlocal1(r);
    CAMLlocalN(args, 4);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_bjactimesfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, yb, fyb, NVEC_BACKLINK(tmpb));
    args[2] = NVEC_BACKLINK(vb);
    args[3] = NVEC_BACKLINK(Jvb);

    retcode = Int_val (caml_callbackN(*call_bjactimesfn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    CAMLreturnT(int, retcode);
}

static int bjacfn(
    long int neqb,
    realtype t,
    N_Vector y,
    N_Vector yb,
    N_Vector fyb,
    DlsMat jacb,
    void *user_data,
    N_Vector tmp1b,
    N_Vector tmp2b,
    N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    int retcode;
    value *backref = user_data;
    CAML_FN (call_bjacfn);

    args[0] = *backref;
    args[1] = make_jac_arg(t, y, yb, fyb, make_triple_tmp(tmp1b, tmp2b, tmp3b));
    args[2] = caml_alloc_final (2, NULL, 0, 1);
    DLSMAT(args[2]) = jacb;

    retcode = Int_val (caml_callbackN(*call_bjacfn,
                                      sizeof (args) / sizeof (*args),
                                      args));

    CAMLreturnT(int, retcode);
}

static int bbandjacfn(
	long int nb,
	long int mupperb,
	long int mlowerb,
	realtype t,
	N_Vector y,
	N_Vector yb,
	N_Vector fyb,
	DlsMat jacb,
	void *user_data, 	 
	N_Vector tmp1b,
	N_Vector tmp2b,
	N_Vector tmp3b)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_bbandjacfn);

    args[0] = *backref;
    args[1] = caml_alloc_tuple(RECORD_CVODES_ADJ_BANDRANGE_SIZE);
    Store_field(args[1], RECORD_CVODES_ADJ_BANDRANGE_MUPPER, Val_long(mupperb));
    Store_field(args[1], RECORD_CVODES_ADJ_BANDRANGE_MLOWER, Val_long(mlowerb));
    args[2] = make_jac_arg(t, y, yb, fyb, make_triple_tmp(tmp1b, tmp2b, tmp3b));
    args[3] = caml_alloc_final(2, NULL, 0, 1);
    DLSMAT(args[3]) = jacb;

    r = Int_val (caml_callbackN(*call_bbandjacfn,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

/* quadrature interface */

CAMLprim void c_cvodes_quad_init(value vdata, value vq0)
{
    CAMLparam2(vdata, vq0);
    CAMLlocal1(r);
    int flag;
    N_Vector q0 = NVEC_VAL(vq0);
    
    flag = CVodeQuadInit(CVODE_MEM_FROM_ML(vdata), quadrhsfn, q0);
    SCHECK_FLAG("CVodeQuadInit", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quad_reinit(value vdata, value vq0)
{
    CAMLparam2(vdata, vq0);
    CAMLlocal1(r);
    int flag;
    N_Vector q0 = NVEC_VAL(vq0);
    
    flag = CVodeQuadReInit(CVODE_MEM_FROM_ML(vdata), q0);
    SCHECK_FLAG("CVodeQuadReInit", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quad_sv_tolerances(value vdata, value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    N_Vector atol_nv = NVEC_VAL(abstol);

    int flag = CVodeQuadSVtolerances(CVODE_MEM_FROM_ML(vdata),
	    Double_val(reltol), atol_nv);
    SCHECK_FLAG("CVodeQuadSVtolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_quad_get(value vdata, value vyq)
{
    CAMLparam2(vdata, vyq);
    N_Vector yq = NVEC_VAL(vyq);
    realtype tret;

    int flag = CVodeGetQuad(CVODE_MEM_FROM_ML(vdata), &tret, yq);
    SCHECK_FLAG("CVodeGetQuad", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim void c_cvodes_quad_get_dky(value vdata, value vt, value vk,
				   value vdkyq)
{
    CAMLparam4(vdata, vt, vk, vdkyq);
    N_Vector dkyq = NVEC_VAL(vdkyq);

    int flag = CVodeGetQuadDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkyq);
	    
    SCHECK_FLAG("CVodeGetQuadDky", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quad_get_err_weights(value vdata, value veqweight)
{
    CAMLparam2(vdata, veqweight);
    N_Vector eqweight = NVEC_VAL(veqweight);

    int flag = CVodeGetQuadErrWeights(CVODE_MEM_FROM_ML(vdata), eqweight);
    SCHECK_FLAG("CVodeGetQuadErrWeights", flag);

    CAMLreturn0;
}

/* sensitivity interface */

CAMLprim void c_cvodes_sens_sv_tolerances(value vdata, value reltol,
					 value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = nvector_table_to_array(abstol);

    int flag = CVodeSensSVtolerances(CVODE_MEM_FROM_ML(vdata),
	    Double_val(reltol), atol_nv);
    free_nvector_array(atol_nv); 
    SCHECK_FLAG("CVodeSensSVtolerances", flag);

    CAMLreturn0;
}

static int decode_sens_method(value vmethod)
{
    switch (Int_val(vmethod)) {
    case VARIANT_CVODES_SENS_METHOD_SIMULTANEOUS:
	return CV_SIMULTANEOUS;

    case VARIANT_CVODES_SENS_METHOD_STAGGERED:
	return CV_STAGGERED;

    case VARIANT_CVODES_SENS_METHOD_STAGGERED1:
	return CV_STAGGERED1;

    default:
	caml_failwith("Illegal sens method.");
    }
}

CAMLprim void c_cvodes_sens_init(value vdata, value vmethod, value vrhsfn,
				value vys0)
{
    CAMLparam4(vdata, vmethod, vrhsfn, vys0);
    int ns = Wosize_val (vys0); /* vys0 : nvector array */
    N_Vector *ys0 = nvector_table_to_array(vys0);

    int flag = CVodeSensInit(CVODE_MEM_FROM_ML(vdata), ns,
			     decode_sens_method(vmethod),
			     ((Bool_val(vrhsfn)) ? sensrhsfn : NULL),
			     ys0);
    free_nvector_array(ys0); 
    SCHECK_FLAG("CVodeSensInit", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_init_1(value vdata, value vmethod, value vrhsfn,
				  value vys0)
{
    CAMLparam4(vdata, vmethod, vrhsfn, vys0);
    int ns = Wosize_val (vys0); /* vys0 : nvector array */
    N_Vector *ys0 = nvector_table_to_array(vys0);

    int flag = CVodeSensInit1(CVODE_MEM_FROM_ML(vdata), ns,
			      decode_sens_method(vmethod),
			      ((Bool_val(vrhsfn)) ? sensrhsfn1 : NULL),
			      ys0);
    free_nvector_array(ys0); 
    SCHECK_FLAG("CVodeSensInit", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_reinit(value vdata, value vmethod, value vs0)
{
    CAMLparam3(vdata, vmethod, vs0);
    CAMLlocal1(r);
    int flag;
    N_Vector *s0 = nvector_table_to_array(vs0);
    
    flag = CVodeSensReInit(CVODE_MEM_FROM_ML(vdata),
			   decode_sens_method(vmethod),
			   s0);
    free_nvector_array(s0);
    SCHECK_FLAG("CVodeSensReInit", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_sens_get(value vdata, value vys)
{
    CAMLparam2(vdata, vys);
    N_Vector *ys = nvector_table_to_array(vys);
    realtype tret;

    int flag = CVodeGetSens(CVODE_MEM_FROM_ML(vdata), &tret, ys);
    free_nvector_array(ys);
    SCHECK_FLAG("CVodeGetSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim void c_cvodes_sens_get_dky(value vdata, value vt, value vk, value vdkys)
{
    CAMLparam4(vdata, vt, vk, vdkys);
    N_Vector *dkys = nvector_table_to_array(vdkys);

    int flag = CVodeGetSensDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), dkys);
    free_nvector_array(dkys);
    SCHECK_FLAG("CVodeGetSensDky", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_sens_get1(value vdata, value vis, value vys)
{
    CAMLparam3(vdata, vis, vys);
    N_Vector ys = NVEC_VAL(vys);
    realtype tret;

    int flag = CVodeGetSens1(CVODE_MEM_FROM_ML(vdata), &tret, Int_val(vis), ys);
    SCHECK_FLAG("CVodeGetSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim void c_cvodes_sens_get_dky1(value vdata, value vt, value vk,
				     value vis, value vdkys)
{
    CAMLparam5(vdata, vt, vk, vis, vdkys);
    N_Vector dkys = NVEC_VAL(vdkys);

    int flag = CVodeGetSensDky1(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
	    Int_val(vk), Int_val(vis), dkys);
	    
    SCHECK_FLAG("CVodeGetSensDky1", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_get_err_weights(value vdata, value vesweight)
{
    CAMLparam2(vdata, vesweight);
    N_Vector *esweight = nvector_table_to_array(vesweight);

    int flag = CVodeGetSensErrWeights(CVODE_MEM_FROM_ML(vdata), esweight);
    free_nvector_array(esweight);
    SCHECK_FLAG("CVodeGetSensErrWeights", flag);

    CAMLreturn0;
}

/* sensitivity/quadrature interface */

CAMLprim void c_cvodes_quadsens_init(value vdata, value vrhsfn, value vyqs0)
{
    CAMLparam3(vdata, vrhsfn, vyqs0);
    N_Vector *yqs0 = nvector_table_to_array(vyqs0);

    int flag = CVodeQuadSensInit(CVODE_MEM_FROM_ML(vdata),
				 Bool_val (vrhsfn) ? quadsensrhsfn : NULL,
				 yqs0);
    free_nvector_array(yqs0); 
    SCHECK_FLAG("CVodeQuadSensInit", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_reinit(value vdata, value vyqs0)
{
    CAMLparam2(vdata, vyqs0);
    N_Vector *yqs0 = nvector_table_to_array(vyqs0);

    int flag = CVodeQuadSensReInit(CVODE_MEM_FROM_ML(vdata), yqs0);
    free_nvector_array(yqs0); 
    SCHECK_FLAG("CVodeQuadSensReInit", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_sv_tolerances(value vdata, value reltol,
					     value abstol)
{
    CAMLparam3(vdata, reltol, abstol);
    N_Vector *atol_nv = nvector_table_to_array(abstol);

    int flag = CVodeQuadSensSVtolerances(CVODE_MEM_FROM_ML(vdata),
	    Double_val(reltol), atol_nv);
    free_nvector_array(atol_nv); 
    SCHECK_FLAG("CVodeQuadSensSVtolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_quadsens_get(value vdata, value vyqs)
{
    CAMLparam2(vdata, vyqs);
    N_Vector *yqs = nvector_table_to_array(vyqs);
    realtype tret;

    int flag = CVodeGetQuadSens(CVODE_MEM_FROM_ML(vdata), &tret, yqs);
    free_nvector_array(yqs); 
    SCHECK_FLAG("CVodeGetQuadSens", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim value c_cvodes_quadsens_get1(value vdata, value vis, value vyqs)
{
    CAMLparam3(vdata, vis, vyqs);
    N_Vector yqs = NVEC_VAL(vyqs);
    realtype tret;

    int flag = CVodeGetQuadSens1(CVODE_MEM_FROM_ML(vdata), &tret,
			         Int_val(vis), yqs);
    SCHECK_FLAG("CVodeGetQuadSens1", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim void c_cvodes_quadsens_get_dky(value vdata, value vt, value vk,
				       value vdkyqs)
{
    CAMLparam4(vdata, vt, vk, vdkyqs);
    N_Vector *dkyqs = nvector_table_to_array(vdkyqs);

    int flag = CVodeGetQuadSensDky(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
				   Int_val(vk), dkyqs);
    free_nvector_array(dkyqs); 
    SCHECK_FLAG("CVodeGetQuadSensDky", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_get_dky1(value vdata, value vt, value vk,
					value vis, value vdkyqs)
{
    CAMLparam5(vdata, vt, vk, vis, vdkyqs);
    N_Vector dkyqs = NVEC_VAL(vdkyqs);

    int flag = CVodeGetQuadSensDky1(CVODE_MEM_FROM_ML(vdata), Double_val(vt),
				    Int_val(vk), Int_val(vis), dkyqs);
    SCHECK_FLAG("CVodeGetQuadSensDky1", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_get_err_weights(value vdata, value veqweights)
{
    CAMLparam2(vdata, veqweights);
    N_Vector *eqweights = nvector_table_to_array(veqweights);

    int flag = CVodeGetQuadSensErrWeights(CVODE_MEM_FROM_ML(vdata), eqweights);
    free_nvector_array(eqweights); 
    SCHECK_FLAG("CVodeGetQuadSensErrWeights", flag);

    CAMLreturn0;
}

/* adjoint interface */

static value forward_solver(value vdata, value vtout, value vyret, int onestep)
{
    CAMLparam3(vdata, vtout, vyret);
    CAMLlocal1(ret);
    N_Vector yret = NVEC_VAL(vyret);
    realtype tret;
    int ncheck;
    enum cvode_solver_result_tag solver_result;

    int flag = CVodeF(CVODE_MEM_FROM_ML(vdata), Double_val(vtout), yret,
		      &tret, onestep ? CV_ONE_STEP : CV_NORMAL, &ncheck);
    switch (flag) {
    case CV_SUCCESS:
	solver_result = VARIANT_CVODE_SOLVER_RESULT_CONTINUE;
	break;

    case CV_TSTOP_RETURN:
	solver_result = VARIANT_CVODE_SOLVER_RESULT_STOPTIMEREACHED;
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
	SCHECK_FLAG ("CVodeF", flag);
    }

    /* Hmm...should this go in the production code or not?  */
    if (Is_block (Field (vdata, RECORD_CVODE_SESSION_EXN_TEMP)))
	abort ();

    ret = caml_alloc_tuple(3);
    Store_field(ret, 0, caml_copy_double(tret));
    Store_field(ret, 1, Val_int(ncheck));
    Store_field(ret, 2, Val_int(solver_result));

    CAMLreturn(ret);
}

CAMLprim value c_cvodes_adj_forward_normal(value vdata, value vtout, value vyret)
{
    CAMLparam3(vdata, vtout, vyret);
    CAMLreturn(forward_solver(vdata, vtout, vyret, 0));
}

CAMLprim value c_cvodes_adj_forward_one_step(value vdata, value vtout,
					    value vyret)
{
    CAMLparam3(vdata, vtout, vyret);
    CAMLreturn(forward_solver(vdata, vtout, vyret, 1));
}

CAMLprim void c_cvodes_adj_sv_tolerances(value vparent, value vwhich,
					value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = CVodeSVtolerancesB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				  Double_val(vreltol), atol_nv);
    SCHECK_FLAG("CVodeSStolerancesB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_set_preconditioner(value vparent,
						   value vwhich,
						   value vset_precsetup,
						   value vset_jac)
{
    CAMLparam4(vparent, vwhich, vset_precsetup, vset_jac);
    int flag;
    void *mem = CVODE_MEM_FROM_ML(vparent);
    int which = Int_val(vwhich);
    CVSpilsPrecSetupFnB bsetup = Bool_val(vset_precsetup) ? bprecsetupfn : NULL;

    flag = CVSpilsSetPreconditionerB(mem, which, bsetup, bprecsolvefn);
    SCHECK_FLAG ("CVSpilsSetPreconditionerB", flag);
    if (Bool_val(vset_jac)) {
	flag = CVSpilsSetJacTimesVecFnB(mem, which, bjactimesfn);
	SCHECK_FLAG ("CVSpilsSetJacTimesVecFnB", flag);
    }

    CAMLreturn0;
}

/* Dense and Band can only be used with serial NVectors.  */
CAMLprim void c_cvodes_adj_dls_dense(value vparent, value vwhich,
				    value vnb, value vset_jac)
{
    CAMLparam3(vparent, vwhich, vset_jac);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    long nbeqs = Long_val(vnb);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVDenseB (cvode_mem, which, nbeqs);
    SCHECK_FLAG ("CVDenseB", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetDenseJacFnB(cvode_mem, which, bjacfn);
	SCHECK_FLAG("CVDlsSetDenseJacFnB", flag);
    }
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_dls_lapack_dense(value vparent, value vwhich,
					   value vnb, value vset_jac)
{
    CAMLparam3 (vparent, vwhich, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    long nbeqs = Long_val(vnb);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVLapackDenseB (cvode_mem, which, nbeqs);
    SCHECK_FLAG ("CVLapackDenseB", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetDenseJacFnB (cvode_mem, which, bjacfn);
	SCHECK_FLAG("CVDlsSetDenseJacFnB", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_dls_set_dense_jac_fn(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);
    int flag = CVDlsSetDenseJacFnB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				   bjacfn);
    SCHECK_FLAG("CVDlsSetDenseJacFnB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_dls_clear_dense_jac_fn(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);
    int flag = CVDlsSetDenseJacFnB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				   NULL);
    SCHECK_FLAG("CVDlsSetDenseJacFnB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_dls_band (value vparent_which, value vnb,
				    value vmupper, value vmlower,
				    value vset_jac)
{
    CAMLparam5(vparent_which, vnb, vmupper, vmlower, vset_jac);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(vnb);
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVBandB (cvode_mem, which, nbeqs,
		    Long_val (vmupper), Long_val (vmlower));
    SCHECK_FLAG ("CVBandB", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetBandJacFnB(cvode_mem, which, bbandjacfn);
	SCHECK_FLAG("CVDlsSetBandJacFnB", flag);
    }
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_dls_lapack_band (value vparent_which, value vnb,
					   value vmupper, value vmlower,
					   value vset_jac)
{
    CAMLparam5(vparent_which, vnb, vmupper, vmlower, vset_jac);
#if SUNDIALS_BLAS_LAPACK
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which, 0));
    long nbeqs = Long_val(vnb);
    int which = Int_val(Field(vparent_which, 1));
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVLapackBandB (cvode_mem, which, nbeqs,
			  Long_val (vmupper), Long_val (vmlower));
    SCHECK_FLAG ("CVLapackBandB", flag);
    if (Bool_val (vset_jac)) {
	flag = CVDlsSetBandJacFnB(cvode_mem, which, bbandjacfn);
	SCHECK_FLAG("CVDlsSetBandJacFnB", flag);
    }
#else
    caml_failwith("Lapack solvers are not available.");
#endif
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_dls_set_band_jac_fn(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);
    int flag = CVDlsSetBandJacFnB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				  bbandjacfn);
    SCHECK_FLAG("CVDlsSetBandJacFnB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_dls_clear_band_jac_fn(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);
    int flag = CVDlsSetBandJacFnB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				  NULL);
    SCHECK_FLAG("CVDlsSetBandJacFnB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_banded_spgmr(value vparent_which_vnb,
					      value vmupper, value vmlower,
					      value vmaxl, value vtype)
{
    CAMLparam5 (vparent_which_vnb, vmupper, vmlower, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which_vnb, 0));
    int which = Int_val(Field(vparent_which_vnb, 1));
    long neqs = Long_val(Field(vparent_which_vnb, 2));
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpgmrB (cvode_mem, which, spils_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpgmrB", flag);
    flag = CVBandPrecInitB (cvode_mem, which, neqs,
			    Long_val (vmupper), Long_val (vmlower));
    SCHECK_FLAG ("CVBandPrecInitB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_banded_spbcg(value vparent_which_vnb,
					      value vmupper, value vmlower,
					      value vmaxl, value vtype)
{
    CAMLparam5 (vparent_which_vnb, vmupper, vmlower, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which_vnb, 0));
    int which = Int_val(Field(vparent_which_vnb, 1));
    long nbeqs = Long_val(Field(vparent_which_vnb, 2));
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpbcgB (cvode_mem, which, spils_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpbcgB", flag);
    flag = CVBandPrecInitB (cvode_mem, which, nbeqs,
			    Long_val (vmupper), Long_val (vmlower));
    SCHECK_FLAG ("CVBandPrecInitB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_banded_sptfqmr(value vparent_which_vnb,
					        value vmupper, value vmlower,
					        value vmaxl, value vtype)
{
    CAMLparam5 (vparent_which_vnb, vmupper, vmlower, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparent_which_vnb, 0));
    int which = Int_val(Field(vparent_which_vnb, 1));
    long nbeqs = Long_val(Field(vparent_which_vnb, 2));
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSptfqmrB (cvode_mem, which, spils_precond_type (vtype),
		       Int_val (vmaxl));
    SCHECK_FLAG ("CVSptfqmrB", flag);
    flag = CVBandPrecInitB (cvode_mem, which, nbeqs,
			    Long_val (vmupper), Long_val (vmlower));
    SCHECK_FLAG ("CVBandPrecInitB", flag);
    CAMLreturn0;
}

CAMLprim value c_cvodes_adj_init_backward(value vparent, value weakref,
					 value vargs, value vwithsens)
{
    CAMLparam4(vparent, weakref, vargs, vwithsens);
    CAMLlocal1(r);
    CAMLlocal2(vlmm, viter);
    int flag, lmm_c, iter_c, which;
    void *parent = CVODE_MEM_FROM_ML(vparent);

    int lmm = Field(vargs, 0);
    int iter = Field(vargs, 1);
    realtype tb0 = Double_val(Field(vargs, 2));
    N_Vector initial_nv = NVEC_VAL(Field(vargs, 3));

    switch (Int_val(lmm)) {
    case VARIANT_CVODE_LMM_ADAMS:
	lmm_c = CV_ADAMS;
	break;

    case VARIANT_CVODE_LMM_BDF:
	lmm_c = CV_BDF;
	break;

    default:
	caml_failwith("Illegal lmm value.");
    }

    if (Is_block(iter)) {
	iter_c = CV_NEWTON;
    } else {
	iter_c = CV_FUNCTIONAL;
    }

    flag = CVodeCreateB(parent, lmm_c, iter_c, &which);
    if (flag != CV_SUCCESS) {
	SCHECK_FLAG("CVodeCreateB", flag);
    }

    if (Bool_val(vwithsens)) {
	flag = CVodeInitBS(parent, which, brhsfn1, tb0, initial_nv);
	if (flag != CV_SUCCESS) {
	    SCHECK_FLAG("CVodeInitBS", flag);
	}
    } else {
	flag = CVodeInitB(parent, which, brhsfn, tb0, initial_nv);
	if (flag != CV_SUCCESS) {
	    SCHECK_FLAG("CVodeInitB", flag);
	}
    }

    value *backref;
    backref = malloc (sizeof (*backref));
    if (backref == NULL) {
	caml_raise_out_of_memory();
    }
    *backref = weakref;
    caml_register_generational_global_root (backref);
    CVodeSetUserDataB (parent, which, backref);

    r = caml_alloc_tuple (4);
    Store_field (r, 0, (value)CVodeGetAdjCVodeBmem(parent, which));
    Store_field (r, 1, Val_int(which));
    Store_field (r, 2, (value)backref);
    Store_field (r, 3, Val_long (0)); /* no err_file = NULL */

    CAMLreturn(r);
}

CAMLprim void c_cvodes_adj_reinit(value vparent, value vwhich,
				 value vtb0, value vyb0)
{
    CAMLparam4(vparent, vwhich, vtb0, vyb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yb0 = NVEC_VAL(vyb0);
    
    flag = CVodeReInitB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			Double_val(vtb0), yb0);
    SCHECK_FLAG("CVodeReInitB", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_adj_get(value vparent, value vwhich, value vyb)
{
    CAMLparam3(vparent, vwhich, vyb);
    N_Vector yq = NVEC_VAL(vyb);
    realtype tret;

    int flag = CVodeGetB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			 &tret, yq);
    SCHECK_FLAG("CVodeGetB", flag);

    CAMLreturn(caml_copy_double(tret));
}

/* adjoint/quadrature interface */

CAMLprim void c_cvodes_adjquad_initb(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = CVodeQuadInitB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			  bquadrhsfn, yqb0);
    SCHECK_FLAG("CVodeQuadInitB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adjquad_initbs(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = CVodeQuadInitBS(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			   bquadrhsfn1, yqb0);
    SCHECK_FLAG("CVodeQuadInitBS", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adjquad_reinit(value vparent, value vwhich, value vyqb0)
{
    CAMLparam3(vparent, vwhich, vyqb0);
    CAMLlocal1(r);
    int flag;
    N_Vector yqb0 = NVEC_VAL(vyqb0);
    
    flag = CVodeQuadReInitB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich), yqb0);
    SCHECK_FLAG("CVodeQuadReInitB", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_adjquad_get(value vparent, value vwhich, value vyqb)
{
    CAMLparam3(vparent, vwhich, vyqb);
    CAMLlocal1(r);
    N_Vector yqb = NVEC_VAL(vyqb);
    realtype tret;

    int flag = CVodeGetQuadB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			     &tret, yqb);
    SCHECK_FLAG("CVodeGetQuadB", flag);

    CAMLreturn(caml_copy_double(tret));
}

CAMLprim void c_cvodes_adjquad_sv_tolerances(value vparent, value vwhich,
					    value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);
    N_Vector atol_nv = NVEC_VAL(vabstol);

    int flag = CVodeQuadSVtolerancesB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				      Double_val(vreltol), atol_nv);
    SCHECK_FLAG("CVodeQuadSVtolerancesB", flag);

    CAMLreturn0;
}

/* Interface without nvectors */

#define MAX_ERRMSG_LEN 256

void cvodes_ml_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	    || flag == CV_ROOT_RETURN
	    || flag == CV_TSTOP_RETURN
	    || flag == CV_WARNING) return;

    switch (flag) {
	case CV_TOO_MUCH_WORK:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchWork"));

	case CV_TOO_MUCH_ACC:
	    caml_raise_constant(*caml_named_value("cvode_TooMuchAccuracy"));

	case CV_ERR_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ErrFailure"));

	case CV_CONV_FAILURE:
	    caml_raise_constant(*caml_named_value("cvode_ConvergenceFailure"));

	/* * */

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

	/* * */

	case CV_ILL_INPUT:
	    caml_raise_constant(*caml_named_value("cvode_IllInput"));
	
	case CV_BAD_K:
	    caml_raise_constant(*caml_named_value("cvode_BadK"));

	case CV_BAD_T:
	    caml_raise_constant(*caml_named_value("cvode_BadT"));

	case CV_BAD_DKY:
	    caml_raise_constant(*caml_named_value("cvode_BadDky"));

	case CV_TOO_CLOSE:
	    caml_raise_constant(*caml_named_value("cvode_TooClose"));

	/* Quadrature */

        case CV_NO_QUAD:
	    caml_raise_constant(*caml_named_value("cvodes_QuadNotInitialized"));

        case CV_QRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_QuadRhsFuncFailure"));

        case CV_FIRST_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstQuadRhsFuncErr"));

        case CV_REPTD_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedQuadRhsFuncErr"));

        case CV_UNREC_QRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableQuadRhsFuncErr"));

	/* Sensitivity */

        case CV_NO_SENS:
	    caml_raise_constant(*caml_named_value("cvodes_SensNotInitialized"));

        case CV_SRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_SensRhsFuncFailure"));

        case CV_FIRST_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstSensRhsFuncErr"));

        case CV_REPTD_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedSensRhsFuncErr"));

        case CV_UNREC_SRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableSensRhsFuncErr"));

        case CV_BAD_IS:
	    caml_raise_constant(*caml_named_value("cvodes_BadIS"));

	/* Sensitivity > Quadrature */

        case CV_NO_QUADSENS:
	    caml_raise_constant(*caml_named_value("cvodes_QuadSensNotInitialized"));

        case CV_QSRHSFUNC_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_QuadSensRhsFuncFailure"));

        case CV_FIRST_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_FirstQuadSensRhsFuncErr"));

        case CV_REPTD_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_RepeatedQuadSensRhsFuncErr"));

        case CV_UNREC_QSRHSFUNC_ERR:
	    caml_raise_constant(*caml_named_value("cvodes_UnrecoverableQuadSensRhsFuncErr"));

	/* Adjoint */

        case CV_NO_ADJ:
	    caml_raise_constant(*caml_named_value("cvodes_AdjointNotInitialized"));

        case CV_NO_FWD:
	    caml_raise_constant(*caml_named_value("cvodes_NoForwardCall"));

        case CV_NO_BCK:
	    caml_raise_constant(*caml_named_value("cvodes_NoBackwardProblem"));

        case CV_BAD_TB0:
	    caml_raise_constant(*caml_named_value("cvodes_BadFinalTime"));

        case CV_REIFWD_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_ForwardReinitializationFailed"));

        case CV_FWD_FAIL:
	    caml_raise_constant(*caml_named_value("cvodes_ForwardFailed"));

	default:
	    /* e.g. CVDIAG_MEM_NULL, CVDIAG_ILL_INPUT, CVDIAG_MEM_FAIL */
	    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
		    CVodeGetReturnFlagName(flag));
	    caml_failwith(exmsg);
    }
}

/* basic quadrature interface */

CAMLprim void c_cvodes_quad_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = CVodeSetQuadErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadErrCon", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quad_ss_tolerances(value vdata,
					  value reltol,
					  value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeQuadSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), Double_val(abstol));
    SCHECK_FLAG("CVodeQuadSStolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_quad_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quad_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quad_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfqevals;
    long int nqetfails;

    flag = CVodeGetQuadStats(CVODE_MEM_FROM_ML(vdata), &nfqevals,
	    &nqetfails);
    SCHECK_FLAG("CVodeGetQuadStats", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(nfqevals));
    Store_field(r, 1, Val_long(nqetfails));

    CAMLreturn(r);
}

/* sensitivity interface */

CAMLprim void c_cvodes_sens_set_err_con(value vdata, value verrcons)
{
    CAMLparam2(vdata, verrcons);
    int flag;
    
    flag = CVodeSetSensErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrcons));
    SCHECK_FLAG("CVodeSetSensErrCon", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_ss_tolerances(value vdata,
					  value reltol,
					  value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeSensSStolerances(CVODE_MEM_FROM_ML(vdata),
		 Double_val(reltol), REAL_ARRAY(abstol));
    SCHECK_FLAG("CVodeSensSStolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSensEEtolerances(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSensEEtolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_set_params(value vdata, value vparams)
{
    CAMLparam2(vdata, vparams);
    CAMLlocal3(vp, vpbar, vplist);

    realtype *p = NULL;
    realtype *pbar = NULL;
    int *plist = NULL;
    int i, ns;

    vp = Field(vparams, RECORD_CVODES_SENS_PARAMS_PVALS);
    vpbar = Field(vparams, RECORD_CVODES_SENS_PARAMS_PBAR);
    vplist = Field(vparams, RECORD_CVODES_SENS_PARAMS_PLIST);

    if (vp != Val_none) p = REAL_ARRAY(Some_val(vp));
    if (vpbar != Val_none) pbar = REAL_ARRAY(Some_val(vpbar));

    if (vplist != Val_none) {
	vplist = Some_val(vplist);
	ns = Wosize_val (vplist); /* vplist : int array */
	plist = calloc(ns, sizeof(int));

	for (i=0; i < ns; ++i) {
	    plist[i] = Int_val(Field(vplist, i));
	}
    }

    int flag = CVodeSetSensParams(CVODE_MEM_FROM_ML(vdata), p, pbar, plist);
    if (plist != NULL) free(plist);
    SCHECK_FLAG("CVodeSetSensParams", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_toggle_off(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSensToggleOff(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSensToggleOff", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_set_dq_method(value vdata, value vdqtype,
					  value vdqrhomax)
{
    CAMLparam3(vdata, vdqtype, vdqrhomax);
    int dqtype;

    switch (Int_val(vdqtype)) {
    case VARIANT_CVODES_SENS_DQ_METHOD_CENTERED:
	dqtype = CV_CENTERED;
	break;

    case VARIANT_CVODES_SENS_DQ_METHOD_FORWARD:
	dqtype = CV_FORWARD;
	break;

    default:
	caml_failwith("Illegal dq method.");
    }

    int flag = CVodeSetSensDQMethod(CVODE_MEM_FROM_ML(vdata), dqtype,
				    Double_val(vdqrhomax));
    SCHECK_FLAG("CVodeSetSensMaxNonlinIters", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_set_max_nonlin_iters(value vdata, value vmaxcors)
{
    CAMLparam2(vdata, vmaxcors);

    int flag = CVodeSetSensMaxNonlinIters(CVODE_MEM_FROM_ML(vdata),
	    Int_val(vmaxcors));
    SCHECK_FLAG("CVodeSetSensMaxNonlinIters", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_sens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_rhs_evals_sens(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetNumRhsEvalsSens(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetNumRhsEvalsSens", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_lin_solv_setups(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumLinSolvSetups(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumLinSolvSetups", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfsevals;
    long int nfevalss;
    long int nsetfails;
    long int nlinsetupss;

    flag = CVodeGetSensStats(CVODE_MEM_FROM_ML(vdata), &nfsevals,
	    &nfevalss, &nsetfails, &nlinsetupss);
    SCHECK_FLAG("CVodeGetSensStats", flag);

    r = caml_alloc_tuple(RECORD_CVODES_SENS_STATS_SIZE);
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_RHS_EVALS, Val_long(nfsevals));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_SENS_EVALS, Val_long(nfevalss));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_ERR_TEST_FAILS,
							   Val_long(nsetfails));
    Store_field(r, RECORD_CVODES_SENS_STATS_NUM_LIN_SOLV_SETUPS,
							 Val_long(nlinsetupss));

    CAMLreturn(r);
}

CAMLprim value c_cvodes_sens_get_num_nonlin_solv_iters(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumNonlinSolvIters(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumNonlinSolvIters", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_num_nonlin_solv_conv_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetSensNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetSensNumNonlinSolvConvFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_sens_get_nonlin_solv_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nsniters;
    long int nsncfails;

    flag = CVodeGetSensNonlinSolvStats(CVODE_MEM_FROM_ML(vdata), &nsniters,
	    &nsncfails);
    SCHECK_FLAG("CVodeGetSensNonlinSolvStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nsniters));
    Store_field(r, 1, Val_long(nsncfails));

    CAMLreturn(r);
}

CAMLprim void c_cvodes_sens_get_num_stgr_nonlin_solv_iters(value vdata,
							  value vr)
{
    CAMLparam2(vdata, vr);

    int flag;

    flag = CVodeGetStgrSensNumNonlinSolvIters(CVODE_MEM_FROM_ML(vdata),
					      LONG_ARRAY(vr));
    SCHECK_FLAG("CVodeGetStgrSensNumNonlinSolvIters", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_sens_get_num_stgr_nonlin_solv_conv_fails(value vdata,
							        value vr)
{
    CAMLparam2(vdata, vr);

    int flag = CVodeGetStgrSensNumNonlinSolvConvFails(CVODE_MEM_FROM_ML(vdata),
						      LONG_ARRAY(vr));
    SCHECK_FLAG("CVodeGetStgrSensNumNonlinSolvConvFails", flag);

    CAMLreturn0;
}

/* sensitivity/quadrature interface */

CAMLprim void c_cvodes_quadsens_set_err_con(value vdata, value verrconq)
{
    CAMLparam2(vdata, verrconq);
    int flag;
    
    flag = CVodeSetQuadSensErrCon(CVODE_MEM_FROM_ML(vdata), Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadSensErrCon", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_ss_tolerances(value vdata,
					      value reltol,
					      value abstol)
{
    CAMLparam3(vdata, reltol, abstol);

    int flag = CVodeQuadSensSStolerances(CVODE_MEM_FROM_ML(vdata),
					 Double_val(reltol),
					 REAL_ARRAY(abstol));
    SCHECK_FLAG("CVodeQuadSensSStolerances", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_quadsens_ee_tolerances(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeQuadSensEEtolerances(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeQuadSensEEtolerances", flag);

    CAMLreturn0;
}

CAMLprim value c_cvodes_quadsens_get_num_rhs_evals(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadSensNumRhsEvals(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadSensNumRhsEvals", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quadsens_get_num_err_test_fails(value vdata)
{
    CAMLparam1(vdata);

    int flag;
    long int v;

    flag = CVodeGetQuadSensNumErrTestFails(CVODE_MEM_FROM_ML(vdata), &v);
    SCHECK_FLAG("CVodeGetQuadSensNumErrTestFails", flag);

    CAMLreturn(Val_long(v));
}

CAMLprim value c_cvodes_quadsens_get_stats(value vdata)
{
    CAMLparam1(vdata);
    CAMLlocal1(r);

    int flag;
    long int nfqsevals;
    long int nqsetfails;

    flag = CVodeGetQuadSensStats(CVODE_MEM_FROM_ML(vdata), &nfqsevals,
				 &nqsetfails);
    SCHECK_FLAG("CVodeGetQuadSensStats", flag);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_long(nfqsevals));
    Store_field(r, 1, Val_long(nqsetfails));

    CAMLreturn(r);
}

/* adjoint interface */

CAMLprim void c_cvodes_adj_init(value vdata, value vnd, value vinterptype)
{
    CAMLparam3(vdata, vnd, vinterptype);
    int interptype;

    switch(Int_val(vinterptype)) {
    case VARIANT_CVODES_ADJ_INTERPOLATION_POLYNOMIAL:
	interptype = CV_POLYNOMIAL;
	break;

    case VARIANT_CVODES_ADJ_INTERPOLATION_HERMITE:
	interptype = CV_HERMITE;
	break;

    default:
	caml_failwith("Illegal interpolation value.");
    }

    int flag = CVodeAdjInit(CVODE_MEM_FROM_ML(vdata), Long_val(vnd),
			    interptype);
    SCHECK_FLAG("CVodeAdjInit", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_ss_tolerances(value vparent, value vwhich,
					 value vreltol, value vabstol)
{
    CAMLparam4(vparent, vwhich, vreltol, vabstol);

    int flag = CVodeSStolerancesB(CVODE_MEM_FROM_ML(vparent),
		 Int_val(vwhich), Double_val(vreltol), Double_val(vabstol));
    SCHECK_FLAG("CVodeSStolerancesB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_diag(value vparent, value vwhich)
{
    CAMLparam2(vparent, vwhich);

    int flag = CVDiagB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich));
    SCHECK_FLAG("CVDiagB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_spgmr (value vparent, value vwhich,
				        value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpgmrB (cvode_mem, which, spils_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpgmrB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_spbcg (value vparent, value vwhich,
				        value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSpbcgB (cvode_mem, which, spils_precond_type (vtype),
		     Int_val (vmaxl));
    SCHECK_FLAG ("CVSpbcgB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_sptfqmr (value vparent, value vwhich,
				          value vmaxl, value vtype)
{
    CAMLparam4 (vparent, vwhich, vmaxl, vtype);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int which = Int_val(vwhich);
    int flag;

    flag = CVodeSetIterTypeB (cvode_mem, which, CV_NEWTON);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    flag = CVSptfqmrB (cvode_mem, which, spils_precond_type (vtype),
		       Int_val (vmaxl));
    SCHECK_FLAG ("CVSptfqmrB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_functional (value vparent, value vwhich)
{
    CAMLparam2 (vparent, vwhich);
    int flag = CVodeSetIterTypeB (CVODE_MEM_FROM_ML (vparent), Int_val(vwhich),
				  CV_FUNCTIONAL);
    SCHECK_FLAG ("CVodeSetIterTypeB", flag);
    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_bsession_finalize(value vdata)
{
    if (CVODE_MEM_FROM_ML(vdata) != NULL) {
	value *backref = CVODE_BACKREF_FROM_ML(vdata);
	// NB: CVodeFree() is *not* called: parents free-up backward problems
	caml_remove_generational_global_root (backref);
	free (backref);
    }
}

CAMLprim void c_cvodes_adj_backward_normal(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = CVodeB(CVODE_MEM_FROM_ML(vdata), Double_val(vtbout), CV_NORMAL);
    SCHECK_FLAG("CVodeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_backward_one_step(value vdata, value vtbout)
{
    CAMLparam2(vdata, vtbout);

    int flag = CVodeB(CVODE_MEM_FROM_ML(vdata), Double_val(vtbout),
		      CV_ONE_STEP);
    SCHECK_FLAG("CVodeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_no_sensitivity(value vdata)
{
    CAMLparam1(vdata);

    int flag = CVodeSetAdjNoSensi(CVODE_MEM_FROM_ML(vdata));
    SCHECK_FLAG("CVodeSetAdjNoSensi", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_max_ord(value vparent, value vwhich,
				       value vmaxord)
{
    CAMLparam3(vparent, vwhich, vmaxord);

    int flag = CVodeSetMaxOrdB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(vmaxord));
    SCHECK_FLAG("CVodeSetMaxOrdB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_max_num_steps(value vparent, value vwhich,
					     value vmxsteps)
{
    CAMLparam3(vparent, vwhich, vmxsteps);

    int flag = CVodeSetMaxNumStepsB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				    Long_val(vmxsteps));
    SCHECK_FLAG("CVodeSetMaxNumStepsB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_init_step(value vparent, value vwhich, value vhin)
{
    CAMLparam3(vparent, vwhich, vhin);

    int flag = CVodeSetInitStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			         Double_val(vhin));
    SCHECK_FLAG("CVodeSetInitStepB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_min_step(value vparent, value vwhich, value vhmin)
{
    CAMLparam3(vparent, vwhich, vhmin);

    int flag = CVodeSetMinStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			        Double_val(vhmin));
    SCHECK_FLAG("CVodeSetMinStepB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_max_step(value vparent, value vwhich, value vhmax)
{
    CAMLparam3(vparent, vwhich, vhmax);

    int flag = CVodeSetMaxStepB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			        Double_val(vhmax));
    SCHECK_FLAG("CVodeSetMaxStepB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_set_stab_lim_det(value vparent, value vwhich,
					    value vstldet)
{
    CAMLparam3(vparent, vwhich, vstldet);

    int flag = CVodeSetStabLimDetB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			           Bool_val(vstldet));
    SCHECK_FLAG("CVodeSetStabLimDetB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_set_prec_type(value vparent, value vwhich,
					       value vptype)
{
    CAMLparam3(vparent, vwhich, vptype);

    int flag = CVSpilsSetPrecTypeB(CVODE_MEM_FROM_ML(vparent),
				   Int_val(vwhich), spils_precond_type(vptype));
    SCHECK_FLAG("CVSpilsSetPrecTypeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_set_gs_type(value vparent, value vwhich,
				             value vgstype)
{
    CAMLparam3(vparent, vwhich, vgstype);

    int flag = CVSpilsSetGSTypeB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				 spils_gs_type(vgstype));
    SCHECK_FLAG("CVSpilsSetGSTypeB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_set_eps_lin(value vparent, value vwhich,
				             value eplifac)
{
    CAMLparam3(vparent, vwhich, eplifac);

    int flag = CVSpilsSetEpsLinB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
				 Double_val(eplifac));
    SCHECK_FLAG("CVSpilsSetEpsLinB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adj_spils_set_maxl(value vparent, value vwhich,
					  value maxl)
{
    CAMLparam3(vparent, vwhich, maxl);

    int flag = CVSpilsSetMaxlB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Int_val(maxl));
    SCHECK_FLAG("CVSpilsSetMaxlB", flag);

    CAMLreturn0;
}

/* adjoint/quadrature interface */

CAMLprim void c_cvodes_adjquad_set_err_con(value vparent, value vwhich,
					   value verrconq)
{
    CAMLparam3(vparent, vwhich, verrconq);
    int flag;
    
    flag = CVodeSetQuadErrConB(CVODE_MEM_FROM_ML(vparent), Int_val(vwhich),
			       Bool_val(verrconq));
    SCHECK_FLAG("CVodeSetQuadErrConB", flag);

    CAMLreturn0;
}

CAMLprim void c_cvodes_adjquad_ss_tolerances(value vparent, value vwhich,
					     value reltol, value abstol)
{
    CAMLparam4(vparent, vwhich, reltol, abstol);

    int flag = CVodeQuadSStolerancesB(CVODE_MEM_FROM_ML(vparent),
				      Int_val(vwhich), Double_val(reltol),
				      Double_val(abstol));
    SCHECK_FLAG("CVodeQuadSStolerancesB", flag);

    CAMLreturn0;
}

