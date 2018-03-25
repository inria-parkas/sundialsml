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

#include "../sundials/sundials_ml.h"

#include <sundials/sundials_config.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_spgmr.h>
#if SUNDIALS_LIB_VERSION >= 260
#include <sundials/sundials_spfgmr.h>
#include <sundials/sundials_pcg.h>
#endif
#include <sundials/sundials_spbcgs.h>
#include <sundials/sundials_sptfqmr.h>
#include <sundials/sundials_nvector.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include "../lsolvers/spils_ml.h"
#include "../nvectors/nvector_ml.h"

// we only use the size in caml_alloc_final
#define NVECTOR_APPROXSIZE(v) (1)

#define CVTYPE(fname) c_nvec_spils_ ## fname
#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define CVTYPESTR(fname) QUOTE(CVTYPE(fname))

CAMLprim value c_spils_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (SPILS, exns);
    CAMLreturn (Val_unit);
}

CAMLprim value c_spils_modified_gs(value vv, value vh, value vk, value vp)
{
    CAMLparam4(vv, vh, vk, vp);

    int p = Int_val(vp);
    int k = Int_val(vk);
    int i;
    realtype new_vk_norm;
    N_Vector* v;

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < k) || (hn < k))
	caml_invalid_argument("Spils.modified_gs: h is too small.");
    if (Wosize_val (vv) < k)
	caml_invalid_argument("Spils.modified_gs: v is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_out_of_memory();
    for (i = 0; i <= p; ++i) {
	v[i] = NVEC_VAL(Field(vv, k - p + i));
    }

    ModifiedGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm);

    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value c_spils_classical_gs(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal3(vv, vh, vs);

    int k = Int_val(Field(vargs, 2));
    int p = Int_val(Field(vargs, 3));
    N_Vector temp;

    vv = Field(vargs, 0);
    vh = Field(vargs, 1);
    vs = Field(vargs, 5);

    int i;
    realtype new_vk_norm;
    N_Vector* v;

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < k) || (hn < k))
	caml_invalid_argument("Spils.classical_gs: h is too small.");
    if (Wosize_val (vv) < k)
	caml_invalid_argument("Spils.classical_gs: v is too small.");
    if (ARRAY1_LEN(vs) < k)
	caml_invalid_argument("Spils.classical_gs: s is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_out_of_memory();
    for (i = 0; i <= p; ++i) {
	v[i] = NVEC_VAL(Field(vv, k - p + i));
    }

    temp = NVEC_VAL(Field(vargs, 4));
    ClassicalGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm,
	        temp, REAL_ARRAY(vs));

    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

static int atimes_f(void *a_fn, N_Vector v, N_Vector z)
{
    CAMLparam0();
    CAMLlocal2(vv, vz);
    int r = 0;

    vv = NVEC_BACKLINK(v);
    vz = NVEC_BACKLINK(z);

    // the data payloads inside vv and vz are only valid during this call,
    // afterward that memory goes back to Sundials. These bigarrays must not
    // be retained by a_fn! If it wants a permanent copy, then it has to
    // make it manually.

    /* NB: Don't trigger GC while processing this return value!  */
    value res = caml_callback2_exn((value)a_fn, vv, vz);
    if (Is_exception_result(res)) {
	res = Extract_exception(res);
	r = (Field(res, 0) == SUNDIALS_EXN_TAG(RecoverableFailure))?1:-1;
    }

    CAMLreturnT(int, r);
}

#if SUNDIALS_LIB_VERSION >= 300
static int psolve_f(void *p_fn, N_Vector r, N_Vector z, realtype tol, int lr)
#else
static int psolve_f(void *p_fn, N_Vector r, N_Vector z, int lr)
#endif
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int fr = 0;

    args[0] = NVEC_BACKLINK(r);	     /* vr */
    args[1] = NVEC_BACKLINK(z);	     /* vz */
#if SUNDIALS_LIB_VERSION >= 300
    args[2] = caml_copy_double(tol); /* tol */
#else
    args[2] = caml_copy_double(0.0); /* tol */
#endif
    args[3] = Val_bool(lr == 1);     /* lr */

    // the data payloads inside vr and vz are only valid during this call,
    // afterward that memory goes back to Sundials. These bigarrays must not
    // be retained by a_fn! If it wants a permanent copy, then it has to
    // make it manually.

    /* NB: Don't trigger GC while processing this return value!  */
    value res = caml_callbackN_exn((value)p_fn, 4, args);
    if (Is_exception_result(res)) {
	res = Extract_exception(res);
	fr = (Field(res, 0) == SUNDIALS_EXN_TAG(RecoverableFailure))?1:-1;
    }

    CAMLreturnT(int, fr);
}

/* SPGMR */

#define SPGMR_SESSION(v) (*((SpgmrMem*)Data_custom_val(v)))

static void spgmr_finalize(value vs)
{
    SpgmrMem s = SPGMR_SESSION(vs);
    SpgmrFree(s);
}

CAMLprim value c_spils_spgmr_make(value vl_max, value vvec_tmpl)
{
    CAMLparam2(vl_max, vvec_tmpl);
    CAMLlocal1(vs);
    SpgmrMem s;
    N_Vector vec_tmpl;
    mlsize_t approx_size = sizeof(SpgmrMemRec)
	+ (Int_val(vl_max) + 3) * NVECTOR_APPROXSIZE(vvec_tmpl);

    vec_tmpl = NVEC_VAL(vvec_tmpl);
    s = SpgmrMalloc(Int_val(vl_max), vec_tmpl);

    if (s == NULL) caml_raise_out_of_memory();
    vs = caml_alloc_final(1, &spgmr_finalize, approx_size, approx_size * 20);
    SPGMR_SESSION(vs) = s;

    CAMLreturn(vs);
}

CAMLprim value c_spils_spgmr_solve(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal4(vr, vpsolve, vs1, vs2);
    int flag;
    int success;
    N_Vector x = NVEC_VAL(Field(vargs, 1));
    N_Vector b = NVEC_VAL(Field(vargs, 2));
    N_Vector s1 = NULL;
    N_Vector s2 = NULL;
    realtype res_norm = 0.0;
    int nli = 0;
    int nps = 0;
    void *p_data = NULL;

    vpsolve = Field(vargs, 10);
    if (vpsolve != Val_none) p_data = (void *)Field(vpsolve, 0);

    vs1 = Field(vargs, 7);
    if (vs1 != Val_none) s1 = NVEC_VAL(Field(vs1, 0));

    vs2 = Field(vargs, 8);
    if (vs2 != Val_none) s2 = NVEC_VAL(Field(vs2, 0));

    flag = SpgmrSolve(SPGMR_SESSION(Field(vargs, 0)),
		      (void *)Field(vargs, 9), // adata = user atimes
		      x,
		      b,
		      spils_precond_type(Field(vargs, 3)), // pretype
		      spils_gs_type(Field(vargs, 4)),	   // gstype
		      Double_val(Field(vargs, 5)),         // delta
		      Int_val(Field(vargs, 6)),            // max_restarts
		      p_data,				   // p_data = user psolve
		      s1,
		      s2,
		      atimes_f,				   // atimes wrapper
		      psolve_f,                            // psolve wrapper
		      &res_norm,
		      &nli,
		      &nps);

    switch(flag) {
      case SPGMR_SUCCESS:
	  success = 1;
          break;

      case SPGMR_RES_REDUCED:
	  success = 0;
	  break;

      case SPGMR_CONV_FAIL:
	caml_raise_constant(SPILS_EXN(ConvFailure));

      case SPGMR_QRFACT_FAIL:
	caml_raise_constant(SPILS_EXN(QRfactFailure));

      case SPGMR_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(1));

      case SPGMR_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(0));

      case SPGMR_ATIMES_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(1));

      case SPGMR_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(0));

      case SPGMR_PSET_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(1));

      case SPGMR_PSET_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(0));

      case SPGMR_GS_FAIL:
	caml_raise_constant(SPILS_EXN(GSFailure));

      case SPGMR_QRSOL_FAIL:
	caml_raise_constant(SPILS_EXN(QRSolFailure));

      default:
	caml_failwith("spmgr_solve: unexpected failure");
    }

    vr = caml_alloc_tuple(4);
    Store_field(vr, 0, Val_bool(success));
    Store_field(vr, 1, caml_copy_double(res_norm));
    Store_field(vr, 2, Val_int(nli));
    Store_field(vr, 3, Val_int(nps));

    CAMLreturn(vr);
}

/* SPFGMR */

#if SUNDIALS_LIB_VERSION >= 260

#define SPFGMR_SESSION(v) (*((SpfgmrMem*)Data_custom_val(v)))

static void spfgmr_finalize(value vs)
{
    SpfgmrMem s = SPFGMR_SESSION(vs);
    SpfgmrFree(s);
}

CAMLprim value c_spils_spfgmr_make(value vl_max, value vvec_tmpl)
{
    CAMLparam2(vl_max, vvec_tmpl);
    CAMLlocal1(vs);
    SpfgmrMem s;
    N_Vector vec_tmpl;
    mlsize_t approx_size = sizeof(SpfgmrMemRec)
	+ (Int_val(vl_max) + 3) * NVECTOR_APPROXSIZE(vvec_tmpl);

    vec_tmpl = NVEC_VAL(vvec_tmpl);
    s = SpfgmrMalloc(Int_val(vl_max), vec_tmpl);

    if (s == NULL) caml_raise_out_of_memory();
    vs = caml_alloc_final(1, &spfgmr_finalize, approx_size, approx_size * 20);
    SPFGMR_SESSION(vs) = s;

    CAMLreturn(vs);
}

CAMLprim value c_spils_spfgmr_solve(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal4(vr, vpsolve, vs1, vs2);
    int flag;
    int success;
    N_Vector x = NVEC_VAL(Field(vargs, 1));
    N_Vector b = NVEC_VAL(Field(vargs, 2));
    N_Vector s1 = NULL;
    N_Vector s2 = NULL;
    realtype res_norm = 0.0;
    int nli = 0;
    int nps = 0;
    void *p_data = NULL;

    vpsolve = Field(vargs, 10);
    if (vpsolve != Val_none) p_data = (void *)Field(vpsolve, 0);

    vs1 = Field(vargs, 7);
    if (vs1 != Val_none) s1 = NVEC_VAL(Field(vs1, 0));

    vs2 = Field(vargs, 8);
    if (vs2 != Val_none) s2 = NVEC_VAL(Field(vs2, 0));

    flag = SpfgmrSolve(SPFGMR_SESSION(Field(vargs, 0)),
		       (void *)Field(vargs, 9), // adata = user atimes
		       x,
		       b,
		       spils_precond_type(Field(vargs, 3)), // pretype
		       spils_gs_type(Field(vargs, 4)),	    // gstype
		       Double_val(Field(vargs, 5)),         // delta
		       Int_val(Field(vargs, 6)),            // max_restarts
		       Int_val(Field(vargs, 11)),	    // maxit
		       p_data,				    // p_data = user psolve
		       s1,
		       s2,
		       atimes_f,			    // atimes wrapper
		       psolve_f,                            // psolve wrapper
		       &res_norm,
		       &nli,
		       &nps);

    switch(flag) {
      case SPFGMR_SUCCESS:
	  success = 1;
          break;

      case SPFGMR_RES_REDUCED:
	  success = 0;
	  break;

      case SPFGMR_CONV_FAIL:
	caml_raise_constant(SPILS_EXN(ConvFailure));

      case SPFGMR_QRFACT_FAIL:
	caml_raise_constant(SPILS_EXN(QRfactFailure));

      case SPFGMR_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(1));

      case SPFGMR_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(0));

      case SPFGMR_ATIMES_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(1));

      case SPFGMR_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(0));

      case SPFGMR_PSET_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(1));

      case SPFGMR_PSET_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(0));

      case SPFGMR_GS_FAIL:
	caml_raise_constant(SPILS_EXN(GSFailure));

      case SPFGMR_QRSOL_FAIL:
	caml_raise_constant(SPILS_EXN(QRSolFailure));

      default:
	caml_failwith("spfmgr_solve: unexpected failure");
    }

    vr = caml_alloc_tuple(4);
    Store_field(vr, 0, Val_bool(success));
    Store_field(vr, 1, caml_copy_double(res_norm));
    Store_field(vr, 2, Val_int(nli));
    Store_field(vr, 3, Val_int(nps));

    CAMLreturn(vr);
}

#else
CAMLprim value c_spils_spfgmr_make(value vl_max, value vvec_tmpl)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}

CAMLprim value c_spils_spfgmr_solve(value vargs)
{
    // Never called anyway since it is not possible to make an SPFGMR.t.
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

/* SPBCG */

#define SPBCG_SESSION(v) (*((SpbcgMem*)Data_custom_val(v)))

static void spbcg_finalize(value vs)
{
    SpbcgMem s = SPBCG_SESSION(vs);
    SpbcgFree(s);
}

CAMLprim value c_spils_spbcg_make(value vl_max, value vvec_tmpl)
{
    CAMLparam2(vl_max, vvec_tmpl);
    CAMLlocal1(vs);
    SpbcgMem s;
    N_Vector vec_tmpl;
    mlsize_t approx_size = sizeof(SpbcgMemRec)
	+ (Int_val(vl_max) + 3) * NVECTOR_APPROXSIZE(vvec_tmpl);

    vec_tmpl = NVEC_VAL(vvec_tmpl);
    s = SpbcgMalloc(Int_val(vl_max), vec_tmpl);

    if (s == NULL) caml_raise_out_of_memory();
    vs = caml_alloc_final(1, &spbcg_finalize, approx_size, approx_size * 20);
    SPBCG_SESSION(vs) = s;

    CAMLreturn(vs);
}

CAMLprim value c_spils_spbcg_solve(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal4(vr, vpsolve, vsx, vsb);
    int flag;
    int success;
    N_Vector x = NVEC_VAL(Field(vargs, 1));
    N_Vector b = NVEC_VAL(Field(vargs, 2));
    N_Vector sx = NULL;
    N_Vector sb = NULL;
    realtype res_norm = 0.0;
    int nli = 0;
    int nps = 0;
    void *p_data = NULL;

    vpsolve = Field(vargs, 8);
    if (vpsolve != Val_none) p_data = (void *)Field(vpsolve, 0);

    vsx = Field(vargs, 5);
    if (vsx != Val_none) sx = NVEC_VAL(Field(vsx, 0));

    vsb = Field(vargs, 6);
    if (vsb != Val_none) sb = NVEC_VAL(Field(vsb, 0));

    flag = SpbcgSolve(SPBCG_SESSION(Field(vargs, 0)),
		      (void *)Field(vargs, 7),		   // adata = user atimes
		      x,
		      b,
		      spils_precond_type(Field(vargs, 3)), // pretype
		      Double_val(Field(vargs, 4)),         // delta
		      p_data,				   // p_data = user psolve
		      sx,
		      sb,
		      atimes_f,				   // atimes wrapper
		      psolve_f,                            // psolve wrapper
		      &res_norm,
		      &nli,
		      &nps);

    switch(flag) {
      case SPBCG_SUCCESS:
	  success = 1;
          break;

      case SPBCG_RES_REDUCED:
	  success = 0;
	  break;

      case SPBCG_CONV_FAIL:
	caml_raise_constant(SPILS_EXN(ConvFailure));

      case SPBCG_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(1));

      case SPBCG_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(0));

      case SPBCG_ATIMES_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(1));

      case SPBCG_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(0));

      case SPBCG_PSET_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(1));

      case SPBCG_PSET_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(0));

      default:
	caml_failwith("spmgr_solve: unexpected failure");
    }

    vr = caml_alloc_tuple(4);
    Store_field(vr, 0, Val_bool(success));
    Store_field(vr, 1, caml_copy_double(res_norm));
    Store_field(vr, 2, Val_int(nli));
    Store_field(vr, 3, Val_int(nps));

    CAMLreturn(vr);
}

/* SPTFQMR */

#define SPTFQMR_SESSION(v) (*((SptfqmrMem*)Data_custom_val(v)))

static void sptfqmr_finalize(value vs)
{
    SptfqmrMem s = SPTFQMR_SESSION(vs);
    SptfqmrFree(s);
}

CAMLprim value c_spils_sptfqmr_make(value vl_max, value vvec_tmpl)
{
    CAMLparam2(vl_max, vvec_tmpl);
    CAMLlocal1(vs);
    SptfqmrMem s;
    N_Vector vec_tmpl;
    mlsize_t approx_size = sizeof(SptfqmrMemRec)
	+ (Int_val(vl_max) + 3) * NVECTOR_APPROXSIZE(vvec_tmpl);

    vec_tmpl = NVEC_VAL(vvec_tmpl);
    s = SptfqmrMalloc(Int_val(vl_max), vec_tmpl);

    if (s == NULL) caml_raise_out_of_memory();
    vs = caml_alloc_final(1, &sptfqmr_finalize, approx_size, approx_size * 20);
    SPTFQMR_SESSION(vs) = s;

    CAMLreturn(vs);
}

CAMLprim value c_spils_sptfqmr_solve(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal4(vr, vpsolve, vsx, vsb);
    int flag;
    int success;
    N_Vector x = NVEC_VAL(Field(vargs, 1));
    N_Vector b = NVEC_VAL(Field(vargs, 2));
    N_Vector sx = NULL;
    N_Vector sb = NULL;
    realtype res_norm = 0.0;
    int nli = 0;
    int nps = 0;
    void *p_data = NULL;

    vpsolve = Field(vargs, 8);
    if (vpsolve != Val_none) p_data = (void *)Field(vpsolve, 0);

    vsx = Field(vargs, 5);
    if (vsx != Val_none) sx = NVEC_VAL(Field(vsx, 0));

    vsb = Field(vargs, 6);
    if (vsb != Val_none) sb = NVEC_VAL(Field(vsb, 0));

    flag = SptfqmrSolve(SPTFQMR_SESSION(Field(vargs, 0)),
		      (void *)Field(vargs, 7),		   // adata = user atimes
		      x,
		      b,
		      spils_precond_type(Field(vargs, 3)), // pretype
		      Double_val(Field(vargs, 4)),         // delta
		      p_data,				   // p_data = user psolve
		      sx,
		      sb,
		      atimes_f,				   // atimes wrapper
		      psolve_f,                            // psolve wrapper
		      &res_norm,
		      &nli,
		      &nps);

    switch(flag) {
      case SPTFQMR_SUCCESS:
	  success = 1;
          break;

      case SPTFQMR_RES_REDUCED:
	  success = 0;
	  break;

      case SPTFQMR_CONV_FAIL:
	  caml_raise_constant(SPILS_EXN(ConvFailure));

      case SPTFQMR_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(1));

      case SPTFQMR_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(0));

      case SPTFQMR_ATIMES_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(1));

      case SPTFQMR_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(0));

      case SPTFQMR_PSET_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(1));

      case SPTFQMR_PSET_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(0));

      default:
	caml_failwith("spmgr_solve: unexpected failure");
    }

    vr = caml_alloc_tuple(4);
    Store_field(vr, 0, Val_bool(success));
    Store_field(vr, 1, caml_copy_double(res_norm));
    Store_field(vr, 2, Val_int(nli));
    Store_field(vr, 3, Val_int(nps));

    CAMLreturn(vr);
}

/* PCG */

#if SUNDIALS_LIB_VERSION >= 260

#define PCG_SESSION(v) (*((PcgMem*)Data_custom_val(v)))

static void pcg_finalize(value vs)
{
    PcgMem s = PCG_SESSION(vs);
    PcgFree(s);
}

CAMLprim value c_spils_pcg_make(value vl_max, value vvec_tmpl)
{
    CAMLparam2(vl_max, vvec_tmpl);
    CAMLlocal1(vs);
    PcgMem s;
    N_Vector vec_tmpl;
    mlsize_t approx_size = sizeof(PcgMemRec)
	+ (Int_val(vl_max) + 3) * NVECTOR_APPROXSIZE(vvec_tmpl);

    vec_tmpl = NVEC_VAL(vvec_tmpl);
    s = PcgMalloc(Int_val(vl_max), vec_tmpl);

    if (s == NULL) caml_raise_out_of_memory();
    vs = caml_alloc_final(1, &pcg_finalize, approx_size, approx_size * 20);
    PCG_SESSION(vs) = s;

    CAMLreturn(vs);
}

CAMLprim value c_spils_pcg_solve(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal2(vr, vpsolve);
    int flag;
    int success;
    N_Vector x = NVEC_VAL(Field(vargs, 1));
    N_Vector b = NVEC_VAL(Field(vargs, 2));
    N_Vector w = NVEC_VAL(Field(vargs, 3));
    realtype res_norm = 0.0;
    int nli = 0;
    int nps = 0;
    void *p_data = NULL;

    vpsolve = Field(vargs, 7);
    if (vpsolve != Val_none) p_data = (void *)Field(vpsolve, 0);

    flag = PcgSolve(PCG_SESSION(Field(vargs, 0)),
		    (void *)Field(vargs, 6), // adata = user atimes
		    x,
		    b,
		    spils_precond_type(Field(vargs, 4)), // pretype
		    Double_val(Field(vargs, 5)),         // delta
		    p_data,				 // p_data = user psolve
		    w,
		    atimes_f,			         // atimes wrapper
		    psolve_f,                            // psolve wrapper
		    &res_norm,
		    &nli,
		    &nps);

    switch(flag) {
      case PCG_SUCCESS:
	  success = 1;
          break;

      case PCG_RES_REDUCED:
	  success = 0;
	  break;

      case PCG_CONV_FAIL:
	caml_raise_constant(SPILS_EXN(ConvFailure));

      case PCG_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(1));

      case PCG_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSolveFailure), Val_bool(0));

      case PCG_ATIMES_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(1));

      case PCG_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(ATimesFailure), Val_bool(0));

      case PCG_PSET_FAIL_REC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(1));

      case PCG_PSET_FAIL_UNREC:
	  caml_raise_with_arg(SPILS_EXN_TAG(PSetFailure), Val_bool(0));

      default:
	caml_failwith("pcg_solve: unexpected failure");
    }

    vr = caml_alloc_tuple(4);
    Store_field(vr, 0, Val_bool(success));
    Store_field(vr, 1, caml_copy_double(res_norm));
    Store_field(vr, 2, Val_int(nli));
    Store_field(vr, 3, Val_int(nps));

    CAMLreturn(vr);
}

#else
CAMLprim value c_spils_pcg_make(value vl_max, value vvec_tmpl)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}

CAMLprim value c_spils_pcg_solve(value vargs)
{
    // Never called anyway since it is not possible to make an PCG.t.
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

int spils_precond_type(value vptype)
{
    CAMLparam1(vptype);

    int ptype = -1;
    switch (Int_val(vptype)) {
    case VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_NONE:
	ptype = PREC_NONE;
	break;

    case VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_LEFT:
	ptype = PREC_LEFT;
	break;

    case VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_RIGHT:
	ptype = PREC_RIGHT;
	break;

    case VARIANT_SPILS_PRECONDITIONING_TYPE_PREC_TYPE_BOTH:
	ptype = PREC_BOTH;
	break;
    }

    CAMLreturnT(int, ptype);
}

int spils_gs_type(value vgstype)
{
    CAMLparam1(vgstype);

    int gstype = -1;
    switch (Int_val(vgstype)) {
    case VARIANT_SPILS_GRAMSCHMIDT_TYPE_MODIFIEDGS:
	gstype = MODIFIED_GS;
	break;

    case VARIANT_SPILS_GRAMSCHMIDT_TYPE_CLASSICALGS:
	gstype = CLASSICAL_GS;
	break;
    }

    CAMLreturnT(int, gstype);
}

CAMLprim value c_spils_qr_fact(value vh, value vq, value vnewjob)
{
    CAMLparam3(vh, vq, vnewjob);
    int r;
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hn = bh->dim[0];

#if SUNDIALS_ML_SAFE == 1
    intnat hm = bh->dim[1];

    if ((hm < hn + 1))
	caml_invalid_argument("Spils.qr_fact: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * hn)
	caml_invalid_argument("Spils.qr_fact: q is too small.");
#endif

    r = QRfact(hn, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), Bool_val(vnewjob));

    if (r != 0) {
	caml_raise_with_arg(SPILS_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

CAMLprim value c_spils_qr_sol(value vh, value vq, value vb)
{
    CAMLparam3(vh, vq, vb);
    int r;
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hn = bh->dim[0];

#if SUNDIALS_ML_SAFE == 1
    intnat hm = bh->dim[1];

    if (hm < hn + 1)
	caml_invalid_argument("Spils.qr_sol: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * hn)
	caml_invalid_argument("Spils.qr_sol: q is too small.");
    if (ARRAY1_LEN(vb) < hn + 1)
	caml_invalid_argument("Spils.qr_sol: b is too small.");
#endif

    r = QRsol(hn, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), REAL_ARRAY(vb));

    if (r != 0) {
	caml_raise_with_arg(SPILS_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

