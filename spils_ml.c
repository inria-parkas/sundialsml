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

#include <sundials/sundials_config.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_spgmr.h>
#include <sundials/sundials_spbcgs.h>
#include <sundials/sundials_sptfqmr.h>
#include <sundials/sundials_nvector.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

#include "sundials_ml.h"
#include "spils_ml.h"
#include "nvector_ml.h"

// we only use the size in caml_alloc_final
#define NVECTOR_APPROXSIZE(v) (1)

#define CVTYPE(fname) c_nvec_spils_ ## fname
#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define CVTYPESTR(fname) QUOTE(CVTYPE(fname))

CAMLprim value c_spils_modified_gs(value vv, value vh, value vk, value vp)
{
    CAMLparam4(vv, vh, vk, vp);

    int p = Int_val(vp);
    int k = Int_val(vk);
    int i;
    realtype new_vk_norm;
    N_Vector* v;

#if CHECK_MATRIX_ACCESS == 1
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
	caml_raise_constant(*caml_named_value("spils_MemoryRequestFailure"));
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

#if CHECK_MATRIX_ACCESS == 1
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
	caml_raise_constant(*caml_named_value("spils_MemoryRequestFailure"));
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
    CAMLlocal3(vv, vz, res);
    int r = 0;

    vv = NVEC_BACKLINK(v);
    vz = NVEC_BACKLINK(z);

    // the data payloads inside vv and vz are only valid during this call,
    // afterward that memory goes back to Sundials. These bigarrays must not
    // be retained by a_fn! If it wants a permanent copy, then it has to
    // make it manually.
    res = caml_callback2_exn((value)a_fn, vv, vz);
    if (Is_exception_result(res)) {
	res = Extract_exception(res);
	r = (Field(res, 0) == *caml_named_value(CVTYPESTR(RecoverableFailure)))?1:-1;
    }

    CAMLreturnT(int, r);
}

static int psolve_f(void *p_fn, N_Vector r, N_Vector z, int lr)
{
    CAMLparam0();
    CAMLlocal3(vr, vz, res);
    int fr = 0;

    vr = NVEC_BACKLINK(r);
    vz = NVEC_BACKLINK(z);

    // the data payloads inside vr and vz are only valid during this call,
    // afterward that memory goes back to Sundials. These bigarrays must not
    // be retained by a_fn! If it wants a permanent copy, then it has to
    // make it manually.
    res = caml_callback3_exn((value)p_fn, vr, vz, Val_bool(lr == 1));
    if (Is_exception_result(res)) {
	res = Extract_exception(res);
	fr = (Field(res, 0) == *caml_named_value(CVTYPESTR(RecoverableFailure)))?1:-1;
    }

    CAMLreturnT(int, fr);
}

/* SPGMR */

#define SPGMR_SESSION(v) (*((SpgmrMem*)Data_custom_val(v)))

static CAMLprim void spgmr_finalize(value vs)
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
    vs = caml_alloc_final(2, &spgmr_finalize, approx_size, approx_size * 20);
    Store_field(vs, 1, (value)s);

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
	caml_raise_constant(*caml_named_value("spils_ConvFailure"));

      case SPGMR_QRFACT_FAIL:
	caml_raise_constant(*caml_named_value("spils_QRfactFailure"));

      case SPGMR_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_PSolveFailure"), Val_bool(1));

      case SPGMR_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_PSolveFailure"), Val_bool(0));

      case SPGMR_ATIMES_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_ATimesFailure"), Val_bool(1));

      case SPGMR_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_ATimesFailure"), Val_bool(0));

      case SPGMR_PSET_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_PSetFailure"), Val_bool(1));

      case SPGMR_PSET_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_PSetFailure"), Val_bool(0));

      case SPGMR_GS_FAIL:
	caml_raise_constant(*caml_named_value("spils_GSFailure"));

      case SPGMR_QRSOL_FAIL:
	caml_raise_constant(*caml_named_value("spils_QRSolFailure"));

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

/* SPBCG */

#define SPBCG_SESSION(v) (*((SpbcgMem*)Data_custom_val(v)))

static CAMLprim void spbcg_finalize(value vs)
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
    vs = caml_alloc_final(2, &spbcg_finalize, approx_size, approx_size * 20);
    Store_field(vs, 1, (value)s);

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
	caml_raise_constant(*caml_named_value("spils_ConvFailure"));

      case SPBCG_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_PSolveFailure"), Val_bool(1));

      case SPBCG_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_PSolveFailure"), Val_bool(0));

      case SPBCG_ATIMES_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_ATimesFailure"), Val_bool(1));

      case SPBCG_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_ATimesFailure"), Val_bool(0));

      case SPBCG_PSET_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_PSetFailure"), Val_bool(1));

      case SPBCG_PSET_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_PSetFailure"), Val_bool(0));

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

static CAMLprim void sptfqmr_finalize(value vs)
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
    vs = caml_alloc_final(2, &sptfqmr_finalize, approx_size, approx_size * 20);
    Store_field(vs, 1, (value)s);

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
	caml_raise_constant(*caml_named_value("spils_ConvFailure"));

      case SPTFQMR_PSOLVE_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_PSolveFailure"), Val_bool(1));

      case SPTFQMR_PSOLVE_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_PSolveFailure"), Val_bool(0));

      case SPTFQMR_ATIMES_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_ATimesFailure"), Val_bool(1));

      case SPTFQMR_ATIMES_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_ATimesFailure"), Val_bool(0));

      case SPTFQMR_PSET_FAIL_REC:
	  caml_raise_with_arg(*caml_named_value("spils_PSetFailure"), Val_bool(1));

      case SPTFQMR_PSET_FAIL_UNREC:
	  caml_raise_with_arg(*caml_named_value("spils_PSetFailure"), Val_bool(0));

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


int spils_precond_type(value vptype)
{
    CAMLparam1(vptype);

    int ptype;
    switch (Int_val(vptype)) {
    case VARIANT_SPILS_PRECONDITIONING_TYPE_PRECNONE:
	ptype = PREC_NONE;
	break;

    case VARIANT_SPILS_PRECONDITIONING_TYPE_PRECLEFT:
	ptype = PREC_LEFT;
	break;

    case VARIANT_SPILS_PRECONDITIONING_TYPE_PRECRIGHT:
	ptype = PREC_RIGHT;
	break;

    case VARIANT_SPILS_PRECONDITIONING_TYPE_PRECBOTH:
	ptype = PREC_BOTH;
	break;
    }

    CAMLreturn(ptype);
}

int spils_gs_type(value vgstype)
{
    CAMLparam1(vgstype);

    int gstype;
    switch (Int_val(vgstype)) {
    case VARIANT_SPILS_GRAMSCHMIDT_TYPE_MODIFIEDGS:
	gstype = MODIFIED_GS;
	break;

    case VARIANT_SPILS_GRAMSCHMIDT_TYPE_CLASSICALGS:
	gstype = CLASSICAL_GS;
	break;
    }

    CAMLreturn(gstype);
}

CAMLprim value c_spils_qr_fact(value vn, value vh, value vq, value vnewjob)
{
    CAMLparam4(vn, vh, vq, vnewjob);
    int r;
    int n = Int_val(vn);

#if CHECK_MATRIX_ACCESS == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < n + 1) || (hn < n))
	caml_invalid_argument("Spils.qr_fact: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * n)
	caml_invalid_argument("Spils.qr_fact: q is too small.");
#endif

    r = QRfact(n, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), Bool_val(vnewjob));

    CAMLreturn(Val_int(r));
}

CAMLprim value c_spils_qr_sol(value vn, value vh, value vq, value vb)
{
    CAMLparam4(vn, vh, vq, vb);
    int r;
    int n = Int_val(vn);

#if CHECK_MATRIX_ACCESS == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < n + 1) || (hn < n))
	caml_invalid_argument("Spils.qr_sol: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * n)
	caml_invalid_argument("Spils.qr_sol: q is too small.");
    if (ARRAY1_LEN(vb) < n + 1)
	caml_invalid_argument("Spils.qr_sol: b is too small.");
#endif

    r = QRsol(n, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), REAL_ARRAY(vb));

    CAMLreturn(Val_int(r));
}

