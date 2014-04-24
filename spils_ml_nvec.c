/***********************************************************************
 *                                                                     *
 *     OCaml interface to Sundials (serial) CVODE and IDA solvers      *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) * 
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

#include <sundials/sundials_config.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_spgmr.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

#include "spils_ml.h"
#include "sundials_ml.h"
#include "nvector_ml.h"

// Call with SPILS_ML_BIGARRAYS to compile for the Serial NVector to
// Bigarray interface code.

#ifdef SPILS_ML_BIGARRAYS

#define CVTYPE(fname) c_ba_spils_ ## fname
#include <nvector/nvector_serial.h>

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0

#define NVECTORIZE_VAL(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDVAL(nv) N_VDestroy(nv)

#define NVECTOR_APPROXSIZE(ba) (Caml_ba_array_val(ba)->dim[0])

#else

#define CVTYPE(fname) c_nvec_spils_ ## fname
#include <sundials/sundials_nvector.h>

#define WRAP_NVECTOR(v) NVEC_DATA(v)
#define RELINQUISH_WRAPPEDNV(v) {}

#define NVECTORIZE_VAL(v) NVEC_VAL(v)
#define RELINQUISH_NVECTORIZEDVAL(nv) {}

// we only use the size in caml_alloc_final
#define NVECTOR_APPROXSIZE(v) (1)

#endif

CAMLprim value CVTYPE(modified_gs)(value vv, value vh, value vk, value vp)
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
    if (caml_array_length(vv) < k)
	caml_invalid_argument("Spils.modified_gs: v is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_constant(*caml_named_value("spils_MemoryRequestFailure"));
    for (i = 0; i <= p; ++i) {
	v[i] = NVECTORIZE_VAL(Field(vv, k - p + i));
    }

    ModifiedGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm);

    for (i = 0; i <= p; ++i) {
	RELINQUISH_NVECTORIZEDVAL(v[i]);
    }
    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value CVTYPE(classical_gs)(value vargs)
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
    if (caml_array_length(vv) < k)
	caml_invalid_argument("Spils.classical_gs: v is too small.");
    if (ARRAY1_LEN(vs) < k)
	caml_invalid_argument("Spils.classical_gs: s is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_constant(*caml_named_value("spils_MemoryRequestFailure"));
    for (i = 0; i <= p; ++i) {
	v[i] = NVECTORIZE_VAL(Field(vv, k - p + i));
    }

    temp = NVECTORIZE_VAL(Field(vargs, 4));
    ClassicalGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm,
	        temp, REAL_ARRAY(vs));
    RELINQUISH_NVECTORIZEDVAL(temp);

    for (i = 0; i <= p; ++i) {
	RELINQUISH_NVECTORIZEDVAL(v[i]);
    }
    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

static int atimes_f(void *a_fn, N_Vector v, N_Vector z)
{
    CAMLparam0();
    CAMLlocal3(vv, vz, res);
    int r = 0;

    vv = WRAP_NVECTOR(v);
    vz = WRAP_NVECTOR(z);

    // the data payloads inside vv and vz are only valid during this call,
    // afterward that memory goes back to Sundials. These bigarrays must not
    // be retained by a_fn! If it wants a permanent copy, then it has to
    // make it manually.
    res = caml_callback2((value)a_fn, vv, vz);
    if (Is_exception_result(res)) r = 1;

    RELINQUISH_WRAPPEDNV(vv);
    RELINQUISH_WRAPPEDNV(vz);

    CAMLreturnT(int, r);
}

static int psolve_f(void *p_fn, N_Vector r, N_Vector z, int lr)
{
    CAMLparam0();
    CAMLlocal3(vr, vz, res);
    int fr = 0;

    vr = WRAP_NVECTOR(r);
    vz = WRAP_NVECTOR(z);

    // the data payloads inside vr and vz are only valid during this call,
    // afterward that memory goes back to Sundials. These bigarrays must not
    // be retained by a_fn! If it wants a permanent copy, then it has to
    // make it manually.
    res = caml_callback3((value)p_fn, vr, vz, Val_bool(lr == 1));
    if (Is_exception_result(res)) fr = 1;

    RELINQUISH_WRAPPEDNV(vr);
    RELINQUISH_WRAPPEDNV(vz);

    CAMLreturnT(int, fr);
}

/* SPGMR */

#define SPGMR_SESSION(v) ((SpgmrMem)Data_custom_val(v))

CAMLprim void spgmr_finalize(value vs)
{
    SpgmrMem s = SPGMR_SESSION(vs);
    SpgmrFree(s);
}

CAMLprim value CVTYPE(spgmr_make)(value vl_max, value vvec_tmpl)
{
    CAMLparam2(vl_max, vvec_tmpl);
    CAMLlocal1(vs);
    SpgmrMem s;
    N_Vector vec_tmpl;
    mlsize_t approx_size = sizeof(SpgmrMemRec)
	+ (Int_val(vl_max) + 3) * NVECTOR_APPROXSIZE(vvec_tmpl);

    vec_tmpl = NVECTORIZE_VAL(vvec_tmpl);
    s = SpgmrMalloc(Int_val(vl_max), vec_tmpl);
    RELINQUISH_NVECTORIZEDVAL(vec_tmpl);

    if (s == NULL) caml_raise_out_of_memory();
    vs = caml_alloc_final(1, &spgmr_finalize, approx_size, approx_size * 20);

    CAMLreturn(vs);
}

CAMLprim value CVTYPE(spgmr_solve)(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal4(vr, vpsolve, vs1, vs2);
    int flag;
    int success;
    N_Vector x = NVECTORIZE_VAL(Field(vargs, 1));
    N_Vector b = NVECTORIZE_VAL(Field(vargs, 2));
    N_Vector s1 = NULL;
    N_Vector s2 = NULL;
    realtype res_norm = 0.0;
    int nli = 0;
    int nps = 0;
    void *p_data = NULL;

    vpsolve = Field(vargs, 10);
    if (vpsolve != Val_none) p_data = (void *)Field(vpsolve, 0);

    vs1 = Field(vargs, 7);
    if (vs1 != Val_none) s1 = NVECTORIZE_VAL(Field(vs1, 0));

    vs2 = Field(vargs, 8);
    if (vs2 != Val_none) s2 = NVECTORIZE_VAL(Field(vs2, 0));

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

    RELINQUISH_NVECTORIZEDVAL(x);
    RELINQUISH_NVECTORIZEDVAL(b);
    if (s1 != NULL) RELINQUISH_NVECTORIZEDVAL(s1);
    if (s2 != NULL) RELINQUISH_NVECTORIZEDVAL(s2);

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

// TODO: Adapt the above for SPBCG
// TODO: Adapt the above for SPTFQMR

