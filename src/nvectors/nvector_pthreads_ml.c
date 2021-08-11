/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2015 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include "../sundials/sundials_ml.h"
#include "nvector_ml.h"
#include "nvector_pthreads_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include <nvector/nvector_pthreads.h>

/* Adapted from sundials-2.6.1/src/nvec_pthreads/nvector_pthreads.c:
   N_VCloneEmpty_Pthreads */
static N_Vector clone_pthreads(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);

    N_Vector v;
    N_VectorContent_Pthreads content;

    if (w == NULL) CAMLreturnT(N_Vector, NULL);
    w_payload = NVEC_BACKLINK(w);
    struct caml_ba_array *w_ba = Caml_ba_array_val(w_payload);

    /* Create vector (we need not copy the data) */
    v_payload = caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim);

    v = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Pthreads), v_payload);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_Pthreads) v->content;

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

    /* Create content */
    content->length      = NV_LENGTH_PT(w);
    content->num_threads = NV_NUM_THREADS_PT(w);
    content->own_data    = 0;
    content->data        = Caml_ba_data_val(v_payload);

    CAMLreturnT(N_Vector, v);
}

/* Clone an "any" nvector by unwrapping and wrapping the RA payload. */
/* Adapted from sundials-2.6.1/src/nvec_pthreads/nvector_pthreads.c:
   N_VCloneEmpty_Pthreads */
static N_Vector clone_any_pthreads(N_Vector w)
{
    CAMLparam0();
    CAMLlocal4(v_wrapped, v_payload, w_wrapped, w_payload);

    N_Vector v;
    N_VectorContent_Pthreads content;

    if (w == NULL) CAMLreturnT(N_Vector, NULL);
    w_wrapped = NVEC_BACKLINK(w);
    w_payload = Field(w_wrapped, 1);

    struct caml_ba_array *w_ba = Caml_ba_array_val(w_payload);

    /* Create vector (we need not copy the data) */
    v_payload = caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim);

    v_wrapped = caml_alloc_tuple(2);
    Store_field(v_wrapped, 0, Field(w_wrapped, 0)); // RA constructor
    Store_field(v_wrapped, 1, v_payload);

    v = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Pthreads), v_wrapped);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_Pthreads) v->content;

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

    /* Create content */
    content->length   = NV_LENGTH_PT(w);
    content->num_threads = NV_NUM_THREADS_PT(w);
    content->own_data = 0;
    content->data     = Caml_ba_data_val(v_payload);

    CAMLreturnT(N_Vector, v);
}

/* Creation from OCaml.  */
/* Adapted from sundials-2.6.1/src/nvec_pthreads/nvector_pthreads.c:
   N_VNewEmpty_Pthreads */
CAMLprim value sunml_nvec_wrap_pthreads(value nthreads,
				        value payload,
				        value checkfn, value clonefn)
{
    CAMLparam4(nthreads, payload, checkfn, clonefn);
    CAMLlocal1(vnvec);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Pthreads content;
    long int length = (Caml_ba_array_val(payload))->dim[0];

    /* Create vector */
    nv = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Pthreads), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_Pthreads) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_pthreads;		    /* ours */
    ops->nvcloneempty      = NULL;
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = sunml_free_cnvec;
#if 270 <= SUNDIALS_LIB_VERSION
    ops->nvgetvectorid	   = N_VGetVectorID_Pthreads;
#endif

    ops->nvspace           = N_VSpace_Pthreads;		    /* theirs */
    ops->nvgetarraypointer = N_VGetArrayPointer_Pthreads;
    ops->nvsetarraypointer = N_VSetArrayPointer_Pthreads;
    ops->nvlinearsum       = N_VLinearSum_Pthreads;
    ops->nvconst           = N_VConst_Pthreads;
    ops->nvprod            = N_VProd_Pthreads;
    ops->nvdiv             = N_VDiv_Pthreads;
    ops->nvscale           = N_VScale_Pthreads;
    ops->nvabs             = N_VAbs_Pthreads;
    ops->nvinv             = N_VInv_Pthreads;
    ops->nvaddconst        = N_VAddConst_Pthreads;
    ops->nvdotprod         = N_VDotProd_Pthreads;
    ops->nvmaxnorm         = N_VMaxNorm_Pthreads;
    ops->nvwrmsnormmask    = N_VWrmsNormMask_Pthreads;
    ops->nvwrmsnorm        = N_VWrmsNorm_Pthreads;
    ops->nvmin             = N_VMin_Pthreads;
    ops->nvwl2norm         = N_VWL2Norm_Pthreads;
    ops->nvl1norm          = N_VL1Norm_Pthreads;
    ops->nvcompare         = N_VCompare_Pthreads;
    ops->nvinvtest         = N_VInvTest_Pthreads;
    ops->nvconstrmask      = N_VConstrMask_Pthreads;
    ops->nvminquotient     = N_VMinQuotient_Pthreads;

#if 400 <= SUNDIALS_LIB_VERSION
    /* fused vector operations (optional, NULL means disabled by default) */
    ops->nvlinearcombination = NULL;
    ops->nvscaleaddmulti     = NULL;
    ops->nvdotprodmulti      = NULL;

    /* vector array operations (optional, NULL means disabled by default) */
    ops->nvlinearsumvectorarray         = NULL;
    ops->nvscalevectorarray             = NULL;
    ops->nvconstvectorarray             = NULL;
    ops->nvwrmsnormvectorarray          = NULL;
    ops->nvwrmsnormmaskvectorarray      = NULL;
    ops->nvscaleaddmultivectorarray     = NULL;
    ops->nvlinearcombinationvectorarray = NULL;
#endif

#if 500 <= SUNDIALS_LIB_VERSION
    ops->nvgetlength	    = N_VGetLength_Pthreads;
    ops->nvgetcommunicator  = NULL;

    ops->nvdotprodlocal     = N_VDotProd_Pthreads;
    ops->nvmaxnormlocal     = N_VMaxNorm_Pthreads;
    ops->nvminlocal         = N_VMin_Pthreads;
    ops->nvl1normlocal      = N_VL1Norm_Pthreads;
    ops->nvinvtestlocal     = N_VInvTest_Pthreads;
    ops->nvconstrmasklocal  = N_VConstrMask_Pthreads;
    ops->nvminquotientlocal = N_VMinQuotient_Pthreads;
    ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Pthreads;
    ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Pthreads;
#endif

    /* Create content */
    content->length      = length;
    content->num_threads = Int_val(nthreads);
    content->own_data    = 0;
    content->data        = Caml_ba_data_val(payload);

    vnvec = NVEC_ALLOC();
    Store_field(vnvec, NVEC_PAYLOAD, payload);
    Store_field(vnvec, NVEC_CPTR,
		sunml_alloc_caml_nvec(nv, sunml_finalize_caml_nvec));
    Store_field(vnvec, NVEC_CHECK, checkfn);
    Store_field(vnvec, NVEC_CLONE, clonefn);

    CAMLreturn(vnvec);
}

/* The "any"-version of a pthreads nvector is created by modifying the
   standard one in two ways:
   1. The payload field is wrapped in the RA constructor.
   2. The nvclone operation is overridden to implement the wrapping operation
      (the current clone_empty_pthreads does not manipulate the backlink). */
CAMLprim value sunml_nvec_anywrap_pthreads(value extconstr,
					   value nthreads, value payload,
					   value checkfn, value clonefn)
{
    CAMLparam5(extconstr, nthreads, payload, checkfn, clonefn);
    CAMLlocal2(vnv, vwrapped);
    N_Vector nv;
    N_Vector_Ops ops;

    vnv = sunml_nvec_wrap_pthreads(nthreads, payload, checkfn, clonefn);
    nv = NVEC_VAL(vnv);
    ops = (N_Vector_Ops) nv->ops;

    ops->nvclone = clone_any_pthreads;

    vwrapped = caml_alloc_tuple(2);
    Store_field(vwrapped, 0, extconstr);
    Store_field(vwrapped, 1, NVEC_BACKLINK(nv));

    Store_field(vnv, 0, vwrapped);
    NVEC_BACKLINK(nv) = vwrapped;

    CAMLreturn(vnv);
}

CAMLprim value sunml_nvec_pthreads_num_threads(value va)
{
    CAMLparam1(va);
    int num_threads = NV_NUM_THREADS_PT(NVEC_VAL(va));

    CAMLreturn(Val_int(num_threads));
}

/** Interface to underlying pthreads nvector functions */

CAMLprim value sunml_nvec_pthreads_n_vlinearsum(value va, value vx, value vb,
					     value vy, value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(y) != NV_LENGTH_PT(x) || NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vlinearsum");
#endif

    N_VLinearSum_Pthreads(Double_val(va), x, Double_val(vb), y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vconst(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst_Pthreads(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vprod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(y) != NV_LENGTH_PT(x) || NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vprod");
#endif

    N_VProd_Pthreads(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(y) != NV_LENGTH_PT(x) || NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vdiv");
#endif

    N_VDiv_Pthreads(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vscale");
#endif

    N_VScale_Pthreads(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vabs");
#endif

    N_VAbs_Pthreads(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vinv");
#endif

    N_VInv_Pthreads(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vaddconst");
#endif

    N_VAddConst_Pthreads(x, Double_val(vb), z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(y) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vdotprod");
#endif

    realtype r = N_VDotProd_Pthreads(x, y);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vmaxnorm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMaxNorm_Pthreads(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vwrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(w) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vwrmsnorm");
#endif

    realtype r = N_VWrmsNorm_Pthreads(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(w) != NV_LENGTH_PT(x) || NV_LENGTH_PT(w) != NV_LENGTH_PT(id))
	caml_invalid_argument("Nvector_pthreads.n_vwrmsnormmask");
#endif

    realtype r = N_VWrmsNormMask_Pthreads(x, w, id);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vmin(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMin_Pthreads(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vwl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(w) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vwl2norm");
#endif

    realtype r = N_VWL2Norm_Pthreads(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vl1norm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VL1Norm_Pthreads(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vcompare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vcompare");
#endif

    N_VCompare_Pthreads(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(z) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vinvtest");
#endif

    booleantype r = N_VInvTest_Pthreads(x, z);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_pthreads_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    N_Vector c = NVEC_VAL(vc);
    N_Vector x = NVEC_VAL(vx);
    N_Vector m = NVEC_VAL(vm);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(x) != NV_LENGTH_PT(c) || NV_LENGTH_PT(m) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vconstrmask");
#endif

    booleantype r = N_VConstrMask_Pthreads(c, x, m);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_pthreads_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    N_Vector num = NVEC_VAL(vnum);
    N_Vector denom = NVEC_VAL(vdenom);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(num) != NV_LENGTH_PT(denom))
	caml_invalid_argument("Nvector_pthreads.n_vminquotient");
#endif

    realtype r = N_VMinQuotient_Pthreads(num, denom);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vspace(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
    sundials_ml_index lrw, liw;

    N_VSpace_Pthreads(NVEC_VAL(vx), &lrw, &liw);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_index(lrw));
    Store_field(r, 1, Val_index(liw));

    CAMLreturn(r);
}

CAMLprim value sunml_nvec_pthreads_n_vgetlength(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
#if 500 <= SUNDIALS_LIB_VERSION
    r = Val_int(N_VGetLength_Pthreads(NVEC_VAL(vx)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

/* fused vector operations */

CAMLprim value sunml_nvec_pthreads_n_vlinearcombination(value vac, value vax,
						   value vz)
{
    CAMLparam3(vac, vax, vz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector z = NVEC_VAL(vz);
    N_Vector *ax;
    int nvec = sunml_arrays_of_nvectors(&ax, 1, vax);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vac) < nvec
	    || NV_LENGTH_PT(ax[0]) != NV_LENGTH_PT(z))
	caml_invalid_argument("Nvector_pthreads.n_vlinearcombination");
#endif

    N_VLinearCombination_Pthreads(nvec, ac, ax, z);
    free(ax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vscaleaddmulti(value vac, value vx, value vay,
					            value vaz)
{
    CAMLparam4(vac, vx, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vay, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vac) < nvec
	    || NV_LENGTH_PT(a[0][0]) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vscaleaddmulti");
#endif

    N_VScaleAddMulti_Pthreads(nvec, ac, x, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vdotprodmulti(value vx, value vay, value vad)
{
    CAMLparam3(vx, vay, vad);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ad = REAL_ARRAY(vad);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *ay;
    int nvec = sunml_arrays_of_nvectors(&ay, 1, vay);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vad) < nvec
		|| NV_LENGTH_PT(ay[0]) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vdotprodmulti");
#endif

    N_VDotProdMulti_Pthreads(nvec, x, ay, ad);
    free(ay);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}


/* vector array operations */

CAMLprim value sunml_nvec_pthreads_n_vlinearsumvectorarray(value va, value vax,
							   value vb, value vay,
							   value vaz)
{
    CAMLparam5(va, vax, vb, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *a[3];
    int nvec = sunml_arrays_of_nvectors(a, 3, vax, vay, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_pthreads.n_vlinearsumvectorarray");
#endif

    N_VLinearSumVectorArray_Pthreads(nvec, Double_val(va), a[0],
					 Double_val(vb), a[1], a[2]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vscalevectorarray(value vac, value vax,
						       value vaz)
{
    CAMLparam3(vac, vax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vac) < nvec)
	caml_invalid_argument("Nvector_pthreads.n_vscalevectorarray");
#endif

    N_VScaleVectorArray_Pthreads(nvec, ac, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vconstvectorarray(value vc, value vaz)
{
    CAMLparam2(vc, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *az;
    int nvec = sunml_arrays_of_nvectors(&az, 1, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_pthreads.n_vconstvectorarray");
#endif

    N_VConstVectorArray_Pthreads(nvec, Double_val(vc), az);
    free(az);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vwrmsnormvectorarray(value vax, value vaw,
							  value van)
{
    CAMLparam3(vax, vaw, van);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *an = REAL_ARRAY(van);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(van) < nvec)
	caml_invalid_argument("Nvector_pthreads.n_vconstvectorarray");
#endif

    N_VWrmsNormVectorArray_Pthreads(nvec, a[0], a[1], an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vwrmsnormmaskvectorarray(value vax,
						value vaw, value vi, value van)
{
    CAMLparam4(vax, vaw, vi, van);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *an = REAL_ARRAY(van);
    N_Vector i = NVEC_VAL(vi);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(van) < nvec
		|| NV_LENGTH_PT(i) != NV_LENGTH_PT(a[0][0]))
	caml_invalid_argument("Nvector_pthreads.n_vconstvectorarray");
#endif

    N_VWrmsNormMaskVectorArray_Pthreads(nvec, a[0], a[1], i, an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vscaleaddmultivectorarray(value vaa,
					    value vax, value vaay, value vaaz)
{
    CAMLparam4(vaa, vax, vaay, vaaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *aa = REAL_ARRAY(vaa);
    N_Vector *ax = NULL;
    int nvec = sunml_arrays_of_nvectors(&ax, 1, vax);
    N_Vector **ayz[2] = { NULL };
    int nvec2, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec2, ayz, 2, vaay, vaaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || !nsum || nvec2 != nvec || ARRAY1_LEN(vaa) < nsum) {
	if (ax != NULL) free(ax);
	if (*ayz != NULL) free(*ayz);
	caml_invalid_argument("Nvector_pthreads.n_vscaleaddmultivectorarray");
    }
#endif

    N_VScaleAddMultiVectorArray_Pthreads(nvec, nsum, aa, ax, ayz[0], ayz[1]);
    free(ax);
    free(*ayz);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_pthreads_n_vlinearcombinationvectorarray(value vac,
							value vaax, value vaz)
{
    CAMLparam3(vac, vaax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *az;
    int nvecz __attribute__((unused))
	= sunml_arrays_of_nvectors(&az, 1, vaz);
    N_Vector **aax;
    int nvec, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec, &aax, 1, vaax);

#if SUNDIALS_ML_SAFE == 1
    if (!nvecz || !nsum || nvec > nvecz || ARRAY1_LEN(vac) < nsum) {
	if (az != NULL) free(az);
	if (aax != NULL) free(aax);
	caml_invalid_argument("Nvector_pthreads.n_vlinearcombinationvectorarray");
    }
#endif

    N_VLinearCombinationVectorArray_Pthreads(nvec, nsum, ac, aax, az);
    free(aax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

/** Reduce operations for pthreads nvectors */

CAMLprim value sunml_nvec_pthreads_n_vwsqrsumlocal(value vx, value vw)
{
    CAMLparam2(vx, vw);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(w) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vwsqrsumlocal");
#endif

    r = N_VWSqrSumLocal_Pthreads(x, w);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_pthreads_n_vwsqrsummasklocal(value vx, value vw,
						       value vid)
{
    CAMLparam3(vx, vw, vid);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_PT(w) != NV_LENGTH_PT(x)
	    || NV_LENGTH_PT(id) != NV_LENGTH_PT(x))
	caml_invalid_argument("Nvector_pthreads.n_vwsqrsummasklocal");
#endif

    r = N_VWSqrSumMaskLocal_Pthreads(x, w, id);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

/** Selectively activate fused and array operations for serial nvectors */

CAMLprim value sunml_nvec_pthreads_enablefusedops(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableFusedOps_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablelinearcombination(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombination_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablescaleaddmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMulti_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enabledotprodmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableDotProdMulti_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablelinearsumvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearSumVectorArray_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablescalevectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleVectorArray_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enableconstvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableConstVectorArray_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablewrmsnormvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormVectorArray_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablewrmsnormmaskvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormMaskVectorArray_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablescaleaddmultivectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMultiVectorArray_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_pthreads_enablelinearcombinationvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombinationVectorArray_Pthreads(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

