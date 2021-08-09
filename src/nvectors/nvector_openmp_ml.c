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
#include "../nvectors/nvector_ml.h"
#include "../nvectors/nvector_openmp_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include <nvector/nvector_openmp.h>

/* Adapted from sundials-2.6.1/src/nvec_pthreads/nvector_openmp.c:
   N_VCloneEmpty_OpenMP */
static N_Vector clone_openmp(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);

    N_Vector v;
    N_VectorContent_OpenMP content;

    if (w == NULL) CAMLreturnT(N_Vector, NULL);
    w_payload = NVEC_BACKLINK(w);
    struct caml_ba_array *w_ba = Caml_ba_array_val(w_payload);

    /* Create vector (we need not copy the data) */
    v_payload = caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim);

    v = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_OpenMP), v_payload);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_OpenMP) v->content;

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

    /* Create content */
    content->length      = NV_LENGTH_OMP(w);
    content->num_threads = NV_NUM_THREADS_OMP(w);
    content->own_data    = 0;
    content->data        = Caml_ba_data_val(v_payload);

    CAMLreturnT(N_Vector, v);
}

/* Clone an "any" nvector by unwrapping and wrapping the RA payload. */
/* Adapted from sundials-2.6.1/src/nvec_pthreads/nvector_openmp.c:
   N_VCloneEmpty_OpenMP */
static N_Vector clone_any_openmp(N_Vector w)
{
    CAMLparam0();
    CAMLlocal4(v_wrapped, v_payload, w_wrapped, w_payload);

    N_Vector v;
    N_VectorContent_OpenMP content;

    if (w == NULL) CAMLreturnT(N_Vector, NULL);
    w_wrapped = NVEC_BACKLINK(w);
    w_payload = Field(w_wrapped, 1);

    struct caml_ba_array *w_ba = Caml_ba_array_val(w_payload);

    /* Create vector (we need not copy the data) */
    v_payload = caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim);

    v_wrapped = caml_alloc_tuple(2);
    Store_field(v_wrapped, 0, Field(w_wrapped, 0)); // RA constructor
    Store_field(v_wrapped, 1, v_payload);

    v = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_OpenMP), v_wrapped);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_OpenMP) v->content;

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

    /* Create content */
    content->length   = NV_LENGTH_OMP(w);
    content->num_threads = NV_NUM_THREADS_OMP(w);
    content->own_data = 0;
    content->data     = Caml_ba_data_val(v_payload);

    CAMLreturnT(N_Vector, v);
}

/*
 * N_VCloneEmpty is used in Sundials as a "light-weight" way to wrap
 * array data for use in calculations with N_Vectors. At the time of
 * writing (Sundials 3.1.0), this feature is only used in the *DenseDQJac
 * routines (e.g., cvDlsDenseDQJac) and the cloned N_Vectors are never
 * passed back into OCaml. So, we do not bother to reproduce the backlink.
 *
 * If the use of this feature is generalized, it will be necessary to
 * duplicate the backlink, and to override the N_VSetArrayPointer function
 * to update the backlink appropriately. Currently, N_VSetArrayPointer is
 * only used in the *DenseDQJac functions and in the *_bbdpre
 * implementations (which only use N_Vectors locally and never pass them
 * back into OCaml).
 */
static N_Vector clone_empty_openmp(N_Vector w)
{
    N_Vector v;

    v = N_VCloneEmpty_OpenMP(w);
    v->ops->nvdestroy = N_VDestroy_OpenMP;
    v->ops->nvclone = NULL;
    v->ops->nvcloneempty = NULL;

    return v;
}

/* Creation from OCaml.  */
/* Adapted from sundials-2.6.1/src/nvec_openmp/nvector_openmp.c:
   N_VNewEmpty_OpenMP */
CAMLprim value sunml_nvec_wrap_openmp(value nthreads,
				   value payload, value checkfn)
{
    CAMLparam3(nthreads, payload, checkfn);
    CAMLlocal1(vnvec);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_OpenMP content;
    long int length = (Caml_ba_array_val(payload))->dim[0];

    /* Create vector */
    nv = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_OpenMP), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_OpenMP) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_openmp;		    /* ours */
    ops->nvcloneempty      = clone_empty_openmp;	    /* ours */
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = sunml_free_cnvec;
#if 270 <= SUNDIALS_LIB_VERSION
    ops->nvgetvectorid	   = N_VGetVectorID_OpenMP;
#endif

    ops->nvspace           = N_VSpace_OpenMP;		    /* theirs */
    ops->nvgetarraypointer = N_VGetArrayPointer_OpenMP;
    ops->nvsetarraypointer = N_VSetArrayPointer_OpenMP;
    ops->nvlinearsum       = N_VLinearSum_OpenMP;
    ops->nvconst           = N_VConst_OpenMP;
    ops->nvprod            = N_VProd_OpenMP;
    ops->nvdiv             = N_VDiv_OpenMP;
    ops->nvscale           = N_VScale_OpenMP;
    ops->nvabs             = N_VAbs_OpenMP;
    ops->nvinv             = N_VInv_OpenMP;
    ops->nvaddconst        = N_VAddConst_OpenMP;
    ops->nvdotprod         = N_VDotProd_OpenMP;
    ops->nvmaxnorm         = N_VMaxNorm_OpenMP;
    ops->nvwrmsnormmask    = N_VWrmsNormMask_OpenMP;
    ops->nvwrmsnorm        = N_VWrmsNorm_OpenMP;
    ops->nvmin             = N_VMin_OpenMP;
    ops->nvwl2norm         = N_VWL2Norm_OpenMP;
    ops->nvl1norm          = N_VL1Norm_OpenMP;
    ops->nvcompare         = N_VCompare_OpenMP;
    ops->nvinvtest         = N_VInvTest_OpenMP;
    ops->nvconstrmask      = N_VConstrMask_OpenMP;
    ops->nvminquotient     = N_VMinQuotient_OpenMP;

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
    ops->nvgetlength	    = N_VGetLength_OpenMP;
    ops->nvgetcommunicator  = NULL;

    ops->nvdotprodlocal     = N_VDotProd_OpenMP;
    ops->nvmaxnormlocal     = N_VMaxNorm_OpenMP;
    ops->nvminlocal         = N_VMin_OpenMP;
    ops->nvl1normlocal      = N_VL1Norm_OpenMP;
    ops->nvinvtestlocal     = N_VInvTest_OpenMP;
    ops->nvconstrmasklocal  = N_VConstrMask_OpenMP;
    ops->nvminquotientlocal = N_VMinQuotient_OpenMP;
    ops->nvwsqrsumlocal     = N_VWSqrSumLocal_OpenMP;
    ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_OpenMP;
#endif

    /* Create content */
    content->length      = length;
    content->num_threads = Int_val(nthreads);
    content->own_data    = 0;
    content->data        = Caml_ba_data_val(payload);

    vnvec = caml_alloc_tuple(3);
    Store_field(vnvec, 0, payload);
    Store_field(vnvec, 1, sunml_alloc_caml_nvec(nv, sunml_finalize_caml_nvec));
    Store_field(vnvec, 2, checkfn);

    CAMLreturn(vnvec);
}

/* The "any"-version of an openmp nvector is created by modifying the
   standard one in two ways:
   1. The payload field is wrapped in the RA constructor.
   2. The nvclone operation is overridden to implement the wrapping operation
      (the current clone_empty_serial does not manipulate the backlink). */
CAMLprim value sunml_nvec_anywrap_openmp(value extconstr,
					 value nthreads,
					 value payload, value checkfn)
{
    CAMLparam4(extconstr, nthreads, payload, checkfn);
    CAMLlocal2(vnv, vwrapped);
    N_Vector nv;
    N_Vector_Ops ops;

    vnv = sunml_nvec_wrap_openmp(nthreads, payload, checkfn);
    nv = NVEC_VAL(vnv);
    ops = (N_Vector_Ops) nv->ops;

    ops->nvclone = clone_any_openmp;

    vwrapped = caml_alloc_tuple(2);
    Store_field(vwrapped, 0, extconstr);
    Store_field(vwrapped, 1, NVEC_BACKLINK(nv));

    Store_field(vnv, 0, vwrapped);
    NVEC_BACKLINK(nv) = vwrapped;

    CAMLreturn(vnv);
}

CAMLprim value sunml_nvec_openmp_num_threads(value va)
{
    CAMLparam1(va);
    int num_threads = NV_NUM_THREADS_OMP(NVEC_VAL(va));

    CAMLreturn(Val_int(num_threads));
}

/** Interface to underlying openmp nvector functions */

CAMLprim value sunml_nvec_openmp_n_vlinearsum(value va, value vx, value vb,
					      value vy, value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(y) != NV_LENGTH_OMP(x)
	    || NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vlinearsum");
#endif

    N_VLinearSum_OpenMP(Double_val(va), x, Double_val(vb), y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vconst(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst_OpenMP(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vprod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(y) != NV_LENGTH_OMP(x)
	    || NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vprod");
#endif

    N_VProd_OpenMP(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(y) != NV_LENGTH_OMP(x)
	    || NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vdiv");
#endif

    N_VDiv_OpenMP(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vscale");
#endif

    N_VScale_OpenMP(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vabs");
#endif

    N_VAbs_OpenMP(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vinv");
#endif

    N_VInv_OpenMP(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vaddconst");
#endif

    N_VAddConst_OpenMP(x, Double_val(vb), z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(y) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vdotprod");
#endif

    realtype r = N_VDotProd_OpenMP(x, y);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vmaxnorm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMaxNorm_OpenMP(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vwrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(w) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vwrmsnorm");
#endif

    realtype r = N_VWrmsNorm_OpenMP(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(w) != NV_LENGTH_OMP(x)
	    || NV_LENGTH_OMP(w) != NV_LENGTH_OMP(id))
	caml_invalid_argument("Nvector_openmp.n_vwrmsnormmask");
#endif

    realtype r = N_VWrmsNormMask_OpenMP(x, w, id);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vmin(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMin_OpenMP(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vwl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(w) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vwl2norm");
#endif

    realtype r = N_VWL2Norm_OpenMP(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vl1norm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VL1Norm_OpenMP(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vcompare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vcompare");
#endif

    N_VCompare_OpenMP(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(z) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vinvtest");
#endif

    booleantype r = N_VInvTest_OpenMP(x, z);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_openmp_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    N_Vector c = NVEC_VAL(vc);
    N_Vector x = NVEC_VAL(vx);
    N_Vector m = NVEC_VAL(vm);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(x) != NV_LENGTH_OMP(c)
	    || NV_LENGTH_OMP(m) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vconstrmask");
#endif

    booleantype r = N_VConstrMask_OpenMP(c, x, m);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_openmp_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    N_Vector num = NVEC_VAL(vnum);
    N_Vector denom = NVEC_VAL(vdenom);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(num) != NV_LENGTH_OMP(denom))
	caml_invalid_argument("Nvector_openmp.n_vminquotient");
#endif

    realtype r = N_VMinQuotient_OpenMP(num, denom);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vspace(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
    sundials_ml_index lrw, liw;

    N_VSpace_OpenMP(NVEC_VAL(vx), &lrw, &liw);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_index(lrw));
    Store_field(r, 1, Val_index(liw));

    CAMLreturn(r);
}

CAMLprim value sunml_nvec_openmp_n_vgetlength(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
#if 500 <= SUNDIALS_LIB_VERSION
    r = Val_int(N_VGetLength_OpenMP(NVEC_VAL(vx)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

/* fused vector operations */

CAMLprim value sunml_nvec_openmp_n_vlinearcombination(value vac, value vax,
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
	    || NV_LENGTH_OMP(ax[0]) != NV_LENGTH_OMP(z))
	caml_invalid_argument("Nvector_serial.n_vlinearcombination");
#endif

    N_VLinearCombination_OpenMP(nvec, ac, ax, z);
    free(ax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vscaleaddmulti(value vac, value vx, value vay,
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
	    || NV_LENGTH_OMP(a[0][0]) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_serial.n_vscaleaddmulti");
#endif

    N_VScaleAddMulti_OpenMP(nvec, ac, x, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vdotprodmulti(value vx, value vay, value vad)
{
    CAMLparam3(vx, vay, vad);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ad = REAL_ARRAY(vad);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *ay;
    int nvec = sunml_arrays_of_nvectors(&ay, 1, vay);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vad) < nvec
		|| NV_LENGTH_OMP(ay[0]) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_serial.n_vdotprodmulti");
#endif

    N_VDotProdMulti_OpenMP(nvec, x, ay, ad);
    free(ay);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}


/* vector array operations */

CAMLprim value sunml_nvec_openmp_n_vlinearsumvectorarray(value va, value vax,
							 value vb, value vay,
							 value vaz)
{
    CAMLparam5(va, vax, vb, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *a[3];
    int nvec = sunml_arrays_of_nvectors(a, 3, vax, vay, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_openmp.n_vlinearsumvectorarray");
#endif

    N_VLinearSumVectorArray_OpenMP(nvec, Double_val(va), a[0],
					 Double_val(vb), a[1], a[2]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vscalevectorarray(value vac, value vax,
						     value vaz)
{
    CAMLparam3(vac, vax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vac) < nvec)
	caml_invalid_argument("Nvector_openmp.n_vscalevectorarray");
#endif

    N_VScaleVectorArray_OpenMP(nvec, ac, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vconstvectorarray(value vc, value vaz)
{
    CAMLparam2(vc, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *az;
    int nvec = sunml_arrays_of_nvectors(&az, 1, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_openmp.n_vconstvectorarray");
#endif

    N_VConstVectorArray_OpenMP(nvec, Double_val(vc), az);
    free(az);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vwrmsnormvectorarray(value vax, value vaw,
						        value van)
{
    CAMLparam3(vax, vaw, van);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *an = REAL_ARRAY(van);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(van) < nvec)
	caml_invalid_argument("Nvector_openmp.n_vconstvectorarray");
#endif

    N_VWrmsNormVectorArray_OpenMP(nvec, a[0], a[1], an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vwrmsnormmaskvectorarray(value vax,
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
		|| NV_LENGTH_OMP(i) != NV_LENGTH_OMP(a[0][0]))
	caml_invalid_argument("Nvector_openmp.n_vconstvectorarray");
#endif

    N_VWrmsNormMaskVectorArray_OpenMP(nvec, a[0], a[1], i, an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vscaleaddmultivectorarray(value vaa,
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
	caml_invalid_argument("Nvector_openmp.n_vscaleaddmultivectorarray");
    }
#endif

    N_VScaleAddMultiVectorArray_OpenMP(nvec, nsum, aa, ax, ayz[0], ayz[1]);
    free(ax);
    free(*ayz);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vlinearcombinationvectorarray(value vac,
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
	caml_invalid_argument("Nvector_openmp.n_vlinearcombinationvectorarray");
    }
#endif

    N_VLinearCombinationVectorArray_OpenMP(nvec, nsum, ac, aax, az);
    free(aax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

/** Reduce operations for openmp nvectors */

CAMLprim value sunml_nvec_openmp_n_vwsqrsumlocal(value vx, value vw)
{
    CAMLparam2(vx, vw);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(w) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vwsqrsumlocal");
#endif

    r = N_VWSqrSumLocal_OpenMP(x, w);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vwsqrsummasklocal(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_OMP(w) != NV_LENGTH_OMP(x)
	    || NV_LENGTH_OMP(id) != NV_LENGTH_OMP(x))
	caml_invalid_argument("Nvector_openmp.n_vwsqrsummasklocal");
#endif

    r = N_VWSqrSumMaskLocal_OpenMP(x, w, id);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

/** Selectively activate fused and array operations for serial nvectors */

CAMLprim value sunml_nvec_openmp_enablefusedops(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableFusedOps_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablelinearcombination(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombination_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablescaleaddmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMulti_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enabledotprodmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableDotProdMulti_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablelinearsumvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearSumVectorArray_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablescalevectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleVectorArray_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enableconstvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableConstVectorArray_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablewrmsnormvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormVectorArray_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablewrmsnormmaskvectorarray(value vx,
							       value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormMaskVectorArray_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablescaleaddmultivectorarray(value vx,
								value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMultiVectorArray_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_enablelinearcombinationvectorarray(value vx,
								    value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombinationVectorArray_OpenMP(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

