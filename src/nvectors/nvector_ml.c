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

#include "../nvectors/nvector_ml.h"
#include "../sundials/sundials_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include <math.h>		/* for nan() */
#include <stdio.h>
#include <stdarg.h>

#include <nvector/nvector_serial.h>

/** Generic nvector functions and macros */

N_Vector sunml_alloc_cnvec(size_t content_size, value backlink)
{
    N_Vector nv;

    /* Alloc memory in C heap */
    nv = (N_Vector)calloc(1, sizeof(struct cnvec));
    if (nv == NULL) return NULL;

    nv->ops = (N_Vector_Ops) calloc(1, sizeof(struct _generic_N_Vector_Ops));
    if (nv->ops == NULL) { free(nv); return(NULL); }

    if (content_size != 0) {
	nv->content = (void *) calloc(1, content_size);
	if (nv->content == NULL) { free(nv->ops); free(nv); return(NULL); }
    }

    NVEC_BACKLINK(nv) = backlink;
    caml_register_generational_global_root(&NVEC_BACKLINK(nv));

    return nv;
}

static mlsize_t nvec_rough_size =
    sizeof(struct _generic_N_Vector)
    + sizeof(value)
    + sizeof(struct _generic_N_Vector_Ops)
    + 4 * sizeof(void *);

CAMLprim value sunml_alloc_caml_nvec(N_Vector nv, value finalizer)
{
    CAMLparam1(finalizer);
    CAMLlocal1(r);

    r = caml_alloc_final(1, custom_finalize_default, 1, 30);
    caml_callback2(*caml_named_value("mlfinalise_register"), finalizer, r);
    NVEC_CVAL(r) = nv;

    CAMLreturn(r);
}

void sunml_free_cnvec(N_Vector nv)
{
    caml_remove_generational_global_root(&NVEC_BACKLINK(nv));
    if (nv->content != NULL) free(nv->content);
    free(nv->ops);
    free(nv);
}

void sunml_finalize_caml_nvec (value vnv)
{
    sunml_free_cnvec (NVEC_CVAL (vnv));
}

void sunml_clone_cnvec_ops(N_Vector dst, N_Vector src)
{
#if 500 <= SUNDIALS_LIB_VERSION
    N_VCopyOps(src, dst); // argument order change is intentional

#else
    N_Vector_Ops ops = (N_Vector_Ops) dst->ops;

    ops->nvclone           = src->ops->nvclone;
    ops->nvcloneempty      = src->ops->nvcloneempty;
    ops->nvdestroy         = src->ops->nvdestroy;
#if 270 <= SUNDIALS_LIB_VERSION
    ops->nvgetvectorid	   = src->ops->nvgetvectorid;
#endif
    ops->nvspace           = src->ops->nvspace;
    ops->nvgetarraypointer = src->ops->nvgetarraypointer;
    ops->nvsetarraypointer = src->ops->nvsetarraypointer;
    ops->nvlinearsum       = src->ops->nvlinearsum;
    ops->nvconst           = src->ops->nvconst;  
    ops->nvprod            = src->ops->nvprod;   
    ops->nvdiv             = src->ops->nvdiv;
    ops->nvscale           = src->ops->nvscale; 
    ops->nvabs             = src->ops->nvabs;
    ops->nvinv             = src->ops->nvinv;
    ops->nvaddconst        = src->ops->nvaddconst;
    ops->nvdotprod         = src->ops->nvdotprod;
    ops->nvmaxnorm         = src->ops->nvmaxnorm;
    ops->nvwrmsnormmask    = src->ops->nvwrmsnormmask;
    ops->nvwrmsnorm        = src->ops->nvwrmsnorm;
    ops->nvmin             = src->ops->nvmin;
    ops->nvwl2norm         = src->ops->nvwl2norm;
    ops->nvl1norm          = src->ops->nvl1norm;
    ops->nvcompare         = src->ops->nvcompare;    
    ops->nvinvtest         = src->ops->nvinvtest;
    ops->nvconstrmask      = src->ops->nvconstrmask;
    ops->nvminquotient     = src->ops->nvminquotient;
#if 400 <= SUNDIALS_LIB_VERSION
    /* fused vector operations */
    ops->nvlinearcombination = src->ops->nvlinearcombination;
    ops->nvscaleaddmulti     = src->ops->nvscaleaddmulti;
    ops->nvdotprodmulti      = src->ops->nvdotprodmulti;

    /* vector array operations */
    ops->nvlinearsumvectorarray         = src->ops->nvlinearsumvectorarray;
    ops->nvscalevectorarray             = src->ops->nvscalevectorarray;
    ops->nvconstvectorarray             = src->ops->nvconstvectorarray;
    ops->nvwrmsnormvectorarray          = src->ops->nvwrmsnormvectorarray;
    ops->nvwrmsnormmaskvectorarray      = src->ops->nvwrmsnormmaskvectorarray;
    ops->nvscaleaddmultivectorarray     = src->ops->nvscaleaddmultivectorarray;
    ops->nvlinearcombinationvectorarray = src->ops->nvlinearcombinationvectorarray;
#endif

#endif
}

N_Vector *sunml_nvector_array_alloc(value vtable)
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

void sunml_nvector_array_free(N_Vector *nvarr)
{
    free(nvarr);
}

CAMLprim value sunml_nvec_get_id(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(vr);
#if 290 <= SUNDIALS_LIB_VERSION
    N_Vector_ID id;

    id = N_VGetVectorID(NVEC_VAL(vx));
    switch (id)
    {
	case SUNDIALS_NVEC_SERIAL:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_SERIAL);
	    break;
	case SUNDIALS_NVEC_PARALLEL:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_PARALLEL);
	    break;
	case SUNDIALS_NVEC_OPENMP:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_OPENMP);
	    break;
	case SUNDIALS_NVEC_PTHREADS:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_PTHREADS);
	    break;
	case SUNDIALS_NVEC_PARHYP:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_PARHYP);
	    break;
	case SUNDIALS_NVEC_PETSC:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_PETSC);
	    break;
	case SUNDIALS_NVEC_CUDA:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_CUDA);
	    break;
	case SUNDIALS_NVEC_RAJA:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_RAJA);
	    break;
#if 400 <= SUNDIALS_LIB_VERSION
	case SUNDIALS_NVEC_OPENMPDEV:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_OPENMPDEV);
	    break;
#endif
#if 500 <= SUNDIALS_LIB_VERSION
	case SUNDIALS_NVEC_TRILINOS:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_TRILINOS);
	    break;
	case SUNDIALS_NVEC_MANYVECTOR:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_MANYVECTOR);
	    break;
	case SUNDIALS_NVEC_MPIMANYVECTOR:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_MPIMANYVECTOR);
	    break;
	case SUNDIALS_NVEC_MPIPLUSX:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_MPIPLUSX);
	    break;
#endif
#if 560 <= SUNDIALS_LIB_VERSION
	case SUNDIALS_NVEC_HIP:
	  vr = Val_int(VARIANT_NVECTOR_ID_TAG_HIP);
	  break;
#endif
#if 570 <= SUNDIALS_LIB_VERSION
	case SUNDIALS_NVEC_SYCL:
	  vr = Val_int(VARIANT_NVECTOR_ID_TAG_SYCL);
	  break;
#endif
#if 640 <= SUNDIALS_LIB_VERSION
	case SUNDIALS_NVEC_KOKKOS:
	  vr = Val_int(VARIANT_NVECTOR_ID_TAG_KOKKOS);
	  break;
#endif
	case SUNDIALS_NVEC_CUSTOM:
	default:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_CUSTOM);
	    break;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

CAMLprim value sunml_nvec_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (NVEC, exns);
    CAMLreturn (Val_unit);
}

/** Serial nvectors * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Creation from Sundials/C.  */
/* Adapted from sundials-2.5.0/src/nvec_ser/nvector_serial.c:
   N_VCloneEmpty_Serial */
static N_Vector clone_serial(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);

    N_Vector v;
    N_VectorContent_Serial content;

    if (w == NULL) CAMLreturnT(N_Vector, NULL);
    w_payload = NVEC_BACKLINK(w);
    struct caml_ba_array *w_ba = Caml_ba_array_val(w_payload);

    /* Create vector (we need not copy the data) */
    v_payload = caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim);

    v = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Serial), v_payload);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_Serial) v->content;

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

#if 600 <= SUNDIALS_LIB_VERSION
    v->sunctx = w->sunctx;
#endif

    /* Create content */
    content->length   = NV_LENGTH_S(w);
    content->own_data = 0;
    content->data     = Caml_ba_data_val(v_payload);

#if SUNDIALS_ML_SAFE == 1
    sundials_ml_index i;
    for (i = 0; i < content->length; ++i)
	content->data[i] = 0.0;
#endif

    CAMLreturnT(N_Vector, v);
}

/* Clone an "any" nvector by unwrapping and wrapping the RA payload. */
/* Adapted from sundials-2.5.0/src/nvec_ser/nvector_serial.c:
   N_VCloneEmpty_Serial */
static N_Vector clone_any_serial(N_Vector w)
{
    CAMLparam0();
    CAMLlocal4(v_wrapped, v_payload, w_wrapped, w_payload);

    N_Vector v;
    N_VectorContent_Serial content;

    if (w == NULL) CAMLreturnT(N_Vector, NULL);
    w_wrapped = NVEC_BACKLINK(w);
    w_payload = Field(w_wrapped, 1);

    struct caml_ba_array *w_ba = Caml_ba_array_val(w_payload);

    /* Create vector (we need not copy the data) */
    v_payload = caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim);

    v_wrapped = caml_alloc_tuple(2);
    Store_field(v_wrapped, 0, Field(w_wrapped, 0)); // RA constructor
    Store_field(v_wrapped, 1, v_payload);

    v = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Serial), v_wrapped);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_Serial) v->content;

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

#if 600 <= SUNDIALS_LIB_VERSION
    v->sunctx = w->sunctx;
#endif

    /* Create content */
    content->length   = NV_LENGTH_S(w);
    content->own_data = 0;
    content->data     = Caml_ba_data_val(v_payload);

#if SUNDIALS_ML_SAFE == 1
    sundials_ml_index i;
    for (i = 0; i < content->length; ++i)
	content->data[i] = 0.0;
#endif

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
static N_Vector clone_empty_serial(N_Vector w)
{
    N_Vector v;

    v = N_VCloneEmpty_Serial(w);
    v->ops->nvdestroy = N_VDestroy_Serial;
    v->ops->nvclone = NULL;
    v->ops->nvcloneempty = NULL;

    return v;
}

/* Creation from OCaml.  */
/* Adapted from sundials-2.5.0/src/nvec_ser/nvector_serial.c:
   N_VNewEmpty_Serial */
CAMLprim value sunml_nvec_wrap_serial(value payload,
				      value checkfn, value clonefn,
				      value context)
{
    CAMLparam4(payload, checkfn, clonefn, context);
    CAMLlocal2(vnvec, vnv_cptr);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Serial content;
    long int length = (Caml_ba_array_val(payload))->dim[0];

    /* Create vector */
    nv = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Serial), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    vnv_cptr = sunml_alloc_caml_nvec(nv, *caml_named_value("sunml_finalize_caml_nvec"));
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_Serial) nv->content;

#if 600 <= SUNDIALS_LIB_VERSION
    nv->sunctx = ML_CONTEXT(context);
#endif

    /* Create vector operation structure */
    ops->nvclone           = clone_serial;		    /* ours */
    ops->nvcloneempty      = clone_empty_serial;	    /* ours */
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = sunml_free_cnvec;
#if 270 <= SUNDIALS_LIB_VERSION
    ops->nvgetvectorid	   = N_VGetVectorID_Serial;
#endif

    ops->nvspace           = N_VSpace_Serial;		    /* theirs */
    ops->nvgetarraypointer = N_VGetArrayPointer_Serial;
    ops->nvsetarraypointer = N_VSetArrayPointer_Serial;
    ops->nvlinearsum       = N_VLinearSum_Serial;
    ops->nvconst           = N_VConst_Serial;
    ops->nvprod            = N_VProd_Serial;
    ops->nvdiv             = N_VDiv_Serial;
    ops->nvscale           = N_VScale_Serial;
    ops->nvabs             = N_VAbs_Serial;
    ops->nvinv             = N_VInv_Serial;
    ops->nvaddconst        = N_VAddConst_Serial;
    ops->nvdotprod         = N_VDotProd_Serial;
    ops->nvmaxnorm         = N_VMaxNorm_Serial;
    ops->nvwrmsnormmask    = N_VWrmsNormMask_Serial;
    ops->nvwrmsnorm        = N_VWrmsNorm_Serial;
    ops->nvmin             = N_VMin_Serial;
    ops->nvwl2norm         = N_VWL2Norm_Serial;
    ops->nvl1norm          = N_VL1Norm_Serial;
    ops->nvcompare         = N_VCompare_Serial;
    ops->nvinvtest         = N_VInvTest_Serial;
    ops->nvconstrmask      = N_VConstrMask_Serial;
    ops->nvminquotient     = N_VMinQuotient_Serial;

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
    ops->nvgetlength	    = N_VGetLength_Serial;
#if 650 <= SUNDIALS_LIB_VERSION
    ops->nvgetlocallength   = N_VGetLength_Serial;
#endif
#if 530 <= SUNDIALS_LIB_VERSION
    ops->nvprint	   = N_VPrint_Serial;
    ops->nvprintfile	   = N_VPrintFile_Serial;
#endif
    ops->nvgetcommunicator  = NULL;

    ops->nvdotprodlocal     = N_VDotProd_Serial;
    ops->nvmaxnormlocal     = N_VMaxNorm_Serial;
    ops->nvminlocal         = N_VMin_Serial;
    ops->nvl1normlocal      = N_VL1Norm_Serial;
    ops->nvinvtestlocal     = N_VInvTest_Serial;
    ops->nvconstrmasklocal  = N_VConstrMask_Serial;
    ops->nvminquotientlocal = N_VMinQuotient_Serial;
    ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Serial;
    ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Serial;
#endif
#if 600 <= SUNDIALS_LIB_VERSION
    /* single buffer reduction operations */
    ops->nvdotprodmultilocal = N_VDotProdMulti_Serial;
    ops->nvdotprodmultiallreduce = NULL;
#endif

    /* Create content */
    content->length   = length;
    content->own_data = 0;
    content->data     = Caml_ba_data_val(payload);

    vnvec = NVEC_ALLOC();
    Store_field(vnvec, NVEC_PAYLOAD, payload);
    Store_field(vnvec, NVEC_CPTR, vnv_cptr);
    Store_field(vnvec, NVEC_CHECK, checkfn);
    Store_field(vnvec, NVEC_CLONE, clonefn);
    Store_field(vnvec, NVEC_CONTEXT, context);

    CAMLreturn(vnvec);
}

/* The "any"-version of a serial nvector is created by modifying the
   standard one in two ways:
   1. The payload field is wrapped in the RA constructor.
   2. The nvclone operation is overridden to implement the wrapping operation
      (the current clone_empty_serial does not manipulate the backlink). */
CAMLprim value sunml_nvec_anywrap_serial(value extconstr,
					 value payload,
					 value checkfn, value clonefn,
					 value context)
{
    CAMLparam5(extconstr, payload, checkfn, clonefn, context);
    CAMLlocal2(vnv, vwrapped);
    N_Vector nv;
    N_Vector_Ops ops;

    vnv = sunml_nvec_wrap_serial(payload, checkfn, clonefn, context);
    nv = NVEC_VAL(vnv);
    ops = (N_Vector_Ops) nv->ops;

    ops->nvclone = clone_any_serial;

    vwrapped = caml_alloc_tuple(2);
    Store_field(vwrapped, 0, extconstr);
    Store_field(vwrapped, 1, payload);

    Store_field(vnv, NVEC_PAYLOAD, vwrapped);
    caml_modify_generational_global_root(&NVEC_BACKLINK(nv), vwrapped);

    CAMLreturn(vnv);
}

/** Custom nvectors * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define CNVEC_OP_TABLE(nvec)  ((nvec)->content)

#define GET_OP(nvec, x) (Field((value)CNVEC_OP_TABLE(nvec), x))

#define HAS_OP(ops, x)	     (Field(ops, x) != Val_none)
#define IS_SOME_OP(nvec, x)  (HAS_OP(CNVEC_OP_TABLE(nvec), x))
#define GET_SOME_OP(nvec, x) (Field(Field(CNVEC_OP_TABLE(nvec), x), 0))

static void free_custom_cnvec(N_Vector v)
{
    caml_remove_generational_global_root((value *)&CNVEC_OP_TABLE(v));
    v->content = NULL;
    sunml_free_cnvec(v);
}

CAMLprim value finalize_custom_caml_nvec(value vnv)
{
    CAMLparam1(vnv);
    free_custom_cnvec (NVEC_CVAL(vnv));
    CAMLreturn(Val_unit);
}

#if SUNDIALS_LIB_VERSION >= 270
static N_Vector_ID getvectorid_custom(N_Vector v)
{
      return SUNDIALS_NVEC_CUSTOM;
}
#endif

// Custom operations
static N_Vector callml_vclone(N_Vector w);
static void callml_vspace(N_Vector v, sundials_ml_index *lrw, sundials_ml_index *liw);
#if 500 <= SUNDIALS_LIB_VERSION
static sundials_ml_index callml_vgetlength(N_Vector v);
#endif
#if 650 <= SUNDIALS_LIB_VERSION
static sundials_ml_index callml_vgetlocallength(N_Vector v);
#endif
#if 530 <= SUNDIALS_LIB_VERSION
static void callml_vprint(N_Vector v);
static void callml_vprintfile(N_Vector v, FILE *logfile);
#endif
static void callml_vlinearsum(sunrealtype a, N_Vector x, sunrealtype b,
			      N_Vector y, N_Vector z);
static void callml_vconst(sunrealtype c, N_Vector z);
static void callml_vprod(N_Vector x, N_Vector y, N_Vector z);
static void callml_vdiv(N_Vector x, N_Vector y, N_Vector z);
static void callml_vscale(sunrealtype c, N_Vector x, N_Vector z);
static void callml_vabs(N_Vector x, N_Vector z);
static void callml_vinv(N_Vector x, N_Vector z);
static void callml_vaddconst(N_Vector x, sunrealtype b, N_Vector z);
static sunrealtype callml_vdotprod(N_Vector x, N_Vector y);
static sunrealtype callml_vmaxnorm(N_Vector x);
static sunrealtype callml_vwrmsnorm(N_Vector x, N_Vector w);
static sunrealtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id);
static sunrealtype callml_vmin(N_Vector x);
static sunrealtype callml_vwl2norm(N_Vector x, N_Vector w);
static sunrealtype callml_vl1norm(N_Vector x);
static void callml_vcompare(sunrealtype c, N_Vector x, N_Vector z);
static sunbooleantype callml_vinvtest(N_Vector x, N_Vector z);
static sunbooleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m);
static sunrealtype callml_vminquotient(N_Vector num, N_Vector denom);

#if 500 <= SUNDIALS_LIB_VERSION
static void *callml_vgetcommunicator(N_Vector x);
#endif

/* Custom fused vector operations */
#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vlinearcombination(int nvec, sunrealtype* c,
				     N_Vector* V, N_Vector z);
static int callml_vscaleaddmulti(int nvec, sunrealtype* a,
			            N_Vector x, N_Vector* Y, N_Vector* Z);
static int callml_vdotprodmulti(int nvec, N_Vector x, N_Vector *Y,
				   sunrealtype* dotprods);
#endif

/* Custom vector array operations */
#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vlinearsumvectorarray(int nvec, sunrealtype a, N_Vector* X,
					sunrealtype b, N_Vector* Y, N_Vector* Z);
static int callml_vscalevectorarray(int nvec, sunrealtype* c, N_Vector* X,
				    N_Vector* Z);
static int callml_vconstvectorarray(int nvecs, sunrealtype c, N_Vector* Z);
static int callml_vwrmsnormvectorarray(int nvecs, N_Vector* X,
				       N_Vector* W, sunrealtype* nrm);
static int callml_vwrmsnormmaskvectorarray(int nvecs, N_Vector* X, N_Vector* W,
					   N_Vector id, sunrealtype* nrm);
static int callml_vscaleaddmultivectorarray(int nvec, int nsum, sunrealtype* a,
					    N_Vector* X, N_Vector** Y,
					    N_Vector** Z);
static int callml_vlinearcombinationvectorarray(int nvec, int nsum, sunrealtype* c,
						N_Vector** X, N_Vector* Z);
#endif

#if 600 <= SUNDIALS_LIB_VERSION
static int callml_nvdotprodmultilocal(int nvec, N_Vector x, N_Vector *Y,
				      sunrealtype *d);
static int callml_nvdotprodmultiallreduce(int nvec, N_Vector x, sunrealtype *d);
#endif

/* Custom vector reduction operators */
#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vdotprodlocal(N_Vector x, N_Vector t);
static sunrealtype callml_vmaxnormlocal(N_Vector x);
static sunrealtype callml_vminlocal(N_Vector x);
static sunrealtype callml_vl1normlocal(N_Vector x);
static sunbooleantype callml_vinvtestlocal(N_Vector x, N_Vector z);
static sunbooleantype callml_vconstrmasklocal(N_Vector c, N_Vector x, N_Vector m);
static sunrealtype callml_vminquotientlocal(N_Vector n, N_Vector d);
static sunrealtype callml_vwsqrsumlocal(N_Vector x, N_Vector w);
static sunrealtype callml_vwsqrsummasklocal(N_Vector x, N_Vector w, N_Vector id);
#endif

/* Creation from OCaml. */
CAMLprim value sunml_nvec_wrap_custom(value mlops, value payload,
				      value checkfn, value clonefn,
				      value context)
{
    CAMLparam5(mlops, payload, checkfn, clonefn, context);
    CAMLlocal2(vcnvec, vnv_cptr);

    N_Vector nv;
    N_Vector_Ops ops;

    /* Create vector */
    nv = sunml_alloc_cnvec(0, payload);
    if (nv == NULL) caml_raise_out_of_memory();
    vnv_cptr = sunml_alloc_caml_nvec(nv, *caml_named_value("finalize_custom_caml_nvec"));
    ops = (N_Vector_Ops) nv->ops;

#if 600 <= SUNDIALS_LIB_VERSION
    nv->sunctx = ML_CONTEXT(context);
#endif

    /* Create vector operation structure */
    ops->nvclone           = callml_vclone;
    ops->nvcloneempty      = NULL;
    ops->nvdestroy         = free_custom_cnvec;
#if 270 <= SUNDIALS_LIB_VERSION
    ops->nvgetvectorid	   = getvectorid_custom;
#endif

    ops->nvspace = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVSPACE))
	ops->nvspace = callml_vspace;

    ops->nvlinearsum       = callml_vlinearsum;
    ops->nvconst           = callml_vconst;
    ops->nvprod            = callml_vprod;
    ops->nvdiv             = callml_vdiv;
    ops->nvscale           = callml_vscale;
    ops->nvabs             = callml_vabs;
    ops->nvinv             = callml_vinv;
    ops->nvaddconst        = callml_vaddconst;
    ops->nvmaxnorm         = callml_vmaxnorm;
    ops->nvwrmsnorm        = callml_vwrmsnorm;
    ops->nvmin             = callml_vmin;
    ops->nvdotprod	   = callml_vdotprod;
    ops->nvcompare	   = callml_vcompare;
    ops->nvinvtest	   = callml_vinvtest;

    ops->nvgetarraypointer = NULL;
    ops->nvsetarraypointer = NULL;

    ops->nvwl2norm = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVWL2NORM))
	ops->nvwl2norm = callml_vwl2norm;

    ops->nvl1norm = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVL1NORM))
	ops->nvl1norm = callml_vl1norm;

    ops->nvwrmsnormmask = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVWRMSNORMMASK))
	ops->nvwrmsnormmask = callml_vwrmsnormmask;

    ops->nvconstrmask = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVCONSTRMASK))
	ops->nvconstrmask = callml_vconstrmask;

    ops->nvminquotient = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVMINQUOTIENT))
	ops->nvminquotient = callml_vminquotient;

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
    ops->nvgetlength	    = callml_vgetlength;

#if 650 <= SUNDIALS_LIB_VERSION
    ops->nvgetlocallength   = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVGETLOCALLENGTH))
	ops->nvgetlocallength = callml_vgetlocallength;
#endif

    ops->nvgetcommunicator  = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVGETCOMMUNICATOR))
	ops->nvgetcommunicator = callml_vgetcommunicator;

    ops->nvdotprodlocal     = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVDOTPROD_LOCAL))
	ops->nvdotprodlocal = callml_vdotprodlocal;

    ops->nvmaxnormlocal     = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVMAXNORM_LOCAL))
	ops->nvmaxnormlocal = callml_vmaxnormlocal;

    ops->nvminlocal         = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVMIN_LOCAL))
	ops->nvminlocal = callml_vminlocal;

    ops->nvl1normlocal      = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVL1NORM_LOCAL))
	ops->nvl1normlocal = callml_vl1normlocal;

    ops->nvinvtestlocal     = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVINVTEST_LOCAL))
	ops->nvinvtestlocal = callml_vinvtestlocal;

    ops->nvconstrmasklocal  = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVCONSTRMASK_LOCAL))
	ops->nvconstrmasklocal = callml_vconstrmasklocal;

    ops->nvminquotientlocal = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVMINQUOTIENT_LOCAL))
	ops->nvminquotientlocal = callml_vminquotientlocal;

    ops->nvwsqrsumlocal     = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVWSQRSUM_LOCAL))
	ops->nvwsqrsumlocal = callml_vwsqrsumlocal;

    ops->nvwsqrsummasklocal = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVWSQRSUMMASK_LOCAL))
	ops->nvwsqrsummasklocal = callml_vwsqrsummasklocal;
#endif
#if 600 <= SUNDIALS_LIB_VERSION
    /* single buffer reduction operations */
    ops->nvdotprodmultilocal = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVDOTPRODMULTI_LOCAL))
	ops->nvdotprodmultilocal = callml_nvdotprodmultilocal;

    ops->nvdotprodmultiallreduce = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVDOTPRODMULTI_ALLREDUCE))
	ops->nvdotprodmultiallreduce = callml_nvdotprodmultiallreduce;
#endif
#if 530 <= SUNDIALS_LIB_VERSION
    ops->nvprint	   = NULL;
    ops->nvprintfile	   = NULL;

    if (HAS_OP(mlops, NVECTOR_OPS_NVPRINTFILE)) {
	ops->nvprint = callml_vprint;
	ops->nvprintfile = callml_vprintfile;
    }
#endif

    /* Create content */
    nv->content = (void *)mlops;
    caml_register_generational_global_root((value *)&CNVEC_OP_TABLE(nv));

    vcnvec = NVEC_ALLOC();
    Store_field(vcnvec, NVEC_PAYLOAD, payload);
    Store_field(vcnvec, NVEC_CPTR, vnv_cptr);
    Store_field(vcnvec, NVEC_CHECK, checkfn);
    Store_field(vcnvec, NVEC_CLONE, clonefn);
    Store_field(vcnvec, NVEC_CONTEXT, context);

    CAMLreturn(vcnvec);
}

/* Creation from Sundials/C. */
N_Vector callml_vclone(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);
    N_Vector v;

    if (w == NULL) CAMLreturnT(N_Vector, NULL);
    w_payload = NVEC_BACKLINK(w);

    /* Create vector */

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (GET_OP(w, NVECTOR_OPS_NVCLONE), w_payload);

    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined clone");
	CAMLreturnT (N_Vector, NULL);
    }
    v_payload = r;
    /* Done processing r.  Now it's OK to trigger GC.  */

    v = sunml_alloc_cnvec(0, v_payload);
    if (v == NULL)
	CAMLreturnT (N_Vector, NULL);

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

#if 600 <= SUNDIALS_LIB_VERSION
    v->sunctx = w->sunctx;
#endif

    /* Create content */
    v->content = (void *) CNVEC_OP_TABLE(w);
    caml_register_generational_global_root((value *)&CNVEC_OP_TABLE(v));

    CAMLreturnT(N_Vector, v);
}

static void callml_vspace(N_Vector v, sundials_ml_index *lrw, sundials_ml_index *liw)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(v, NVECTOR_OPS_NVSPACE);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(v));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined space");
	fputs ("Sundials/ML has no sensible value to return to Sundials, "
	       "and incorrect values risk memory corruption.  Abort.", stderr);
	fflush (stderr);
	abort ();
    } else {
	*lrw = Index_val(Field(r, 0));
	*liw = Index_val(Field(r, 1)) + (nvec_rough_size / sizeof(int));
    }

    CAMLreturn0;
}

#if 500 <= SUNDIALS_LIB_VERSION
static sundials_ml_index callml_vgetlength(N_Vector v)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(v, NVECTOR_OPS_NVGETLENGTH);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(v));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined getlength");
	fputs ("Sundials/ML has no sensible value to return to Sundials, "
	       "and incorrect values risk memory corruption.  Abort.", stderr);
	fflush (stderr);
	abort ();
    }

    CAMLreturnT(sundials_ml_index, Int_val(r));
}
#endif

#if 650 <= SUNDIALS_LIB_VERSION
static sundials_ml_index callml_vgetlocallength(N_Vector v)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(v, NVECTOR_OPS_NVGETLOCALLENGTH);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(v));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined getlocallength");
	fputs ("Sundials/ML has no sensible value to return to Sundials, "
	       "and incorrect values risk memory corruption.  Abort.", stderr);
	fflush (stderr);
	abort ();
    }

    CAMLreturnT(sundials_ml_index, Int_val(r));
}
#endif

#if 530 <= SUNDIALS_LIB_VERSION
static void callml_vprint(N_Vector v)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(v, NVECTOR_OPS_NVPRINTFILE);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(v), Val_none);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined printfile");
	fputs ("Sundials/ML has no sensible value to return to Sundials, "
	       "and incorrect values risk memory corruption.  Abort.", stderr);
	fflush (stderr);
	abort ();
    }

    CAMLreturn0;
}

static void callml_vprintfile(N_Vector v, FILE* logfile)
{
    CAMLparam0();
    CAMLlocal2(mlop, vologfile);
    mlop = GET_SOME_OP(v, NVECTOR_OPS_NVPRINTFILE);
    vologfile = sunml_sundials_wrap_file(logfile);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(v), vologfile);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined printfile");
	fputs ("Sundials/ML has no sensible value to return to Sundials, "
	       "and incorrect values risk memory corruption.  Abort.", stderr);
	fflush (stderr);
	abort ();
    }

    CAMLreturn0;
}
#endif

static void callml_vlinearsum(sunrealtype a, N_Vector x, sunrealtype b,
			      N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    CAMLlocalN(args, 5);

    mlop = GET_OP(x, NVECTOR_OPS_NVLINEARSUM);

    args[0] = caml_copy_double(a);
    args[1] = NVEC_BACKLINK(x);
    args[2] = caml_copy_double(b);
    args[3] = NVEC_BACKLINK(y);
    args[4] = NVEC_BACKLINK(z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 5, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined linearsum");
    }

    CAMLreturn0;
}

static void callml_vconst(sunrealtype c, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vc);

    vc = caml_copy_double (c);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(GET_OP(z, NVECTOR_OPS_NVCONST),
				 vc, NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined const");

    CAMLreturn0;
}

static void callml_vprod(N_Vector x, N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVPROD);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(x),
				  NVEC_BACKLINK(y), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined prod");

    CAMLreturn0;
}

static void callml_vdiv(N_Vector x, N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVDIV);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(mlop, NVEC_BACKLINK(x),
				 NVEC_BACKLINK(y), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined div");

    CAMLreturn0;
}

static void callml_vscale(sunrealtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vc);

    vc = caml_copy_double(c);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(GET_OP(x, NVECTOR_OPS_NVSCALE), vc,
				 NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined scale");

    CAMLreturn0;
}

static void callml_vabs(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVABS);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined abs");

    CAMLreturn0;
}

static void callml_vinv(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVINV);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined inv");

    CAMLreturn0;
}

static void callml_vaddconst(N_Vector x, sunrealtype b, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vb);

    vb = caml_copy_double(b);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (GET_OP(x, NVECTOR_OPS_NVADDCONST),
				  NVEC_BACKLINK(x), vb, NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined addconst");

    CAMLreturn0;
}

static sunrealtype callml_vdotprod(N_Vector x, N_Vector y)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVDOTPROD);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(y));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined dotprod");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

static sunrealtype callml_vmaxnorm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVMAXNORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined maxnorm");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

static sunrealtype callml_vwrmsnorm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVWRMSNORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(w));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined wrmsnorm");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

static sunrealtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMMASK);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(x),
				  NVEC_BACKLINK(w), NVEC_BACKLINK(id));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined wrmsnormmask");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

static sunrealtype callml_vmin(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVMIN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined min");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

static sunrealtype callml_vwl2norm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWL2NORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(w));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined wl2norm");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

static sunrealtype callml_vl1norm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVL1NORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined l1norm");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

static void callml_vcompare(sunrealtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vc);

    vc = caml_copy_double(c);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (GET_OP(x, NVECTOR_OPS_NVCOMPARE), vc,
				  NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined compare");

    CAMLreturn0;
}

static sunbooleantype callml_vinvtest(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVINVTEST);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined invtest");
	CAMLreturnT(sunbooleantype, 0);
    }

    CAMLreturnT(sunbooleantype, Bool_val(r));
}

static sunbooleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVCONSTRMASK);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(c),
				  NVEC_BACKLINK(x), NVEC_BACKLINK(m));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined constrmask");
	CAMLreturnT(sunbooleantype, 0);
    }

    CAMLreturnT(sunbooleantype, Bool_val(r));
}

static sunrealtype callml_vminquotient(N_Vector num, N_Vector denom)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(num, NVECTOR_OPS_NVMINQUOTIENT);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(num),
				  NVEC_BACKLINK(denom));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined minquotient");
	CAMLreturnT(sunrealtype, nan(""));
    }

    CAMLreturnT(sunrealtype, Double_val(r));
}

void sunml_nvectors_into_array(int n, value vy, N_Vector *y)
{
    int i;
    for (i = 0; i < n; ++i) {
	Store_field(vy, i, NVEC_BACKLINK(y[i]));
    }
}

// NB: Normally, we should worry about relinquishing the elements of vy
// after we are finished using them (so as not to block the GC), but we
// instead make the assumption that these elements come from 'within'
// Sundials and thus that they would anyway not be GC-ed.
#if 400 <= SUNDIALS_LIB_VERSION
value sunml_wrap_to_nvector_table(int n, N_Vector *y)
{
    CAMLparam0();
    CAMLlocal1(vy);
    int i;

    vy = caml_alloc_tuple(n);

    for (i = 0; i < n; ++i) {
	Store_field(vy, i, NVEC_BACKLINK(y[i]));
    }

    CAMLreturn (vy);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
value sunml_wrap_to_nvector_tables(int n1, int n2, N_Vector **yy)
{
    CAMLparam0();
    CAMLlocal1(vyy);
    int i;

    vyy = caml_alloc_tuple(n1);

    for (i = 0; i < n1; ++i) {
	Store_field(vyy, i, sunml_wrap_to_nvector_table(n2, yy[i]));
    }

    CAMLreturn (vyy);
}
#endif

/* fused vector operations */
#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vlinearcombination(int nvec, sunrealtype* c,
				     N_Vector* V, N_Vector z)
{
    CAMLparam0();
    CAMLlocal3(mlop, vc, vv);
    intnat n = nvec;

    mlop = GET_SOME_OP(z, NVECTOR_OPS_NVLINEARCOMBINATION);
    vv = sunml_wrap_to_nvector_table(nvec, V);
    vc = caml_ba_alloc(BIGARRAY_FLOAT, 1, c, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vc, vv, NVEC_BACKLINK(z));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined linearcombination");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vscaleaddmulti(int nvec, sunrealtype* a,
			         N_Vector x, N_Vector* Y, N_Vector* Z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    CAMLlocalN(args, 4);
    intnat n = nvec;

    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVSCALEADDMULTI);
    args[0] = caml_ba_alloc(BIGARRAY_FLOAT, 1, a, &n);
    args[1] = NVEC_BACKLINK(x);
    args[2] = sunml_wrap_to_nvector_table(nvec, Y);
    args[3] = sunml_wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 4, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined scaleaddmulti");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vdotprodmulti(int nvec, N_Vector x, N_Vector *Y,
				sunrealtype* dotprods)
{
    CAMLparam0();
    CAMLlocal3(mlop, vy, vdotprods);
    intnat n = nvec;

    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVDOTPRODMULTI);
    vy = sunml_wrap_to_nvector_table(nvec, Y);
    vdotprods = caml_ba_alloc(BIGARRAY_FLOAT, 1, dotprods, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(x), vy, vdotprods);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined dotprodmulti");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

/* vector array operations */
#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vlinearsumvectorarray(int nvec, sunrealtype a, N_Vector* X,
					   sunrealtype b, N_Vector* Y,
                                           N_Vector* Z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    CAMLlocalN(args, 5);

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVLINEARSUMVECTORARRAY);
    args[0] = caml_copy_double(a);
    args[1] = sunml_wrap_to_nvector_table(nvec, X);
    args[2] = caml_copy_double(b);
    args[3] = sunml_wrap_to_nvector_table(nvec, Y);
    args[4] = sunml_wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 5, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined linearsumvectorarray");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vscalevectorarray(int nvec, sunrealtype* c, N_Vector* X, N_Vector* Z)
{
    CAMLparam0();
    CAMLlocal4(mlop, vc, vx, vz);
    intnat n = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVSCALEVECTORARRAY);
    vc = caml_ba_alloc(BIGARRAY_FLOAT, 1, c, &n);
    vx = sunml_wrap_to_nvector_table(nvec, X);
    vz = sunml_wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vc, vx, vz);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined scalevectorarray");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vconstvectorarray(int nvec, sunrealtype c, N_Vector* Z)
{
    CAMLparam0();
    CAMLlocal3(mlop, vc, vz);

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(Z[0], NVECTOR_OPS_NVCONSTVECTORARRAY);
    vc = caml_copy_double(c);
    vz = sunml_wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, vc, vz);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined constvectorarray");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vwrmsnormvectorarray(int nvec, N_Vector* X,
				       N_Vector* W, sunrealtype* nrm)
{
    CAMLparam0();
    CAMLlocal4(mlop, vx, vw, vnrm);
    intnat n = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVWRMSNORMVECTORARRAY);
    vx = sunml_wrap_to_nvector_table(nvec, X);
    vw = sunml_wrap_to_nvector_table(nvec, W);
    vnrm = caml_ba_alloc(BIGARRAY_FLOAT, 1, nrm, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vx, vw, vnrm);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined wrmsnormvectorarray");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vwrmsnormmaskvectorarray(int nvec, N_Vector* X, N_Vector* W,
					   N_Vector id, sunrealtype* nrm)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    CAMLlocalN(args, 4);
    intnat n = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVWRMSNORMMASKVECTORARRAY);
    args[0] = sunml_wrap_to_nvector_table(nvec, X);
    args[1] = sunml_wrap_to_nvector_table(nvec, W);
    args[2] = NVEC_BACKLINK(id);
    args[3] = caml_ba_alloc(BIGARRAY_FLOAT, 1, nrm, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 4, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
				    "user-defined wrmsnormmaskvectorarray");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vscaleaddmultivectorarray(int nvec, int nsum, sunrealtype* a,
					       N_Vector* X, N_Vector** Y,
					       N_Vector** Z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    CAMLlocalN(args, 4);
    intnat n = nsum;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVSCALEADDMULTIVECTORARRAY);
    args[0] = caml_ba_alloc(BIGARRAY_FLOAT, 1, a, &n);
    args[1] = sunml_wrap_to_nvector_table(nvec, X);
    args[2] = sunml_wrap_to_nvector_tables(nsum, nvec, Y);
    args[3] = sunml_wrap_to_nvector_tables(nsum, nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 4, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
				    "user-defined scaleaddmultivectorarray");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int callml_vlinearcombinationvectorarray(int nvec, int nsum, sunrealtype* c,
						N_Vector** X, N_Vector* Z)
{
    CAMLparam0();
    CAMLlocal4(mlop, vc, vxx, vz);
    intnat n2 = nsum;

    if ((nvec <= 0) || (nsum <= 0)) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(Z[0], NVECTOR_OPS_NVLINEARCOMBINATIONVECTORARRAY);
    vc = caml_ba_alloc(BIGARRAY_FLOAT, 1, c, &n2);
    vxx = sunml_wrap_to_nvector_tables(nsum, nvec, X);
    vz = sunml_wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vc, vxx, vz);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
			    "user-defined linearcombinationvectorarray");
	CAMLreturnT(int, 0);
    }
    CAMLreturnT(int, 1);
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
/* Must correspond with camlmpi.h after replacing MPI_Comm* by void* */
#define Comm_val_addr(comm) ((void *) &Field(comm, 1))

static void *callml_vgetcommunicator(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVGETCOMMUNICATOR);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined getcommunicator");

    CAMLreturnT(void *, Comm_val_addr(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vdotprodlocal(N_Vector x, N_Vector t)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVDOTPROD_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(t));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined dotprodlocal");

    CAMLreturnT(sunrealtype, Double_val(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vmaxnormlocal(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVMAXNORM_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined maxnormlocal");

    CAMLreturnT(sunrealtype, Double_val(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vminlocal(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVMIN_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined minlocal");

    CAMLreturnT(sunrealtype, Double_val(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vl1normlocal(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVL1NORM_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined l1normlocal");

    CAMLreturnT(sunrealtype, Double_val(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunbooleantype callml_vinvtestlocal(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVINVTEST_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined invtestlocal");

    CAMLreturnT(sunbooleantype, Bool_val(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunbooleantype callml_vconstrmasklocal(N_Vector c, N_Vector x, N_Vector m)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVCONSTRMASK_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(c),
					NVEC_BACKLINK(x),
					NVEC_BACKLINK(m));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined constrmasklocal");

    CAMLreturnT(sunbooleantype, Bool_val(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vminquotientlocal(N_Vector n, N_Vector d)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(n, NVECTOR_OPS_NVMINQUOTIENT_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(n), NVEC_BACKLINK(d));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined minquotientlocal");

    CAMLreturnT(sunrealtype, caml_copy_double(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vwsqrsumlocal(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWSQRSUM_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(w));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined wsqrsumlocal");

    CAMLreturnT(sunrealtype, Bool_val(r));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static sunrealtype callml_vwsqrsummasklocal(N_Vector x, N_Vector w, N_Vector id)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWSQRSUMMASK_LOCAL);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(x),
					NVEC_BACKLINK(w),
					NVEC_BACKLINK(id));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined wsqrsummasklocal");

    CAMLreturnT(sunrealtype, Bool_val(r));
}
#endif

#if 600 <= SUNDIALS_LIB_VERSION
static int callml_nvdotprodmultilocal(int nvec, N_Vector x, N_Vector *Y,
				      sunrealtype *d)
{
    CAMLparam0();
    CAMLlocal3(mlop, vy, vd);
    intnat nv = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVDOTPRODMULTI_LOCAL);
    vy = sunml_wrap_to_nvector_table(nvec, Y);
    vd = caml_ba_alloc(BIGARRAY_FLOAT, 1, d, &nv);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(x), vy, vd);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined dotprodmultilocal");
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, 1);
}

static int callml_nvdotprodmultiallreduce(int nvec, N_Vector x, sunrealtype *d)
{
    CAMLparam0();
    CAMLlocal2(mlop, vd);
    intnat nv = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVDOTPRODMULTI_ALLREDUCE);
    vd = caml_ba_alloc(BIGARRAY_FLOAT, 1, d, &nv);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), vd);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined dotprodallreduce");
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, 1);
}
#endif

/** Interface to underlying serial nvector functions */

CAMLprim value sunml_nvec_ser_linearsum(value va, value vx, value vb,
					   value vy, value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_VLinearSum_Serial(Double_val(va), NVEC_VAL(vx),
			Double_val(vb), NVEC_VAL(vy),
			NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_const(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst_Serial(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_prod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VProd_Serial(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_div(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VDiv_Serial(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_scale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VScale_Serial(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_abs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VAbs_Serial(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_inv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VInv_Serial(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_addconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_VAddConst_Serial(NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_dotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    sunrealtype r = N_VDotProd_Serial(NVEC_VAL(vx), NVEC_VAL(vy));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_maxnorm(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = N_VMaxNorm_Serial(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_wrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    sunrealtype r = N_VWrmsNorm_Serial(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_wrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    sunrealtype r = N_VWrmsNormMask_Serial(NVEC_VAL(vx), NVEC_VAL(vw),
					NVEC_VAL(vid));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_min(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = N_VMin_Serial(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_wl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    sunrealtype r = N_VWL2Norm_Serial(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_l1norm(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = N_VL1Norm_Serial(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_compare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VCompare_Serial(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_invtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    sunbooleantype r = N_VInvTest_Serial(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_ser_constrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    sunbooleantype r = N_VConstrMask_Serial(NVEC_VAL(vc),
					 NVEC_VAL(vx),
					 NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_ser_minquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    sunrealtype r = N_VMinQuotient_Serial(NVEC_VAL(vnum), NVEC_VAL(vdenom));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_space(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
    sundials_ml_index lrw, liw;

    N_VSpace_Serial(NVEC_VAL(vx), &lrw, &liw);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_index(lrw));
    Store_field(r, 1, Val_index(liw));

    CAMLreturn(r);
}

CAMLprim value sunml_nvec_ser_getlength(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
#if 500 <= SUNDIALS_LIB_VERSION
    r = Val_int(N_VGetLength_Serial(NVEC_VAL(vx)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_nvec_ser_print_file(value vx, value volog)
{
    CAMLparam2(vx, volog);
#if 270 <= SUNDIALS_LIB_VERSION
    if (volog == Val_none) {
	N_VPrint_Serial(NVEC_VAL(vx));
    } else {
#if 310 <= SUNDIALS_LIB_VERSION
	N_VPrintFile_Serial(NVEC_VAL(vx), ML_CFILE(Some_val(volog)));
#else
	caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* fused vector operations */

#if 400 <= SUNDIALS_LIB_VERSION
/* arrays_of_nvectors(rr, 3, vax, vay, vaz) creates three contiguous
   arrays of N_Vectors (a single call to calloc), checking that they
   all have the same length, and that all the nvectors they contain
   are of the same length, and storing the corresponding pointers
   in rr[0] for vax, rr[1] for vay, and rr[2] for vaz. Returns the
   length of vax/vay/vaz on success, or 0 on fail. After a successful
   call (only), the caller must free the allocated memory with free(*r). */
int sunml_arrays_of_nvectors(N_Vector *r[], int n, ...)
{
    va_list valist;
    int i, j;
    int nvecs;
    N_Vector *px = NULL;
    value vax;

    va_start(valist, n);
    vax = va_arg(valist, value);
    nvecs = Wosize_val(vax);

    px = calloc(n * nvecs, sizeof(N_Vector));
    if (px == NULL) goto error;

    j = 0;
    while (1) {
	r[j] = px;

	for (i=0; i < nvecs; ++i, ++px) {
	    *px = NVEC_VAL(Field(vax, i));
	}
	++j;
	if (j == n) break;

	vax = va_arg(valist, value);
    }

    va_end(valist);
    return nvecs;

error:
    if (px != NULL) free(r[0]);
    for (j=0; j < n; ++j) r[0] = NULL;
    va_end(valist);
    return 0;
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
/* sunml_arrays_of_nvectors2(nr, nc, vv, 3, vaax, vaay, vaaz) creates three
   contiguous two-dimensional arrays of N_Vectors (a single call to malloc),
   checking that they all have the same dimensions, and that all the nvectors
   they contain are of the same length, and storing the corresponding pointers
   in vv[0] for vaax, vv[1] for vaay, and vv[2] for vaaz. It sets nrows and
   ncols to the first and second dimensions of the vv[i]'s on success, or both
   to zero on fail. After a successful call (only), the caller must free the
   allocated memory with free(*r). */
void sunml_arrays_of_nvectors2(int* nrows, int *ncols, N_Vector **vv[],
			       int n, ...)
{
    va_list valist;
    int j;
    N_Vector **px = NULL;
    N_Vector *pd;
    value vaax, vax;
    int nr, nc, r, c;

    va_start(valist, n);
    vaax = va_arg(valist, value);
    *nrows = nr = Wosize_val(vaax);
    if (nr < 1) goto error;

    vax = Field(vaax, 0);
    *ncols = nc = Wosize_val(vax);
    if (nc < 1) goto error;

    px = malloc(n * ((nr * sizeof(N_Vector *))
		     + (nr * nc * sizeof(N_Vector))));
    if (px == NULL) goto error;

    // the row pointers start at px; the column data starts at pd
    pd = (N_Vector *)(px + (n * nr));

    j = 0;
    while (1) {
	vv[j] = px;

	for (r=0; r < nr; ++r, ++px) {
	    *px = pd;
	    vax = Field(vaax, r);
	    for (c=0; c < nc; ++c, ++pd) {
		*pd = NVEC_VAL(Field(vax, c));
	    }
	}

	++j;
	if (j == n) break;

	vaax = va_arg(valist, value);
    }

    va_end(valist);
    return;

error:
    if (px != NULL) free(vv[0]);
    for (j=0; j < n; ++j) vv[j] = NULL;
    *nrows = *ncols = 0;
    va_end(valist);
    return;
}
#endif

CAMLprim value sunml_nvec_ser_linearcombination(value vac, value vax,
						   value vz)
{
    CAMLparam3(vac, vax, vz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector z = NVEC_VAL(vz);
    N_Vector *ax;
    int nvec = sunml_arrays_of_nvectors(&ax, 1, vax);
    if (!nvec) caml_raise_out_of_memory();

    N_VLinearCombination_Serial(nvec, ac, ax, z);
    free(ax);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_scaleaddmulti(value vac, value vx, value vay,
					       value vaz)
{
    CAMLparam4(vac, vx, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vay, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VScaleAddMulti_Serial(nvec, ac, x, a[0], a[1]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_dotprodmulti(value vx, value vay, value vad)
{
    CAMLparam3(vx, vay, vad);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ad = REAL_ARRAY(vad);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *ay;
    int nvec = sunml_arrays_of_nvectors(&ay, 1, vay);
    if (!nvec) caml_raise_out_of_memory();

    N_VDotProdMulti_Serial(nvec, x, ay, ad);
    free(ay);
#endif
    CAMLreturn(Val_unit);
}

/* vector array operations */

CAMLprim value sunml_nvec_ser_linearsumvectorarray(value va, value vax,
						      value vb, value vay,
						      value vaz)
{
    CAMLparam5(va, vax, vb, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *a[3];
    int nvec = sunml_arrays_of_nvectors(a, 3, vax, vay, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VLinearSumVectorArray_Serial(nvec, Double_val(va), a[0],
					 Double_val(vb), a[1], a[2]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_scalevectorarray(value vac, value vax,
						  value vaz)
{
    CAMLparam3(vac, vax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VScaleVectorArray_Serial(nvec, ac, a[0], a[1]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_constvectorarray(value vc, value vaz)
{
    CAMLparam2(vc, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *az;
    int nvec = sunml_arrays_of_nvectors(&az, 1, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VConstVectorArray_Serial(nvec, Double_val(vc), az);
    free(az);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_wrmsnormvectorarray(value vax, value vaw,
						     value van)
{
    CAMLparam3(vax, vaw, van);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *an = REAL_ARRAY(van);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);
    if (!nvec) caml_raise_out_of_memory();

    N_VWrmsNormVectorArray_Serial(nvec, a[0], a[1], an);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_wrmsnormmaskvectorarray(value vax,
						value vaw, value vi, value van)
{
    CAMLparam4(vax, vaw, vi, van);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *an = REAL_ARRAY(van);
    N_Vector i = NVEC_VAL(vi);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);
    if (!nvec) caml_raise_out_of_memory();

    N_VWrmsNormMaskVectorArray_Serial(nvec, a[0], a[1], i, an);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_scaleaddmultivectorarray(value vaa,
					    value vax, value vaay, value vaaz)
{
    CAMLparam4(vaa, vax, vaay, vaaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *aa = REAL_ARRAY(vaa);
    N_Vector *ax = NULL;
    int nvec = sunml_arrays_of_nvectors(&ax, 1, vax);
    if (!nvec) caml_raise_out_of_memory();
    N_Vector **ayz[2] = { NULL };
    int nvec2, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec2, ayz, 2, vaay, vaaz);
    if (!nsum) caml_raise_out_of_memory();

    N_VScaleAddMultiVectorArray_Serial(nvec, nsum, aa, ax, ayz[0], ayz[1]);
    free(ax);
    free(*ayz);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_linearcombinationvectorarray(value vac,
							value vaax, value vaz)
{
    CAMLparam3(vac, vaax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector *az;
    int nvecz __attribute__((unused))
	= sunml_arrays_of_nvectors(&az, 1, vaz);
    N_Vector **aax;
    int nvec, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec, &aax, 1, vaax);
    if (!nsum) caml_raise_out_of_memory();

    N_VLinearCombinationVectorArray_Serial(nvec, nsum, ac, aax, az);
    free(aax);
#endif
    CAMLreturn(Val_unit);
}

/** Reduce operations for serial nvectors */

CAMLprim value sunml_nvec_ser_wsqrsumlocal(value vx, value vw)
{
    CAMLparam2(vx, vw);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

    r = N_VWSqrSumLocal_Serial(x, w);
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_wsqrsummasklocal(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

    r = N_VWSqrSumMaskLocal_Serial(x, w, id);
#endif
    CAMLreturn(caml_copy_double(r));
}

/** Selectively activate fused and array operations for serial nvectors */

CAMLprim value sunml_nvec_ser_enablefusedops(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableFusedOps_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablelinearcombination(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombination_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablescaleaddmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMulti_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enabledotprodmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableDotProdMulti_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablelinearsumvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearSumVectorArray_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablescalevectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleVectorArray_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enableconstvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableConstVectorArray_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablewrmsnormvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormVectorArray_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablewrmsnormmaskvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormMaskVectorArray_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablescaleaddmultivectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMultiVectorArray_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_enablelinearcombinationvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombinationVectorArray_Serial(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/** Selectively activate fused and array operations for custom nvectors */

CAMLprim value sunml_nvec_custom_enablefusedops(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (Bool_val(vv)) {
	r =    IS_SOME_OP(x, NVECTOR_OPS_NVLINEARCOMBINATION)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVSCALEADDMULTI)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVDOTPRODMULTI)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVLINEARSUMVECTORARRAY)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVSCALEVECTORARRAY)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVCONSTVECTORARRAY)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMVECTORARRAY)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMMASKVECTORARRAY)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVSCALEADDMULTIVECTORARRAY)
	    && IS_SOME_OP(x, NVECTOR_OPS_NVLINEARCOMBINATIONVECTORARRAY);

	if (r) {
	    x->ops->nvlinearcombination = callml_vlinearcombination;
	    x->ops->nvscaleaddmulti = callml_vscaleaddmulti;
	    x->ops->nvdotprodmulti = callml_vdotprodmulti;
	    x->ops->nvlinearsumvectorarray = callml_vlinearsumvectorarray;
	    x->ops->nvscalevectorarray = callml_vscalevectorarray;
	    x->ops->nvconstvectorarray = callml_vconstvectorarray;
	    x->ops->nvwrmsnormvectorarray = callml_vwrmsnormvectorarray;
	    x->ops->nvwrmsnormmaskvectorarray = callml_vwrmsnormmaskvectorarray;
	    x->ops->nvscaleaddmultivectorarray =
		callml_vscaleaddmultivectorarray;
	    x->ops->nvlinearcombinationvectorarray =
		callml_vlinearcombinationvectorarray;
	}
    } else {
	x->ops->nvlinearcombination = NULL;
	x->ops->nvscaleaddmulti     = NULL;
	x->ops->nvdotprodmulti      = NULL;
	x->ops->nvlinearsumvectorarray         = NULL;
	x->ops->nvscalevectorarray             = NULL;
	x->ops->nvconstvectorarray             = NULL;
	x->ops->nvwrmsnormvectorarray          = NULL;
	x->ops->nvwrmsnormmaskvectorarray      = NULL;
	x->ops->nvscaleaddmultivectorarray     = NULL;
	x->ops->nvlinearcombinationvectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablelinearcombination(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVLINEARCOMBINATION)) {
	    x->ops->nvlinearcombination = callml_vlinearcombination;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvlinearcombination = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablescaleaddmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVSCALEADDMULTI)) {
	    x->ops->nvscaleaddmulti = callml_vscaleaddmulti;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvscaleaddmulti = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enabledotprodmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVDOTPRODMULTI)) {
	    x->ops->nvdotprodmulti = callml_vdotprodmulti;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvdotprodmulti = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablelinearsumvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVLINEARSUMVECTORARRAY)) {
	    x->ops->nvlinearsumvectorarray = callml_vlinearsumvectorarray;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvlinearsumvectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablescalevectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVSCALEVECTORARRAY)) {
	    x->ops->nvscalevectorarray = callml_vscalevectorarray;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvscalevectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enableconstvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVCONSTVECTORARRAY)) {
	    x->ops->nvconstvectorarray = callml_vconstvectorarray;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvconstvectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablewrmsnormvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMVECTORARRAY)) {
	    x->ops->nvwrmsnormvectorarray = callml_vwrmsnormvectorarray;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvwrmsnormvectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablewrmsnormmaskvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMMASKVECTORARRAY)) {
	    x->ops->nvwrmsnormmaskvectorarray = callml_vwrmsnormmaskvectorarray;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvwrmsnormmaskvectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablescaleaddmultivectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVSCALEADDMULTIVECTORARRAY)) {
	    x->ops->nvscaleaddmultivectorarray =
		callml_vscaleaddmultivectorarray;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvscaleaddmultivectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

CAMLprim value sunml_nvec_custom_enablelinearcombinationvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
    sunbooleantype r = 1;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    if (Bool_val(vv)) {
	if (IS_SOME_OP(x, NVECTOR_OPS_NVLINEARCOMBINATIONVECTORARRAY)) {
	    x->ops->nvlinearcombinationvectorarray =
		callml_vlinearcombinationvectorarray;
	} else {
	    r = 0;
	}
    } else {
	x->ops->nvlinearcombinationvectorarray = NULL;
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_bool(r));
}

/** Test for nvector features */

CAMLprim value sunml_nvec_has_linearcombination(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvlinearcombination != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_scaleaddmulti(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvscaleaddmulti != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_dotprodmulti(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvdotprodmulti != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_linearsumvectorarray(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvlinearsumvectorarray != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_scalevectorarray(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvscalevectorarray != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_constvectorarray(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvconstvectorarray != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_wrmsnormvectorarray(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvwrmsnormvectorarray != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_wrmsnormmaskvectorarray(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvwrmsnormmaskvectorarray != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_scaleaddmultivectorarray(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvscaleaddmultivectorarray != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_linearcombinationvectorarray(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvlinearcombinationvectorarray != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_dotprodlocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvdotprodlocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_maxnormlocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvmaxnormlocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_minlocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvminlocal != NULL);
#endif

    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_l1normlocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvl1normlocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_invtestlocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvinvtestlocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_constrmasklocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvconstrmasklocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_minquotientlocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvminquotientlocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_wsqrsumlocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvwsqrsumlocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_wsqrsummasklocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvwsqrsummasklocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_dotprodmultilocal(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 600 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvdotprodmultilocal != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_has_dotprodmultiallreduce(value vx)
{
    CAMLparam1(vx);
    sunbooleantype r = 0;
#if 600 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    r = (x->ops->nvdotprodmultiallreduce != NULL);
#endif
    CAMLreturn(Val_bool(r));
}

/** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
 * Interface to underlying generic nvector functions */

CAMLprim value sunml_nvec_linearsum(value va, value vx, value vb,
				    value vy, value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_VLinearSum(Double_val(va), NVEC_VAL(vx), Double_val(vb),
		 NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_const(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_prod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VProd(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_div(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VDiv(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_scale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VScale(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_abs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VAbs(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_inv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VInv(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_addconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_VAddConst(NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_dotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    sunrealtype r = N_VDotProd(NVEC_VAL(vx), NVEC_VAL(vy));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_maxnorm(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = N_VMaxNorm(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_wrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    sunrealtype r = N_VWrmsNorm(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_wrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvwrmsnormmask == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    sunrealtype r = N_VWrmsNormMask(x, NVEC_VAL(vw), NVEC_VAL(vid));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_min(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = N_VMin(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_wl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvwl2norm == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    sunrealtype r = N_VWL2Norm(x, NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_l1norm(value vx)
{
    CAMLparam1(vx);
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvl1norm == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    sunrealtype r = N_VL1Norm(x);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_compare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VCompare(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_invtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    sunbooleantype r = N_VInvTest(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_constrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    N_Vector c = NVEC_VAL(vc);

    if (c->ops->nvconstrmask == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    sunbooleantype r = N_VConstrMask(c, NVEC_VAL(vx), NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_minquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    N_Vector num = NVEC_VAL(vnum);

    if (num->ops->nvminquotient == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    sunrealtype r = N_VMinQuotient(num, NVEC_VAL(vdenom));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_space(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
    sundials_ml_index lrw, liw;
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvspace == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));
    N_VSpace(x, &lrw, &liw);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_index(lrw));
    Store_field(r, 1, Val_index(liw));

    CAMLreturn(r);
}

CAMLprim value sunml_nvec_getlength(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
#if 500 <= SUNDIALS_LIB_VERSION
    r = Val_int(N_VGetLength(NVEC_VAL(vx)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_nvec_getlocallength(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
#if 650 <= SUNDIALS_LIB_VERSION
    r = Val_int(N_VGetLocalLength(NVEC_VAL(vx)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_nvec_print_file(value vx, value volog)
{
    CAMLparam2(vx, volog);
#if 530 <= SUNDIALS_LIB_VERSION
    if (volog == Val_none) {
	N_VPrint(NVEC_VAL(vx));
    } else {
	N_VPrintFile(NVEC_VAL(vx), ML_CFILE(Some_val(volog)));
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

/* fused vector operations */

CAMLprim value sunml_nvec_linearcombination(value vac, value vax, value vz)
{
    CAMLparam3(vac, vax, vz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector z = NVEC_VAL(vz);
    N_Vector *ax;
    int nvec = sunml_arrays_of_nvectors(&ax, 1, vax);
    if (!nvec) caml_raise_out_of_memory();

    N_VLinearCombination(nvec, ac, ax, z);
    free(ax);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_scaleaddmulti(value vac, value vx, value vay,
				        value vaz)
{
    CAMLparam4(vac, vx, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vay, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VScaleAddMulti(nvec, ac, x, a[0], a[1]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_dotprodmulti(value vx, value vay, value vad)
{
    CAMLparam3(vx, vay, vad);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ad = REAL_ARRAY(vad);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *ay;
    int nvec = sunml_arrays_of_nvectors(&ay, 1, vay);
    if (!nvec) caml_raise_out_of_memory();

    N_VDotProdMulti(nvec, x, ay, ad);
    free(ay);
#endif
    CAMLreturn(Val_unit);
}

/* vector array operations */

CAMLprim value sunml_nvec_linearsumvectorarray(value va, value vax,
					       value vb, value vay, value vaz)
{
    CAMLparam5(va, vax, vb, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *a[3];
    int nvec = sunml_arrays_of_nvectors(a, 3, vax, vay, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VLinearSumVectorArray(nvec, Double_val(va), a[0],
					 Double_val(vb), a[1], a[2]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_scalevectorarray(value vac, value vax, value vaz)
{
    CAMLparam3(vac, vax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VScaleVectorArray(nvec, ac, a[0], a[1]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_constvectorarray(value vc, value vaz)
{
    CAMLparam2(vc, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *az;
    int nvec = sunml_arrays_of_nvectors(&az, 1, vaz);
    if (!nvec) caml_raise_out_of_memory();

    N_VConstVectorArray(nvec, Double_val(vc), az);
    free(az);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_wrmsnormvectorarray(value vax, value vaw, value van)
{
    CAMLparam3(vax, vaw, van);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *an = REAL_ARRAY(van);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);
    if (!nvec) caml_raise_out_of_memory();

    N_VWrmsNormVectorArray(nvec, a[0], a[1], an);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_wrmsnormmaskvectorarray(value vax, value vaw,
						  value vi, value van)
{
    CAMLparam4(vax, vaw, vi, van);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *an = REAL_ARRAY(van);
    N_Vector i = NVEC_VAL(vi);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);
    if (!nvec) caml_raise_out_of_memory();

    N_VWrmsNormMaskVectorArray(nvec, a[0], a[1], i, an);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_scaleaddmultivectorarray(value vaa, value vax,
						   value vaay, value vaaz)
{
    CAMLparam4(vaa, vax, vaay, vaaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *aa = REAL_ARRAY(vaa);
    N_Vector *ax = NULL;
    int nvec = sunml_arrays_of_nvectors(&ax, 1, vax);
    if (!nvec) caml_raise_out_of_memory();
    N_Vector **ayz[2] = { NULL };
    int nvec2, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec2, ayz, 2, vaay, vaaz);
    if (!nsum) caml_raise_out_of_memory();

    N_VScaleAddMultiVectorArray(nvec, nsum, aa, ax, ayz[0], ayz[1]);
    free(ax);
    free(*ayz);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_linearcombinationvectorarray(value vac,
							value vaax, value vaz)
{
    CAMLparam3(vac, vaax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    sunrealtype *ac = REAL_ARRAY(vac);
    N_Vector *az;
    int nvecz __attribute__((unused))
	= sunml_arrays_of_nvectors(&az, 1, vaz);
    N_Vector **aax;
    int nvec, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec, &aax, 1, vaax);
    if (!nsum) caml_raise_out_of_memory();

    N_VLinearCombinationVectorArray(nvec, nsum, ac, aax, az);
    free(aax);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_dotprodlocal(value vx, value vy)
{
    CAMLparam2(vx, vy);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvdotprodlocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VDotProdLocal(x, NVEC_VAL(vy));
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_maxnormlocal(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvmaxnormlocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VMaxNormLocal(x);
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_minlocal(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvminlocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VMinLocal(x);
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_l1normlocal(value vx)
{
    CAMLparam1(vx);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvl1normlocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VL1NormLocal(x);
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_invtestlocal(value vx, value vz)
{
    CAMLparam2(vx, vz);
    sunbooleantype r;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvinvtestlocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VInvTestLocal(x, NVEC_VAL(vz));
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_constrmasklocal(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    sunbooleantype r;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector c = NVEC_VAL(vc);

    if (c->ops->nvconstrmasklocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VConstrMaskLocal(c, NVEC_VAL(vx), NVEC_VAL(vm));
#endif
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_minquotientlocal(value vn, value vd)
{
    CAMLparam2(vn, vd);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector n = NVEC_VAL(vn);

    if (n->ops->nvminquotientlocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VMinQuotientLocal(n, NVEC_VAL(vd));
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_wsqrsumlocal(value vx, value vw)
{
    CAMLparam2(vx, vw);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvwsqrsumlocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VWSqrSumLocal(x, NVEC_VAL(vw));
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_wsqrsummasklocal(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    sunrealtype r = 0.0;
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    if (x->ops->nvwsqrsummasklocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    r = N_VWSqrSumMaskLocal(x, NVEC_VAL(vw), NVEC_VAL(vid));
#endif
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_dotprodmultilocal(value vx, value vay, value vd)
{
    CAMLparam3(vx, vay, vd);
#if 600 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    sunrealtype *d = REAL_ARRAY(vd);
    N_Vector *ay = NULL;
    int nvec;

    if (x->ops->nvdotprodmultilocal == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    nvec = sunml_arrays_of_nvectors(&ay, 1, vay);
    N_VDotProdMultiLocal(nvec, NVEC_VAL(vx), ay, d);
    free(ay);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_dotprodmultiallreduce(value vx, value vd)
{
    CAMLparam2(vx, vd);
#if 600 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    sunrealtype *d = REAL_ARRAY(vd);
    int nvec_total = (Caml_ba_array_val(vd))->dim[0];

    if (x->ops->nvdotprodmultiallreduce == NULL)
	caml_raise_constant(NVEC_EXN(OperationNotProvided));

    N_VDotProdMultiAllReduce(nvec_total, NVEC_VAL(vx), d);
#endif
    CAMLreturn(Val_unit);
}

