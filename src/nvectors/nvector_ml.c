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
    nv = (N_Vector)malloc(sizeof(struct cnvec));
    if (nv == NULL) return NULL;

    nv->ops = NULL;
    nv->ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
    if (nv->ops == NULL) { free(nv); return(NULL); }

    nv->content = NULL;
    if (content_size != 0) {
	nv->content = (void *) malloc(content_size);
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

CAMLprim value sunml_alloc_caml_nvec(N_Vector nv, void (*finalizer)(value))
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_final(1, finalizer, 1, 30);
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
    N_Vector_Ops ops = (N_Vector_Ops) dst->ops;

    ops->nvclone           = src->ops->nvclone;
    ops->nvcloneempty      = src->ops->nvcloneempty;
    ops->nvdestroy         = src->ops->nvdestroy;
#if SUNDIALS_LIB_VERSION >= 270
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
#if SUNDIALS_LIB_VERSION >= 400
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
#if SUNDIALS_LIB_VERSION >= 290
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
	case SUNDIALS_NVEC_OPENMPDEV:
	    vr = Val_int(VARIANT_NVECTOR_ID_TAG_OPENMPDEV);
	    break;
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

    /* Create content */
    content->length   = NV_LENGTH_S(w);
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
CAMLprim value sunml_nvec_wrap_serial(value payload, value checkfn)
{
    CAMLparam2(payload, checkfn);
    CAMLlocal1(vnvec);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Serial content;
    long int length = (Caml_ba_array_val(payload))->dim[0];

    /* Create vector */
    nv = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Serial), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_Serial) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_serial;		    /* ours */
    ops->nvcloneempty      = clone_empty_serial;	    /* ours */
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = sunml_free_cnvec;
#if SUNDIALS_LIB_VERSION >= 270
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

#if SUNDIALS_LIB_VERSION >= 400
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

    /* Create content */
    content->length   = length;
    content->own_data = 0;
    content->data     = Caml_ba_data_val(payload);

    vnvec = caml_alloc_tuple(3);
    Store_field(vnvec, 0, payload);
    Store_field(vnvec, 1, sunml_alloc_caml_nvec(nv, sunml_finalize_caml_nvec));
    Store_field(vnvec, 2, checkfn);

    CAMLreturn(vnvec);
}

/** Custom nvectors * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define CNVEC_OP_TABLE(nvec)  ((nvec)->content)

#define GET_OP(nvec, x) (Field((value)CNVEC_OP_TABLE(nvec), x))

#define HAS_OP(ops, x)	     (Field(ops, x) != Val_int(0))
#define IS_SOME_OP(nvec, x)  (HAS_OP(CNVEC_OP_TABLE(nvec), x))
#define GET_SOME_OP(nvec, x) (Field(Field(CNVEC_OP_TABLE(nvec), x), 0))

static void free_custom_cnvec(N_Vector v)
{
    caml_remove_generational_global_root((value *)&CNVEC_OP_TABLE(v));
    v->content = NULL;
    sunml_free_cnvec(v);
}

static void finalize_custom_caml_nvec(value vnv)
{
    free_custom_cnvec (NVEC_CVAL(vnv));
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
static void callml_vlinearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
static void callml_vconst(realtype c, N_Vector z);
static void callml_vprod(N_Vector x, N_Vector y, N_Vector z);
static void callml_vdiv(N_Vector x, N_Vector y, N_Vector z);
static void callml_vscale(realtype c, N_Vector x, N_Vector z);
static void callml_vabs(N_Vector x, N_Vector z);
static void callml_vinv(N_Vector x, N_Vector z);
static void callml_vaddconst(N_Vector x, realtype b, N_Vector z);
static realtype callml_vdotprod(N_Vector x, N_Vector y);
static realtype callml_vmaxnorm(N_Vector x);
static realtype callml_vwrmsnorm(N_Vector x, N_Vector w);
static realtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id);
static realtype callml_vmin(N_Vector x);
static realtype callml_vwl2norm(N_Vector x, N_Vector w);
static realtype callml_vl1norm(N_Vector x);
static void callml_vcompare(realtype c, N_Vector x, N_Vector z);
static booleantype callml_vinvtest(N_Vector x, N_Vector z);
static booleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m);
static realtype callml_vminquotient(N_Vector num, N_Vector denom);

/* Custom fused vector operations */
static int callml_vlinearcombination(int nvec, realtype* c,
				     N_Vector* V, N_Vector z);
static int callml_vscaleaddmulti(int nvec, realtype* a,
			            N_Vector x, N_Vector* Y, N_Vector* Z);
static int callml_vdotprodmulti(int nvec, N_Vector x, N_Vector *Y,
				   realtype* dotprods);

/* Custom vector array operations */
static int callml_vlinearsumvectorarray(int nvec, realtype a, N_Vector* X,
					realtype b, N_Vector* Y, N_Vector* Z);
static int callml_vscalevectorarray(int nvec, realtype* c, N_Vector* X,
				    N_Vector* Z);
static int callml_vconstvectorarray(int nvecs, realtype c, N_Vector* Z);
static int callml_vwrmsnormvectorarray(int nvecs, N_Vector* X,
				       N_Vector* W, realtype* nrm);
static int callml_vwrmsnormmaskvectorarray(int nvecs, N_Vector* X, N_Vector* W,
					   N_Vector id, realtype* nrm);
static int callml_vscaleaddmultivectorarray(int nvec, int nsum, realtype* a,
					    N_Vector* X, N_Vector** Y,
					    N_Vector** Z);
static int callml_vlinearcombinationvectorarray(int nvec, int nsum, realtype* c,
						N_Vector** X, N_Vector* Z);

/* Creation from OCaml. */
CAMLprim value sunml_nvec_wrap_custom(value mlops, value payload, value checkfn)
{
    CAMLparam3(mlops, payload, checkfn);
    CAMLlocal1(vcnvec);

    N_Vector nv;
    N_Vector_Ops ops;

    /* Create vector */
    nv = sunml_alloc_cnvec(0, payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;

    /* Create vector operation structure */
    ops->nvclone           = callml_vclone;
    ops->nvcloneempty      = NULL;
    ops->nvdestroy         = free_custom_cnvec;
#if SUNDIALS_LIB_VERSION >= 270
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

#if SUNDIALS_LIB_VERSION >= 400
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

    /* Create content */
    nv->content = (void *)mlops;
    caml_register_generational_global_root((value *)&CNVEC_OP_TABLE(nv));

    vcnvec = caml_alloc_tuple(3);
    Store_field(vcnvec, 0, payload);
    Store_field(vcnvec, 1, sunml_alloc_caml_nvec(nv, finalize_custom_caml_nvec));
    Store_field(vcnvec, 2, checkfn);

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
					"user-defined n_vclone");
	CAMLreturnT (N_Vector, NULL);
    }
    v_payload = r;
    /* Done processing r.  Now it's OK to trigger GC.  */

    v = sunml_alloc_cnvec(0, v_payload);
    if (v == NULL)
	CAMLreturnT (N_Vector, NULL);

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

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
					"user-defined n_vspace");
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

static void callml_vlinearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
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
					"user-defined n_vlinearsum");
    }

    CAMLreturn0;
}

static void callml_vconst(realtype c, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vc);

    vc = caml_copy_double (c);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(GET_OP(z, NVECTOR_OPS_NVCONST),
				 vc, NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vconst");

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
					"user-defined n_vprod");

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
					"user-defined n_vdiv");

    CAMLreturn0;
}

static void callml_vscale(realtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vc);

    vc = caml_copy_double(c);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(GET_OP(x, NVECTOR_OPS_NVSCALE), vc,
				 NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vscale");

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
					"user-defined n_vabs");

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
					"user-defined n_vinv");

    CAMLreturn0;
}

static void callml_vaddconst(N_Vector x, realtype b, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vb);

    vb = caml_copy_double(b);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (GET_OP(x, NVECTOR_OPS_NVADDCONST),
				  NVEC_BACKLINK(x), vb, NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vaddconst");

    CAMLreturn0;
}

static realtype callml_vdotprod(N_Vector x, N_Vector y)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVDOTPROD);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(y));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vdotprod");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

static realtype callml_vmaxnorm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVMAXNORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vmaxnorm");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

static realtype callml_vwrmsnorm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVWRMSNORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(w));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vwrmsnorm");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

static realtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMMASK);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(x),
				  NVEC_BACKLINK(w), NVEC_BACKLINK(id));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vwrmsnormmask");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

static realtype callml_vmin(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVMIN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vmin");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

static realtype callml_vwl2norm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWL2NORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(w));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vwl2norm");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

static realtype callml_vl1norm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVL1NORM);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (mlop, NVEC_BACKLINK(x));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vl1norm");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

static void callml_vcompare(realtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(vc);

    vc = caml_copy_double(c);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (GET_OP(x, NVECTOR_OPS_NVCOMPARE), vc,
				  NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r))
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vcompare");

    CAMLreturn0;
}

static booleantype callml_vinvtest(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVINVTEST);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vinvtest");
	CAMLreturnT(booleantype, 0);
    }

    CAMLreturnT(booleantype, Bool_val(r));
}

static booleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVCONSTRMASK);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(c),
				  NVEC_BACKLINK(x), NVEC_BACKLINK(m));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vconstrmask");
	CAMLreturnT(booleantype, 0);
    }

    CAMLreturnT(booleantype, Bool_val(r));
}

static realtype callml_vminquotient(N_Vector num, N_Vector denom)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(num, NVECTOR_OPS_NVMINQUOTIENT);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, NVEC_BACKLINK(num),
				  NVEC_BACKLINK(denom));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vminquotient");
	CAMLreturnT(realtype, nan(""));
    }

    CAMLreturnT(realtype, Double_val(r));
}

// NB: Normally, we should worry about relinquishing the elements of vy
// after we are finished using them (so as not to block the GC), but we
// instead make the assumption that these elements come from 'within'
// Sundials and thus that they would anyway not be GC-ed.
static value wrap_to_nvector_table(int n, N_Vector *y)
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

static value wrap_to_nvector_tables(int n1, int n2, N_Vector **yy)
{
    CAMLparam0();
    CAMLlocal1(vyy);
    int i;

    vyy = caml_alloc_tuple(n1);

    for (i = 0; i < n1; ++i) {
	Store_field(vyy, i, wrap_to_nvector_table(n2, yy[i]));
    }

    CAMLreturn (vyy);
}

/* fused vector operations */
static int callml_vlinearcombination(int nvec, realtype* c,
				     N_Vector* V, N_Vector z)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal3(mlop, vc, vv);
    intnat n = nvec;

    mlop = GET_SOME_OP(z, NVECTOR_OPS_NVLINEARCOMBINATION);
    vv = wrap_to_nvector_table(nvec, V);
    vc = caml_ba_alloc(BIGARRAY_FLOAT, 1, c, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vc, vv, NVEC_BACKLINK(z));
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vlinearcombination");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vscaleaddmulti(int nvec, realtype* a,
			         N_Vector x, N_Vector* Y, N_Vector* Z)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(mlop);
    CAMLlocalN(args, 4);
    intnat n = nvec;

    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVSCALEADDMULTI);
    args[0] = caml_ba_alloc(BIGARRAY_FLOAT, 1, a, &n);
    args[1] = NVEC_BACKLINK(x);
    args[2] = wrap_to_nvector_table(nvec, Y);
    args[3] = wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 4, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vscaleaddmulti");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vdotprodmulti(int nvec, N_Vector x, N_Vector *Y,
				realtype* dotprods)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal3(mlop, vy, vdotprods);
    intnat n = nvec;

    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVDOTPRODMULTI);
    vy = wrap_to_nvector_table(nvec, Y);
    vdotprods = caml_ba_alloc(BIGARRAY_FLOAT, 1, dotprods, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, NVEC_BACKLINK(x), vy, vdotprods);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vdotprodmulti");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}


/* vector array operations */
static int callml_vlinearsumvectorarray(int nvec, realtype a, N_Vector* X,
					   realtype b, N_Vector* Y,
                                           N_Vector* Z)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(mlop);
    CAMLlocalN(args, 5);

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVLINEARSUMVECTORARRAY);
    args[0] = caml_copy_double(a);
    args[1] = wrap_to_nvector_table(nvec, X);
    args[2] = caml_copy_double(b);
    args[3] = wrap_to_nvector_table(nvec, Y);
    args[4] = wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 5, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vlinearsumvectorarray");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vscalevectorarray(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal4(mlop, vc, vx, vz);
    intnat n = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVSCALEVECTORARRAY);
    vc = caml_ba_alloc(BIGARRAY_FLOAT, 1, c, &n);
    vx = wrap_to_nvector_table(nvec, X);
    vz = wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vc, vx, vz);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vscalevectorarray");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vconstvectorarray(int nvec, realtype c, N_Vector* Z)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal3(mlop, vc, vz);

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(Z[0], NVECTOR_OPS_NVCONSTVECTORARRAY);
    vc = caml_copy_double(c);
    vz = wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (mlop, vc, vz);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vconstvectorarray");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vwrmsnormvectorarray(int nvec, N_Vector* X,
				       N_Vector* W, realtype* nrm)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal4(mlop, vx, vw, vnrm);
    intnat n = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVWRMSNORMVECTORARRAY);
    vx = wrap_to_nvector_table(nvec, X);
    vw = wrap_to_nvector_table(nvec, W);
    vnrm = caml_ba_alloc(BIGARRAY_FLOAT, 1, nrm, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vx, vw, vnrm);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined n_vwrmsnormvectorarray");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vwrmsnormmaskvectorarray(int nvec, N_Vector* X, N_Vector* W,
					   N_Vector id, realtype* nrm)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(mlop);
    CAMLlocalN(args, 4);
    intnat n = nvec;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVWRMSNORMMASKVECTORARRAY);
    args[0] = wrap_to_nvector_table(nvec, X);
    args[1] = wrap_to_nvector_table(nvec, W);
    args[2] = NVEC_BACKLINK(id);
    args[3] = caml_ba_alloc(BIGARRAY_FLOAT, 1, nrm, &n);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 4, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
				    "user-defined n_vwrmsnormmaskvectorarray");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vscaleaddmultivectorarray(int nvec, int nsum, realtype* a,
					       N_Vector* X, N_Vector** Y,
					       N_Vector** Z)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(mlop);
    CAMLlocalN(args, 4);
    intnat n = nsum;

    if (nvec <= 0) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(X[0], NVECTOR_OPS_NVSCALEADDMULTIVECTORARRAY);
    args[0] = caml_ba_alloc(BIGARRAY_FLOAT, 1, a, &n);
    args[1] = wrap_to_nvector_table(nvec, X);
    args[2] = wrap_to_nvector_tables(nsum, nvec, Y);
    args[3] = wrap_to_nvector_tables(nsum, nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (mlop, 4, args);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
				    "user-defined n_vscaleaddmultivectorarray");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

static int callml_vlinearcombinationvectorarray(int nvec, int nsum, realtype* c,
						N_Vector** X, N_Vector* Z)
{
    CAMLparam0();
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal4(mlop, vc, vxx, vz);
    intnat n2 = nsum;

    if ((nvec <= 0) || (nsum <= 0)) CAMLreturnT(int, 1);

    mlop = GET_SOME_OP(Z[0], NVECTOR_OPS_NVLINEARCOMBINATIONVECTORARRAY);
    vc = caml_ba_alloc(BIGARRAY_FLOAT, 1, c, &n2);
    vxx = wrap_to_nvector_tables(nsum, nvec, X);
    vz = wrap_to_nvector_table(nvec, Z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (mlop, vc, vxx, vz);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
			    "user-defined n_vlinearcombinationvectorarray");
	CAMLreturnT(int, 0);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturnT(int, 1);
}

/** Interface to underlying serial nvector functions */

CAMLprim value sunml_nvec_ser_n_vlinearsum(value va, value vx, value vb, value vy,
				       value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(y) != NV_LENGTH_S(x) || NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vlinearsum");
#endif

    N_VLinearSum_Serial(Double_val(va), x, Double_val(vb), y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vconst(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst_Serial(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vprod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(y) != NV_LENGTH_S(x) || NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vprod");
#endif

    N_VProd_Serial(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(y) != NV_LENGTH_S(x) || NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vdiv");
#endif

    N_VDiv_Serial(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vscale");
#endif

    N_VScale_Serial(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vabs");
#endif

    N_VAbs_Serial(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vinv");
#endif

    N_VInv_Serial(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vaddconst");
#endif

    N_VAddConst_Serial(x, Double_val(vb), z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(y) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vdotprod");
#endif

    realtype r = N_VDotProd_Serial(x, y);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vmaxnorm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMaxNorm_Serial(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vwrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(w) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vwrmsnorm");
#endif

    realtype r = N_VWrmsNorm_Serial(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(w) != NV_LENGTH_S(x) || NV_LENGTH_S(w) != NV_LENGTH_S(id))
	caml_invalid_argument("Nvector_serial.n_vwrmsnormmask");
#endif

    realtype r = N_VWrmsNormMask_Serial(x, w, id);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vmin(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMin_Serial(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vwl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(w) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vwl2norm");
#endif

    realtype r = N_VWL2Norm_Serial(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vl1norm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VL1Norm_Serial(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vcompare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vcompare");
#endif

    N_VCompare_Serial(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(z) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vinvtest");
#endif

    booleantype r = N_VInvTest_Serial(x, z);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_ser_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    N_Vector c = NVEC_VAL(vc);
    N_Vector x = NVEC_VAL(vx);
    N_Vector m = NVEC_VAL(vm);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(x) != NV_LENGTH_S(c) || NV_LENGTH_S(m) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vconstrmask");
#endif

    booleantype r = N_VConstrMask_Serial(c, x, m);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_ser_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    N_Vector num = NVEC_VAL(vnum);
    N_Vector denom = NVEC_VAL(vdenom);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LENGTH_S(num) != NV_LENGTH_S(denom))
	caml_invalid_argument("Nvector_serial.n_vminquotient");
#endif

    realtype r = N_VMinQuotient_Serial(num, denom);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vspace(value vx)
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
#if SUNDIALS_ML_SAFE == 1
    int nx = -1;
#endif

    va_start(valist, n);
    vax = va_arg(valist, value);
    nvecs = Wosize_val(vax);
#if SUNDIALS_ML_SAFE == 1
    if (nvecs < 1) goto error;
#endif

    px = calloc(n * nvecs, sizeof(N_Vector));
    if (px == NULL) goto error;

    j = 0;
    while (1) {
	r[j] = px;

	for (i=0; i < nvecs; ++i, ++px) {
	    *px = NVEC_VAL(Field(vax, i));
#if SUNDIALS_ML_SAFE == 1
	    if (NV_LENGTH_S((*px)) != nx) {
		if (nx != -1) goto error;
		nx = NV_LENGTH_S((*px));
	    }
#endif
	}
	++j;
	if (j == n) break;

	vax = va_arg(valist, value);
#if SUNDIALS_ML_SAFE == 1
	if (Wosize_val(vax) != nvecs) goto error;
#endif
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
#if SUNDIALS_ML_SAFE == 1
    int nx = -1;
#endif
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
#if SUNDIALS_ML_SAFE == 1
	    if (Wosize_val(vax) != nc) goto error;
#endif
	    for (c=0; c < nc; ++c, ++pd) {
		*pd = NVEC_VAL(Field(vax, c));
#if SUNDIALS_ML_SAFE == 1
		if (NV_LENGTH_S((*pd)) != nx) {
		    if (nx != -1) goto error;
		    nx = NV_LENGTH_S((*pd));
		}
#endif
	    }
	}

	++j;
	if (j == n) break;

	vaax = va_arg(valist, value);
#if SUNDIALS_ML_SAFE == 1
	if (Wosize_val(vaax) != nr) goto error;
#endif
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

CAMLprim value sunml_nvec_ser_n_vlinearcombination(value vac, value vax,
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
	    || NV_LENGTH_S(ax[0]) != NV_LENGTH_S(z))
	caml_invalid_argument("Nvector_serial.n_vlinearcombination");
#endif

    N_VLinearCombination_Serial(nvec, ac, ax, z);
    free(ax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vscaleaddmulti(value vac, value vx, value vay,
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
	    || NV_LENGTH_S(a[0][0]) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vscaleaddmulti");
#endif

    N_VScaleAddMulti_Serial(nvec, ac, x, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vdotprodmulti(value vx, value vay, value vad)
{
    CAMLparam3(vx, vay, vad);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ad = REAL_ARRAY(vad);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *ay;
    int nvec = sunml_arrays_of_nvectors(&ay, 1, vay);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vad) < nvec
		|| NV_LENGTH_S(ay[0]) != NV_LENGTH_S(x))
	caml_invalid_argument("Nvector_serial.n_vdotprodmulti");
#endif

    N_VDotProdMulti_Serial(nvec, x, ay, ad);
    free(ay);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}


/* vector array operations */

CAMLprim value sunml_nvec_ser_n_vlinearsumvectorarray(value va, value vax,
						      value vb, value vay,
						      value vaz)
{
    CAMLparam5(va, vax, vb, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *a[3];
    int nvec = sunml_arrays_of_nvectors(a, 3, vax, vay, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_serial.n_vlinearsumvectorarray");
#endif

    N_VLinearSumVectorArray_Serial(nvec, Double_val(va), a[0],
					 Double_val(vb), a[1], a[2]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vscalevectorarray(value vac, value vax,
						  value vaz)
{
    CAMLparam3(vac, vax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vac) < nvec)
	caml_invalid_argument("Nvector_serial.n_vscalevectorarray");
#endif

    N_VScaleVectorArray_Serial(nvec, ac, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vconstvectorarray(value vc, value vaz)
{
    CAMLparam2(vc, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *az;
    int nvec = sunml_arrays_of_nvectors(&az, 1, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_serial.n_vconstvectorarray");
#endif

    N_VConstVectorArray_Serial(nvec, Double_val(vc), az);
    free(az);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vwrmsnormvectorarray(value vax, value vaw,
						     value van)
{
    CAMLparam3(vax, vaw, van);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *an = REAL_ARRAY(van);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(van) < nvec)
	caml_invalid_argument("Nvector_serial.n_vconstvectorarray");
#endif

    N_VWrmsNormVectorArray_Serial(nvec, a[0], a[1], an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vwrmsnormmaskvectorarray(value vax,
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
		|| NV_LENGTH_S(i) != NV_LENGTH_S(a[0][0]))
	caml_invalid_argument("Nvector_serial.n_vconstvectorarray");
#endif

    N_VWrmsNormMaskVectorArray_Serial(nvec, a[0], a[1], i, an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vscaleaddmultivectorarray(value vaa,
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
	caml_invalid_argument("Nvector_serial.n_vscaleaddmultivectorarray");
    }
#endif

    N_VScaleAddMultiVectorArray_Serial(nvec, nsum, aa, ax, ayz[0], ayz[1]);
    free(ax);
    free(*ayz);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vlinearcombinationvectorarray(value vac,
							value vaax, value vaz)
{
    CAMLparam3(vac, vaax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *az;
    int nvecz = sunml_arrays_of_nvectors(&az, 1, vaz);
    N_Vector **aax;
    int nvec, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec, &aax, 1, vaax);

#if SUNDIALS_ML_SAFE == 1
    if (!nvecz || !nsum || nvec > nvecz || ARRAY1_LEN(vac) < nsum) {
	if (az != NULL) free(az);
	if (aax != NULL) free(aax);
	caml_invalid_argument("Nvector_serial.n_vlinearcombinationvectorarray");
    }
#endif

    N_VLinearCombinationVectorArray_Serial(nvec, nsum, ac, aax, az);
    free(aax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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
    booleantype r = 1;
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

