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

/** Interface to underlying serial nvector functions */

CAMLprim value sunml_nvec_ser_n_vlinearsum(value va, value vx, value vb, value vy,
				       value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_VLinearSum_Serial(Double_val(va), NVEC_VAL(vx), Double_val(vb),
			NVEC_VAL(vy), NVEC_VAL(vz));
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
    N_VProd_Serial(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VDiv_Serial(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VScale_Serial(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VAbs_Serial(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VInv_Serial(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_VAddConst_Serial(NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    realtype r = N_VDotProd_Serial(NVEC_VAL(vx), NVEC_VAL(vy));
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
    realtype r = N_VWrmsNorm_Serial(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_ser_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    realtype r = N_VWrmsNormMask_Serial(NVEC_VAL(vx), NVEC_VAL(vw),
					NVEC_VAL(vid));
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
    realtype r = N_VWL2Norm_Serial(NVEC_VAL(vx), NVEC_VAL(vw));
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
    N_VCompare_Serial(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_ser_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    booleantype r = N_VInvTest_Serial(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_ser_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    booleantype r = N_VConstrMask_Serial(NVEC_VAL(vc), NVEC_VAL(vx),
					 NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_ser_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    realtype r = N_VMinQuotient_Serial(NVEC_VAL(vnum), NVEC_VAL(vdenom));
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

