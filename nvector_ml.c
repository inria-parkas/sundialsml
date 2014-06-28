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

#include "nvector_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include <nvector/nvector_serial.h>

#ifdef CVODE_ML_WITHMPI
#include <nvector/nvector_parallel.h>

/* Must correspond with camlmpi.h */
#define Comm_val(comm) (*((MPI_Comm *) &Field(comm, 1)))
#endif

/** Generic nvector functions and macros */

static N_Vector alloc_cnvec(size_t content_size, value backlink)
{
    N_Vector nv;

    /* Alloc memory in C heap */
    nv = (N_Vector)malloc(sizeof(struct cnvec));
    if (nv == NULL) return NULL;

    nv->ops = NULL;
    nv->ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
    if (nv->ops == NULL) { free(nv); return(NULL); }

    nv->content = NULL;
    nv->content = (void *) malloc(content_size);
    if (nv->content == NULL) { free(nv->ops); free(nv); return(NULL); }

    NVEC_BACKLINK(nv) = backlink;
    caml_register_global_root(&NVEC_BACKLINK(nv));

    return nv;
}

static void free_cnvec(N_Vector nv)
{
    caml_remove_global_root(&NVEC_BACKLINK(nv));
    free(nv->content);
    free(nv->ops);
    free(nv);
}

static mlsize_t nvec_rough_size =
    sizeof(struct _generic_N_Vector)
    + sizeof(value)
    + sizeof(struct _generic_N_Vector_Ops)
    + 4 * sizeof(void *);

static CAMLprim value val_cnvec(N_Vector nv, void (*finalizer)(value))
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_alloc_final(2, finalizer, nvec_rough_size, nvec_rough_size * 30);
    Store_field(r, 1, (value)nv);

    CAMLreturn(r);
}

static CAMLprim void finalize_cnvec(value vnv)
{
    free_cnvec(NVEC_CVAL(vnv));
}

static void clone_ops(N_Vector dst, N_Vector src)
{
    N_Vector_Ops ops = (N_Vector_Ops) dst->ops;

    ops->nvclone           = src->ops->nvclone;
    ops->nvcloneempty      = src->ops->nvcloneempty;
    ops->nvdestroy         = src->ops->nvdestroy;
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

/** Serial nvectors * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Adapted from sundials-2.5.0/src/nvec_ser/nvector_serial.c:
   N_VCloneEmpty_Serial */
static N_Vector clone_serial(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);

    N_Vector v;
    N_VectorContent_Serial content;

    if (w == NULL) return(NULL);
    w_payload = NVEC_BACKLINK(w);
    struct caml_ba_array *w_ba = Caml_ba_array_val(w_payload);

    /* Create vector (we need not copy the data) */
    v_payload = caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim);
    
    v = alloc_cnvec(sizeof(struct _N_VectorContent_Serial), v_payload);
    if (v == NULL) caml_raise_out_of_memory();
    content = (N_VectorContent_Serial) v->content;

    /* Create vector operation structure */
    clone_ops(v, w);

    /* Create content */
    content->length   = NV_LENGTH_S(w);
    content->own_data = 0;
    content->data     = Caml_ba_data_val(v_payload);

    CAMLreturnT(N_Vector, v);
}

/* Adapted from sundials-2.5.0/src/nvec_ser/nvector_serial.c:
   N_VNewEmpty_Serial */
CAMLprim value ml_nvec_wrap_serial(value payload)
{
    CAMLparam1(payload);
    CAMLlocal1(vnvec);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Serial content;
    long int length = (Caml_ba_array_val(payload))->dim[0];

    /* Create vector */
    nv = alloc_cnvec(sizeof(struct _N_VectorContent_Serial), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_Serial) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_serial;		    /* ours */
    ops->nvcloneempty      = NULL;
    ops->nvdestroy         = free_cnvec;

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

    vnvec = caml_alloc_tuple(2);
    Store_field(vnvec, 0, payload);
    Store_field(vnvec, 1, val_cnvec(nv, finalize_cnvec));

    CAMLreturn(vnvec);
}

/** Parallel nvectors * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef CVODE_ML_WITHMPI

/* Adapted from sundials-2.5.0/src/nvec_par/nvector_parallel.c:
   N_VCloneEmpty_Parallel */
static N_Vector clone_parallel(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);

    N_Vector v;
    N_VectorContent_Parallel content;

    if (w == NULL) return(NULL);
    w_payload = NVEC_BACKLINK(w);
    struct caml_ba_array *w_ba = Caml_ba_array_val(Field(w_payload, 0));

    /* Create vector (we need not copy the data) */
    v_payload = caml_alloc_tuple(3);
    Store_field(v_payload, 0,
		caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim));
    Store_field(v_payload, 1, Field(w_payload, 1))
    Store_field(v_payload, 2, Field(w_payload, 2))
    
    v = alloc_cnvec(sizeof(struct _N_VectorContent_Parallel), v_payload);
    if (v == NULL) caml_raise_out_of_memory();
    content = (N_VectorContent_Parallel) v->content;

    /* Create vector operation structure */
    clone_ops(v, w);

    /* Attach lengths and communicator */
    content->local_length  = NV_LOCLENGTH_P(w);
    content->global_length = NV_GLOBLENGTH_P(w);
    content->comm          = NV_COMM_P(w);
    content->own_data      = 0;
    content->data          = Caml_ba_data_val(Field(v_payload, 0));

    CAMLreturnT(N_Vector, v);
}

/* Adapted from sundials-2.5.0/src/nvec_par/nvector_parallel.c:
   N_VNewEmpty_Parallel */
CAMLprim value ml_nvec_wrap_parallel(value payload)
{
    CAMLparam0(payload);
    CAMLlocal2(vnvec, vlocalba);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Parallel content;
    MPI_Comm comm,
    long int n, nsum, local_length, global_length;

    vlocalba      = Field(payload, 0);
    local_length  = (Caml_ba_array_val(vlocalba))->dim[0];
    global_length = Long_val(Field(payload, 1));
    comm          = Comm_val(Field(payload, 2));

    /* Compute global length as sum of local lengths */
    n = local_length;
    MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
    if (Nsum != global_length)
        caml_raise_constant(
	    *caml_named_value("nvector_parallel_IncorrectGlobalSize"));

    /* Create vector */
    nv = alloc_cnvec(sizeof(struct _N_VectorContent_Parallel), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_Parallel) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_parallel;		    /* ours */
    ops->nvcloneempty      = NULL;
    ops->nvdestroy         = free_cnvec;

    ops->nvspace           = N_VSpace_Parallel;		    /* theirs */
    ops->nvgetarraypointer = N_VGetArrayPointer_Parallel;
    ops->nvsetarraypointer = N_VSetArrayPointer_Parallel;
    ops->nvlinearsum       = N_VLinearSum_Parallel;
    ops->nvconst           = N_VConst_Parallel;
    ops->nvprod            = N_VProd_Parallel;
    ops->nvdiv             = N_VDiv_Parallel;
    ops->nvscale           = N_VScale_Parallel;
    ops->nvabs             = N_VAbs_Parallel;
    ops->nvinv             = N_VInv_Parallel;
    ops->nvaddconst        = N_VAddConst_Parallel;
    ops->nvdotprod         = N_VDotProd_Parallel;
    ops->nvmaxnorm         = N_VMaxNorm_Parallel;
    ops->nvwrmsnormmask    = N_VWrmsNormMask_Parallel;
    ops->nvwrmsnorm        = N_VWrmsNorm_Parallel;
    ops->nvmin             = N_VMin_Parallel;
    ops->nvwl2norm         = N_VWL2Norm_Parallel;
    ops->nvl1norm          = N_VL1Norm_Parallel;
    ops->nvcompare         = N_VCompare_Parallel;
    ops->nvinvtest         = N_VInvTest_Parallel;
    ops->nvconstrmask      = N_VConstrMask_Parallel;
    ops->nvminquotient     = N_VMinQuotient_Parallel;

    /* Attach lengths and communicator */
    content->local_length  = local_length;
    content->global_length = global_length;
    content->comm          = comm
    content->own_data      = 0;
    content->data          = Caml_ba_data_val(vlocalba);

    vnvec = caml_alloc_tuple(2);
    Store_field(vnvec, 0, payload);
    Store_field(vnvec, 1, val_cnvec(nv, finalize_cnvec));

    CAMLreturn(vnvec);
}

#endif

/** Custom nvectors * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define CNVEC_OP_TABLE(nvec)  (*(value*)((nvec)->content))

#define GET_OP(nvec, x) (Field(CNVEC_OP_TABLE(nvec), x))

#define HAS_OP(ops, x)	     (Field(ops, x) != Val_int(0))
#define IS_SOME_OP(nvec, x)  (HAS_OP(CNVEC_OP_TABLE(nvec), x))
#define GET_SOME_OP(nvec, x) (Field(Field(CNVEC_OP_TABLE(nvec), x), 0))

CAMLprim void callml_vdestroy(N_Vector v)
{
    CAMLparam0();
    CAMLlocal1(mlop);

    if (IS_SOME_OP(v, NVECTOR_OPS_NVDESTROY)) {
	mlop = GET_SOME_OP(v, NVECTOR_OPS_NVDESTROY);
	caml_callback(mlop, NVEC_BACKLINK(v));
    }

    caml_remove_generational_global_root(&CNVEC_OP_TABLE(v));
    free_cnvec(v);

    CAMLreturn0;
}

static CAMLprim void finalize_custom_cnvec(value vnv)
{
    callml_vdestroy(NVEC_CVAL(vnv));
}

CAMLprim value ml_nvec_wrap_custom(value mlops, value payload)
{
    CAMLparam2(mlops, payload);
    CAMLlocal1(vcnvec);

    N_Vector nv;
    N_Vector_Ops ops;

    /* Create vector */
    nv = alloc_cnvec(0, payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;

    /* Create vector operation structure */
    ops->nvclone           = callml_vclone;
    ops->nvcloneempty      = NULL;
    ops->nvdestroy         = callml_vdestroy;

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

    ops->nvgetarraypointer = NULL;
    ops->nvsetarraypointer = NULL;

    ops->nvdotprod = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVDOTPROD))
	ops->nvdotprod = callml_vdotprod;

    ops->nvcompare = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVCOMPARE))
	ops->nvcompare = callml_vcompare;

    ops->nvinvtest = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVINVTEST))
	ops->nvinvtest = callml_vinvtest;

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
    caml_register_generational_global_root(&CNVEC_OP_TABLE(nv));

    vcnvec = caml_alloc_tuple(2);
    Store_field(vcnvec, 0, payload);
    Store_field(vcnvec, 1, val_cnvec(nv, finalize_custom_cnvec));

    CAMLreturn(vcnvec);
}

CAMLprim N_Vector callml_vclone(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);

    N_Vector v;

    if (w == NULL) return(NULL);
    w_payload = NVEC_BACKLINK(w);

    /* Create vector */
    v_payload = caml_callback(GET_OP(w, NVECTOR_OPS_NVCLONE), w_payload);

    v = alloc_cnvec(0, v_payload);
    if (v == NULL) caml_raise_out_of_memory();

    /* Create vector operation structure */
    clone_ops(v, w);

    /* Create content */
    v->content = (void *) CNVEC_OP_TABLE(w);
    caml_register_generational_global_root(&CNVEC_OP_TABLE(v));

    CAMLreturnT(N_Vector, v);
}

CAMLprim void callml_vspace(N_Vector v, long int *lrw, long int *liw)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(v, NVECTOR_OPS_NVSPACE);

    r = caml_callback(mlop, NVEC_BACKLINK(v));

    *lrw = Long_val(Field(r, 0));
    *liw = Long_val(Field(r, 1)) + (nvec_rough_size / sizeof(int));

    CAMLreturn0;
}

CAMLprim void callml_vlinearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
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

    caml_callbackN(mlop, 5, args);

    CAMLreturn0;
}

CAMLprim void callml_vconst(realtype c, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(z, NVECTOR_OPS_NVCONST);

    caml_callback2(mlop, caml_copy_double(c), NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim void callml_vprod(N_Vector x, N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVPROD);

    caml_callback3(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(y), NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim void callml_vdiv(N_Vector x, N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVDIV);

    caml_callback3(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(y), NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim void callml_vscale(realtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVSCALE);

    caml_callback3(mlop, caml_copy_double(c),
	    NVEC_BACKLINK(x), NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim void callml_vabs(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVABS);

    caml_callback2(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim void callml_vinv(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVINV);

    caml_callback2(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim void callml_vaddconst(N_Vector x, realtype b, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVADDCONST);

    caml_callback3(mlop, NVEC_BACKLINK(x), caml_copy_double(b),
	    NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim realtype callml_vdotprod(N_Vector x, N_Vector y)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVDOTPROD);

    r = caml_callback2(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(y));

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vmaxnorm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(x, NVECTOR_OPS_NVMAXNORM);

    r = caml_callback(mlop, NVEC_BACKLINK(x));

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vwrmsnorm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(x, NVECTOR_OPS_NVWRMSNORM);

    r = caml_callback2(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(w));

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMMASK);

    r = caml_callback3(mlop, NVEC_BACKLINK(x),
	    NVEC_BACKLINK(w), NVEC_BACKLINK(id));

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vmin(N_Vector x)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(x, NVECTOR_OPS_NVMIN);

    r = caml_callback(mlop, NVEC_BACKLINK(x));

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vwl2norm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWL2NORM);

    r = caml_callback2(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(w));

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vl1norm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVL1NORM);

    r = caml_callback(mlop, NVEC_BACKLINK(x));

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim void callml_vcompare(realtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVCOMPARE);

    caml_callback3(mlop, caml_copy_double(c),
	    NVEC_BACKLINK(x), NVEC_BACKLINK(z));

    CAMLreturn0;
}

CAMLprim booleantype callml_vinvtest(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVINVTEST);

    r = caml_callback2(mlop, NVEC_BACKLINK(x), NVEC_BACKLINK(z));

    CAMLreturnT(booleantype, Bool_val(r));
}

CAMLprim booleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVCONSTRMASK);

    r = caml_callback3(mlop, NVEC_BACKLINK(c),
	    NVEC_BACKLINK(x), NVEC_BACKLINK(m));

    CAMLreturnT(booleantype, Bool_val(r));
}

CAMLprim realtype callml_vminquotient(N_Vector num, N_Vector denom)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(num, NVECTOR_OPS_NVMINQUOTIENT);

    r = caml_callback2(mlop, NVEC_BACKLINK(num), NVEC_BACKLINK(denom));

    CAMLreturnT(realtype, Double_val(r));
}

