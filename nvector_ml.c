/***********************************************************************
 *                                                                     *
 *               OCaml interface to (serial) Sundials                  *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

/* OCaml interface for custom NVectors. */

/*
    ML NVectors can be created in two different ways:

    1. From OCaml by Nvector_array.make or Nvector_array.wrap.

       on creation:
	   ml_nvec_new is invoked (it allocates memory in the C heap)
       on deletion:
	   finalize_nec is invoked to free memory
	   (it, in turn, invokes callml_vdestroy)

    2. From C (Sundials) by NV_Clone.

       on creation:
           callml_vclone is invoked (it allocates memory in the C heap)
       on deletion:
	    callml_vdestroy is invoked to free memory

    The underlying N_Vector is allocated in the C heap (rather than inside a
    custom value) so that it will not be moved by the garbage collector.
    
    An ML NVector (in the C heap) points to two OCaml values:
	content.data
	content.callback

    These are registered as global roots so that the garbage collector does
    not destroy the table of operations or the underlying data. These
    registrations also ensure that the values are updated correctly if the
    underlying data is moved (compacted).
 */

#include "nvector_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>

enum nvector_ops_tag {
  NVECTOR_OPS_NVCLONE = 0,
  NVECTOR_OPS_NVDESTROY,
  NVECTOR_OPS_NVSPACE,
  NVECTOR_OPS_NVLINEARSUM,
  NVECTOR_OPS_NVCONST,
  NVECTOR_OPS_NVPROD,
  NVECTOR_OPS_NVDIV,
  NVECTOR_OPS_NVSCALE,
  NVECTOR_OPS_NVABS,
  NVECTOR_OPS_NVINV,
  NVECTOR_OPS_NVADDCONST,
  NVECTOR_OPS_NVMAXNORM,
  NVECTOR_OPS_NVWRMSNORM,
  NVECTOR_OPS_NVMIN,

  NVECTOR_OPS_NVDOTPROD,
  NVECTOR_OPS_NVCOMPARE,
  NVECTOR_OPS_NVINVTEST,

  NVECTOR_OPS_NVWL2NORM,
  NVECTOR_OPS_NVL1NORM,
  NVECTOR_OPS_NVWRMSNORMMASK,
  NVECTOR_OPS_NVCONSTRMASK,
  NVECTOR_OPS_NVMINQUOTIENT,
  NVECTOR_OPS_SIZE
};

#define GET_OP(nvec, x) (Field(NVEC_CONTENT(nvec)->callbacks, x))
#define GET_DATA(nvec) (NVEC_CONTENT(nvec)->data)
#define EMPTY_DATA Val_int(0)

#define HAS_OP(ops, x) (Field(ops, x) != Val_int(0))
#define IS_SOME_OP(nvec, x) (HAS_OP(NVEC_CONTENT(nvec)->callbacks, x))
#define GET_SOME_OP(nvec, x) (Field(Field(NVEC_CONTENT(nvec)->callbacks, x), 0))

static mlsize_t nvec_rough_size =
    sizeof(struct _generic_N_Vector_Ops) + sizeof(struct _ml_nvec_content);

CAMLprim void callml_vdestroy(N_Vector v)
{
    CAMLparam0();
    CAMLlocal1(mlop);

    if (GET_DATA(v) != EMPTY_DATA) {
	if (IS_SOME_OP(v, NVECTOR_OPS_NVDESTROY)) {
	    mlop = GET_SOME_OP(v, NVECTOR_OPS_NVDESTROY);
	    caml_callback(mlop, GET_DATA(v));
	}
    }

    caml_remove_global_root((&NVEC_CONTENT(v)->data));
    caml_remove_generational_global_root((&NVEC_CONTENT(v)->callbacks));
    free(v);

    CAMLreturn0;
}

static void finalize_nvec(value nvec)
{
    callml_vdestroy(NVEC_VAL(nvec));
}

static N_Vector malloc_nvec()
{
    char* mem = NULL; // pointer arithmetic in bytes
    N_Vector nv;

    /* Alloc memory in C heap */
    mem = malloc(sizeof(*nv) + sizeof(*nv->ops) + sizeof(ml_nvec_content));
    if (mem == NULL) return NULL;

    nv = (N_Vector) mem;
    nv->ops     = (N_Vector_Ops) (mem + sizeof(*nv));
    nv->content = (ml_nvec_content) (mem + sizeof(*nv) + sizeof(*nv->ops));

    return nv;
}

CAMLprim value ml_nvec_new(value mlops, value data)
{
    CAMLparam2(mlops, data);
    CAMLlocal1(rv);
    ml_nvec_content content = NULL;

    N_Vector nv = malloc_nvec();
    if (nv == NULL) caml_raise_out_of_memory();
    content = NVEC_CONTENT(nv);

    rv = caml_alloc_final(2, finalize_nvec,
			  nvec_rough_size, nvec_rough_size * 50);
    Store_field(rv, 1, (value)nv);

    /* Create vector operation structure */
    nv->ops->nvclone           = callml_vclone;
    nv->ops->nvcloneempty      = callml_vcloneempty;
    nv->ops->nvdestroy         = callml_vdestroy;

    nv->ops->nvspace = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVSPACE))
	nv->ops->nvspace = callml_vspace;

    nv->ops->nvlinearsum       = callml_vlinearsum;
    nv->ops->nvconst           = callml_vconst;
    nv->ops->nvprod            = callml_vprod;
    nv->ops->nvdiv             = callml_vdiv;
    nv->ops->nvscale           = callml_vscale;
    nv->ops->nvabs             = callml_vabs;
    nv->ops->nvinv             = callml_vinv;
    nv->ops->nvaddconst        = callml_vaddconst;
    nv->ops->nvmaxnorm         = callml_vmaxnorm;
    nv->ops->nvwrmsnorm        = callml_vwrmsnorm;
    nv->ops->nvmin             = callml_vmin;

    nv->ops->nvgetarraypointer = NULL;
    nv->ops->nvsetarraypointer = NULL;

    nv->ops->nvdotprod = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVDOTPROD))
	nv->ops->nvdotprod = callml_vdotprod;

    nv->ops->nvcompare = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVCOMPARE))
	nv->ops->nvcompare = callml_vcompare;

    nv->ops->nvinvtest = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVINVTEST))
	nv->ops->nvinvtest = callml_vinvtest;

    nv->ops->nvwl2norm = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVWL2NORM))
	nv->ops->nvwl2norm = callml_vwl2norm;

    nv->ops->nvl1norm = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVL1NORM))
	nv->ops->nvl1norm = callml_vl1norm;

    nv->ops->nvwrmsnormmask = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVWRMSNORMMASK))
	nv->ops->nvwrmsnormmask = callml_vwrmsnormmask;

    nv->ops->nvconstrmask = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVCONSTRMASK))
	nv->ops->nvconstrmask = callml_vconstrmask;

    nv->ops->nvminquotient = NULL;
    if (HAS_OP(mlops, NVECTOR_OPS_NVMINQUOTIENT))
	nv->ops->nvminquotient = callml_vminquotient;

    /* Create content */
    content->callbacks = mlops;
    content->data      = data;
    caml_register_generational_global_root(&content->callbacks);
    caml_register_global_root(&content->data);

    CAMLreturn(rv);
}

CAMLprim value ml_nvec_data(N_Vector v)
{
    CAMLparam0();
    CAMLreturn(GET_DATA(v));
}

CAMLprim N_Vector callml_vcloneempty(N_Vector w)
{
    /* Create vector */
    N_Vector v = malloc_nvec();
    if (v == NULL) return NULL;

    ml_nvec_content content = NVEC_CONTENT(v);

    /* Create vector operation structure */
    v->ops->nvclone           = w->ops->nvclone;
    v->ops->nvcloneempty      = w->ops->nvcloneempty;
    v->ops->nvdestroy         = w->ops->nvdestroy;
    v->ops->nvspace           = w->ops->nvspace;
    v->ops->nvgetarraypointer = w->ops->nvgetarraypointer;
    v->ops->nvsetarraypointer = w->ops->nvsetarraypointer;
    v->ops->nvlinearsum       = w->ops->nvlinearsum;
    v->ops->nvconst           = w->ops->nvconst;  
    v->ops->nvprod            = w->ops->nvprod;   
    v->ops->nvdiv             = w->ops->nvdiv;
    v->ops->nvscale           = w->ops->nvscale; 
    v->ops->nvabs             = w->ops->nvabs;
    v->ops->nvinv             = w->ops->nvinv;
    v->ops->nvaddconst        = w->ops->nvaddconst;
    v->ops->nvdotprod         = w->ops->nvdotprod;
    v->ops->nvmaxnorm         = w->ops->nvmaxnorm;
    v->ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
    v->ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
    v->ops->nvmin             = w->ops->nvmin;
    v->ops->nvwl2norm         = w->ops->nvwl2norm;
    v->ops->nvl1norm          = w->ops->nvl1norm;
    v->ops->nvcompare         = w->ops->nvcompare;    
    v->ops->nvinvtest         = w->ops->nvinvtest;
    v->ops->nvconstrmask      = w->ops->nvconstrmask;
    v->ops->nvminquotient     = w->ops->nvminquotient;

    /* Create content */
    content->callbacks = NVEC_CONTENT(w)->callbacks;
    content->data = EMPTY_DATA;
    caml_register_generational_global_root(&content->callbacks);
    caml_register_global_root(&content->data);

    return v;
}

CAMLprim N_Vector callml_vclone(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    N_Vector v = NULL;

    mlop = GET_OP(w, NVECTOR_OPS_NVCLONE);

    r = caml_callback(mlop, GET_DATA(w));

    v = callml_vcloneempty(w);
    if (v != NULL) {
	NVEC_CONTENT(v)->data = r;
    }

    CAMLreturnT(N_Vector, v);
}

CAMLprim void callml_vspace(N_Vector v, long int *lrw, long int *liw)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(v, NVECTOR_OPS_NVSPACE);

    r = caml_callback(mlop, GET_DATA(v));

    *lrw = Long_val(Field(r, 0));
    *liw = Long_val(Field(r, 1));

    CAMLreturn0;
}

CAMLprim void callml_vlinearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    CAMLlocalN(args, 5);

    mlop = GET_OP(x, NVECTOR_OPS_NVLINEARSUM);

    args[0] = caml_copy_double(a);
    args[1] = NVEC_CONTENT(x)->data;
    args[2] = caml_copy_double(b);
    args[3] = NVEC_CONTENT(y)->data;
    args[4] = NVEC_CONTENT(z)->data;

    caml_callbackN(mlop, 5, args);

    CAMLreturn0;
}

CAMLprim void callml_vconst(realtype c, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(z, NVECTOR_OPS_NVCONST);

    caml_callback2(mlop, caml_copy_double(c), NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim void callml_vprod(N_Vector x, N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVPROD);

    caml_callback3(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(y)->data,
	    NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim void callml_vdiv(N_Vector x, N_Vector y, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVDIV);

    caml_callback3(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(y)->data,
	    NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim void callml_vscale(realtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVSCALE);

    caml_callback3(mlop, caml_copy_double(c), NVEC_CONTENT(x)->data, 
	    NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim void callml_vabs(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVABS);

    caml_callback2(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim void callml_vinv(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVINV);

    caml_callback2(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim void callml_vaddconst(N_Vector x, realtype b, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_OP(x, NVECTOR_OPS_NVADDCONST);

    caml_callback3(mlop, NVEC_CONTENT(x)->data, caml_copy_double(b),
	    NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim realtype callml_vdotprod(N_Vector x, N_Vector y)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVDOTPROD);

    r = caml_callback2(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(y)->data);

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vmaxnorm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(x, NVECTOR_OPS_NVMAXNORM);

    r = caml_callback(mlop, NVEC_CONTENT(x)->data);

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vwrmsnorm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(x, NVECTOR_OPS_NVWRMSNORM);

    r = caml_callback2(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(w)->data);

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWRMSNORMMASK);

    r = caml_callback3(mlop, NVEC_CONTENT(x)->data,
	    NVEC_CONTENT(w)->data, NVEC_CONTENT(id)->data);

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vmin(N_Vector x)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_OP(x, NVECTOR_OPS_NVMIN);

    r = caml_callback(mlop, NVEC_CONTENT(x)->data);

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vwl2norm(N_Vector x, N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVWL2NORM);

    r = caml_callback2(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(w)->data);

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim realtype callml_vl1norm(N_Vector x)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVL1NORM);

    r = caml_callback(mlop, NVEC_CONTENT(x)->data);

    CAMLreturnT(realtype, Double_val(r));
}

CAMLprim void callml_vcompare(realtype c, N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal1(mlop);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVCOMPARE);

    caml_callback3(mlop, caml_copy_double(c),
	    NVEC_CONTENT(x)->data, NVEC_CONTENT(z)->data);

    CAMLreturn0;
}

CAMLprim booleantype callml_vinvtest(N_Vector x, N_Vector z)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVINVTEST);

    r = caml_callback2(mlop, NVEC_CONTENT(x)->data, NVEC_CONTENT(z)->data);

    CAMLreturnT(booleantype, Bool_val(r));
}

CAMLprim booleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(x, NVECTOR_OPS_NVCONSTRMASK);

    r = caml_callback3(mlop, NVEC_CONTENT(c)->data,
	    NVEC_CONTENT(x)->data, NVEC_CONTENT(m)->data);

    CAMLreturnT(booleantype, Bool_val(r));
}

CAMLprim realtype callml_vminquotient(N_Vector num, N_Vector denom)
{
    CAMLparam0();
    CAMLlocal2(mlop, r);
    mlop = GET_SOME_OP(num, NVECTOR_OPS_NVMINQUOTIENT);

    r = caml_callback2(mlop, NVEC_CONTENT(num)->data,
	    NVEC_CONTENT(denom)->data);

    CAMLreturnT(realtype, Double_val(r));
}

