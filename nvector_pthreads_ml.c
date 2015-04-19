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

#include "sundials_ml.h"
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

    v = alloc_cnvec(sizeof(struct _N_VectorContent_Pthreads), v_payload);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_Pthreads) v->content;

    /* Create vector operation structure */
    clone_cnvec_ops(v, w);

    /* Create content */
    content->length      = NV_LENGTH_PT(w);
    content->num_threads = NV_NUM_THREADS_PT(w);
    content->own_data    = 0;
    content->data        = Caml_ba_data_val(v_payload);

    CAMLreturnT(N_Vector, v);
}

/* Creation from OCaml.  */
/* Adapted from sundials-2.6.1/src/nvec_pthreads/nvector_pthreads.c:
   N_VNewEmpty_Pthreads */
CAMLprim value ml_nvec_wrap_pthreads(value nthreads,
				     value payload, value checkfn)
{
    CAMLparam3(nthreads, payload, checkfn);
    CAMLlocal1(vnvec);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Pthreads content;
    long int length = (Caml_ba_array_val(payload))->dim[0];

    /* Create vector */
    nv = alloc_cnvec(sizeof(struct _N_VectorContent_Pthreads), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_Pthreads) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_pthreads;		    /* ours */
    ops->nvcloneempty      = NULL;
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = free_cnvec;

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

    /* Create content */
    content->length      = length;
    content->num_threads = Int_val(nthreads);
    content->own_data    = 0;
    content->data        = Caml_ba_data_val(payload);

    vnvec = caml_alloc_tuple(3);
    Store_field(vnvec, 0, payload);
    Store_field(vnvec, 1, alloc_caml_nvec(nv, finalize_caml_nvec));
    Store_field(vnvec, 2, checkfn);

    CAMLreturn(vnvec);
}

CAMLprim value ml_nvec_pthreads_num_threads(value va)
{
    CAMLparam1(va);
    int num_threads = NV_NUM_THREADS_PT(NVEC_VAL(va));

    CAMLreturn(Val_int(num_threads));
}

