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
#include "nvector_openmp_ml.h"

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

    v = alloc_cnvec(sizeof(struct _N_VectorContent_OpenMP), v_payload);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);

    content = (N_VectorContent_OpenMP) v->content;

    /* Create vector operation structure */
    clone_cnvec_ops(v, w);

    /* Create content */
    content->length      = NV_LENGTH_OMP(w);
    content->num_threads = NV_NUM_THREADS_OMP(w);
    content->own_data    = 0;
    content->data        = Caml_ba_data_val(v_payload);

    CAMLreturnT(N_Vector, v);
}

/* Creation from OCaml.  */
/* Adapted from sundials-2.6.1/src/nvec_openmp/nvector_openmp.c:
   N_VNewEmpty_OpenMP */
CAMLprim value ml_nvec_wrap_openmp(value nthreads,
				   value payload, value checkfn)
{
    CAMLparam3(nthreads, payload, checkfn);
    CAMLlocal1(vnvec);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_OpenMP content;
    long int length = (Caml_ba_array_val(payload))->dim[0];

    /* Create vector */
    nv = alloc_cnvec(sizeof(struct _N_VectorContent_OpenMP), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_OpenMP) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_openmp;		    /* ours */
    ops->nvcloneempty      = NULL;
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = free_cnvec;

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

CAMLprim value ml_nvec_openmp_num_threads(value va)
{
    CAMLparam1(va);
    int num_threads = NV_NUM_THREADS_OMP(NVEC_VAL(va));

    CAMLreturn(Val_int(num_threads));
}

