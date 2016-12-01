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
#if SUNDIALS_LIB_VERSION >= 270
    ops->nvgetvectorid	   = SUNDIALS_NVEC_PTHREADS;
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

/** Interface to underlying pthreads nvector functions */

CAMLprim value ml_nvec_pthreads_n_vlinearsum(value va, value vx, value vb,
					     value vy, value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_VLinearSum_Pthreads(Double_val(va), NVEC_VAL(vx), Double_val(vb),
			  NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vconst(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst_Pthreads(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vprod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VProd_Pthreads(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VDiv_Pthreads(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VScale_Pthreads(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VAbs_Pthreads(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VInv_Pthreads(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_VAddConst_Pthreads(NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    realtype r = N_VDotProd_Pthreads(NVEC_VAL(vx), NVEC_VAL(vy));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_pthreads_n_vmaxnorm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMaxNorm_Pthreads(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_pthreads_n_vwrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    realtype r = N_VWrmsNorm_Pthreads(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_pthreads_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    realtype r = N_VWrmsNormMask_Pthreads(NVEC_VAL(vx), NVEC_VAL(vw),
					  NVEC_VAL(vid));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_pthreads_n_vmin(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMin_Pthreads(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_pthreads_n_vwl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    realtype r = N_VWL2Norm_Pthreads(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_pthreads_n_vl1norm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VL1Norm_Pthreads(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_pthreads_n_vcompare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VCompare_Pthreads(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value ml_nvec_pthreads_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    booleantype r = N_VInvTest_Pthreads(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
}

CAMLprim value ml_nvec_pthreads_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    booleantype r = N_VConstrMask_Pthreads(NVEC_VAL(vc), NVEC_VAL(vx),
					   NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
}

CAMLprim value ml_nvec_pthreads_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    realtype r = N_VMinQuotient_Pthreads(NVEC_VAL(vnum), NVEC_VAL(vdenom));
    CAMLreturn(caml_copy_double(r));
}

