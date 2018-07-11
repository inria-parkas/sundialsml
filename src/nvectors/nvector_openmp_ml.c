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
    nv = alloc_cnvec(sizeof(struct _N_VectorContent_OpenMP), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_OpenMP) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_openmp;		    /* ours */
    ops->nvcloneempty      = NULL;
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = free_cnvec;
#if SUNDIALS_LIB_VERSION >= 270
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

CAMLprim value sunml_nvec_openmp_num_threads(value va)
{
    CAMLparam1(va);
    int num_threads = NV_NUM_THREADS_OMP(NVEC_VAL(va));

    CAMLreturn(Val_int(num_threads));
}

/** Interface to underlying openmp nvector functions */

CAMLprim value sunml_nvec_openmp_n_vlinearsum(value va, value vx, value vb, value vy,
				       value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_VLinearSum_OpenMP(Double_val(va), NVEC_VAL(vx), Double_val(vb),
			NVEC_VAL(vy), NVEC_VAL(vz));
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
    N_VProd_OpenMP(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VDiv_OpenMP(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VScale_OpenMP(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VAbs_OpenMP(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VInv_OpenMP(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_VAddConst_OpenMP(NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    realtype r = N_VDotProd_OpenMP(NVEC_VAL(vx), NVEC_VAL(vy));
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
    realtype r = N_VWrmsNorm_OpenMP(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_openmp_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    realtype r = N_VWrmsNormMask_OpenMP(NVEC_VAL(vx), NVEC_VAL(vw),
					NVEC_VAL(vid));
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
    realtype r = N_VWL2Norm_OpenMP(NVEC_VAL(vx), NVEC_VAL(vw));
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
    N_VCompare_OpenMP(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_openmp_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    booleantype r = N_VInvTest_OpenMP(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_openmp_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    booleantype r = N_VConstrMask_OpenMP(NVEC_VAL(vc), NVEC_VAL(vx),
					 NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_openmp_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    realtype r = N_VMinQuotient_OpenMP(NVEC_VAL(vnum), NVEC_VAL(vdenom));
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

