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
#include "nvector_parallel_ml.h"

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include <nvector/nvector_parallel.h>

/* Must correspond with camlmpi.h */
#define Comm_val(comm) (*((MPI_Comm *) &Field(comm, 1)))

/** Parallel nvectors * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
    Store_field(v_payload, 1, Field(w_payload, 1));
    Store_field(v_payload, 2, Field(w_payload, 2));
    
    v = alloc_cnvec(sizeof(struct _N_VectorContent_Parallel), v_payload);
    if (v == NULL) caml_raise_out_of_memory();
    content = (N_VectorContent_Parallel) v->content;

    /* Create vector operation structure */
    clone_cnvec_ops(v, w);

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
    CAMLparam1(payload);
    CAMLlocal2(vnvec, vlocalba);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Parallel content;
    MPI_Comm comm;
    long int n, nsum, local_length, global_length;

    vlocalba      = Field(payload, 0);
    local_length  = (Caml_ba_array_val(vlocalba))->dim[0];
    global_length = Long_val(Field(payload, 1));
    comm          = Comm_val(Field(payload, 2));

    /* Compute global length as sum of local lengths */
    n = local_length;
    MPI_Allreduce(&n, &nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
    if (nsum != global_length)
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
    content->comm          = comm;
    content->own_data      = 0;
    content->data          = Caml_ba_data_val(vlocalba);

    vnvec = caml_alloc_tuple(2);
    Store_field(vnvec, 0, payload);
    Store_field(vnvec, 1, val_cnvec(nv, finalize_cnvec));

    CAMLreturn(vnvec);
}

CAMLprim void ml_nvec_par_n_vlinearsum(value va, value vx, value vb, value vy,
				       value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_VLinearSum_Parallel(Double_val(va), NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim void ml_nvec_par_n_vconst(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst_Parallel(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim void ml_nvec_par_n_vprod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VProd_Parallel(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim void ml_nvec_par_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_VDiv_Parallel(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim void ml_nvec_par_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VScale_Parallel(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim void ml_nvec_par_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VAbs_Parallel(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim void ml_nvec_par_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_VInv_Parallel(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim void ml_nvec_par_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_VAddConst_Parallel(NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim value ml_nvec_par_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    realtype r = N_VDotProd_Parallel(NVEC_VAL(vx), NVEC_VAL(vy));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_par_n_vmaxnorm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMaxNorm_Parallel(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_par_n_vwrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    realtype r = N_VWrmsNorm_Parallel(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_par_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    realtype r = N_VWrmsNormMask_Parallel(NVEC_VAL(vx), NVEC_VAL(vw), NVEC_VAL(vid));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_par_n_vmin(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMin_Parallel(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_par_n_vwl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    realtype r = N_VWL2Norm_Parallel(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value ml_nvec_par_n_vl1norm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VL1Norm_Parallel(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim void ml_nvec_par_n_vcompare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_VCompare_Parallel(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn0;
}

CAMLprim value ml_nvec_par_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    booleantype r = N_VInvTest_Parallel(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
}

CAMLprim value ml_nvec_par_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    booleantype r = N_VConstrMask_Parallel(NVEC_VAL(vc), NVEC_VAL(vx), NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
}

CAMLprim value ml_nvec_par_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    realtype r = N_VMinQuotient_Parallel(NVEC_VAL(vnum), NVEC_VAL(vdenom));
    CAMLreturn(caml_copy_double(r));
}

