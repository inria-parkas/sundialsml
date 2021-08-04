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

#include "../sundials/sundials_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../nvectors/nvector_parallel_ml.h"

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

CAMLprim value sunml_nvector_parallel_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (NVECTOR_PARALLEL, exns);
    CAMLreturn (Val_unit);
}


/* Adapted from sundials-2.5.0/src/nvec_par/nvector_parallel.c:
   N_VCloneEmpty_Parallel */
static N_Vector clone_parallel(N_Vector w)
{
    CAMLparam0();
    CAMLlocal2(v_payload, w_payload);

    N_Vector v;
    N_VectorContent_Parallel content;

    if (w == NULL) CAMLreturnT (N_Vector, NULL);
    w_payload = NVEC_BACKLINK(w);
    struct caml_ba_array *w_ba = Caml_ba_array_val(Field(w_payload, 0));

    /* Create vector (we need not copy the data) */
    v_payload = caml_alloc_tuple(3);
    Store_field(v_payload, 0,
		caml_ba_alloc(w_ba->flags, w_ba->num_dims, NULL, w_ba->dim));
    Store_field(v_payload, 1, Field(w_payload, 1));
    Store_field(v_payload, 2, Field(w_payload, 2));
    
    v = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Parallel), v_payload);
    if (v == NULL) CAMLreturnT (N_Vector, NULL);
    content = (N_VectorContent_Parallel) v->content;

    /* Create vector operation structure */
    sunml_clone_cnvec_ops(v, w);

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
CAMLprim value sunml_nvec_wrap_parallel(value payload, value checkfn)
{
    CAMLparam2(payload, checkfn);
    CAMLlocal2(vnvec, vlocalba);

    N_Vector nv;
    N_Vector_Ops ops;
    N_VectorContent_Parallel content;
    MPI_Comm comm;
    sundials_ml_index local_length, global_length;

    vlocalba      = Field(payload, 0);
    local_length  = (Caml_ba_array_val(vlocalba))->dim[0];
    global_length = Index_val(Field(payload, 1));
    comm          = Comm_val(Field(payload, 2));

#if SUNDIALS_ML_SAFE
    {
    /* Compute global length as sum of local lengths */
    sundials_ml_index nsum;
    sundials_ml_index n = local_length;
    MPI_Allreduce(&n, &nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
    if (nsum != global_length)
        caml_raise_constant(NVECTOR_PARALLEL_EXN (IncorrectGlobalSize));
    }
#endif

    /* Create vector */
    nv = sunml_alloc_cnvec(sizeof(struct _N_VectorContent_Parallel), payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (N_VectorContent_Parallel) nv->content;

    /* Create vector operation structure */
    ops->nvclone           = clone_parallel;		    /* ours */
    ops->nvcloneempty      = NULL;
    ops->nvdestroy         = sunml_free_cnvec;
#if 270 <= SUNDIALS_LIB_VERSION
    ops->nvgetvectorid	   = N_VGetVectorID_Parallel;
#endif

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
    ops->nvgetlength	    = N_VGetLength_Parallel;
    ops->nvgetcommunicator  = N_VGetCommunicator_Parallel;

    ops->nvdotprodlocal     = N_VDotProdLocal_Parallel;
    ops->nvmaxnormlocal     = N_VMaxNormLocal_Parallel;
    ops->nvminlocal         = N_VMinLocal_Parallel;
    ops->nvl1normlocal      = N_VL1NormLocal_Parallel;
    ops->nvinvtestlocal     = N_VInvTestLocal_Parallel;
    ops->nvconstrmasklocal  = N_VConstrMaskLocal_Parallel;
    ops->nvminquotientlocal = N_VMinQuotientLocal_Parallel;
    ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Parallel;
    ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Parallel;
#endif

    /* Attach lengths and communicator */
    content->local_length  = local_length;
    content->global_length = global_length;
    content->comm          = comm;
    content->own_data      = 0;
    content->data          = Caml_ba_data_val(vlocalba);

    vnvec = caml_alloc_tuple(3);
    Store_field(vnvec, 0, payload);
    Store_field(vnvec, 1, sunml_alloc_caml_nvec(nv, sunml_finalize_caml_nvec));
    Store_field(vnvec, 2, checkfn);

    CAMLreturn(vnvec);
}

CAMLprim value sunml_nvec_par_n_vlinearsum(value va, value vx, value vb, value vy,
					value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(y) != NV_LOCLENGTH_P(x)
	    || NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vlinearsum");
#endif

    N_VLinearSum_Parallel(Double_val(va), x, Double_val(vb), y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vconst(value vc, value vz)
{
    CAMLparam2(vc, vz);
    N_VConst_Parallel(Double_val(vc), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vprod(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(y) != NV_LOCLENGTH_P(x)
	    || NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vprod");
#endif

    N_VProd_Parallel(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vdiv(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(y) != NV_LOCLENGTH_P(x)
	    || NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vdiv");
#endif

    N_VDiv_Parallel(x, y, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vscale(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vscale");
#endif

    N_VScale_Parallel(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vabs(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vabs");
#endif

    N_VAbs_Parallel(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vinv(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vinv");
#endif

    N_VInv_Parallel(x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vaddconst(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vaddconst");
#endif

    N_VAddConst_Parallel(x, Double_val(vb), z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vdotprod(value vx, value vy)
{
    CAMLparam2(vx, vy);
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(y) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vdotprod");
#endif

    realtype r = N_VDotProd_Parallel(x, y);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vmaxnorm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMaxNorm_Parallel(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vwrmsnorm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(w) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vwrmsnorm");
#endif

    realtype r = N_VWrmsNorm_Parallel(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vwrmsnormmask(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(w) != NV_LOCLENGTH_P(x)
	    || NV_LOCLENGTH_P(w) != NV_LOCLENGTH_P(id))
	caml_invalid_argument("Nvector_parallel.n_vwrmsnormmask");
#endif
    realtype r = N_VWrmsNormMask_Parallel(x, w, id);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vmin(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VMin_Parallel(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vwl2norm(value vx, value vw)
{
    CAMLparam2(vx, vw);
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(w) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vwl2norm");
#endif

    realtype r = N_VWL2Norm_Parallel(x, w);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vl1norm(value vx)
{
    CAMLparam1(vx);
    realtype r = N_VL1Norm_Parallel(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vcompare(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vcompare");
#endif
    N_VCompare_Parallel(Double_val(vc), x, z);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_n_vinvtest(value vx, value vz)
{
    CAMLparam2(vx, vz);
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vinvtest");
#endif

    booleantype r = N_VInvTest_Parallel(x, z);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_par_n_vconstrmask(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    N_Vector c = NVEC_VAL(vc);
    N_Vector x = NVEC_VAL(vx);
    N_Vector m = NVEC_VAL(vm);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(x) != NV_LOCLENGTH_P(c)
	    || NV_LOCLENGTH_P(m) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vconstrmask");
#endif

    booleantype r = N_VConstrMask_Parallel(c, x, m);
    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_par_n_vminquotient(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
    N_Vector num = NVEC_VAL(vnum);
    N_Vector denom = NVEC_VAL(vdenom);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(num) != NV_LOCLENGTH_P(denom))
	caml_invalid_argument("Nvector_parallel.n_vminquotient");
#endif

    realtype r = N_VMinQuotient_Parallel(num, denom);
    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vspace(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
    sundials_ml_index lrw, liw;

    N_VSpace_Parallel(NVEC_VAL(vx), &lrw, &liw);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_index(lrw));
    Store_field(r, 1, Val_index(liw));

    CAMLreturn(r);
}

CAMLprim value sunml_nvec_par_n_vgetlength(value vx)
{
    CAMLparam1(vx);
    CAMLlocal1(r);
#if 500 <= SUNDIALS_LIB_VERSION
    r = Val_int(N_VGetLength_Parallel(NVEC_VAL(vx)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

/* fused vector operations */

CAMLprim value sunml_nvec_par_n_vlinearcombination(value vac, value vax,
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
	    || NV_LOCLENGTH_P(ax[0]) != NV_LOCLENGTH_P(z))
	caml_invalid_argument("Nvector_parallel.n_vlinearcombination");
#endif

    N_VLinearCombination_Parallel(nvec, ac, ax, z);
    free(ax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vscaleaddmulti(value vac, value vx, value vay,
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
	    || NV_LOCLENGTH_P(a[0][0]) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vscaleaddmulti");
#endif

    N_VScaleAddMulti_Parallel(nvec, ac, x, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vdotprodmulti(value vx, value vay, value vad)
{
    CAMLparam3(vx, vay, vad);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ad = REAL_ARRAY(vad);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *ay;
    int nvec = sunml_arrays_of_nvectors(&ay, 1, vay);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vad) < nvec
		|| NV_LOCLENGTH_P(ay[0]) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vdotprodmulti");
#endif

    N_VDotProdMulti_Parallel(nvec, x, ay, ad);
    free(ay);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}


/* vector array operations */

CAMLprim value sunml_nvec_par_n_vlinearsumvectorarray(value va, value vax,
						      value vb, value vay,
						      value vaz)
{
    CAMLparam5(va, vax, vb, vay, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *a[3];
    int nvec = sunml_arrays_of_nvectors(a, 3, vax, vay, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_parallel.n_vlinearsumvectorarray");
#endif

    N_VLinearSumVectorArray_Parallel(nvec, Double_val(va), a[0],
					 Double_val(vb), a[1], a[2]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vscalevectorarray(value vac, value vax,
						  value vaz)
{
    CAMLparam3(vac, vax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(vac) < nvec)
	caml_invalid_argument("Nvector_parallel.n_vscalevectorarray");
#endif

    N_VScaleVectorArray_Parallel(nvec, ac, a[0], a[1]);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vconstvectorarray(value vc, value vaz)
{
    CAMLparam2(vc, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector *az;
    int nvec = sunml_arrays_of_nvectors(&az, 1, vaz);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec) caml_invalid_argument("Nvector_parallel.n_vconstvectorarray");
#endif

    N_VConstVectorArray_Parallel(nvec, Double_val(vc), az);
    free(az);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vwrmsnormvectorarray(value vax, value vaw,
						     value van)
{
    CAMLparam3(vax, vaw, van);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *an = REAL_ARRAY(van);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);

#if SUNDIALS_ML_SAFE == 1
    if (!nvec || ARRAY1_LEN(van) < nvec)
	caml_invalid_argument("Nvector_parallel.n_vconstvectorarray");
#endif

    N_VWrmsNormVectorArray_Parallel(nvec, a[0], a[1], an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vwrmsnormmaskvectorarray(value vax,
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
		|| NV_LOCLENGTH_P(i) != NV_LOCLENGTH_P(a[0][0]))
	caml_invalid_argument("Nvector_parallel.n_vconstvectorarray");
#endif

    N_VWrmsNormMaskVectorArray_Parallel(nvec, a[0], a[1], i, an);
    free(*a);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vscaleaddmultivectorarray(value vaa,
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
	caml_invalid_argument("Nvector_parallel.n_vscaleaddmultivectorarray");
    }
#endif

    N_VScaleAddMultiVectorArray_Parallel(nvec, nsum, aa, ax, ayz[0], ayz[1]);
    free(ax);
    free(*ayz);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_nvec_par_n_vlinearcombinationvectorarray(value vac,
							value vaax, value vaz)
{
    CAMLparam3(vac, vaax, vaz);
#if 400 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *az;
    int nvecz __attribute__((unused))
	= sunml_arrays_of_nvectors(&az, 1, vaz);
    N_Vector **aax;
    int nvec, nsum;
    sunml_arrays_of_nvectors2(&nsum, &nvec, &aax, 1, vaax);

#if SUNDIALS_ML_SAFE == 1
    if (!nvecz || !nsum || nvec > nvecz || ARRAY1_LEN(vac) < nsum) {
	if (az != NULL) free(az);
	if (aax != NULL) free(aax);
	caml_invalid_argument("Nvector_parallel.n_vlinearcombinationvectorarray");
    }
#endif

    N_VLinearCombinationVectorArray_Parallel(nvec, nsum, ac, aax, az);
    free(aax);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

/** Reduce operations for parallel nvectors */

CAMLprim value sunml_nvec_par_n_vdotprodlocal(value vx, value vy)
{
    CAMLparam2(vx, vy);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector y = NVEC_VAL(vy);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(y) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vdotprodlocal");
#endif

    r = N_VDotProdLocal_Parallel(x, y);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vmaxnormlocal(value vx)
{
    CAMLparam1(vx);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    r = N_VMaxNormLocal_Parallel(x);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vminlocal(value vx)
{
    CAMLparam1(vx);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    r = N_VMinLocal_Parallel(x);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vl1normlocal(value vx)
{
    CAMLparam1(vx);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);

    r = N_VL1NormLocal_Parallel(x);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vinvtestlocal(value vx, value vz)
{
    CAMLparam2(vx, vz);
    booleantype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector z = NVEC_VAL(vz);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(z) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vinvtestlocal");
#endif

    r = N_VInvTestLocal_Parallel(x, z);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_par_n_vconstrmasklocal(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
    booleantype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector c = NVEC_VAL(vc);
    N_Vector x = NVEC_VAL(vx);
    N_Vector m = NVEC_VAL(vm);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(x) != NV_LOCLENGTH_P(c)
	    || NV_LOCLENGTH_P(m) != NV_LOCLENGTH_P(c))
	caml_invalid_argument("Nvector_parallel.n_vconstrmasklocal");
#endif

    r = N_VConstrMaskLocal_Parallel(c, x, m);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(Val_bool(r));
}

CAMLprim value sunml_nvec_par_n_vminquotientlocal(value vn, value vd)
{
    CAMLparam2(vn, vd);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector n = NVEC_VAL(vn);
    N_Vector d = NVEC_VAL(vd);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(n) != NV_LOCLENGTH_P(d))
	caml_invalid_argument("Nvector_parallel.n_vminquotientlocal");
#endif

    r = N_VMinQuotientLocal_Parallel(n, d);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vwsqrsumlocal(value vx, value vw)
{
    CAMLparam2(vx, vw);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(w) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vwsqrsumlocal");
#endif

    r = N_VWSqrSumLocal_Parallel(x, w);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

CAMLprim value sunml_nvec_par_n_vwsqrsummasklocal(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
    realtype r;

#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector x = NVEC_VAL(vx);
    N_Vector w = NVEC_VAL(vw);
    N_Vector id = NVEC_VAL(vid);

#if SUNDIALS_ML_SAFE == 1
    if (NV_LOCLENGTH_P(w) != NV_LOCLENGTH_P(x)
	    || NV_LOCLENGTH_P(id) != NV_LOCLENGTH_P(x))
	caml_invalid_argument("Nvector_parallel.n_vwsqrsummasklocal");
#endif

    r = N_VWSqrSumMaskLocal_Parallel(x, w, id);

#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn(caml_copy_double(r));
}

/** Selectively activate fused and array operations for serial nvectors */

CAMLprim value sunml_nvec_par_enablefusedops(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableFusedOps_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablelinearcombination(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombination_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablescaleaddmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMulti_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enabledotprodmulti(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableDotProdMulti_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablelinearsumvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearSumVectorArray_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablescalevectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleVectorArray_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enableconstvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableConstVectorArray_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablewrmsnormvectorarray(value vx, value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormVectorArray_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablewrmsnormmaskvectorarray(value vx,
							       value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableWrmsNormMaskVectorArray_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablescaleaddmultivectorarray(value vx,
								value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableScaleAddMultiVectorArray_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nvec_par_enablelinearcombinationvectorarray(value vx,
								    value vv)
{
    CAMLparam2(vx, vv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_VEnableLinearCombinationVectorArray_Parallel(NVEC_VAL(vx), Bool_val(vv));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn (Val_unit);
}

