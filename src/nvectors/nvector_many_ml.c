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

#if 500 < SUNDIALS_LIB_VERSION

/* Macro to handle separate MPI-aware/unaware installations */
#ifdef MANYVECTOR_BUILD_WITH_MPI
#include <nvector/nvector_mpimanyvector.h>
#include <nvector/nvector_mpiplusx.h>

#define MVAPPEND(fun) fun##_MPIManyVector
#define SUNML_NVEC(fun) sunml_nvec_##fun##_mpimany
#define SUNML_NVEC_OP(fun) sunml_nvec_mpimany_##fun

/* Must correspond with camlmpi.h */
#define Comm_val(comm) (*((MPI_Comm *) &Field(comm, 1)))

#else
#include <nvector/nvector_manyvector.h>

#define MVAPPEND(fun) fun##_ManyVector
#define SUNML_NVEC(fun) sunml_nvec_##fun##_many
#define SUNML_NVEC_OP(fun) sunml_nvec_many_##fun

#endif


#endif

#if 500 <= SUNDIALS_LIB_VERSION
enum do_clone_mode {
    CLONE_NORMAL,
    CLONE_ANY,
    CLONE_X_NORMAL, // MPIPLUSX
    CLONE_X_ANY,    // MPIPLUSX
};

static void free_cnvec_many(N_Vector nv)
{
    MVAPPEND(N_VectorContent) content =
	(MVAPPEND(N_VectorContent))nv->content;

    caml_remove_generational_global_root(&NVEC_BACKLINK(nv));
    if (content != NULL) {
	// We do not have to N_VDestroy the component subvectors; this will
	// happen when the OCaml payload array and its contents is garbage
	// collected.
	free(content->subvec_array);
	free(content);
    }
    free(nv->ops);
    free(nv);
}

static void finalize_caml_nvec_many (value vnv)
{
    free_cnvec_many(NVEC_CVAL (vnv));
}

/* Creation from Sundials/C.  */
/* Adapted from sundials-5.7.0/src/nvector/manyvector/nvector_manyvector.c:
   ManyVectorClone */
static N_Vector do_clone_many(N_Vector src, enum do_clone_mode clonemode)
{
    CAMLparam0();

    CAMLlocal4(dstpayload, dstarray, srcpayload, srcarray);
    CAMLlocal2(dstwrapped, srcwrapped);
#ifdef MANYVECTOR_BUILD_WITH_MPI
    CAMLlocal1(vcomm);
#endif
    static const value *pnvector_clone = NULL;
    sunindextype i;

    N_Vector dst = NULL;
    MVAPPEND(N_VectorContent) dstcontent = NULL;
    MVAPPEND(N_VectorContent) srccontent =
	(MVAPPEND(N_VectorContent))src->content;
    N_Vector *dst_subvec_array = NULL;

    if (src== NULL) CAMLreturnT(N_Vector, NULL);

    if (!pnvector_clone)
	pnvector_clone = caml_named_value("Nvector.clone");

    dst_subvec_array =
	(N_Vector *) calloc(srccontent->num_subvectors, sizeof(N_Vector));
    if (dst_subvec_array == NULL) CAMLreturnT (N_Vector, NULL);

    /* Allocate the new N_Vector and set the C-side fields. */
    dst = sunml_alloc_cnvec(sizeof *dstcontent, Val_unit);
    if (dst == NULL) {
	free(dst_subvec_array);
	CAMLreturnT (N_Vector, NULL);
    }
    dstcontent = dst->content;
    dstcontent->num_subvectors = srccontent->num_subvectors;
    dstcontent->global_length = srccontent->global_length;
    dstcontent->own_data = SUNFALSE; // wait for garbage collection
    dstcontent->subvec_array = dst_subvec_array;
    sunml_clone_cnvec_ops(dst, src);

    /* Duplicate the source array and the nvectors that it contains. */
    switch (clonemode) {
    case CLONE_ANY:
    case CLONE_X_ANY:
	srcwrapped = NVEC_BACKLINK(src);
	srcpayload = Field(srcwrapped, 1); // field 0 of extensible variant
	break;
    
    case CLONE_NORMAL:
    case CLONE_X_NORMAL:
    default:
	srcpayload = NVEC_BACKLINK(src);
	break;
    }
    srcarray = Field(srcpayload, 0);
#ifdef MANYVECTOR_BUILD_WITH_MPI
    vcomm = Field(srcpayload, 1);
#endif

    switch (clonemode) {
    case CLONE_NORMAL:
    case CLONE_ANY:
	dstarray = caml_alloc_tuple(Wosize_val(srcarray));

	/* Clone vectors into the subvector array */
	for (i=0; i < dstcontent->num_subvectors; i++) {
	    /* NB: Don't trigger GC while processing this return value!  */
	    value r = caml_callback_exn(*pnvector_clone, Field(srcarray, i));
	    if (Is_exception_result (r)) {
		sunml_free_cnvec(dst);
		free(dst_subvec_array);
		sunml_warn_discarded_exn (Extract_exception (r),
						"many vector n_vclone");
		CAMLreturnT (N_Vector, NULL);
	    }

	    dstcontent->subvec_array[i] = NVEC_VAL(r);
	    Store_field(dstarray, i, r);
	}
	break;

    case CLONE_X_NORMAL:
    case CLONE_X_ANY:
    default: {
	/* Clone vectors into the subvector array */
	value r = caml_callback_exn(*pnvector_clone, srcarray);
	if (Is_exception_result (r)) {
	    sunml_free_cnvec(dst);
	    free(dst_subvec_array);
	    sunml_warn_discarded_exn (Extract_exception (r),
					    "mpiplusx vector n_vclone");
	    CAMLreturnT (N_Vector, NULL);
	}
	dstarray = r;
	dstcontent->subvec_array[0] = NVEC_VAL(r);
	break; }
    }

    /* Set the OCaml payload of the new nvector. */
#ifdef MANYVECTOR_BUILD_WITH_MPI
    dstpayload = caml_alloc_tuple(3);
    Store_field(dstpayload, 0, dstarray);
    Store_field(dstpayload, 1, vcomm);
    Store_field(dstpayload, 2, Field(srcpayload, 2));
#else
    dstpayload = caml_alloc_tuple(2);
    Store_field(dstpayload, 0, dstarray);
    Store_field(dstpayload, 1, Field(srcpayload, 1));
#endif

    switch (clonemode) {
    case CLONE_ANY:
    case CLONE_X_ANY:
	dstwrapped = caml_alloc_tuple(2);
	Store_field(dstwrapped, 0, Field(srcwrapped, 0)); // Many constructor
	Store_field(dstwrapped, 1, dstpayload);

	NVEC_BACKLINK(dst) = dstwrapped;
	break;

    case CLONE_NORMAL:
    case CLONE_X_NORMAL:
    default:
	NVEC_BACKLINK(dst) = dstpayload;
	break;
    }

    CAMLreturnT(N_Vector, dst);
}

static N_Vector clone_many(N_Vector src)
{
    return do_clone_many(src, CLONE_NORMAL);
}

/* Clone an "any" nvector by unwrapping and wrapping the Many payload. */
static N_Vector clone_any_many(N_Vector src)
{
    return do_clone_many(src, CLONE_ANY);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
static N_Vector clone_mpiplusx(N_Vector src)
{
    return do_clone_many(src, CLONE_X_NORMAL);
}

/* Clone an "any" nvector by unwrapping and wrapping the MpiPlusX payload. */
static N_Vector clone_any_mpiplusx(N_Vector src)
{
    return do_clone_many(src, CLONE_X_ANY);
}
#endif

/*
 * N_VCloneEmpty is used in Sundials as a "light-weight" way to wrap
 * array data for use in calculations with N_Vectors. At the time of
 * writing (Sundials 5.7.0), this feature is only used in the *DenseDQJac
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
static N_Vector clone_empty_many(N_Vector w)
{
    N_Vector v;

    v = MVAPPEND(N_VCloneEmpty)(w);
    v->ops->nvdestroy = MVAPPEND(N_VDestroy);
    v->ops->nvclone = NULL;
    v->ops->nvcloneempty = NULL;

    return v;
}
#endif

/* Creation from OCaml.  */
/* Adapted from sundials-5.7.0/src/nvector/manyvector/nvector_manyvector.c:
   ManyVectorClone */
#if 500 <= SUNDIALS_LIB_VERSION
static value do_wrap(value payload,
		     value checkfn, value clonefn,
		     booleantype mpiplusx)
{
    CAMLparam3(payload, checkfn, clonefn);
    CAMLlocal3(vnvec, vnvs, vglen);
#ifdef MANYVECTOR_BUILD_WITH_MPI
    CAMLlocal1(vcomm);
#endif
    N_Vector nv;
    N_Vector_Ops ops;
    MVAPPEND(N_VectorContent) content = NULL;
    sunindextype num_subvectors, i;
    N_Vector *subvec_array;

    /* Create vector */
    nv = sunml_alloc_cnvec(sizeof *content, payload);
    if (nv == NULL) caml_raise_out_of_memory();
    ops = (N_Vector_Ops) nv->ops;
    content = (MVAPPEND(N_VectorContent)) nv->content;

    /* Create C-side array and mirror OCaml array */
#ifdef MANYVECTOR_BUILD_WITH_MPI
    vnvs = Field(payload, 0);
    vcomm = Field(payload, 1);
    vglen = Field(payload, 2);

    if (mpiplusx) {
	payload = caml_alloc_tuple(2);
	Store_field(payload, 0, vnvs);
	Store_field(payload, 1, vcomm);
    }
#else
    vnvs = Field(payload, 0);
    vglen = Field(payload, 1);
#endif
    num_subvectors = (mpiplusx ? 1 : Wosize_val(vnvs));
    subvec_array = calloc(num_subvectors, sizeof(N_Vector));
    if (subvec_array == NULL) caml_raise_out_of_memory();
    if (mpiplusx) {
	subvec_array[0] = NVEC_VAL(vnvs);
    } else {
	for (i = 0; i < num_subvectors; ++i) {
	    subvec_array[i] = NVEC_VAL(Field(vnvs, i));
	}
    }

    /* Create vector operation structure */
    ops->nvcloneempty      = clone_empty_many;		    /* ours */
    /* This is registered but only ever called for C-allocated clones. */
    ops->nvdestroy         = free_cnvec_many;		    /* ours */

    ops->nvspace           = MVAPPEND(N_VSpace);
#ifdef MANYVECTOR_BUILD_WITH_MPI
    ops->nvclone           = mpiplusx			    /* ours */
				? clone_mpiplusx
				: clone_many;
    ops->nvgetvectorid	   = mpiplusx
				? N_VGetVectorID_MPIPlusX
				: MVAPPEND(N_VGetVectorID); /* theirs */
    ops->nvgetcommunicator = N_VGetCommunicator_MPIManyVector;
    ops->nvgetarraypointer = mpiplusx
				? N_VGetArrayPointer_MPIPlusX
				: NULL;
    ops->nvsetarraypointer = mpiplusx
				? N_VSetArrayPointer_MPIPlusX
				: NULL;
#else
    ops->nvclone           = clone_many;		    /* ours */
    ops->nvgetvectorid	   = MVAPPEND(N_VGetVectorID);	    /* theirs */
    ops->nvgetcommunicator = NULL;
    ops->nvgetarraypointer = NULL;
    ops->nvsetarraypointer = NULL;
#endif
    ops->nvgetlength	   = MVAPPEND(N_VGetLength);
    ops->nvlinearsum       = MVAPPEND(N_VLinearSum);
    ops->nvconst           = MVAPPEND(N_VConst);
    ops->nvprod            = MVAPPEND(N_VProd);
    ops->nvdiv             = MVAPPEND(N_VDiv);
    ops->nvscale           = MVAPPEND(N_VScale);
    ops->nvabs             = MVAPPEND(N_VAbs);
    ops->nvinv             = MVAPPEND(N_VInv);
    ops->nvaddconst        = MVAPPEND(N_VAddConst);
    ops->nvwrmsnormmask    = MVAPPEND(N_VWrmsNormMask);
    ops->nvwrmsnorm        = MVAPPEND(N_VWrmsNorm);
    ops->nvwl2norm         = MVAPPEND(N_VWL2Norm);
    ops->nvcompare         = MVAPPEND(N_VCompare);
#ifdef MANYVECTOR_BUILD_WITH_MPI
    ops->nvdotprod         = MVAPPEND(N_VDotProd);
    ops->nvmaxnorm         = MVAPPEND(N_VMaxNorm);
    ops->nvmin             = MVAPPEND(N_VMin);
    ops->nvl1norm          = MVAPPEND(N_VL1Norm);
    ops->nvinvtest         = MVAPPEND(N_VInvTest);
    ops->nvconstrmask      = MVAPPEND(N_VConstrMask);
    ops->nvminquotient     = MVAPPEND(N_VMinQuotient);
#else
    ops->nvdotprod         = MVAPPEND(N_VDotProdLocal);
    ops->nvmaxnorm         = MVAPPEND(N_VMaxNormLocal);
    ops->nvmin             = MVAPPEND(N_VMinLocal);
    ops->nvl1norm          = MVAPPEND(N_VL1NormLocal);
    ops->nvinvtest         = MVAPPEND(N_VInvTestLocal);
    ops->nvconstrmask      = MVAPPEND(N_VConstrMaskLocal);
    ops->nvminquotient     = MVAPPEND(N_VMinQuotientLocal);
#endif
    /* fused vector operations */
    ops->nvlinearcombination = MVAPPEND(N_VLinearCombination);
    ops->nvscaleaddmulti     = MVAPPEND(N_VScaleAddMulti);
    ops->nvdotprodmulti      = MVAPPEND(N_VDotProdMulti);
    /* vector array operations */
    ops->nvlinearsumvectorarray         = NULL;
    ops->nvscalevectorarray             = NULL;
    ops->nvconstvectorarray             = NULL;
    ops->nvwrmsnormvectorarray          = MVAPPEND(N_VWrmsNormVectorArray);
    ops->nvwrmsnormmaskvectorarray      = MVAPPEND(N_VWrmsNormMaskVectorArray);
    ops->nvscaleaddmultivectorarray     = NULL;
    ops->nvlinearcombinationvectorarray = NULL;
    /* local operations */
    ops->nvdotprodlocal     = MVAPPEND(N_VDotProdLocal);
    ops->nvmaxnormlocal     = MVAPPEND(N_VMaxNormLocal);
    ops->nvminlocal         = MVAPPEND(N_VMinLocal);
    ops->nvl1normlocal      = MVAPPEND(N_VL1NormLocal);
    ops->nvinvtestlocal     = MVAPPEND(N_VInvTestLocal);
    ops->nvconstrmasklocal  = MVAPPEND(N_VConstrMaskLocal);
    ops->nvminquotientlocal = MVAPPEND(N_VMinQuotientLocal);
    ops->nvwsqrsumlocal     = MVAPPEND(N_VWSqrSumLocal);
    ops->nvwsqrsummasklocal = MVAPPEND(N_VWSqrSumMaskLocal);

    /* Create content */
    content->num_subvectors = num_subvectors;
    content->own_data = SUNFALSE;
    content->subvec_array = subvec_array;
    content->global_length = Int_val(vglen);
#ifdef MANYVECTOR_BUILD_WITH_MPI
    content->comm = Comm_val(vcomm);
#endif

    vnvec = NVEC_ALLOC();
    Store_field(vnvec, NVEC_PAYLOAD, payload);
    Store_field(vnvec, NVEC_CPTR,
		sunml_alloc_caml_nvec(nv, finalize_caml_nvec_many));
    Store_field(vnvec, NVEC_CHECK, checkfn);
    Store_field(vnvec, NVEC_CLONE, clonefn);

    CAMLreturn(vnvec);
}
#endif

CAMLprim value SUNML_NVEC(wrap)(value payload, value checkfn, value clonefn)
{
#if 500 <= SUNDIALS_LIB_VERSION
    return do_wrap(payload, checkfn, clonefn, SUNFALSE);
#else
    return Val_unit;
#endif
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
CAMLprim value sunml_nvec_wrap_mpiplusx(value payload,
					value checkfn, value clonefn)
{
#if 500 <= SUNDIALS_LIB_VERSION
    return do_wrap(payload, checkfn, clonefn, SUNTRUE);
#else
    return Val_unit;
#endif
}
#endif

/* The "any"-version of a many-vector nvector is created by modifying the
   standard one in two ways:
   1. The payload field is wrapped in the RA constructor.
   2. The nvclone operation is overridden to implement the wrapping operation
      (the current clone_empty_many does not manipulate the backlink). */
#if 500 <= SUNDIALS_LIB_VERSION
static value do_anywrap(value extconstr, value payload,
			value checkfn, value clonefn,
			booleantype mpiplusx)
{
    CAMLparam4(extconstr, payload, checkfn, clonefn);
    CAMLlocal2(vnv, vwrapped);
    N_Vector nv;
    N_Vector_Ops ops;

    vnv = do_wrap(payload, checkfn, clonefn, mpiplusx);
    nv = NVEC_VAL(vnv);
    ops = (N_Vector_Ops) nv->ops;

#ifdef MANYVECTOR_BUILD_WITH_MPI
    ops->nvclone = mpiplusx
		    ? clone_any_mpiplusx
		    : clone_any_many;
#else
    ops->nvclone = clone_any_many;
#endif

    vwrapped = caml_alloc_tuple(2);
    Store_field(vwrapped, 0, extconstr);
    Store_field(vwrapped, 1, NVEC_BACKLINK(nv));

    Store_field(vnv, 0, vwrapped);
    NVEC_BACKLINK(nv) = vwrapped;

    CAMLreturn(vnv);
}
#endif

CAMLprim value SUNML_NVEC(anywrap)(value extconstr, value payload,
			           value checkfn, value clonefn)
{
#if 500 <= SUNDIALS_LIB_VERSION
    return do_anywrap(extconstr, payload, checkfn, clonefn, SUNFALSE);
#else
    return Val_unit;
#endif
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
CAMLprim value sunml_nvec_anywrap_mpiplusx(value extconstr, value payload,
					   value checkfn, value clonefn)
{
#if 500 <= SUNDIALS_LIB_VERSION
    return do_anywrap(extconstr, payload, checkfn, clonefn, SUNTRUE);
#else
    return Val_unit;
#endif
}
#endif

/* * * * Operations * * * */

CAMLprim value SUNML_NVEC_OP(n_vlinearsum)(value va, value vx, value vb,
					   value vy, value vz)
{
    CAMLparam5(va, vx, vb, vy, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VLinearSum)(Double_val(va), NVEC_VAL(vx),
			   Double_val(vb), NVEC_VAL(vy),
			   NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vconst)(value vc, value vz)
{
    CAMLparam2(vc, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VConst)(Double_val(vc), NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vprod)(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VProd)(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vdiv)(value vx, value vy, value vz)
{
    CAMLparam3(vx, vy, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VDiv)(NVEC_VAL(vx), NVEC_VAL(vy), NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vscale)(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VScale)(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vabs)(value vx, value vz)
{
    CAMLparam2(vx, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VAbs)(NVEC_VAL(vx), NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vinv)(value vx, value vz)
{
    CAMLparam2(vx, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VInv)(NVEC_VAL(vx), NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vaddconst)(value vx, value vb, value vz)
{
    CAMLparam3(vx, vb, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VAddConst)(NVEC_VAL(vx), Double_val(vb), NVEC_VAL(vz));
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vwrmsnorm)(value vx, value vw)
{
    CAMLparam2(vx, vw);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VWrmsNorm)(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vwrmsnormmask)(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VWrmsNormMask)(NVEC_VAL(vx), NVEC_VAL(vw),
					   NVEC_VAL(vid));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
CAMLprim value SUNML_NVEC_OP(n_vdotprod)(value vx, value vy)
{
    CAMLparam2(vx, vy);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VDotProd)(NVEC_VAL(vx), NVEC_VAL(vy));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vmaxnorm)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VMaxNorm)(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vmin)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VMin)(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}
#endif

CAMLprim value SUNML_NVEC_OP(n_vwl2norm)(value vx, value vw)
{
    CAMLparam2(vx, vw);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VWL2Norm)(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
CAMLprim value SUNML_NVEC_OP(n_vl1norm)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VL1Norm)(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}
#endif

CAMLprim value SUNML_NVEC_OP(n_vcompare)(value vc, value vx, value vz)
{
    CAMLparam3(vc, vx, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    MVAPPEND(N_VCompare)(Double_val(vc), NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn (Val_unit);
#endif
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
CAMLprim value SUNML_NVEC_OP(n_vinvtest)(value vx, value vz)
{
    CAMLparam2(vx, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    booleantype r = MVAPPEND(N_VInvTest)(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
#else
    CAMLreturn (Val_unit);
#endif
}
#endif

#ifdef MANYVECTOR_BUILD_WITH_MPI
CAMLprim value SUNML_NVEC_OP(n_vconstrmask)(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
#if 500 <= SUNDIALS_LIB_VERSION
    booleantype r = MVAPPEND(N_VConstrMask)(NVEC_VAL(vc), NVEC_VAL(vx),
					    NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
#else
    CAMLreturn (Val_unit);
#endif
}
#endif

#ifdef MANYVECTOR_BUILD_WITH_MPI
CAMLprim value SUNML_NVEC_OP(n_vminquotient)(value vnum, value vdenom)
{
    CAMLparam2(vnum, vdenom);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VMinQuotient)(NVEC_VAL(vnum), NVEC_VAL(vdenom));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}
#endif

CAMLprim value SUNML_NVEC_OP(n_vspace)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(r);
    sundials_ml_index lrw, liw;

    MVAPPEND(N_VSpace)(NVEC_VAL(vx), &lrw, &liw);

    r = caml_alloc_tuple(2);
    Store_field(r, 0, Val_index(lrw));
    Store_field(r, 1, Val_index(liw));

    CAMLreturn(r);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vgetlength)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(r);
    r = Val_int(MVAPPEND(N_VGetLength)(NVEC_VAL(vx)));
    CAMLreturn(r);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vlinearcombination)(value vac, value vax,
						   value vz)
{
    CAMLparam3(vac, vax, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector z = NVEC_VAL(vz);
    N_Vector *ax;
    int nvec = sunml_arrays_of_nvectors(&ax, 1, vax);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VLinearCombination)(nvec, ac, ax, z);
    free(ax);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vscaleaddmulti)(value vac, value vx, value vay,
					       value vaz)
{
    CAMLparam4(vac, vx, vay, vaz);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vay, vaz);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VScaleAddMulti)(nvec, ac, x, a[0], a[1]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vdotprodmulti)(value vx, value vay, value vad)
{
    CAMLparam3(vx, vay, vad);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype *ad = REAL_ARRAY(vad);
    N_Vector x = NVEC_VAL(vx);
    N_Vector *ay;
    int nvec = sunml_arrays_of_nvectors(&ay, 1, vay);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VDotProdMulti)(nvec, x, ay, ad);
    free(ay);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vlinearsumvectorarray)(value va, value vax,
						      value vb, value vay,
						      value vaz)
{
    CAMLparam5(va, vax, vb, vay, vaz);
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector *a[3];
    int nvec = sunml_arrays_of_nvectors(a, 3, vax, vay, vaz);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VLinearSumVectorArray)(nvec, Double_val(va), a[0],
				      Double_val(vb), a[1], a[2]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vscalevectorarray)(value vac, value vax,
						  value vaz)
{
    CAMLparam3(vac, vax, vaz);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype *ac = REAL_ARRAY(vac);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaz);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VScaleVectorArray)(nvec, ac, a[0], a[1]);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vconstvectorarray)(value vc, value vaz)
{
    CAMLparam2(vc, vaz);
#if 500 <= SUNDIALS_LIB_VERSION
    N_Vector *az;
    int nvec = sunml_arrays_of_nvectors(&az, 1, vaz);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VConstVectorArray)(nvec, Double_val(vc), az);
    free(az);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vwrmsnormvectorarray)(value vax, value vaw,
						     value van)
{
    CAMLparam3(vax, vaw, van);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype *an = REAL_ARRAY(van);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VWrmsNormVectorArray)(nvec, a[0], a[1], an);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value SUNML_NVEC_OP(n_vwrmsnormmaskvectorarray)(value vax, value vaw,
							 value vi, value van)
{
    CAMLparam4(vax, vaw, vi, van);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype *an = REAL_ARRAY(van);
    N_Vector i = NVEC_VAL(vi);
    N_Vector *a[2];
    int nvec = sunml_arrays_of_nvectors(a, 2, vax, vaw);
    if (!nvec) caml_raise_out_of_memory();

    MVAPPEND(N_VWrmsNormMaskVectorArray)(nvec, a[0], a[1], i, an);
    free(*a);
#endif
    CAMLreturn(Val_unit);
}

/* CAMLprim value SUNML_NVEC_OP(n_vscaleaddmultivectorarray)(value vaa, value vax,
   							  value vaay, value vaaz) */

/* CAMLprim value SUNML_NVEC_OP(n_vlinearcombinationvectorarray)(value vac,
   							      value vaax,
   							      value vaz) */

CAMLprim value SUNML_NVEC_OP(n_vdotprodlocal)(value vx, value vw)
{
    CAMLparam2(vx, vw);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VDotProdLocal)(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vmaxnormlocal)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VMaxNormLocal)(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vminlocal)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VMinLocal)(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vl1normlocal)(value vx)
{
    CAMLparam1(vx);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VL1NormLocal)(NVEC_VAL(vx));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vinvtestlocal)(value vx, value vz)
{
    CAMLparam2(vx, vz);
#if 500 <= SUNDIALS_LIB_VERSION
    booleantype r = MVAPPEND(N_VInvTestLocal)(NVEC_VAL(vx), NVEC_VAL(vz));
    CAMLreturn(Val_bool(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vconstrmasklocal)(value vc, value vx, value vm)
{
    CAMLparam3(vc, vx, vm);
#if 500 <= SUNDIALS_LIB_VERSION
    booleantype r = MVAPPEND(N_VConstrMaskLocal)(NVEC_VAL(vc), NVEC_VAL(vx),
						 NVEC_VAL(vm));
    CAMLreturn(Val_bool(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vminquotientlocal)(value vn, value vd)
{
    CAMLparam2(vn, vd);
#if 500 <= SUNDIALS_LIB_VERSION
    booleantype r = MVAPPEND(N_VMinQuotientLocal)(NVEC_VAL(vn), NVEC_VAL(vd));
    CAMLreturn(Val_bool(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vwsqrsumlocal)(value vx, value vw)
{
    CAMLparam2(vx, vw);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VWSqrSumLocal)(NVEC_VAL(vx), NVEC_VAL(vw));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value SUNML_NVEC_OP(n_vwsqrsummasklocal)(value vx, value vw, value vid)
{
    CAMLparam3(vx, vw, vid);
#if 500 <= SUNDIALS_LIB_VERSION
    realtype r = MVAPPEND(N_VWSqrSumMaskLocal)(NVEC_VAL(vx), NVEC_VAL(vw),
					       NVEC_VAL(vid));
    CAMLreturn(caml_copy_double(r));
#else
    CAMLreturn (Val_unit);
#endif
}

