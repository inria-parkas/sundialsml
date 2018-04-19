/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2018 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include "../config.h"

#include "../sundials/sundials_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../lsolvers/lsolver_ml.h"
#include "../lsolvers/matrix_ml.h"

#if SUNDIALS_LIB_VERSION >= 300
#include <sundials/sundials_linearsolver.h>

#include <stdio.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>
#include <sundials/sundials_math.h>

#ifdef SUNDIALS_ML_KLU
#include <sunlinsol/sunlinsol_klu.h>
#endif

#ifdef SUNDIALS_ML_SUPERLUMT
#include <sunlinsol/sunlinsol_superlumt.h>
#endif

#ifdef SUNDIALS_ML_LAPACK
#include <sunlinsol/sunlinsol_lapackband.h>
#include <sunlinsol/sunlinsol_lapackdense.h>
#endif

#endif

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

CAMLprim void ml_lsolver_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (LSOLVER, exns);
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Linear Solver cptrs
 */
#if SUNDIALS_LIB_VERSION >= 300
static void finalize_lsolver(value vls)
{
    SUNLinearSolver ls = LSOLVER_VAL(vls);
    if (ls) SUNLinSolFree(ls);
}

static value alloc_lsolver(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(vcptr);

    vcptr = caml_alloc_final(1, &finalize_lsolver, 1, 20);
    LSOLVER_VAL(vcptr) = ls;

    CAMLreturn(vcptr);
}
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Direct
 */

CAMLprim value ml_lsolver_dense(value vnvec, value vdmat)
{
    CAMLparam2(vnvec, vdmat);
#if SUNDIALS_LIB_VERSION >= 300
    SUNMatrix dmat = MAT_VAL(vdmat);
    SUNLinearSolver ls = SUNDenseLinearSolver(NVEC_VAL(vnvec), dmat);

    if (ls == NULL) {
	if (SUNDenseMatrix_Rows(dmat) != SUNDenseMatrix_Columns(dmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_lapack_dense(value vnvec, value vdmat)
{
    CAMLparam2(vnvec, vdmat);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_LAPACK
    SUNMatrix dmat = MAT_VAL(vdmat);
    SUNLinearSolver ls = SUNLapackDense(NVEC_VAL(vnvec), dmat);

    if (ls == NULL) {
	if (SUNDenseMatrix_Rows(bmat) != SUNDenseMatrix_Columns(bmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_band(value vnvec, value vbmat)
{
    CAMLparam2(vnvec, vbmat);
#if SUNDIALS_LIB_VERSION >= 300
    SUNMatrix bmat = MAT_VAL(vbmat);
    SUNLinearSolver ls = SUNBandLinearSolver(NVEC_VAL(vnvec), bmat);

    if (ls == NULL) {
	if (SUNBandMatrix_Rows(bmat) != SUNBandMatrix_Columns(bmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNBandMatrix_StoredUpperBandwidth(bmat) <
	    SUNMIN(SUNBandMatrix_Rows(bmat) - 1,
		   SUNBandMatrix_LowerBandwidth(bmat)
		   + SUNBandMatrix_UpperBandwidth(bmat)))
	    caml_raise_constant(LSOLVER_EXN(InsufficientStorageUpperBandwidth));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_lapack_band(value vnvec, value vbmat)
{
    CAMLparam2(vnvec, vbmat);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_LAPACK
    SUNMatrix bmat = MAT_VAL(vbmat);
    SUNLinearSolver ls = SUNLapackBand(NVEC_VAL(vnvec), bmat);

    if (ls == NULL) {
	if (SUNBandMatrix_Rows(bmat) != SUNBandMatrix_Columns(bmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNBandMatrix_StoredUpperBandwidth(bmat) <
	    SUNMIN(SUNBandMatrix_Rows(bmat) - 1,
		   SUNBandMatrix_LowerBandwidth(bmat)
		   + SUNBandMatrix_UpperBandwidth(bmat)))
	    caml_raise_constant(LSOLVER_EXN(InsufficientStorageUpperBandwidth));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_klu(value vnvec, value vsmat)
{
    CAMLparam2(vnvec, vsmat);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_KLU
    SUNMatrix smat = MAT_VAL(vsmat);
    SUNLinearSolver ls = SUNKLU(NVEC_VAL(vnvec), smat);

    if (ls == NULL) {
	if (SUNSparseMatrix_Rows(smat) != SUNSparseMatrix_Columns(smat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim void ml_lsolver_klu_reinit(value vcptr, value vsmat)
{
    CAMLparam2(vcptr, vsmat);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_KLU
    // reinit is done at ML level; nnz arg is ignore when reinit_type = 2
    SUNKLUReInit(LSOLVER_VAL(vcptr), MAT_VAL(vsmat), 0, 2);
#endif
    CAMLreturn0;
}

CAMLprim void ml_lsolver_klu_set_ordering(value vcptr, value vordering)
{
    CAMLparam2(vcptr, vordering);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_KLU
    // ignore return value
    SUNKLUSetOrdering(LSOLVER_VAL(vcptr), Int_val(vordering));
#endif
    CAMLreturn0;
}

CAMLprim value ml_lsolver_superlumt(value vnvec, value vsmat, value vnthreads)
{
    CAMLparam2(vnvec, vsmat);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_SUPERLUMT
    SUNMatrix smat = MAT_VAL(vsmat);
    SUNLinearSolver ls = SUNSuperLUMT(NVEC_VAL(vnvec), smat,
	    Int_val(vnthreads));

    if (ls == NULL) {
	if (SUNSparseMatrix_Rows(smat) != SUNSparseMatrix_Columns(smat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim void ml_lsolver_superlumt_set_ordering(value vcptr, value vordering)
{
    CAMLparam2(vcptr, vordering);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_SUPERLUMT
    SUNSuperLUMTSetOrdering(LSOLVER_VAL(vcptr), Int_val(vordering));
#endif
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Iterative
 */

int lsolver_precond_type(value vptype)
{
    switch (Int_val(vptype)) {
    case VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_NONE:
	return PREC_NONE;

    case VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_LEFT:
	return PREC_LEFT;

    case VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_RIGHT:
	return PREC_RIGHT;

    case VARIANT_LSOLVER_PRECONDITIONING_TYPE_PREC_TYPE_BOTH:
	return PREC_BOTH;
    }

    return -1;
}

int lsolver_gs_type(value vgstype)
{
    switch (Int_val(vgstype)) {
    case VARIANT_LSOLVER_GRAMSCHMIDT_TYPE_MODIFIEDGS:
	return MODIFIED_GS;

    case VARIANT_LSOLVER_GRAMSCHMIDT_TYPE_CLASSICALGS:
	return CLASSICAL_GS;
    }

    return -1;
}

CAMLprim void ml_lsolver_set_prec_type(value vcptr, value vsolver,
	value vpretype, value vdocheck)
{
    CAMLparam4(vcptr, vsolver, vpretype, vdocheck);

#if SUNDIALS_LIB_VERSION >= 300
    int old_pretype = PREC_NONE;
    int pretype = lsolver_precond_type(vpretype);
    SUNLinearSolver lsolv = LSOLVER_VAL(vcptr);

    if (Bool_val(vdocheck)) {
	switch (Int_val(vsolver)) {
	    case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPFGMR:
		old_pretype =
		    ((SUNLinearSolverContent_SPFGMR)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
		old_pretype =
		    ((SUNLinearSolverContent_SPGMR)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPBCGS:
		old_pretype =
		    ((SUNLinearSolverContent_SPBCGS)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR:
		old_pretype =
		    ((SUNLinearSolverContent_SPTFQMR)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG:
		old_pretype =
		    ((SUNLinearSolverContent_PCG)(lsolv->content))->pretype;
		break;
	}

	if ((old_pretype == PREC_NONE) && (pretype != PREC_NONE))
	    caml_raise_constant(LSOLVER_EXN(IllegalPrecType));
    }

    // ignore returned values
    switch (Int_val(vsolver)) {
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPFGMR:
	    SUNSPFGMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
	    SUNSPGMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPBCGS:
	    SUNSPBCGSSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR:
	    SUNSPTFQMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG:
	    SUNPCGSetPrecType(lsolv, pretype);
	    break;
    }
#endif

    CAMLreturn0;
}

CAMLprim void ml_lsolver_set_maxl(value vcptr, value vsolver, value vmaxl)
{
    CAMLparam3(vcptr, vsolver, vmaxl);

#if SUNDIALS_LIB_VERSION >= 300
    switch (Int_val(vsolver)) {
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPBCGS:
	    SUNSPBCGSSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR:
	    SUNSPTFQMRSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG:
	    SUNPCGSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPFGMR:
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
	    break;
    }
#endif

    CAMLreturn0;
}

CAMLprim void ml_lsolver_set_gs_type(value vcptr, value vsolver, value vgst)
{
    CAMLparam3(vcptr, vsolver, vgst);

#if SUNDIALS_LIB_VERSION >= 300
    switch (Int_val(vsolver)) {
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPFGMR:
	    SUNSPFGMRSetGSType(LSOLVER_VAL(vcptr), lsolver_gs_type(vgst));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
	    SUNSPGMRSetGSType(LSOLVER_VAL(vcptr), lsolver_gs_type(vgst));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPBCGS:
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR:
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG:
	    break;
    }
#endif

    CAMLreturn0;
}

CAMLprim void ml_lsolver_set_max_restarts(value vcptr, value vsolver,
	value vmaxr)
{
    CAMLparam3(vcptr, vsolver, vmaxr);

#if SUNDIALS_LIB_VERSION >= 300
    switch (Int_val(vsolver)) {
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPFGMR:
	    SUNSPFGMRSetMaxRestarts(LSOLVER_VAL(vcptr), Int_val(vmaxr));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
	    SUNSPGMRSetMaxRestarts(LSOLVER_VAL(vcptr), Int_val(vmaxr));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPBCGS:
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR:
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG:
	    break;
    }
#endif

    CAMLreturn0;
}

CAMLprim value ml_lsolver_spbcgs(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if SUNDIALS_LIB_VERSION >= 300
    SUNLinearSolver ls = SUNSPBCGS(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_spfgmr(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if SUNDIALS_LIB_VERSION >= 300
    SUNLinearSolver ls = SUNSPFGMR(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_spgmr(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if SUNDIALS_LIB_VERSION >= 300
    SUNLinearSolver ls = SUNSPGMR(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_sptfqmr(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if SUNDIALS_LIB_VERSION >= 300
    SUNLinearSolver ls = SUNSPTFQMR(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value ml_lsolver_pcg(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if SUNDIALS_LIB_VERSION >= 300
    SUNLinearSolver ls = SUNPCG(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Custom
 */

#if SUNDIALS_LIB_VERSION >= 300

#define LSOLV_OP_TABLE(ls)  ((ls)->content)
#define GET_OP(ls, x) (Field((value)LSOLV_OP_TABLE(ls), RECORD_LSOLVER_OPS_ ## x))

#define CHECK_EXCEPTION(result)						 \
    (Is_exception_result (result)					 \
     ? lsolver_translate_exception (result = Extract_exception (result)) \
     : result)

#define CHECK_EXCEPTION_SUCCESS(result)					 \
    (Is_exception_result (result)					 \
     ? lsolver_translate_exception (result = Extract_exception (result)) \
     : SUNLS_SUCCESS)

static int lsolver_translate_exception(value vexn)
{
    CAMLparam1(vexn);
    CAMLlocal1(vtag);
    int r;

    vtag = Field(vexn, 0);

    if (vtag == SUNDIALS_EXN_TAG(RecoverableFailure)) {
	r = 100;

    } else if (vtag == LSOLVER_EXN_TAG(ATimesFailure)) {
	r = Bool_val(Field(vexn, 1)) ? SUNLS_ATIMES_FAIL_REC
				     : SUNLS_ATIMES_FAIL_UNREC;

    } else if (vtag == LSOLVER_EXN_TAG(PSetFailure)) {
	r = Bool_val(Field(vexn, 1)) ? SUNLS_PSET_FAIL_REC
				     : SUNLS_PSET_FAIL_UNREC;

    } else if (vtag == LSOLVER_EXN_TAG(PSolveFailure)) {
	r = Bool_val(Field(vexn, 1)) ? SUNLS_PSOLVE_FAIL_UNREC
				     : SUNLS_PSOLVE_FAIL_REC;

    } else if (vtag == LSOLVER_EXN_TAG(GSFailure)) {
	r = SUNLS_GS_FAIL;

    } else if (vtag == LSOLVER_EXN_TAG(QRSolFailure)) {
	r = SUNLS_QRSOL_FAIL;

    } else if (vtag == LSOLVER_EXN_TAG(ResReduced)) {
	r = SUNLS_RES_REDUCED;

    } else if (vtag == LSOLVER_EXN_TAG(ConvFailure)) {
	r = SUNLS_CONV_FAIL;

    } else if (vtag == LSOLVER_EXN_TAG(QRfactFailure)) {
	r = SUNLS_QRFACT_FAIL;

    } else if (vtag == LSOLVER_EXN_TAG(LUfactFailure)) {
	r = SUNLS_LUFACT_FAIL;

    } else if (vtag == LSOLVER_EXN_TAG(PackageFailure)) {
	r = Bool_val(Field(vexn, 1)) ? SUNLS_PACKAGE_FAIL_REC
				     : SUNLS_PACKAGE_FAIL_UNREC;

    } else if (vtag == LSOLVER_EXN_TAG(InvalidArgument)) {
	r = SUNLS_ILL_INPUT;

    } else {
	r = -100;
    }

    CAMLreturnT(int, r);
}

static void lsolver_raise_error(int r, const char *msg)
{
    CAMLparam0();
    CAMLlocal1(varg);

    switch (r) {
	case SUNLS_ILL_INPUT:
	    caml_invalid_argument(msg);

	case SUNLS_ATIMES_FAIL_UNREC:
	    caml_raise_with_arg(LSOLVER_EXN(ATimesFailure), Val_false);

	case SUNLS_PSET_FAIL_UNREC:
	    caml_raise_with_arg(LSOLVER_EXN(PSetFailure), Val_false);

	case SUNLS_PSOLVE_FAIL_UNREC:
	    caml_raise_with_arg(LSOLVER_EXN(PSolveFailure), Val_false);

	case SUNLS_PACKAGE_FAIL_UNREC:
	    caml_raise_with_arg(LSOLVER_EXN(PackageFailure), Val_false);

	case SUNLS_GS_FAIL:
	    caml_raise_constant(LSOLVER_EXN(GSFailure));

	case SUNLS_QRSOL_FAIL:
	    caml_raise_constant(LSOLVER_EXN(QRSolFailure));

	case SUNLS_RES_REDUCED:
	    caml_raise_constant(LSOLVER_EXN(ResReduced));

	case SUNLS_CONV_FAIL:
	    caml_raise_constant(LSOLVER_EXN(ConvFailure));

	case SUNLS_ATIMES_FAIL_REC:
	    caml_raise_with_arg(LSOLVER_EXN(ATimesFailure), Val_true);

	case SUNLS_PSET_FAIL_REC:
	    caml_raise_with_arg(LSOLVER_EXN(PSetFailure), Val_true);

	case SUNLS_PSOLVE_FAIL_REC:
	    caml_raise_with_arg(LSOLVER_EXN(PSolveFailure), Val_true);

	case SUNLS_PACKAGE_FAIL_REC:
	    caml_raise_with_arg(LSOLVER_EXN(PackageFailure), Val_true);

	case SUNLS_QRFACT_FAIL:
	    caml_raise_constant(LSOLVER_EXN(QRfactFailure));

	case SUNLS_LUFACT_FAIL:
	    caml_raise_constant(LSOLVER_EXN(LUfactFailure));

	default:
	    varg = caml_alloc_tuple(2);
	    Store_field(varg, 0, caml_copy_string(msg));
	    Store_field(varg, 1, Val_int(r));
	    caml_raise_with_arg(LSOLVER_EXN(InternalFailure), varg);
    }

    CAMLreturn0;
}

SUNLinearSolver_Type callml_custom_gettype_direct(SUNLinearSolver ls)
{
    return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_Type callml_custom_gettype_iterative(SUNLinearSolver ls)
{
    return SUNLINEARSOLVER_ITERATIVE;
}

struct atimes_with_data {
    ATimesFn atimes_func;
    void *atimes_data;
};
#define ATIMES_WITH_DATA(v) ((struct atimes_with_data *)Data_custom_val(v))

CAMLprim value ml_lsolver_call_atimes(value vcptr, value vv, value vz)
{
    CAMLparam3(vcptr, vv, vz);
    int r;

    r = ATIMES_WITH_DATA(vcptr)->atimes_func(
	    ATIMES_WITH_DATA(vcptr)->atimes_data, NVEC_VAL(vv), NVEC_VAL(vz));
    if (r != 0) lsolver_raise_error(r, "atimes");

    CAMLreturn(Val_unit);
}

int callml_custom_setatimes(SUNLinearSolver ls, void* A_data, ATimesFn ATimes)
{
    CAMLparam0();
    CAMLlocal2(vcptr, r);

    vcptr = caml_alloc_final(
		(sizeof(struct atimes_with_data) + sizeof(value) - 1)
		    / sizeof(value),
		NULL, 0, 1);
    ATIMES_WITH_DATA(vcptr)->atimes_func = ATimes;
    ATIMES_WITH_DATA(vcptr)->atimes_data = A_data;

    r = caml_callback_exn(GET_OP(ls, SET_ATIMES), vcptr);

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

struct precond_with_data {
    PSetupFn psetup_func;
    PSolveFn psolve_func;
    void *precond_data;
};
#define PRECOND_WITH_DATA(v) ((struct precond_with_data *)Data_custom_val(v))

CAMLprim value ml_lsolver_call_psetup(value vcptr)
{
    CAMLparam1(vcptr);
    int r;
    PSetupFn psetup = PRECOND_WITH_DATA(vcptr)->psetup_func;

    if (psetup != NULL) {
	r = psetup(PRECOND_WITH_DATA(vcptr)->precond_data);
	if (r != 0) lsolver_raise_error(r, "psetup");
    }

    CAMLreturn(Val_unit);
}

CAMLprim value ml_lsolver_call_psolve(value vcptr, value vr, value vz,
				      value vtol, value vlr)
{
    CAMLparam5(vcptr, vr, vz, vtol, vlr);
    int r;
    PSolveFn psolve = PRECOND_WITH_DATA(vcptr)->psolve_func;

    if (psolve != NULL) {
	r = psolve(PRECOND_WITH_DATA(vcptr)->precond_data,
		   NVEC_VAL(vr),
		   NVEC_VAL(vz),
		   Double_val(vtol),
		   Bool_val(vlr) ? 1 : 2);
	if (r != 0) lsolver_raise_error(r, "psolve");
    }

    CAMLreturn(Val_unit);
}

int callml_custom_setpreconditioner(SUNLinearSolver ls, void* P_data,
				    PSetupFn Pset, PSolveFn Psol)
{
    CAMLparam0();
    CAMLlocal2(vcptr, r);

    vcptr = caml_alloc_final(
		(sizeof(struct precond_with_data) + sizeof(value) - 1)
		    / sizeof(value),
		NULL, 0, 1);
    PRECOND_WITH_DATA(vcptr)->psetup_func = Pset;
    PRECOND_WITH_DATA(vcptr)->psolve_func = Psol;
    PRECOND_WITH_DATA(vcptr)->precond_data = P_data;

    r = caml_callback3_exn(GET_OP(ls, SET_PRECONDITIONER),
	    vcptr, Val_bool(Pset != NULL), Val_bool(Psol != NULL));

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

int callml_custom_setscalingvectors(SUNLinearSolver ls, N_Vector s1, N_Vector s2)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback2_exn(GET_OP(ls, SET_SCALING_VECTORS),
	    (s1 == NULL) ? Val_none : Some_val(NVEC_VAL(s1)),
	    (s2 == NULL) ? Val_none : Some_val(NVEC_VAL(s2)));

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

int callml_custom_initialize(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, INIT), Val_unit);

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

int callml_custom_setup(SUNLinearSolver ls, SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, SETUP),
	    (A == NULL) ? Val_unit : MAT_BACKLINK(A));

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

int callml_custom_solve(SUNLinearSolver ls, SUNMatrix A, N_Vector x,
                        N_Vector b, realtype tol)
{
    CAMLparam0();
    CAMLlocal1(r);
    CAMLlocalN(args, 4);

    Store_field(args, 0, (A == NULL) ? Val_unit : MAT_BACKLINK(A));
    Store_field(args, 1, NVEC_BACKLINK(x));
    Store_field(args, 2, NVEC_BACKLINK(b));
    Store_field(args, 3, caml_copy_double(tol));

    r = caml_callbackN_exn(GET_OP(ls, SOLVE), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

int callml_custom_numiters(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_NUM_ITERS), Val_unit);
    if (Is_exception_result (r)) {
	sundials_ml_warn_discarded_exn (Extract_exception (r),
					"user-defined num iters handler");
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, Int_val(r));
}

realtype callml_custom_resnorm(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_RES_NORM), Val_unit);
    if (Is_exception_result (r)) {
	sundials_ml_warn_discarded_exn (Extract_exception (r),
					"user-defined res norm handler");
	CAMLreturnT(realtype, 0.0);
    }

    CAMLreturnT(realtype, Double_val(r));
}

N_Vector callml_custom_resid(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_RES_ID), Val_unit);
    if (Is_exception_result (r)) {
	sundials_ml_warn_discarded_exn (Extract_exception (r),
					"user-defined res id handler");
	CAMLreturnT(N_Vector, NULL);
    }

    CAMLreturnT(N_Vector, NVEC_VAL(r));
}

int callml_custom_space(SUNLinearSolver ls, long int *lenrwLS, long int *leniwLS)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_WORK_SPACE), Val_unit);
    if (Is_exception_result (r)) {
	r = Extract_exception (r);
	lenrwLS = 0;
	leniwLS = 0;
	CAMLreturnT(int, lsolver_translate_exception(r));
    }

    *lenrwLS = Long_val(Field(r, 0));
    *leniwLS = Long_val(Field(r, 1));

    CAMLreturnT(int, SUNLS_SUCCESS);
}

int callml_custom_free(SUNLinearSolver ls)
{
    caml_remove_generational_global_root((value *)&(ls->content));
    if (ls->ops != NULL) free(ls->ops);
    free(ls);

    return(SUNLS_SUCCESS);
}

#else // SUNDIALS_LIB_VERSION < 300

CAMLprim value ml_lsolver_call_atimes(value vcptr, value vv, value vz)
{
    CAMLparam3(vcptr, vv, vz);
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    CAMLreturn(Val_unit);
}

CAMLprim value ml_lsolver_call_psetup(value vcptr)
{
    CAMLparam1(vcptr);
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    CAMLreturn(Val_unit);
}

CAMLprim value ml_lsolver_call_psolve(value vcptr, value vr, value vz,
	value vtol, value vlr)
{
    CAMLparam5(vcptr, vr, vz, vtol, vlr);
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    CAMLreturn(Val_unit);
}
#endif

CAMLprim value ml_lsolver_make_custom(value vid, value vops, value vhasops)
{
    CAMLparam3(vid, vops, vhasops);
#if SUNDIALS_LIB_VERSION >= 300
    SUNLinearSolver ls;
    SUNLinearSolver_Ops ops;

    ls = (SUNLinearSolver)malloc(sizeof *ls);
    if (ls == NULL) caml_raise_out_of_memory();

    ops = (SUNLinearSolver_Ops) malloc(
	    sizeof(struct _generic_SUNLinearSolver_Ops));
    if (ops == NULL) {
	free(ls);
	caml_raise_out_of_memory();
    }

    /* Attach operations */
    ops->gettype           = (Int_val(vid) == 0)
				? callml_custom_gettype_direct
				: callml_custom_gettype_iterative;
    ops->initialize	   = callml_custom_initialize;
    ops->setup             = callml_custom_setup;
    ops->solve             = callml_custom_solve;
    ops->lastflag          = NULL;
    ops->free              = callml_custom_free;
    ops->setatimes         =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_SET_ATIMES))
	? callml_custom_setatimes : NULL;
    ops->setpreconditioner =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_SET_PRECONDITIONER))
	? callml_custom_setpreconditioner : NULL;
    ops->setscalingvectors =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_SET_SCALING_VECTORS))
	? callml_custom_setscalingvectors : NULL;
    ops->numiters          =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_GET_NUM_ITERS))
	? callml_custom_numiters : NULL;
    ops->resnorm           =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_GET_RES_NORM))
	? callml_custom_resnorm : NULL;
    ops->resid             =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_GET_RES_ID))
	? callml_custom_resid : NULL;
    ops->space             =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_GET_WORK_SPACE))
	? callml_custom_space : NULL;

    ls->ops = ops;
    ls->content = (void *)vops;
    caml_register_generational_global_root((void *)&(ls->content));

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Low-level iterative algorithms on arrays
 */

CAMLprim value c_spils_modified_gs(value vv, value vh, value vk, value vp)
{
    CAMLparam4(vv, vh, vk, vp);

    int p = Int_val(vp);
    int k = Int_val(vk);
    int i;
    realtype new_vk_norm;
    N_Vector* v;

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < k) || (hn < k))
	caml_invalid_argument("Spils.modified_gs: h is too small.");
    if (Wosize_val (vv) < k)
	caml_invalid_argument("Spils.modified_gs: v is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_out_of_memory();
    for (i = 0; i <= p; ++i) {
	v[i] = NVEC_VAL(Field(vv, k - p + i));
    }

    ModifiedGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm);

    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value c_spils_classical_gs(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal3(vv, vh, vs);

    int k = Int_val(Field(vargs, 2));
    int p = Int_val(Field(vargs, 3));
    N_Vector temp;

    vv = Field(vargs, 0);
    vh = Field(vargs, 1);
    vs = Field(vargs, 5);

    int i;
    realtype new_vk_norm;
    N_Vector* v;

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if ((hm < k) || (hn < k))
	caml_invalid_argument("Spils.classical_gs: h is too small.");
    if (Wosize_val (vv) < k)
	caml_invalid_argument("Spils.classical_gs: v is too small.");
    if (ARRAY1_LEN(vs) < k)
	caml_invalid_argument("Spils.classical_gs: s is too small.");
#endif

    v = calloc(p + 1, sizeof(N_Vector));

    if (v == NULL)
	caml_raise_out_of_memory();
    for (i = 0; i <= p; ++i) {
	v[i] = NVEC_VAL(Field(vv, k - p + i));
    }

    temp = NVEC_VAL(Field(vargs, 4));
    ClassicalGS(v, ARRAY2_ACOLS(vh), p + 1, p, &new_vk_norm,
	        temp, REAL_ARRAY(vs));

    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value c_spils_qr_fact(value vh, value vq, value vnewjob)
{
    CAMLparam3(vh, vq, vnewjob);
    int r;
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hn = bh->dim[0];

#if SUNDIALS_ML_SAFE == 1
    intnat hm = bh->dim[1];

    if ((hm < hn + 1))
	caml_invalid_argument("Spils.qr_fact: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * hn)
	caml_invalid_argument("Spils.qr_fact: q is too small.");
#endif

    r = QRfact(hn, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), Bool_val(vnewjob));

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

CAMLprim value c_spils_qr_sol(value vh, value vq, value vb)
{
    CAMLparam3(vh, vq, vb);
    int r;
    struct caml_ba_array *bh = ARRAY2_DATA(vh);
    intnat hn = bh->dim[0];

#if SUNDIALS_ML_SAFE == 1
    intnat hm = bh->dim[1];

    if (hm < hn + 1)
	caml_invalid_argument("Spils.qr_sol: h is too small.");
    if (ARRAY1_LEN(vq) < 2 * hn)
	caml_invalid_argument("Spils.qr_sol: q is too small.");
    if (ARRAY1_LEN(vb) < hn + 1)
	caml_invalid_argument("Spils.qr_sol: b is too small.");
#endif

    r = QRsol(hn, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), REAL_ARRAY(vb));

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

