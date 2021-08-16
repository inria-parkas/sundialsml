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
#include "../lsolvers/sundials_linearsolver_ml.h"
#include "../lsolvers/sundials_matrix_ml.h"

#if 300 <= SUNDIALS_LIB_VERSION
#include <sundials/sundials_linearsolver.h>

#include <assert.h>

#include <stdio.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>

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

#include <nvector/nvector_serial.h>
#endif

#include <sundials/sundials_math.h>
#include <sundials/sundials_iterative.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#define MAX_ERRMSG_LEN 256

CAMLprim void sunml_lsolver_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (LSOLVER, exns);
#if 300 <= SUNDIALS_LIB_VERSION
    assert((int)VARIANT_LSOLVER_TYPE_DIRECT
	    == SUNLINEARSOLVER_DIRECT);
    assert((int)VARIANT_LSOLVER_TYPE_ITERATIVE
	    == SUNLINEARSOLVER_ITERATIVE);
#endif
#if 400 <= SUNDIALS_LIB_VERSION
    assert((int)VARIANT_LSOLVER_TYPE_ITERATIVE_MATRIX
	    == SUNLINEARSOLVER_MATRIX_ITERATIVE);
#endif
    CAMLreturn0;
}

value sunml_lsolver_exception_from_flag(int linflag)
{
    CAMLparam0();
    CAMLlocal2(vro, vr);

    if (linflag > 0) {
	vr = caml_alloc_small(1, 0);
	Field(vr, 0) = LSOLVER_EXN(ZeroInDiagonal);
	Field(vr, 1) = Val_int(linflag);
	Store_some(vro, vr);
    } else {
	switch (linflag) {
#if 300 <= SUNDIALS_LIB_VERSION
	    case SUNLS_ATIMES_FAIL_UNREC:
		vr = caml_alloc_small(1, 0);
		Field(vr, 0) = LSOLVER_EXN(ATimesFailure);
		Field(vr, 1) = Val_bool(0);
		Store_some(vro, vr);
		break;

	    case SUNLS_PSET_FAIL_UNREC:
		vr = caml_alloc_small(1, 0);
		Field(vr, 0) = LSOLVER_EXN(PSetFailure);
		Field(vr, 1) = Val_bool(0);
		Store_some(vro, vr);
		break;

	    case SUNLS_PSOLVE_FAIL_UNREC:
		vr = caml_alloc_small(1, 0);
		Field(vr, 0) = LSOLVER_EXN(PSolveFailure);
		Field(vr, 1) = Val_bool(0);
		Store_some(vro, vr);
		break;

	    case SUNLS_PACKAGE_FAIL_UNREC:
		vr = caml_alloc_small(1, 0);
		Field(vr, 0) = LSOLVER_EXN(PackageFailure);
		Field(vr, 1) = Val_bool(0);
		Store_some(vro, vr);
		break;

	    case SUNLS_GS_FAIL:
		Store_some(vro, LSOLVER_EXN(GSFailure));
		break;

	    case SUNLS_QRSOL_FAIL:
		Store_some(vro, LSOLVER_EXN(QRSolFailure));
		break;
#endif
#if 400 <= SUNDIALS_LIB_VERSION
	    case SUNLS_VECTOROP_ERR:
		Store_some(vro, LSOLVER_EXN(VectorOpError));
		break;
#endif
	    default:
		vro = Val_none;
	}
    }

    CAMLreturn(vro);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Linear Solver cptrs
 */
#if 300 <= SUNDIALS_LIB_VERSION
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

CAMLprim value sunml_lsolver_dense(value vnvec, value vdmat)
{
    CAMLparam2(vnvec, vdmat);
#if 300 <= SUNDIALS_LIB_VERSION
    SUNMatrix dmat = MAT_VAL(vdmat);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_Dense(NVEC_VAL(vnvec), dmat);
#else
    SUNLinearSolver ls = SUNDenseLinearSolver(NVEC_VAL(vnvec), dmat);
#endif

    if (ls == NULL) {
	if (SUNDenseMatrix_Rows(dmat) != SUNDenseMatrix_Columns(dmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNDenseMatrix_Rows(dmat) != NV_LENGTH_S(NVEC_VAL(vnvec)))
	    caml_raise_constant(LSOLVER_EXN(MatrixVectorMismatch));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_lapack_dense(value vnvec, value vdmat)
{
    CAMLparam2(vnvec, vdmat);
#if 300 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_LAPACK
    SUNMatrix dmat = MAT_VAL(vdmat);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_LapackDense(NVEC_VAL(vnvec), dmat);
#else
    SUNLinearSolver ls = SUNLapackDense(NVEC_VAL(vnvec), dmat);
#endif

    if (ls == NULL) {
	if (SUNDenseMatrix_Rows(dmat) != SUNDenseMatrix_Columns(dmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNDenseMatrix_Rows(dmat) != NV_LENGTH_S(NVEC_VAL(vnvec)))
	    caml_raise_constant(LSOLVER_EXN(MatrixVectorMismatch));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_band(value vnvec, value vbmat)
{
    CAMLparam2(vnvec, vbmat);
#if 300 <= SUNDIALS_LIB_VERSION
    SUNMatrix bmat = MAT_VAL(vbmat);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_Band(NVEC_VAL(vnvec), bmat);
#else
    SUNLinearSolver ls = SUNBandLinearSolver(NVEC_VAL(vnvec), bmat);
#endif

    if (ls == NULL) {
	if (SUNBandMatrix_Rows(bmat) != SUNBandMatrix_Columns(bmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNBandMatrix_StoredUpperBandwidth(bmat) <
	    SUNMIN(SUNBandMatrix_Rows(bmat) - 1,
		   SUNBandMatrix_LowerBandwidth(bmat)
		   + SUNBandMatrix_UpperBandwidth(bmat)))
	    caml_raise_constant(LSOLVER_EXN(InsufficientStorageUpperBandwidth));

	if (SUNBandMatrix_Rows(bmat) != NV_LENGTH_S(NVEC_VAL(vnvec)))
	    caml_raise_constant(LSOLVER_EXN(MatrixVectorMismatch));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_lapack_band(value vnvec, value vbmat)
{
    CAMLparam2(vnvec, vbmat);
#if 300 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_LAPACK
    SUNMatrix bmat = MAT_VAL(vbmat);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_LapackBand(NVEC_VAL(vnvec), bmat);
#else
    SUNLinearSolver ls = SUNLapackBand(NVEC_VAL(vnvec), bmat);
#endif

    if (ls == NULL) {
	if (SUNBandMatrix_Rows(bmat) != SUNBandMatrix_Columns(bmat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNBandMatrix_StoredUpperBandwidth(bmat) <
	    SUNMIN(SUNBandMatrix_Rows(bmat) - 1,
		   SUNBandMatrix_LowerBandwidth(bmat)
		   + SUNBandMatrix_UpperBandwidth(bmat)))
	    caml_raise_constant(LSOLVER_EXN(InsufficientStorageUpperBandwidth));

	if (SUNBandMatrix_Rows(bmat) != NV_LENGTH_S(NVEC_VAL(vnvec)))
	    caml_raise_constant(LSOLVER_EXN(MatrixVectorMismatch));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_klu(value vnvec, value vsmat)
{
    CAMLparam2(vnvec, vsmat);
#if 300 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_KLU
    SUNMatrix smat = MAT_VAL(vsmat);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_KLU(NVEC_VAL(vnvec), smat);
#else
    SUNLinearSolver ls = SUNKLU(NVEC_VAL(vnvec), smat);
#endif

    if (ls == NULL) {
	if (SUNSparseMatrix_Rows(smat) != SUNSparseMatrix_Columns(smat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNBandMatrix_Rows(smat) != NV_LENGTH_S(NVEC_VAL(vnvec)))
	    caml_raise_constant(LSOLVER_EXN(MatrixVectorMismatch));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim void sunml_lsolver_klu_reinit(value vcptr, value vsmat)
{
    CAMLparam2(vcptr, vsmat);
#if   400 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_KLU
    SUNLinSol_KLUReInit(LSOLVER_VAL(vcptr), MAT_VAL(vsmat),
			0, SUNKLU_REINIT_PARTIAL);
#elif 312 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_KLU
    // reinit is done at ML level; nnz arg is ignored on partial reinit
    SUNKLUReInit(LSOLVER_VAL(vcptr), MAT_VAL(vsmat), 0, SUNKLU_REINIT_PARTIAL);
#elif 300 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_KLU
    // reinit is done at ML level; nnz arg is ignored when reinit_type = 2
    SUNKLUReInit(LSOLVER_VAL(vcptr), MAT_VAL(vsmat), 0, 2);
#endif
    CAMLreturn0;
}

CAMLprim void sunml_lsolver_klu_set_ordering(value vcptr, value vordering)
{
    CAMLparam2(vcptr, vordering);
#if   400 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_KLU
    // ignore return value
    SUNLinSol_KLUSetOrdering(LSOLVER_VAL(vcptr), Int_val(vordering));
#elif 300 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_KLU
    // ignore return value
    SUNKLUSetOrdering(LSOLVER_VAL(vcptr), Int_val(vordering));
#endif
    CAMLreturn0;
}

CAMLprim value sunml_lsolver_superlumt(value vnvec, value vsmat, value vnthreads)
{
    CAMLparam2(vnvec, vsmat);
#if 300 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_SUPERLUMT
    SUNMatrix smat = MAT_VAL(vsmat);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_SuperLUMT(NVEC_VAL(vnvec), smat,
	    Int_val(vnthreads));
#else
    SUNLinearSolver ls = SUNSuperLUMT(NVEC_VAL(vnvec), smat,
	    Int_val(vnthreads));
#endif

    if (ls == NULL) {
	if (SUNSparseMatrix_Rows(smat) != SUNSparseMatrix_Columns(smat))
	    caml_raise_constant(LSOLVER_EXN(MatrixNotSquare));

	if (SUNBandMatrix_Rows(smat) != NV_LENGTH_S(NVEC_VAL(vnvec)))
	    caml_raise_constant(LSOLVER_EXN(MatrixVectorMismatch));

	caml_raise_out_of_memory();
    }

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim void sunml_lsolver_superlumt_set_ordering(value vcptr, value vordering)
{
    CAMLparam2(vcptr, vordering);
#if   400 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_SUPERLUMT
    SUNLinSol_SuperLUMTSetOrdering(LSOLVER_VAL(vcptr), Int_val(vordering));
#elif 300 <= SUNDIALS_LIB_VERSION && defined SUNDIALS_ML_SUPERLUMT
    SUNSuperLUMTSetOrdering(LSOLVER_VAL(vcptr), Int_val(vordering));
#endif
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Iterative
 */

int sunml_lsolver_precond_type(value vptype)
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

int sunml_lsolver_gs_type(value vgstype)
{
    switch (Int_val(vgstype)) {
    case VARIANT_LSOLVER_GRAMSCHMIDT_TYPE_MODIFIEDGS:
	return MODIFIED_GS;

    case VARIANT_LSOLVER_GRAMSCHMIDT_TYPE_CLASSICALGS:
	return CLASSICAL_GS;
    }

    return -1;
}

CAMLprim void sunml_lsolver_set_prec_type(value vcptr, value vsolver,
	value vpretype, value vdocheck)
{
    CAMLparam4(vcptr, vsolver, vpretype, vdocheck);

#if 300 <= SUNDIALS_LIB_VERSION
    const char* interrmsg = "internal error in sunml_lsolver_set_prec_type";
    int old_pretype = PREC_NONE;
    int pretype = sunml_lsolver_precond_type(vpretype);
    SUNLinearSolver lsolv = LSOLVER_VAL(vcptr);

    if (Bool_val(vdocheck)) {
	switch (Int_val(vsolver)) {
	    case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
		old_pretype =
		    ((SUNLinearSolverContent_SPFGMR)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
		old_pretype =
		    ((SUNLinearSolverContent_SPGMR)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
		old_pretype =
		    ((SUNLinearSolverContent_SPBCGS)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
		old_pretype =
		    ((SUNLinearSolverContent_SPTFQMR)(lsolv->content))->pretype;
		break;

	    case VARIANT_LSOLVER_SOLVER_DATA_PCG:
		old_pretype =
		    ((SUNLinearSolverContent_PCG)(lsolv->content))->pretype;
		break;

	    default:
		caml_failwith(interrmsg);
	}

	if ((old_pretype == PREC_NONE) && (pretype != PREC_NONE))
	    caml_raise_constant(LSOLVER_EXN(IllegalPrecType));
    }

    // ignore returned values
    switch (Int_val(vsolver)) {
#if 400 <= SUNDIALS_LIB_VERSION
	case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	    SUNLinSol_SPFGMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	    SUNLinSol_SPGMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	    SUNLinSol_SPBCGSSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	    SUNLinSol_SPTFQMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	    SUNLinSol_PCGSetPrecType(lsolv, pretype);
	    break;
#else
	case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	    SUNSPFGMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	    SUNSPGMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	    SUNSPBCGSSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	    SUNSPTFQMRSetPrecType(lsolv, pretype);
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	    SUNPCGSetPrecType(lsolv, pretype);
	    break;
#endif
	default:
	    caml_failwith(interrmsg);
    }
#endif

    CAMLreturn0;
}

CAMLprim void sunml_lsolver_set_maxl(value vcptr, value vsolver, value vmaxl)
{
    CAMLparam3(vcptr, vsolver, vmaxl);

#if 300 <= SUNDIALS_LIB_VERSION
    switch (Int_val(vsolver)) {
#if 400 <= SUNDIALS_LIB_VERSION
	case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	    SUNLinSol_SPBCGSSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	    SUNLinSol_SPTFQMRSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	    SUNLinSol_PCGSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;
#else
	case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	    SUNSPBCGSSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	    SUNSPTFQMRSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	    SUNPCGSetMaxl(LSOLVER_VAL(vcptr), Int_val(vmaxl));
	    break;
#endif

	case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	    break;

	default:
	    caml_failwith("internal error in sunml_lsolver_set_maxl");
    }
#endif

    CAMLreturn0;
}

CAMLprim void sunml_lsolver_set_gs_type(value vcptr, value vsolver, value vgst)
{
    CAMLparam3(vcptr, vsolver, vgst);

#if 300 <= SUNDIALS_LIB_VERSION
    switch (Int_val(vsolver)) {
#if 400 <= SUNDIALS_LIB_VERSION
	case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	    SUNLinSol_SPFGMRSetGSType(LSOLVER_VAL(vcptr),
				      sunml_lsolver_gs_type(vgst));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	    SUNLinSol_SPGMRSetGSType(LSOLVER_VAL(vcptr),
				     sunml_lsolver_gs_type(vgst));
	    break;
#else
	case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	    SUNSPFGMRSetGSType(LSOLVER_VAL(vcptr), sunml_lsolver_gs_type(vgst));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	    SUNSPGMRSetGSType(LSOLVER_VAL(vcptr), sunml_lsolver_gs_type(vgst));
	    break;
#endif

	case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	    break;
	default:
	    caml_failwith("internal error in sunml_lsolver_set_gs_type");
    }
#endif

    CAMLreturn0;
}

CAMLprim void sunml_lsolver_set_max_restarts(value vcptr, value vsolver,
	value vmaxr)
{
    CAMLparam3(vcptr, vsolver, vmaxr);

#if 300 <= SUNDIALS_LIB_VERSION
    switch (Int_val(vsolver)) {
#if 400 <= SUNDIALS_LIB_VERSION
	case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	    SUNLinSol_SPFGMRSetMaxRestarts(LSOLVER_VAL(vcptr), Int_val(vmaxr));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	    SUNLinSol_SPGMRSetMaxRestarts(LSOLVER_VAL(vcptr), Int_val(vmaxr));
	    break;
#else
	case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	    SUNSPFGMRSetMaxRestarts(LSOLVER_VAL(vcptr), Int_val(vmaxr));
	    break;

	case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	    SUNSPGMRSetMaxRestarts(LSOLVER_VAL(vcptr), Int_val(vmaxr));
	    break;
#endif

	case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	    break;
	default:
	    caml_failwith("internal error in sunml_lsolver_set_max_restarts");
    }
#endif

    CAMLreturn0;
}

CAMLprim value sunml_lsolver_spbcgs(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if 300 <= SUNDIALS_LIB_VERSION
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_SPBCGS(NVEC_VAL(vnvec),
					  PREC_NONE, Int_val(vmaxl));
#else
    SUNLinearSolver ls = SUNSPBCGS(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
#endif
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_spfgmr(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if 300 <= SUNDIALS_LIB_VERSION
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_SPFGMR(NVEC_VAL(vnvec),
					  PREC_NONE, Int_val(vmaxl));
#else
    SUNLinearSolver ls = SUNSPFGMR(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
#endif
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_spgmr(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if 300 <= SUNDIALS_LIB_VERSION
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_SPGMR(NVEC_VAL(vnvec),
					 PREC_NONE, Int_val(vmaxl));
#else
    SUNLinearSolver ls = SUNSPGMR(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
#endif
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_sptfqmr(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if 300 <= SUNDIALS_LIB_VERSION
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_SPTFQMR(NVEC_VAL(vnvec), PREC_NONE,
					   Int_val(vmaxl));
#else
    SUNLinearSolver ls = SUNSPTFQMR(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
#endif
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

CAMLprim value sunml_lsolver_pcg(value vmaxl, value vnvec)
{
    CAMLparam2(vmaxl, vnvec);
#if 300 <= SUNDIALS_LIB_VERSION
#if 400 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls = SUNLinSol_PCG(NVEC_VAL(vnvec),
				       PREC_NONE, Int_val(vmaxl));
#else
    SUNLinearSolver ls = SUNPCG(NVEC_VAL(vnvec), PREC_NONE, Int_val(vmaxl));
#endif
    if (ls == NULL) caml_raise_out_of_memory();

    CAMLreturn(alloc_lsolver(ls));
#else
    CAMLreturn(Val_unit);
#endif
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Custom
 */

#if 300 <= SUNDIALS_LIB_VERSION

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

    } else if (vtag == LSOLVER_EXN_TAG(VectorOpError)) {
#if 400 <= SUNDIALS_LIB_VERSION
	r = SUNLS_VECTOROP_ERR;
#else
	r = -100;
#endif

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

static SUNLinearSolver_Type callml_custom_gettype_direct(SUNLinearSolver ls)
{
    return SUNLINEARSOLVER_DIRECT;
}

static SUNLinearSolver_Type callml_custom_gettype_iterative(SUNLinearSolver ls)
{
    return SUNLINEARSOLVER_ITERATIVE;
}

#if 400 <= SUNDIALS_LIB_VERSION
static SUNLinearSolver_Type callml_custom_gettype_matrix_iterative(
							    SUNLinearSolver ls)
{
    return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}
#endif

struct atimes_with_data {
    ATimesFn atimes_func;
    void *atimes_data;
};
#define ATIMES_WITH_DATA(v) ((struct atimes_with_data *)Data_custom_val(v))

CAMLprim value sunml_lsolver_call_atimes(value vcptr, value vv, value vz)
{
    CAMLparam3(vcptr, vv, vz);
    int r;

    r = ATIMES_WITH_DATA(vcptr)->atimes_func(
	    ATIMES_WITH_DATA(vcptr)->atimes_data, NVEC_VAL(vv), NVEC_VAL(vz));
    if (r != 0)
	caml_raise_with_arg(LSOLVER_EXN(ATimesFailure), Val_bool(r > 0));

    CAMLreturn(Val_unit);
}

static int callml_custom_setatimes(SUNLinearSolver ls, void* A_data,
				   ATimesFn ATimes)
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

CAMLprim value sunml_lsolver_call_psetup(value vcptr)
{
    CAMLparam1(vcptr);
    int r;
    PSetupFn psetup = PRECOND_WITH_DATA(vcptr)->psetup_func;

    if (psetup != NULL) {
	r = psetup(PRECOND_WITH_DATA(vcptr)->precond_data);
	if (r != 0)
	    caml_raise_with_arg(LSOLVER_EXN(PSetFailure), Val_bool(r > 0));
    }

    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_call_psolve(value vcptr, value vr, value vz,
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
	if (r != 0)
	    caml_raise_with_arg(LSOLVER_EXN(PSolveFailure), Val_bool(r > 0));
    }

    CAMLreturn(Val_unit);
}

static int callml_custom_setpreconditioner(SUNLinearSolver ls, void* P_data,
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

static int callml_custom_setscalingvectors(SUNLinearSolver ls,
					   N_Vector s1, N_Vector s2)
{
    CAMLparam0();
    CAMLlocal3(r, ss1, ss2);

    ss1 = Val_none;
    if (s1 != NULL) Store_some(ss1, NVEC_BACKLINK(s1));
    ss2 = Val_none;
    if (s2 != NULL) Store_some(ss2, NVEC_BACKLINK(s2));
    r = caml_callback2_exn(GET_OP(ls, SET_SCALING_VECTORS), ss1, ss2);

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

static int callml_custom_initialize(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, INIT), Val_unit);

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

static int callml_custom_setup(SUNLinearSolver ls, SUNMatrix A)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, SETUP),
	    (A == NULL) ? Val_unit : MAT_BACKLINK(A));

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

static int callml_custom_solve(SUNLinearSolver ls, SUNMatrix A, N_Vector x,
                               N_Vector b, realtype tol)
{
    CAMLparam0();
    CAMLlocal1(r);
    CAMLlocalN(args, 4);

    args[0] = (A == NULL) ? Val_unit : MAT_BACKLINK(A);
    args[1] = NVEC_BACKLINK(x);
    args[2] = NVEC_BACKLINK(b);
    args[3] = caml_copy_double(tol);

    r = caml_callbackN_exn(GET_OP(ls, SOLVE), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION_SUCCESS(r));
}

#if 500 <= SUNDIALS_LIB_VERSION
static SUNLinearSolver_ID callml_custom_get_id(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);
    SUNLinearSolver_ID id;

    r = caml_callback_exn(GET_OP(ls, GET_ID), Val_unit);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined get id handler");
	id = SUNLINEARSOLVER_CUSTOM;
    } else {
	switch (Int_val(r)) {
	case VARIANT_LSOLVER_ID_BAND:
	    id = SUNLINEARSOLVER_BAND;
	    break;

	case VARIANT_LSOLVER_ID_DENSE:
	    id = SUNLINEARSOLVER_DENSE;
	    break;

	case VARIANT_LSOLVER_ID_KLU:
	    id = SUNLINEARSOLVER_KLU;
	    break;

	case VARIANT_LSOLVER_ID_LAPACKBAND:
	    id = SUNLINEARSOLVER_LAPACKBAND;
	    break;

	case VARIANT_LSOLVER_ID_LAPACKDENSE:
	    id = SUNLINEARSOLVER_LAPACKDENSE;
	    break;

	case VARIANT_LSOLVER_ID_PCG:
	    id = SUNLINEARSOLVER_PCG;
	    break;

	case VARIANT_LSOLVER_ID_SPBCGS:
	    id = SUNLINEARSOLVER_SPBCGS;
	    break;

	case VARIANT_LSOLVER_ID_SPFGMR:
	    id = SUNLINEARSOLVER_SPFGMR;
	    break;

	case VARIANT_LSOLVER_ID_SPGMR:
	    id = SUNLINEARSOLVER_SPGMR;
	    break;

	case VARIANT_LSOLVER_ID_SPTFQMR:
	    id = SUNLINEARSOLVER_SPTFQMR;
	    break;

	case VARIANT_LSOLVER_ID_SUPERLUDIST:
	    id = SUNLINEARSOLVER_SUPERLUDIST;
	    break;

	case VARIANT_LSOLVER_ID_SUPERLUMT:
	    id = SUNLINEARSOLVER_SUPERLUMT;
	    break;

	case VARIANT_LSOLVER_ID_CUSTOM:
	default:
	    id = SUNLINEARSOLVER_CUSTOM;
	    break;
	}
    }

    CAMLreturnT(SUNLinearSolver_ID, id);
}
#endif

static int callml_custom_numiters(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_NUM_ITERS), Val_unit);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined num iters handler");
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, Int_val(r));
}

static realtype callml_custom_resnorm(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_RES_NORM), Val_unit);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined res norm handler");
	CAMLreturnT(realtype, 0.0);
    }

    CAMLreturnT(realtype, Double_val(r));
}

static N_Vector callml_custom_resid(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_RES_ID), Val_unit);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined res id handler");
	CAMLreturnT(N_Vector, NULL);
    }

    CAMLreturnT(N_Vector, NVEC_VAL(r));
}

#if 500 <= SUNDIALS_LIB_VERSION
static sunindextype callml_custom_lastflag(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(r);

    r = caml_callback_exn(GET_OP(ls, GET_LAST_FLAG), Val_unit);
    if (Is_exception_result (r)) {
	sunml_warn_discarded_exn (Extract_exception (r),
					"user-defined last flag handler");
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, Int_val(r));
}
#endif

static int callml_custom_space(SUNLinearSolver ls,
			       long int *lenrwLS, long int *leniwLS)
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

static int callml_custom_free(SUNLinearSolver ls)
{
    caml_remove_generational_global_root((value *)&(ls->content));
    if (ls->ops != NULL) free(ls->ops);
    free(ls);

    return(SUNLS_SUCCESS);
}

#else // SUNDIALS_LIB_VERSION < 300

CAMLprim value sunml_lsolver_call_atimes(value vcptr, value vv, value vz)
{
    CAMLparam3(vcptr, vv, vz);
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_call_psetup(value vcptr)
{
    CAMLparam1(vcptr);
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_call_psolve(value vcptr, value vr, value vz,
	value vtol, value vlr)
{
    CAMLparam5(vcptr, vr, vz, vtol, vlr);
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
    CAMLreturn(Val_unit);
}
#endif

CAMLprim value sunml_lsolver_make_custom(value vlstype, value vlsid,
					 value vops, value vhasops)
{
    CAMLparam4(vlstype, vlsid, vops, vhasops);
#if 300 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver ls;
    SUNLinearSolver_Ops ops;
    SUNLinearSolver_Type (*callml_get_type)(SUNLinearSolver ls) = NULL;

    ls = (SUNLinearSolver)malloc(sizeof *ls);
    if (ls == NULL) caml_raise_out_of_memory();

    ops = (SUNLinearSolver_Ops) malloc(
	    sizeof(struct _generic_SUNLinearSolver_Ops));
    if (ops == NULL) {
	free(ls);
	caml_raise_out_of_memory();
    }

    switch (Int_val(vlstype)) {
	case SUNLINEARSOLVER_DIRECT:
	    callml_get_type = callml_custom_gettype_direct;
	    break;
	case SUNLINEARSOLVER_ITERATIVE:
	    callml_get_type = callml_custom_gettype_iterative;
	    break;
#if 400 <= SUNDIALS_LIB_VERSION
	case SUNLINEARSOLVER_MATRIX_ITERATIVE:
	    callml_get_type = callml_custom_gettype_matrix_iterative;
	    break;
#endif
	default:
	    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
	    break;
    }

    /* Attach operations */
    ops->gettype           = callml_get_type;

#if 500 <= SUNDIALS_LIB_VERSION
    ops->getid		   = callml_custom_get_id;
    ops->initialize        =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_INIT))
	? callml_custom_initialize : NULL;
    ops->setup             =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_SETUP))
	? callml_custom_setup : NULL;
#else
    ops->initialize	   = callml_custom_initialize;
    ops->setup             = callml_custom_setup;
#endif

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
#if 500 <= SUNDIALS_LIB_VERSION
    ops->lastflag          =
	Bool_val(Field(vhasops, RECORD_LSOLVER_HASOPS_GET_LAST_FLAG))
	? callml_custom_lastflag : NULL;
#endif
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

CAMLprim value sunml_spils_modified_gs(value vv, value vh, value vk, value vp)
{
    CAMLparam4(vv, vh, vk, vp);

    int p = Int_val(vp);
    int k = Int_val(vk);
    int i;
    int i0 = SUNMAX(k-p, 0);
    realtype new_vk_norm;
    N_Vector* v;

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_BA(vh);
    intnat hn = bh->dim[0];
    intnat hm = bh->dim[1];

    if (hn < k + 1)
	caml_invalid_argument("modified_gs: h is too small (dim1 < k + 1).");
    if (hm < k)
	caml_invalid_argument("modified_gs: h is too small (dim2 < k).");
    if (Wosize_val (vv) < k + 1)
	caml_invalid_argument("modified_gs: v is too small (< k + 1).");
#endif

    v = calloc(k + 1, sizeof(N_Vector));

    if (v == NULL) caml_raise_out_of_memory();

    for (i = i0; i <= k; ++i)
	v[i] = NVEC_VAL(Field(vv, i));

    ModifiedGS(v, ARRAY2_ACOLS(vh), k, p, &new_vk_norm);

    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value sunml_spils_classical_gs(value vargs)
{
    CAMLparam1(vargs);
    CAMLlocal4(vv, vh, vs, vtemp);

    int k = Int_val(Field(vargs, 2));
    int p = Int_val(Field(vargs, 3));
    N_Vector *temp;
    int i;
    int i0 = SUNMAX(k-p, 0);
    realtype new_vk_norm;
    N_Vector* v;

    vv = Field(vargs, 0);
    vh = Field(vargs, 1);
    vs = Field(vargs, 4);
    vtemp = Field(vargs, 5);

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_BA(vh);
    intnat hn = bh->dim[0];
    intnat hm = bh->dim[1];

    if (hn < k)
	caml_invalid_argument("classical_gs: h is too small (< k).");
    if (hm < k)
	caml_invalid_argument("classical_gs: h is too small (< k).");

    if (ARRAY1_LEN(vs) < k + 1)
	caml_invalid_argument("classical_gs: s is too small (< k + 1).");
#endif

    temp = sunml_nvector_array_alloc(vtemp);
    if (temp == NULL) caml_raise_out_of_memory();

    v = calloc(k + 1, sizeof(N_Vector));
    if (v == NULL) {
	sunml_nvector_array_free(temp);
	caml_raise_out_of_memory();
    }

    for (i = i0; i <= k; ++i)
	v[i] = NVEC_VAL(Field(vv, i));

#if 400 <= SUNDIALS_LIB_VERSION
    ClassicalGS(v, ARRAY2_ACOLS(vh), k, p, &new_vk_norm,
	        REAL_ARRAY(vs), temp);
#else
    ClassicalGS(v, ARRAY2_ACOLS(vh), k, p, &new_vk_norm,
	        temp[0], REAL_ARRAY(vs));
#endif

    sunml_nvector_array_free(temp);
    free(v);

    CAMLreturn(caml_copy_double(new_vk_norm));
}

CAMLprim value sunml_spils_qr_fact(value vn, value vh, value vq, value vnewjob)
{
    CAMLparam4(vn, vh, vq, vnewjob);
    int r;
    int n = Int_val(vn);

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_BA(vh);
    intnat hn = bh->dim[0];
    intnat hm = bh->dim[1];

    if (hn < n + 1)
	caml_invalid_argument("qr_fact: h is too small (< n + 1).");
    if (hm < n)
	caml_invalid_argument("qr_fact: h is too small (< n).");
    if (ARRAY1_LEN(vq) < 2 * n)
	caml_invalid_argument("qr_fact: q is too small (< 2n).");
#endif

    r = QRfact(n, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), Bool_val(vnewjob));

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_spils_qr_sol(value vn, value vh, value vq, value vb)
{
    CAMLparam4(vn, vh, vq, vb);
    int r;
    int n = Int_val(vn);

#if SUNDIALS_ML_SAFE == 1
    struct caml_ba_array *bh = ARRAY2_BA(vh);
    intnat hm = bh->dim[1];
    intnat hn = bh->dim[0];

    if (hn < n + 1)
	caml_invalid_argument("qr_sol: h is too small (< n + 1).");
    if (hm < n)
	caml_invalid_argument("qr_sol: h is too small (< n).");
    if (ARRAY1_LEN(vq) < 2 * n)
	caml_invalid_argument("qr_sol: q is too small (< 2n).");
    if (ARRAY1_LEN(vb) < n + 1)
	caml_invalid_argument("qr_sol: b is too small (< n + 1).");
#endif

    r = QRsol(n, ARRAY2_ACOLS(vh), REAL_ARRAY(vq), REAL_ARRAY(vb));

    if (r != 0) {
	caml_raise_with_arg(MATRIX_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }

    CAMLreturn (Val_unit);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Invoking linear solver functions from OCaml
 */

#if 300 <= SUNDIALS_LIB_VERSION
static int ocaml_atimes(void *callback_croot, N_Vector v, N_Vector z)
{
    CAMLparam0();
    CAMLlocalN (args, 2);

    value *croot = VPTRCROOT(callback_croot);

    args[0] = NVEC_BACKLINK(v);
    args[1] = NVEC_BACKLINK(z);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn
	(Field(*croot, RECORD_LSOLVER_OCAML_CALLBACKS_ATIMES), 2, args);

    CAMLreturnT(int, CHECK_EXCEPTION(r));
}

static int ocaml_psetup(void *callback_croot)
{
    CAMLparam0();

    value *croot = VPTRCROOT(callback_croot);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn
	(Field(*croot, RECORD_LSOLVER_OCAML_CALLBACKS_PSETUP), Val_unit);

    CAMLreturnT(int, CHECK_EXCEPTION(r));
}

int ocaml_psolve(void *callback_croot,
		 N_Vector w, N_Vector z, realtype tol, int lr)
{
    CAMLparam0();
    CAMLlocalN (args, 4);

    value *croot = VPTRCROOT(callback_croot);

    args[0] = NVEC_BACKLINK(w);
    args[1] = NVEC_BACKLINK(z);
    args[2] = caml_copy_double(tol);
    args[3] = (lr == 1) ? Val_true : Val_false;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn
	(Field(*croot, RECORD_LSOLVER_OCAML_CALLBACKS_PSOLVE), 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION(r));
}
#endif

#if 300 <= SUNDIALS_LIB_VERSION
static void sunml_lsolver_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == SUNLS_SUCCESS) return;

    switch (flag) {
#if 400 <= SUNDIALS_LIB_VERSION
	case SUNLS_ILL_INPUT:
	    caml_invalid_argument(call);

	case SUNLS_MEM_FAIL:
	    caml_raise_out_of_memory();

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

	case SUNLS_VECTOROP_ERR:
	    caml_raise_constant(LSOLVER_EXN(VectorOpError));

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
#endif
	default:
	    if (flag > 0) {
		caml_raise_constant(SUNDIALS_EXN(RecoverableFailure));
	    } else {
		snprintf(exmsg, MAX_ERRMSG_LEN, "%s: unexpected error", call);
		caml_failwith(exmsg);
	    }
    }
}

#define CHECK_FLAG(call, flag) if (flag != SUNLS_SUCCESS) \
				 sunml_lsolver_check_flag(call, flag)

#endif

CAMLprim value sunml_lsolver_get_type(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(r);
#if 300 <= SUNDIALS_LIB_VERSION
    r = Val_int(SUNLinSolGetType(LSOLVER_VAL(vcptr)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_lsolver_get_id(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(r);
#if 500 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver_ID id = SUNLinSolGetID(LSOLVER_VAL(vcptr));
    enum lsolver_linear_solver_id_tag result;
    static char exmsg[MAX_ERRMSG_LEN] = "";

    switch(id) {
	case SUNLINEARSOLVER_BAND:
            result = VARIANT_LSOLVER_ID_BAND;
            break;

	case SUNLINEARSOLVER_DENSE:
            result = VARIANT_LSOLVER_ID_DENSE;
            break;

	case SUNLINEARSOLVER_KLU:
            result = VARIANT_LSOLVER_ID_KLU;
            break;

	case SUNLINEARSOLVER_LAPACKBAND:
            result = VARIANT_LSOLVER_ID_LAPACKBAND;
            break;

	case SUNLINEARSOLVER_LAPACKDENSE:
            result = VARIANT_LSOLVER_ID_LAPACKDENSE;
            break;

	case SUNLINEARSOLVER_PCG:
            result = VARIANT_LSOLVER_ID_PCG;
            break;

	case SUNLINEARSOLVER_SPBCGS:
            result = VARIANT_LSOLVER_ID_SPBCGS;
            break;

	case SUNLINEARSOLVER_SPFGMR:
            result = VARIANT_LSOLVER_ID_SPFGMR;
            break;

	case SUNLINEARSOLVER_SPGMR:
            result = VARIANT_LSOLVER_ID_SPGMR;
            break;

	case SUNLINEARSOLVER_SPTFQMR:
            result = VARIANT_LSOLVER_ID_SPTFQMR;
            break;

	case SUNLINEARSOLVER_SUPERLUDIST:
            result = VARIANT_LSOLVER_ID_SUPERLUDIST;
            break;

	case SUNLINEARSOLVER_SUPERLUMT:
            result = VARIANT_LSOLVER_ID_SUPERLUMT;
            break;

	case SUNLINEARSOLVER_CUSTOM:
            result = VARIANT_LSOLVER_ID_CUSTOM;
            break;

	case SUNLINEARSOLVER_CUSOLVERSP_BATCHQR:
	    snprintf(exmsg, MAX_ERRMSG_LEN,
		"get_id: SUNLINEARSOLVER_CUSOLVERSP_BATCHQR is not supported");
	    caml_failwith(exmsg);

	default:
	    snprintf(exmsg, MAX_ERRMSG_LEN, "get_id: unknown id (%d)", id);
	    caml_failwith(exmsg);
    }

    r = Val_int(result);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_lsolver_set_atimes(value vcptr, value vcroot)
{
    CAMLparam2(vcroot, vcroot);
#if 300 <= SUNDIALS_LIB_VERSION
    int flag = SUNLinSolSetATimes(LSOLVER_VAL(vcptr), VPTRCROOT(vcroot),
				  ocaml_atimes);
    CHECK_FLAG("SUNLinSolSetATimes", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_set_preconditioner(value vcptr, value vcroot)
{
    CAMLparam2(vcptr, vcroot);
#if 300 <= SUNDIALS_LIB_VERSION
    int flag = SUNLinSolSetPreconditioner(LSOLVER_VAL(vcptr), VPTRCROOT(vcroot),
					  ocaml_psetup, ocaml_psolve);
    CHECK_FLAG("SUNLinSolSetATimes", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_set_scaling_vectors(value vcptr,
						 value vs1, value vs2)
{
    CAMLparam3(vcptr, vs1, vs2);
#if 300 <= SUNDIALS_LIB_VERSION
    int flag = SUNLinSolSetScalingVectors(LSOLVER_VAL(vcptr),
					  NVEC_VAL(vs1), NVEC_VAL(vs2));
    CHECK_FLAG("SUNLinSolSetScalingVectors", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_initialize(value vcptr)
{
    CAMLparam1(vcptr);
#if 300 <= SUNDIALS_LIB_VERSION
    int flag = SUNLinSolInitialize(LSOLVER_VAL(vcptr));
    CHECK_FLAG("SUNLinSolInitialize", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_setup(value vcptr, value vm)
{
    CAMLparam1(vcptr);
#if 300 <= SUNDIALS_LIB_VERSION
    int flag = SUNLinSolSetup(LSOLVER_VAL(vcptr), MAT_VAL(vm));
    CHECK_FLAG("SUNLinSolSetup", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_solve(value vcptr, value va, value vx,
				   value vb, value vtol)
{
    CAMLparam5(vcptr, va, vx, vb, vtol);
#if 300 <= SUNDIALS_LIB_VERSION
    int flag = SUNLinSolSolve(LSOLVER_VAL(vcptr), MAT_VAL(va), NVEC_VAL(vx),
			      NVEC_VAL(vb), Double_val(vtol));
    CHECK_FLAG("SUNLinSolSolve", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(Val_unit);
}

CAMLprim value sunml_lsolver_iters(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(r);
#if 300 <= SUNDIALS_LIB_VERSION
    r = Val_int(SUNLinSolNumIters(LSOLVER_VAL(vcptr)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_lsolver_res_norm(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(r);
#if 300 <= SUNDIALS_LIB_VERSION
    r = caml_copy_double(SUNLinSolResNorm(LSOLVER_VAL(vcptr)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_lsolver_res_id(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(r);
#if 300 <= SUNDIALS_LIB_VERSION
    N_Vector resid = SUNLinSolResid(LSOLVER_VAL(vcptr));
    r = NVEC_BACKLINK(resid);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_lsolver_last_flag(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(r);
#if 500 <= SUNDIALS_LIB_VERSION
    r = Val_int(SUNLinSolLastFlag(LSOLVER_VAL(vcptr)));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim value sunml_lsolver_space(value vcptr)
{
    CAMLparam1(vcptr);
    CAMLlocal1(r);
#if 300 <= SUNDIALS_LIB_VERSION
    long int lenrwLS, leniwLS;
    int flag = SUNLinSolSpace(LSOLVER_VAL(vcptr), &lenrwLS, &leniwLS);
    CHECK_FLAG("SUNLinSolSpace", flag);

    r = caml_alloc_tuple(3);
    Store_field(r, 0, Val_int(lenrwLS));
    Store_field(r, 1, Val_int(leniwLS));
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(r);
}

CAMLprim void sunml_lsolver_set_info_file(value vcptr, value vsolver,
					  value vfile)
{
    CAMLparam3(vcptr, vsolver, vfile);
#if 530 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver lsolv = LSOLVER_VAL(vcptr);
    const char* interrmsg = "internal error in sunml_lsolver_set_info_file";
    FILE *file = ML_CFILE(vfile);
    int flag;

    switch (Int_val(vsolver)) {
    case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	flag = SUNLinSolSetInfoFile_SPFGMR(lsolv, file);
	CHECK_FLAG("SUNLinSolSetInfofile_SPFGMR", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	flag = SUNLinSolSetInfoFile_SPGMR(lsolv, file);
	CHECK_FLAG("SUNLinSolSetInfofile_SPGMR", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	flag = SUNLinSolSetInfoFile_SPBCGS(lsolv, file);
	CHECK_FLAG("SUNLinSolSetInfofile_SPBCGS", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	flag = SUNLinSolSetInfoFile_SPTFQMR(lsolv, file);
	CHECK_FLAG("SUNLinSolSetInfofile_SPTFQMR", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	flag = SUNLinSolSetInfoFile_PCG(lsolv, file);
	CHECK_FLAG("SUNLinSolSetInfofile_PCG", flag);
	break;

    default:
	caml_failwith(interrmsg);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_lsolver_set_print_level(value vcptr, value vsolver,
					    value vlevel)
{
    CAMLparam3(vcptr, vsolver, vlevel);
#if 530 <= SUNDIALS_LIB_VERSION
    SUNLinearSolver lsolv = LSOLVER_VAL(vcptr);
    const char* interrmsg = "internal error in sunml_lsolver_set_print_level";
    int level = Int_val(vlevel);
    int flag;

    switch (Int_val(vsolver)) {
    case VARIANT_LSOLVER_SOLVER_DATA_SPFGMR:
	flag = SUNLinSolSetPrintLevel_SPFGMR(lsolv, level);
	CHECK_FLAG("SUNLinSolSetInfofile_SPFGMR", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_SPGMR:
	flag = SUNLinSolSetPrintLevel_SPGMR(lsolv, level);
	CHECK_FLAG("SUNLinSolSetInfofile_SPGMR", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_SPBCGS:
	flag = SUNLinSolSetPrintLevel_SPBCGS(lsolv, level);
	CHECK_FLAG("SUNLinSolSetInfofile_SPBCGS", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_SPTFQMR:
	flag = SUNLinSolSetPrintLevel_SPTFQMR(lsolv, level);
	CHECK_FLAG("SUNLinSolSetInfofile_SPTFQMR", flag);
	break;

    case VARIANT_LSOLVER_SOLVER_DATA_PCG:
	flag = SUNLinSolSetPrintLevel_PCG(lsolv, level);
	CHECK_FLAG("SUNLinSolSetInfofile_PCG", flag);
	break;

    default:
	caml_failwith(interrmsg);
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

