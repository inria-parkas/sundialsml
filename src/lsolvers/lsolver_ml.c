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
    SUNLinearSolver ls = LSOLVER(vls);
    if (ls) SUNLinSolFree(ls);
}

static value alloc_lsolver(SUNLinearSolver ls)
{
    CAMLparam0();
    CAMLlocal1(vcptr);

    vcptr = caml_alloc_final(1, &finalize_lsolver, 1, 20);
    LSOLVER_CPTR(vcptr) = ls;

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
	    LSOLVER_EXN(MatrixNotSquare);

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
	    LSOLVER_EXN(MatrixNotSquare);

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
	    LSOLVER_EXN(MatrixNotSquare);

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
	    LSOLVER_EXN(MatrixNotSquare);

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
	    LSOLVER_EXN(MatrixNotSquare);

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
    SUNKLUReInit(LSOLVER(vcptr), MAT_VAL(vsmat), 0, 2);
#endif
    CAMLreturn0;
}

CAMLprim void ml_lsolver_klu_set_ordering(value vcptr, value vordering)
{
    CAMLparam2(vcptr, vordering);
#if SUNDIALS_LIB_VERSION >= 300 && defined SUNDIALS_ML_KLU
    // ignore return value
    SUNKLUSetOrdering(LSOLVER(vcptr), Int_val(vordering));
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
	    LSOLVER_EXN(MatrixNotSquare);

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
    SUNSuperLUMTSetOrdering(LSOLVER(vcptr), Int_val(vordering));
#endif
    CAMLreturn0;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Iterative
 */

static int lsolver_precond_type(value vptype)
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

static int lsolver_gs_type(value vgstype)
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
	value vpretype)
{
    CAMLparam2(vcptr, vpretype);

#if SUNDIALS_LIB_VERSION >= 300
    int pretype = lsolver_precond_type(vpretype);

    // ignore returned values
    switch (Int_val(vsolver)) {
	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPFGMR:
	    SUNSPFGMRSetPrecType(LSOLVER(vcptr), pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
	    SUNSPGMRSetPrecType(LSOLVER(vcptr), pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPBCGS:
	    SUNSPBCGSSetPrecType(LSOLVER(vcptr), pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR:
	    SUNSPTFQMRSetPrecType(LSOLVER(vcptr), pretype);
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG:
	    SUNPCGSetPrecType(LSOLVER(vcptr), pretype);
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
	    SUNSPBCGSSetMaxl(LSOLVER(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPTFQMR:
	    SUNSPTFQMRSetMaxl(LSOLVER(vcptr), Int_val(vmaxl));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_PCG:
	    SUNPCGSetMaxl(LSOLVER(vcptr), Int_val(vmaxl));
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
	    SUNSPFGMRSetGSType(LSOLVER(vcptr), lsolver_gs_type(vgst));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
	    SUNSPGMRSetGSType(LSOLVER(vcptr), lsolver_gs_type(vgst));
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
	    SUNSPFGMRSetMaxRestarts(LSOLVER(vcptr), Int_val(vmaxr));
	    break;

	case VARIANT_LSOLVER_ITERATIVE_SOLVER_SPGMR:
	    SUNSPGMRSetMaxRestarts(LSOLVER(vcptr), Int_val(vmaxr));
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

