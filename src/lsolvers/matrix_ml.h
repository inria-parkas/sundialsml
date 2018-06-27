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

#ifndef _MATRIX_ML_H__
#define _MATRIX_ML_H__

#include <caml/mlvalues.h>

#if SUNDIALS_LIB_VERSION >= 300
#include <sundials/sundials_matrix.h>
#else
#include <sundials/sundials_dense.h>
#include <sundials/sundials_sparse.h>
#endif

/*  OCaml matrix interface : Sundials >= 3.0.0

    The implementation works over three layers:

	 OCaml			 C
	-------			---
    1.  Bigarray		realtype* / sunindextype*

	The underlying storage is managed by Bigarrays, it is freed when
	the corresponding bigarray is collected.

    2.	matrix_content		SUNMatrixContent_Dense/_Band/_Sparse

	The matrix_content has three fields:

	- payload: refers to the underlying bigarrays (layer 1), which prevent
		   them being collected while the matrix_content exists.

	- rawptr:  a custom value holding a pointer to the associated
		   SUNMatrixContent_* object. Its finalizer unregisters
	           the associated global root and frees the allocated memory
		   (but not the array memory which is freed by layer 1 as
		    soon as there are no other references from OCaml -- since
		    the global root no longer exists).

        - valid:   this field is unused for Sundials >= 3.0.0.

	There is no need to extend the SUNMatrixContent_* records with a
	"backlink" because such values are never passed to callbacks
	(a SUNMatrix is passed).

	SUNMatrixContent_* records are never manipulated directly by Sundials,
	which always works through SUNMatrix objects. This means, in
	particular, that there is always a corresponding matrix_content on
	the OCaml side (possibly held through a Matrix.t).

    3.	Matrix.t		SUNMatrix

	The Matrix.t has four fields:

	- payload: refers to the underlying matrix_content (layer 2), which
		   prevent them being collected while the Matrix.t exists.

	- rawptr:  a custom value holding a pointer to the associated
		   SUNMatrix object. Its finalizer unregisters the associated
		   global root and frees the allocated memory (but not the
		   content memory which is freed by layer 2 as soon as there
		   are no other references from OCaml -- since the global root
		   no longer exists).

	- id:	   this field is not used from C.

	- mat_ops: for custom nvectors (only), this field is used to implement
		   callbacks.
	
	The SUNMatrix record is extended with a "backlink" global
	root to the payload value using the technique described in
	nvectors/nvector_ml.h.

	For "created" Matrix.t's (see nvector/nvector_ml.h), there is always
	a corresponding Matrix.t on the OCaml side. The SUNMatrix is freed
	when the associated Matrix.t is collected and the "backlink" global
	root is unregistered at this time.

	For "cloned" SUNMatrix's there is no Matrix.t on the OCaml side, but
	there is a matrix_content since it is required for callbacks. The
	"backlink" global root is unregistered when the SUNMatrix destroy
	operation is invoked.

    Custom Matrix.t's
    -----------------
    The payload is the value being wrapped. The content field is set to point
    to a Value containing the OCaml callback table, it is also registered as a
    global root. The user must ensure the callbacks do not hold references to
    any particular matrix, for otherwise that matrix is never reclaimed.

    Sundials creation/cloning of SUNMatrix's
    ----------------------------------------
    Care must be taken with the SUNBandMatrix/SUNDenseMatrix/SUNSparseMatrix
    creation functions since, unlike the clone operation, they do not create
    objects that respect the above design (there is no backlink and the clone
    and destroy operations are not overridden with the Sundials/ML versions).

    SUNBandMatrix is called from the *_bandpre and *_bbdpre preconditioners,
    but the resulting SUNMatrix's are only used internally and never passed
    via callbacks into OCaml or 'mixed' in operations with SUNMatrix's created
    by Sundials/ML. They thus do not concern the interface code.

    SUNSparseMatrix is called from SUNSparseFromDenseMatrix,
    SUNSparseFromBandMatrix, and SUNKLUReInit. These latter three functions
    are not used within Sundials, we reimplement them to function properly
    with OCaml.

    SUNBandMatrix is also called from SUNMatScaleAdd_Band (via the internal
    function SMScaleAddNew_Band) to reallocate the storage within the A matrix
    if necessary. SUNSparseMatrix is also called from the SUNMatScaleAdd_Sparse
    and SUNMatScaleAddI_Sparse functions for similar reasons. Sundials/ML
    overrides these operations to ensure the correct construction of OCaml
    compatible SUNMatrix values.
 */


/*  OCaml matrix interface : Sundials <= 2.7.0

    Sundials/ML aims to work with older versions of the Sundials library
    (which may still be installed by certain package managers, notably on
     Debian stable, which we take as a reference and which evolves slowly).
    This requires working around the differences between Sundials versions
    2.x and 3.x.

    The basic idea is
    
    (i) to implement Matrix.t entirely from the OCaml side, since SUNMatrix
	does not exist in Sundials <= 2.7.0;

    (ii) to use matrix_content as a substitute for Dls.Dense.t, Dls.Band.t,
         and Sls.t, with payload holding the OCaml values, and rawptr being
	 a DlsMat or SlsMat as in earlier versions of the interface;

	 NB: The payload field of matrix_content is unused, since we cannot
	 know whether or not the underlyig array has been reallocated.
	 
    (iii) not to control array storage using bigarrays, but instead to use
	  the valid field to indicate that the underlying storage has been
	  freed (we have no choice anyway since Sundials controls the
	  lifetimes directly and there is no "operator" indirection as with
	  SUNMatrix's and Nvector's).

    In other words, layer 3 is "emulated" in OCaml only, and layers 2 and 1
    are implemented exactly as in the Dls and Sls modules. Note that the
    types in these modules cannot be used directly because that would induce
    an incompatability at the OCaml level (Dls.Dense.t versus Matrix.Dense.t).
    The underlying C structures are different: DlsMat for Dls.Dense.t and
    SUNLinearSolverContent_Dense for Matrix.Dense.t.
 */

enum mat_matrix_index {
    RECORD_MAT_MATRIX_PAYLOAD   = 0,
    RECORD_MAT_MATRIX_RAWPTR,
    RECORD_MAT_MATRIX_ID,
    RECORD_MAT_MATRIX_MATOPS,
    RECORD_MAT_MATRIX_SIZE /* This has to come last. */
};

#if SUNDIALS_LIB_VERSION >= 300

/* Map a matrix_content.rawptr to a (void *) to a MAT_CONTENT_*_TYPE */
#define MAT_CONTENT(v) (*(void **)Data_custom_val(v))

/* Map a matrix_content.rawptr to a pointer to a MAT_CONTENT_*_TYPE */
#define MAT_CONTENT_DENSE(v)  (*(SUNMatrixContent_Dense *)Data_custom_val(v))
#define MAT_CONTENT_BAND(v)   (*(SUNMatrixContent_Band *)Data_custom_val(v))
#define MAT_CONTENT_SPARSE(v) (*(SUNMatrixContent_Sparse *)Data_custom_val(v))

#define MAT_CONTENT_DENSE_TYPE SUNMatrixContent_Dense
#define MAT_CONTENT_BAND_TYPE SUNMatrixContent_Band
#define MAT_CONTENT_SPARSE_TYPE SUNMatrixContent_Sparse

struct csmat {
    struct _generic_SUNMatrix smat;
    value backlink;
};

// Return the OCaml version of the sunmatrix payload
#define MAT_BACKLINK(smat) (((struct csmat *)smat)->backlink)

/* Map a Matrix.rawptr to a SUNMatrix */
#define MAT_CVAL(v) (*(SUNMatrix *)Data_custom_val(v))
// MAT_VAL turns an OCaml Matrix.t into a c-sunmatrix
#define MAT_VAL(v) (MAT_CVAL(Field(v, RECORD_MAT_MATRIX_RAWPTR)))

#else // SUNDIALS_LIB_VERSION < 300
#define DLSMAT(v) (*(DlsMat *)Data_custom_val(v))
#define SLSMAT(v) (*(SlsMat *)Data_custom_val(v))

#define MAT_CONTENT_DENSE(v)  DLSMAT(v)
#define MAT_CONTENT_BAND(v)   DLSMAT(v)
#define MAT_CONTENT_SPARSE(v) SLSMAT(v)
#define MAT_CONTENT_DENSE_TYPE DlsMat
#define MAT_CONTENT_BAND_TYPE DlsMat
#define MAT_CONTENT_SPARSE_TYPE SlsMat

CAMLprim value c_matrix_dense_wrap(DlsMat a);
CAMLprim value c_matrix_band_wrap(DlsMat a);

CAMLprim value c_matrix_sparse_wrap(SlsMat a);
CAMLprim value ml_matrix_sparse_rewrap(value vm);

#endif

// Convert Sundials CSC_MAT=0 and CSR_MAT=1 constants into
// Matrix.Sparse.sformat values
#define MAT_TO_SFORMAT(x) (Val_int(x))
#define MAT_FROM_SFORMAT(x) (Int_val(x))

enum mat_matrix_id_tag {
    MATRIX_ID_DENSE = 0,
    MATRIX_ID_BAND,
    MATRIX_ID_SPARSE,
    MATRIX_ID_CUSTOM,
};

enum mat_matrix_content_index {
    RECORD_MAT_MATRIXCONTENT_PAYLOAD   = 0,
    RECORD_MAT_MATRIXCONTENT_RAWPTR,
    RECORD_MAT_MATRIXCONTENT_VALID,
    RECORD_MAT_MATRIXCONTENT_SIZE /* This has to come last. */
};

enum mat_matrix_ops_index {
    RECORD_MAT_MATRIXOPS_CLONE   = 0,
    RECORD_MAT_MATRIXOPS_ZERO,
    RECORD_MAT_MATRIXOPS_COPY,
    RECORD_MAT_MATRIXOPS_SCALE_ADD,
    RECORD_MAT_MATRIXOPS_SCALE_ADDI,
    RECORD_MAT_MATRIXOPS_MATVEC,
    RECORD_MAT_MATRIXOPS_SPACE,
    RECORD_MAT_MATRIXOPS_SIZE /* This has to come last. */
};

enum mat_band_data_index {
    RECORD_MAT_BANDDATA_DATA = 0,
    RECORD_MAT_BANDDATA_DIMS,
    RECORD_MAT_BANDDATA_SIZE /* This has to come last. */
};

enum mat_band_dimensions_index {
    RECORD_MAT_BANDDIMENSIONS_N   = 0,
    RECORD_MAT_BANDDIMENSIONS_MU,
    RECORD_MAT_BANDDIMENSIONS_SMU,
    RECORD_MAT_BANDDIMENSIONS_ML,
    RECORD_MAT_BANDDIMENSIONS_SIZE /* This has to come last. */
};

enum mat_sparse_data_index {
    RECORD_MAT_SPARSEDATA_IDXVALS   = 0,
    RECORD_MAT_SPARSEDATA_IDXPTRS,
    RECORD_MAT_SPARSEDATA_DATA,
    RECORD_MAT_SPARSEDATA_SFORMAT,
    RECORD_MAT_SPARSEDATA_SIZE /* This has to come last. */
};

#define MAT_UNWRAP(v) (Field(v, RECORD_MAT_MATRIX_PAYLOAD))

/* This enum must list exceptions in the same order as the call to
 * c_init_module in matrix.ml.  */
enum mat_exn_index {
    MATRIX_EXN_Invalidated = 0,
    MATRIX_EXN_IncompatibleArguments,
    MATRIX_EXN_ZeroDiagonalElement,
    MATRIX_EXN_SET_SIZE
};

#define MATRIX_EXN(name)     REGISTERED_EXN(MATRIX, name)
#define MATRIX_EXN_TAG(name) REGISTERED_EXN_TAG(MATRIX, name)

#endif

