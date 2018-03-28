(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Generic matrices.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNMatrix module
    @since 3.0.0 *)

(** Generic operations that all matrix types must implement. Failure is
    signalled by raising an exception.

    @nocvode <node> Description of the SUNMatrix module *)
type ('m, 'd, 'k) matrix_ops = {
  m_clone     : 'm -> 'm;
  (** Create a new, distinct matrix from an existing one without
      copying the contents of the original matrix. *)

  m_zero      : 'm -> unit;
  (** Set all elements to zero. *)

  m_copy      : 'm -> 'm -> unit;
  (** [m_copy a b] copies the contents of [a] into [b]. *)

  m_scale_add  : float -> 'm -> 'm -> unit;
  (** [m_scale_add c a b] calculates $a = ca + b$. *)

  m_scale_addi : float -> 'm -> unit;
  (** [m_scale_addi c a] calculates $a = ca + I$. *)

  m_matvec    : 'm -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> unit;
  (** [m_matvec a x y] calculates $y = ax$. *)

  m_space     : 'm -> int * int;
  (** [lrw, liw = m_space a] returns the number of realtype words [lrw] and
      integer words [liw] required to store the matrix [a]. *)
}

(** Raised in {!Sundials.sundials_version} < 3.0.0 on an attempt to
    access a value that has become invalid. Such values refer to matrices
    that no longer exist in the underlying library. Values never
    become invalid in {!Sundials.sundials_version} >= 3.0.0. *)
exception Invalidated

(** Raised if matrix operation arguments are mutually incompatible. *)
exception IncompatibleArguments

(** Dense matrices *)
module Dense :
  sig
    (* TODO: check documentation links *)
    (** A dense matrix. Values of this type are typically passed to linear
        solver callback functions (like {!Cvode.dense_jac_fn}, 
        {!Ida.dense_jac_fn}, and {!Kinsol.dense_jac_fn}).

        @nocvode <node> The SUNMatrix_Dense implementation *)
    type t

    (** {4 Basic access} *)

    (** [make m n x] returns an [m] by [n] dense matrix with elements set
        to [x].

        @nocvode <node> SUNDenseMatrix *)
    val make : int -> int -> float -> t

    (** [create m n] returns an uninitialized [m] by [n] dense matrix.

        @nocvode <node> SUNDenseMatrix *)
    val create : int -> int -> t

    (** [m, n = size a] returns the numbers of columns [m] and rows [n]
        of [a].

        @nocvode <node> SM_ROWS_D
        @nocvode <node> SM_COLUMNS_D
        @nocvode <node> SUNDenseMatrix_Rows
        @nocvode <node> SUNDenseMatrix_Columns *)
    val size : t -> int * int

    (* TOPLEVEL-PRINTER: Matrix.Dense.pp *)
    (** Pretty-print a dense matrix using the
        {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
    val pp : Format.formatter -> t -> unit

    (** Pretty-print a dense matrix using the
        {{:OCAML_DOC_ROOT(Format.html)} Format} module.
        The defaults are: [start="\["], [stop="\]"], [sep=";"],
        [indent=4], [itemsep=" "], and
        [item=fun f r c->Format.fprintf f "(%2d,%2d)=% -15e" r c] (see
        {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}).
        The [indent] argument specifies the indent for wrapped rows. *)
    val ppi : ?start:string -> ?stop:string -> ?sep:string
              -> ?indent:int -> ?itemsep:string
              -> ?item:(Format.formatter -> int -> int -> float -> unit)
              -> Format.formatter -> t -> unit

    (** [get a i j] returns the value at row [i] and column [j] of [a].

        @nocvode <node> SM_ELEMENT_D *)
    val get : t -> int -> int -> float

    (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v].

        @nocvode <node> SM_ELEMENT_D *)
    val set : t -> int -> int -> float -> unit

    (** [update a i j f] sets the value at row [i] and column [j] of [a]
        to [f v].

        @nocvode <node> SM_ELEMENT_D *)
    val update : t -> int -> int -> (float -> float) -> unit

    (** Direct access to the underlying storage array, which is accessed
        column first (unlike in {!get}).

        NB: For {!Sundials.sundials_version} < 3.0.0, this access is
        potentially unsafe and {b must} only be used when the underlying
        storage is valid, which will be the case in callbacks.

        @nocvode <node> SM_CONTENT_D *)
    val unwrap : t -> Sundials.real_array2

    (** {4 Operations} *)

    (** Operations on dense matrices. *)
    val ops : (t, Nvector_serial.data, [>Nvector_serial.kind] as 'k) matrix_ops

    (** [scale_add c A B] calculates $A = cA + B$.

        @nocvode <node> SUNMatScaleAdd
        @nocvode <node> SUNMatScaleAdd_Dense *)
    val scale_add : float -> t -> t -> unit

    (** [scale_addi c A] calculates $A = cA + I$.

        @nocvode <node> SUNMatScaleAddI
        @nocvode <node> SUNMatScaleAddI_Dense *)
    val scale_addi : float -> t -> unit

    (** Compute the matrix-vector product $y = Ax$.
    
        @nocvode <node> SUNMatMatvec
        @nocvode <node> SUNMatMatvec_Dense *)
    val matvec :
      t -> x:'nk Nvector_serial.any -> y:'nk Nvector_serial.any -> unit

    (** Fills a matrix with zeros.

        @nocvode <node> SUNMatZero
        @nocvode <node> SUNMatZero_Dense *)
    val set_to_zero : t -> unit

    (** [blit src dst] copies the contents of [src] into [dst]. Both
        must have the same size.

        @nocvode <node> SUNMatCopy
        @nocvode <node> SUNMatCopy_Dense *)
    val blit : t -> t -> unit

    (** [lrw, liw = space a] returns the storage requirements of [a] as
        [lrw] realtype words and [liw] integer words.
    
        @nocvode <node> SUNMatSpace
        @nocvode <node> SUNMatSpace_Dense *)
    val space : t -> int * int

    (** {4 Low-level details} *)

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. Has no effect when
        {!Sundials.sundials_version} >= 3.0.0. *)
    val invalidate : t -> unit
  end

(** Banded matrices *)
module Band :
  sig
    (* TODO: check documentation links *)
    (** A band matrix. Values of this type are typically passed to linear
        solver callback functions (like {!Cvode.dense_jac_fn}, 
        {!Ida.dense_jac_fn}, and {!Kinsol.dense_jac_fn}).

        @nocvode <node> The SUNMatrix_Band implementation *)
    type t

    (** Band matrix dimensions. If the result will not be LU factored then
        {% $\mathtt{smu} = \mathtt{mu}$ %}, otherwise
        {% $\mathtt{smu} = \min(\mathtt{n}-1, \mathtt{mu} + \mathtt{ml})$ %}.
        The extra space is used to store U. *)
    type dimensions = {
        n   : int;  (** Matrix size: [n] by [n]. *)
        mu  : int;  (** Upper bandwidth. *)
        smu : int;  (** Storage upper bandwidth. *)
        ml  : int;  (** Lower bandwidth. *)
      }

    (** {4 Basic access} *)

    (** Returns a band matrix with the given {!dimensions} and all elements
        initialized to the given value.

        @nocvode <node> SUNBandMatrix *)
    val make : dimensions -> float -> t

    (** Returns an uninitialized band matrix with the given {!dimensions}.

        @nocvode <node> SUNBandMatrix *)
    val create : dimensions -> t

    (** [m, n = size a] returns the numbers of columns [m] and rows [n]
        of [a].

        @nocvode <node> SM_ROWS_B
        @nocvode <node> SM_COLUMNS_B
        @nocvode <node> SUNBandMatrix_Rows
        @nocvode <node> SUNBandMatrix_Columns *)
    val size : t -> int * int

    (** Returns the dimensions of a band matrix.

        @nocvode <node> SM_COLUMNS_B
        @nocvode <node> SM_UBAND_B
        @nocvode <node> SM_SUBAND_B
        @nocvode <node> SM_LBAND_B
        @nocvode <node> SUNBandMatrix_Columns
        @nocvode <node> SUNBandMatrix_UpperBandwidth
        @nocvode <node> SUNBandMatrix_StoredUpperBandwidth
        @nocvode <node> SUNBandMatrix_LowerBandwidth *)
    val dims : t -> dimensions

    (* TOPLEVEL-PRINTER: Matrix.Band.pp *)
    (** Pretty-print a band matrix using the
        {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
    val pp : Format.formatter -> t -> unit

    (** Pretty-print a band matrix using the
        {{:OCAML_DOC_ROOT(Format.html)} Format} module.
        The defaults are: [start="\["], [stop="\]"], [sep=";"],
        [indent=4], [itemsep=" "], and
        [item=fun f r c->Format.fprintf f "(%2d,%2d)=% -15e" r c] (see
        {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}).
        The [indent] argument specifies the indent for wrapped rows. *)
    val ppi : ?start:string -> ?stop:string -> ?sep:string
              -> ?indent:int -> ?itemsep:string
              -> ?item:(Format.formatter -> int -> int -> float -> unit)
              -> Format.formatter -> t -> unit

    (** [get a i j] returns the value at row [i] and column [j] of [a].
        Only rows and columns satisfying
        {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
        {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

        @nocvode <node> SM_ELEMENT_B *)
    val get : t -> int -> int -> float

    (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v].
        Only rows and columns satisfying
        {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
        {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

        @nocvode <node> SM_ELEMENT_B *)
    val set : t -> int -> int -> float -> unit

    (** [update a i j f] sets the value at row [i] and column [j] of [a]
        to [f v]. Only rows and columns satisfying
        {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
        {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

        @nocvode <node> SM_ELEMENT_B *)
    val update : t -> int -> int -> (float -> float) -> unit

    (** Direct access to the underlying storage array, which is accessed
        column first (unlike in {!get}).

        NB: The {!scale_add} operation, invoked either directly or from within
        a solver, will replace the underlying storage of its first matrix
        argument if the second matrix has a strictly larger bandwidth.
        Similarly, the {!blit} operation, invoked either directly or from
        within a solver, will replace the underlying storage of its second
        matrix argument if the first matrix has a strictly larger bandwidth.
        In both cases, any previously 'unwrapped' array is no longer
        associated with the matrix storage.

        NB: For {!Sundials.sundials_version} < 3.0.0, this access is
        potentially unsafe and {b must} only be used when the underlying
        storage is valid, which will be the case in callbacks.

        @nocvode <node> SM_CONTENT_B *)
    val unwrap : t -> Sundials.real_array2

    (** {4 Operations} *)

    (** Operations on dense matrices. *)

    val ops : (t, Nvector_serial.data, [>Nvector_serial.kind]) matrix_ops

    (** [scale_add c A B] calculates $A = cA + B$.

        NB: This operation, invoked either directly or from within a solver,
        will replace the underlying storage of its first matrix argument if
        the second matrix has a strictly larger bandwidth.

        @nocvode <node> SUNMatScaleAdd
        @nocvode <node> SUNMatScaleAdd_Band *)
    val scale_add : float -> t -> t -> unit

    (** [scale_addi c A] calculates $A = cA + I$.

        @nocvode <node> SUNMatScaleAddI
        @nocvode <node> SUNMatScaleAddI_Band *)
    val scale_addi : float -> t -> unit

    (** Compute the matrix-vector product $y = Ax$.
    
        @nocvode <node> SUNMatMatvec
        @nocvode <node> SUNMatMatvec_Band *)
    val matvec :
      t -> x:'nk Nvector_serial.any -> y:'nk Nvector_serial.any -> unit

    (** Fills a matrix with zeros.

        @nocvode <node> SUNMatZero
        @nocvode <node> SUNMatZero_Band *)
    val set_to_zero : t -> unit

    (** [blit src dst] copies the contents of [src] into [dst]. Both
        must have the same size.

        NB: This operation, invoked either directly or from within a solver,
        will replace the underlying storage of its second matrix argument if
        the first matrix has a strictly larger bandwidth.

        @nocvode <node> SUNMatCopy
        @nocvode <node> SUNMatCopy_Band *)
    val blit : t -> t -> unit

    (** [lrw, liw = space a] returns the storage requirements of [a] as
        [lrw] realtype words and [liw] integer words.
    
        @nocvode <node> SUNMatSpace
        @nocvode <node> SUNMatSpace_Band *)
    val space : t -> int * int

    (** {4 Low-level details} *)

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. Has no effect when
        {!Sundials.sundials_version} >= 3.0.0. *)
    val invalidate : t -> unit
  end

(** Sparse matrices *)
module Sparse :
  sig
    type csc  (* Compressed-sparse-column format ([CSC_MAT]). *)
    type csr  (* Compressed-sparse-row format ([CSR_MAT]). *)

    (* TODO: check documentation links *)
    (** A spare matrix. Values of this type are typically passed to linear
        solver callback functions (like {!Cvode.dense_jac_fn}, 
        {!Ida.dense_jac_fn}, and {!Kinsol.dense_jac_fn}).

        @nocvode <node> The SUNMatrix_Sparse implementation *)
    type 's t
    type t_csc = csc t
    type t_csr = csr t

    (** Array of row or column indices *)
    type index_array =
      (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t

    (** {4 Basic access} *)

    (** [make_csc m n nnz] returns an [m] by [n] sparse matrix in
        compressed-sparse-column format with a potential
        for [nnz] non-zero elements. All elements are initially zero.

        @nocvode <node> SUNSparseMatrix *)
    val make_csc : int -> int -> int -> csc t

    (** [make_csr m n nnz] returns an [m] by [n] sparse matrix in
        compressed-sparse-row format with a potential
        for [nnz] non-zero elements. All elements are initially zero.

        @since 2.7.0
        @raise Sundials.NotImplementedBySundialsVersion CSR format not available.
        @nocvode <node> SUNSparseMatrix *)
    val make_csr : int -> int -> int -> csr t

    (** Creates a sparse matrix in compressed-sparse-column format from a dense
        matrix by copying all values of magnitude greater than the given
        tolerance.

        @nocvode <node> SUNSparseFromDenseMatrix *)
    val csc_from_dense : float -> Dense.t -> csc t

    (** Creates a sparse matrix in compressed-sparse-row format from a dense
        matrix by copying all values of magnitude greater than the given
        tolerance.

        @since 2.7.0
        @raise Sundials.NotImplementedBySundialsVersion CSR format not available.
        @nocvode <node> SUNSparseFromDenseMatrix *)
    val csr_from_dense : float -> Dense.t -> csr t

    (** Creates a sparse matrix in compressed-sparse-column format from a band
        matrix by copying all values of magnitude greater than the given
        tolerance.

        @nocvode <node> SUNSparseFromBandMatrix *)
    val csc_from_band : float -> Band.t -> csc t

    (** Creates a sparse matrix in compressed-sparse-row format from a band
        matrix by copying all values of magnitude greater than the given
        tolerance.

        @since 2.7.0
        @raise Sundials.NotImplementedBySundialsVersion CSR format not available.
        @nocvode <node> SUNSparseFromBandMatrix *)
    val csr_from_band : float -> Band.t -> csr t

    (** [m, n = size a] returns the numbers of columns [m] and rows [n]
        of [a].

        @nocvode <node> SM_ROWS_S
        @nocvode <node> SM_COLUMNS_S
        @nocvode <node> SUNSparseMatrix_Rows
        @nocvode <node> SUNSparseMatrix_Columns *)
    val size : 's t -> int * int

    (** [nnz, np = dims m] returns the allocated number of nonzeros [nnz] and
        of the number [np] of columns (for csc) or rows (for csr) in the
        matrix [m].

        @nocvode <node> SM_NNZ_S
        @nocvode <node> SM_NP_S
        @nocvode <node> SUNSparseMatrix_NNZ
        @nocvode <node> SUNSparseMatrix_NP *)
    val dims : 's t -> int * int

    (* TOPLEVEL-PRINTER: Matrix.Sparse.pp *)
    (** Pretty-print a sparse matrix using the
        {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
    val pp : Format.formatter -> 's t -> unit

    (** Pretty-print a sparse matrix using the
        {{:OCAML_DOC_ROOT(Format.html)} Format} module.
        The defaults are: [start="\["], [stop="\]"], [sep=";"],
        [indent=4], [itemsep=" "],
        [rowcol=fun f i->Format.fprintf f "%2d: " i], and
        [item=fun f r c->Format.fprintf f "(%2d,%2d)=% -15e" r c] (see
        {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}).
        The [indent] argument specifies the indent for wrapped rows. *)
    val ppi : ?start:string -> ?stop:string -> ?sep:string
              -> ?indent:int -> ?itemsep:string
              -> ?rowcol:(Format.formatter -> int -> unit)
              -> ?item:(Format.formatter -> int -> float -> unit)
              -> Format.formatter -> 's t -> unit

    (** [set_col a j idx] sets the data index of column [j] to [idx]. *)
    val set_col : csc t -> int -> int -> unit

    (** [get_col a j] returns the data index of column [j]. *)
    val get_col : csc t -> int -> int

    (** [set_row a j idx] sets the data index of row [j] to [idx]. *)
    val set_row : csr t -> int -> int -> unit

    (** [get_row a j] returns the data index of row [j]. *)
    val get_row : csr t -> int -> int

    (** [set a idx i v] sets the [idx]th row/column to [i] and its value
        to [v]. *)
    val set : 'f t -> int -> int -> float -> unit

    (** [r, v = get a idx] returns the row/column [r] and value [v] at the
        [idx]th position. *)
    val get : 'f t -> int -> int * float

    (** Direct access to the underlying sparse storage arrays.
        In the call [vals, ptrs, data = unwrap m] ,
        - [vals] contains the row (if [csc]) or column (if [csr]) indices of
                 nonzero matrix entries in [data],
        - [ptrs] contains the indices of the columnn (if [csc]) or
                 row (if [csr]) entries in [data] and [vals], and
        - [data] contains the values of the nonzero entries.

        NB: The {!scale_add}, {!scale_addi}, {!blit}, and {!resize} functions,
        invoked either directly or from within a solver, may replace the
        underlying storage. In these cases, any previously 'unwrapped' arrays
        are no longer associated with the matrix storage.

        NB: For {!Sundials.sundials_version} < 3.0.0, this access is
        potentially unsafe and {b must} only be used when the underlying
        storage is valid, which will be the case in callbacks unless the
        {!scale_add}, {!scale_addi}, {!blit}, and {!resize} functions are
        used.

        @nocvode <node> SM_INDEXVALS_S
        @nocvode <node> SM_INDEXPTRS_S
        @nocvode <node> SM_DATA_S
    *)
    val unwrap : 's t -> index_array * index_array * Sundials.RealArray.t

    (** Reallocates the underlying arrays to the given number of non-zero
        elements, or otherwise to the current number of non-zero elements .
        
        NB: The {!resize} operation may replace the underlying storage of the
        matrix argument. In this case, any previously 'unwrapped' array is no
        longer associated with the matrix storage.

        @nocvode <node> SUNSparseMatrix_Reallocate *)
    val resize : ?nnz:int -> 's t -> unit

    (** {4 Operations} *)

    (** Operations on dense matrices. *)
    val ops : ('s t, Nvector_serial.data, [>Nvector_serial.kind]) matrix_ops

    (** [scale_add c A B] calculates $A = cA + B$.

        NB: The {!scale_add} operation, invoked either directly or from within
        a solver, may replace the underlying storage of its first matrix
        argument does not contain the sparsity of the second matrix argument.
        In this case, any previously 'unwrapped' array is no longer associated
        with the matrix storage.

        @nocvode <node> SUNMatScaleAdd
        @nocvode <node> SUNMatScaleAdd_Sparse *)
    val scale_add    : float -> 's t -> 's t -> unit

    (** [scale_addi c A] calculates $A = cA + I$.

        NB: The {!scale_add} operation, invoked either directly or from within
        a solver, may replace the underlying storage of its matrix argument
        if it does not already contain a complete diagonal. In this
        case, any previously 'unwrapped' array is no longer associated with
        the matrix storage.

        @nocvode <node> SUNMatScaleAddI
        @nocvode <node> SUNMatScaleAddI_Sparse *)
    val scale_addi   : float -> 's t -> unit

    (** Compute the matrix-vector product $y = Ax$.
    
        @nocvode <node> SUNMatMatvec
        @nocvode <node> SUNMatMatvec_Sparse *)
    val matvec :
      's t -> x:'nk Nvector_serial.any -> y:'nk Nvector_serial.any -> unit

    (** Fills a matrix with zeros.

        @nocvode <node> SUNMatZero
        @nocvode <node> SUNMatZero_Sparse *)
    val set_to_zero : 's t -> unit

    (** [blit src dst] copies the contents of [src] into [dst]. Both
        must have the same size.

        NB: This operation, invoked either directly or from within a solver,
        may replace the underlying storage of its second matrix argument if it
        does not contain the sparsity of the first matrix argument.
        In this case, any previously 'unwrapped' array is no longer associated
        with the matrix storage.

        @nocvode <node> SUNMatCopy
        @nocvode <node> SUNMatCopy_Sparse *)
    val blit : 's t -> 's t -> unit

    (** [lrw, liw = space a] returns the storage requirements of [a] as
        [lrw] realtype words and [liw] integer words.
    
        @nocvode <node> SUNMatSpace
        @nocvode <node> SUNMatSpace_Sparse *)
    val space : 's t -> int * int

    (** {4 Low-level details} *)

    (** [set_rowval a idx i] sets the [idx]th row to [i]. *)
    val set_rowval : csc t -> int -> int -> unit

    (** [r = get_rowval a idx] returns the row [r] at the [idx]th position. *)
    val get_rowval : csc t -> int -> int

    (** [set_colval a idx i] sets the [idx]th column to [i]. *)
    val set_colval : csr t -> int -> int -> unit

    (** [c = get_colval a idx] returns the column [c] at the [idx]th
        position. *)
    val get_colval : csr t -> int -> int

    (** [set_data a idx v] sets the value of the [idx]th row [v]. *)
    val set_data : 'f t -> int -> float -> unit

    (** [v = get_data a idx] returns the value [v] at the [idx]th position. *)
    val get_data : 'f t -> int -> float

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. Has no effect when
        {!Sundials.sundials_version} >= 3.0.0. *)
    val invalidate : 's t -> unit
  end

type csc = Sparse.csc
type csr = Sparse.csr

(** Distinguishes a library-supplied matrix from a custom one. *)
type standard

(** Distinguishes a user-supplied matrix from a standard one. *)
type custom

(** A generic matrix with a payload of type ['m]. The ['k] type argument tracks
    whether the matrix is {!standard} or {!custom}. The ['nd] and ['nk] type
    arguments track the compatiblity of the {{!matrix_ops}m_matvec} vector
    parameters.
 
    @nocvode <node> SUNMatrix *)
type ('k, 'm, 'nd, 'nk) t

type 'nk dense =
  (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type 'nk band =
  (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type ('s, 'nk) sparse =
  (standard, 's Sparse.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

(** [make m n x] returns an [m] by [n] (dense) matrix with elements set
    to [x].

    @nocvode <node> SUNDenseMatrix *)
val make_dense : int -> int -> float -> 'nk dense

(** Creates a (dense) matrix by wrapping an existing dense matrix. The two
    values share the same underlying storage.

    @nocvode <node> SUNDenseMatrix *)
val wrap_dense : Dense.t -> 'nk dense

(** [make dimensions x] returns a (band) matrix with the given
    {!dimensions} and all elements initialized to [x].

    @nocvode <node> SUNBandMatrix *)
val make_band : Band.dimensions -> float -> 'nk band

(** Creates a (band) matrix by wrapping an existing band matrix. The two
    values share the same underlying storage.

    @nocvode <node> SUNBandMatrix *)
val wrap_band : Band.t -> 'nk band

(** [make m n nnz] returns an [m] by [n] (sparse compressed-sparse-column)
    matrix with a potential for [nnz] non-zero elements. All elements are
    initially zero.

    @nocvode <node> SUNSparseMatrix *)
val make_sparse_csc : int -> int -> int -> (csc, 'nk) sparse

(** [make m n nnz] returns an [m] by [n] (sparse compressed-sparse-row)
    matrix with a potential for [nnz] non-zero elements. All elements are
    initially zero.

    @nocvode <node> SUNSparseMatrix *)
val make_sparse_csr : int -> int -> int -> (csr, 'nk) sparse

(** Creates a (sparse) matrix by wrapping an existing sparse matrix. The two
    values share the same underlying storage.

    @nocvode <node> SUNSparseMatrix *)
val wrap_sparse : 's Sparse.t -> ('s, 'nk) sparse

(** Wrap a custom matrix value.

    @nocvode <node> Description of the SUNMatrix module *)
val wrap_custom : ('m, 'nd, 'nk) matrix_ops -> 'm -> (custom, 'm, 'nd, 'nk) t

(** Matrix internal type identifiers.
 
    @nocvode <node> SUNMatrix_ID *)
type id =
  | Dense
  | Band
  | Sparse
  | Custom

(** Return a record of matrix operations. *)
val get_ops : ('k, 'm, 'nd, 'nk) t -> ('m, 'nd, 'nk) matrix_ops

(** Return the internal type identifier of a matrix.

    @nocvode <node> SUNMatGetID *)
val get_id : ('k, 'm, 'nd, 'nk) t -> id

(** Direct access to the underlying storage array, which is accessed
    column first (unlike in {!get}).

    @nocvode <node> SM_CONTENT_B *)
val unwrap : ('k, 'm, 'nd, 'nk) t -> 'm

(** {4 Operations} *)

(** [scale_add c A B] calculates $A = cA + B$.

    @nocvode <node> SUNMatScaleAdd *)
val scale_add : float -> ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit

(** [scale_addi c A] calculates $A = cA + I$.

    @nocvode <node> SUNMatScaleAddI *)
val scale_addi : float -> ('k, 'm, 'nd, 'nk) t -> unit

(** Compute the matrix-vector product $y = Ax$.

    @nocvode <node> SUNMatMatvec *)
val matvec :
  ('k, 'm, 'nd, 'nk) t -> x:('nd, 'nk) Nvector.t
                       -> y:('nd, 'nk) Nvector.t -> unit

(** Fills a matrix with zeros.

    @nocvode <node> SUNMatZero *)
val set_to_zero : ('k, 'm, 'nd, 'nk) t -> unit

(** [blit src dst] copies the contents of [src] into [dst]. Both
    must have the same size.

    @nocvode <node> SUNMatCopy *)
val blit : ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit

(** [lrw, liw = space a] returns the storage requirements of [a] as
    [lrw] realtype words and [liw] integer words.

    @nocvode <node> SUNMatSpace *)
val space : ('k, 'm, 'nd, 'nk) t -> int * int

