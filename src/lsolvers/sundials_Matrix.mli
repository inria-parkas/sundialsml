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

    @matrix <index.html> Matrix Data Structures
    @since 3.0.0 *)

open Sundials

(** {2:shared Shared definitions} *)

(** Generic operations that all matrix types must implement. Failure is
    signalled by raising an exception.

    There are two type variables:
    - ['m] represents the matrix data manipulated by the operations;
    - ['d] is the type of the vector arguments to [m_matvec].

    @matrix <SUNMatrix_API_link.html#description-of-the-sunmatrix-modules> Description of the SUNMatrix modules *)
type ('m, 'd) matrix_ops = { (* {{{ *)
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

  m_matvec_setup : ('m -> unit) option;
  (** Performs any setup necessary to perform a matrix-vector product. *)

  m_matvec    : 'm -> 'd -> 'd -> unit;
  (** [m_matvec a x y] calculates $y = ax$. *)

  m_space     : 'm -> int * int;
  (** [lrw, liw = m_space a] returns the number of realtype words [lrw] and
      integer words [liw] required to store the matrix [a]. *)
} (* }}} *)

(** Raised in {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0 on an attempt to
    access a value that has become invalid. Such values refer to matrices
    that no longer exist in the underlying library. Values never
    become invalid in {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0. *)
exception Invalidated

(** Raised if matrix operation arguments are mutually incompatible. *)
exception IncompatibleArguments

(** Raised if a zero diagonal element is found during factorization using a
    low-level routine like
    {{!Sundials_LinearSolver.Iterative.Algorithms.qr_fact}LinearSolver.Iterative.Algorithms.qr_fact},
    {{!Sundials_LinearSolver.Iterative.Algorithms.qr_sol}LinearSolver.Iterative.Algorithms.qr_sol},
    or {!ArrayDense.getrf}.
    The argument gives the equation number (from 1). *)
exception ZeroDiagonalElement of int

(** {2:content Matrix content} *)

(** Dense matrices

    @matrix_data <SUNMatrix_links.html#the-sunmatrix-dense-module> The SUNMATRIX_DENSE module *)
module Dense : sig (* {{{ *)

  (** A dense matrix. Values of this type are typically passed to linear
      solver callback functions (like {!Cvode.Dls.jac_fn},
      {!Ida.Dls.jac_fn}, and {!Kinsol.Dls.jac_fn}). *)
  type t

  (** {3:dense_basic Basic access} *)

  (** [make m n x] returns an [m] by [n] dense matrix with elements set
      to [x].

      @matrix_data SUNDenseMatrix *)
  val make : int -> int -> float -> t

  (** [create m n] returns an uninitialized [m] by [n] dense matrix.

      @matrix_data SUNDenseMatrix *)
  val create : int -> int -> t

  (** [m, n = size a] returns the numbers of rows [m] and columns [n]
      of [a].

      @matrix_data SM_ROWS_D
      @matrix_data SM_COLUMNS_D
      @matrix_data SUNDenseMatrix_Rows
      @matrix_data SUNDenseMatrix_Columns *)
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
            -> unit
            -> Format.formatter -> t -> unit

  (** [get a i j] returns the value at row [i] and column [j] of [a].

      @matrix_data SM_ELEMENT_D *)
  val get : t -> int -> int -> float

  (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v].

      @matrix_data SM_ELEMENT_D *)
  val set : t -> int -> int -> float -> unit

  (** [update a i j f] sets the value at row [i] and column [j] of [a]
      to [f v].

      @matrix_data SM_ELEMENT_D *)
  val update : t -> int -> int -> (float -> float) -> unit

  (** Direct access to the underlying storage array, which is accessed
      column first (unlike in {!get}).

      NB: For {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0, this access is
      potentially unsafe and {b must} only be used when the underlying
      storage is valid, which will be the case in callbacks.

      @matrix_data SM_CONTENT_D *)
  val unwrap : t -> RealArray2.data

  (** {3:dense_ops Operations} *)

  (** Operations on dense matrices. *)
  val ops : (t, Nvector_serial.data) matrix_ops

  (** [scale_add c A B] calculates $A = cA + B$.

      @matrix SUNMatScaleAdd (SUNMatScaleAdd_Dense) *)
  val scale_add : float -> t -> t -> unit

  (** [scale_addi c A] calculates $A = cA + I$.

      @matrix SUNMatScaleAddI (SUNMatScaleAddI_Dense) *)
  val scale_addi : float -> t -> unit

  (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

      @matrix SUNMatMatvec (SUNMatMatvec_Dense) *)
  val matvec : t -> RealArray.t -> RealArray.t -> unit

  (** Fills a matrix with zeros.

      @matrix SUNMatZero (SUNMatZero_Dense) *)
  val set_to_zero : t -> unit

  (** [blit ~src ~dst] copies the contents of [src] into [dst]. Both
      must have the same size.

      @matrix SUNMatCopy (SUNMatCopy_Dense) *)
  val blit : src:t -> dst:t -> unit

  (** [lrw, liw = space a] returns the storage requirements of [a] as
      [lrw] realtype words and [liw] integer words.

      @matrix SUNMatSpace (SUNMatSpace_Dense) *)
  val space : t -> int * int

  (** {3:dense_lowlevel Low-level details} *)

  (** Called internally when the corresponding value in the underlying
      library ceases to exist. Has no effect when
      {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0. *)
  val invalidate : t -> unit

end (* }}} *)

(** Banded matrices

    @matrix <SUNMatrix_links.html#the-sunmatrix-band-module> The SUNMATRIX_BAND Module *)
module Band : sig (* {{{ *)

  (** A band matrix. Values of this type are typically passed to linear
      solver callback functions (like {!Cvode.Dls.jac_fn},
      {!Ida.Dls.jac_fn}, and {!Kinsol.Dls.jac_fn}). *)
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

  (** {3:band_basic Basic access} *)

  (** Returns a band matrix with the given {!dimensions} and all elements
      initialized to the given value.

      @matrix_data SUNBandMatrixStorage *)
  val make : dimensions -> float -> t

  (** Returns an uninitialized band matrix with the given {!dimensions}.

      @matrix_data SUNBandMatrixStorage *)
  val create : dimensions -> t

  (** [m, n = size a] returns the numbers of rows [m] and columns [n] of [a].
      
      NB: [m] and [n] are always equal for band matrices.

      @matrix_data SM_ROWS_B
      @matrix_data SM_COLUMNS_B
      @matrix_data SUNBandMatrix_Rows
      @matrix_data SUNBandMatrix_Columns *)
  val size : t -> int * int

  (** Returns the dimensions of a band matrix.

      @matrix_data SM_COLUMNS_B
      @matrix_data SM_UBAND_B
      @matrix_data SM_SUBAND_B
      @matrix_data SM_LBAND_B
      @matrix_data SUNBandMatrix_Columns
      @matrix_data SUNBandMatrix_UpperBandwidth
      @matrix_data SUNBandMatrix_StoredUpperBandwidth
      @matrix_data SUNBandMatrix_LowerBandwidth *)
  val dims : t -> dimensions

  (* TOPLEVEL-PRINTER: Matrix.Band.pp *)
  (** Pretty-print a band matrix using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
  val pp : Format.formatter -> t -> unit

  (** Pretty-print a band matrix using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module.
      The defaults are: [start="\["], [stop="\]"], [sep=";"],
      [indent=4], [itemsep=" "], [empty="           ~           "] and
      [item=fun f r c->Format.fprintf f "(%2d,%2d)=% -15e" r c] (see
      {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}).
      The [indent] argument specifies the indent for wrapped rows. *)
  val ppi : ?start:string -> ?stop:string -> ?sep:string
            -> ?indent:int -> ?itemsep:string -> ?empty:string
            -> ?item:(Format.formatter -> int -> int -> float -> unit)
            -> unit
            -> Format.formatter -> t -> unit

  (** [get a i j] returns the value at row [i] and column [j] of [a].
      Only rows and columns satisfying
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

      @matrix_data SM_ELEMENT_B *)
  val get : t -> int -> int -> float

  (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v].
      Only rows and columns satisfying
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

      @matrix_data SM_ELEMENT_B *)
  val set : t -> int -> int -> float -> unit

  (** [update a i j f] sets the value at row [i] and column [j] of [a]
      to [f v]. Only rows and columns satisfying
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

      @matrix_data SM_ELEMENT_B *)
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

      NB: For {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0, this access is
      potentially unsafe and {b must} only be used when the underlying
      storage is valid, which will be the case in callbacks.

      @matrix_data SM_CONTENT_B *)
  val unwrap : t -> RealArray2.data

  (** {3:band_ops Operations} *)

  (** Operations on band matrices. *)
  val ops : (t, Nvector_serial.data) matrix_ops

  (** [scale_add c A B] calculates $A = cA + B$.

      NB: This operation, invoked either directly or from within a solver,
      will replace the underlying storage of its first matrix argument if
      the second matrix has a strictly larger bandwidth.

      @matrix SUNMatScaleAdd (SUNMatScaleAdd_Band) *)
  val scale_add : float -> t -> t -> unit

  (** [scale_addi c A] calculates $A = cA + I$.

      @matrix SUNMatScaleAddI (SUNMatScaleAddI_Band) *)
  val scale_addi : float -> t -> unit

  (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

      @matrix SUNMatMatvec (SUNMatMatvec_Band) *)
  val matvec : t -> RealArray.t -> RealArray.t -> unit

  (** Fills a matrix with zeros.

      @matrix SUNMatZero (SUNMatZero_Band) *)
  val set_to_zero : t -> unit

  (** [blit ~src ~dst] copies the contents of [src] into [dst]. Both
      must have the same size.

      NB: This operation, invoked either directly or from within a solver,
      will replace the underlying storage of its second matrix argument if
      the first matrix has a strictly larger bandwidth.

      @matrix SUNMatCopy (SUNMatCopy_Band) *)
  val blit : src:t -> dst:t -> unit

  (** [lrw, liw = space a] returns the storage requirements of [a] as
      [lrw] realtype words and [liw] integer words.

      @matrix SUNMatSpace (SUNMatSpace_Band) *)
  val space : t -> int * int

  (** {3:band_lowlevel Low-level details} *)

  (** Called internally when the corresponding value in the underlying
      library ceases to exist. Has no effect when
      {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0. *)
  val invalidate : t -> unit

end (* }}} *)

(** Sparse matrices

    @matrix <SUNMatrix_links.html#the-sunmatrix-sparse-module> The SUNMATRIX_SPARSE Module *)
module Sparse : sig (* {{{ *)

  type csc (** Compressed-sparse-column format. *)

  type csr (** Compressed-sparse-row format. *)

  (** Matrix storage formats. *)
  type _ sformat =
    | CSC : csc sformat (** Compressed-sparse-column format ([CSC_MAT]). *)
    | CSR : csr sformat (** Compressed-sparse-row format ([CSR_MAT]). *)

  (** A sparse matrix. Values of this type are typically passed to linear
      solver callback functions (like {!Cvode.Dls.jac_fn},
      {!Ida.Dls.jac_fn}, and {!Kinsol.Dls.jac_fn}).

      The type argument specifies the storage format, either
      {{!csc}compressed-sparse-column format} or
      {{!csr}compressed-sparse-row format}. *)
  type 's t

  (** Array of row or column indices *)
  type index_array =
    (Index.t, Index.index_elt, Bigarray.c_layout)
    Bigarray.Array1.t

  (** {3:sparse_basic Basic access} *)

  (** [make fmt m n nnz] returns an [m] by [n] sparse matrix in the
      specified format with a potential for [nnz] non-zero elements.
      All elements are initially zero.

      The {{!sformat}CSR} format is only available from Sundials 2.7.0
      onwards.

      @matrix_data SUNSparseMatrix *)
  val make : 's sformat -> int -> int -> int -> 's t

  (** Creates a sparse matrix in in the specified format from a dense matrix
      by copying all values of magnitude greater than the given tolerance.

      The {{!sformat}CSR} format is only available from Sundials 2.7.0
      onwards.

      @matrix_data SUNSparseFromDenseMatrix *)
  val from_dense : 's sformat -> float -> Dense.t -> 's t

  (** Creates a sparse matrix in the specified format from a band matrix by
      copying all values of magnitude greater than the given tolerance.

      The {{!sformat}CSR} format is only available from Sundials 2.7.0
      onwards.

      @matrix_data SUNSparseFromBandMatrix *)
  val from_band : 's sformat -> float -> Band.t -> 's t

  (** Return the matrix format. *)
  val sformat : 's t -> 's sformat

  (** Returns true iff the matrix format is {{!Sparse.sformat}CSC}.
      It is essentially a version of {!sformat} with less typing
      complications. *)
  val is_csc : 's t -> bool

  (** [m, n = size a] returns the numbers of rows [m] and columns [n]
      of [a].

      @matrix_data SM_ROWS_S
      @matrix_data SM_COLUMNS_S
      @matrix_data SUNSparseMatrix_Rows
      @matrix_data SUNSparseMatrix_Columns *)
  val size : 's t -> int * int

  (** [nnz, np = dims m] returns the allocated number of nonzeros [nnz] and
      of the number [np] of columns (for csc) or rows (for csr) in the
      matrix [m].

      @matrix_data SM_NNZ_S
      @matrix_data SM_NP_S
      @matrix_data SUNSparseMatrix_NNZ
      @matrix_data SUNSparseMatrix_NP *)
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
            -> unit
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
      - [ptrs] contains the indices of the column (if [csc]) or
               row (if [csr]) entries in [data] and [vals], and
      - [data] contains the values of the nonzero entries.

      NB: The {!scale_add}, {!scale_addi}, {!blit}, and {!resize} functions,
      invoked either directly or from within a solver, may replace the
      underlying storage. In these cases, any previously 'unwrapped' arrays
      are no longer associated with the matrix storage.

      NB: For {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0, this access is
      potentially unsafe and {b must} only be used when the underlying
      storage is valid, which will be the case in callbacks unless the
      {!scale_add}, {!scale_addi}, {!blit}, and {!resize} functions are
      used.

      @matrix_data SM_INDEXVALS_S
      @matrix_data SM_INDEXPTRS_S
      @matrix_data SM_DATA_S
  *)
  val unwrap : 's t -> index_array * index_array * RealArray.t

  (** Reallocates the underlying arrays to the given number of non-zero
      elements, or otherwise to the current number of non-zero elements .

      NB: The {!resize} operation may replace the underlying storage of the
      matrix argument. In this case, any previously 'unwrapped' array is no
      longer associated with the matrix storage.

      @matrix_data SUNSparseMatrix_Realloc
      @matrix_data SUNSparseMatrix_Reallocate *)
  val resize : ?nnz:int -> 's t -> unit

  (** {3:sparse_ops Operations} *)

  (** Operations on sparse matrices. *)
  val ops : ('s t, Nvector_serial.data) matrix_ops

  (** [scale_add c A B] calculates $A = cA + B$.

      NB: The {!scale_add} operation, invoked either directly or from within
      a solver, may replace the underlying storage of its first matrix
      argument does not contain the sparsity of the second matrix argument.
      In this case, any previously 'unwrapped' array is no longer associated
      with the matrix storage.

      @matrix SUNMatScaleAdd (SUNMatScaleAdd_Sparse) *)
  val scale_add    : float -> 's t -> 's t -> unit

  (** [scale_addi c A] calculates $A = cA + I$.

      NB: The {!scale_add} operation, invoked either directly or from within
      a solver, may replace the underlying storage of its matrix argument
      if it does not already contain a complete diagonal. In this
      case, any previously 'unwrapped' array is no longer associated with
      the matrix storage.

      @matrix SUNMatScaleAddI (SUNMatScaleAddI_Sparse) *)
  val scale_addi   : float -> 's t -> unit

  (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

      @matrix SUNMatMatvec (SUNMatMatvec_Sparse) *)
  val matvec : 's t -> RealArray.t -> RealArray.t -> unit

  (** Fills a matrix with zeros.

      @matrix SUNMatZero (SUNMatZero_Sparse) *)
  val set_to_zero : 's t -> unit

  (** [blit ~src ~dst] copies the contents of [src] into [dst]. Both
      must have the same size.

      NB: This operation, invoked either directly or from within a solver,
      may replace the underlying storage of its second matrix argument if it
      does not contain the sparsity of the first matrix argument.
      In this case, any previously 'unwrapped' array is no longer associated
      with the matrix storage.

      @matrix SUNMatCopy (SUNMatCopy_Sparse) *)
  val blit : src:'s t -> dst:'s t -> unit

  (** Create a new sparse matrix in {{!sformat}CSR} format from the contents
      of an existing one in {{!sformat}CSC} format.

      @since 5.2.0
      @matrix_data SUNSparseMatrix_ToCSR *)
  val copy_to_csr : csc t -> csr t

  (** Create a new sparse matrix in {{!sformat}CSC} format from the contents
      of an existing one in {{!sformat}CSR} format.

      @since 5.2.0
      @matrix_data SUNSparseMatrix_ToCSC *)
  val copy_to_csc : csr t -> csc t

  (** [lrw, liw = space a] returns the storage requirements of [a] as
      [lrw] realtype words and [liw] integer words.

      @matrix SUNMatSpace (SUNMatSpace_Sparse) *)
  val space : 's t -> int * int

  (** {3:sparse_lowlevel Low-level details} *)

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
      {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0. *)
  val invalidate : 's t -> unit

end (* }}} *)

(** {2:array Arrays as matrices} *)

(** General purpose dense matrix operations on arrays.
    @nodoc The DENSE Module *)
module ArrayDense : sig (* {{{ *)

  (** A dense matrix accessible directly through a
      {{:OCAML_DOC_ROOT(Bigarray.Array2.html)} Bigarray}.

      @nodoc SUNDlsMat_newDenseMat *)
  type t = RealArray2.t

  (** {3:arraydense_basic Basic access} *)

  (** [make m n x] returns an [m] by [n] array dense matrix with elements set
      to [x].

      @nodoc SUNDlsMat_newDenseMat *)
  val make : int -> int -> float -> t

  (** [create m n] returns an uninitialized [m] by [n] array dense matrix.

       @nodoc SUNDlsMat_newDenseMat *)
  val create : int -> int -> t

  (** [m, n = size a] returns the numbers of rows [m] and columns [n]
      of [a]. *)
  val size : t -> int * int

  (** Pretty-print an array dense matrix using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
  val pp : Format.formatter -> t -> unit

  (** [get a i j] returns the value at row [i] and column [j] of [a]. *)
  val get : t -> int -> int -> float

  (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v]. *)
  val set : t -> int -> int -> float -> unit

  (** [update a i j f] sets the value at row [i] and column [j] of [a]
      to [f v]. *)
  val update : t -> int -> int -> (float -> float) -> unit

  (** Direct access to the underlying storage array, which is accessed
      column first (unlike in {!get}). *)
  val unwrap : t -> RealArray2.data

  (** {3:arraydense_ops Operations} *)

  (** Operations on array-based dense matrices. *)
  val ops : (t, RealArray.t) matrix_ops

  (** [scale_add c A B] calculates $A = cA + B$. *)
  val scale_add : float -> t -> t -> unit

  (** [scale_addi c A] calculates $A = cA + I$. *)
  val scale_addi : float -> t -> unit

  (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

      @nodoc denseMatvec
      @since 2.6.0 *)
  val matvec : t -> RealArray.t -> RealArray.t -> unit

  (** Fills the matrix with zeros.

      @nodoc SUNDlsMat_SetToZero *)
  val set_to_zero    : t -> unit

  (** [blit ~src ~dst] copies the contents of [src] into [dst]. Both
      must have the same size.

      @nodoc SUNDlsMat_denseCopy *)
  val blit  : src:t -> dst:t -> unit

  (** [lrw, liw = space a] returns the storage requirements of [a] as
      [lrw] realtype words and [liw] integer words. *)
  val space : t -> int * int

  (** {3:arraydense_calcs Calculations} *)

  (** Increments a square matrix by the identity matrix.

      @nodoc SUNDlsMat_denseAddIdentity *)
  val add_identity : t -> unit

  (** Multiplies each element by a constant.

      @nodoc SUNDlsMat_denseScale *)
  val scale : float -> t -> unit

  (** [getrf a p] performs the LU factorization of the square matrix [a] with
      partial pivoting according to [p]. The values in [a] are overwritten
      with those of the calculated L and U matrices. The diagonal belongs to
      U. The diagonal of L is all 1s. Multiplying L by U gives a permutation
      of [a], according to the values of [p]: [p.{k} = j] means that rows [k]
      and [j] were swapped (in order, where [p.{0}] swaps against the
      original matrix [a]).

      @nodoc SUNDlsMat_denseGETRF
      @raise ZeroDiagonalElement Zero found in matrix diagonal *)
  val getrf : t -> LintArray.t -> unit

  (** [getrs a p b] finds the solution of [ax = b] using an LU factorization
      found by {!getrf}. Both [p] and [b] must have the same number of rows
      as [a].

      @nodoc SUNDlsMat_denseGETRS *)
  val getrs : t -> LintArray.t -> RealArray.t -> unit

  (** Like {!getrs} but stores [b] starting at a given offset. *)
  val getrs'
        : t -> LintArray.t -> RealArray.t -> int -> unit

  (** Performs Cholesky factorization of a real symmetric positive matrix.

      @nodoc SUNDlsMat_densePOTRF *)
  val potrf : t -> unit

  (** [potrs a b] finds the solution of [ax = b] using the Cholesky
      factorization found by {!potrf}. [a] must be an n by n matrix and [b]
      must be of length n.

      @nodoc SUNDlsMat_densePOTRS *)
  val potrs : t -> RealArray.t -> unit

  (** [geqrf a beta work] performs the QR factorization of [a]. [a] must be
      an [m] by [n] matrix, where [m >= n]. The [beta] vector must have
      length [n]. The [work] vector must have length [m].

      @nodoc SUNDlsMat_denseGEQRF *)
  val geqrf : t -> RealArray.t -> RealArray.t -> unit

  (** [ormqr q beta v w work] computes the product {% w = qv %}. [Q] is
      an [m] by [n] matrix calculated using {!geqrf} with [m >= n],
      [beta] has length [n], [v] has length [n], [w] has length [m], and
      [work] has length [m].

      @param q       matrix modified by {!geqrf}
      @param beta    vector passed to {!geqrf}
      @param v       vector multiplier
      @param w       result vector
      @param work    temporary vector used in the calculation
      @nodoc SUNDlsMat_denseORMQR *)
  val ormqr :
    a:t -> beta:RealArray.t -> v:RealArray.t
      -> w:RealArray.t -> work:RealArray.t -> unit

end (* }}} *)

(** General-purpose band matrix operations on arrays.
    @nodoc The BAND Module *)
module ArrayBand : sig (* {{{ *)

  type smu = int (** Storage upper-bandwidth. *)

  type mu  = int (** Upper-bandwidth. *)

  type ml  = int (** Lower-bandwidth. *)

  (** A band matrix accessible directly through a
      {{:OCAML_DOC_ROOT(Bigarray.Array2.html)} Bigarray}.

      The layout of these arrays are characterized by the {e storage upper
      bandwidth} [smu]. Given an array [a], the first dimension indexes
      the diagonals, with the main diagonal ({% $i = j$ %}) at [a.{smu, _}].
      The value in the [i]th row and [j]th column provided
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} is at
      [a.{i - j + smu, j}].

      @nodoc SUNDlsMat_newBandMat *)
  type t = RealArray2.t * (smu * mu * ml)

  (** {3:arrayband_basic Basic access} *)

  (** [make (smu, mu, ml) n v] returns an [n] by [n] band matrix with
      storage upper bandwidth [smu], upper bandwidth [sm],
      lower half-bandwidth [ml], and all elements initialized to [v].

      If the result will not be LU factored then
      {% $\mathtt{smu} = \mathtt{mu}$ %}, otherwise
      {% $\mathtt{smu} = \min(\mathtt{n}-1, \mathtt{mu} + \mathtt{ml})$ %}.
      The extra space is used to store U after a call to {!gbtrf}.

      @nodoc SUNDlsMat_newBandMat *)
  val make : smu * mu * ml -> int -> float -> t

  (** [create smu ml n] returns an uninitialized [n] by [n] band matrix with
      storage upper bandwidth [smu] and lower half-bandwidth [ml].

      @nodoc SUNDlsMat_newBandMat *)
  val create : smu * mu * ml -> int -> t

  (** [m, n = size a] returns the numbers of rows [m] and columns [n] of [a].

      NB: [m] and [n] are always equal for band matrices. *)
  val size : t -> int * int

  (** Returns the dimensions of an array band matrix. *)
  val dims : t -> smu * mu * ml

  (** Pretty-print a band matrix using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
  val pp : Format.formatter -> t -> unit

  (** Pretty-print an array band matrix using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module.
      The defaults are: [start="\["], [stop="\]"], [sep=";"],
      [indent=4], [itemsep=" "], [empty="           ~           "] and
      [item=fun f r c->Format.fprintf f "(%2d,%2d)=% -15e" r c] (see
      {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}).
      The [indent] argument specifies the indent for wrapped rows. *)
  val ppi : ?start:string -> ?stop:string -> ?sep:string
            -> ?indent:int -> ?itemsep:string -> ?empty:string
            -> ?item:(Format.formatter -> int -> int -> float -> unit)
            -> unit
            -> Format.formatter -> t -> unit

  (** [get a i j] returns the value at row [i] and column [j] of [a].
      Only rows and columns satisfying
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid. *)
  val get : t -> int -> int -> float

  (** [set a i j v] sets the value at row [i] and column [j] of [a]
      to [v]. Only rows and columns satisfying
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid. *)
  val set : t -> int -> int -> float -> unit

  (** [update a i j f] sets the value at row [i] and column [j] of [a]
      to [f v]. Only rows and columns satisfying
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid. *)
  val update : t -> int -> int -> (float -> float) -> unit

  (** Direct access to the underlying storage array, which is accessed
      column first (unlike in {!get}). *)
  val unwrap : t -> RealArray2.data

  (** {3:arrayband_ops Operations} *)

  (** Operations on array-based band matrices. *)
  val ops : (t, RealArray.t) matrix_ops

  (** [scale_add c a b] calculates $A = cA + B$.

      NB: Unlike the {!Band.scale_add} operation, this operation raises an
      exception if [b] has a greater bandwidth than [a], i.e., it never
      resizes [a]. *)
  val scale_add : float -> t -> t -> unit

  (** [scale_addi ml c A] calculates $A = cA + I$. *)
  val scale_addi : float -> t -> unit

  (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

      @nodoc SUNDlsMat_bandMatvec *)
  val matvec : t -> RealArray.t -> RealArray.t -> unit

  (** Fills the matrix with zeros.

      @nodoc SUNDlsMat_SetToZero *)
  val set_to_zero : t -> unit

  (** [blit ~src ~dst] copies the contents of [src] into [dst].

      @nodoc SUNDlsMat_bandCopy *)
  val blit : src:t -> dst:t -> unit

  (** [lrw, liw = space a] returns the storage requirements of [a] as
      [lrw] realtype words and [liw] integer words. *)
  val space : t -> int * int

  (** {3:arrayband_calcs Calculations} *)

  (** Increment a square matrix by the identity matrix.

      @nodoc SUNDlsMat_bandAddIdentity *)
  val add_identity : t -> unit

  (** [scale c a] multiplies each element of the band matrix [a] by [c].

      @nodoc SUNDlsMat_bandScale *)
  val scale : float -> t -> unit

  (** [gbtrf a p] performs the LU factorization of [a] with partial pivoting
      according to [p]. The values in [a] are overwritten with those of the
      calculated L and U matrices. The diagonal belongs to U. The diagonal
      of L is all 1s. U may occupy elements up to bandwidth [smu]
      (rather than to [mu]).

      @nodoc SUNDlsMat_bandGBTRF *)
  val gbtrf : t -> LintArray.t -> unit

  (** [gbtrs a p b] finds the solution of [ax = b] using LU factorization.
      Both [p] and [b] must have the same number of rows as [a].

      @nodoc SUNDlsMat_bandGBTRS *)
  val gbtrs : t -> LintArray.t -> RealArray.t -> unit

end (* }}} *)

(** {2:generic Generic matrices} *)

(** Distinguishes a library-supplied matrix from a custom one. *)
type standard

(** Distinguishes a user-supplied matrix from a standard one. *)
type custom

(** A generic matrix with a payload of type ['m]. The ['k] type argument tracks
    whether the matrix is {!standard} or {!custom}. The ['nd] and ['nk] type
    arguments track the compatiblity of the {{!matrix_ops}m_matvec} vector
    parameters.

    @matrix SUNMatrix *)
type ('k, 'm, 'nd, 'nk) t

(** Generic matrix with Dense content. *)
type 'nk dense =
  (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

(** By default, [dense n] returns an [n] by [n] dense matrix with all elements
    initialized to [0.0]. Optional arguments allow specifying the number of rows
    ([m]) and the initial value ([i]).

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @matrix_data SUNDenseMatrix *)
val dense : ?context:Context.t -> ?m:int -> ?i:float -> int -> 'nk dense

(** Creates a (dense) matrix by wrapping an existing dense matrix. The two
    values share the same underlying storage.

    @matrix_data SUNDenseMatrix *)
val wrap_dense : ?context:Context.t -> Dense.t -> 'nk dense

(** Generic matrix with Band content. *)
type 'nk band =
  (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

(** By default, [band n] returns an [n] by [n] band matrix with all bandwidths
    equal to 2 and all values initialized to [0.0].
    Optional arguments allow specifying the upper bandwidth [mu], the lower
    bandwidth [ml], the storage upper bandwidth [smu] and the initial
    values [i]. If [mu] is given but not [ml], then [ml] is set to [mu].
    If [mu] is given but not [smu], then [smu] is set to [mu+ml].

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @matrix_data SUNBandMatrix
    @matrix_data SUNBandMatrixStorage *)
val band :
  ?context:Context.t
  -> ?mu:int
  -> ?smu:int
  -> ?ml:int
  -> ?i:float
  -> int
  -> 'nk band

(** Creates a (band) matrix by wrapping an existing band matrix. The two
    values share the same underlying storage.

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @matrix_data SUNBandMatrix *)
val wrap_band : ?context:Context.t -> Band.t -> 'nk band

(** Generic matrix with Sparse content. *)
type ('s, 'nk) sparse =
  (standard, 's Sparse.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

(** By default, [sparse_csc n] returns an [n] by [n] sparse matrix in
    {{!Sparse.sformat}CSC} format with the capacity for [n / 10] non-zero
    elements and all elements initialized to [0.0]. Optional arguments allow
    specifying the number of rows [m], and the number of non-zero
    elements [nnz].

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @matrix_data SUNSparseMatrix *)
val sparse_csc :
     ?context:Context.t
  -> ?m:int
  -> ?nnz:int
  -> int
  -> (Sparse.csc, 'nk) sparse

(** As for {!sparse_csc} but the returned matrix is in {{!Sparse.sformat}CSR}
    format.

    The {{!Sparse.sformat}CSR} format is only available from Sundials 2.7.0
    onwards.

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @matrix_data SUNSparseMatrix *)
val sparse_csr :
     ?context:Context.t
  -> ?m:int
  -> ?nnz:int
  -> int
  -> (Sparse.csr, 'nk) sparse

(** Creates a (sparse) matrix by wrapping an existing sparse matrix. The two
    values share the same underlying storage.

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @matrix_data SUNSparseMatrix *)
val wrap_sparse : ?context:Context.t -> 's Sparse.t -> ('s, 'nk) sparse

(** Generic matrix with array-based dense content. *)
type 'nk arraydense = (custom, ArrayDense.t, RealArray.t, 'nk) t

(** By default, [arraydense n] returns an [n] by [n] dense matrix with all
    elements initialized to [0.0]. Optional arguments allow specifying the
    number of rows ([m]) and the initial value ([i]).

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @nodoc SUNDlsMat_newDenseMat *)
val arraydense :
     ?context:Context.t
  -> ?m:int
  -> ?i:float
  -> int
  -> 'nk arraydense

(** Creates an (array-based dense) matrix by wrapping an existing array-based
    dense matrix. The two values share the same underlying storage.

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @nodoc SUNDlsMat_newDenseMat *)
val wrap_arraydense : ?context:Context.t -> ArrayDense.t -> 'nk arraydense

(** Generic matrix with array-based band content. *)
type 'nk arrayband = (custom, ArrayBand.t, RealArray.t, 'nk) t

(** By default, [band n] returns an [n] by [n] band matrix with all bandwidths
    equal to 2 and all values initialized to [0.0].
    Optional arguments allow specifying the upper bandwidth [mu], the lower
    bandwidth [ml], the storage upper bandwidth [smu] and the initial
    values [i]. If [mu] is given but not [smu], then [smu] is set to [mu],
    otherwise if [smu] is given but not [mu], then [mu] is set to [smu].
    If [ml] is not given, then it is set to [mu].

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @nodoc SUNDlsMat_newBandMat *)
val arrayband :
     ?context:Context.t
  -> ?mu:int
  -> ?smu:int
  -> ?ml:int
  -> ?i:float
  -> int
  -> 'nk arrayband

(** Creates an (array-based band) matrix by wrapping an existing array-based
    band matrix. The two values share the same underlying storage.

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument. *)
val wrap_arrayband :?context:Context.t ->  ArrayBand.t -> 'nk arrayband

(** Wrap a custom matrix value.

    By default, the matrix is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing
    an optional [context] argument.

    @matrix <SUNMatrix_API_link.html#description-of-the-sunmatrix-operations> Description of the SUNMATRIX operations *)
val wrap_custom :
     ('m, 'nd) matrix_ops
  -> ?context:Context.t
  -> 'm
  -> (custom, 'm, 'nd, 'nk) t

(** Matrix internal type identifiers.

    @matrix SUNMatGetID *)
type (_,_,_,_) id =
  | Dense : (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind]) id
  | Band  : (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind]) id
  | Sparse : (standard, 's Sparse.t, Nvector_serial.data, [>Nvector_serial.kind]) id
  | Custom : (custom, 'm, 'nd, 'nk) id
  | ArrayDense : (custom, ArrayDense.t, RealArray.t, 'nk) id
  | ArrayBand  : (custom, ArrayBand.t, RealArray.t, 'nk) id

(** Return a record of matrix operations. *)
val get_ops : ('k, 'm, 'nd, 'nk) t -> ('m, 'nd) matrix_ops

(** Return the internal type identifier of a matrix.

    @matrix SUNMatGetID *)
val get_id : ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) id

(** Direct access to the underlying storage array, which is accessed
    column first (unlike in {!Dense.get}, {!Band.get}, and {!Sparse.get}). *)
val unwrap : ('k, 'm, 'nd, 'nk) t -> 'm

(** {3:generic_ops Operations} *)

(** [scale_add c A B] calculates $A = cA + B$.

    @matrix SUNMatScaleAdd *)
val scale_add : float -> ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit

(** [scale_addi c A] calculates $A = cA + I$.

    @matrix SUNMatScaleAddI *)
val scale_addi : float -> ('k, 'm, 'nd, 'nk) t -> unit

(** Perform any setup required before a matrix-vector product.

    @matrix SUNMatMatvecSetup
    @since 5.0.0 *)
val matvec_setup : ('k, 'm, 'nd, 'nk) t -> unit

(** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

    @matrix SUNMatMatvec *)
val matvec :
  ('k, 'm, 'nd, 'nk) t -> ('nd, 'nk) Nvector.t
                       -> ('nd, 'nk) Nvector.t -> unit

(** Fills a matrix with zeros.

    @matrix SUNMatZero *)
val set_to_zero : ('k, 'm, 'nd, 'nk) t -> unit

(** [blit ~src ~dst] copies the contents of [src] into [dst]. Both
    must have the same size.

    @matrix SUNMatCopy *)
val blit : src:('k, 'm, 'nd, 'nk) t -> dst:('k, 'm, 'nd, 'nk) t -> unit

(** [lrw, liw = space a] returns the storage requirements of [a] as
    [lrw] realtype words and [liw] integer words.

    @matrix SUNMatSpace *)
val space : ('k, 'm, 'nd, 'nk) t -> int * int

(** Prints a dense matrix to the given log file.

    NB: Not supported in {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

    @matrix_data SUNDenseMatrix_Print *)
val print_dense : 'nk dense -> Logfile.t -> unit

(** Prints a band matrix to the given log file.

    NB: Not supported in {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

    @matrix_data SUNBandMatrix_Print *)
val print_band : 'nk band -> Logfile.t -> unit

(** Prints a sparse matrix to the given log file.

    NB: Not supported in {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

    @matrix_data SUNSparseMatrix_Print *)
val print_sparse : ('s, 'nk) sparse -> Logfile.t -> unit

(** Pretty-print a generic matrix using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module.
    For {{!id}Custom} matrices, it simply prints {e <custom matrix>}. *)
val pp : Format.formatter -> ('k, 'm, 'nd, 'nk) t -> unit

