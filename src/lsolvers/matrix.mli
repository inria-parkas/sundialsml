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

(** {2:shared Shared definitions} *)

(** Generic operations that all matrix types must implement. Failure is
    signalled by raising an exception.

    @nocvode <node> Description of the SUNMatrix module *)
type ('m, 'd, 'k) matrix_ops = { (* {{{ *)
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

  m_matvec    : 'm -> 'd -> 'd -> unit;
  (** [m_matvec a x y] calculates $y = ax$. *)

  m_space     : 'm -> int * int;
  (** [lrw, liw = m_space a] returns the number of realtype words [lrw] and
      integer words [liw] required to store the matrix [a]. *)
} (* }}} *)

(** Raised in {!Sundials.sundials_version} < 3.0.0 on an attempt to
    access a value that has become invalid. Such values refer to matrices
    that no longer exist in the underlying library. Values never
    become invalid in {!Sundials.sundials_version} >= 3.0.0. *)
exception Invalidated

(** Raised if matrix operation arguments are mutually incompatible. *)
exception IncompatibleArguments

(** Raised if a zero diagonal element is found during factorization using a
    low-level routine like {!Lsolver.Iterative.Algorithms.qr_fact},
    {!Lsolver.Iterative.Algorithms.qr_sol},
    or {!Matrix.ArrayDense.getrf}.
    The argument gives the equation number (from 1). *)
exception ZeroDiagonalElement of int

(** {2:content Matrix content} *)

(** Dense matrices *)
module Dense : (* {{{ *)
  sig
    (** A dense matrix. Values of this type are typically passed to linear
        solver callback functions (like {!Cvode.Dls.jac_fn},
        {!Ida.Dls.jac_fn}, and {!Kinsol.Dls.jac_fn}).

        @nocvode <node> The SUNMatrix_Dense implementation *)
    type t

    (** {3:basic Basic access} *)

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

    (** {3:ops Operations} *)

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

    (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

        @nocvode <node> SUNMatMatvec
        @nocvode <node> SUNMatMatvec_Dense *)
    val matvec : t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

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

    (** {3:lowlevel Low-level details} *)

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. Has no effect when
        {!Sundials.sundials_version} >= 3.0.0. *)
    val invalidate : t -> unit
  end (* }}} *)

(** Banded matrices *)
module Band : (* {{{ *)
  sig
    (** A band matrix. Values of this type are typically passed to linear
        solver callback functions (like {!Cvode.Dls.jac_fn},
        {!Ida.Dls.jac_fn}, and {!Kinsol.Dls.jac_fn}).

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

    (** {3:basic Basic access} *)

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

    (** {3:ops Operations} *)

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

    (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

        @nocvode <node> SUNMatMatvec
        @nocvode <node> SUNMatMatvec_Band *)
    val matvec : t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

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

    (** {3:lowlevel Low-level details} *)

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. Has no effect when
        {!Sundials.sundials_version} >= 3.0.0. *)
    val invalidate : t -> unit
  end (* }}} *)

(** Sparse matrices *)
module Sparse : (* {{{ *)
  sig
    type csc (** Compressed-sparse-column format. *)
    type csr (** Compressed-sparse-row format. *)

    (* Matrix storage formats. *)
    type _ sformat =
      | CSC : csc sformat (** Compressed-sparse-column format ([CSC_MAT]). *)
      | CSR : csr sformat (** Compressed-sparse-row format ([CSR_MAT]). *)

    (** A sparse matrix. Values of this type are typically passed to linear
        solver callback functions (like {!Cvode.Dls.jac_fn},
        {!Ida.Dls.jac_fn}, and {!Kinsol.Dls.jac_fn}).

        @nocvode <node> The SUNMatrix_Sparse implementation *)
    type 's t

    (** Array of row or column indices *)
    type index_array =
      (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t

    (** {3:basic Basic access} *)

    (** [make fmt m n nnz] returns an [m] by [n] sparse matrix in the
        specified format with a potential for [nnz] non-zero elements.
        All elements are initially zero.

        The {{!sformat}CSR} format is only available from Sundials 2.7.0
        onwards.

        @nocvode <node> SUNSparseMatrix *)
    val make : 's sformat -> int -> int -> int -> 's t

    (** Creates a sparse matrix in in the specified format from a dense matrix
        by copying all values of magnitude greater than the given tolerance.

        The {{!sformat}CSR} format is only available from Sundials 2.7.0
        onwards.

        @nocvode <node> SUNSparseFromDenseMatrix *)
    val from_dense : 's sformat -> float -> Dense.t -> 's t

    (** Creates a sparse matrix in the specified format from a band matrix by
        copying all values of magnitude greater than the given tolerance.

        The {{!sformat}CSR} format is only available from Sundials 2.7.0
        onwards.

        @nocvode <node> SUNSparseFromBandMatrix *)
    val from_band : 's sformat -> float -> Band.t -> 's t

    (** Return the matrix format. *)
    val sformat : 's t -> 's sformat

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

    (** {3:ops Operations} *)

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

    (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

        @nocvode <node> SUNMatMatvec
        @nocvode <node> SUNMatMatvec_Sparse *)
    val matvec : 's t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

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

    (** {3:lowlevel Low-level details} *)

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
  end (* }}} *)

(** {2:array Arrays as matrices} *)

(** General purpose dense matrix operations on arrays.
    @cvode <node9#ss:dense> The DENSE Module *)
module ArrayDense : sig (* {{{ *)

  (** A dense matrix accessible directly through a
      {{:OCAML_DOC_ROOT(Bigarray.Array2.html)} Bigarray}.

      @cvode <node9#ss:dense> Small dense matrices
      @cvode <node9#ss:dense> newDenseMat *)
  type t = Sundials.RealArray2.t

  (** {4 Basic access} *)

  (** [make m n x] returns an [m] by [n] dense matrix with elements set
      to [x].

      @cvode <node9#ss:dense> newDenseMat *)
  val make : int -> int -> float -> t

  (** [create m n] returns an uninitialized [m] by [n] dense matrix.

       @cvode <node9#ss:dense> newDenseMat *)
  val create : int -> int -> t

  (** [get a i j] returns the value at row [i] and column [j] of [a]. *)
  val get : t -> int -> int -> float

  (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v]. *)
  val set : t -> int -> int -> float -> unit

  (** [update a i j f] sets the value at row [i] and column [j] of [a]
      to [f v]. *)
  val update : t -> int -> int -> (float -> float) -> unit

  (** Fills the matrix with zeros.

      @cvode <node9#ss:dense> setToZero *)
  val set_to_zero    : t -> unit

  (** {4 Calculations} *)

  (** Increments a square matrix by the identity matrix.

      @cvode <node9#ss:dense> denseAddIdentity *)
  val add_identity : t -> unit

  (** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

      @nocvode <node9#ss:dense> denseMatvec
      @since 2.6.0 *)
  val matvec : t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

  (** [blit src dst] copies the contents of [src] into [dst]. Both
      must have the same size.

      @cvode <node9#ss:dense> denseCopy *)
  val blit  : t -> t -> unit

  (** Multiplies each element by a constant.

      @cvode <node9#ss:dense> denseScale *)
  val scale : float -> t -> unit

  (** [getrf a p] performs the LU factorization of the square matrix [a] with
      partial pivoting according to [p]. The values in [a] are overwritten
      with those of the calculated L and U matrices. The diagonal belongs to
      U. The diagonal of L is all 1s. Multiplying L by U gives a permutation
      of [a], according to the values of [p]: [p.{k} = j] means that rows [k]
      and [j] were swapped (in order, where [p.{0}] swaps against the
      original matrix [a]).

      @cvode <node9#ss:dense> denseGETRF
      @raise ZeroDiagonalElement Zero found in matrix diagonal *)
  val getrf : t -> Sundials.LintArray.t -> unit

  (** [getrs a p b] finds the solution of [ax = b] using an LU factorization
      found by {!getrf}. Both [p] and [b] must have the same number of rows
      as [a].

      @cvode <node9#ss:dense> denseGETRS *)
  val getrs : t -> Sundials.LintArray.t -> Sundials.RealArray.t -> unit

  (** Like {!getrs} but stores [b] starting at a given offset. *)
  val getrs'
        : t -> Sundials.LintArray.t -> Sundials.RealArray.t -> int -> unit

  (** Performs Cholesky factorization of a real symmetric positive matrix.

      @cvode <node9#ss:dense> densePOTRF *)
  val potrf : t -> unit

  (** [potrs a b] finds the solution of [ax = b] using the Cholesky
      factorization found by {!potrf}. [a] must be an n by n matrix and [b]
      must be of length n.

      @cvode <node9#ss:dense> densePOTRS *)
  val potrs : t -> Sundials.RealArray.t -> unit

  (** [geqrf a beta work] performs the QR factorization of [a]. [a] must be
      an [m] by [n] matrix, where [m >= n]. The [beta] vector must have
      length [n]. The [work] vector must have length [m].

      @cvode <node9#ss:dense> denseGEQRF *)
  val geqrf : t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

  (** [ormqr q beta v w work] computes the product {% w = qv %}. [Q] is
      an [m] by [n] matrix calculated using {!geqrf} with [m >= n],
      [beta] has length [n], [v] has length [n], [w] has length [m], and
      [work] has length [m].

      @param q       matrix modified by {!geqrf}
      @param beta    vector passed to {!geqrf}
      @param v       vector multiplier
      @param w       result vector
      @param work    temporary vector used in the calculation
      @cvode <node9#ss:dense> denseORMQR *)
  val ormqr :
    a:t -> beta:Sundials.RealArray.t -> v:Sundials.RealArray.t
      -> w:Sundials.RealArray.t -> work:Sundials.RealArray.t -> unit

end (* }}} *)

(** General-purpose band matrix operations on arrays.
    @cvode <node9#ss:band> The BAND Module *)
module ArrayBand : sig (* {{{ *)

  (** A band matrix accessible directly through a
      {{:OCAML_DOC_ROOT(Bigarray.Array2.html)} Bigarray}.

      The layout of these arrays are characterized by the {e storage upper
      bandwidth} [smu]. Given an array [a], the first dimension indexes
      the diagonals, with the main diagonal ({% $i = j$ %}) at [a.{smu, _}].
      The value in the [i]th row and [j]th column provided
      {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} is at
      [a.{i - j + smu, j}].

      @cvode <node9#ss:band> newBandMat *)
  type t = Sundials.RealArray2.t

  type smu = int (** Storage upper-bandwidth. *)

  type mu = int  (** Upper-bandwidth. *)
  type ml = int  (** Lower-bandwidth. *)

  (** {4 Basic access} *)

  (** [create n smu ml] returns an [n] by [n] band matrix with
      storage upper bandwidth [smu] and lower half-bandwidth [ml].
      If the result will not be LU factored then
      {% $\mathtt{smu} = \mathtt{mu}$ %}, otherwise
      {% $\mathtt{smu} = \min(\mathtt{n}-1, \mathtt{mu} + \mathtt{ml})$ %}.
      The extra space is used to store U after a call to {!gbtrf}.

      @cvode <node9#ss:band> newBandMat *)
  val create : int -> smu -> ml -> t

  (** [get a smu i j] returns the value at row [i] and column [j] of [a].
      [smu] is the storage upper bandwidth. Only rows and columns
      satisfying {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid. *)
  val get : t -> smu -> int -> int -> float

  (** [set a smu i j v] sets the value at row [i] and column [j] of [a]
      to [v]. [smu] is the storage upper bandwidth. Only rows and columns
      satisfying {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid. *)
  val set : t -> smu -> int -> int -> float -> unit

  (** [update a smu i j f] sets the value at row [i] and column [j] of [a]
      to [f v]. [smu] is the storage upper bandwidth. Only rows and columns
      satisfying {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
      {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid. *)
  val update : t -> smu -> int -> int -> (float -> float) -> unit

  (** {4 Calculations} *)

  (** Increment a square matrix by the identity matrix.

      @cvode <node9#ss:band> bandAddIdentity *)
  val add_identity : t -> smu -> unit

  (** The call [matvec a smu mu ml x y] computes the matrix-vector product
      $y = Ax$.

      @nocvode <node9#ss:band> bandMatvec
      @since 2.6.0 *)
  val matvec : t -> smu -> mu -> ml
                -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

  (** [blit src dst src_smu dst_smu copy_mu copy_ml] copies the contents
      of [src] into [dst]. The storage upper bandwidths of [src] and [dst]
      are, respectively, [src_smu] and [dst_smu]. The bandwidth to copy is
      given by [copy_mu] and [copy_ml]. Both matrices must have the
      same size.

      @cvode <node9#ss:band> bandCopy *)
  val blit : t -> t -> smu -> smu -> int -> int -> unit

  (** [scale c a smu mu ml] multiplies each element of the band matrix
      [a] by [c].
      The parameters [smu], [mu], and [ml] give, respectively,
      the storage upper, upper, and lower bandwidths.

      @cvode <node9#ss:band> bandScale *)
  val scale : float -> t -> smu -> mu -> ml -> unit

  (** [gbtrf a smu mu ml p] performs the LU factorization of [a] with
      partial pivoting according to [p]. The parameters [smu], [mu], and
      [ml] give, respectively, the storage upper, upper, and lower
      bandwidths. The values in [a] are overwritten with those of the
      calculated L and U matrices. The diagonal belongs to U. The diagonal
      of L is all 1s. U may occupy elements up to bandwidth [smu]
      (rather than to [mu]).

      @cvode <node9#ss:band> bandGBTRF *)
  val gbtrf : t -> smu -> mu -> ml -> Sundials.LintArray.t -> unit

  (** [gbtrs a smu ml p b] finds the solution of [ax = b] using LU
      factorization. The parameters [smu] and [ml] give, respectively, the
      storage upper and lower bandwidths. Both [p] and [b] must have the same
      number of rows as [a].

      @cvode <node9#ss:band> bandGBTRS *)
  val gbtrs : t -> smu -> ml -> Sundials.LintArray.t ->
    Sundials.RealArray.t -> unit

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

    @nocvode <node> SUNMatrix *)
type ('k, 'm, 'nd, 'nk) t

(** Generic matrix with Dense content. *)
type 'nk dense =
  (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

(** Generic matrix with Band content. *)
type 'nk band =
  (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

(** Generic matrix with Sparse content. *)
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
    {!Band.dimensions} and all elements initialized to [x].

    @nocvode <node> SUNBandMatrix *)
val make_band : Band.dimensions -> float -> 'nk band

(** Creates a (band) matrix by wrapping an existing band matrix. The two
    values share the same underlying storage.

    @nocvode <node> SUNBandMatrix *)
val wrap_band : Band.t -> 'nk band

(** [make m n nnz] returns an [m] by [n] matrix in the specified format with
    a potential for [nnz] non-zero elements. All elements are initially zero.

    The {{!Sparse.sformat}CSR} format is only available from Sundials 2.7.0
    onwards.

    @nocvode <node> SUNSparseMatrix *)
val make_sparse : 's Sparse.sformat -> int -> int -> int -> ('s, 'nk) sparse

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
    column first (unlike in {!Dense.get}, {!Band.get}, and {!Sparse.get}).

    @nocvode <node> SM_CONTENT_B *)
val unwrap : ('k, 'm, 'nd, 'nk) t -> 'm

(** {3 Operations} *)

(** [scale_add c A B] calculates $A = cA + B$.

    @nocvode <node> SUNMatScaleAdd *)
val scale_add : float -> ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit

(** [scale_addi c A] calculates $A = cA + I$.

    @nocvode <node> SUNMatScaleAddI *)
val scale_addi : float -> ('k, 'm, 'nd, 'nk) t -> unit

(** The call [matvec a x y] computes the matrix-vector product $y = Ax$.

    @nocvode <node> SUNMatMatvec *)
val matvec :
  ('k, 'm, 'nd, 'nk) t -> ('nd, 'nk) Nvector.t
                       -> ('nd, 'nk) Nvector.t -> unit

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

