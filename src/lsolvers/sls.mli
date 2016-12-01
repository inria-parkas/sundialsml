(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Sparse Linear Solvers.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node9#s:sls>  The SLS Module *)

(** Sparse matrices as passed to callback functions.
    @nocvode <node9#ss:sparse> The SLS Module *)
module SparseMatrix :
  sig
    type csc  (* Compressed-sparse-column format ([CSC_MAT]). *)
    type csr  (* Compressed-sparse-row format ([CSR_MAT]). *)

    (** A sparse matrix. Values of this type are passed
        to linear solver callback functions, in which
        case their lifetimes are determined by the underlying library and they
        should only be used within the callback to avoid the
        {{!Invalidated}Invalidated} exception.

        The type argument ['sformat] specifies the storage format.
        It is either {csc} or {csr}.

        @nocvode <node9#s:sls>  Type SlsMat *)
    type 'sformat t = 'sformat Sls_impl.t
    type t_csc = csc t
    type t_csr = csr t

    (** {4 Basic access} *)

    (** [make m n nnz] returns an [m] by [n] sparse matrix in
        compressed-sparse-column format with a potential
        for [nnz] non-zero elements. All elements are initially zero.

        @nocvode <node9#ss:sparse> SparseNewMat *)
    val make_csc : int -> int -> int -> csc t

    (** [create m n nnz] returns an uninitialized [m] by [n] sparse matrix
        in compressed-sparse-column format with a potential for [nnz]
        non-zero elements.

        @cvode <node9#ss:sparse> SparseNewMat *)
    val create_csc : int -> int -> int -> csc t

    (** [make m n nnz] returns an [m] by [n] sparse matrix in
        compressed-sparse-row format with a potential
        for [nnz] non-zero elements. All elements are initially zero.

        @since 2.7.0
        @raise Sundials.NotImplementedBySundialsVersion CSR format not available.
        @nocvode <node9#ss:sparse> SparseNewMat *)
    val make_csr : int -> int -> int -> csr t

    (** [create m n nnz] returns an uninitialized [m] by [n] sparse matrix in
        compressed-sparse-row format with a potential for [nnz] non-zero
        elements.

        @since 2.7.0
        @raise Sundials.NotImplementedBySundialsVersion CSR format not available.
        @cvode <node9#ss:sparse> SparseNewMat *)
    val create_csr : int -> int -> int -> csr t

    (** [m, n, nnz = size a] returns the numbers of columns [m], rows [n], and
        the maximum number of non-zero elements of [a]. *)
    val size  : 'f t -> int * int * int

    (** Prints a sparse matrix to stdout.

        In versions of Sundials prior to 2.7.0, the only valid log file
        is {Sundials.Logfile.stdout}.

        @raise Sundials.NotImplementedBySundialsVersion log file not supported.
        @nocvode <node9#ss:sparse> SparsePrintMat *)
    val print : Sundials.Logfile.t -> 'f t -> unit

    (* TOPLEVEL-PRINTER: Sls.SparseMatrix.pp *)
    (** Pretty-print a sparse matrix using the
        {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
    val pp : Format.formatter -> 'f t -> unit

    (** Create a compressed-sparse-column matrix from a dense one.

        @nocvode <node9#ss:sparse> SparseFromDenseMat *)
    val csc_from_dense : Dls.DenseMatrix.t -> csc t

    (** Create a compressed-sparse-row matrix from a dense one.

        @nocvode <node9#ss:sparse> SparseFromDenseMat *)
    val csr_from_dense : Dls.DenseMatrix.t -> csr t

    (** Create a sparse matrix from a banded one.

        @nocvode <node9#ss:sparse> SparseFromDenseMat *)
    val csc_from_band : Dls.BandMatrix.t -> csc t

    (** Create a sparse matrix from a banded one.

        @nocvode <node9#ss:sparse> SparseFromDenseMat *)
    val csr_from_band : Dls.BandMatrix.t -> csr t

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

    (** Fills a matrix with zeros.

        @nocvode <node9#ss:sparse> SparseSetMatToZero *)
    val set_to_zero    : 'f t -> unit

    (** Reallocates enoughs space for the given number of non-zero values. *)
    val realloc        : 'f t -> int -> unit

    (** {4 Calculations} *)

    (** Increments a square matrix by the identity matrix.

        @nocvode <node9#ss:sparse> SparseAddIdentityMat *)
    val add_identity   : 'f t -> unit

    (** [blit src dst] copies the contents of [src] into [dst].

        @nocvode <node9#ss:sparse> SparseCopyMat *)
    val blit : 'f t -> 'f t -> unit

    (** Multiplies each element by a constant.

        @nocvode <node9#ss:sparse> SparseScaleMat *)
    val scale    : float -> 'f t -> unit

    (** Adds two matrices.

        @nocvode <node9#ss:sparse> SparseAddMat *)
    val add    : 'f t -> 'f t -> unit

    (** [matvec a x y] computes the matrix-vector product [y = A*x].
        If [a] has dimensions [m] by [n], then [x] must be of size [n],
        and [y] must be of size [m].
        
        @nocvode <node9#ss:sparse> SparseMatvec *)
    val matvec : 'f t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

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

    (** Raised on an attempt to access a value that has become invalid. Such
        values refer to matrices that no longer exist in the underlying
        library. *)
    exception Invalidated

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. *)
    val invalidate : 'f t -> unit
  end

