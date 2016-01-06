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
    (** A sparse matrix. Values of this type are passed
        to linear solver callback functions, in which
        case their lifetimes are determined by the underlying library and they
        should only be used within the callback to avoid the
        {{!Invalidated}Invalidated} exception.

        @nocvode <node9#s:sls>  Type SlsMat *)
    type t

    (** {4 Basic access} *)

    (** [make m n nnz] returns an [m] by [n] sparse matrix with a potential
        for [nnz] non-zero elements. All elements are initially zero.

        @nocvode <node9#ss:sparse> NewSparseMat *)
    val make : int -> int -> int -> t

    (** [create m n nnz] returns an uninitialized [m] by [n] sparse matrix
        with a potential for [nnz] non-zero elements.

        @cvode <node9#ss:dense> NewDenseMat *)
    val create : int -> int -> int -> t

    (** [m, n, nnz = size a] returns the numbers of columns [m], rows [n], and
        the maximum number of non-zero elements of [a]. *)
    val size  : t -> int * int * int

    (** Prints a sparse matrix to stdout.

        @nocvode <node9#ss:sparse> PrintSparseMat *)
    val print : t -> unit

    (** Create a sparse matrix from a dense one.

        @nocvode <node9#ss:sparse> SlsConvertDls *)
    val from_dense : Dls.DenseMatrix.t -> t

    (** Create a sparse matrix from a banded one.

        @nocvode <node9#ss:sparse> SlsConvertDls *)
    val from_band : Dls.BandMatrix.t -> t

    (** [set_col a j idx] sets the data index of column [j] to [idx]. *)
    val set_col : t -> int -> int -> unit

    (** [get_col a j] returns the data index of column [j]. *)
    val get_col : t -> int -> int

    (** [set a idx i v] sets the [idx]th row to [i] and its value to [v]. *)
    val set : t -> int -> int -> float -> unit

    (** [r, v = get a idx] returns the row [r] and value [v] at the [idx]th
        position. *)
    val get : t -> int -> int * float

    (** Fills a matrix with zeros.

        @nocvode <node9#ss:sparse> SlsSetToZero *)
    val set_to_zero    : t -> unit

    (** Reallocates enoughs space for the given number of non-zero values. *)
    val realloc        : t -> int -> unit

    (** {4 Calculations} *)

    (** Increments a square matrix by the identity matrix.

        @nocvode <node9#ss:sparse> AddIdentitySparseMat *)
    val add_identity   : t -> unit

    (** [blit src dst] copies the contents of [src] into [dst].

        @nocvode <node9#ss:sparse> CopySparseMat *)
    val blit : t -> t -> unit

    (** Multiplies each element by a constant.

        @nocvode <node9#ss:sparse> ScaleSparseMat *)
    val scale    : float -> t -> unit

    (** Multiplies each element by a constant.

        @nocvode <node9#ss:sparse> SlsAddMat *)
    val add    : t -> t -> unit

    (** [matvec a x y] computes the matrix-vector product [y = A*x].
        If [a] has dimensions [m] by [n], then [x] must be of size [n],
        and [y] must be of size [m]. *)
    val matvec : t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

    (** {4 Low-level details} *)

    (** [set_rowval a idx i] sets the [idx]th row to [i]. *)
    val set_rowval : t -> int -> int -> unit

    (** [r = get_rowval a idx] returns the row [r] at the [idx]th position. *)
    val get_rowval : t -> int -> int

    (** [set_data a idx v] sets the value of the [idx]th row [v]. *)
    val set_data : t -> int -> float -> unit

    (** [v = get_data a idx] returns the value [v] at the [idx]th position. *)
    val get_data : t -> int -> float

    (** Raised on an attempt to access a value that has become invalid. Such
        values refer to matrices that no longer exist in the underlying
        library. *)
    exception Invalidated

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. *)
    val invalidate : t -> unit
  end

