(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(** Direct Linear Solvers.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)

    @cvode <node9#s:dls>  The DLS Modules *)

(** Raised by {!DenseMatrix.getrf} and {!ArrayDenseMatrix.getrf} on a zero
    diagonal element during factorization. The argument gives the column
    index (from 1).

    @cvode <node9#ss:dense> DenseGETRF/denseGETRF *)
exception ZeroDiagonalElement of int

(** {2 Dense matrices} *)

(** Dense matrices as passed to callback functions.
    @cvode <node9#ss:dense> The DENSE Module *)
module DenseMatrix :
  sig
    (** An abstract dense matrix. Values of this type are typically passed
        to linear solver callback functions (like {!Cvode.Dls.dense_jac_fn}, 
        {!Ida.Dls.dense_jac_fn}, and {!Kinsol.Dls.dense_jac_fn}), in which
        case their lifetimes are determined by the underlying library and they
        should only be used within the callback to avoid the
        {{!Invalidated}Invalidated} exception.

        @cvode <node9#s:dls>  Type DlsMat *)
    type t

    (** {4 Basic access} *)

    (** [make m n x] returns an [m] by [n] dense matrix with elements set
        to [x].

        @cvode <node9#ss:dense> NewDenseMat *)
    val make : int -> int -> float -> t

    (** [create m n] returns an uninitialized [m] by [n] dense matrix.

        @cvode <node9#ss:dense> NewDenseMat *)
    val create : int -> int -> t

    (** [m, n = size a] returns the numbers of columns [m] and rows [n]
        of [a]. *)
    val size  : t -> int * int

    (** Prints a dense matrix to stdout.

        @cvode <node9#ss:dense> PrintMat *)
    val print : t -> unit

    (** [get a i j] returns the value at row [i] and column [j] of [a].

        @cvode <node9#s:dls> DENSE_ELEM *)
    val get : t -> int -> int -> float

    (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v].

        @cvode <node9#s:dls> DENSE_ELEM *)
    val set : t -> int -> int -> float -> unit

    (** Fills a matrix with zeros.

        @cvode <node9#ss:dense> SetToZero *)
    val set_to_zero    : t -> unit

    (** {4 Calculations} *)

    (** Increments a square matrix by the identity matrix.

        @cvode <node9#ss:dense> AddIdentity *)
    val add_identity   : t -> unit

    (** [blit src dst] copies the contents of [src] into [dst]. Both
        must have the same size.

        @cvode <node9#ss:dense> DenseCopy *)
    val blit : t -> t -> unit

    (** Multiplies each element by a constant.

        @cvode <node9#ss:dense> DenseScale *)
    val scale    : float -> t -> unit

    (** [getrf a p] performs the LU factorization of the square matrix [a] with
        partial pivoting according to [p]. The values in [a] are overwritten
        with those of the calculated L and U matrices. The diagonal belongs to
        U. The diagonal of L is all 1s. Multiplying L by U gives a permutation
        of [a], according to the values of [p]: [p.{k} = j] means that rows [k]
        and [j] were swapped (in order, where [p.{0}] swaps against the
        original matrix [a]).

        @cvode <node9#ss:dense> DenseGETRF
        @raise ZeroDiagonalElement Zero found in matrix diagonal *)
    val getrf    : t -> Sundials.LintArray.t -> unit

    (** [getrs a p b] finds the solution of [ax = b] using an LU factorization
        found by {!getrf}. Both [p] and [b] must have the same number of
        rows as [a].

        @cvode <node9#ss:dense> DenseGETRS *)
    val getrs    : t -> Sundials.LintArray.t -> Sundials.RealArray.t -> unit

    (** Performs Cholesky factorization of a real symmetric positive matrix.

        @cvode <node9#ss:dense> DensePOTRF *)
    val potrf    : t -> unit

    (** [potrs a b] finds the solution of [ax = b] using the Cholesky
        factorization found by {!potrf}. [a] must be an n by n matrix and [b]
        must be of length n.

        @cvode <node9#ss:dense> DensePOTRS *)
    val potrs    : t -> Sundials.RealArray.t -> unit

    (** [geqrf a beta work] performs the QR factorization of [a]. [a] must be
        an [m] by [n] matrix, where [m >= n]. The [beta] vector must have
        length [n]. The [work] vector must have length [m].

        @cvode <node9#ss:dense> DenseGEQRF *)
    val geqrf    : t -> Sundials.RealArray.t -> Sundials.RealArray.t -> unit

    (** [ormqr q beta v w work] computes the product {% w = qv %}. [Q] is
        an [m] by [n] matrix calculated using {!geqrf} with [m >= n],
        [beta] has length [n], [v] has length [n], [w] has length [m], and
        [work] has length [m].

        @param q       matrix modified by {!geqrf}
        @param beta    vector passed to {!geqrf}
        @param v       vector multiplier
        @param w       result vector
        @param work    temporary vector used in the calculation
        @cvode <node9#ss:dense> DenseORMQR *)
    val ormqr :
      a:t -> beta:Sundials.RealArray.t -> v:Sundials.RealArray.t
        -> w:Sundials.RealArray.t -> work:Sundials.RealArray.t -> unit

    (** {4 Low-level details} *)

    (** Raised on an attempt to access a value that has become invalid. Such
        values refer to matrices that no longer exist in the underlying
        library. *)
    exception Invalidated

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. *)
    val invalidate : t -> unit

    (** Potentially unsafe access to the underlying storage. This array
        {b must} only be used when the underlying storage is valid, which
        will be the case in callbacks. The array is accessed column first
        (unlike in {!get}). *)
    val unsafe_unwrap
      : t -> (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
  end

(** General purpose dense matrix operations on arrays.
    @cvode <node9#ss:dense> The DENSE Module *)
module ArrayDenseMatrix :
  sig
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

    (** Fills the matrix with zeros.
     
        @cvode <node9#ss:dense> setToZero *)
    val set_to_zero    : t -> unit

    (** {4 Calculations} *)

    (** Increments a square matrix by the identity matrix.

        @cvode <node9#ss:dense> denseAddIdentity *)
    val add_identity : t -> unit

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
  end

(** {2 Band matrices} *)

(** Band matrices as passed to callback functions.
    @cvode <node9#ss:band> The BAND Module *)
module BandMatrix :
  sig
    (** An abstract band matrix. Values of this type are typically passed
        to linear solver callback functions (like {!Cvode.Dls.band_jac_fn}, 
        {!Ida.Dls.band_jac_fn}, and {!Kinsol.Dls.band_jac_fn}), in which case
        their lifetimes are determined by the underlying library and they
        should only be used within the callback to avoid the
        {{!Invalidated}Invalidated} exception.

        @cvode <node9#s:dls>  Type DlsMat *)
    type t

    (** Band matrix dimensions. If the result will not be LU factored then
        {% $\mathtt{smu} = \mathtt{mu}$ %}, otherwise
        {% $\mathtt{smu} = \min(\mathtt{n}-1, \mathtt{mu} + \mathtt{ml})$ %}.
        The extra space is used to store U after a call to {!gbtrf}. *)
    type dimensions = {
        n   : int;  (** Matrix size: [n] by [n]. *)
        mu  : int;  (** Upper bandwidth. *)
        smu : int;  (** Storage upper bandwidth. *)
        ml  : int;  (** Lower bandwidth. *)
      }

    (** {4 Basic access} *)

    (** Returns a band matrix with the given {!dimensions} and all elements
        initialized to the given value.

        @cvode <node9#ss:band> NewBandMat *)
    val make : dimensions -> float -> t

    (** Returns an uninitialized band matrix with the given {!dimensions}.

        @cvode <node9#ss:band> NewBandMat *)
    val create : dimensions -> t

    (** Returns the dimensions of a band matrix. *)
    val size  : t -> dimensions

    (** Prints a band matrix to stdout.

        @cvode <node9#ss:band> PrintMat *)
    val print : t -> unit

    (** [get a i j] returns the value at row [i] and column [j] of [a].
        Only rows and columns satisfying
        {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
        {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

        @cvode <node9#s:dls> BAND_ELEM *)
    val get : t -> int -> int -> float

    (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v].
        Only rows and columns satisfying
        {% $\mathtt{i} \leq \mathtt{j} + \mathtt{ml}$ %} and
        {% $\mathtt{j} \leq \mathtt{i} + \mathtt{smu}$ %} are valid.

        @cvode <node9#s:dls> BAND_ELEM *)
    val set : t -> int -> int -> float -> unit

    (** {4 Calculations} *)

    (** Fills a matrix with zeros.

        @cvode <node9#ss:band> SetToZero *)
    val set_to_zero    : t -> unit

    (** Increment a square matrix by the identity matrix.

        @cvode <node9#ss:band> AddIdentity *)
    val add_identity   : t -> unit

    (** [blit src dst copymu copyml] copies the contents of [src] into [dst].
        The bandwidth to copy is given by [copy_mu] and [copy_ml]. Both
        matrices must have the same size.

        @cvode <node9#ss:band> BandCopy *)
    val blit : t -> t -> int -> int -> unit

    (** Multiplies each element by a constant.

        @cvode <node9#ss:band> BandScale *)
    val scale : float -> t -> unit

    (** [gbtrf a p] performs the LU factorization of [a] with partial pivoting
        according to [p]. The values in [a] are overwritten with those of the
        calculated L and U matrices. The diagonal belongs to U. The diagonal
        of L is all 1s. U may occupy elements up to bandwidth [smu] (rather
        than to [mu]).

        @cvode <node9#ss:band> BandGBTRF *)
    val gbtrf : t -> Sundials.LintArray.t -> unit

    (** [gbtrs a p b] finds the solution of [ax = b] using an LU factorization
        found by {!gbtrf}. Both [p] and [b] must have the same number of rows
        as [a].

        @cvode <node9#ss:band> BandGBTRS *)
    val gbtrs : t -> Sundials.LintArray.t -> Sundials.RealArray.t -> unit

    (** {4 Low-level details} *)

    (** Raised on an attempt to access a value that has become invalid. Such
        values refer to matrices that no longer exist in the underlying
        library. *)
    exception Invalidated

    (** Called internally when the corresponding value in the underlying
        library ceases to exist. *)
    val invalidate : t -> unit

    (** Potentially unsafe access to the underlying storage. This array
        {b must} only be used when the underlying storage is valid, which
        will be the case in callbacks. The array is accessed column first
        (unlike in {!get}). *)
    val unsafe_unwrap
      : t -> (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
  end

(** General-purpose band matrix operations on arrays.
    @cvode <node9#ss:band> The BAND Module *)
module ArrayBandMatrix :
  sig
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

    type smu = int
    type mu = int
    type ml = int

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

    (** {4 Calculations} *)

    (** Increment a square matrix by the identity matrix.

        @cvode <node9#ss:band> bandAddIdentity *)
    val add_identity : t -> smu -> unit


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
  end

