(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(**
  {2 Direct Linear Solvers.}
  @cvode <node9#s:dls>  The DLS Modules
 *)

(**
 Thrown by the getrf functions if a zero diagonal element is encountered during
 factorization. The argument indicates the column index (from 1).

 @cvode <node9#ss:dense> DenseGETRF/denseGETRF 
 *)
exception ZeroDiagonalElement of int

(** {3 Dense matrices}
    @cvode <node9#ss:dense> The DENSE Module *)

(** Operations for creating and manipulating dense matrices. *)
module DenseMatrix :
  sig
    (**
    This type represents a [DlsMat] returned from a call to
    {!make}.

     @cvode <node9#s:dls>  Type DlsMat
     @cvode <node9#ss:dense> NewDenseMat 
     *)
    type t

    (** {4 Basic access} *)

    (**
     [make m n] returns an [m] by [n]  dense matrix.

     @cvode <node9#ss:dense> NewDenseMat
     *)
    val make  : int -> int -> t

    (**
     [m, n = size a] returns the number of columns ([m]) and rows ([n]) of the
     dense matrix [a].
     *)
    val size  : t -> (int * int)

    (**
     Prints a dense matrix to stdout.

     @cvode <node9#ss:dense> PrintMat
     *)
    val print : t -> unit

    (**
     [get a i j] returns the value at row [i] and column [j] in [a],
     where [0 <= i < m] and [0 <= j < n].

     @cvode <node9#s:dls> DENSE_ELEM
     *)
    val get : t -> int -> int -> float

    (**
     [set a i j v] stores the value [v] at row [i] and column [j] in [a],
     where [0 <= i < m] and [0 <= j < n].

     @cvode <node9#s:dls> DENSE_ELEM
     *)
    val set : t -> int -> int -> float -> unit

    (** {4 Calculations} *)

    (**
     Fills the matrix with zeros.

     @cvode <node9#ss:dense> SetToZero
     *)
    val set_to_zero    : t -> unit

    (**
     Increment a square matrix by the identity matrix.

     @cvode <node9#ss:dense> AddIdentity
     *)
    val add_identity   : t -> unit

    (**
     [copy src dst] copies the contents of one matrix into another.

     @cvode <node9#ss:dense> DenseCopy
     *)
    val copy     : t -> t -> unit

    (**
     [scale c a] multiplies each element of [a] by [c].

     @cvode <node9#ss:dense> DenseScale
     *)
    val scale    : float -> t -> unit

    (**
     [getrf a p] performs the LU factorization of [a] with partial pivoting
     according to [p].

     @cvode <node9#ss:dense> DenseGETRF
     @raise ZeroDiagonalElement Zero found in matrix diagonal
     *)
    val getrf    : t -> Sundials.lint_array -> unit

    (**
     [getrs a p b] finds the solution of [ax = b] using LU factorization.
     [a] must be a square matrix.

     @cvode <node9#ss:dense> DenseGETRS
     *)
    val getrs    : t -> Sundials.lint_array -> Sundials.real_array -> unit

    (**
     Performs Cholesky factorization of a real symmetric positive matrix.

     @cvode <node9#ss:dense> DensePOTRF
     *)
    val potrf    : t -> unit

    (**
     [potrs a b] finds the solution of [ax = b] using Cholesky factorization.

     @cvode <node9#ss:dense> DensePOTRS
     *)
    val potrs    : t -> Sundials.real_array -> unit

    (**
     [geqrf a beta work] performs the QR factorization of m by n matrix, with
     m >= n.

     @cvode <node9#ss:dense> DenseGEQRF
     *)
    val geqrf    : t -> Sundials.real_array -> Sundials.real_array -> unit

    (**
     [ormqr a beta v w work] computes the product [w = Qv], with Q calculated using {!geqrf}.

     @param a       matrix passed to {!geqrf}
     @param beta    vector apssed to {!geqrf}
     @param v       vector multiplier
     @param w       result vector
     @param work    temporary vector used in the calculation
     @cvode <node9#ss:dense> DenseORMQR
     *)
    val ormqr :
      a:t -> beta:Sundials.real_array -> v:Sundials.real_array ->
      w:Sundials.real_array -> work:Sundials.real_array -> unit

  end

(** {3 Array-based dense matrices}
    @cvode <node9#ss:dense> The DENSE Module *)

(** Operations for creating and manipulating dense matrices based on arrays. *)
module ArrayDenseMatrix :
  sig
    (**
     This type represents a two-dimensional matrix returned from a call to
     {!make}, {!Sundials.Realarray2.make}, or {!Sundials.Realarray2.wrap}.

     @cvode <node9#ss:dense> Small dense matrices
     @cvode <node9#ss:dense> newDenseMat 
     *)
    type t = Sundials.Realarray2.t

    (** {4 Basic access} *)

    (**
     [make m n] returns an [m] by [n] dense small matrix.

     @cvode <node9#ss:dense> newDenseMat
     *)
    val make : int -> int -> t

    (**
     [get a i j] returns the value at row [i] and column [j] in the m by
     n matrix [a], where 0 <= [i] < m and 0 <= [j] < n.
     *)
    val get : t -> int -> int -> float

    (**
     [set a i j v] stores the value [v] at row [i] and column [j] in the
     m by n matrix [a], where 0 <= [i] < m and 0 <= [j] < n.
     *)
    val set : t -> int -> int -> float -> unit

    (** {4 Calculations} *)

    (**
     [copy src dst] copies the contents of one matrix into another. Both
     must have the same size.

     @cvode <node9#ss:dense> denseCopy
     *)
    val copy  : t -> t -> unit

    (*
     [scale c a] multiplies each element of the matrix [a] by [c].

     @cvode <node9#ss:dense> denseScale
     *)
    val scale : float -> t -> unit

    (**
     [add_identity a] increments the (square) matrix [a] by the identity matrix.

     @cvode <node9#ss:dense> denseAddIdentity
     *)
    val add_identity : t -> unit

    (**
     [getrf a p] performs the LU factorization of an m by n matrix [a] with
     partial pivoting according to [p] (of length n).

     @cvode <node9#ss:dense> denseGETRF
     @raise ZeroDiagonalElement Zero found in matrix diagonal
     *)
    val getrf : t -> Sundials.lint_array -> unit

    (**
     [getrs a p b] finds the solution of [ax = b] using LU factorization.
     [a] must be an n by n matrix and [b] and [p] must be of length n.

     @cvode <node9#ss:dense> denseGETRS
     *)
    val getrs : t -> Sundials.lint_array -> Sundials.real_array -> unit

    (**
     [potrf a] performs the Cholesky factorization of a real symmetric positive
     n by n matrix.

     @cvode <node9#ss:dense> densePOTRF
     *)
    val potrf : t -> unit

    (**
     [potrs a b] finds the solution of [ax = b] using Cholesky
     factorization. [a] must be an n by n matrix and [b] must be of
     length n.

     @cvode <node9#ss:dense> densePOTRS
     *)
    val potrs : t -> Sundials.real_array -> unit

    (**
     [geqrf a beta work] performs the QR factorization of an
     m by n matrix, where m >= n. The [beta] vector must be of length n.
     The [work] vector must be of length m.

     @cvode <node9#ss:dense> denseGEQRF
     *)
    val geqrf : t -> Sundials.real_array -> Sundials.real_array -> unit

    (**
     [ormqr a beta v w work] computes the product [w = Qv], with Q calculated
     using {!geqrf}, where [a] is an m by n matrix with m >= n, [beta] is of
     length n, [v] is of length n, [w] is of length m, and [work] is of length
     m.

     @param a       matrix passed to {!geqrf}
     @param beta    vector apssed to {!geqrf}
     @param v       vector multiplier
     @param w       result vector
     @param work    temporary vector used in the calculation
     @cvode <node9#ss:dense> denseORMQR
     *)
    val ormqr :
      a:t -> beta:Sundials.real_array ->
      v:Sundials.real_array -> w:Sundials.real_array -> work:Sundials.real_array
      -> unit
  end

(** {3 Banded matrices}
    @cvode <node9#ss:band> The BAND Module *)

(** Operations for creating and manipulating banded matrices. *)
module BandMatrix :
  sig
    (**
    This type represents a [DlsMat] returned from a call to
    {!make}.

     @cvode <node9#s:dls>  Type DlsMat
     @cvode <node9#ss:band> NewBandMat 
     *)
    type t

    (** {4 Basic access} *)

    (**
     [make n mu ml smu] returns an [n] by [n] band matrix of upper
     bandwith [mu] and lower bandwidth [ml].
     - If [smu] = [mu], the result will {b not} be LU factored.
     - Otherwise pass [smu] = min([n]-1, [mu] + [ml]).

     @cvode <node9#ss:band> NewBandMat
     *)
    val make : int -> int -> int -> int -> t

    (**
     [mu, smu, ml = size a] returns the upper bandwidth ([mu]), upper storage
     bandwidth ([smu]), and lower bandwidth ([ml]) of the band matrix [a].
     *)
    val size  : t -> (int * int * int)

    (**
     Prints a band matrix to stdout.

     @cvode <node9#ss:band> PrintMat
     *)
    val print : t -> unit

    (**
     [get a i j] returns the value at row [i] and column [j] of the n by n
     matrix [a],
     where 0 <= [i], [j] <= n - 1 and [j] - mu(A) <= [i] <= [j] + ml(A).

     @cvode <node9#s:dls> BAND_ELEM
     *)
    val get : t -> int -> int -> float

    (**
      [set a i j v] stores the value [v] at row [i] and column [j] of the
      n by n matrix [a], where 0 <= [i], [j] <= n - 1 and
      [j] - mu(A) <= [i] <= [j] + ml(A).

      @cvode <node9#s:dls> BAND_ELEM
     *)
    val set : t -> int -> int -> float -> unit

    (** {4 Calculations} *)

    (**
     Fills the matrix with zeros.

     @cvode <node9#ss:band> SetToZero
     *)
    val set_to_zero    : t -> unit

    (**
     Increment a square matrix by the identity matrix.

     @cvode <node9#ss:band> AddIdentity
     *)
    val add_identity   : t -> unit

    (**
     [copy src dst copymu copyml] copies the submatrix with upper and lower
     bandwidths [copymu] and [copyml] of the n by n band matrix [src] into the n
     by n band matrix [dst].

     @cvode <node9#ss:band> BandCopy
     *)
    val copy : t -> t -> int -> int -> unit

    (**
     [scale c a] multiplies each element of [a] by [c].

     @cvode <node9#ss:band> BandScale
     *)
    val scale : float -> t -> unit

    (**
     [gbtrf a p] performs the LU factorization of [a] with partial pivoting
     according to [p].

     @cvode <node9#ss:band> BandGBTRF
     *)
    val gbtrf : t -> Sundials.lint_array -> unit

    (**
     [gbtrs a p b] finds the solution of [ax = b] using LU factorization.
     [a] must be a square matrix.

     @cvode <node9#ss:band> BandGBTRS
     *)
    val gbtrs : t -> Sundials.lint_array -> Sundials.real_array -> unit

    (** {4 Column access} *)

    (** Access banded matrix columns *)
    module Col :
      sig
        (**
         This type represents a bandmatrix ([DlsMat]) column.

         The underlying array cannot be exposed directly in OCaml as a
         {{:OCAML_DOC_ROOT(Bigarray)} Bigarray} because there is no simple way
         to attach a custom finalize function to such a big array.

         @cvode <node9#s:dls> BAND_COL
         *)
        type c

        (**
         [get_col a j] returns the diagonal element of the j-th column of the n
         by n band matrix [a], where 0 <= [j] <= n - 1.
         The resulting column may be indexed from -mu([a]) to ml([a]).

         @cvode <node9#s:dls> BAND_COL
         *)
        val get_col : t -> int -> c

        (**
         [get c i j] returns the ([i], [j])th entry of the band matrix from
         which the column [c] has already been selected;
         provided that [j] - mu(c) <= [i] <= [j] + ml(c).

         @cvode <node9#s:dls> BAND_COL_ELEM
         *)
        val get : c -> int -> int -> float

        (**
         [set c i j v] stores the value [v] at the ([i], [j])th entry of
         the band matrix from which the column [c] has already been selected;
         provided that [j] - mu(c) <= [i] <= [j] + ml(c).

         @cvode <node9#s:dls> BAND_COL_ELEM
         *)
        val set : c -> int -> int -> float -> unit
      end
  end

(** {3 Array-based banded matrices}
    @cvode <node9#ss:band> The BAND Module *)

(** Operations for creating and manipulating banded matrices based on arrays. *)
module ArrayBandMatrix :
  sig
    (**
     This type represents a two-dimensional matrix returned from a call to
     {!make}, {!Sundials.Realarray2.make}, or {!Sundials.Realarray2.wrap}.

     @cvode <node9#ss:band> newBandMat 
     *)
    type t = Sundials.Realarray2.t

    (** {4 Basic access} *)

    (**
     [make n smu ml] returns an [n] by [n] band matrix with lower
     half-bandwidth [ml] and storage upper bandwidth [smu].

     @cvode <node9#ss:band> newBandMat
     *)
    val make : int -> int -> int -> t

    (**
     [get a smu i j] returns the value at row [i] and column [j] in the 
     band matrix [a] with storage upper bandwidth [smu].
     *)
    val get : t -> int -> int -> int -> float

    (**
     [set a smu i j v] stores the value [v] at row [i] and column [j] in the
     band matrix [a] with storage upper bandwidth [smu].
     *)
    val set : t -> int -> int -> int -> float -> unit

    (** {4 Calculations} *)

    (**
     [copy src dst a_smu b_smu copymu copyml] copies the submatrix with
     upper and lower bandwidths [copymu] and [copyml] of the band
     matrix [src] into the band matrix [dst]. Both matrices must be
     square and of the same size.

     @cvode <node9#ss:band> bandCopy
     *)
    val copy : t -> t -> int -> int -> int -> int -> unit

    (**
     [scale c a mu ml smu] multiplies each element of the (square) band
     matrix [a], of upper bandwidth [mu], lower bandwidth [ml], and
     storage upper bandwidth [smu], by [c].

     @cvode <node9#ss:band> bandScale
     *)
    val scale : float -> t -> int -> int -> int -> unit

    (**
     [add_idenity a smu] increments the (square) band matrix [a] of storage
     upper bandwidth [smu] by the identity matrix.

     @cvode <node9#ss:band> bandAddIdentity
     *)
    val add_identity : t -> int -> unit

    (**
     [gbtrf a mu ml smu p] performs the LU factorization of the (square) band
     matrix [a], having bandwidths [mu] and [ml], and storage bandwidth [smu],
     with partial pivoting according to [p]. If [a] has size n by n then [p]
     must have length n.

     @cvode <node9#ss:band> bandGBTRF
     *)
    val gbtrf : t -> int -> int -> int -> Sundials.lint_array -> unit

    (**
     [gbtrs a smu ml p b] finds the solution of [ax = b] using LU factorization.
     [a] must be a (Square) band matrix having storage bandwidth [smu] and lower
     bandwidth [ml]. If [a] has size n by n then [p] and [b] both must have
     length n.

     @cvode <node9#ss:band> bandGBTRS
     *)
    val gbtrs : t -> int -> int -> Sundials.lint_array ->
      Sundials.real_array -> unit
  end

