(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Direct Linear Solvers.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @cvode <node9#s:dls>  The DLS Modules *)

(** General purpose dense matrix operations on arrays.
    @cvode <node9#ss:dense> The DENSE Module *)
module ArrayDenseMatrix : sig (* {{{ *)

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

  (** Compute the matrix-vector product $y = Ax$.

      @nocvode <node9#ss:dense> denseMatvec
      @since 2.6.0 *)
  val matvec : t -> x:Sundials.RealArray.t -> y:Sundials.RealArray.t -> unit

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
      @raise Lsolver.ZeroDiagonalElement Zero found in matrix diagonal *)
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
module ArrayBandMatrix : sig (* {{{ *)

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

  (** Compute the matrix-vector product $y = Ax$.
  
      @nocvode <node9#ss:band> bandMatvec
      @since 2.6.0 *)
  val matvec : t -> smu -> mu -> ml
                -> x:Sundials.RealArray.t -> y:Sundials.RealArray.t -> unit

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

