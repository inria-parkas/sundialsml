(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(** Generic nvector types and operations.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

(** Represents an nvector of kind ['kind] with underlying data of type ['data].
    The type argument ['kind] is either {!Nvector_serial.kind},
    {!Nvector_parallel.kind}, or {!Nvector_custom.kind}.
    It is needed because some linear solvers make additional
    assumptions about the underlying vector representation.

    @cvode <node7#s:nvector> N_Vector *)
type ('data, 'kind) t

(** An alias for {!t}. *)
type ('data, 'kind) nvector = ('data, 'kind) t

(** [unwrap nv] returns the data underlying the nvector [nv]. *)
val unwrap : ('data, 'kind) t -> 'data

(** Raised when an nvector argument is incompatible with a session. 
    For example, when a solver session was initialized with an nvector
    having 10 elements, and a later call passes an nvector with only
    9 elements. The exact details depend on the nvector instantiation. *)
exception IncompatibleNvector

(** [check v1 v2] checks [v1] and [v2] for compatibility.

    @raise IncompatibleNvector The vectors are not compatible. *)
val check : ('data, 'kind) t -> ('data, 'kind) t -> unit

(** Generic vector operations.

    @cvode <node7> Description of the NVECTOR module. *)

(** Basic operations underlying an nvector. *)
module type NVECTOR_OPS =
  sig
    (** The vector type. *)
    type t

    (** Create a new, distinct vector from an existing one. *)
    val n_vclone        : t -> t

    (** [n_vlinearsum a x b y z] calculates [z = a*x + b*y]. *)
    val n_vlinearsum    : float -> t -> float -> t -> t -> unit

    (** [n_vconst c z] sets all of [z] to [c]. *)
    val n_vconst        : float -> t -> unit

    (** [n_vprod x y z] calculates [z = x * y] (pointwise). *)
    val n_vprod         : t -> t -> t -> unit

    (** [n_vdiv x y z] calculates [z = x / y] (pointwise). *)
    val n_vdiv          : t -> t -> t -> unit

    (** [n_vscale c x z] calculates [z = c *. x]. *)
    val n_vscale        : float -> t -> t -> unit

    (** [n_vabs x z] calculates [z = abs(x)]. *)
    val n_vabs          : t -> t -> unit

    (** [n_vinv x z] calculates [z = 1/x] (pointwise). *)
    val n_vinv          : t -> t -> unit

    (** [n_vaddconst x b z] calculates [z = x + b]. *)
    val n_vaddconst     : t -> float -> t -> unit

    (** [n_vdotprod x y] returns the dot product of [x] and [y]. *)
    val n_vdotprod      : t -> t -> float

    (** [n_vmaxnorm x] returns the maximum absolute value in x. *)
    val n_vmaxnorm      : t -> float

    (** [n_vwrmsnorm x w] returns the weighted root-mean-square norm of [x]
        with weight vector [w]. *)
    val n_vwrmsnorm     : t -> t -> float

    (** [n_vmin x] returns the smallest element in [x]. *)
    val n_vmin          : t -> float

    (** [n_vcompare c x z] calculates
        [z(i) = if abs x(i) >= c then 1 else 0]. *)
    val n_vcompare      : float -> t -> t -> unit

  (** [n_vinvtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
      zero values. This routine returns [true] if all components of [x] are
      nonzero (successful inversion) and [false] otherwise (not all elements
      inverted). *)
    val n_vinvtest      : t -> t -> bool

    (** [n_vwl2norm x w] returns the weighted ([w]) Euclidean l2 norm of [x]. *)
    val n_vwl2norm      : t -> t -> float

    (** [n_vl1norm x] returns the l1 norm of [x]. *)
    val n_vl1norm       : t -> float

    (** [n_vmaxnormmask x w id] returns the weighted root-mean-square norm
        of [x] using only elements where the corresponding [id] is non-zero. *)
    val n_vwrmsnormmask : t -> t -> t -> float

    (** [n_vconstrmask c x m] calculates [m(i) = Pi x(i)] returning the
        conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
        [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)
    val n_vconstrmask   : t -> t -> t -> bool

    (** [n_vminquotient num denom] returns the minimum of [num(i) / denom(i)].
        Zero [denom] elements are skipped. *)
    val n_vminquotient  : t -> t -> float
  end

(** Basic structure of a concrete nvector implementation module. *)
module type NVECTOR =
  sig
    (** Classifies the internal structure of an nvector. *)
    type kind

    (** The data wrapped within an nvector. *)
    type data

    (** An alias for the nvector type. *)
    type t = (data, kind) nvector

    (** Wrap data in an nvector. *)
    val wrap : data -> t

    (** Standard operations over nvectors. *)
    module Ops : NVECTOR_OPS with type t = t

    (** Standard operations over the underlying data. *)
    module DataOps : NVECTOR_OPS with type t = data
  end

