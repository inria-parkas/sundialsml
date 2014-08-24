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

(** Generic vector operations.

    @cvode <node7> Description of the NVECTOR module. *)
module type NVECTOR_OPS =
  sig
    (** Vector type *)
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

    (** [n_vdotprod x y] returns the dot produce of [x] and [y]. *)
    val n_vdotprod      : t -> t -> float

    (** [n_vmaxnorm x] returns the maximum absolute value in x. *)
    val n_vmaxnorm      : t -> float

    (** [n_vwrmsnorm x w] returns the weighted root-mean-square norm of [x]
        with weight vector [w]. *)
    val n_vwrmsnorm     : t -> t -> float

    (** [n_vmaxnormmask x w id] returns the weighted root-mean-square norm
        of [x] using only elements where the corresponding [id] is non-zero. *)
    val n_vwrmsnormmask : t -> t -> t -> float

    (** [n_vmin x] returns the smallest element in [x]. *)
    val n_vmin          : t -> float

    (** [n_vwl2norm x w] returns the weighted ([w]) Euclidean l2 norm of [x]. *)
    val n_vwl2norm      : t -> t -> float

    (** [n_vl1norm x] returns the l1 norm of [x]. *)
    val n_vl1norm       : t -> float

    (** [n_vcompare c x z] calculates
        [z(i) = if abs x(i) >= c then 1 else 0]. *)
    val n_vcompare      : float -> t -> t -> unit

    (** [n_vinvtest x z] calculates [z(i) = 1 / x(i)]. *)
    val n_vinvtest      : t -> t -> bool

    (** [n_vconstrmask c x m] calculates [m(i) = Pi x(i)] returning the
        conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
        [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)
    val n_vconstrmask   : t -> t -> t -> bool

    (** [n_vminquotient num denom] returns the minimum of [num(i) / denom(i)].
        Zero [denom] elements are skipped. *)
    val n_vminquotient  : t -> t -> float
  end

