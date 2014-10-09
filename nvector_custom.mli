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

(** An interface for creating custom nvectors in OCaml.
 
    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @cvode <node7#s:nvector> The NVECTOR Module
 *)

(** Represents an nvector manipulated by operations written in OCaml. Note that
    such operations entail the additional runtime cost of an OCaml callback. *)
type kind

(** The type scheme of custom nvectors. *)
type 'a t = ('a, kind) Nvector.t

(**
   The set of operations required to define an nvector. Some operations
   are optional; default values are either provided by the OCaml interface
   or the Sundials library.

   @cvode <node7#s:nvector> _generic_N_Vector_Ops
 *)
type 'a nvector_ops = {

  n_vclone           : 'a -> 'a;
  (** Creates a new, distinct vector from an existing one without
      necessarily copying the contents of the original vector. *)

  n_vdestroy         : ('a -> unit) option;
  (** Destroys the N_Vector v and frees memory allocated for its internal
      data. *)

  n_vspace           : ('a -> int * int) option;
  (** Returns storage requirements for one nvector [(lrw, liw)], where
      [lrw] is the number of realtype words and [liw] is the number of
      integer words . *)

  n_vlinearsum       : float -> 'a -> float -> 'a -> 'a -> unit;
  (** [n_vlinearsum a x b y z] calculates [z = ax + by]. *)

  n_vconst           : float -> 'a -> unit;
  (** [n_vconst c z] sets all of [z] to [c]. *)

  n_vprod            : 'a -> 'a -> 'a -> unit;
  (** [n_vprod x y z] calculates [z = x * y] (pointwise). *)

  n_vdiv             : 'a -> 'a -> 'a -> unit;
  (** [n_vdiv x y z] calculates [z = x / y] (pointwise). *)

  n_vscale           : float -> 'a -> 'a -> unit;
  (** [n_vscale c x z] calculates [z = c *. x]. *)

  n_vabs             : 'a -> 'a -> unit;
  (** [n_vabs x z] calculates [z = abs(x)]. *)

  n_vinv             : 'a -> 'a -> unit;
  (** [n_vinv x z] calculates [z = 1/x] (pointwise). *)

  n_vaddconst        : 'a -> float -> 'a -> unit;
  (** [n_vaddconst x b z] calculates [z = x + b]. *)

  n_vmaxnorm         : 'a -> float;
  (** [n_vmaxnorm x] returns the maximum absolute value in x. *)

  n_vwrmsnorm        : 'a -> 'a -> float;
  (** [n_vwrmsnorm x w] returns the weighted root-mean-square norm of [x]
      with weight vector [w]. *)

  n_vmin             : 'a -> float;
  (** [n_vmin x] returns the smallest element in [x]. *)

  n_vdotprod         : 'a -> 'a -> float;
  (** [n_vdotprod x y] returns the dot product of [x] and [y]. *)

  n_vcompare         : float -> 'a -> 'a -> unit;
  (** [n_vcompare c x z] calculates [z(i) = if abs x(i) >= c then 1 else 0]. *)

  n_vinvtest         : 'a -> 'a -> bool;
  (** [n_vinvtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
      zero values. This routine must return [true] if all components of [x] are
      nonzero (successful inversion) and [false] otherwise (not all elements
      inverted). *)

  n_vwl2norm         : ('a -> 'a -> float) option;
  (** [m = n_vwl2norm x w] returns the weighted Euclidean l_2 norm of [x]
      with weight vector [w], i.e.,
      {% $m = \sqrt{\sum_{i=0}^{n-1}(\mathtt{x}_i\mathtt{w}_i)}$ %}. *)

  n_vl1norm          : ('a -> float) option;
    (** [n_vl1norm x] returns the l1 norm of [x], i.e.,
         {% $\sum_{i=0}^{n-1}\lvert\mathtt{x}_i\rvert$ %}. *)

  n_vwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
  (** [n_vmaxnormmask x w id] returns the weighted root-mean-square norm
      of [x] using only elements where the corresponding [id] is non-zero. *)

  n_vconstrmask      : ('a -> 'a -> 'a -> bool) option;
  (** [n_vconstrmask c x m] calculates [m(i) = Pi x(i)] returning the
      conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
      [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)

  n_vminquotient     : ('a -> 'a -> float) option;
  (** [n_vminquotient num denom] returns the minimum of [num(i) / denom(i)].
      Zero [denom] elements are skipped. If no such quotients are found,
      then {!Sundials.big_real} is returned. *)
}

(** [make_wrap ops] takes a set of operations on the data
    type ['a] and yields a function for lifting values of type ['a]
    into ['a] nvectors which can be passed to a solver. *)
val make_wrap  : 'a nvector_ops -> 'a -> 'a t

(** [add_tracing p ops] modifies a set of {!nvector_ops} so that
    a message, prefixed by [p], is printed each time an operation
    is called. This function is intended to help debug sets of
    vector operations. *)
val add_tracing     : string -> 'a nvector_ops -> 'a nvector_ops

(** Thrown for operations not provided to {!MakeOps} *)
exception OperationNotSupported

(** Turn a set of {!nvector_ops} into an nvector module. *)
module MakeOps : functor (A : sig
    type data
    val ops : data nvector_ops
  end) -> Nvector.NVECTOR with type data = A.data
                           and type kind = kind

