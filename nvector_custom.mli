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

(** An interface for creating custom nvectors in OCaml.
 
    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvode <node7#s:nvector> The NVECTOR Module *)

(** Represents an nvector manipulated by operations written in OCaml. Note that
    such operations entail the additional runtime cost of an OCaml callback. *)
type kind

(** The type scheme of custom nvectors. *)
type 'd t = ('d, kind) Nvector.t

(** The set of operations required to define an nvector. Some operations
    are optional; default values are either provided by the OCaml interface
    or the Sundials library.

    Custom nvectors and their operations must
    {ul
      {- Catch all exceptions since they and their backtraces cannot be
         propagated reliably. Backtraces can be produced manually by
         wrapping an operating in [try...with] and using
         {{:OCAML_DOC_ROOT(Printexc.html)} Printexc}.}
      {- Not create indirect reference loops back to the enclosing nvector. Such
         loops will not be properly garbage collected.}}

    @cvode <node7#s:nvector> _generic_N_Vector_Ops *)
type 'd nvector_ops = {
  n_vcheck           : 'd -> 'd ->  bool;
  (** Returns [true] if the vectors are compatible. See {!Nvector.check}. *)

  n_vclone           : 'd -> 'd;
  (** Creates a new, distinct vector from an existing one without
      necessarily copying the contents of the original vector. *)

  n_vspace           : ('d -> int * int) option;
  (** Returns storage requirements for one nvector [(lrw, liw)], where
      [lrw] is the number of realtype words and [liw] is the number of
      integer words . *)

  n_vlinearsum       : float -> 'd -> float -> 'd -> 'd -> unit;
  (** [n_vlinearsum a x b y z] calculates [z = ax + by]. *)

  n_vconst           : float -> 'd -> unit;
  (** [n_vconst c z] sets all of [z] to [c]. *)

  n_vprod            : 'd -> 'd -> 'd -> unit;
  (** [n_vprod x y z] calculates [z = x * y] (pointwise). *)

  n_vdiv             : 'd -> 'd -> 'd -> unit;
  (** [n_vdiv x y z] calculates [z = x / y] (pointwise). *)

  n_vscale           : float -> 'd -> 'd -> unit;
  (** [n_vscale c x z] calculates [z = c *. x]. *)

  n_vabs             : 'd -> 'd -> unit;
  (** [n_vabs x z] calculates [z = abs(x)]. *)

  n_vinv             : 'd -> 'd -> unit;
  (** [n_vinv x z] calculates [z = 1/x] (pointwise). *)

  n_vaddconst        : 'd -> float -> 'd -> unit;
  (** [n_vaddconst x b z] calculates [z = x + b]. *)

  n_vmaxnorm         : 'd -> float;
  (** [n_vmaxnorm x] returns the maximum absolute value in x. *)

  n_vwrmsnorm        : 'd -> 'd -> float;
  (** [n_vwrmsnorm x w] returns the weighted root-mean-square norm of [x]
      with weight vector [w]. *)

  n_vmin             : 'd -> float;
  (** [n_vmin x] returns the smallest element in [x]. *)

  n_vdotprod         : 'd -> 'd -> float;
  (** [n_vdotprod x y] returns the dot product of [x] and [y]. *)

  n_vcompare         : float -> 'd -> 'd -> unit;
  (** [n_vcompare c x z] calculates [z(i) = if abs x(i) >= c then 1 else 0]. *)

  n_vinvtest         : 'd -> 'd -> bool;
  (** [n_vinvtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
      zero values. This routine must return [true] if all components of [x] are
      nonzero (successful inversion) and [false] otherwise (not all elements
      inverted). *)

  n_vwl2norm         : ('d -> 'd -> float) option;
  (** [m = n_vwl2norm x w] returns the weighted Euclidean l_2 norm of [x]
      with weight vector [w], i.e.,
      {% $m = \sqrt{\sum_{i=0}^{n-1}(\mathtt{x}_i\mathtt{w}_i)}$ %}. *)

  n_vl1norm          : ('d -> float) option;
    (** [n_vl1norm x] returns the l1 norm of [x], i.e.,
         {% $\sum_{i=0}^{n-1}\lvert\mathtt{x}_i\rvert$ %}. *)

  n_vwrmsnormmask    : ('d -> 'd -> 'd -> float) option;
  (** [n_vmaxnormmask x w id] returns the weighted root-mean-square norm
      of [x] using only elements where the corresponding [id] is non-zero. *)

  n_vconstrmask      : ('d -> 'd -> 'd -> bool) option;
  (** [n_vconstrmask c x m] calculates [m(i) = Pi x(i)] returning the
      conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
      [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)

  n_vminquotient     : ('d -> 'd -> float) option;
  (** [n_vminquotient num denom] returns the minimum of [num(i) / denom(i)].
      Zero [denom] elements are skipped. If no such quotients are found,
      then {!Sundials.big_real} is returned. *)
}

(** [make_wrap ops] takes a set of operations on the data
    type ['d] and yields a function for lifting values of type ['d]
    into ['d] nvectors which can be passed to a solver. *)
val make_wrap  : 'd nvector_ops -> 'd -> 'd t

(** [add_tracing p ops] modifies a set of {!nvector_ops} so that
    a message, prefixed by [p], is printed each time an operation
    is called. This function is intended to help debug sets of
    vector operations. *)
val add_tracing     : string -> 'd nvector_ops -> 'd nvector_ops

(** Thrown for operations not provided to {!MakeOps} *)
exception OperationNotSupported

(** Turn a set of {!nvector_ops} into an nvector module. *)
module MakeOps : functor (A : sig
    type data
    val ops : data nvector_ops
  end) -> Nvector.NVECTOR with type data = A.data
                           and type kind = kind

