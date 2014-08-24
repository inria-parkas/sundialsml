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

(***********************************************************************)
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Interface to Custom imperative and persistent NVectors.
 
    The Sundials solvers are written in a data-independent manner. They all
    operate on generic vectors through a set of operations defined by the
    particular nvector implementation. The OCaml interface replicates this
    design choice.

    The OCaml interface provides two ways to define custom nvectors: either
    through a set of (imperative) operations on an underlying mutable data type
    (default), or a set of (functional) operations on an underlying {!Immutable}
    data type. A set of {!Immutable} operations is converted into a set of
    mutable ones, which can be used with a solver, by
    {!Immutable.from_immutable} which introduces a reference in the underlying
    data type and modifies the set of operations to update it in-place.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @cvode <node7#s:nvector> The NVECTOR Module
 *)

(** Represents an nvector manipulated by operations written in OCaml. Note that
    such operations entail the additional runtime cost of a lookup and callback
    into OCaml. *)
type kind
type 'a t = ('a, kind) Sundials.nvector

(** {2 Custom imperative nvectors} *)

(**
   The set of operations required to define an nvector. Some operations
   are optional; default values are either provided by the OCaml interface
   or the Sundials library.

   @cvode <node7#s:nvector> _generic_N_Vector_Ops
 *)
type 'a nvector_ops = {

  n_vclone           : 'a -> 'a;
  (** Creates a new nvector of the same type as an existing vector.
      It need not copy the existing vector, but rather just allocate
      storage for the new vector. *)

  n_vdestroy         : ('a -> unit) option;
  (** Destroys the N_Vector v and frees memory allocated for its internal
      data. *)

  n_vspace           : ('a -> int * int) option;
  (** Returns storage requirements for one nvector. In the pair [(lrw,
      liw)], [lrw] is the number of realtype words required, and [liw] is
      the number of integer words required. *)

  n_vlinearsum       : float -> 'a -> float -> 'a -> 'a -> unit;
  (** [n_vlinearsum a x b y z] performs the operation
      [z] = [a]*[x] + [b]*[y], where [a] and [b] are scalars and
      [x] and [y] are vectors:
      [z]({i i}) = [a]*[x]({i i}) + [b]*[y]({i i}). *)

  n_vconst           : float -> 'a -> unit;
  (** [n_vconst c z] sets all components of [z] to [c]: [z]({i i}) = c. *)

  n_vprod            : 'a -> 'a -> 'a -> unit;
  (** [n_vprod x y z] sets [z] to be the component-wise product of [x]
      and [y]: [z]({i i}) = x({i i}) * y({i i}). *)

  n_vdiv             : 'a -> 'a -> 'a -> unit;
  (** [n_vdiv x y z] sets [z] to be the component-wise ratio of [x] and
      [y]: [z]({i i}) = [x]({i i}) / [y]({i i}).
      This function should only be called with a [y] whose components
      are all guaranteed to be nonzero. *)

  n_vscale           : float -> 'a -> 'a -> unit;
  (** [n_vscale c x z] scales [x] by [c] and returns the result in [z]:
      [z]({i i}) = [c] * [x]({i i}). *)

  n_vabs             : 'a -> 'a -> unit;
  (** [n_vabs x z] sets the components of [z] to the absolute values of
      the components of [x]: [z]({i i}) = |[x]({i i})|. *)

  n_vinv             : 'a -> 'a -> unit;
  (** [n_vinv x z] sets the components of [z] to be the inverses of the
     components of [x]: [z]({i i}) = 1.0 / [x]({i i}).
     This function should only be called with an [x] whose
     components are all guaranteed to be nonzero. *)

  n_vaddconst        : 'a -> float -> 'a -> unit;
  (** [n_vaddconst x b z] adds [b] to all components of [x] and
      returns the result in [z]: [z]({i i}) = [x]({i i}) + [b]. *)

  n_vmaxnorm         : 'a -> float;
  (** [m = n_vmaxnorm x] returns the maximum norm of [x]:
      [m] = |[x]({i i})|. *)

  n_vwrmsnorm        : 'a -> 'a -> float;
  (** [m = n_vwrmsnorm x w] returns the weighted root-mean-square norm of
      [x] with weight vector [w]:
        [m] = sqroot( sum( sqr([x]({i i}) * w({i i})) ) / n ). *)

  n_vmin             : 'a -> float;
  (** [m = n_vmin(x)] returns the smallest element of [x],
      [m] = min([x]({i i})). *)

  n_vdotprod         : ('a -> 'a -> float) option;
  (** [d = n_vdotprod x y] returns the ordinary dot product of [x] and [y],
      [d] = sum([x]({i i}) * [y]({i i})). *)

  n_vcompare         : (float -> 'a -> 'a -> unit) option;
  (** [n_vcompare c x z] compares the components of [x] to the scalar [c]
      and returns [z] such that [z]({i i}) = 1.0 if |[x]({i i})| >= [c]
      and [z]({i i}) = 0.0 otherwise. *)

  n_vinvtest         : ('a -> 'a -> bool) option;
  (** [t = n_vinvtest x z] sets the components of [z] to be the inverses of
      the components of [x], with prior testing for zero values:
        [z]({i i}) = 1.0 / [x]({i i}).
      This routine returns [true] if all components of [x] are nonzero
      (successful inversion) and [false] otherwise. *)

  n_vwl2norm         : ('a -> 'a -> float) option;
  (** [m = n_vwl2norm x w] returns the weighted Euclidean l_2 norm of [x]
      with weight vector [w]
      [m] = sqroot( sum( sqr([x]({i i}) * w({i i})) ) ). *)

  n_vl1norm          : ('a -> float) option;
  (** [m = n_vl1norm x] returns the l_1 norm of [x]:
      [m] = sum(|x({i i})|) *)

  n_vwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
  (** [m = n_vwrmsnormmask x w id] returns the weighted root mean square
      norm of [x] with weight vector [w] built using only the elements
      of [x] corresponding to nonzero elements of [id]:
        [m] = sqroot( sum( sqr([x]({i i}) * w({i i})
        * sign(id({i i}))) ) / n ).
   *)

  n_vconstrmask      : ('a -> 'a -> 'a -> bool) option;
  (** [t = n_vconstrmask c x m] performs the following constraint tests:
        - [x]({i i}) > 0 if [c]({i i})  = 2,
        - [x]({i i}) >= 0 if [c]({i i}) = 1,
        - [x]({i i}) <= 0 if [c]({i i}) = -1,
        - [x]({i i}) < 0 if [c]({i i})  = -2.

      There is no constraint on [x]({i i}) if [c]({i i}) = 0.
      This routine returns [false] if any element failed the constraint
      test, and [true] if all passed. It also sets a mask vector [m],
      with elements equal to 1.0 where the constraint test failed,
      and 0.0 where the test passed. This routine is used only for
      constraint checking. *)

  n_vminquotient     : ('a -> 'a -> float) option;
  (** [minq = n_vminquotient num denom] returns the minimum of the
      quotients obtained by term-wise dividing [num]({i i}) by
      [denom]({i i}). A zero element in [denom] will be skipped.
      If no such quotients are found, then {!Sundials.big_real}
      is returned. *)
}

(** [make_nvector ops] takes a set of operations on the data
    type ['a] and yields a function for lifting values of type ['a]
    into ['a] nvectors which can be passed to a solver. *)
val make  : 'a nvector_ops -> 'a -> 'a t

(** [add_tracing p ops] modifies a set of {!nvector_ops} so that
    a message, prefixed by [p], is printed each time an operation
    is called. This function is intended to help debug sets of
    vector operations. *)
val add_tracing     : string -> 'a nvector_ops -> 'a nvector_ops

(** {2 Custom persistent nvectors} *)

module Immutable :
sig
  type 'a mutable_nvector_ops = 'a nvector_ops

  type 'a nvector_ops = {
    n_vclone           : 'a -> 'a;
    (** Creates a new nvector of the same type as an existing vector.
        It need not copy the existing vector, but rather just allocate
        storage for the new vector. *)

    n_vdestroy         : ('a -> unit) option;
    (** Destroys the N_Vector v and frees memory allocated for its internal
        data. *)

    n_vspace           : ('a -> int * int) option;
    (** Returns storage requirements for one nvector. In the pair [(lrw,
        liw)], [lrw] is the number of realtype words required, and [liw] is
        the number of integer words required. *)

    n_vlinearsum       : float -> 'a -> float -> 'a -> 'a;
    (** [z = n_vlinearsum a x b y] performs the operation
        [z] = [a]*[x] + [b]*[y], where [a] and [b] are scalars and
        [x] and [y] are vectors:
        [z]({i i}) = [a]*[x]({i i}) + [b]*[y]({i i}). *)

    n_vconst           : float -> 'a;
    (** [z = n_vconst c] returns a [z] with all components equal to [c]:
        [z]({i i}) = c. *)

    n_vprod            : 'a -> 'a -> 'a;
    (** [z = n_vprod x y] returns the component-wise product of [x]
        and [y]: [z]({i i}) = x({i i}) * y({i i}). *)

    n_vdiv             : 'a -> 'a -> 'a;
    (** [z = n_vdiv x y] returns the component-wise ratio of [x] and
        [y]: [z]({i i}) = [x]({i i}) / [y]({i i}).
        This function should only be called with a [y] whose components
        are all guaranteed to be nonzero. *)

    n_vscale           : float -> 'a -> 'a;
    (** [z = n_vscale c x] returns the result of scaling [x] by [c]:
        [z]({i i}) = [c] * [x]({i i}). *)

    n_vabs             : 'a -> 'a;
    (** [z = n_vabs x] returns the absolute values of the components of [x]:
        [z]({i i}) = |[x]({i i})|. *)

    n_vinv             : 'a -> 'a;
    (** [z = n_vinv x] returns the inverses of the components of [x]:
        [z]({i i}) = 1.0 / [x]({i i}).
        This function should only be called with an [x] whose
        components are all guaranteed to be nonzero. *)

    n_vaddconst        : 'a -> float -> 'a;
    (** [z = n_vaddconst x b] returns the result of adding [b] to all
        components of [x]: [z]({i i}) = [x]({i i}) + [b]. *)

    n_vmaxnorm         : 'a -> float;
    (** [m = n_vmaxnorm x] returns the maximum norm of [x]:
        [m] = |[x]({i i})|. *)

    n_vwrmsnorm        : 'a -> 'a -> float;
    (** [m = n_vwrmsnorm x w] returns the weighted root-mean-square norm of
        [x] with weight vector [w]:
          [m] = sqroot( sum( sqr([x]({i i}) * w({i i})) ) / n ). *)

    n_vmin             : 'a -> float;
    (** [m = n_vmin(x)] returns the smallest element of [x],
        [m] = min([x]({i i})). *)

    n_vdotprod         : ('a -> 'a -> float) option;
    (** [d = n_vdotprod x y] returns the ordinary dot product of [x] and [y],
        [d] = sum([x]({i i}) * [y]({i i})). *)

    n_vcompare         : (float -> 'a -> 'a) option;
    (** [z = n_vcompare c x] compares the components of [x] to the scalar [c]
        and returns [z] such that [z]({i i}) = 1.0 if |[x]({i i})| >= [c]
        and [z]({i i}) = 0.0 otherwise. *)

    n_vinvtest         : ('a -> 'a -> bool) option;
    (** [t = n_vinvtest x z] sets the components of [z] to be the inverses of
        the components of [x], with prior testing for zero values:
          [z]({i i}) = 1.0 / [x]({i i}).
        This routine returns [true] if all components of [x] are nonzero
        (successful inversion) and [false] otherwise. *)

    n_vwl2norm         : ('a -> 'a -> float) option;
    (** [m = n_vwl2norm x w] returns the weighted Euclidean l_2 norm of [x]
        with weight vector [w]
        [m] = sqroot( sum( sqr([x]({i i}) * w({i i})) ) ). *)

    n_vl1norm          : ('a -> float) option;
    (** [m = n_vl1norm x] returns the l_1 norm of [x]:
        [m] = sum(|x({i i})|) *)

    n_vwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
    (** [m = n_vwrmsnormmask x w id] returns the weighted root mean square
        norm of [x] with weight vector [w] built using only the elements
        of [x] corresponding to nonzero elements of [id]:
          [m] = sqroot( sum( sqr([x]({i i}) * w({i i})
          * sign(id({i i}))) ) / n ).
     *)

    n_vconstrmask      : ('a -> 'a -> 'a -> bool) option;
    (** [t = n_vconstrmask c x m] performs the following constraint tests:
          - [x]({i i}) > 0 if [c]({i i})  = 2,
          - [x]({i i}) >= 0 if [c]({i i}) = 1,
          - [x]({i i}) <= 0 if [c]({i i}) = -1,
          - [x]({i i}) < 0 if [c]({i i})  = -2.

        There is no constraint on [x]({i i}) if [c]({i i}) = 0.
        This routine returns [false] if any element failed the constraint
        test, and [true] if all passed. It also sets a mask vector [m],
        with elements equal to 1.0 where the constraint test failed,
        and 0.0 where the test passed. This routine is used only for
        constraint checking. *)

    n_vminquotient     : ('a -> 'a -> float) option;
    (** [minq = n_vminquotient num denom] returns the minimum of the
        quotients obtained by term-wise dividing [num]({i i}) by
        [denom]({i i}). A zero element in [denom] will be skipped.
        If no such quotients are found, then {!Sundials.big_real}
        is returned. *)
  }

  (** Transforms a set of vector operations on an immutable data type into a
      set of vector operations on references to that data type. *)
  val from_immutable  : 'a nvector_ops -> 'a ref mutable_nvector_ops

  (** [make_nvector ops] takes a set of operations on the data
      type 'a and yields a function for lifting values of type ['a ref]
      into ['a ref] nvectors which can be passed to a solver. *)
  val make            : 'a nvector_ops -> 'a ref -> 'a ref t

  (** Extracts the data from an nvector. *)
  val unwrap          : 'a ref t -> 'a
end

