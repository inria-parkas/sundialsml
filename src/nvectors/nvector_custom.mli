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

(** Represents an MPI communicator without introducing any unwanted
    dependencies on MPI. See {!Nvector_parallel.hide_communicator}. *)
type communicator

(** The set of operations required to define an nvector. Some operations
    are optional; default values are either provided by the OCaml interface
    or the Sundials library.

    Custom nvectors' payloads and operations must not hold any
    reference that loops back to the enclosing nvector. Such loops
    will not be properly garbage collected.

    Nvector operations are not allowed to throw exceptions, as they
    cannot be propagated reliably.  Uncaught exceptions will be
    discarded with a warning.

    Note that the fused and array custom nvector operations currently
    reallocate fresh arrays at each call. There is thus a tradeoff between
    the speed advantages of providing a single callback that handles many
    values at once but allocates more heap memory and multiple callbacks.

    @cvode <node7#s:nvector> _generic_N_Vector_Ops *)
type 'd nvector_ops = { (* {{{ *)
  check           : 'd -> 'd ->  bool;
  (** Returns [true] if the vectors are compatible. See {!Nvector.check}. *)

  clone           : 'd -> 'd;
  (** Creates a new, distinct vector from an existing one without
      necessarily copying the contents of the original vector. *)

  space           : ('d -> int * int) option;
  (** Returns storage requirements for one nvector [(lrw, liw)], where
      [lrw] is the number of realtype words and [liw] is the number of
      integer words . *)

  getlength       : 'd -> int;
  (** Returns the number of "active" entries. This value is cumulative
      across all processes in a parallel environment. *)

  print           : ('d -> Sundials.Logfile.t option -> unit) option;
  (** Print to the given logfile (stdout, by default). *)

  linearsum       : float -> 'd -> float -> 'd -> 'd -> unit;
  (** [linearsum a x b y z] calculates [z = ax + by]. *)

  const           : float -> 'd -> unit;
  (** [const c z] sets all of [z] to [c]. *)

  prod            : 'd -> 'd -> 'd -> unit;
  (** [prod x y z] calculates [z = x * y] (pointwise). *)

  div             : 'd -> 'd -> 'd -> unit;
  (** [div x y z] calculates [z = x / y] (pointwise). *)

  scale           : float -> 'd -> 'd -> unit;
  (** [scale c x z] calculates [z = c *. x]. *)

  abs             : 'd -> 'd -> unit;
  (** [abs x z] calculates [z = abs(x)]. *)

  inv             : 'd -> 'd -> unit;
  (** [inv x z] calculates [z = 1/x] (pointwise). *)

  addconst        : 'd -> float -> 'd -> unit;
  (** [addconst x b z] calculates [z = x + b]. *)

  maxnorm         : 'd -> float;
  (** [maxnorm x] returns the maximum absolute value in x. *)

  wrmsnorm        : 'd -> 'd -> float;
  (** [wrmsnorm x w] returns the weighted root-mean-square norm of [x]
      with weight vector [w]. *)

  min             : 'd -> float;
  (** [min x] returns the smallest element in [x]. *)

  dotprod         : 'd -> 'd -> float;
  (** [dotprod x y] returns the dot product of [x] and [y]. *)

  compare         : float -> 'd -> 'd -> unit;
  (** [compare c x z] calculates [z(i) = if abs x(i) >= c then 1 else 0]. *)

  invtest         : 'd -> 'd -> bool;
  (** [invtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
      zero values. This routine must return [true] if all components of [x] are
      nonzero (successful inversion) and [false] otherwise (not all elements
      inverted). *)

  wl2norm         : ('d -> 'd -> float) option;
  (** [m = wl2norm x w] returns the weighted Euclidean l_2 norm of [x]
      with weight vector [w], i.e.,
      {% $m = \sqrt{\sum_{i=0}^{n-1}(\mathtt{x}_i\mathtt{w}_i)}$ %}. *)

  l1norm          : ('d -> float) option;
    (** [l1norm x] returns the l1 norm of [x], i.e.,
         {% $\sum_{i=0}^{n-1}\lvert\mathtt{x}_i\rvert$ %}. *)

  wrmsnormmask    : ('d -> 'd -> 'd -> float) option;
  (** [maxnormmask x w id] returns the weighted root-mean-square norm
      of [x] using only elements where the corresponding [id] is non-zero. *)

  constrmask      : ('d -> 'd -> 'd -> bool) option;
  (** [constrmask c x m] calculates [m(i) = Pi x(i)] returning the
      conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
      [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)

  minquotient     : ('d -> 'd -> float) option;
  (** [minquotient num denom] returns the minimum of [num(i) / denom(i)].
      Zero [denom] elements are skipped. If no such quotients are found,
      then {{!Sundials_Config.big_real}Config.big_real} is returned. *)

  getcommunicator : ('d -> communicator) option;
  (** Returns the MPI communicator associated with an nvector. *)

  (* optional fused vector operations *)

  linearcombination :
    (Sundials.RealArray.t -> 'd array -> 'd -> unit) option;
  (** [linearcombination c x z] calculates
      [z(i) = c(0)*x(0)(i) + ... + c(nv-1)*x(nv-1)(i)] for the [nv] elements
      of [c] and [x], where [i] ranges over the nvector elements. *)

  scaleaddmulti :
    (Sundials.RealArray.t -> 'd -> 'd array -> 'd array -> unit) option;
  (** [scaleaddmulti c x y z] calculates
      [z(j)(i) = c(j)*x(i) + y(j)(i)], where [j] ranges over the array
      elements, and [i] ranges over the nvector elements. *)

  dotprodmulti :
    ('d -> 'd array -> Sundials.RealArray.t -> unit) option;
  (** [dotprodmulti x y d] calculates
      [d(j) = x(0)*y(j)(0) + ... + x(n-1)*y(j)(n-1)] for the [n] elements in
      the nvectors and where [j] ranges over the array elements. *)

  (* vector array operations *)

  linearsumvectorarray :
    (float -> 'd array -> float -> 'd array -> 'd array -> unit) option;
  (** [linearsumvectorarray a x b y z] calculates
      [z(j)(i) = a*x(j)(i) + b*y(j)(i)], where [j] ranges over the array
      elements and [i] ranges over the nvector elements. *)

  scalevectorarray :
    (Sundials.RealArray.t -> 'd array -> 'd array -> unit) option;
  (** [scalevectorarray c x z] calculates [z(j)(i) = c(j)*x(j)(i)],
      where [j] ranges over the array elements and [i] ranges over the nvector
      elements. *)

  constvectorarray :
    (float -> 'd array -> unit) option;
  (** [constvectorarray c x] sets [z(j)(i) = c],
      where [j] ranges over the array elements and [i] ranges over the
      nvector elements. *)

  wrmsnormvectorarray :
    ('d array -> 'd array -> Sundials.RealArray.t -> unit) option;
  (** [wrmsnormvectorarray x w m] calculates
      [m(j) = sqrt(((x(j)(0)*w(j)(0))^2 + ... + (x(j)(n-1)*w(j)(n-1))^2)/n)]
      for the [n] elements in the nvectors and where [j] ranges over the array
      elements. *)

  wrmsnormmaskvectorarray :
    ('d array -> 'd array -> 'd -> Sundials.RealArray.t -> unit) option;
  (** [wrmsnormvectorarray x w id m] calculates
      [m(j) = sqrt(((x(j)(0)*w(j)(0)*H(id(o)))^2 + ... + (x(j)(n-1)*w(j)(n-1)*H(id(n-1))^2)/n)]
      for the [n] elements in the nvectors, where [j] ranges over the array
      elements, and where [H(x) = if x > 0 then 1. else 0]. *)

  scaleaddmultivectorarray :
    (Sundials.RealArray.t -> 'd array -> 'd array array -> 'd array array ->
      unit) option;
  (** [scaleaddmultivectorarray a x yy zz] calculates
      [zz(j)(k)(i) = a(k)*x(k)(i) + yy(j)(k)(i)] where [j] and [k] range over
      the arrays and [i] ranges over the nvector elements. *)

  linearcombinationvectorarray :
    (Sundials.RealArray.t -> 'd array array -> 'd array -> unit) option;
  (** [linearcombinationvectorarray c xx z] calculates
      [z(k)(i) = c(0)*x(0)(k)(i) + ... + c(ns)*x(ns)(k)(i)] where [k] ranges
      over the array elements, [ns] is the number of arrays in [xx], and
      [i] ranges over the nvector elements. *)

  (* optional reduction operations *)

  dotprod_local      : ('d -> 'd -> float) option;
  (** Perform [dotprod] on task-local elements. *)

  maxnorm_local      : ('d -> float) option;
  (** Perform [maxnorm] on task-local elements. *)

  min_local          : ('d -> float) option;
  (** Returns the smallest task-local element. *)

  l1norm_local       : ('d -> float) option;
  (** Perform [l1norm] on task-local elements. *)

  invtest_local      : ('d -> 'd -> bool) option;
  (** Perform [invtest] on task-local elements. *)

  constrmask_local   : ('d -> 'd -> 'd -> bool) option;
  (** Perform [constrmask] on task-local elements. *)

  minquotient_local  : ('d -> 'd -> float) option;
  (** Perform [minquotient] on task-local elements. *)

  wsqrsum_local      : ('d -> 'd -> float) option;
  (** [wsqrsum x w] calculates the weighted squared sum of [x] with
      weight vector [w]. *)

  wsqrsummask_local  : ('d -> 'd -> 'd -> float) option;
  (** [wsqrsummask x w id] calculates the weighted squared sum of [x]
      with weight vector [w] for the elements where [id] is positive. *)
} (* }}} *)

(** Instantiation of custom nvectors.
    [make_wrap ops] takes set a set of operations on the data
    type ['d] and yields a function for lifting values of type ['d]
    into ['d] nvectors which can be passed to a solver. *)
val make_wrap  : 'd nvector_ops -> ?with_fused_ops:bool -> 'd -> 'd t

(** Add tracing to custom operations.
    [add_tracing p ops] modifies a set of {!nvector_ops} so that
    a message, prefixed by [p], is printed each time an operation
    is called. This function is intended to help debug sets of
    vector operations. *)
val add_tracing     : string -> 'd nvector_ops -> 'd nvector_ops

(** Selectively enable or disable fused and array operations.
    The [with_fused_ops] argument enables or disables all such operations.

    @since 4.0.0
    @raise OperationNotSupported The requested functionality has not been provided.
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val enable :
     ?with_fused_ops                       : bool
  -> ?with_linear_combination              : bool
  -> ?with_scale_add_multi                 : bool
  -> ?with_dot_prod_multi                  : bool
  -> ?with_linear_sum_vector_array         : bool
  -> ?with_scale_vector_array              : bool
  -> ?with_const_vector_array              : bool
  -> ?with_wrms_norm_vector_array          : bool
  -> ?with_wrms_norm_mask_vector_array     : bool
  -> ?with_scale_add_multi_vector_array    : bool
  -> ?with_linear_combination_vector_array : bool
  -> 'd t
  -> unit

(** Turn a set of {!nvector_ops} into an nvector module. *)
module MakeOps : functor (A : sig
    type data
    val ops : data nvector_ops
  end) -> Nvector.NVECTOR with type data = A.data
                           and type kind = kind

(** A generic nvector interface to custom nvectors.

    Create custom nvectors using the generic nvector interface where the
    payload is wrapped with a constructor from {!Nvector.gdata}. *)
module Any : sig (* {{{ *)

  (** Adapt a set of nvector operations so that they work with a payload of
      type {!Nvector.gdata}. It is better to manually implement the fused and
      array operations to avoid the creation of intermediate arrays.

      The [project] function should raise {!Nvector.BadGenericType} if applied
      to the wrong constructor. *)
  val convert_ops :
       inject:('d -> Nvector.gdata)
    -> project:(Nvector.gdata -> 'd)
    -> 'd nvector_ops
    -> Nvector.gdata nvector_ops

  (** Instantiation of custom nvectors.
      [make_wrap ops] takes set a set of operations on the data
      type {!Nvector.gdata} and yields a function for lifting values of
      type {!Nvector.gdata} into generic nvectors which can be passed to
      a solver.

      The optional arguments permit to enable fused and array operations for
      a given nvector (they are disabled by default).

      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 2.9.0 *)
  val make_wrap :
       Nvector.gdata nvector_ops
    -> inject:('d -> Nvector.gdata)
    -> ?with_fused_ops                       : bool
    -> ?with_linear_combination              : bool
    -> ?with_scale_add_multi                 : bool
    -> ?with_dot_prod_multi                  : bool
    -> ?with_linear_sum_vector_array         : bool
    -> ?with_scale_vector_array              : bool
    -> ?with_const_vector_array              : bool
    -> ?with_wrms_norm_vector_array          : bool
    -> ?with_wrms_norm_mask_vector_array     : bool
    -> ?with_scale_add_multi_vector_array    : bool
    -> ?with_linear_combination_vector_array : bool
    -> 'd
    -> Nvector.any

end (* }}} *)

