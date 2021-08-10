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

(** Generic nvector types and operations.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

(** Represents an nvector of kind ['kind] with underlying data of type ['data].
    The type argument ['kind] is either {!Nvector_serial.kind},
    {!Nvector_parallel.kind}, {!Nvector_custom.kind}, {!Nvector_openmp.kind},
    or {!Nvector_pthreads.kind}.
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

(** Vector type identifiers. *)
type nvector_id =
    Serial
  | Parallel
  | OpenMP
  | Pthreads
  | ParHyp
  | PETSc
  | CUDA
  | RAJA
  | OpenMPdev
  | Trilinos
  | ManyVector
  | MpiManyVector
  | MpiPlusX
  | Custom

(** Returns the vector type identifier.

    @since 2.9.0 *)
val get_id : ('data, 'kind) t -> nvector_id

(** {2:hasops Test vector operations}

    These functions test whether an nvector is equipped with a given
    operation. *)
(* {{{ *)

external has_n_vlinearcombination            : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vlinearcombination" [@@noalloc]
external has_n_vscaleaddmulti                : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vscaleaddmulti" [@@noalloc]
external has_n_vdotprodmulti                 : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vdotprodmulti" [@@noalloc]
external has_n_vlinearsumvectorarray         : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vlinearsumvectorarray" [@@noalloc]
external has_n_vscalevectorarray             : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vscalevectorarray" [@@noalloc]
external has_n_vconstvectorarray             : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vconstvectorarray" [@@noalloc]
external has_n_vwrmsnormvectorarray          : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwrmsnormvectorarray" [@@noalloc]
external has_n_vwrmsnormmaskvectorarray      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwrmsnormmaskvectorarray" [@@noalloc]
external has_n_vscaleaddmultivectorarray     : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vscaleaddmultivectorarray" [@@noalloc]
external has_n_vlinearcombinationvectorarray : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vlinearcombinationvectorarray" [@@noalloc]

module Local : sig
  external has_n_vdotprod      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vdotprodlocal" [@@noalloc]
  external has_n_vmaxnorm      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vmaxnormlocal" [@@noalloc]
  external has_n_vmin          : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vminlocal" [@@noalloc]
  external has_n_vl1norm       : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vl1normlocal" [@@noalloc]
  external has_n_vinvtest      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vinvtestlocal" [@@noalloc]
  external has_n_vconstrmask   : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vconstrmasklocal" [@@noalloc]
  external has_n_vminquotient  : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vminquotientlocal" [@@noalloc]
  external has_n_vwsqrsum      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwsqrsumlocal" [@@noalloc]
  external has_n_vwsqrsummask  : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwsqrsummasklocal" [@@noalloc]
end

(* }}} *)

(** {2:vecops Vector operations}

    @cvode <node7> Description of the NVECTOR module. *)

open Sundials

(** Basic operations underlying an nvector. *)
module type NVECTOR_OPS =
  sig (* {{{ *)
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

    (** [lrw, liw = n_vspace c] returns the number of realtype words [lrw] and
        integer words [liw] required to store [c]. *)
    val n_vspace : t -> int * int

    (** Returns the number of "active" entries. This value is cumulative
        across all processes in a parallel environment.

       @since 5.0.0 *)
    val n_vgetlength : t -> int

    (* Fused and array operations *)

    (** [n_vlinearcombination c x z] calculates
        {% $z_i = \sum_{j=0}^{n_v-1} c_j (x_j)_i$ %}.
        The sum is over the {% $n_v$ %} elements of [c] and [x]. *)
    val n_vlinearcombination : Sundials.RealArray.t -> t array -> t -> unit

    (** [n_vscaleaddmulti c x y z] scales [x] and adds it to the
        {% $n_v$ %} vectors in [y].
        That is,
        {% $\forall j=0,\ldots,n_v-1, (z_j)_i = c_j x_i + (y_j)_i}$ %}. *)
    val n_vscaleaddmulti
      : Sundials.RealArray.t -> t -> t array -> t array -> unit

    (** [n_vdotprodmulti x y d] calculates the dot product of [x] with
        the {% $n_v$ %} elements of [y].
        That is, {% $\forall j=0,\ldots,n_v-1,
                       d_j = \sum_{i=0}^{n-1} x_i (y_j)_i$ %}. *)
    val n_vdotprodmulti      : t -> t array -> Sundials.RealArray.t -> unit

    (* {2: Vector array operations} *)

    (** [n_vlinearsumvectorarray a x b y z] computes the linear sum of the
        {% $n_v$ %} elements of [x] and [y].
        That is,
        {% $\forall j=0,\ldots,n_v-1, (z_j)_i = a (x_j)_i + b (y_j)_i$ %}. *)
    val n_vlinearsumvectorarray
      : float -> t array -> float -> t array -> t array -> unit

    (** [n_vscalevectorarray c x z] scales each of the {% $n_v$ %} vectors
        of [x].
        That is, {% $\forall j=0,\ldots,n_v-1, (z_j)_i = c_j (x_j)_i$ %}. *)
    val n_vscalevectorarray
      : Sundials.RealArray.t -> t array -> t array -> unit

    (** [n_vconstvectorarray c x] sets all elements of the {% $n_v$ %} nvectors
        in [x] to [c].
        That is, {% $\forall j=0,\ldots,n_v, (z_j)_i = c$ %}. *)
    val n_vconstvectorarray
      : float -> t array -> unit

    (** [n_vwrmsnormvectorarray x w m] computes the weighted root mean
        square norm of the {% $n_v$ %} vectors in [x] and [w].
        That is,
        {% $\forall j=0,\ldots,n_v,
            m_j = \left( \frac{1}{n}
            \sum_{i=0}^{n-1} ((x_j)_i (w_j)_i)^2
            \right)^\frac{1}{2}$ %},
        where {% $n$ %} is the number of elements in each nvector. *)
    val n_vwrmsnormvectorarray
      : t array -> t array -> Sundials.RealArray.t -> unit

    (** [n_vwrmsnormmaskvectorarray x w id m] computes the weighted root mean
        square norm of the {% $n_v$ %} vectors in [x] and [w].
        That is,
        {% $\forall j=0,\ldots,n_v,
            m_j = \left( \frac{1}{n}
            \sum_{i=0}^{n-1} ((x_j)_i (w_j)_i H(\mathit{id}_i))^2
            \right)^\frac{1}{2}$ %},
        where {% $H(x) = \begin{cases} 1 & \text{if } x > 0 \\
                                       0 & \text{otherwise}
                         \end{cases}$ %} and
        {% $n$ %} is the number of elements in each nvector. *)
    val n_vwrmsnormmaskvectorarray
      : t array -> t array -> t -> Sundials.RealArray.t -> unit

    (** [n_vscaleaddmultivectorarray a x yy zz] scales and adds {% $n_v$ %}
        vectors in [x] across the {% $n_s$ %} vector arrays in [yy].
        That is, {% $\forall j=0,\ldots,n_s-1,
                     \forall k=0,\ldots,n_v-1,
          (\mathit{zz}_{j,k})_i = a_k (x_k)_i + (\mathit{yy}_{j,k})_i$ %}. *)
    val n_vscaleaddmultivectorarray
      : Sundials.RealArray.t -> t array -> t array array -> t array array -> unit

    (** [n_vlinearcombinationvectorarray c xx z] computes the linear
        combinations of {% $n_s$ %} vector arrays containing {% $n_v$ %}
        vectors.
        That is, {% $\forall k=0,\ldots,n_v-1,
                      (z_k)_i = \sum_{j=0}^{n_s-1} c_j (x_{j,k})_i$ %}. *)
    val n_vlinearcombinationvectorarray
      : Sundials.RealArray.t -> t array array -> t array -> unit

    (* {2: Vector local reduction operations} *)

    (** Compute the task-local portions of certain operations. *)
    module Local : sig (* {{{ *)
      (** [n_vdotprod x y] returns the dot product of [x] and [y]. *)
      val n_vdotprod      : t -> t -> float

      (** [n_vmaxnorm x] returns the maximum absolute value in x. *)
      val n_vmaxnorm      : t -> float

      (** [n_vmin x] returns the smallest element in [x]. *)
      val n_vmin          : t -> float

      (** [n_vl1norm x] returns the l1 norm of [x]. *)
      val n_vl1norm       : t -> float

      (** [n_vinvtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
          zero values. This routine returns [true] if all components of [x] are
          nonzero (successful inversion) and [false] otherwise (not all elements
          inverted). *)
      val n_vinvtest      : t -> t -> bool

      (** [n_vconstrmask c x m] calculates [m(i) = Pi x(i)] returning the
          conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
          [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)
      val n_vconstrmask   : t -> t -> t -> bool

      (** [n_vminquotient num denom] returns the minimum of [num(i) / denom(i)].
          Zero [denom] elements are skipped. *)
      val n_vminquotient  : t -> t -> float

      (** [n_vwsqrsum x w] calculates the weighted squared sum of [x] with
          weight vector [w].
          That is, {% $s = \sum_{i=0}^{n_\mathit{local} - 1}(x_i w_i)^2$ %}. *)
      val n_vwsqrsum      : t -> t -> float

      (** [n_vwsqrsummask x w id] calculates the weighted squared sum of [x]
          with weight vector [w] for the elements where [id] is positive.
          That is, {% $m =
          \sum_{i=0}^{n_\mathit{local} - 1}(x_i w_i H(\math{id}_i))^2$ %}
          where {% $H(\alpha) = \begin{cases}
                                  1 & \alpha > 0 \\
                                  0 & \alpha \le 0
                                \end{cases} $ %}. *)
      val n_vwsqrsummask  : t -> t -> t -> float
    end (* }}} *)

  end (* }}} *)

(** Basic structure of a concrete nvector implementation module. *)
module type NVECTOR =
  sig (* {{{ *)

    (** Classifies the internal structure of an nvector. *)
    type kind

    (** The data wrapped within an nvector. *)
    type data

    (** An alias for the nvector type. *)
    type t = (data, kind) nvector

    (** Wrap data in an nvector. *)
    val wrap : ?with_fused_ops:bool -> data -> t

    (** Selectively enable or disable fused and array operations. *)
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
      -> t
      -> unit

    (** Standard operations over nvectors. *)
    module Ops : NVECTOR_OPS with type t = t

    (** Standard operations over the underlying data. *)
    module DataOps : NVECTOR_OPS with type t = data

  end (* }}} *)

(** {2:genvec Generic nvectors}

    The two arguments of the standard nvector type encode, respectively, the
    type of the data payload that is manipulated from OCaml and the type of
    the underlying implementation. This encoding has two advantages. First,
    callback functions have direct access to the data payload. Second,
    incorrect combinations of nvectors can be rejected at compile time
    (although dynamic checks are still required on payload dimensions). There
    are two main disadvantages to this encoding. First, type signatures and
    error messages become more complicated. Second, it is not possible to form
    lists or arrays of heterogeneous nvectors, despite the object-oriented
    style of the underlying implementation.

    Generic nvectors are an alternative interface that inverses the advantages
    and disadvantages of standard nvectors. The data payload is now accessed
    through a constructor that otherwise hides the underlying type. Since all
    such nvectors have the same values for the two type arguments, they can be
    stored together in lists and arrays regardless of their underlying
    payloads and implementations. The price to pay is additional run-time
    checks and the fact that certain errors can no longer be detected
    statically. *)

(** Represents generic nvector data. This type is extensible so that other
    nvector modules can support the generic form. *)
type gdata = ..

(** The wrapper for {{!Sundials.RealArray.t}RealArray}s. *)
type gdata += RA of Sundials.RealArray.t

(** Represents an nvector whose data must be accessed through a constructor
    in {!gdata}. *)
type gkind

(** The type of a generic nvector. *)
type any = (gdata, gkind) t

module Any : sig (* {{{ *)
  type t = any

  external has_n_vlinearcombination            : t -> bool
      = "sunml_nvec_has_n_vlinearcombination" [@@noalloc]
  external has_n_vscaleaddmulti                : t -> bool
      = "sunml_nvec_has_n_vscaleaddmulti" [@@noalloc]
  external has_n_vdotprodmulti                 : t -> bool
      = "sunml_nvec_has_n_vdotprodmulti" [@@noalloc]
  external has_n_vlinearsumvectorarray         : t -> bool
      = "sunml_nvec_has_n_vlinearsumvectorarray" [@@noalloc]
  external has_n_vscalevectorarray             : t -> bool
      = "sunml_nvec_has_n_vscalevectorarray" [@@noalloc]
  external has_n_vconstvectorarray             : t -> bool
      = "sunml_nvec_has_n_vconstvectorarray" [@@noalloc]
  external has_n_vwrmsnormvectorarray          : t -> bool
      = "sunml_nvec_has_n_vwrmsnormvectorarray" [@@noalloc]
  external has_n_vwrmsnormmaskvectorarray      : t -> bool
      = "sunml_nvec_has_n_vwrmsnormmaskvectorarray" [@@noalloc]
  external has_n_vscaleaddmultivectorarray     : t -> bool
      = "sunml_nvec_has_n_vscaleaddmultivectorarray" [@@noalloc]
  external has_n_vlinearcombinationvectorarray : t -> bool
      = "sunml_nvec_has_n_vlinearcombinationvectorarray" [@@noalloc]

  module Local : sig
    external has_n_vdotprod      : t -> bool
      = "sunml_nvec_has_n_vdotprodlocal" [@@noalloc]
    external has_n_vmaxnorm      : t -> bool
      = "sunml_nvec_has_n_vmaxnormlocal" [@@noalloc]
    external has_n_vmin          : t -> bool
      = "sunml_nvec_has_n_vminlocal" [@@noalloc]
    external has_n_vl1norm       : t -> bool
      = "sunml_nvec_has_n_vl1normlocal" [@@noalloc]
    external has_n_vinvtest      : t -> bool
      = "sunml_nvec_has_n_vinvtestlocal" [@@noalloc]
    external has_n_vconstrmask   : t -> bool
      = "sunml_nvec_has_n_vconstrmasklocal" [@@noalloc]
    external has_n_vminquotient  : t -> bool
      = "sunml_nvec_has_n_vminquotientlocal" [@@noalloc]
    external has_n_vwsqrsum      : t -> bool
      = "sunml_nvec_has_n_vwsqrsumlocal" [@@noalloc]
    external has_n_vwsqrsummask  : t -> bool
      = "sunml_nvec_has_n_vwsqrsummasklocal" [@@noalloc]
  end
end (* }}} *)

(** A {!gdata} value did not have the expected wrapper. *)
exception BadGenericType

(** Operations on generic nvectors. Note that generic nvectors cannot be
    cloned from OCaml.

    @raise Invalid_argument If {!Ops.nv_clone} is called. *)
module Ops : NVECTOR_OPS with type t = any

