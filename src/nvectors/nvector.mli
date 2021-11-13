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

(** The type of any nvector that can be used as a serial nvector. *)
type 'k serial = (Sundials.RealArray.t, [>`Serial] as 'k) t

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

(** Clone an nvector. Cloning duplicates the payload and preserves the
    status of fused and array operations. *)
val clone : ('data, 'kind) t -> ('data, 'kind) t

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

(** {2:vecops Vector operations}

    @cvode <node7> Description of the NVECTOR module. *)

open Sundials

(** Basic operations underlying an nvector. *)
module type NVECTOR_OPS =
  sig (* {{{ *)
    (** The vector type. *)
    type t

    (** Create a new, distinct vector from an existing one. *)
    val clone        : t -> t

    (** [linearsum a x b y z] calculates [z = a*x + b*y]. *)
    val linearsum    : float -> t -> float -> t -> t -> unit

    (** [const c z] sets all of [z] to [c]. *)
    val const        : float -> t -> unit

    (** [prod x y z] calculates [z = x * y] (pointwise). *)
    val prod         : t -> t -> t -> unit

    (** [div x y z] calculates [z = x / y] (pointwise). *)
    val div          : t -> t -> t -> unit

    (** [scale c x z] calculates [z = c *. x]. *)
    val scale        : float -> t -> t -> unit

    (** [abs x z] calculates [z = abs(x)]. *)
    val abs          : t -> t -> unit

    (** [inv x z] calculates [z = 1/x] (pointwise). *)
    val inv          : t -> t -> unit

    (** [addconst x b z] calculates [z = x + b]. *)
    val addconst     : t -> float -> t -> unit

    (** [dotprod x y] returns the dot product of [x] and [y]. *)
    val dotprod      : t -> t -> float

    (** [maxnorm x] returns the maximum absolute value in x. *)
    val maxnorm      : t -> float

    (** [wrmsnorm x w] returns the weighted root-mean-square norm of [x]
        with weight vector [w]. *)
    val wrmsnorm     : t -> t -> float

    (** [min x] returns the smallest element in [x]. *)
    val min          : t -> float

    (** [compare c x z] calculates
        [z(i) = if abs x(i) >= c then 1 else 0]. *)
    val compare      : float -> t -> t -> unit

    (** [invtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
        zero values. This routine returns [true] if all components of [x] are
        nonzero (successful inversion) and [false] otherwise (not all elements
        inverted). *)
    val invtest      : t -> t -> bool

    (** [wl2norm x w] returns the weighted ([w]) Euclidean l2 norm of [x]. *)
    val wl2norm      : t -> t -> float

    (** [l1norm x] returns the l1 norm of [x]. *)
    val l1norm       : t -> float

    (** [maxnormmask x w id] returns the weighted root-mean-square norm
        of [x] using only elements where the corresponding [id] is non-zero. *)
    val wrmsnormmask : t -> t -> t -> float

    (** [constrmask c x m] calculates [m(i) = Pi x(i)] returning the
        conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
        [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)
    val constrmask   : t -> t -> t -> bool

    (** [minquotient num denom] returns the minimum of [num(i) / denom(i)].
        Zero [denom] elements are skipped. *)
    val minquotient  : t -> t -> float

    (** [lrw, liw = space c] returns the number of realtype words [lrw] and
        integer words [liw] required to store [c]. *)
    val space : t -> int * int

    (** Returns the number of "active" entries. This value is cumulative
        across all processes in a parallel environment.

       @since 5.0.0 *)
    val getlength : t -> int

    (** Prints to the given logfile (stdout, by default).

        @since 5.3.0 *)
    val print : ?logfile:Logfile.t -> t -> unit

    (* Fused and array operations *)

    (** [linearcombination c x z] calculates
        {% $z_i = \sum_{j=0}^{n_v-1} c_j (x_j)_i$ %}.
        The sum is over the {% $n_v$ %} elements of [c] and [x]. *)
    val linearcombination : Sundials.RealArray.t -> t array -> t -> unit

    (** [scaleaddmulti c x y z] scales [x] and adds it to the
        {% $n_v$ %} vectors in [y].
        That is,
        {% $\forall j=0,\ldots,n_v-1, (z_j)_i = c_j x_i + (y_j)_i}$ %}. *)
    val scaleaddmulti
      : Sundials.RealArray.t -> t -> t array -> t array -> unit

    (** [dotprodmulti x y d] calculates the dot product of [x] with
        the {% $n_v$ %} elements of [y].
        That is, {% $\forall j=0,\ldots,n_v-1,
                       d_j = \sum_{i=0}^{n-1} x_i (y_j)_i$ %}. *)
    val dotprodmulti      : t -> t array -> Sundials.RealArray.t -> unit

    (* {2: Vector array operations} *)

    (** [linearsumvectorarray a x b y z] computes the linear sum of the
        {% $n_v$ %} elements of [x] and [y].
        That is,
        {% $\forall j=0,\ldots,n_v-1, (z_j)_i = a (x_j)_i + b (y_j)_i$ %}. *)
    val linearsumvectorarray
      : float -> t array -> float -> t array -> t array -> unit

    (** [scalevectorarray c x z] scales each of the {% $n_v$ %} vectors
        of [x].
        That is, {% $\forall j=0,\ldots,n_v-1, (z_j)_i = c_j (x_j)_i$ %}. *)
    val scalevectorarray
      : Sundials.RealArray.t -> t array -> t array -> unit

    (** [constvectorarray c x] sets all elements of the {% $n_v$ %} nvectors
        in [x] to [c].
        That is, {% $\forall j=0,\ldots,n_v, (z_j)_i = c$ %}. *)
    val constvectorarray
      : float -> t array -> unit

    (** [wrmsnormvectorarray x w m] computes the weighted root mean
        square norm of the {% $n_v$ %} vectors in [x] and [w].
        That is,
        {% $\forall j=0,\ldots,n_v,
            m_j = \left( \frac{1}{n}
            \sum_{i=0}^{n-1} ((x_j)_i (w_j)_i)^2
            \right)^\frac{1}{2}$ %},
        where {% $n$ %} is the number of elements in each nvector. *)
    val wrmsnormvectorarray
      : t array -> t array -> Sundials.RealArray.t -> unit

    (** [wrmsnormmaskvectorarray x w id m] computes the weighted root mean
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
    val wrmsnormmaskvectorarray
      : t array -> t array -> t -> Sundials.RealArray.t -> unit

    (** [scaleaddmultivectorarray a x yy zz] scales and adds {% $n_v$ %}
        vectors in [x] across the {% $n_s$ %} vector arrays in [yy].
        That is, {% $\forall j=0,\ldots,n_s-1,
                     \forall k=0,\ldots,n_v-1,
          (\mathit{zz}_{j,k})_i = a_k (x_k)_i + (\mathit{yy}_{j,k})_i$ %}. *)
    val scaleaddmultivectorarray
      : Sundials.RealArray.t -> t array -> t array array -> t array array -> unit

    (** [linearcombinationvectorarray c xx z] computes the linear
        combinations of {% $n_s$ %} vector arrays containing {% $n_v$ %}
        vectors.
        That is, {% $\forall k=0,\ldots,n_v-1,
                      (z_k)_i = \sum_{j=0}^{n_s-1} c_j (x_{j,k})_i$ %}. *)
    val linearcombinationvectorarray
      : Sundials.RealArray.t -> t array array -> t array -> unit

    (* {2: Vector local reduction operations} *)

    (** Compute the task-local portions of certain operations. *)
    module Local : sig (* {{{ *)
      (** [dotprod x y] returns the dot product of [x] and [y]. *)
      val dotprod      : t -> t -> float

      (** [maxnorm x] returns the maximum absolute value in x. *)
      val maxnorm      : t -> float

      (** [min x] returns the smallest element in [x]. *)
      val min          : t -> float

      (** [l1norm x] returns the l1 norm of [x]. *)
      val l1norm       : t -> float

      (** [invtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
          zero values. This routine returns [true] if all components of [x] are
          nonzero (successful inversion) and [false] otherwise (not all elements
          inverted). *)
      val invtest      : t -> t -> bool

      (** [constrmask c x m] calculates [m(i) = Pi x(i)] returning the
          conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
          [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)
      val constrmask   : t -> t -> t -> bool

      (** [minquotient num denom] returns the minimum of [num(i) / denom(i)].
          Zero [denom] elements are skipped. *)
      val minquotient  : t -> t -> float

      (** [wsqrsum x w] calculates the weighted squared sum of [x] with
          weight vector [w].
          That is, {% $s = \sum_{i=0}^{n_\mathit{local} - 1}(x_i w_i)^2$ %}. *)
      val wsqrsum      : t -> t -> float

      (** [wsqrsummask x w id] calculates the weighted squared sum of [x]
          with weight vector [w] for the elements where [id] is positive.
          That is, {% $m =
          \sum_{i=0}^{n_\mathit{local} - 1}(x_i w_i H(\math{id}_i))^2$ %}
          where {% $H(\alpha) = \begin{cases}
                                  1 & \alpha > 0 \\
                                  0 & \alpha \le 0
                                \end{cases} $ %}. *)
      val wsqrsummask  : t -> t -> t -> float
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

(** Operations on any type of nvector. *)
module Ops : sig (* {{{ *)

    (** Create a new, distinct vector from an existing one. *)
    val clone        : ('d, 'k) t -> ('d, 'k) t

    (** [linearsum a x b y z] calculates [z = a*x + b*y]. *)
    val linearsum
      : float -> ('d, 'k) t -> float -> ('d, 'k) t -> ('d, 'k) t -> unit

    (** [const c z] sets all of [z] to [c]. *)
    val const        : float -> ('d, 'k) t -> unit

    (** [prod x y z] calculates [z = x * y] (pointwise). *)
    val prod         : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> unit

    (** [div x y z] calculates [z = x / y] (pointwise). *)
    val div          : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> unit

    (** [scale c x z] calculates [z = c *. x]. *)
    val scale        : float -> ('d, 'k) t -> ('d, 'k) t -> unit

    (** [abs x z] calculates [z = abs(x)]. *)
    val abs          : ('d, 'k) t -> ('d, 'k) t -> unit

    (** [inv x z] calculates [z = 1/x] (pointwise). *)
    val inv          : ('d, 'k) t -> ('d, 'k) t -> unit

    (** [addconst x b z] calculates [z = x + b]. *)
    val addconst     : ('d, 'k) t -> float -> ('d, 'k) t -> unit

    (** [dotprod x y] returns the dot product of [x] and [y]. *)
    val dotprod      : ('d, 'k) t -> ('d, 'k) t -> float

    (** [maxnorm x] returns the maximum absolute value in x. *)
    val maxnorm      : ('d, 'k) t -> float

    (** [wrmsnorm x w] returns the weighted root-mean-square norm of [x]
        with weight vector [w]. *)
    val wrmsnorm     : ('d, 'k) t -> ('d, 'k) t -> float

    (** [min x] returns the smallest element in [x]. *)
    val min          : ('d, 'k) t -> float

    (** [compare c x z] calculates
        [z(i) = if abs x(i) >= c then 1 else 0]. *)
    val compare      : float -> ('d, 'k) t -> ('d, 'k) t -> unit

    (** [invtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
        zero values. This routine returns [true] if all components of [x] are
        nonzero (successful inversion) and [false] otherwise (not all elements
        inverted). *)
    val invtest      : ('d, 'k) t -> ('d, 'k) t -> bool

    (** [wl2norm x w] returns the weighted ([w]) Euclidean l2 norm of [x]. *)
    val wl2norm      : ('d, 'k) t -> ('d, 'k) t -> float

    (** [l1norm x] returns the l1 norm of [x]. *)
    val l1norm       : ('d, 'k) t -> float

    (** [maxnormmask x w id] returns the weighted root-mean-square norm
        of [x] using only elements where the corresponding [id] is non-zero. *)
    val wrmsnormmask : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> float

    (** [constrmask c x m] calculates [m(i) = Pi x(i)] returning the
        conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
        [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)
    val constrmask   : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> bool

    (** [minquotient num denom] returns the minimum of [num(i) / denom(i)].
        Zero [denom] elements are skipped. *)
    val minquotient  : ('d, 'k) t -> ('d, 'k) t -> float

    (** [lrw, liw = space c] returns the number of realtype words [lrw] and
        integer words [liw] required to store [c]. *)
    val space : ('d, 'k) t -> int * int

    (** Returns the number of "active" entries. This value is cumulative
        across all processes in a parallel environment.

       @since 5.0.0 *)
    val getlength : ('d, 'k) t -> int

    (** Prints to the given logfile (stdout, by default).

        @since 5.3.0 *)
    val print : ?logfile:Logfile.t -> ('d, 'k) t -> unit

    (* Fused and array operations *)

    (** [linearcombination c x z] calculates
        {% $z_i = \sum_{j=0}^{n_v-1} c_j (x_j)_i$ %}.
        The sum is over the {% $n_v$ %} elements of [c] and [x]. *)
    val linearcombination
      : Sundials.RealArray.t -> ('d, 'k) t array -> ('d, 'k) t -> unit

    (** [scaleaddmulti c x y z] scales [x] and adds it to the
        {% $n_v$ %} vectors in [y].
        That is,
        {% $\forall j=0,\ldots,n_v-1, (z_j)_i = c_j x_i + (y_j)_i}$ %}. *)
    val scaleaddmulti
      :    Sundials.RealArray.t
        -> ('d, 'k) t
        -> ('d, 'k) t array
        -> ('d, 'k) t array
        -> unit

    (** [dotprodmulti x y d] calculates the dot product of [x] with
        the {% $n_v$ %} elements of [y].
        That is, {% $\forall j=0,\ldots,n_v-1,
                       d_j = \sum_{i=0}^{n-1} x_i (y_j)_i$ %}. *)
    val dotprodmulti
      : ('d, 'k) t -> ('d, 'k) t array -> Sundials.RealArray.t -> unit

    (* {2: Vector array operations} *)

    (** [linearsumvectorarray a x b y z] computes the linear sum of the
        {% $n_v$ %} elements of [x] and [y].
        That is,
        {% $\forall j=0,\ldots,n_v-1, (z_j)_i = a (x_j)_i + b (y_j)_i$ %}. *)
    val linearsumvectorarray
      :    float
        -> ('d, 'k) t array
        -> float
        -> ('d, 'k) t array
        -> ('d, 'k) t array
        -> unit

    (** [scalevectorarray c x z] scales each of the {% $n_v$ %} vectors
        of [x].
        That is, {% $\forall j=0,\ldots,n_v-1, (z_j)_i = c_j (x_j)_i$ %}. *)
    val scalevectorarray
      : Sundials.RealArray.t -> ('d, 'k) t array -> ('d, 'k) t array -> unit

    (** [constvectorarray c x] sets all elements of the {% $n_v$ %} nvectors
        in [x] to [c].
        That is, {% $\forall j=0,\ldots,n_v, (z_j)_i = c$ %}. *)
    val constvectorarray
      : float -> ('d, 'k) t array -> unit

    (** [wrmsnormvectorarray x w m] computes the weighted root mean
        square norm of the {% $n_v$ %} vectors in [x] and [w].
        That is,
        {% $\forall j=0,\ldots,n_v,
            m_j = \left( \frac{1}{n}
            \sum_{i=0}^{n-1} ((x_j)_i (w_j)_i)^2
            \right)^\frac{1}{2}$ %},
        where {% $n$ %} is the number of elements in each nvector. *)
    val wrmsnormvectorarray
      : ('d, 'k) t array -> ('d, 'k) t array -> Sundials.RealArray.t -> unit

    (** [wrmsnormmaskvectorarray x w id m] computes the weighted root mean
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
    val wrmsnormmaskvectorarray
      :    ('d, 'k) t array
        -> ('d, 'k) t array
        -> ('d, 'k) t
        -> Sundials.RealArray.t
        -> unit

    (** [scaleaddmultivectorarray a x yy zz] scales and adds {% $n_v$ %}
        vectors in [x] across the {% $n_s$ %} vector arrays in [yy].
        That is, {% $\forall j=0,\ldots,n_s-1,
                     \forall k=0,\ldots,n_v-1,
          (\mathit{zz}_{j,k})_i = a_k (x_k)_i + (\mathit{yy}_{j,k})_i$ %}. *)
    val scaleaddmultivectorarray
      :    Sundials.RealArray.t
        -> ('d, 'k) t array
        -> ('d, 'k) t array array
        -> ('d, 'k) t array array
        -> unit

    (** [linearcombinationvectorarray c xx z] computes the linear
        combinations of {% $n_s$ %} vector arrays containing {% $n_v$ %}
        vectors.
        That is, {% $\forall k=0,\ldots,n_v-1,
                      (z_k)_i = \sum_{j=0}^{n_s-1} c_j (x_{j,k})_i$ %}. *)
    val linearcombinationvectorarray
      :    Sundials.RealArray.t
        -> ('d, 'k) t array array
        -> ('d, 'k) t array
        -> unit

    (** {3:hasnvecop Availability of nvector array operations} *)

    (** Indicates whether an nvector supports {!linearcombination}. *)
    external has_linearcombination            : ('d, 'k) t -> bool
        = "sunml_nvec_has_linearcombination" [@@noalloc]

    (** Indicates whether an nvector supports {!scaleaddmulti}. *)
    external has_scaleaddmulti                : ('d, 'k) t -> bool
        = "sunml_nvec_has_scaleaddmulti" [@@noalloc]

    (** Indicates whether an nvector supports {!dotprodmulti}. *)
    external has_dotprodmulti                 : ('d, 'k) t -> bool
        = "sunml_nvec_has_dotprodmulti" [@@noalloc]

    (** Indicates whether an nvector supports {!linearsumvectorarray}. *)
    external has_linearsumvectorarray         : ('d, 'k) t -> bool
        = "sunml_nvec_has_linearsumvectorarray" [@@noalloc]

    (** Indicates whether an nvector supports {!scalevectorarray}. *)
    external has_scalevectorarray             : ('d, 'k) t -> bool
        = "sunml_nvec_has_scalevectorarray" [@@noalloc]

    (** Indicates whether an nvector supports {!constvectorarray}. *)
    external has_constvectorarray             : ('d, 'k) t -> bool
        = "sunml_nvec_has_constvectorarray" [@@noalloc]

    (** Indicates whether an nvector supports {!wrmsnormvectorarray}. *)
    external has_wrmsnormvectorarray          : ('d, 'k) t -> bool
        = "sunml_nvec_has_wrmsnormvectorarray" [@@noalloc]

    (** Indicates whether an nvector supports {!wrmsnormmaskvectorarray}. *)
    external has_wrmsnormmaskvectorarray      : ('d, 'k) t -> bool
        = "sunml_nvec_has_wrmsnormmaskvectorarray" [@@noalloc]

    (** Indicates whether an nvector supports {!scaleaddmultivectorarray}. *)
    external has_scaleaddmultivectorarray     : ('d, 'k) t -> bool
        = "sunml_nvec_has_scaleaddmultivectorarray" [@@noalloc]

    (** Indicates whether an nvector supports {!linearcombinationvectorarray}. *)
    external has_linearcombinationvectorarray : ('d, 'k) t -> bool
        = "sunml_nvec_has_linearcombinationvectorarray" [@@noalloc]

    (* {2: Vector local reduction operations} *)

    (** Compute the task-local portions of certain operations. *)
    module Local : sig (* {{{ *)
      (** [dotprod x y] returns the dot product of [x] and [y]. *)
      val dotprod      : ('d, 'k) t -> ('d, 'k) t -> float

      (** [maxnorm x] returns the maximum absolute value in x. *)
      val maxnorm      : ('d, 'k) t -> float

      (** [min x] returns the smallest element in [x]. *)
      val min          : ('d, 'k) t -> float

      (** [l1norm x] returns the l1 norm of [x]. *)
      val l1norm       : ('d, 'k) t -> float

      (** [invtest x z] calculates [z(i) = 1 / x(i)] with prior testing for
          zero values. This routine returns [true] if all components of [x] are
          nonzero (successful inversion) and [false] otherwise (not all elements
          inverted). *)
      val invtest      : ('d, 'k) t -> ('d, 'k) t -> bool

      (** [constrmask c x m] calculates [m(i) = Pi x(i)] returning the
          conjunction. The value of [Pi] depends on [c(i)]: [2: x(i) > 0],
          [1: x(i) >= 0], [0: true], [-1: x(i) <= 0], and [-2: x(i) < 0]. *)
      val constrmask   : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> bool

      (** [minquotient num denom] returns the minimum of [num(i) / denom(i)].
          Zero [denom] elements are skipped. *)
      val minquotient  : ('d, 'k) t -> ('d, 'k) t -> float

      (** [wsqrsum x w] calculates the weighted squared sum of [x] with
          weight vector [w].
          That is, {% $s = \sum_{i=0}^{n_\mathit{local} - 1}(x_i w_i)^2$ %}. *)
      val wsqrsum      : ('d, 'k) t -> ('d, 'k) t -> float

      (** [wsqrsummask x w id] calculates the weighted squared sum of [x]
          with weight vector [w] for the elements where [id] is positive.
          That is, {% $m =
          \sum_{i=0}^{n_\mathit{local} - 1}(x_i w_i H(\math{id}_i))^2$ %}
          where {% $H(\alpha) = \begin{cases}
                                  1 & \alpha > 0 \\
                                  0 & \alpha \le 0
                                \end{cases} $ %}. *)
      val wsqrsummask  : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> float

      (** {3:hasnveclop Availability of nvector local operations} *)

      (** Indicates whether an nvector supports a local {!dotprod}. *)
      external has_dotprod      : ('d, 'k) t -> bool
        = "sunml_nvec_has_dotprodlocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!maxnorm}. *)
      external has_maxnorm      : ('d, 'k) t -> bool
        = "sunml_nvec_has_maxnormlocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!min}. *)
      external has_min          : ('d, 'k) t -> bool
        = "sunml_nvec_has_minlocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!l1norm}. *)
      external has_l1norm       : ('d, 'k) t -> bool
        = "sunml_nvec_has_l1normlocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!invtest}. *)
      external has_invtest      : ('d, 'k) t -> bool
        = "sunml_nvec_has_invtestlocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!constrmask}. *)
      external has_constrmask   : ('d, 'k) t -> bool
        = "sunml_nvec_has_constrmasklocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!minquotient}. *)
      external has_minquotient  : ('d, 'k) t -> bool
        = "sunml_nvec_has_minquotientlocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!wsqrsum}. *)
      external has_wsqrsum      : ('d, 'k) t -> bool
        = "sunml_nvec_has_wsqrsumlocal" [@@noalloc]

      (** Indicates whether an nvector supports a local {!wsqrsummask}. *)
      external has_wsqrsummask  : ('d, 'k) t -> bool
        = "sunml_nvec_has_wsqrsummasklocal" [@@noalloc]

    end (* }}} *)

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

(** Generic wrapper for {{!Sundials_RealArray.t}RealArray}. *)
type gdata += RA of Sundials.RealArray.t

(** Represents an nvector whose data must be accessed through a constructor
    in {!gdata}. *)
type gkind

(** The type of a generic nvector. *)
type any = (gdata, gkind) t

(** A {!gdata} value did not have the expected wrapper. *)
exception BadGenericType

(** The requested operation is not provided by the given nvector. *)
exception OperationNotProvided

(** Lifts a set of operations on a type [t] to a set of operations on a
    generic nvector payload. *)
module MakeDataOps :
  functor (Ops : sig
    include NVECTOR_OPS

    val unwrap : gdata -> t
    val wrap   : t -> gdata
  end) -> NVECTOR_OPS with type t = gdata

