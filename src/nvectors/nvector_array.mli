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

(** A custom nvector based on float arrays.

  @version VERSION()
  @author Timothy Bourke (Inria/ENS)
  @author Jun Inoue (Inria/ENS)
  @author Marc Pouzet (UPMC/ENS/Inria)
 *)

(** An abstract set of functions for working manipulating nvectors
    where the underlying data structure is an array of [float]s.  *)
module type ARRAY_NVECTOR =
  sig
    (** Type of the underlying array. *)
    type data

    (** Array nvectors are custom nvectors. *)
    type kind = Nvector_custom.kind

    (** An alias for the nvector type. *)
    type t = data Nvector_custom.t

    (** The set of nvector operations on an array. *)
    val array_nvec_ops  : data Nvector_custom.nvector_ops

    (** [make n x] creates an nvector containing an array
        of [n] elements, each of which is equal to [x]. *)
    val make            : int -> float -> t

    (** Lifts an array to an nvector. *)
    val wrap            : ?with_fused_ops:bool -> data -> t

    (** Returns the array underlying an nvector. *)
    val unwrap          : t -> data

    (** Selectively enable or disable fused and array operations.
        The [with_fused_ops] argument enables or disables all such operations.

        @since 4.0.0
        @cvode <node5> N_VEnableFusedOps_Serial
        @cvode <node5> N_VEnableLinearCombination_Serial
        @cvode <node5> N_VEnableScaleAddMulti_Serial
        @cvode <node5> N_VEnableDotProdMulti_Serial
        @cvode <node5> N_VEnableLinearSumVectorArray_Serial
        @cvode <node5> N_VEnableScaleVectorArray_Serial
        @cvode <node5> N_VEnableConstVectorArray_Serial
        @cvode <node5> N_VEnableWrmsNormVectorArray_Serial
        @cvode <node5> N_VEnableWrmsNormMaskVectorArray_Serial
        @cvode <node5> N_VEnableScaleAddMultiVectorArray_Serial
        @cvode <node5> N_VEnableLinearCombinationVectorArray_Serial
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
      -> t
      -> unit

    (** Standard operations over array nvectors. *)
    module Ops : Nvector.NVECTOR_OPS with type t = t

    (** Standard operations over the underlying array. *)
    module DataOps : Nvector.NVECTOR_OPS with type t = data
  end

(** Produce a custom nvector from basic operations on
    an underlying array. *)
module Make : functor (A : sig
      type data
      val get       : data -> int -> float
      val set       : data -> int -> float -> unit
      val fill      : data -> float -> unit
      val make      : int -> float -> data
      val clone     : data -> data
      val length    : data -> int
    end) -> ARRAY_NVECTOR with type data = A.data

(** Nvector on {{:OCAML_DOC_ROOT(Array.html)} Array}s of [float]s. *)
include ARRAY_NVECTOR with type data = float array

