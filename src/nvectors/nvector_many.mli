(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2021 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** The standard many-vector nvectors of Sundials.

    @nvector <NVector_links.html#the-nvector-manyvector-module> The NVECTOR_MANYVECTOR Module
    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @since 5.0.0 *)

open Sundials

(** The data in underlying nvectors is exposed as an array of wrapped values
    paired with the sum of their lengths. *)
type data = Nvector.any ROArray.t * int

(** Represents the internal layout of a many-vector nvector. *)
type kind

(** The type of many-vector nvectors. *)
type t = (data, kind) Nvector.t

(** Generic wrapper for {!data}. *)
type Nvector.gdata += Many of data

(** Creates a many-vector nvector from an array of generic nvectors.

    @nvector N_VNew_ManyVector
    @since 5.0.0 *)
val wrap : ?context:Context.t -> Nvector.any ROArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> data

(** Returns the sum of the lengths of the component nvectors. *)
val length : t -> int

(** Returns the number of subectors in the array. *)
val num_subvectors : t -> int

(** Selectively enable or disable fused and array operations.
    The [with_fused_ops] argument enables or disables all such operations.

    @nvector N_VEnableFusedOps_ManyVector
    @nvector N_VEnableLinearCombination_ManyVector
    @nvector N_VEnableScaleAddMulti_ManyVector
    @nvector N_VEnableDotProdMulti_ManyVector
    @nvector N_VEnableLinearSumVectorArray_ManyVector
    @nvector N_VEnableScaleVectorArray_ManyVector
    @nvector N_VEnableConstVectorArray_ManyVector
    @nvector N_VEnableWrmsNormVectorArray_ManyVector
    @nvector N_VEnableWrmsNormMaskVectorArray_ManyVector
    @nvector N_VEnableDotProdMultiLocal_ManyVector *)
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
  -> ?with_dot_prod_multi_local            : bool
  -> t
  -> unit

(** Underlying nvector operations on many-vector nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations implemented in OCaml on many-vector payloads. *)
module DataOps : Nvector.NVECTOR_OPS with type t = data

(** A generic nvector interface to many-vector nvectors.

    Create many-vector nvectors using the generic nvector interface where the
    payload is wrapped with the {{!Nvector.gdata}Many} constructor. *)
module Any : sig (* {{{ *)

  (** Creates a generic nvector from an array of generic nvectors.

      @nvector N_VNew_ManyVector
      @since 5.0.0 *)
  val wrap : ?context:Context.t -> Nvector.any ROArray.t -> Nvector.any

  (** Returns the payload of the generic vector if it was constructed with
      {{!Nvector.gdata}Many}, otherwise raises {!Nvector.BadGenericType}. *)
  val unwrap : Nvector.any -> data

  (** Selectively enable or disable fused and array operations.
      The [with_fused_ops] argument enables or disables all such operations.

      @nvector N_VEnableFusedOps_ManyVector
      @nvector N_VEnableLinearCombination_ManyVector
      @nvector N_VEnableScaleAddMulti_ManyVector
      @nvector N_VEnableDotProdMulti_ManyVector
      @nvector N_VEnableLinearSumVectorArray_ManyVector
      @nvector N_VEnableScaleVectorArray_ManyVector
      @nvector N_VEnableConstVectorArray_ManyVector
      @nvector N_VEnableWrmsNormVectorArray_ManyVector
      @nvector N_VEnableWrmsNormMaskVectorArray_ManyVector
      @nvector N_VEnableDotProdMultiLocal_ManyVector
      @raise Nvector.BadGenericType If not called on a many nvector *)
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
    -> ?with_dot_prod_multi_local            : bool
    -> Nvector.any
    -> unit

end (* }}} *)

