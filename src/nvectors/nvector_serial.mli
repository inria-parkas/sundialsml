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

(** Standard serial nvectors of Sundials.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvode <node7#ss:nvec_ser> NVECTOR_SERIAL *)

open Sundials

(** Serial nvectors are based on
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarrays} of floats. *)
type data = RealArray.t

(** Represents any nvector that can be treated as a serial nvector. That is,
    any nvector whose underlying elements can be accessed as an array
    locally. *)
type kind = [`Serial]

(** The type of serial nvectors. *)
type t = (data, kind) Nvector.t

(** The type of any nvector that can be used as a serial nvector. *)
type 'k any = (data, [>kind] as 'k) Nvector.t

(** [make n iv] creates a new serial nvector with [n] elements, each initialized
    to [iv].

    The optional argument enables the fused and array operations for a given
    nvector (they are disabled by default).

    @cvode <node5> N_VEnableFusedOps_Serial
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val make : ?with_fused_ops:bool -> int -> float -> t

(** [wrap a] creates a new serial nvector over the elements of [a].

    The optional arguments permit to enable all the fused and array operations
    for a given nvector (they are disabled by default).

    @cvode <node5> N_VEnableFusedOps_Serial
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val wrap : ?with_fused_ops:bool -> RealArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> RealArray.t

(* TOPLEVEL-PRINTER: Nvector_serial.pp *)
(** Pretty-print a serial nvector using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

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

(** Underlying nvector operations on serial nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations on {{!Sundials.RealArray}RealArray}s implemented in OCaml. *)
module DataOps : Nvector.NVECTOR_OPS with type t = RealArray.t

(** {2:genvec Generic nvector interface}

    Create serial nvectors using the generic nvector interface where the
    payload is wrapped with the {{!Nvector.gdata}RA} constructor. *)
module Any : sig (* {{{ *)

  (** [make n iv] creates a new serial nvector with [n] elements,
      each initialized to [iv].

      The optional arguments permit to enable fused and array operations for
      a given nvector (they are disabled by default).

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
      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 2.9.0 *)
  val make :
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
    -> int
    -> float
    -> Nvector.any

  (** [wrap a] creates a new serial nvector over the elements of [a].

      The optional arguments permit to enable fused and array operations for
      a given nvector (they are disabled by default).

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
      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 2.9.0 *)
  val wrap :
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
    -> RealArray.t
    -> Nvector.any

end (* }}} *)

