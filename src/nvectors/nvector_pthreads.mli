(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** The Pthreads nvectors of Sundials (requires pthreads).

    @nvector <NVector_links.html#the-nvector-pthreads-module> The NVECTOR_PTHREADS Module
    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** Pthreads nvectors are based on
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarrays} of floats. *)
type data = RealArray.t

(** Represents the internal layout of a Pthreads nvector.
    Pthreads nvectors can usually be used wherever serial nvectors can. *)
type kind = [`Pthreads|Nvector_serial.kind]

(** The type of Pthreads nvectors. *)
type t = (data, kind) Nvector.t

(** [make nthreads n iv] creates a new Pthreads nvector with [nthreads]
    threads and [n] elements inialized to [iv].

    The optional argument enables the fused and array operations for a given
    nvector (they are disabled by default).

    @nvector N_VNew_Pthreads
    @nvector N_VEnableFusedOps_Pthreads
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val make :
     ?context:Context.t
  -> ?with_fused_ops:bool
  -> int
  -> int
  -> float
  -> t

(** [wrap nthreads a] creates a new Pthreads nvector with [nthreads] threads
    over the elements of [a].

    The optional arguments permit to enable all the fused and array operations
    for a given nvector (they are disabled by default).

    @nvector N_VMake_Pthreads
    @nvector N_VEnableFusedOps_Pthreads
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val wrap :
     ?context:Context.t
  -> ?with_fused_ops:bool
  -> int
  -> RealArray.t
  -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> RealArray.t

(* TOPLEVEL-PRINTER: Nvector_pthreads.pp *)
(** Pretty-print a Pthreads nvector using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

(** Returns the number of threads used within a Pthreads nvector. *)
val num_threads : t -> int

(** Selectively enable or disable fused and array operations.
    The [with_fused_ops] argument enables or disables all such operations.

    @nvector N_VEnableFusedOps_Pthreads
    @nvector N_VEnableLinearCombination_Pthreads
    @nvector N_VEnableScaleAddMulti_Pthreads
    @nvector N_VEnableDotProdMulti_Pthreads
    @nvector N_VEnableLinearSumVectorArray_Pthreads
    @nvector N_VEnableScaleVectorArray_Pthreads
    @nvector N_VEnableConstVectorArray_Pthreads
    @nvector N_VEnableWrmsNormVectorArray_Pthreads
    @nvector N_VEnableWrmsNormMaskVectorArray_Pthreads
    @nvector N_VEnableScaleAddMultiVectorArray_Pthreads
    @nvector N_VEnableLinearCombinationVectorArray_Pthreads
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
    @since 4.0.0 *)
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

(** Underlyling nvector operations on Pthreads nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** A generic nvector interface to Pthreads nvectors.

    Create Pthreads nvectors using the generic nvector interface where the
    payload is wrapped with the {{!Nvector.gdata}RA} constructor. *)
module Any : sig (* {{{ *)

  (** [make nthreads n iv] creates a new Pthreads nvector with [nthreads]
      threads and [n] elements inialized to [iv].

      The optional argument enables the fused and array operations for a given
      nvector (they are disabled by default).

      @nvector N_VNew_Pthreads
      @nvector N_VEnableFusedOps_Pthreads
      @nvector N_VEnableLinearCombination_Pthreads
      @nvector N_VEnableScaleAddMulti_Pthreads
      @nvector N_VEnableDotProdMulti_Pthreads
      @nvector N_VEnableLinearSumVectorArray_Pthreads
      @nvector N_VEnableScaleVectorArray_Pthreads
      @nvector N_VEnableConstVectorArray_Pthreads
      @nvector N_VEnableWrmsNormVectorArray_Pthreads
      @nvector N_VEnableWrmsNormMaskVectorArray_Pthreads
      @nvector N_VEnableScaleAddMultiVectorArray_Pthreads
      @nvector N_VEnableLinearCombinationVectorArray_Pthreads
      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 5.0.0 *)
  val make :
       ?context                              : Context.t
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
    -> int
    -> int
    -> float
    -> Nvector.any

  (** [wrap nthreads a] creates a new Pthreads nvector with [nthreads] threads
      over the elements of [a].

      The optional arguments permit to enable all the fused and array operations
      for a given nvector (they are disabled by default).

      @nvector N_VMake_Pthreads
      @nvector N_VEnableFusedOps_Pthreads
      @nvector N_VEnableLinearCombination_Pthreads
      @nvector N_VEnableScaleAddMulti_Pthreads
      @nvector N_VEnableDotProdMulti_Pthreads
      @nvector N_VEnableLinearSumVectorArray_Pthreads
      @nvector N_VEnableScaleVectorArray_Pthreads
      @nvector N_VEnableConstVectorArray_Pthreads
      @nvector N_VEnableWrmsNormVectorArray_Pthreads
      @nvector N_VEnableWrmsNormMaskVectorArray_Pthreads
      @nvector N_VEnableScaleAddMultiVectorArray_Pthreads
      @nvector N_VEnableLinearCombinationVectorArray_Pthreads
      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 5.0.0 *)
  val wrap :
       ?context                              : Context.t
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
    -> int
    -> RealArray.t
    -> Nvector.any

  (** Returns the payload of the generic vector if it was constructed with
      {{!Nvector.gdata}RA} and an array of reals, otherwise
      raises {!Nvector.BadGenericType}. *)
  val unwrap : Nvector.any -> RealArray.t

  (** Selectively enable or disable fused and array operations.
      The [with_fused_ops] argument enables or disables all such operations.

      @nvector N_VEnableFusedOps_Pthreads
      @nvector N_VEnableLinearCombination_Pthreads
      @nvector N_VEnableScaleAddMulti_Pthreads
      @nvector N_VEnableDotProdMulti_Pthreads
      @nvector N_VEnableLinearSumVectorArray_Pthreads
      @nvector N_VEnableScaleVectorArray_Pthreads
      @nvector N_VEnableConstVectorArray_Pthreads
      @nvector N_VEnableWrmsNormVectorArray_Pthreads
      @nvector N_VEnableWrmsNormMaskVectorArray_Pthreads
      @nvector N_VEnableScaleAddMultiVectorArray_Pthreads
      @nvector N_VEnableLinearCombinationVectorArray_Pthreads
      @raise Nvector.BadGenericType If not called on a pthreads nvector
      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 4.0.0 *)
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
    -> Nvector.any
    -> unit

end (* }}} *)

