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

(** The OpenMP nvectors of Sundials (requires OpenMP).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvode <node7#ss:nvec_par> NVECTOR_PARALLEL *)

open Sundials

(** OpenMP nvectors are based on
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarrays} of floats. *)
type data = RealArray.t

(** Represents the internal layout of an OpenMP nvector.
    OpenMP nvectors can usually be used wherever serial nvectors can. *)
type kind = [`OpenMP|Nvector_serial.kind]

(** The type of OpenMP nvectors. *)
type t = (data, kind) Nvector.t

(** [make nthreads n iv] creates a new OpenMP nvector with [nthreads]
    threads and [n] elements inialized to [iv].

    The optional arguments permit to enable all the fused and array operations
    for a given nvector (they are disabled by default).

    @cvode <node5> N_VEnableFusedOps_OpenMP
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val make : ?with_fused_ops:bool -> int -> int -> float -> t

(** [wrap nthreads a] creates a new OpenMP nvector with [nthreads] threads
    over the elements of [a].

    The optional arguments permit to enable all the fused and array operations
    for a given nvector (they are disabled by default).

    @cvode <node5> N_VEnableFusedOps_OpenMP
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val wrap : ?with_fused_ops:bool -> int -> RealArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> RealArray.t

(* TOPLEVEL-PRINTER: Nvector_openmp.pp *)
(** Pretty-print an OpenMP nvector using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

(** Returns the number of threads used within an OpenMP nvector. *)
val num_threads : t -> int

(** Selectively enable or disable fused and array operations.
    The [with_fused_ops] argument enables or disables all such operations.

    @since 4.0.0
    @cvode <node5> N_VEnableFusedOps_OpenMP
    @cvode <node5> N_VEnableLinearCombination_OpenMP
    @cvode <node5> N_VEnableScaleAddMulti_OpenMP
    @cvode <node5> N_VEnableDotProdMulti_OpenMP
    @cvode <node5> N_VEnableLinearSumVectorArray_OpenMP
    @cvode <node5> N_VEnableScaleVectorArray_OpenMP
    @cvode <node5> N_VEnableConstVectorArray_OpenMP
    @cvode <node5> N_VEnableWrmsNormVectorArray_OpenMP
    @cvode <node5> N_VEnableWrmsNormMaskVectorArray_OpenMP
    @cvode <node5> N_VEnableScaleAddMultiVectorArray_OpenMP
    @cvode <node5> N_VEnableLinearCombinationVectorArray_OpenMP
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

(** Underlying nvector operations on OpenMP nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** {2:genvec Generic nvector interface}

    Create OpenMP nvectors using the generic nvector interface where the
    payload is wrapped with the {{!Nvector.gdata}RA} constructor. *)
module Any : sig (* {{{ *)

  (** [make nthreads n iv] creates a new OpenMP nvector with [nthreads]
      threads and [n] elements inialized to [iv].

      The optional arguments permit to enable all the fused and array operations
      for a given nvector (they are disabled by default).

      @cvode <node5> N_VEnableFusedOps_OpenMP
      @cvode <node5> N_VEnableLinearCombination_OpenMP
      @cvode <node5> N_VEnableScaleAddMulti_OpenMP
      @cvode <node5> N_VEnableDotProdMulti_OpenMP
      @cvode <node5> N_VEnableLinearSumVectorArray_OpenMP
      @cvode <node5> N_VEnableScaleVectorArray_OpenMP
      @cvode <node5> N_VEnableConstVectorArray_OpenMP
      @cvode <node5> N_VEnableWrmsNormVectorArray_OpenMP
      @cvode <node5> N_VEnableWrmsNormMaskVectorArray_OpenMP
      @cvode <node5> N_VEnableScaleAddMultiVectorArray_OpenMP
      @cvode <node5> N_VEnableLinearCombinationVectorArray_OpenMP
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
    -> int
    -> float
    -> Nvector.any

  (** [wrap nthreads a] creates a new OpenMP nvector with [nthreads] threads
      over the elements of [a].

      The optional arguments permit to enable all the fused and array operations
      for a given nvector (they are disabled by default).

      @cvode <node5> N_VEnableFusedOps_OpenMP
      @cvode <node5> N_VEnableLinearCombination_OpenMP
      @cvode <node5> N_VEnableScaleAddMulti_OpenMP
      @cvode <node5> N_VEnableDotProdMulti_OpenMP
      @cvode <node5> N_VEnableLinearSumVectorArray_OpenMP
      @cvode <node5> N_VEnableScaleVectorArray_OpenMP
      @cvode <node5> N_VEnableConstVectorArray_OpenMP
      @cvode <node5> N_VEnableWrmsNormVectorArray_OpenMP
      @cvode <node5> N_VEnableWrmsNormMaskVectorArray_OpenMP
      @cvode <node5> N_VEnableScaleAddMultiVectorArray_OpenMP
      @cvode <node5> N_VEnableLinearCombinationVectorArray_OpenMP
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
    -> int
    -> RealArray.t
    -> Nvector.any

end (* }}} *)

