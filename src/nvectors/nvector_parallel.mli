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

(** The standard parallel nvectors of Sundials (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvode <node7#ss:nvec_par> NVECTOR_PARALLEL *)

open Sundials

(** Parallel nvectors carry triples of a local
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarray} of floats,
    a global length, and an MPI communicator. *)
type data = RealArray.t * int * Mpi.communicator

(** Represents the internal layout of a parallel nvector. *)
type kind

(** The type of parallel nvectors. *)
type t = (data, kind) Nvector.t

(** Generic wrapper for {!data}. *)
type Nvector.gdata += Par of data

(** Raised by make if the given global length is not consistent with the sum of
    local lengths across all parallel instances. *)
exception IncorrectGlobalSize

(** [make nl ng c iv] creates a new parallel nvector with [nl] local elements,
    that is part of a global array with [ng] elements. The local elements are
    initialized to [iv], and communications occur on [c].

    The optional argument enables the fused and array operations for a given
    nvector (they are disabled by default).

    @cvode <node5> N_VEnableFusedOps_Parallel
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val make : ?with_fused_ops:bool -> int -> int -> Mpi.communicator -> float -> t

(** Creates an nvector with a distinct underlying array but that shares the
    original global size and communicator. *)
val clone : t -> t

(** [wrap a] creates a new parallel nvector from [a].

    The optional arguments permit to enable all the fused and array operations
    for a given nvector (they are disabled by default).

    @cvode <node5> N_VEnableFusedOps_Parallel
    @raise Config.NotImplementedBySundialsVersion Fused and array operations not available. *)
val wrap : ?with_fused_ops:bool -> data -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> data

(* TOPLEVEL-PRINTER: Nvector_parallel.pp *)
(** Pretty-print the local portion of a parallel nvector using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

(** [local_array nv] returns the local array [a] underlying the parallel
    nvector [nv]. *)
val local_array : t -> RealArray.t

(** Returns the number of local elements for a parallel nvector. *)
val local_length : t -> int

(** Returns the number of global elements for a parallel nvector. *)
val global_length : t -> int

(** Returns the communicator used for the parallel nvector. *)
val communicator : t -> Mpi.communicator

(** Return the communicator associated with any nvector.

    @nocvode <node> N_VGetCommunicator
    @since 5.0.0 *)
val get_communicator : ('d, 'k) Nvector.t -> Mpi.communicator option

(** Hides an MPI communicator for use in custom nvector functions. *)
val hide_communicator : Mpi.communicator -> Nvector_custom.communicator

(** Selectively enable or disable fused and array operations.
    The [with_fused_ops] argument enables or disables all such operations.

    @since 4.0.0
    @cvode <node5> N_VEnableFusedOps_Parallel
    @cvode <node5> N_VEnableLinearCombination_Parallel
    @cvode <node5> N_VEnableScaleAddMulti_Parallel
    @cvode <node5> N_VEnableDotProdMulti_Parallel
    @cvode <node5> N_VEnableLinearSumVectorArray_Parallel
    @cvode <node5> N_VEnableScaleVectorArray_Parallel
    @cvode <node5> N_VEnableConstVectorArray_Parallel
    @cvode <node5> N_VEnableWrmsNormVectorArray_Parallel
    @cvode <node5> N_VEnableWrmsNormMaskVectorArray_Parallel
    @cvode <node5> N_VEnableScaleAddMultiVectorArray_Parallel
    @cvode <node5> N_VEnableLinearCombinationVectorArray_Parallel
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

(** Produce a set of parallel {!Nvector.NVECTOR_OPS} from basic
    operations on an underlying array. *)
module MakeOps : functor (A : sig
      type local_data
      val get       : local_data -> int -> float
      val set       : local_data -> int -> float -> unit
      val fill      : local_data -> float -> unit
      val make      : int -> float -> local_data
      val clone     : local_data -> local_data
      val length    : local_data -> int
    end) -> Nvector.NVECTOR_OPS
            with type t = A.local_data * int * Mpi.communicator

(** Underlying nvector operations on parallel nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations on {!data} implemented in OCaml. *)
module DataOps : Nvector.NVECTOR_OPS with type t = data

(** A generic nvector interface to parallel nvectors.

    Create parallel nvectors using the generic nvector interface where the
    payload is wrapped with the {{!Nvector.gdata}Par} constructor. *)
module Any : sig (* {{{ *)

  (** [make nl ng c iv] creates a new parallel nvector with [nl] local elements,
      that is part of a global array with [ng] elements. The local elements are
      initialized to [iv], and communications occur on [c].

      The optional argument enables the fused and array operations for a given
      nvector (they are disabled by default).

      @cvode <node5> N_VEnableFusedOps_Parallel
      @cvode <node5> N_VEnableLinearCombination_Parallel
      @cvode <node5> N_VEnableScaleAddMulti_Parallel
      @cvode <node5> N_VEnableDotProdMulti_Parallel
      @cvode <node5> N_VEnableLinearSumVectorArray_Parallel
      @cvode <node5> N_VEnableScaleVectorArray_Parallel
      @cvode <node5> N_VEnableConstVectorArray_Parallel
      @cvode <node5> N_VEnableWrmsNormVectorArray_Parallel
      @cvode <node5> N_VEnableWrmsNormMaskVectorArray_Parallel
      @cvode <node5> N_VEnableScaleAddMultiVectorArray_Parallel
      @cvode <node5> N_VEnableLinearCombinationVectorArray_Parallel
      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 5.0.0 *)
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
    -> Mpi.communicator
    -> float
    -> Nvector.any

  (** [wrap a] creates a new parallel nvector from [a].

      The optional arguments permit to enable all the fused and array operations
      for a given nvector (they are disabled by default).

      @cvode <node5> N_VEnableFusedOps_Parallel
      @cvode <node5> N_VEnableLinearCombination_Parallel
      @cvode <node5> N_VEnableScaleAddMulti_Parallel
      @cvode <node5> N_VEnableDotProdMulti_Parallel
      @cvode <node5> N_VEnableLinearSumVectorArray_Parallel
      @cvode <node5> N_VEnableScaleVectorArray_Parallel
      @cvode <node5> N_VEnableConstVectorArray_Parallel
      @cvode <node5> N_VEnableWrmsNormVectorArray_Parallel
      @cvode <node5> N_VEnableWrmsNormMaskVectorArray_Parallel
      @cvode <node5> N_VEnableScaleAddMultiVectorArray_Parallel
      @cvode <node5> N_VEnableLinearCombinationVectorArray_Parallel
      @raise Config.NotImplementedBySundialsVersion Fused and array operations not available.
      @since 5.0.0 *)
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
    -> data
    -> Nvector.any

  (** Returns the payload of the generic vector if it was constructed with
      {{!Nvector.gdata}Par} and a parallel payload, otherwise
      raises {!Nvector.BadGenericType}. *)
  val unwrap : Nvector.any -> data

  (** Selectively enable or disable fused and array operations.
      The [with_fused_ops] argument enables or disables all such operations.

      @since 4.0.0
      @cvode <node5> N_VEnableFusedOps_Parallel
      @cvode <node5> N_VEnableLinearCombination_Parallel
      @cvode <node5> N_VEnableScaleAddMulti_Parallel
      @cvode <node5> N_VEnableDotProdMulti_Parallel
      @cvode <node5> N_VEnableLinearSumVectorArray_Parallel
      @cvode <node5> N_VEnableScaleVectorArray_Parallel
      @cvode <node5> N_VEnableConstVectorArray_Parallel
      @cvode <node5> N_VEnableWrmsNormVectorArray_Parallel
      @cvode <node5> N_VEnableWrmsNormMaskVectorArray_Parallel
      @cvode <node5> N_VEnableScaleAddMultiVectorArray_Parallel
      @cvode <node5> N_VEnableLinearCombinationVectorArray_Parallel
      @raise Nvector.BadGenericType If not called on a parallel nvector
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
    -> Nvector.any
    -> unit
end (* }}} *)

