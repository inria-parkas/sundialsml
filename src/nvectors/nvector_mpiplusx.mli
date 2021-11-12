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

(** The standard mpiplusx nvectors of Sundials. They couple an nvector with an
    mpi communicator.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @nocvode <node> NVECTOR_MPIPLUSX
    @since 5.0.0 *)

open Sundials

(** The data in underlying nvectors is exposed as an array of wrapped values
    together with an MPI communicator. *)
type data = Nvector.any * Mpi.communicator

(** Represents the internal layout of an mpiplusx nvector. *)
type kind

(** The type of mpiplusx nvectors. *)
type t = (data, kind) Nvector.t

(** A generic data wrapper for {!data}. *)
type Nvector.gdata += MpiPlusX of data

(** Creates an mpiplusx nvector from an mpi communicator and a generic nvector.

    An optional argument permits to enable all the fused and array operations
    for a given nvector (they are disabled by default).

    @since 5.0.0
    @cvode <node> N_VMake_MPIPlusX
    @cvode <node5> N_VEnableFusedOps_MPIManyVector *)
val wrap : ?with_fused_ops:bool -> Mpi.communicator -> Nvector.any -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> data

(** Returns the communicator used for the nvector.
    See also: {!Nvector_parallel.get_communicator} and
    {!Nvector_parallel.hide_communicator}. *)
val communicator : t -> Mpi.communicator

(** Selectively enable or disable fused and array operations.
    The [with_fused_ops] argument enables or disables all such operations.

    @cvode <node5> N_VEnableFusedOps_MPIManyVector
    @cvode <node5> N_VEnableLinearCombination_MPIManyVector
    @cvode <node5> N_VEnableScaleAddMulti_MPIManyVector
    @cvode <node5> N_VEnableDotProdMulti_MPIManyVector
    @cvode <node5> N_VEnableLinearSumVectorArray_MPIManyVector
    @cvode <node5> N_VEnableScaleVectorArray_MPIManyVector
    @cvode <node5> N_VEnableConstVectorArray_MPIManyVector
    @cvode <node5> N_VEnableWrmsNormVectorArray_MPIManyVector
    @cvode <node5> N_VEnableWrmsNormMaskVectorArray_MPIManyVector *)
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
  -> t
  -> unit

(** Underlying nvector operations on mpiplusx nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations implemented in OCaml on mpiplusx payloads. *)
module DataOps : Nvector.NVECTOR_OPS with type t = data

(** A generic nvector interface to mpiplusx nvectors.

    Create mpiplusx nvectors using the generic nvector interface where
    the payload is wrapped with the {{!Nvector.gdata}MpiPlusX} constructor. *)
module Any : sig

  (** Creates a generic nvector from an mpi communicator and a generic
      nvector.

      An optional argument permits to enable all the fused and array
      operations for a given nvector (they are disabled by default).

      @since 5.0.0
      @cvode <node> N_VMake_MPIPlusX
      @cvode <node5> N_VEnableFusedOps_MPIManyVector *)
  val wrap
    : ?with_fused_ops:bool -> Mpi.communicator -> Nvector.any -> Nvector.any

  (** Returns the payload of the generic vector if it was constructed with
      {{!Nvector.gdata}MpiPlusX}, otherwise raises {!Nvector.BadGenericType}. *)
  val unwrap : Nvector.any -> data

  (** Selectively enable or disable fused and array operations.
      The [with_fused_ops] argument enables or disables all such operations.

      @cvode <node5> N_VEnableFusedOps_MPIManyVector
      @cvode <node5> N_VEnableLinearCombination_MPIManyVector
      @cvode <node5> N_VEnableScaleAddMulti_MPIManyVector
      @cvode <node5> N_VEnableDotProdMulti_MPIManyVector
      @cvode <node5> N_VEnableLinearSumVectorArray_MPIManyVector
      @cvode <node5> N_VEnableScaleVectorArray_MPIManyVector
      @cvode <node5> N_VEnableConstVectorArray_MPIManyVector
      @cvode <node5> N_VEnableWrmsNormVectorArray_MPIManyVector
      @cvode <node5> N_VEnableWrmsNormMaskVectorArray_MPIManyVector
      @raise Nvector.BadGenericType If not called on an MPIMany nvector *)
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
    -> Nvector.any
    -> unit
end

