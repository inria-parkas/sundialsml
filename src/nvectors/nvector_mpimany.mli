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

(** The standard mpimany-vector nvectors of Sundials.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @nocvode <node> NVECTOR_MPIMANYVECTOR
    @since 5.0.0 *)

open Sundials

(** The data in underlying nvectors is exposed as an array of wrapped values
    together with an MPI communicator and the sum of their lengths. *)
type data = Nvector.any ROArray.t * Mpi.communicator * int

(** Represents the internal layout of an mpimany-vector nvector. *)
type kind

(** The type of mpimany-vector nvectors. *)
type t = (data, kind) Nvector.t

(** Generic wrapper for {!data}. *)
type Nvector.gdata += MpiMany of data

(** Creates a mpimany-vector nvector from an array of generic nvectors.
    If the communicator argument is not given then all nvectors in the array
    that have a communicator must have the same communicator which becomes the
    communicator of the new nvector (at least one element must have a
    communicator). If the communicator argument is given, each array element
    may use any or no communicator.

    An optional argument permits to enable all the fused and array operations
    for a given nvector (they are disabled by default).

    @since 5.0.0
    @cvode <node> N_VNew_MPIManyVector
    @cvode <node> N_VMake_MPIManyVector
    @cvode <node5> N_VEnableFusedOps_MPIManyVector
    @raise Invalid_arg if an mpi communicator is not specified or found. *)
val wrap
  : ?with_fused_ops:bool -> ?comm:Mpi.communicator -> Nvector.any ROArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> data

(** Returns the sum of the lengths of the component nvectors. *)
val length : t -> int

(** Returns the number of subectors in the array. *)
val num_subvectors : t -> int

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

(** Underlying nvector operations on mpimany-vector nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations implemented in OCaml on mpimany-vector payloads. *)
module DataOps : Nvector.NVECTOR_OPS with type t = data

(** A generic nvector interface to mpimany-vector nvectors.

    Create mpimany-vector nvectors using the generic nvector interface where
    the payload is wrapped with the {{!Nvector.gdata}MpiMany} constructor. *)
module Any : sig (* {{{ *)

  (** Creates a generic nvector from an array of generic nvectors.
      If the communicator argument is not given then all nvectors in the array
      that have a communicator must have the same communicator which becomes
      the communicator of the new nvector (at least one element must have a
      communicator). If the communicator argument is given, each array element
      may use any or no communicator.

      An optional argument permits to enable all the fused and array operations
      for a given nvector (they are disabled by default).

      @since 5.0.0
      @cvode <node> N_VNew_MPIManyVector
      @cvode <node> N_VMake_MPIManyVector
      @cvode <node> N_VEnableFusedOps_MPIManyVector
      @raise Invalid_arg if an mpi communicator is not specified or found. *)
  val wrap
    : ?with_fused_ops:bool -> ?comm:Mpi.communicator -> Nvector.any ROArray.t -> Nvector.any

  (** Returns the payload of the generic vector if it was constructed with
      {{!Nvector.gdata}MpiMany}, otherwise raises {!Nvector.BadGenericType}. *)
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

end (* }}} *)

