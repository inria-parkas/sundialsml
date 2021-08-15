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

(** A generic data wrapper for {!data}. *)
type Nvector.gdata += MpiMany of data

(** Creates a mpimany-vector nvector from an array of generic nvectors.
    The parallel array elements must all use the same communicator which
    becomes the communicator of the resulting nvector. Specifying the
    communicator explicitly is useful if none of the array elements has
    one.

    @since 5.0.0
    @raises Invalid_arg if an mpi communicator is not specified or found. *)
val wrap : ?comm:Mpi.communicator -> Nvector.any ROArray.t -> t

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

(** Underlying nvector operations on mpimany-vector nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations implemented in OCaml on mpimany-vector payloads. *)
module DataOps : Nvector.NVECTOR_OPS with type t = data

(** {2:genvec Generic nvector interface}

    Create mpimany-vector nvectors using the generic nvector interface where
    the payload is wrapped with the {{!Nvector.gdata}MpiMany} constructor. *)
module Any : sig

  (** Creates a generic nvector from an array of generic nvectors.
      The parallel array elements must all use the same communicator which
      becomes the communicator of the resulting nvector. Specifying the
      communicator explicitly is useful if none of the array elements has
      one.

      @since 5.0.0
      @raises Invalid_arg if an mpi communicator is not specified or found. *)
  val wrap : ?comm:Mpi.communicator -> Nvector.any ROArray.t -> Nvector.any

end

