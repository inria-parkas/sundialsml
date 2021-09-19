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

    @since 5.0.0 *)
val wrap : Mpi.communicator -> Nvector.any -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> data

(** Returns the communicator used for the nvector.
    See also: {!Nvector_parallel.get_communicator} and
    {!Nvector_parallel.hide_communicator}. *)
val communicator : t -> Mpi.communicator

(** Underlying nvector operations on mpiplusx nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations implemented in OCaml on mpiplusx payloads. *)
module DataOps : Nvector.NVECTOR_OPS with type t = data

(** {2:genvec Generic nvector interface}

    Create mpiplusx nvectors using the generic nvector interface where
    the payload is wrapped with the {{!Nvector.gdata}MpiPlusX} constructor. *)
module Any : sig

  (** Creates a generic nvector from an mpi communicator and a generic
      nvector.

      @since 5.0.0 *)
  val wrap : Mpi.communicator -> Nvector.any -> Nvector.any

  (** Returns the payload of the generic vector if it was constructed with
      {{!Nvector.gdata}MpiPlusX}, otherwise raises {!Nvector.BadGenericType}. *)
  val unwrap : Nvector.any -> data

end

