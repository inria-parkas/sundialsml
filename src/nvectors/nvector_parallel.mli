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

(** Parallel nvectors carry triples of a local
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarray} of floats,
    a global length, and an MPI communicator. *)
type data = Sundials.RealArray.t * int * Mpi.communicator

(** Represents the internal layout of a parallel nvector. *)
type kind

(** The type of parallel nvectors. *)
type t = (data, kind) Nvector.t

(** Raised by make if the given global length is not consistent with the sum of
    local lengths across all parallel instances. *)
exception IncorrectGlobalSize

(** [make nl ng c iv] creates a new parallel nvector with [nl] local elements,
    that is part of a global array with [ng] elements. The local elements are
    initialized to [iv], and communications occur on [c]. *)
val make : int -> int -> Mpi.communicator -> float -> t

(** Creates an nvector with a distinct underlying array but that shares the
    original global size and communicator. *)
val clone : t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> data

(** [local_array nv] returns the local array [a] underlying the parallel
    nvector [nv]. *)
val local_array : t -> Sundials.RealArray.t

(** Returns the number of local elements for a parallel nvector. *)
val local_length : t -> int

(** Returns the number of global elements for a parallel nvector. *)
val global_length : t -> int

(** Returns the communicator used for the parallel nvector. *)
val communicator : t -> Mpi.communicator

(** Produce a set of parallel {!Nvector.NVECTOR_OPS} from basic operations on
    an underlying array. *)
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

