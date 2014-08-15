(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(** Interface to Parallel NVectors manipulated by the standard operations
    provided by Sundials. Requires MPI.
 
    @cvode <node7#ss:nvec_par> NVECTOR_PARALLEL *)
type kind

(** (local array, global length, MPI communicator) *)
type data = Sundials.RealArray.t * int * Mpi.communicator
type t = (data, kind) Sundials.nvector

(** Raised by make if the given global length is not consistent with the sum of
    local lengths across all parallel instances. *)
exception IncorrectGlobalSize

(** [make nl ng c iv] creates a new parallel nvector with [nl] local elements,
    that is part of a global array with [ng] elements. The local elements are
    initialized to [iv], and communications occur on [c]. *)
val make : int -> int -> Mpi.communicator -> float -> t

(** [unwrap nv] returns the local data [a] underlying the parallel nvector
    [nv]. *)
val unwrap : t -> Sundials.RealArray.t

(** Returns the number of global elements for a parallel nvector. *)
val global_length : t -> int

(** Returns the communicator used for by the parallel nvector. *)
val communicator : t -> Mpi.communicator

module Ops : sig
  val n_vclone        : t -> t
  val n_vlinearsum    : float -> t -> float -> t -> t -> unit
  val n_vconst        : float -> t -> unit
  val n_vprod         : t -> t -> t -> unit
  val n_vdiv          : t -> t -> t -> unit
  val n_vscale        : float -> t -> t -> unit
  val n_vabs          : t -> t -> unit
  val n_vinv          : t -> t -> unit
  val n_vaddconst     : t -> float -> t -> unit
  val n_vdotprod      : t -> t -> float
  val n_vmaxnorm      : t -> float
  val n_vwrmsnorm     : t -> t -> float
  val n_vwrmsnormmask : t -> t -> t -> float
  val n_vmin          : t -> float
  val n_vwl2norm      : t -> t -> float
  val n_vl1norm       : t -> float
  val n_vcompare      : float -> t -> t -> unit
  val n_vinvtest      : t -> t -> bool
  val n_vconstrmask   : t -> t -> t -> bool
  val n_vminquotient  : t -> t -> float
end

