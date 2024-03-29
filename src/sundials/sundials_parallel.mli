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

(** Generic definitions for use with MPI.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS) *)

(** Performance profiling using MPI

    @profiler SUNProfiler
    @since 6.0.0 *)
module Profiler : sig

  (** Creates a new profiler with the given name.

      @profiler SUNProfiler_Create *)
  val make : Mpi.communicator -> string -> Sundials.Profiler.t

end

(** Contexts for creating Sundials values using MPI

    Provides an alternative to {!Sundials.Context.make} that allows
    specifying an MPI communciator.

    @context SUNContext
    @since 6.0.0 *)
module Context : sig

  (** Create a new context.

      @context SUNContext_Create *)
  val make :
       ?profiler:Sundials.Profiler.t
    -> Mpi.communicator
    -> Sundials.Context.t

end

