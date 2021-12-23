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

(** Contexts for creating Sundials values using MPI

    Provides an alternative to {!Sundials_context.make} that allows
    specifying an MPI communciator.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)

    @cvode <node> SUNContext
    @since 6.0.0 *)

(** Create a new context.

    @cvode <node> SUNContext_Create *)
val make : Mpi.communicator -> Sundials_Context.t

