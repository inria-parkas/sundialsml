(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2021 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

module Profiler = struct
  module I = Sundials_impl.Profiler

  external make : Mpi.communicator -> string -> I.t
    = "sunml_profiler_make_parallel"

end

module Context = struct

  module I = Sundials_impl.Context

  external c_make : Mpi.communicator -> I.cptr
    = "sunml_context_make_parallel"

  let make ?profiler comm =
    let ctx = { I.cptr = c_make comm; I.profiler = None } in
    (match profiler with
     | Some p -> Sundials.Context.set_profiler ctx p
     | None ->
         if not Sundials_configuration.caliper_enabled
         then Sundials.Context.set_profiler ctx
                (Profiler.make comm "SUNContext Default"));
    ctx

end

