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
    = "sunml_profiler_create_parallel"

end

module Logger = struct
  module I = Sundials_impl.Logger

  external make : Mpi.communicator -> int -> I.t
    = "sunml_logger_create_parallel"

  external get_output_rank : I.t -> int
    = "sunml_logger_get_output_rank"

end

module Context = struct

  module I = Sundials_impl.Context

  external c_make : Mpi.communicator -> I.cptr * Sundials_impl.Logger.t
    = "sunml_context_create_parallel"

  let make ?profiler ?logger comm =
    let cptr, original_logger = c_make comm in
    let ctx = { I.cptr = cptr; I.profiler = None;
                I.logger = original_logger }
    in
    (match profiler with
     | Some p -> Sundials.Context.set_profiler ctx p
     | None ->
         if not Sundials_configuration.caliper_enabled
         then Sundials.Context.set_profiler ctx
                (Profiler.make comm "SUNContext Default"));
    (match logger with Some l -> Sundials.Context.set_logger ctx l | None -> ());
    ctx

end

