(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2020 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(* Global internal definitions *)

external crash : string -> 'a = "sunml_crash"

(* A simple way of sharing values between OCaml and C. *)
module Vptr : sig

  (* A vptr ("value pointer") involves three elements:
     1. val: the OCaml value being shared between OCaml and C
     2. croot : a malloc'ed pointer to val registered as a global root
     3. vptr: pairs val with a custom value wrapping croot

     From OCaml, val is accessed via (unwrap vptr).
     From C, val is accessed via croot (VPTRCROOT(val), see sundials_ml.h).

     When no more references are held to vptr, it may be collected by the GC.
     When this happens, its second component (the custom value) is also
     collected (since there can be no other references to it), and its
     finalizer removes the global root and frees croot.

     The structure referenced by val must not reference vptr as this would
     create an inter-heap cycle (val -> vptr -> croot -> val) and prevent
     garbage collection. The block referenced by val must no longer be
     accessed through croot after garbage collection. This means that the
     OCaml program must reference vptr for as long as C functions may try to
     go through croot.
                                    .
                 OCaml GCed heap    .   C heap (malloc'ed values)
                                    .
                  -------+          .
                         |          .
                        \|/         .
                        +---+---+   .
                   vptr:| * | * |   .
                        +-|-+--\+   .
                         /      \   .
                        /        +------------+
                       /            .         |
                     \|/            .        \|/
                   +----+           .       +---+
               val:| 'a |           . croot:| * | (global root)
                   +----+           .       +-|-+
                     /|\            .         |
                      +-----------------------+
                                    .
  *)
  type 'a vptr

  val make   : 'a -> 'a vptr

  val unwrap : 'a vptr -> 'a

end = struct

  type 'a vcptr
  type 'a vptr = 'a * 'a vcptr

  external make : 'a -> 'a vptr
    = "sunml_make_vptr"

  let unwrap (x, _) = x

end

(* Callbacks into Sundials combine a C function pointer and an OCaml function
   to invoke it. *)
module Callback = struct

  type 'f cfunptr

  type 'f cfun = {
    cptr : 'f cfunptr;
    call : 'f;
  }

  let invoke { call; _ } = call

end

(* Compatiblity modes and Sundials versions *)
module Version = struct

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode2 =
  match Sundials_configuration.sundials_version with
  | 2,_,_ -> true
  | _ -> false

(* "Simulate" Nonlinear Solvers in Sundials < 4.0.0 *)
let in_compat_mode2_3 =
  match Sundials_configuration.sundials_version with
  | 2,_,_ -> true
  | 3,_,_ -> true
  | _ -> false

let lt400 =
  match Sundials_configuration.sundials_version with
  | 2,_,_ -> true
  | 3,_,_ -> true
  | _ -> false

let lt500 =
  let m, _, _ = Sundials_configuration.sundials_version in
  m < 5

let lt530 =
  let m, n, _ = Sundials_configuration.sundials_version in
  m < 5 || (m = 5 && n < 3)

let lt540 =
  let m, n, _ = Sundials_configuration.sundials_version in
  m < 5 || (m = 5 && n < 4)

let lt580 =
  let m, n, _ = Sundials_configuration.sundials_version in
  m < 5 || (m = 5 && n < 8)

let lt600 =
  let m, _, _ = Sundials_configuration.sundials_version in
  m < 6

let lt620 =
  let m, n, _ = Sundials_configuration.sundials_version in
  m < 6 || (m = 6 && n < 2)

let lt640 =
  let m, n, _ = Sundials_configuration.sundials_version in
  m < 6 || (m = 6 && n < 4)

let has_nvector_get_id =
  match Sundials_configuration.sundials_version with
  | 2,n,_ -> n >= 9
  | _ -> true

end

module Logfile = struct
  type t

  external c_stderr : unit -> t
    = "sunml_sundials_stderr"

  external c_stdout : unit -> t
    = "sunml_sundials_stdout"

  external fopen : string -> bool -> t
    = "sunml_sundials_fopen"

  let stderr = c_stderr ()
  let stdout = c_stdout ()

  let openfile ?(trunc=false) fpath = fopen fpath trunc

  external output_string : t -> string -> unit
    = "sunml_sundials_write"

  external output_bytes : t -> bytes -> unit
    = "sunml_sundials_write"

  external flush : t -> unit
    = "sunml_sundials_fflush"

  external close : t -> unit
    = "sunml_sundials_close"
end

module Profiler = struct
  type t

  external make : string -> t
    = "sunml_profiler_make"
end

module Logger = struct
  type t
end

module Context = struct

  type cptr

  type t = {
    cptr : cptr;
    mutable profiler : Profiler.t option;
    mutable logger : Logger.t;
  }

  exception ExternalProfilerInUse

  external c_make : unit -> cptr * Logger.t
    = "sunml_context_make"

  external c_set_profiler : cptr -> Profiler.t -> unit
    = "sunml_context_set_profiler"

  let set_profiler ({ cptr; _ } as context) profiler =
    c_set_profiler cptr profiler;
    context.profiler <- Some profiler

  external c_set_logger : cptr -> Logger.t -> unit
    = "sunml_context_set_logger"

  let set_logger ({ cptr; _ } as context) logger =
    c_set_logger cptr logger;
    context.logger <- logger

  let get_logger { logger; _} = logger

  let make ?profiler ?logger () =
    let cptr, original_logger = c_make () in
    let ctx = { cptr; profiler = None; logger = original_logger } in
    (match profiler with
     | Some p -> set_profiler ctx p
     | None ->
         if not Sundials_configuration.caliper_enabled
         then set_profiler ctx (Profiler.make "SUNContext Default"));
    (match logger with Some l -> set_logger ctx l | None -> ());
    ctx

  let default_context = (Weak.create 1 : t Weak.t)

  let default () =
    match Weak.get default_context 0 with
    | Some c -> c
    | None ->
        let ctx = make () in
        Weak.set default_context 0 (Some ctx);
        ctx

  let get = function
    | None -> default ()
    | Some ctx -> ctx

  let get_profiler { profiler; _ } =
    match profiler with
    | Some p -> p
    | None -> raise ExternalProfilerInUse

end

