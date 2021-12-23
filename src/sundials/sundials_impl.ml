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

let has_nvector_get_id =
  match Sundials_configuration.sundials_version with
  | 2,n,_ -> n >= 9
  | _ -> true

end

module Context = struct

  type cptr

  type t = { cptr : cptr }

  external c_make : unit -> cptr
    = "sunml_context_make"

  let make () = { cptr = c_make () }

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

end

