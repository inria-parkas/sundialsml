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

