(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** A {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} Bigarray} of integers. *)
type t = (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t

(** [make n v] returns an array with [n] elements each set to [v]. *)
val make  : int -> int -> t

(** [create n] returns an uninitialized array with [n] elements. *)
val create  : int -> t

(** An array with no elements. *)
val empty : t

(** Pretty-print an array using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

(** Pretty-print an array using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module.
    The defaults are: [start="\["], [stop="\]"], [sep="; "], and
    [item=fun f->Format.fprintf f "%2d=% 6d"] (see
  {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}). *)
val ppi : ?start:string -> ?stop:string -> ?sep:string
          -> ?item:(Format.formatter -> int -> int -> unit)
          -> unit
          -> Format.formatter -> t -> unit

