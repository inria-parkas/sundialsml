(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

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

