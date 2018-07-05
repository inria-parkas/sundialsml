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
  = "c_sundials_stderr"

external c_stdout : unit -> t
  = "c_sundials_stdout"

external fopen : string -> bool -> t
  = "c_sundials_fopen"

let stderr = c_stderr ()
let stdout = c_stdout ()

let openfile ?(trunc=false) fpath = fopen fpath trunc

external flush : t -> unit
  = "c_sundials_fflush"

