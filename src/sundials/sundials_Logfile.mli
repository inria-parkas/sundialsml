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

(** A rudimentary interface to C streams for logging in Sundials. *)

(** An open log file. *)
type t

(** The stderr file. *)
val stderr : t

(** The stdout file. *)
val stdout : t

(** Opens the named file. When [trunc] is false, the default, writes are
    appended to the file. When [trunc] is true, the opened file is
    truncated to zero length. Files are closed on garbage collection. *)
val openfile : ?trunc:bool -> string -> t

(** Writes the given string to an open log file. *)
val output_string : t -> string -> unit

(** Writes the given byte sequence to an open log file. *)
val output_bytes : t -> bytes -> unit

(** Flushes the given file. *)
val flush : t -> unit

(** Closes the given file. *)
val close : t -> unit

