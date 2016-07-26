(*---------------------------------------------------------------------*)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2016 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file, unlike most      *)
(*  files constituting this library, is distributed under the Q        *)
(*  Public License, refer to the file LICENSE.                         *)
(*                                                                     *)
(*---------------------------------------------------------------------*)

(** Utility functions for setting up the toplevel.

 @version VERSION()
 @author Timothy Bourke (Inria/ENS)
 @author Jun Inoue (Inria/ENS)
 @author Marc Pouzet (UPMC/ENS/Inria)
 *)

(** Execute a toplevel directive, such as #install_printer.  The
    argument should contain the full directive including the terminating
    ;; token.  *)
val eval_str : string -> bool

(** Perform #install_printer on each function named by the input
    strings.  *)
val install_printers : string list -> unit
