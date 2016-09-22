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


(* This code is taken from
   http://stackoverflow.com/questions/29336244/ocaml-automate-custom-pretty-printer-installing
 *)
let eval_str str =
  let lexbuf = Lexing.from_string str in
  let phrase = !Toploop.parse_toplevel_phrase lexbuf in
  Toploop.execute_phrase false Format.err_formatter phrase
;;

let install_printers ps =
  let install p = not (eval_str (Printf.sprintf "#install_printer %s;;" p)) in
  match List.filter install ps with
  | [] -> ()
  | ps ->
     (prerr_string "Errors installing the following printers:\n";
      List.iter (fun p -> Printf.fprintf stderr "  %s\n" p) ps;
      flush stderr)

let () =
  install_printers
    ["Sundials.RealArray.pp";
     "Sundials.RealArray2.pp";
     "Sundials.LintArray.pp";
     "Sundials.Roots.pp";
     "Sundials.RootDirs.pp";
    ]
