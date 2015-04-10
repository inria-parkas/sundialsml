(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Kinsol_impl
include SlsTypes

(* Must correspond with kin_klu_ordering_tag *)
type ordering =
     Amd
   | ColAmd
   | Natural

external c_klu : serial_session -> int -> int -> unit
  = "c_kin_klu_init"

let klu f nnz session onv =
  (match onv with
   | Some nv -> session.neqs <- Sundials.RealArray.length (Nvector.unwrap nv)
   | None -> ());
  session.ls_callbacks <- SlsKluCallback { jacfn = f; smat = None };
  c_klu session neqs nnz

external c_set_ordering : serial_session -> ordering -> unit
  = "c_kin_klu_set_ordering"

let set_ordering session ordering =
  ls_check_klu session;
  c_set_ordering session ordering

external c_reinit : serial_session -> int -> int -> bool -> unit
  = "c_kin_klu_reinit"

let reinit session n nnz realloc =
  ls_check_klu session;
  c_reinit session n nnz realloc

external c_get_num_jac_evals : serial_session -> int
  = "c_kin_klu_get_num_jac_evals"

let get_num_jac_evals session =
  ls_check_klu session;
  c_get_num_jac_evals session

