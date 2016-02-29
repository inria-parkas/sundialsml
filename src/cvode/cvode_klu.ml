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

open Cvode_impl
include SlsTypes

(* Must correspond with cvode_klu_ordering_tag *)
type ordering =
     Amd
   | ColAmd
   | Natural

external c_klu : 'k serial_session -> int -> int -> unit
  = "c_cvode_klu_init"

let klu f nnz session nv =
  let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
  session.ls_callbacks <- SlsKluCallback { jacfn = f; smat = None };
  session.ls_precfns <- NoPrecFns;
  c_klu session neqs nnz

external c_set_ordering : 'k serial_session -> ordering -> unit
  = "c_cvode_klu_set_ordering"

let set_ordering session ordering =
  ls_check_klu session;
  c_set_ordering session ordering

external c_reinit : 'k serial_session -> int -> int -> bool -> unit
  = "c_cvode_klu_reinit"

let reinit session n nnz realloc =
  ls_check_klu session;
  c_reinit session n nnz realloc

external c_get_num_jac_evals : 'k serial_session -> int
  = "c_cvode_klu_get_num_jac_evals"

let get_num_jac_evals session =
  ls_check_klu session;
  c_get_num_jac_evals session

