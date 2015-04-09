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

(* Must correspond with cvode_superlumt_ordering_tag *)
type ordering =
     Natural
   | MinDegreeProd
   | MinDegreeSum
   | ColAmd

external c_superlumt : serial_session -> int -> int -> int -> unit
  = "c_cvode_superlumt_init"

let superlumt f ~nnz ~nthreads session nv =
  let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
  session.ls_callbacks <- SlsSuperlumtCallback { jacfn = f; smat = None };
  c_superlumt session neqs nnz nthreads

external c_set_ordering : serial_session -> ordering -> unit
  = "c_cvode_superlumt_set_ordering"

let set_ordering session ordering =
  ls_check_superlumt session;
  c_set_ordering session ordering

external c_get_num_jac_evals : serial_session -> int
  = "c_cvode_superlumt_get_num_jac_evals"

let get_num_jac_evals session =
  ls_check_superlumt session;
  c_get_num_jac_evals session

