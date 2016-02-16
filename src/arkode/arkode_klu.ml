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

open Arkode_impl
include SlsTypes

(* Must correspond with arkode_klu_ordering_tag *)
type ordering =
     Amd
   | ColAmd
   | Natural

external c_klu : 'k serial_session -> int -> int -> unit
  = "c_arkode_klu_init"

let klu f nnz session nv =
  let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
  session.ls_callbacks <- SlsKluCallback { jacfn = f; smat = None };
  session.ls_precfns <- NoPrecFns;
  c_klu session neqs nnz

external c_set_ordering : 'k serial_session -> ordering -> unit
  = "c_arkode_klu_set_ordering"

let set_ordering session ordering =
  ls_check_klu session;
  c_set_ordering session ordering

external c_reinit : 'k serial_session -> int -> int -> bool -> unit
  = "c_arkode_klu_reinit"

let reinit session n nnz realloc =
  ls_check_klu session;
  c_reinit session n nnz realloc

external c_get_num_jac_evals : 'k serial_session -> int
  = "c_arkode_klu_get_num_jac_evals"

let get_num_jac_evals session =
  ls_check_klu session;
  c_get_num_jac_evals session

module Mass = struct
  include SlsTypes.MassTypes

  external c_klu : 'k serial_session -> int -> int -> unit
    = "c_arkode_mass_klu_init"

  let klu f nnz session nv =
    let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
    session.mass_callbacks <- SlsKluMassCallback { massfn = f; smat = None };
    session.mass_precfns <- NoMassPrecFns;
    c_klu session neqs nnz

  external c_set_ordering : 'k serial_session -> ordering -> unit
    = "c_arkode_mass_klu_set_ordering"

  let set_ordering session ordering =
    mass_check_klu session;
    c_set_ordering session ordering

  external c_reinit : 'k serial_session -> int -> int -> bool -> unit
    = "c_arkode_mass_klu_reinit"

  let reinit session n nnz realloc =
    mass_check_klu session;
    c_reinit session n nnz realloc

  external c_get_num_evals : 'k serial_session -> int
    = "c_arkode_klu_get_num_mass_evals"

  let get_num_evals session =
    mass_check_klu session;
    c_get_num_jac_evals session
end

