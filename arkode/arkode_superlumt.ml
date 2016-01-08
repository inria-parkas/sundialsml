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

(* Must correspond with arkode_superlumt_ordering_tag *)
type ordering =
     Natural
   | MinDegreeProd
   | MinDegreeSum
   | ColAmd

external c_superlumt : serial_session -> int -> int -> int -> unit
  = "c_arkode_superlumt_init"

let superlumt f ~nnz ~nthreads session nv =
  let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
  session.ls_callbacks <- SlsSuperlumtCallback { jacfn = f; smat = None };
  session.ls_precfns <- NoPrecFns;
  c_superlumt session neqs nnz nthreads

external c_set_ordering : serial_session -> ordering -> unit
  = "c_arkode_superlumt_set_ordering"

let set_ordering session ordering =
  ls_check_superlumt session;
  c_set_ordering session ordering

external c_get_num_jac_evals : serial_session -> int
  = "c_arkode_superlumt_get_num_jac_evals"

let get_num_jac_evals session =
  ls_check_superlumt session;
  c_get_num_jac_evals session

module Mass = struct
  include SlsTypes.MassTypes

  external c_superlumt : serial_session -> int -> int -> int -> unit
    = "c_arkode_mass_superlumt_init"

  let superlumt f ~nnz ~nthreads session nv =
    let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
    session.mass_callbacks
      <- SlsSuperlumtMassCallback { massfn = f; smat = None };
    session.mass_precfns <- NoMassPrecFns;
    c_superlumt session neqs nnz nthreads

  external c_set_ordering : serial_session -> ordering -> unit
    = "c_arkode_mass_superlumt_set_ordering"

  let set_ordering session ordering =
    mass_check_superlumt session;
    c_set_ordering session ordering

  external c_get_num_evals : serial_session -> int
    = "c_arkode_superlumt_get_num_mass_evals"

  let get_num_evals session =
    mass_check_superlumt session;
    c_get_num_evals session
end
