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
include AdjointTypes'.SlsTypes
let tosession = AdjointTypes.tosession

type ordering = Cvode.Sls.Klu.ordering =
     Amd
   | ColAmd
   | Natural

type sparse_jac_fn =
    NoSens of sparse_jac_fn_no_sens
  | WithSens of sparse_jac_fn_with_sens

external c_klub : 'k serial_session -> int -> int -> int -> bool -> unit
  = "c_cvodes_klub_init"

let klu f nnz bs nv =
  let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
  let session = tosession bs in
  let parent, which = AdjointTypes.parent_and_which bs in
  let use_sens = match f with WithSens _ -> true | NoSens _ -> false in
  c_klub parent which neqs nnz use_sens;
  session.ls_precfns <- NoPrecFns;
  match f with
  | NoSens fns ->
      session.ls_callbacks <-
          BSlsKluCallback { jacfn = fns; smat = None }
  | WithSens fbs ->
      session.ls_callbacks <-
          BSlsKluCallbackSens { jacfn_sens = fbs; smat_sens = None }

let set_ordering bs = Cvode.Sls.Klu.set_ordering (tosession bs)
let reinit bs = Cvode.Sls.Klu.reinit (tosession bs)
let get_num_jac_evals bs = Cvode.Sls.Klu.get_num_jac_evals (tosession bs)

