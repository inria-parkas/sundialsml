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

open Ida_impl
include AdjointTypes'.SlsTypes
let tosession = AdjointTypes.tosession

type ordering = Ida_klu.ordering =
     Amd
   | ColAmd
   | Natural

type sparse_jac_fn =
    NoSens of sparse_jac_fn_no_sens
  | WithSens of sparse_jac_fn_with_sens

external c_klub : 'k serial_session -> int -> int -> int -> bool -> unit
  = "c_idas_klub_init"

let klu f nnz bs nv nv' =
  let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
  let session = tosession bs in
  let parent, which = AdjointTypes.parent_and_which bs in
  let use_sens = match f with NoSens _ -> false | WithSens _ -> true in
  c_klub parent which neqs nnz use_sens;
  session.ls_precfns <- NoPrecFns;
  match f with
  | NoSens fns ->
      session.ls_callbacks <-
          BSlsKluCallback { jacfn = fns; smat = None }
  | WithSens fbs ->
      session.ls_callbacks <-
        BSlsKluCallbackSens { jacfn_sens = fbs; smat_sens = None }

let set_ordering bs = Ida_klu.set_ordering (tosession bs)
let reinit bs = Ida_klu.reinit (tosession bs)
let get_num_jac_evals bs = Ida_klu.get_num_jac_evals (tosession bs)

