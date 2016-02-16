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

type ordering = Cvode_superlumt.ordering =
     Natural
   | MinDegreeProd
   | MinDegreeSum
   | ColAmd

type sparse_jac_fn =
    NoSens of sparse_jac_fn_no_sens
  | WithSens of sparse_jac_fn_with_sens

external c_superlumtb : ('k serial_session * int)
                        -> int -> int -> int -> bool -> unit
  = "c_cvodes_superlumtb_init"

let superlumt f ~nnz ~nthreads bs nv =
  let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
  let session = tosession bs in
  let use_sens =
    match f with
    | NoSens fns ->
        (session.ls_callbacks <-
            BSlsSuperlumtCallback { jacfn = fns; smat = None };
         false)
    | WithSens fbs ->
        (session.ls_callbacks <-
            BSlsSuperlumtCallbackSens { jacfn_sens = fbs; smat_sens = None };
         true)
  in
  session.ls_precfns <- NoPrecFns;
  c_superlumtb (AdjointTypes.parent_and_which bs) neqs nnz nthreads use_sens

let set_ordering bs = Cvode_superlumt.set_ordering (tosession bs)
let get_num_jac_evals bs = Cvode_superlumt.get_num_jac_evals (tosession bs)

