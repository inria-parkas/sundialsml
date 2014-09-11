(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

include Cvode_impl
include CvodesBbdTypes

(* These types can't be defined in Cvode_impl because they introduce
   dependence on Mpi.  Some duplication is unavoidable.  *)
type data = Nvector_parallel.data
type kind = Nvector_parallel.kind

type parallel_session = Cvode_bbd.parallel_session
type parallel_bsession = (data, kind) AdjointTypes.bsession
type parallel_linear_solver = (data, kind) AdjointTypes.linear_solver

let tosession = AdjointTypes.tosession

module Impl = CvodesBbdParamTypes
type local_fn = data Impl.local_fn
type comm_fn = data Impl.comm_fn
type callbacks =
  {
    local_fn : local_fn;
    comm_fn : comm_fn option;
  }

let bbd_callbacks { local_fn; comm_fn } =
  { Impl.local_fn = local_fn; Impl.comm_fn = comm_fn }

let call_bbbdlocal session t y yb glocal =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBBDCallback { Impl.local_fn = f } ->
      adjust_retcode session true (f t y yb) glocal
  | _ -> assert false

let call_bbbdcomm session t y yb =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBBDCallback { Impl.comm_fn = Some f } ->
      adjust_retcode session true (f t y) yb
  | _ -> assert false

external c_bbd_prec_initb
    : (parallel_session * int) -> int
      -> Cvode_bbd.bandwidths -> float -> bool -> unit
    = "c_cvodes_bbd_prec_initb"

let parent_and_which s =
  match (tosession s).sensext with
  | BwdSensExt se -> (se.parent, se.which)
  | _ -> failwith "Internal error: bsession invalid"

external c_spils_spgmr
  : ('a, 'k) session -> int -> int -> Spils.preconditioning_type -> unit
  = "c_cvodes_adj_spils_spgmr"

external c_spils_spbcg
  : ('a, 'k) session -> int -> int -> Spils.preconditioning_type -> unit
  = "c_cvodes_adj_spils_spbcg"

external c_spils_sptfqmr
  : ('a, 'k) session -> int -> int -> Spils.preconditioning_type -> unit
  = "c_cvodes_adj_spils_sptfqmr"

let spgmr maxl prec_type bws dqrely cb bs nv =
  let parent, which = parent_and_which bs in
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let ba, _, _ = Sundials.unvec nv in
  let localn   = Sundials.RealArray.length ba in
  c_spils_spgmr parent which maxl prec_type;
  c_bbd_prec_initb (parent, which) localn bws dqrely
                                            (cb.comm_fn <> None);
  (tosession bs).ls_callbacks <- BBBDCallback (bbd_callbacks cb)

let spbcg maxl prec_type bws dqrely cb bs nv =
  let parent, which = parent_and_which bs in
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let ba, _, _ = Sundials.unvec nv in
  let localn   = Sundials.RealArray.length ba in
  c_spils_spbcg parent which maxl prec_type;
  c_bbd_prec_initb (parent, which) localn bws dqrely
                                            (cb.comm_fn <> None);
  (tosession bs).ls_callbacks <- BBBDCallback (bbd_callbacks cb)

let sptfqmr maxl prec_type bws dqrely cb bs nv =
  let parent, which = parent_and_which bs in
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let ba, _, _ = Sundials.unvec nv in
  let localn   = Sundials.RealArray.length ba in
  c_spils_sptfqmr parent which maxl prec_type;
  c_bbd_prec_initb (parent, which) localn bws dqrely
                                            (cb.comm_fn <> None);
  (tosession bs).ls_callbacks <- BBBDCallback (bbd_callbacks cb)

external c_bbd_prec_reinitb
    : parallel_session -> int -> int -> int -> float -> unit
    = "c_cvodes_bbd_prec_reinitb"

let reinit bs mudq mldq dqrely =
  let parent, which = parent_and_which bs in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  c_bbd_prec_reinitb parent which mudq mldq dqrely

let get_work_space bs = Cvode_bbd.get_work_space (tosession bs)
let get_num_gfn_evals bs = Cvode_bbd.get_num_gfn_evals (tosession bs)


(* Let C code know about some of the values in this module.  *)
type fcn = Fcn : 'a -> fcn
external c_init_module : fcn array -> unit =
  "c_cvodes_bbd_init_module"

let _ =
  c_init_module
    (* Functions must be listed in the same order as
       callback_index in cvodes_bbd_ml.c.  *)
    [|Fcn call_bbbdlocal;
      Fcn call_bbbdcomm;
    |]
