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

open Kinsol_impl

type data = Nvector_parallel.data
type parallel_session = (data, Nvector_parallel.kind) session
type parallel_linear_solver = (data, Nvector_parallel.kind) linear_solver

type bandwidths =
  {
    mudq    : int;
    mldq    : int;
    mukeep  : int;
    mlkeep  : int;
  }

type callbacks =
  {
    local_fn : data -> data -> unit;
    comm_fn  : (data -> unit) option;
  }

let bbd_callbacks { local_fn; comm_fn } =
  { Bbd.local_fn = local_fn; Bbd.comm_fn = comm_fn }

let call_bbdlocal session u gval =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBDCallback { Bbd.local_fn = f } ->
      adjust_retcode session true (f u) gval
  | _ -> assert false

let call_bbdcomm session u =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBDCallback { Bbd.comm_fn = Some f } -> adjust_retcode session true f u
  | _ -> assert false

let _ =
  Callback.register "c_kinsol_call_bbdlocal"      call_bbdlocal;
  Callback.register "c_kinsol_call_bbdcomm"       call_bbdcomm

external c_bbd_prec_init
    : parallel_session -> int -> bandwidths -> float -> bool -> unit
    = "c_kinsol_bbd_prec_init"

external c_set_max_restarts : ('a, 'k) session -> int -> unit
    = "c_kinsol_spils_set_max_restarts"

external c_spils_spgmr : ('a, 'k) session -> int -> unit
  = "c_kinsol_spils_spgmr"

external c_spils_spbcg : ('a, 'k) session -> int -> unit
  = "c_kinsol_spils_spbcg"

external c_spils_sptfqmr : ('a, 'k) session -> int -> unit
  = "c_kinsol_spils_sptfqmr"

let spgmr maxl omaxrs bws dqrely cb session onv =
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let localn =
    match onv with
      None -> 0
    | Some nv -> let ba, _, _ = Sundials.unvec nv in
                 Sundials.RealArray.length ba
  in
  c_spils_spgmr session maxl;
  (match omaxrs with
   | None -> ()
   | Some maxrs -> c_set_max_restarts session maxrs);
  c_bbd_prec_init session localn bws dqrely (cb.comm_fn <> None);
  session.ls_callbacks <- BBDCallback (bbd_callbacks cb)

let spbcg maxl bws dqrely cb session onv =
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let localn =
    match onv with
      None -> 0
    | Some nv -> let ba, _, _ = Sundials.unvec nv in
                 Sundials.RealArray.length ba
  in
  c_spils_spbcg session maxl;
  c_bbd_prec_init session localn bws dqrely (cb.comm_fn <> None);
  session.ls_callbacks <- BBDCallback (bbd_callbacks cb)

let sptfqmr maxl bws dqrely cb session onv =
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let localn =
    match onv with
      None -> 0
    | Some nv -> let ba, _, _ = Sundials.unvec nv in
                 Sundials.RealArray.length ba
  in
  c_spils_sptfqmr session maxl;
  c_bbd_prec_init session localn bws dqrely (cb.comm_fn <> None);
  session.ls_callbacks <- BBDCallback (bbd_callbacks cb)

external get_work_space : parallel_session -> int * int
    = "c_kinsol_bbd_get_work_space"

external get_num_gfn_evals : parallel_session -> int
    = "c_kinsol_bbd_get_num_gfn_evals"

