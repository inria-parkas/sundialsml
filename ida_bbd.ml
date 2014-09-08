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

open Ida_impl

type data = Nvector_parallel.data
type parallel_session = (data, Nvector_parallel.kind) Ida.session
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
    local_fn : float -> data -> data -> data -> unit;
    comm_fn  : (float -> data -> data -> unit) option;
  }

let bbd_callbacks { local_fn; comm_fn } =
  { Bbd.local_fn = local_fn; Bbd.comm_fn = comm_fn }

let call_bbdlocal session t y y' gval =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBDCallback { Bbd.local_fn = f } ->
      adjust_retcode session true (f t y y') gval
  | _ -> assert false

let call_bbdcomm session t y y' =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBDCallback { Bbd.comm_fn = Some f } ->
    adjust_retcode session true (f t y) y'
  | _ -> assert false

external c_bbd_prec_init
    : parallel_session -> int -> bandwidths -> float -> bool -> unit
    = "c_ida_bbd_prec_init"

external c_spils_spgmr
  : ('a, 'k) session -> int -> unit
  = "c_ida_spils_spgmr"

external c_spils_spbcg
  : ('a, 'k) session -> int -> unit
  = "c_ida_spils_spbcg"

external c_spils_sptfqmr
  : ('a, 'k) session -> int -> unit
  = "c_ida_spils_sptfqmr"

let spgmr maxl bws dqrely cb session nv nv' =
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let ba, _, _ = Sundials.unvec nv in
  let localn   = Sundials.RealArray.length ba in
  c_spils_spgmr session maxl;
  c_bbd_prec_init session localn bws dqrely (cb.comm_fn <> None);
  session.ls_callbacks <- BBDCallback (bbd_callbacks cb)

let spbcg maxl bws dqrely cb session nv nv' =
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let ba, _, _ = Sundials.unvec nv in
  let localn   = Sundials.RealArray.length ba in
  c_spils_spbcg session maxl;
  c_bbd_prec_init session localn bws dqrely (cb.comm_fn <> None);
  session.ls_callbacks <- BBDCallback (bbd_callbacks cb)

let sptfqmr maxl bws dqrely cb session nv nv' =
  let maxl   = match maxl with None -> 0 | Some ml -> ml in
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  let ba, _, _ = Sundials.unvec nv in
  let localn   = Sundials.RealArray.length ba in
  c_spils_sptfqmr session maxl;
  c_bbd_prec_init session localn bws dqrely (cb.comm_fn <> None);
  session.ls_callbacks <- BBDCallback (bbd_callbacks cb)

external c_bbd_prec_reinit
    : parallel_session -> int -> int -> float -> unit
    = "c_ida_bbd_prec_reinit"

let reinit s mudq mldq dqrely =
  let dqrely = match dqrely with None -> 0.0 | Some v -> v in
  c_bbd_prec_reinit s mudq mldq dqrely

external get_work_space : parallel_session -> int * int
    = "c_ida_bbd_get_work_space"

external get_num_gfn_evals : parallel_session -> int
    = "c_ida_bbd_get_num_gfn_evals"


(* Let C code know about some of the values in this module.  *)
type fcn = Fcn : 'a -> fcn
external c_init_module : fcn array -> unit =
  "c_ida_bbd_init_module"

let _ =
  c_init_module
    (* Functions must be listed in the same order as
       callback_index in ida_bbd_ml.c.  *)
    [|Fcn call_bbdlocal;
      Fcn call_bbdcomm;
    |]
