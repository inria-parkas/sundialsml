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

include Ida_impl
include IdaBbdTypes

type parallel_session =
        (Nvector_parallel.data, Nvector_parallel.kind) Ida.session
type parallel_preconditioner =
        (Nvector_parallel.data, Nvector_parallel.kind) SpilsTypes.preconditioner

module Impl = IdaBbdParamTypes
type local_fn = Nvector_parallel.data Impl.local_fn
type comm_fn = Nvector_parallel.data Impl.comm_fn
type callbacks =
  {
    local_fn : local_fn;
    comm_fn : comm_fn option;
  }

let bbd_callbacks { local_fn; comm_fn } =
  { Impl.local_fn = local_fn; Impl.comm_fn = comm_fn }

let call_bbdlocal session t y y' gval =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBDCallback { Impl.local_fn = f } ->
      adjust_retcode session true (f t y y') gval
  | _ -> assert false

let call_bbdcomm session t y y' =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBDCallback { Impl.comm_fn = Some f } ->
    adjust_retcode session true (f t y) y'
  | _ -> assert false

external c_bbd_prec_init
    : parallel_session -> int -> bandwidths -> float -> bool -> unit
    = "c_ida_bbd_prec_init"

let init_preconditioner dqrely bandwidths callbacks session nv nv' =
  let ba, _, _ = Nvector.unwrap nv in
  let localn   = Sundials.RealArray.length ba in
  c_bbd_prec_init session localn bandwidths dqrely (callbacks.comm_fn <> None);
  session.ls_callbacks <- BBDCallback (bbd_callbacks callbacks)

let prec_left ?(dqrely=0.0) bandwidths ?comm_fn local_fn =
  SpilsTypes.InternalPrecLeft
    (init_preconditioner dqrely bandwidths { local_fn; comm_fn })

external c_bbd_prec_reinit
    : parallel_session -> int -> int -> float -> unit
    = "c_ida_bbd_prec_reinit"

let reinit s ?(dqrely=0.0) mudq mldq =
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
