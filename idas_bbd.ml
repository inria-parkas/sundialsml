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
include IdasBbdTypes

type data = Nvector_parallel.data
type kind = Nvector_parallel.kind

type parallel_session = (data, Nvector_parallel.kind) Idas.session
type parallel_bsession = (data, Nvector_parallel.kind) Idas.Adjoint.bsession
type parallel_preconditioner =
  (data, kind) AdjointTypes.SpilsTypes.preconditioner

let tosession = AdjointTypes.tosession

module Impl = IdasBbdParamTypes
type local_fn = data Impl.local_fn
type comm_fn = data Impl.comm_fn
type callbacks =
  {
    local_fn : local_fn;
    comm_fn : comm_fn option;
  }

let call_bbbdlocal session t y y' yb y'b glocal =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBBDCallback { Impl.local_fn = f } ->
      adjust_retcode session true (f t y y' yb y'b) glocal
  | _ -> assert false

let call_bbbdcomm session t y y' yb y'b =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BBBDCallback { Impl.comm_fn = Some f } ->
      adjust_retcode session true (f t y y' yb) y'b
  | _ -> assert false

let bbd_callbacks { local_fn; comm_fn } =
  { Impl.local_fn = local_fn; Impl.comm_fn = comm_fn }

external c_bbd_prec_initb
    : (parallel_session * int) -> int
      -> Ida_bbd.bandwidths -> float -> bool -> unit
    = "c_idas_bbd_prec_initb"

let parent_and_which s =
  match (tosession s).sensext with
  | BwdSensExt se -> (se.parent, se.which)
  | _ -> failwith "Internal error: bsession invalid"

let init_preconditioner dqrely bandwidths callbacks bs parent which nv nv' =
  let ba, _, _ = Nvector.unwrap nv in
  let localn   = Sundials.RealArray.length ba in
  c_bbd_prec_initb (parent, which) localn bandwidths dqrely
    (callbacks.comm_fn <> None);
  (tosession bs).ls_callbacks <- BBBDCallback (bbd_callbacks callbacks)

let prec_left ?(dqrely=0.0) bandwidths callbacks =
  AdjointTypes.SpilsTypes.InternalPrecLeft
    (init_preconditioner dqrely bandwidths callbacks)

external c_bbd_prec_reinitb
    : parallel_session -> int -> int -> int -> float -> unit
    = "c_idas_bbd_prec_reinitb"

let reinit bs ?(dqrely=0.0) mudq mldq =
  let parent, which = parent_and_which bs in
  c_bbd_prec_reinitb parent which mudq mldq dqrely

let get_work_space bs = Ida_bbd.get_work_space (tosession bs)
let get_num_gfn_evals bs = Ida_bbd.get_num_gfn_evals (tosession bs)


(* Let C code know about some of the values in this module.  *)
type fcn = Fcn : 'a -> fcn
external c_init_module : fcn array -> unit =
  "c_idas_bbd_init_module"

let _ =
  c_init_module
    (* Functions must be listed in the same order as
       callback_index in idas_bbd_ml.c.  *)
    [|Fcn call_bbbdlocal;
      Fcn call_bbbdcomm;
    |]
