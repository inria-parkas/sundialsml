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
