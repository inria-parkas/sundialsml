(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)
open Sundials
include Ida_impl
include IdasBbdTypes

type parallel_session =
        (Nvector_parallel.data, Nvector_parallel.kind) Idas.session
type parallel_bsession =
        (Nvector_parallel.data, Nvector_parallel.kind) Idas.Adjoint.bsession
type parallel_preconditioner =
  (Nvector_parallel.data, Nvector_parallel.kind)
                AdjointTypes.SpilsTypes.preconditioner

let tosession = AdjointTypes.tosession

module Impl = IdasBbdParamTypes
type local_fn = Nvector_parallel.data Impl.local_fn
type comm_fn = Nvector_parallel.data Impl.comm_fn
type precfns =
  {
    local_fn : local_fn;
    comm_fn : comm_fn option;
  }

let bbd_precfns { local_fn; comm_fn } =
  { Impl.local_fn = local_fn; Impl.comm_fn = comm_fn }

external c_bbd_prec_initb
    : (parallel_session * int) -> int
      -> Ida_bbd.bandwidths -> float -> bool -> unit
    = "sunml_idas_bbd_prec_initb"

let parent_and_which s =
  match (tosession s).sensext with
  | BwdSensExt se -> (se.parent, se.which)
  | _ -> failwith "Internal error: bsession invalid"

let init_preconditioner dqrely bandwidths precfns bs parent which nv =
  let ba, _, _ = Nvector.unwrap nv in
  let localn   = RealArray.length ba in
  c_bbd_prec_initb (parent, which) localn bandwidths dqrely
    (precfns.comm_fn <> None);
  (tosession bs).ls_precfns <- BBBDPrecFns (bbd_precfns precfns)

let prec_left ?(dqrely=0.0) bandwidths ?comm local_fn =
  LSI.Iterative.(PrecLeft,
    init_preconditioner dqrely bandwidths { local_fn; comm_fn = comm })

external c_bbd_prec_reinitb
    : parallel_session -> int -> int -> int -> float -> unit
    = "sunml_idas_bbd_prec_reinitb"

let reinit bs ?(dqrely=0.0) mudq mldq =
  ls_check_spils_bbd (tosession bs);
  let parent, which = parent_and_which bs in
  c_bbd_prec_reinitb parent which mudq mldq dqrely

let get_work_space bs = Ida_bbd.get_work_space (tosession bs)
let get_num_gfn_evals bs = Ida_bbd.get_num_gfn_evals (tosession bs)

