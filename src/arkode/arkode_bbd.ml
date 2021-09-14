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
open Sundials
include Arkode_impl
include ArkodeBbdTypes

(* These types can't be defined in Arkode_impl because they introduce
   dependence on Mpi.  Some duplication is unavoidable.  *)
type data = Nvector_parallel.data
type kind = Nvector_parallel.kind

type 'step parallel_session =
  (Nvector_parallel.data, Nvector_parallel.kind, 'step) Arkode_impl.session
  constraint 'step = [<Arkode.arkstep|Arkode.mristep]

type 'step parallel_preconditioner =
  (Nvector_parallel.data, Nvector_parallel.kind, 'step)
    Arkode_impl.SpilsTypes.preconditioner
  constraint 'step = [<Arkode.arkstep|Arkode.mristep]

module Impl = ArkodeBbdParamTypes
type local_fn = data Impl.local_fn
type comm_fn = data Impl.comm_fn
type precfns =
  {
    local_fn : local_fn;
    comm_fn : comm_fn option;
  }

let bbd_precfns { local_fn; comm_fn } =
  { Impl.local_fn = local_fn; Impl.comm_fn = comm_fn }

external c_bbd_prec_init
    : 'step parallel_session -> int -> bandwidths -> float -> bool -> unit
    = "sunml_arkode_bbd_prec_init"

let init_preconditioner dqrely bandwidths precfns session nv =
  let ba, _, _ = Nvector.unwrap nv in
  let localn   = RealArray.length ba in
  c_bbd_prec_init session localn bandwidths dqrely (precfns.comm_fn <> None);
  session.ls_precfns <- BBDPrecFns (bbd_precfns precfns)

let prec_left ?(dqrely=0.0) bandwidths ?comm local_fn =
  LSI.Iterative.(PrecLeft,
    init_preconditioner dqrely bandwidths { local_fn ; comm_fn = comm })

let prec_right ?(dqrely=0.0) bandwidths ?comm local_fn =
  LSI.Iterative.(PrecRight,
    init_preconditioner dqrely bandwidths { local_fn ; comm_fn = comm })

let prec_both ?(dqrely=0.0) bandwidths ?comm local_fn =
  LSI.Iterative.(PrecBoth,
    init_preconditioner dqrely bandwidths { local_fn ; comm_fn = comm })

external c_bbd_prec_reinit
    : 'step parallel_session -> int -> int -> float -> unit
    = "sunml_arkode_bbd_prec_reinit"

let reinit s ?(dqrely=0.0) mudq mldq =
  ls_check_spils_bbd s;
  match s.ls_precfns with
  | BBDPrecFns _ -> c_bbd_prec_reinit s mudq mldq dqrely
  | _ -> raise LinearSolver.InvalidLinearSolver

external get_work_space : 'step parallel_session -> int * int
    = "sunml_arkode_bbd_get_work_space"

let get_work_space s =
  ls_check_spils_bbd s;
  get_work_space s

external get_num_gfn_evals : 'step parallel_session -> int
    = "sunml_arkode_bbd_get_num_gfn_evals"

let get_num_gfn_evals s =
  ls_check_spils_bbd s;
  get_num_gfn_evals s

