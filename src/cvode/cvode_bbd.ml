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

include Cvode_impl
include CvodeBbdTypes

(* These types can't be defined in Cvode_impl because they introduce
   dependence on Mpi.  Some duplication is unavoidable.  *)
type data = Nvector_parallel.data
type kind = Nvector_parallel.kind

type parallel_session = (data, kind) session
type parallel_preconditioner = (data, kind) Cvode.Iterative.preconditioner

module Impl = CvodeBbdParamTypes
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
    : parallel_session -> int -> bandwidths -> float -> bool -> unit
    = "c_cvode_bbd_prec_init"

let init_preconditioner dqrely bandwidths precfns session nv =
  let ba, _, _ = Nvector.unwrap nv in
  let localn   = Sundials.RealArray.length ba in
  c_bbd_prec_init session localn bandwidths dqrely (precfns.comm_fn <> None);
  session.ls_precfns <- BBDPrecFns (bbd_precfns precfns)

let prec_left ?(dqrely=0.0) bandwidths ?comm local_fn =
  Lsolver_impl.Iterative.(PrecLeft,
    init_preconditioner dqrely bandwidths { local_fn ; comm_fn = comm })

let prec_right ?(dqrely=0.0) bandwidths ?comm local_fn =
  Lsolver_impl.Iterative.(PrecRight,
    init_preconditioner dqrely bandwidths { local_fn ; comm_fn = comm })

let prec_both ?(dqrely=0.0) bandwidths ?comm local_fn =
  Lsolver_impl.Iterative.(PrecBoth,
    init_preconditioner dqrely bandwidths { local_fn ; comm_fn = comm })

external c_bbd_prec_reinit
    : parallel_session -> int -> int -> float -> unit
    = "c_cvode_bbd_prec_reinit"

let reinit s ?(dqrely=0.0) mudq mldq =
  ls_check_spils_bbd s;
  match s.ls_precfns with
  | BBDPrecFns _ -> c_bbd_prec_reinit s mudq mldq dqrely
  | _ -> raise Sundials.InvalidLinearSolver

external get_work_space : parallel_session -> int * int
    = "c_cvode_bbd_get_work_space"

let get_work_space s =
  ls_check_spils_bbd s;
  get_work_space s

external get_num_gfn_evals : parallel_session -> int
    = "c_cvode_bbd_get_num_gfn_evals"

let get_num_gfn_evals s =
  ls_check_spils_bbd s;
  get_num_gfn_evals s

