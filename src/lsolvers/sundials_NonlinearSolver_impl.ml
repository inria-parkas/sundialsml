(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2019 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

(* Hack to ensure that Sundials.c_init_module is executed so that the global
   exceptions are properly registered. *)
let e = Sundials.RecoverableFailure

(* Types underlying NonlinearSolver.  *)

(* This module defines the nonlinear_solver types which are manipulated
   (abstractly) via the NonlinearSolver module. They are declared here so
   that they can also be used (concretely) from the <solver> modules.
   To ensure that this type will be opaque outside of Sundials/ML, we simply
   do not install the sundials_NonlinearSolver_impl.cmi file. *)

exception VectorOpError
exception IncorrectUse
exception ExtFail
exception NonlinearSolverInUse

module Senswrapper = struct (* {{{ *)

  type ('d, 'k) t

  external data : ('d, 'k) t -> 'd array
    = "sunml_senswrapper_data"

end (* }}} *)

(* For O/Onls, it doesn't really matter what we pass as the mem argument since
   it is not possible to do anything with it.
   For O/Cnls, the mem argument is provided by the C stubs *)

type user = unit

(* Used as a phantom type tag. Define as a variant type to help GADT
   exhaustiveness checking.
   See: https://github.com/ocaml/ocaml/issues/7028#issuecomment-473051557 *)
type 'a integrator = Integrator of 'a

type ('d, 's) sysfn = 'd -> 'd -> 's -> unit

type ('d, 's) lsetupfn = bool -> 's -> bool

type ('d, 's) lsolvefn = 'd -> 's -> unit

(* Used in sundials_nonlinearsolver_ml.c: nlsolver_convtest_tag *)
type convtest =
  | Success
  | Continue
  | Recover

type ('d, 's) convtestfn = 'd -> 'd -> float -> 'd -> 's -> convtest

type ('d, 'k, 's) cptr

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_callbacks_index *)
type ('d, 'k, 's) callbacks = {
  mutable sysfn      : ('d, 's) sysfn;
  mutable lsetupfn   : ('d, 's) lsetupfn;
  mutable lsolvefn   : ('d, 's) lsolvefn;
  mutable convtestfn : ('d, 's) convtestfn;
}

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_type *)
type nonlinear_solver_type =
  | RootFind
  | FixedPoint

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_solver_tag *)
type ('d, 'k, _) solver =
  | NewtonSolver     : ('d, 'k, 's) solver
  | FixedPointSolver : int (* argument for backwards compatability in Arkode *)
                       -> ('d, 'k, 's) solver
  | CustomSolver     : (('d, 'k) Nvector.t, 's) ops
                       -> ('d, 'k, 's) solver
  | CustomSensSolver : (('d, 'k) Senswrapper.t, 'a integrator) ops
                       -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) solver

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_index *)
and ('d, 'k, 's) nonlinear_solver = {
  rawptr : ('d, 'k, 's) cptr;
  solver : ('d, 'k, 's) solver;
  callbacks : ('d, 'k, 's) callbacks;
  mutable attached : bool;
}

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_ops_index *)
and ('nv, 's) ops = {
  nls_type           : nonlinear_solver_type;
  init               : (unit -> unit) option;

  setup              : ('nv -> 's -> unit) option;
  solve              : 'nv -> 'nv -> 'nv -> float -> bool -> 's -> unit;

  set_sys_fn         : ('nv, 's) sysfn -> unit;
  set_lsetup_fn      : (('nv, 's) lsetupfn -> unit) option;
  set_lsolve_fn      : (('nv, 's) lsolvefn -> unit) option;
  set_convtest_fn    : (('nv, 's) convtestfn -> unit) option;
  set_max_iters      : (int -> unit) option;
  set_info_file      : (Logfile.t -> unit) option;
  set_print_level    : (int -> unit) option;

  get_num_iters      : (unit -> int) option;
  get_cur_iter       : (unit -> int) option;
  get_num_conv_fails : (unit -> int) option;
}

(* Used to reference a nonlinear solver to prevent garbage collection,
   while hiding the fact that the data argument is a senswrapper. *)
type ('d, 'k, 's) nonlinear_solver_hold =
  | NoNLS
  | NLS of ('d, 'k, 's) nonlinear_solver
  | NLS_sens of (('d, 'k) Senswrapper.t, 'k, 's) nonlinear_solver

let attach ({ attached } as s) =
  if attached then raise NonlinearSolverInUse;
  s.attached <- true

let detach s =
  s.attached <- false

let assert_senswrapper_solver { solver } =
  match solver with
  | CustomSolver _ -> raise IncorrectUse
  | CustomSensSolver _ -> ()
  (* The native nl_solvers can handle senswrappers *)
  | FixedPointSolver _ -> ()
  | NewtonSolver     -> ()

let get_type (type d k s) ({ rawptr; solver } : (d, k, s) nonlinear_solver) =
  match solver with
  | FixedPointSolver _            -> FixedPoint
  | NewtonSolver                  -> RootFind
  | CustomSolver { nls_type }     -> nls_type
  | CustomSensSolver { nls_type } -> nls_type

let empty_sysfn _ _ _ =
  Sundials_impl.crash "Internal error: nls_sysn called\n"
let empty_lsetupfn _ _ =
  Sundials_impl.crash "Internal error: nls_lsetupfn called\n"
let empty_lsolvefn _ _ =
  Sundials_impl.crash "Internal error: nls_lsolvefn called\n"
let empty_convtestfn _ _ _ _ =
  Sundials_impl.crash "Internal error: nls_convtestfn called\n"

let empty_callbacks () = {
    sysfn      = empty_sysfn;
    lsetupfn   = empty_lsetupfn;
    lsolvefn   = empty_lsolvefn;
    convtestfn = empty_convtestfn;
  }

