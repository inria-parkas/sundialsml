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

(* Used as a phantom type tag. Define as a variant type to help GADT
   exhaustiveness checking.
   See: https://github.com/ocaml/ocaml/issues/7028#issuecomment-473051557 *)
type 'a integrator = Integrator of 'a

type ('nv, 's) sysfn = 'nv -> 'nv -> 's -> unit

type 's lsetupfn = bool -> 's -> bool

type ('nv, 's) lsolvefn = 'nv -> 's -> unit

(* Used in sundials_nonlinearsolver_ml.c: nlsolver_convtest_tag *)
type convtest =
  | Success
  | Continue
  | Recover

type ('nv, 's) convtestfn' = 'nv -> 'nv -> float -> 'nv -> 's -> convtest

type ('d, 'k, 's, 'v) cptr

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_callbacks_index *)
type ('nv, 's) callbacks = {
  mutable sysfn      : ('nv, 's) sysfn;
  mutable lsetupfn   : 's lsetupfn;
  mutable lsolvefn   : ('nv, 's) lsolvefn;
  mutable convtestfn : ('nv, 's) convtestfn';
}

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_type *)
type nonlinear_solver_type =
  | RootFind
  | FixedPoint

type ('d, 'k, 's, _) solver =
  (* C-NLS used from C or OCaml *)
  | NewtonSolver :
      ('d, 's) callbacks (* C-NLS to C callback: only used for convtestfn
                            C-NLS to OCaml callback: backlink to 'd *)
      -> ('d, 'k, 's, [`Nvec]) solver

  (* C-NLS used from C or OCaml *)
  | NewtonSolverSens :
      (('d, 'k) Senswrapper.t, 's) callbacks
                         (* C-NLS to C callback: only used for convtestfn
                            C-NLS to OCaml callback: backlink to 'd *)
      -> ('d, 'k, 's, [`Sens]) solver

  (* C-NLS used from C or OCaml *)
  | FixedPointSolver :
      ('d, 's) callbacks (* C-NLS to C callback: only used for convtestfn
                            C-NLS to OCaml callback: backlink to 'd *)
      * int (* argument for backwards compatability in Arkode *)
      -> ('d, 'k, 's, [`Nvec]) solver

  (* C-NLS used from C or OCaml *)
  | FixedPointSolverSens :
      (('d, 'k) Senswrapper.t, 's) callbacks
                         (* C-NLS to C callback: only used for convtestfn
                            C-NLS to OCaml callback: backlink to 'd *)
      * int (* argument for backwards compatability in Arkode *)
      -> ('d, 'k, 's, [`Sens]) solver

  (* OCaml-NLS used from C or OCaml *)
  | CustomSolver :
      (('d, 'k) Nvector.t, 's) callbacks
                         (* OCaml-NLS to C callback: needs N_Vector
                            OCaml-NLS to OCaml callback: would prefer 'd... *)
      * (('d, 'k) Nvector.t, 'd, 's, [`Nvec]) ops
      -> ('d, 'k, 's, [`Nvec]) solver

  (* OCaml-NLS used from C only *)
  (* A special effort is required for Senswrappers. They must be passed
     directly to setup and solve so that they can then be passed onto
     to sysfn, lsolvefn, and convtestfn. Passing the payload directly would
     not be useful since a Senswrapper cannot be constructed from OCaml. *)
  | CustomSolverSens :
      (('d, 'k) Senswrapper.t, 's) callbacks
                          (* OCaml-NLS to C callback: needs Senswrapper *)
      * (('d, 'k) Senswrapper.t, ('d, 'k) Senswrapper.t, 's, [`Sens]) ops
      -> ('d, 'k, 's, [`Sens]) solver

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_index *)
and ('d, 'k, 's, 'v) nonlinear_solver = {
  rawptr : ('d, 'k, 's, 'v) cptr;
  solver : ('d, 'k, 's, 'v) solver;
  context : Sundials.Context.t;
  mutable info_file : Logfile.t option;
  mutable attached : bool;
}

and 's convtest_callback =
  { f : 'd1 'k1 't2 'd2 'k2. ('d1, 'k1, 't2, [`Nvec]) nonlinear_solver
      -> (('d2, 'k2) Nvector.t, 's) convtestfn' }
  [@@unboxed]

and 's convtest_callback_sens =
  { f : 'd1 'k1 't2 'd2 'k2. ('d1, 'k1, 't2, [`Sens]) nonlinear_solver
      -> (('d2, 'k2) Senswrapper.t, 's) convtestfn' }
  [@@unboxed]

and ('nv, 's, 'v) convtestfn =
  | CConvTest
    : 's convtest_callback cfun -> ('nv, 's, [`Nvec]) convtestfn
  | CSensConvTest
    : 's convtest_callback_sens cfun -> ('nv, 's, [`Sens]) convtestfn
  | OConvTest of ('nv, 's) convtestfn'

(* Distinguish operations that take an nvector ('nv), since the C function
   (sysfn, lsolvefn, convtestfn) passed from the integrator requires
   one, from those that take the payload directly ('d), since they are invoked
   from C with nvectors. This asymmetry is an unfortunate consequence of the
   implementation of nvectors in Sundials/ML. It occurs since there are call
   chains from C to OCaml to C. *)

(* Accessed from sundials_nonlinearsolver_ml.c: nlsolver_ops_index *)
and ('nv, 'd, 's, 'v) ops = {
  nls_type           : nonlinear_solver_type;
  init               : (unit -> unit) option;

  setup              : ('d -> 's -> unit) option;
  solve              : 'd -> 'd -> 'd -> float -> bool -> 's -> unit;

  set_sys_fn         : ('nv, 's) sysfn -> unit;
  set_lsetup_fn      : ('s lsetupfn -> unit) option;
  set_lsolve_fn      : (('nv, 's) lsolvefn -> unit) option;
  set_convtest_fn    : (('d, 's, 'v) convtestfn -> unit) option;
  set_max_iters      : (int -> unit) option;

  get_num_iters      : (unit -> int) option;
  get_cur_iter       : (unit -> int) option;
  get_num_conv_fails : (unit -> int) option;
}

(* Used to reference a nonlinear solver to prevent garbage collection,
   while hiding the fact that the data argument is a senswrapper. *)
type ('d, 'k, 's) nonlinear_solver_hold =
  | NoNLS
  | NLS of ('d, 'k, 's, [`Nvec]) nonlinear_solver
  | NLS_sens of ('d, 'k, 's, [`Sens]) nonlinear_solver

let attach ({ attached } as s) =
  if attached then raise NonlinearSolverInUse;
  s.attached <- true

let detach s =
  s.attached <- false

let get_type (type d k s v) ({ solver; _ } : (d, k, s, v) nonlinear_solver) =
  match solver with
  | FixedPointSolver _ -> FixedPoint
  | FixedPointSolverSens _ -> FixedPoint
  | NewtonSolver _ -> RootFind
  | NewtonSolverSens _ -> RootFind
  | CustomSolver (_, { nls_type }) -> nls_type
  | CustomSolverSens (_, { nls_type }) -> nls_type

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

