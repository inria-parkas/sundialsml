(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

(* Hack to ensure that Sundials.c_init_module is executed so that the global
   exceptions are properly registered. *)
let e = Sundials.RecoverableFailure

(* Types shared between Kinsol and Kinsol_bbd.  See the notes on
   Cvode_impl about the rationale behind this module.  *)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in kinsol_ml.h (and code in kinsol_ml.c) must also be updated.
 *)

module LSI = Sundials_LinearSolver_impl

(* Dummy callbacks.  These dummes getting called indicates a fatal
   bug.  Rather than raise an exception (which may or may not get
   propagated properly depending on the context), we immediately abort
   the program. *)

type ('a, 'k) nvector = ('a, 'k) Nvector.t

type 'a double = 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_u   : 'a;
    jac_fu  : 'a;
    jac_tmp : 't
  }

type real_array = RealArray.t

type bandrange = { mupper : int; mlower : int; }

module DirectTypes = struct
  type 'm jac_fn =
    (real_array double, real_array) jacobian_arg -> 'm -> unit

  (* These fields are accessed from kinsol_ml.c *)
  type 'm jac_callback =
    {
      jacfn : 'm jac_fn;
      mutable jmat : 'm option
    }

  let no_callback = fun _ _ -> Sundials_impl.crash "no direct callback"
end

module SpilsTypes' = struct
  type 'a solve_arg = { uscale : 'a; fscale : 'a; }

  type 'a prec_solve_fn =
    (unit, 'a) jacobian_arg
    -> 'a solve_arg
    -> 'a
    -> unit

  type 'a prec_setup_fn =
    (unit, 'a) jacobian_arg
    -> 'a solve_arg
    -> unit

  type 'a jac_times_vec_fn = 'a -> 'a -> 'a -> bool -> bool

  type 'a precfns =
    {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }
end

module KinsolBbdParamTypes = struct
  type 'a local_fn = 'a -> 'a -> unit
  type 'a comm_fn = 'a -> unit
  type 'a precfns =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end

module KinsolBbdTypes = struct
  type bandwidths =
    {
      mudq    : int;
      mldq    : int;
      mukeep  : int;
      mlkeep  : int;
    }
end

type kin_mem
type c_weak_ref

type 'a sysfn = 'a -> 'a -> unit
type errh = Util.error_details -> unit
type infoh = Util.error_details -> unit

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.

   The ls_solver field only exists to ensure that the linear solver is not
   garbage collected while still being used by a session.
*)

type ('a, 'k) session = {
  kinsol    : kin_mem;
  backref   : c_weak_ref;
  initvec   : ('a, 'k) Nvector.t;   (* for the set_linear_solver call. *)
  checkvec  : (('a, 'k) Nvector.t -> unit);
  context   : Context.t;

  mutable neqs       : int;    (* only valid for 'kind = serial *)
  mutable exn_temp   : exn option;

  mutable sysfn      : 'a sysfn;
  mutable errh       : errh;
  mutable infoh      : infoh;

  mutable error_file : Logfile.t option;

  mutable ls_solver    : LSI.held_linear_solver;
  mutable ls_callbacks : ('a, 'k) linsolv_callbacks;
  mutable ls_precfns : 'a linsolv_precfns;
}

and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  | DlsDenseCallback
      of Matrix.Dense.t DirectTypes.jac_callback
  | DlsBandCallback
      of Matrix.Band.t  DirectTypes.jac_callback

  | SlsKluCallback
      : ('s Matrix.Sparse.t) DirectTypes.jac_callback
        -> ('a, 'kind) linsolv_callbacks
  | SlsSuperlumtCallback
      : ('s Matrix.Sparse.t) DirectTypes.jac_callback
        -> ('a, 'kind) linsolv_callbacks

  (* Custom *)
  | DirectCustomCallback :
      'm DirectTypes.jac_callback -> ('a, 'kind) linsolv_callbacks

  | SpilsCallback1 of 'a SpilsTypes'.jac_times_vec_fn option
  | SpilsCallback2 of 'a sysfn

and 'a linsolv_precfns =
  | NoPrecFns
  | PrecFns of 'a SpilsTypes'.precfns
  | BBDPrecFns of 'a KinsolBbdParamTypes.precfns

(* Linear solver check functions *)

let ls_check_direct session =
  if Sundials_configuration.safe then
    match session.ls_callbacks with
    | DlsDenseCallback _ | DlsBandCallback _
    | SlsKluCallback _
    | SlsSuperlumtCallback _
    | DirectCustomCallback _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils session =
  if Sundials_configuration.safe then
    match session.ls_callbacks with
    | SpilsCallback1 _ | SpilsCallback2 _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils_bbd session =
  if Sundials_configuration.safe then
    match session.ls_precfns with
    | BBDPrecFns _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

(* Types that depend on session *)
type 'kind serial_session = (Nvector_serial.data, 'kind) session
                            constraint 'kind = [>Nvector_serial.kind]

type ('data, 'kind) linear_solver =
  ('data, 'kind) session
  -> ('data, 'kind) nvector
  -> unit

module SpilsTypes = struct
  include SpilsTypes'

  type ('a, 'k) set_preconditioner =
    ('a, 'k) session
    -> ('a, 'k) nvector
    -> unit

  type ('a, 'k) preconditioner =
    LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner

end

let read_weak_ref x : ('a, 'k) session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

let dummy_sysfn _ _ =
  Sundials_impl.crash "Internal error: dummy_sysfn called\n"
let dummy_errh _ =
  Sundials_impl.crash "Internal error: dummy_errh called\n"
let dummy_infoh _ =
  Sundials_impl.crash "Internal error: dummy_infoh called\n"
