(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_serial.h (and code in cvode_serial.c) must also be updated.
 *)

include Sundials
include Dls

type preconditioning_type =
  | PrecNone
  | PrecLeft
  | PrecRight
  | PrecBoth

type bandrange = { mupper : int; mlower : int }

type sprange = int

type linear_solver =
  | Dense
  | Band of bandrange
  | LapackDense
  | LapackBand of bandrange
  | Spgmr of sprange
  | Spbcg of sprange
  | Sptfqmr of sprange

type iter =
  | Newton of linear_solver
  | Functional

type solver_result =
  | Continue
  | RootsFound
  | StopTimeReached

type root_direction = RootDirs.root_direction

type error_details = {
    error_code : int;
    module_name : string;
    function_name : string;
    error_message : string;
  }

(* Solver exceptions *)
exception IllInput
exception TooClose
exception TooMuchWork
exception TooMuchAccuracy
exception ErrFailure
exception ConvergenceFailure
exception LinearInitFailure
exception LinearSetupFailure
exception LinearSolveFailure
exception ResFuncFailure
exception FirstResFuncFailure
exception RepeatedResFuncErr
exception UnrecoverableResFuncErr
exception RootFuncFailure

(* Initial condition calculator exceptions *)
exception NoRecovery
exception BadEwt

(* get_dky exceptions *)
exception BadK
exception BadT
exception BadDky

let no_roots = (0, (fun _ _ _ _ -> ()))

(* Throw inside the f callback if the derivatives cannot be calculated at
   the given time. *)
exception RecoverableFailure

type integrator_stats = {
    num_steps : int;
    num_res_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

exception StopTimeReached

exception ZeroDiagonalElement of int

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("ida_RecoverableFailure",      RecoverableFailure);

    ("ida_StopTimeReached",         StopTimeReached);
    ("ida_IllInput",                IllInput);
    ("ida_TooClose",                TooClose);
    ("ida_TooMuchWork",             TooMuchWork);
    ("ida_TooMuchAccuracy",         TooMuchAccuracy);
    ("ida_ErrFailure",              ErrFailure);
    ("ida_ConvergenceFailure",      ConvergenceFailure);
    ("ida_LinearInitFailure",       LinearInitFailure);
    ("ida_LinearSetupFailure",      LinearSetupFailure);
    ("ida_LinearSolveFailure",      LinearSolveFailure);
    ("ida_ResFuncFailure",          ResFuncFailure);
    ("ida_FirstResFuncFailure",     FirstResFuncFailure);
    ("ida_RepeatedResFuncErr",      RepeatedResFuncErr);
    ("ida_NoRecovery",              NoRecovery);
    ("ida_BadEwt",                  BadEwt);
    ("ida_UnrecoverableResFuncErr", UnrecoverableResFuncErr);
    ("ida_RootFuncFailure",         RootFuncFailure);

    ("ida_BadK",                    BadK);
    ("ida_BadT",                    BadT);
    ("ida_BadDky",                  BadDky);

    ("ida_ZeroDiagonalElement",     ZeroDiagonalElement 0);
  ]
