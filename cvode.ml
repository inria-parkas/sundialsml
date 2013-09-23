(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_serial.h (and code in cvode_serial.c) must also be updated.
 *)

include Sundials
include Dls

type lmm =
  | Adams
  | BDF

type preconditioning_type =
  | PrecNone
  | PrecLeft
  | PrecRight
  | PrecBoth

type bandrange = { mupper : int; mlower : int }

type sprange = { pretype : preconditioning_type; maxl: int }

type linear_solver =
  | Dense
  | LapackDense
  | Band of bandrange
  | LapackBand of bandrange
  | Diag
  | Spgmr of sprange
  | Spbcg of sprange
  | Sptfqmr of sprange
  | BandedSpgmr of sprange * bandrange
  | BandedSpbcg of sprange * bandrange
  | BandedSptfqmr of sprange * bandrange

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
exception RhsFuncFailure
exception FirstRhsFuncErr
exception RepeatedRhsFuncErr
exception UnrecoverableRhsFuncErr
exception RootFuncFailure

(* get_dky exceptions *)
exception BadK
exception BadT
exception BadDky

let no_roots = (0, (fun _ _ _ -> ()))

(* Throw inside the f callback if the derivatives cannot be calculated at
   the given time. *)
exception RecoverableFailure

type integrator_stats = {
    num_steps : int;
    num_rhs_evals : int;
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
    ("cvode_RecoverableFailure",      RecoverableFailure);

    ("cvode_StopTimeReached",         StopTimeReached);
    ("cvode_IllInput",                IllInput);
    ("cvode_TooClose",                TooClose);
    ("cvode_TooMuchWork",             TooMuchWork);
    ("cvode_TooMuchAccuracy",         TooMuchAccuracy);
    ("cvode_ErrFailure",              ErrFailure);
    ("cvode_ConvergenceFailure",      ConvergenceFailure);
    ("cvode_LinearInitFailure",       LinearInitFailure);
    ("cvode_LinearSetupFailure",      LinearSetupFailure);
    ("cvode_LinearSolveFailure",      LinearSolveFailure);
    ("cvode_RhsFuncFailure",          RhsFuncFailure);
    ("cvode_FirstRhsFuncErr",         FirstRhsFuncErr);
    ("cvode_RepeatedRhsFuncErr",      RepeatedRhsFuncErr);
    ("cvode_UnrecoverableRhsFuncErr", UnrecoverableRhsFuncErr);
    ("cvode_RootFuncFailure",         RootFuncFailure);

    ("cvode_BadK",                    BadK);
    ("cvode_BadT",                    BadT);
    ("cvode_BadDky",                  BadDky);

    ("cvode_ZeroDiagonalElement",     ZeroDiagonalElement 0);
  ]


