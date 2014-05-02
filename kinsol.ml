(***********************************************************************)
(*                                                                     *)
(*               OCaml interface to (serial) Sundials                  *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

exception IllInput                       (* KIN_ILL_INPUT *)
exception LineSearchNonConvergence       (* KIN_LINESEARCH_NONCONV *)
exception MaxIterationsReached           (* KIN_MAXITER_REACHED *)
exception MaxNewtonStepExceeded          (* KIN_MXNEWT_5X_EXCEEDED *)
exception LineSearchBetaConditionFailure (* KIN_LINESEARCH_BCFAIL *)
exception LinearSolverNoRecovery         (* KIN_LINSOLV_NO_RECOVERY *)
exception LinearSolverInitFailure        (* KIN_LINIT_FAIL *)
exception LinearSetupFailure             (* KIN_LSETUP_FAIL *)
exception LinearSolverFailure            (* KIN_LSOLVE_FAIL *)
exception SystemFunctionFailure          (* KIN_SYSFUNC_FAIL *)
exception FirstSystemFunctionFailure     (* KIN_FIRST_SYSFUNC_FAIL *)
exception RepeatedSystemFunctionFailure  (* KIN_REPTD_SYSFUNC_ERR *)

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("kinsol_RecoverableFailure",             Sundials.RecoverableFailure);

    ("kinsol_IllInput",                       IllInput);
    ("kinsol_LineSearchNonConvergence",       LineSearchNonConvergence);
    ("kinsol_MaxIterationsReached",           MaxIterationsReached);
    ("kinsol_MaxNewtonStepExceeded",          MaxNewtonStepExceeded);
    ("kinsol_LineSearchBetaConditionFailure", LineSearchBetaConditionFailure);
    ("kinsol_LinearSolverNoRecovery",         LinearSolverNoRecovery);
    ("kinsol_LinearSolverInitFailure",        LinearSolverInitFailure);
    ("kinsol_LinearSetupFailure",             LinearSetupFailure);
    ("kinsol_LinearSolverFailure",            LinearSolverFailure);
    ("kinsol_SystemFunctionFailure",          SystemFunctionFailure);
    ("kinsol_FirstSystemFunctionFailure",     FirstSystemFunctionFailure);
    ("kinsol_RepeatedSystemFunctionFailure",  RepeatedSystemFunctionFailure);
  ]

