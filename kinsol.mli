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

(***********************************************************************)
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Vector-independent types and values for the KINSOL solver.

 @version VERSION()
 @author Timothy Bourke (Inria)
 @author Jun Inoue (Inria)
 @author Marc Pouzet (LIENS)
 *)

(** {3:exceptions Exceptions} *)

(** @kinsol <node5#sss:kinsol> KIN_ILL_INPUT *)
exception IllInput

(** @kinsol <node5#sss:kinsol> KIN_LINESEARCH_NONCONV *)
exception LineSearchNonConvergence

(** @kinsol <node5#sss:kinsol> KIN_MAXITER_REACHED *)
exception MaxIterationsReached

(** @kinsol <node5#sss:kinsol> KIN_MXNEWT_5X_EXCEEDED *)
exception MaxNewtonStepExceeded

(** @kinsol <node5#sss:kinsol> KIN_LINESEARCH_BCFAIL *)
exception LineSearchBetaConditionFailure

(** @kinsol <node5#sss:kinsol> KIN_LINSOLV_NO_RECOVERY *)
exception LinearSolverNoRecovery

(** @kinsol <node5#sss:kinsol> KIN_LINIT_FAIL *)
exception LinearSolverInitFailure

(** @kinsol <node5#sss:kinsol> KIN_LSETUP_FAIL *)
exception LinearSetupFailure

(** @kinsol <node5#sss:kinsol> KIN_LSOLVE_FAIL *)
exception LinearSolverFailure

(** @kinsol <node5#sss:kinsol> KIN_SYSFUNC_FAIL *)
exception SystemFunctionFailure

(** @kinsol <node5#sss:kinsol> KIN_FIRST_SYSFUNC_FAIL *)
exception FirstSystemFunctionFailure

(** @kinsol <node5#sss:kinsol> KIN_REPTD_SYSFUNC_ERR *)
exception RepeatedSystemFunctionFailure

