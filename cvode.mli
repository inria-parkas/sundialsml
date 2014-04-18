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

(***********************************************************************)
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Vector-independent types and values for the CVODE solver.

 @version VERSION()
 @author Timothy Bourke (Inria)
 @author Jun Inoue (Inria)
 @author Marc Pouzet (LIENS)
 *)

include module type of Sundials
  with type Roots.t = Sundials.Roots.t
  and type Roots.root_event = Sundials.Roots.root_event
  and type RootDirs.t = Sundials.RootDirs.t
  and type RootDirs.root_direction = Sundials.RootDirs.root_direction
  and type solver_result = Sundials.solver_result
  and type error_details = Sundials.error_details

(** {2 General} *)

(**
 Specify a linear multistep method.

 @cvode <node3#ss:ivp_sol> IVP Solution
 @cvode <node5#sss:cvodemalloc> CVodeCreate
 *)
type lmm =
  | Adams   (** Non-stiff systems; Adams-Moulton formulas *)
  | BDF     (** Stiff systems;     Backward Differentiation Formulas *)

(**
 Type of preconditioning for Krylov solvers.
 @cvode <node3#s:preconditioning> Preconditioning
 @cvode <node5#sss:lin_solv_init> CVSpgmr/CVSpbcg/CVSptfqrm
 *)
and preconditioning_type =
  | PrecNone
  | PrecLeft
  | PrecRight
  | PrecBoth

(**
 Values for root directions.
 @cvode <node5#sss:optin_root> CVodeSetRootDirection
 *)
type root_direction = RootDirs.root_direction

(**
 This is a convenience value for signalling to {!Cvode_serial.init} and
 {!Cvode_nvector.init} that there are no roots (zero-crossings) to
 monitor.
 *)
val no_roots : (int * ('a -> 'b -> 'c -> unit))

(** {3 Integrator statistics} *)

(**
 Aggregated integrator statistics.
 @cvode <node5#sss:optout_main> CVodeGetIntegratorStats
 *)
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

(** {3:exceptions Exceptions} *)

(** @cvode <node5#sss:cvode> CV_ILL_INPUT *)
exception IllInput

(** @cvode <node5#sss:cvode> CV_TOO_CLOSE *)
exception TooClose

(** @cvode <node5#sss:cvode> CV_TOO_MUCH_WORK *)
exception TooMuchWork

(** @cvode <node5#sss:cvode> CV_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** @cvode <node5#sss:cvode> CV_ERR_FAIL *)
exception ErrFailure                

(** @cvode <node5#sss:cvode> CV_CONV_FAIL *)
exception ConvergenceFailure        

(** @cvode <node5#sss:cvode> CV_LINIT_FAIL *)
exception LinearInitFailure         

(** @cvode <node5#sss:cvode> CV_LSETUP_FAIL *)
exception LinearSetupFailure        

(** @cvode <node5#sss:cvode> CV_LSOLVE_FAIL *)
exception LinearSolveFailure        

(** @cvode <node5#sss:cvode> CV_RHSFUNC_FAIL *)
exception RhsFuncFailure

(** @cvode <node5#sss:cvode> CV_FIRST_RHSFUNC_ERR *)
exception FirstRhsFuncErr

(** @cvode <node5#sss:cvode> CV_REPTD_RHSFUNC_ERR *)
exception RepeatedRhsFuncErr        

(** @cvode <node5#sss:cvode> CV_UNREC_RHSFUNC_ERR *)
exception UnrecoverableRhsFuncErr   

(** @cvode <node5#sss:cvode> CV_RTFUNC_FAIL *)
exception RootFuncFailure           

exception BadK      (** k is not in the range 0, 1, ..., q_u (CV_BAD_K)
                        @cvode <node5#ss:optional_dky> CVodeGetDky *)

exception BadT      (** t is not in the interval
                        \[t_n - h_u, t_n\] (CV_BAD_T)
                        @cvode <node5#ss:optional_dky> CVodeGetDky *)
exception BadDky    (** invalid dky argument (CV_BAD_DKY)
                        @cvode <node5#ss:optional_dky> CVodeGetDky *)

(**
 This exception may be thrown inside the RHS callback function (f)
 to indicate that one or more derivatives cannot be calculated at
 the given time offset. *)
exception RecoverableFailure


(**
 Thrown by the getrf functions if a zero diagonal element is encountered during
 factorization. The argument indicates the column index (from 1).

 @cvode <node9#ss:dense> DenseGETRF/denseGETRF 
 *)
exception ZeroDiagonalElement of int

include module type of Dls
  with type BandMatrix.t = Dls.BandMatrix.t
  and  type ArrayBandMatrix.t = Dls.ArrayBandMatrix.t
  and  type DenseMatrix.t = Dls.DenseMatrix.t
  and  type ArrayDenseMatrix.t = Dls.ArrayDenseMatrix.t

