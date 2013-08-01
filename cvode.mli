(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*           Timothy Bourke (INRIA) and Marc Pouzet (LIENS)            *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
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
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 *)

include module type of Sundials
  with type Roots.t = Sundials.Roots.t
  and type Roots.root_event = Sundials.Roots.root_event
  and type RootDirs.t = Sundials.RootDirs.t
  and type RootDirs.root_direction = Sundials.RootDirs.root_direction

(** {2 General} *)

(** {3 Solver initialisation} *)

(**
 Specify a linear multistep method.

 @cvode <node3#ss:ivp_sol> IVP Solution
 @cvode <node5#sss:cvodemalloc> CVodeCreate
 *)
type lmm =
  | Adams   (** Non-stiff systems; Adams-Moulton formulas *)
  | BDF     (** Stiff systems;     Backward Differentiation Formulas *)

(**
 Specify a solution method.

 @cvode <node3#ss:ivp_sol> IVP Solution
 @cvode <node5#sss:cvodemalloc> CVodeCreate
 *)
type iter =
  | Newton of linear_solver (** Newton iteration with a given linear solver *)
  | Functional              (** Functional iteration (non-stiff systems only) *)

(**
 Specify a linear solver.

 The Lapack solvers require that both Sundials and the Ocaml interface were
 built to link with a LAPACK library.

 The Banded Krylov solvers imply an additional call to
 {{:CVODE_DOC_ROOT(node5#sss:cvbandpre)} CVBandPrecInit}.

 @cvode <node5#sss:lin_solv_init> Linear Solver Specification
                                                 Functions
 *)
and linear_solver =
  | Dense                                   (** Direct with dense matrix,
                                                see {!Cvode_serial.Dls} and
                                                {!Cvode_nvector.Dls}.*)
  | LapackDense                             (** Direct with dense matrix,
                                                with Lapack,
                                                see {!Cvode_serial.Dls} and
                                                {!Cvode_nvector.Dls}.*)

  | Band of bandrange                       (** Direct with banded matrix,
                                                see {!Cvode_serial.Dls}
                                                and {!Cvode_nvector.Dls}.
                                             *)
  | LapackBand of bandrange                 (** Direct with banded matrix
                                                with Lapack,
                                                see {!Cvode_serial.Dls}
                                                and {!Cvode_nvector.Dls}.
                                             *)

  | Diag                                    (** Diagonal approximation
                                                of the Jacobian,
                                                see {!Cvode_serial.Diag}
                                                and {!Cvode_nvector.Diag}. *)

  | Spgmr of sprange                        (** Krylov Spils solver: SPGMR,
                                                see {!Cvode_serial.Spils}
                                                and {!Cvode_nvector.Spils}. *)
  | Spbcg of sprange                        (** Krylov Spils solver: SPBCG,
                                                see {!Cvode_serial.Spils}
                                                and {!Cvode_nvector.Spils}. *)
  | Sptfqmr of sprange                      (** Krylov Spils solver: SPFQMR,
                                                see {!Cvode_serial.Spils}
                                                and {!Cvode_nvector.Spils}. *)

  | BandedSpgmr of sprange * bandrange      (** Krylov Spils solver
                                                with banded matrix: SPGMR,
                                                see {!Cvode_serial.Spils},
                                                {!Cvode_serial.BandPrec},
                                                {!Cvode_nvector.Spils}, and
                                                {!Cvode_nvector.BandPrec}. *)
  | BandedSpbcg of sprange * bandrange      (** Krylov Spils solver
                                                with banded matrix: SPBCG,
                                                see {!Cvode_serial.Spils},
                                                {!Cvode_serial.BandPrec},
                                                {!Cvode_nvector.Spils}, and
                                                {!Cvode_nvector.BandPrec}. *)
  | BandedSptfqmr of sprange * bandrange    (** Krylov Spils solver
                                                with banded matrix: SPTFQMR,
                                                see {!Cvode_serial.Spils},
                                                {!Cvode_serial.BandPrec},
                                                {!Cvode_nvector.Spils}, and
                                                {!Cvode_nvector.BandPrec}. *)

(**
 @cvode <node5#sss:lin_solve_init> CVBand
 @cvode <node5#sss:cvbandpre> CVBandPrecInit
 *)
and bandrange = {
    mupper : int; (** upper half-bandwidth of the Jacobian approximation. *)
    mlower : int; (** lower half-bandwidth of the Jacobian approximation. *)
  }


(**
 Parameters for Krylov solvers.
 @cvode <node5#sss:lin_solv_init> CVSpgmr/CVSpbcg/CVSptfqrm
 *)
and sprange = {
    pretype : preconditioning_type;
    maxl: int
  }

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

(** {3 Solver results, statistics, and errors} *)

(**
 Possible values returned when a solver step function succeeds.
 Failures are indicated by exceptions.

 @cvode <node5#sss:cvode> CVode
 *)
type solver_result =
  | Continue            (** CV_SUCCESS *)
  | RootsFound          (** CV_ROOT_RETURN *)
  | StopTimeReached     (** CV_TSTOP_RETURN *)

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

(**
 Type of values passed to a registered error handler function.

 @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
 *)
type error_details = {
    error_code : int;
    module_name : string;
    function_name : string;
    error_message : string;
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
  with type Bandmatrix.t = Dls.Bandmatrix.t
  and  type Directbandmatrix.t = Dls.Directbandmatrix.t
  and  type Densematrix.t = Dls.Densematrix.t
  and  type Directdensematrix.t = Dls.Directdensematrix.t
