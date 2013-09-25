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
(*               User Documentation for IDA v2.7.0                     *)
(*         Alan C. Hindmarsh, Radu Serban, and Aaron Collier           *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Vector-independent types and values for the IDA solver.

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

(** {2 General} *)

(** {3 Solver initialisation} *)

(**
 Specify a linear solver.

 The Lapack solvers require that both Sundials and the OCaml interface were
 built to link with a LAPACK library.

 IDA supports direct linear solvers and Krylov solvers, but for the latter
 does not support banded preconditioning.

 @ida <node5#sss:lin_solv_init> Linear Solver Specification
                                                 Functions
 *)
type linear_solver =
  | Dense                                   (** Direct with dense matrix,
                                                see {!Ida_serial.Dls} and
                                                {!Ida_nvector.Dls}.*)
  | Band of bandrange                       (** Direct with banded matrix,
                                                see {!Ida_serial.Dls}
                                                and {!Ida_nvector.Dls}.
                                             *)
  | LapackDense                             (** Direct with dense matrix,
                                                with Lapack,
                                                see {!Ida_serial.Dls} and
                                                {!Ida_nvector.Dls}.*)

  | LapackBand of bandrange                 (** Direct with banded matrix
                                                with Lapack,
                                                see {!Ida_serial.Dls}
                                                and {!Ida_nvector.Dls}.
                                             *)

  | Spgmr of sprange                        (** Krylov Spils solver: SPGMR,
                                                see {!Ida_serial.Spils}
                                                and {!Ida_nvector.Spils}. *)
  | Spbcg of sprange                        (** Krylov Spils solver: SPBCG,
                                                see {!Ida_serial.Spils}
                                                and {!Ida_nvector.Spils}. *)
  | Sptfqmr of sprange                      (** Krylov Spils solver: SPFQMR,
                                                see {!Ida_serial.Spils}
                                                and {!Ida_nvector.Spils}. *)

(**
 @ida <node5#sss:lin_solve_init> IDABand
 @ida <node5#sss:idabandpre> IDABandPrecInit
 *)
and bandrange = {
    mupper : int; (** upper half-bandwidth of the Jacobian approximation. *)
    mlower : int; (** lower half-bandwidth of the Jacobian approximation. *)
  }


(**
 Parameters for Krylov solvers.
 @ida <node5#sss:lin_solv_init> IDASpgmr/IDASpbcg/IDASptfqrm
 *)
and sprange = int

(**
 Values for root directions.
 @ida <node5#sss:optin_root> IDASetRootDirection
 *)
type root_direction = RootDirs.root_direction

(**
 This is a convenience value for signalling to {!Ida_serial.init} and
 {!Ida_nvector.init} that there are no roots (zero-crossings) to
 monitor.
 *)
val no_roots : (int * ('a -> 'b -> 'c -> 'd -> unit))

(** {3 Solver results, statistics, and errors} *)

(**
 Possible values returned when a solver step function succeeds.
 Failures are indicated by exceptions.

 @ida <node5#sss:ida> IDA
 *)
type solver_result =
  | Continue            (** IDA_SUCCESS *)
  | RootsFound          (** IDA_ROOT_RETURN *)
  | StopTimeReached     (** IDA_TSTOP_RETURN *)

(**
 Aggregated integrator statistics.
 @ida <node5#sss:optout_main> IDAGetIntegratorStats
 *)
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

(**
 Type of values passed to a registered error handler function.

 @ida <node5#sss:optin_main> IDASetErrHandlerFn
 *)
type error_details = {
    error_code : int;
    module_name : string;
    function_name : string;
    error_message : string;
  }

(** {3:exceptions Exceptions} *)

(** @ida <node5#sss:ida> IDA_ILL_INPUT *)
exception IllInput

(** @ida <node5#sss:ida> IDA_TOO_CLOSE *)
exception TooClose

(** @ida <node5#sss:ida> IDA_TOO_MUCH_WORK *)
exception TooMuchWork

(** @ida <node5#sss:ida> IDA_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** @ida <node5#sss:ida> IDA_ERR_FAIL *)
exception ErrFailure                

(** @ida <node5#sss:ida> IDA_CONV_FAIL *)
exception ConvergenceFailure        

(** @ida <node5#sss:ida> IDA_LINIT_FAIL *)
exception LinearInitFailure         

(** @ida <node5#sss:ida> IDA_LSETUP_FAIL *)
exception LinearSetupFailure        

(** @ida <node5#sss:ida> IDA_LSOLVE_FAIL *)
exception LinearSolveFailure        

(** @ida <node5#sss:ida> IDA_RES_FAIL *)
exception ResFuncFailure

(** @ida <node5#sss:ida> IDA_FIRST_RES_FAIL *)
exception FirstResFuncFailure       

(** @ida <node5#sss:ida> IDA_REP_RES_ERR *)
exception RepeatedResFuncErr        

(** @ida <node5#sss:ida> IDA_UNREC_RESFUNC_ERR *)
exception UnrecoverableResFuncErr   

(** @ida <node5#sss:ida> IDA_RTFUNC_FAIL *)
exception RootFuncFailure           

exception BadK      (** k is not in the range 0, 1, ..., q_u (IDA_BAD_K)
                        @ida <node5#ss:optional_dky> IDAGetDky *)

exception BadT      (** t is not in the interval
                        \[t_n - h_u, t_n\] (IDA_BAD_T)
                        @ida <node5#ss:optional_dky> IDAGetDky *)
exception BadDky    (** invalid dky argument (IDA_BAD_DKY)
                        @ida <node5#ss:optional_dky> IDAGetDky *)

(**
 This exception may be thrown inside the RES callback function (f)
 to indicate that one or more derivatives cannot be calculated at
 the given time offset. *)
exception RecoverableFailure


(**
 Thrown by the getrf functions if a zero diagonal element is encountered during
 factorization. The argument indicates the column index (from 1).

 @ida <node9#ss:dense> DenseGETRF/denseGETRF 
 *)
exception ZeroDiagonalElement of int

include module type of Dls
  with type Bandmatrix.t = Dls.Bandmatrix.t
  and  type Directbandmatrix.t = Dls.Directbandmatrix.t
  and  type Densematrix.t = Dls.Densematrix.t
  and  type Directdensematrix.t = Dls.Directdensematrix.t
