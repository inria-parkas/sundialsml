(* Aug 2010, Timothy Bourke (INRIA) *)

module Carray :
  sig
    type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

    val kind : (float, Bigarray.float64_elt) Bigarray.kind
    val layout : Bigarray.c_layout Bigarray.layout

    val empty : t
    val create : int -> t
    val of_array : float array -> t
    val fill : t -> float -> unit
    val length : t -> int

    val print_with_time : float -> t -> unit
  end

type val_array = Carray.t
type der_array = Carray.t
type rootval_array = Carray.t

module Roots :
  sig
    type t
    val empty : t
    val create : int -> t
    val print : t -> unit
    val get : t -> int -> bool
    val set : t -> int -> bool -> unit
    val length : t -> int 
    val reset : t -> unit
  end

type lmm =
| Adams
| BDF

type bandrange = {mupper : int; mlower : int}
type sprange = { pretype : int; maxl: int }

type linear_solver =
| Dense
| LapackDense
| Band of bandrange
| LapackBand of bandrange
| Diag
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

type root_direction =
| Increasing
| Decreasing
| IncreasingOrDecreasing

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
exception RhsFuncErr
exception FirstRhsFuncFailure
exception RepeatedRhsFuncErr
exception UnrecoverableRhsFuncErr
exception RootFuncFailure

type session

val no_roots : (int * (float -> val_array -> Roots.t -> int))

(* Throw inside the f callback if the derivatives cannot be calculated at the
   given time. *)
exception RecoverableFailure

val init :
    lmm
    -> iter
    -> (float -> val_array -> der_array -> unit)
    -> (int * (float -> val_array -> rootval_array -> unit))
    -> val_array
    -> session
val nroots : session -> int
val neqs : session -> int

val reinit : session -> float -> val_array -> unit
val set_tolerances : session -> float -> Carray.t -> unit
val get_roots : session -> Roots.t -> unit

val advance : session -> float -> val_array -> float * solver_result
val step : session -> float -> val_array -> float * solver_result
val free : session -> unit

type integrator_stats = {
  steps : int;
  rhs_evals : int;
  linear_solver_setups : int;
  error_test_failures : int;
  last_internal_order : int;
  next_internal_order : int;
  initial_step_size : float;
  last_step_size : float;
  next_step_size : float;
  internal_time : float
}

val integrator_stats : session -> integrator_stats
val last_step_size : session -> float
val next_step_size : session -> float

val print_integrator_stats : session -> unit

(* optional input functions *)

(* path to an error file, whether to truncate (true) or append (false) *)
val set_error_file : session -> string -> bool -> unit 

(* Error handler function *)
val set_error_handler : session -> (error_details -> unit) -> unit 

(* Maximum order for BDF or Adams method *)
val set_max_ord : session -> int -> unit 

(* Maximum no. of internal steps before tout *)
val set_max_num_steps : session -> int -> unit 

(* Maximum no. of warnings for tn + h = tn *)
val set_max_hnil_warns : session -> int -> unit 

(* Flag to activate stability limit detection *)
val set_stability_limit_detection : session -> bool -> unit 

(* Initial step size *)
val set_initial_step_size : session -> float -> unit 

(* Minimum absolute step size *)
val set_min_abs_step_size : session -> float -> unit 

(* Maximum absolute step size *)
val set_max_abs_step_size : session -> float -> unit 

(* Value of tstop *)
val set_stop_time : session -> float -> unit 

(* Maximum no. of error test failures *)
val set_max_error_test_failures : session -> int -> unit 

(* Maximum no. of nonlinear iterations *)
val set_max_nonlinear_iterations : session -> int -> unit 

(* Maximum no. of convergence failures *)
val set_max_convergence_failures : session -> int -> unit 

(* Coefficient in the nonlinear convergence test *)
val set_nonlinear_convergence_coeffficient : session -> float -> unit 

(* Nonlinear iteration type *)
val set_nonlinear_iteration_type : session -> iter -> unit 

(* Direction of zero-crossing *)
val set_root_direction : session -> root_direction array -> unit 
val set_all_root_directions : session -> root_direction -> unit 

(* Disable rootfinding warnings *)
val disable_inactive_root_warnings : session -> unit 

(* TODO:
(* CVDLS linear solvers *)
val CVDlsSetDenseJacFn  : (* Dense Jacobian function  *)
val CVDlsSetBandJacFn  : (* Band Jacobian function  *)

(* CVSPILS linear solvers *)
val CVSpilsSetPreconditioner  : (* Preconditioner functions  *)
val CVSpilsSetJacTimesVecFn  : (* Jacobian-times-vector function  *)
val CVSpilsSetPrecType  : (* Preconditioning type  *)
val CVSpilsSetEpsLin  : (* Ratio between linear and nonlinear tolerances  *)
val CVSpilsSetGSType  : (* Type of Gram-Schmidt orthogonalization(a)  *)
val CVSpilsSetMaxl  : (* Maximum Krylov subspace size(b)  *)
*)

