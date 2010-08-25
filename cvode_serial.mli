(* Aug 2010, Timothy Bourke (INRIA) *)

val kind : (float, Bigarray.float64_elt) Bigarray.kind
val layout : Bigarray.c_layout Bigarray.layout
type c_array =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
val empty : c_array

type val_array = c_array
type der_array = c_array
type rootval_array = c_array

val create : int -> c_array
val of_array : float array -> c_array 
val fill : c_array -> float -> unit

val print_results : float -> c_array -> unit

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

type root_direction =
| Rising
| Falling
| RisingAndFalling

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
val set_tolerances : session -> float -> c_array -> unit
val get_roots : session -> Roots.t -> unit

val advance : session -> float -> val_array -> float * bool
val step : session -> float -> val_array -> float * bool
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

