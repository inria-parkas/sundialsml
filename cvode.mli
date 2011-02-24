(* Aug 2010, Timothy Bourke (INRIA) *)

include module type of Sundials

type lmm =
  | Adams
  | BDF

type preconditioning_type =
  | PrecNone
  | PrecLeft
  | PrecRight
  | PrecBoth

type bandrange = { mupper : int; mlower : int }

val sprange_default_maxl : int
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

(* get_dky exceptions *)
exception BadK
exception BadT
exception BadDky

val no_roots : (int * ('a -> 'b -> 'c -> unit))

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

(* direct linear solvers functions *)

(* Thrown by GETRF routines for a zero diagonal element at the given
   column index. *)
exception ZeroDiagonalElement of int

module Densematrix :
  sig
    type t

    val new_dense_mat  : int * int -> t
    val print_mat      : t -> unit

    val set_to_zero    : t -> unit
    val add_identity   : t -> unit

    val copy     : t -> t -> unit
    val scale    : float -> t -> unit
    val getrf    : t -> int_array -> unit
    val getrs    : t -> int_array -> real_array -> unit
    val potrf    : t -> unit
    val potrs    : t -> real_array -> unit
    val geqrf    : t -> real_array -> real_array -> unit

    type ormqr = {
          beta : real_array;
          vn   : real_array;
          vm   : real_array;
          work : real_array;
        }

    val ormqr : t -> ormqr -> unit

    val get : t -> (int * int) -> float
    val set : t -> (int * int) -> float -> unit

    module Direct :
      sig
        type t

        val new_dense_mat  : int * int -> t

        val get : t -> (int * int) -> float
        val set : t -> (int * int) -> float -> unit

        val copy  : t -> t -> int * int -> unit
        val scale : float -> t -> int * int -> unit
        val add_identity : t -> int -> unit
        val getrf : t -> int * int -> int_array -> unit
        val getrs : t -> int -> int_array -> real_array -> unit
        val potrf : t -> int -> unit
        val potrs : t -> int -> real_array -> unit
        val geqrf : t -> int * int -> real_array -> real_array -> unit
        val ormqr : t -> int * int -> ormqr -> unit
      end
  end

module Bandmatrix :
  sig
    type t

    val new_band_mat : int * int * int * int -> t (* n, mu, ml, smu *)
    val print_mat : t -> unit

    val set_to_zero    : t -> unit
    val add_identity   : t -> unit

    val copy : t -> t -> int -> int -> unit
    val scale : float -> t -> unit
    val gbtrf : t -> int_array -> unit
    val gbtrs : t -> int_array -> real_array -> unit

    val get : t -> (int * int) -> float
    val set : t -> (int * int) -> float -> unit

    module Col :
      sig
        type c

        val get_col : t -> int -> c

        val get : c -> int -> int -> float
        val set : c -> int -> int -> float -> unit
      end

    module Direct :
      sig
        type t

        val new_band_mat : int * int * int -> t (* n smu ml *)

        val get : t -> (int * int) -> float
        val set : t -> (int * int) -> float -> unit

        val copy : t -> t -> int -> int -> int -> int -> int -> unit
               (*  a    b    n     a_smu  b_smu  copymu  copyml *)

        val scale : float -> t -> int -> int -> int -> int -> unit
               (*  c         a    n      mu     ml     smu *)

        val add_identity : t -> int -> int -> unit
               (*          a    n      smu *)

        val gbtrf : t -> int -> int -> int -> int -> int_array -> unit
               (*   a    n      mu     ml     smu    p *)

        val gbtrs
            : t -> int -> int -> int -> int_array -> real_array -> unit
            (*a    n      smu    ml     p            b *)
      end
  end

