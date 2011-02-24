(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*       Timothy Bourke (INRIA Rennes) and Marc Pouzet (LIENS)         *)
(*                                                                     *)
(*  Copyright 2011 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

(** Vector-independent types and values for the CVODE solver.

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 *)

include module type of Sundials

(**
 Specify a linear multistep method in {!Cvode_serial.init} and
 {!Cvode_nvector.init}.

 @see <CVODE_DOC_ROOT(node5#ss:ivp_sol)> IVP Solution
 @see <CVODE_DOC_ROOT(node5#sss:cvodemalloc)> CVodeCreate
 *)
type lmm =
  | Adams   (** Non-stiff systems; Adams-Moulton formulas *)
  | BDF     (** Stiff systems;     Backward Differentiation Formulas *)

(**
 Specify a solution method in {!Cvode_serial.init} and {!Cvode_nvector.init}.

 @see <CVODE_DOC_ROOT(node5#ss:ivp_sol)> IVP Solution
 @see <CVODE_DOC_ROOT(node5#sss:cvodemalloc)> CVodeCreate
 *)
type iter =
  | Newton of linear_solver (** Newton iteration with a given linear solver *)
  | Functional              (** Functional iteration (non-stiff systems only) *)

(**
 Specify a linear solver in {!Cvode_serial.init}, {!Cvode_serial.set_iter_type},
 {!Cvode_nvector.init}, and {!Cvode_nvector.set_iter_type}.

 The Lapack solvers require that both Sundials and the Ocaml interface were
 built to link with a LAPACK library.

 The Banded Krylov solvers imply an additional call to
 {{:CVODE_DOC_ROOT(node5#sss:cvbandpre)} CVBandPrecInit}.

 @see <CVODE_DOC_ROOT(node5#sss:lin_solve_init)> Linear Solver Specification
                                                 Functions
 *)
and linear_solver =
  | Dense                                   (** Direct with dense matrix *)
  | LapackDense                             (** Direct with dense matrix,
                                                with Lapack *)

  | Band of bandrange                       (** Direct with banded matrix *)
  | LapackBand of bandrange                 (** Direct with banded matrix
                                                with Lapack *)

  | Diag                                    (** Diagonal approximation
                                                of the Jacobian *)
  | Spgmr of sprange                        (** Krylov Spils solver: SPGMR *)
  | Spbcg of sprange                        (** Krylov Spils solver: SPBCG *)
  | Sptfqmr of sprange                      (** Krylov Spils solver: SPFQMR *)

  | BandedSpgmr of sprange * bandrange      (** Krylov Spils solver
                                                with banded matrix: SPGMR *)
  | BandedSpbcg of sprange * bandrange      (** Krylov Spils solver
                                                with banded matrix: SPBCG *)
  | BandedSptfqmr of sprange * bandrange    (** Krylov Spils solver
                                                with banded matrix: SPTFQMR *)

and bandrange = {
    mupper : int;
    mlower : int
  }

and sprange = {
    pretype : preconditioning_type;
    maxl: int
  }

and preconditioning_type =
  | PrecNone
  | PrecLeft
  | PrecRight
  | PrecBoth

(**
 Possible values returned when a solver step function ({!Cvode_serial.normal},
 {!Cvode_serial.one_step}, {!Cvode_nvector.normal}, {!Cvode_nvector.one_step})
 succeeds.
 Failures are indicated by exceptions.

 @see <CVODE_DOC_ROOT(node5#sss:cvode)> CVode
 *)
type solver_result =
  | Continue            (** CV_SUCCESS *)
  | RootsFound          (** CV_ROOT_RETURN *)
  | StopTimeReached     (** CV_TSTOP_RETURN *)


(**
 Exceptions that may be thrown when solver step function
 ({!Cvode_serial.normal}, {!Cvode_serial.one_step}, {!Cvode_nvector.normal},
  {!Cvode_nvector.one_step}) fails.

 @see <CVODE_DOC_ROOT(node5#sss:cvode)> CVode
 *)
exception IllInput                  (** CV_ILL_INPUT *)
exception TooClose                  (** CV_TOO_CLOSE *)
exception TooMuchWork               (** CV_TOO_MUCH_WORK *)
exception TooMuchAccuracy           (** CV_TOO_MUCH_ACC *)
exception ErrFailure                (** CV_ERR_FAIL *)
exception ConvergenceFailure        (** CV_CONV_FAIL *)
exception LinearInitFailure         (** CV_LINIT_FAIL *)
exception LinearSetupFailure        (** CV_LSETUP_FAIL *)
exception LinearSolveFailure        (** CV_LSOLVE_FAIL *)
exception RhsFuncErr                (** CV_RHSFUNC_FAIL *)
exception FirstRhsFuncFailure       (** CV_FIRST_RHSFUNC_FAIL *)
exception RepeatedRhsFuncErr        (** CV_REPTD_RHSFUNC_ERR *)
exception UnrecoverableRhsFuncErr   (** CV_UNREC_RHSFUNC_ERR *)
exception RootFuncFailure           (** CV_RTFUNC_FAIL *)

(**
 Type of values passed to a registered error handler function
 ({!Cvode_serial.set_err_handler_fn}, {!Cvode_nvector.set_err_handler_fn}).

 @see <CVODE_DOC_ROOT(node5#sss:optin_main)> CVodeSetErrHandlerFn
 *)
type error_details = {
    error_code : int;
    module_name : string;
    function_name : string;
    error_message : string;
  }


(**
 @see <CVODE_DOC_ROOT(node5#s:)> 
 *)
type root_direction =
  | Increasing
  | Decreasing
  | IncreasingOrDecreasing

(* get_dky exceptions *)
(**
 @see <CVODE_DOC_ROOT(node5#s:)> 
 *)
exception BadK
exception BadT
exception BadDky

(**
 @see <CVODE_DOC_ROOT(node5#s:)> 
 *)
val no_roots : (int * ('a -> 'b -> 'c -> unit))

(**
 @see <CVODE_DOC_ROOT(node5#s:)> 
 *)
(* Throw inside the f callback if the derivatives cannot be calculated at
   the given time. *)
exception RecoverableFailure

(**
 @see <CVODE_DOC_ROOT(node5#s:)> 
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

(** {2 Data structures for direct linear solvers} *)

(**
 Thrown by GETRF routines for a zero diagonal element at the given column index.

 @see <CVODE_DOC_ROOT(node5#s:)> 
 *)
exception ZeroDiagonalElement of int

(** {3 Dense linear solver} *)

module Densematrix :
  sig
    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    type t

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val new_dense_mat  : int * int -> t

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val print_mat      : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val set_to_zero    : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val add_identity   : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val copy     : t -> t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val scale    : float -> t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val getrf    : t -> int_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val getrs    : t -> int_array -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val potrf    : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val potrs    : t -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val geqrf    : t -> real_array -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    type ormqr = {
          beta : real_array;
          vn   : real_array;
          vm   : real_array;
          work : real_array;
        }

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val ormqr : t -> ormqr -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val get : t -> (int * int) -> float

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val set : t -> (int * int) -> float -> unit

    module Direct :
      sig
        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        type t

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val new_dense_mat  : int * int -> t

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val get : t -> (int * int) -> float

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val set : t -> (int * int) -> float -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val copy  : t -> t -> int * int -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val scale : float -> t -> int * int -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val add_identity : t -> int -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val getrf : t -> int * int -> int_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val getrs : t -> int -> int_array -> real_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val potrf : t -> int -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val potrs : t -> int -> real_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val geqrf : t -> int * int -> real_array -> real_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val ormqr : t -> int * int -> ormqr -> unit
      end
  end

(** {3 Banded linear solver} *)

module Bandmatrix :
  sig
    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    type t

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val new_band_mat : int * int * int * int -> t (* n, mu, ml, smu *)

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val print_mat : t -> unit


    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val set_to_zero    : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val add_identity   : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val copy : t -> t -> int -> int -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val scale : float -> t -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val gbtrf : t -> int_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val gbtrs : t -> int_array -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val get : t -> (int * int) -> float

    (**
     @see <CVODE_DOC_ROOT(node5#s:)> 
     *)
    val set : t -> (int * int) -> float -> unit

    module Col :
      sig
        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        type c

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val get_col : t -> int -> c

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val get : c -> int -> int -> float

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val set : c -> int -> int -> float -> unit
      end

    module Direct :
      sig
        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        type t

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val new_band_mat : int * int * int -> t (* n smu ml *)

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val get : t -> (int * int) -> float

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val set : t -> (int * int) -> float -> unit

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val copy : t -> t -> int -> int -> int -> int -> int -> unit
               (*  a    b    n     a_smu  b_smu  copymu  copyml *)

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val scale : float -> t -> int -> int -> int -> int -> unit
               (*  c         a    n      mu     ml     smu *)

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val add_identity : t -> int -> int -> unit
               (*          a    n      smu *)

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val gbtrf : t -> int -> int -> int -> int -> int_array -> unit
               (*   a    n      mu     ml     smu    p *)

        (**
         @see <CVODE_DOC_ROOT(node5#s:)> 
         *)
        val gbtrs
            : t -> int -> int -> int -> int_array -> real_array -> unit
            (*a    n      smu    ml     p            b *)
      end
  end

