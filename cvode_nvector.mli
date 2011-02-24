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

(** CVODE interface for abstract NVectors

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 *)

include module type of Cvode

type 'a nvector = 'a Nvector.nvector
type 'a session

val init :
    lmm
    -> iter
    -> (float -> 'a -> 'a -> unit)
    -> (int * (float -> 'a -> Sundials.Roots.val_array -> unit))
    -> 'a nvector
    -> 'a session

val init' :
    lmm
    -> iter
    -> (float -> 'a -> 'a -> unit)
    -> (int * (float -> 'a -> Sundials.Roots.val_array -> unit))
    -> 'a nvector
    -> float (* start time *)
    -> 'a session

val nroots : 'a session -> int
val neqs : 'a session -> int

val reinit : 'a session -> float -> 'a nvector -> unit

val sv_tolerances : 'a session -> float -> 'a nvector -> unit
val ss_tolerances : 'a session -> float -> float -> unit
val wf_tolerances : 'a session -> ('a -> 'a -> unit) -> unit

val get_root_info : 'a session -> Sundials.Roots.t -> unit

val normal : 'a session -> float -> 'a nvector -> float * solver_result
val one_step : 'a session -> float -> 'a nvector -> float * solver_result

val get_dky : 'a session -> float -> int -> 'a nvector -> unit

val get_work_space          : 'a session -> int * int
val get_num_steps           : 'a session -> int
val get_num_rhs_evals       : 'a session -> int
val get_num_lin_solv_setups : 'a session -> int
val get_num_err_test_fails  : 'a session -> int
val get_last_order          : 'a session -> int
val get_current_order       : 'a session -> int

val get_actual_init_step    : 'a session -> float
val get_last_step           : 'a session -> float
val get_current_step        : 'a session -> float
val get_current_time        : 'a session -> float

val get_integrator_stats    : 'a session -> integrator_stats
val print_integrator_stats  : 'a session -> unit
val last_step_size          : 'a session -> float
val next_step_size          : 'a session -> float

(* optional input functions *)

(* path to an error file, whether to truncate (true) or append (false) *)
val set_error_file : 'a session -> string -> bool -> unit

(* Error handler function *)
val set_err_handler_fn : 'a session -> (error_details -> unit) -> unit
val clear_err_handler_fn : 'a session -> unit

(* Maximum order for BDF or Adams method *)
val set_max_ord : 'a session -> int -> unit

(* Maximum no. of internal steps before tout *)
val set_max_num_steps : 'a session -> int -> unit

(* Maximum no. of warnings for tn + h = tn *)
val set_max_hnil_warns : 'a session -> int -> unit

(* Flag to activate stability limit detection *)
val set_stab_lim_det : 'a session -> bool -> unit

(* Initial step size *)
val set_init_step : 'a session -> float -> unit

(* Minimum absolute step size *)
val set_min_step : 'a session -> float -> unit

(* Maximum absolute step size *)
val set_max_step : 'a session -> float -> unit

(* Value of tstop *)
val set_stop_time : 'a session -> float -> unit

(* Maximum no. of error test failures *)
val set_max_err_test_fails : 'a session -> int -> unit

(* Maximum no. of nonlinear iterations *)
val set_max_nonlin_iters : 'a session -> int -> unit

(* Maximum no. of convergence failures *)
val set_max_conv_fails : 'a session -> int -> unit

(* Coefficient in the nonlinear convergence test *)
val set_nonlin_conv_coef : 'a session -> float -> unit

(* Nonlinear iteration type *)
val set_iter_type : 'a session -> iter -> unit

(* Direction of zero-crossing *)
val set_root_direction : 'a session -> root_direction array -> unit
val set_all_root_directions : 'a session -> root_direction -> unit

(* Disable rootfinding warnings *)
val set_no_inactive_root_warn : 'a session -> unit

(* Optional output functions *)

(* No. of order reductions due to stability limit detection *)
val get_num_stab_lim_order_reds : 'a session -> int

(* Suggested factor for tolerance scaling *)
val get_tol_scale_factor : 'a session -> float

(* Error weight vector for state variables *)
val get_err_weights : 'a session -> 'a nvector -> unit

(* Estimated local error vector *)
val get_est_local_errors : 'a session -> 'a nvector -> unit

(* No. of nonlinear solver iterations  *)
val get_num_nonlin_solv_iters : 'a session -> int

(* No. of nonlinear convergence failures  *)
val get_num_nonlin_solv_conv_fails : 'a session -> int

(* No. of calls to user root function  *)
val get_num_g_evals : 'a session -> int

(* direct linear solvers functions *)

type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;
    jac_y   : 'a;
    jac_fy  : 'a;
    jac_tmp : 't
  }

type 'a triple_tmp = 'a * 'a * 'a

module Dls :
  sig
    val set_dense_jac_fn :
         'a session
      -> (('a triple_tmp, 'a) jacobian_arg -> Densematrix.t -> unit)
      -> unit

    val clear_dense_jac_fn : 'a session -> unit

    val set_band_jac_fn :
         'a session
      -> (('a triple_tmp, 'a) jacobian_arg -> int -> int -> Bandmatrix.t
          -> unit)
      -> unit

    val clear_band_jac_fn : 'a session -> unit

    val get_work_space : 'a session -> int * int

    (* No. of Jacobian evaluations *)
    val get_num_jac_evals : 'a session -> int

    (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
    val get_num_rhs_evals : 'a session -> int
  end

module Diag :
  sig
    val get_work_space : 'a session -> int * int

    (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
    val get_num_rhs_evals : 'a session -> int
  end

module BandPrec :
  sig
    val get_work_space : 'a session -> int * int

    (* No. of r.h.s. calls for finite diff. banded Jacobian evals. *)
    val get_num_rhs_evals : 'a session -> int
  end

(* iterative linear solvers *)
module Spils :
  sig
    type 'a solve_arg =
      {
        rhs   : 'a;
        gamma : float;
        delta : float;
        left  : bool; (* true: left, false: right *)
      }

    type 'a single_tmp = 'a nvector

    type gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    val set_preconditioner :
      'a session
      -> (('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool)
      -> (('a single_tmp, 'a) jacobian_arg -> 'a solve_arg -> 'a nvector
          -> unit)
      -> unit

    val set_jac_times_vec_fn :
      'a session
      -> (('a single_tmp, 'a) jacobian_arg
          -> 'a (* v *)
          -> 'a (* Jv *)
          -> unit)
      -> unit
    val clear_jac_times_vec_fn : 'a session -> unit

    val set_prec_type : 'a session -> preconditioning_type -> unit

    val set_gs_type :
      'a session -> gramschmidt_type -> unit

    val set_eps_lin : 'a session -> float -> unit

    val set_maxl : 'a session -> int -> unit

    val get_work_space       : 'a session -> int * int
    val get_num_prec_evals   : 'a session -> int
    val get_num_prec_solves  : 'a session -> int
    val get_num_lin_iters    : 'a session -> int
    val get_num_conv_fails   : 'a session -> int
    val get_num_jtimes_evals : 'a session -> int
    val get_num_rhs_evals    : 'a session -> int
  end

