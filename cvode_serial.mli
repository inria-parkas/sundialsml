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

(** CVODE interface for serial NVectors

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 *)

include module type of Cvode

type nvec = Sundials.Carray.t
type val_array = Sundials.Carray.t
type der_array = Sundials.Carray.t

type session

val init :
    lmm
    -> iter
    -> (float -> val_array -> der_array -> unit)
    -> (int * (float -> val_array -> Sundials.Roots.val_array -> unit))
    -> val_array
    -> session

val init' :
    lmm
    -> iter
    -> (float -> val_array -> der_array -> unit)
    -> (int * (float -> val_array -> Sundials.Roots.val_array -> unit))
    -> val_array
    -> float (* start time *)
    -> session

val nroots : session -> int
val neqs : session -> int

val reinit : session -> float -> val_array -> unit

val sv_tolerances : session -> float -> nvec -> unit
val ss_tolerances : session -> float -> float -> unit
val wf_tolerances : session -> (val_array -> nvec -> unit) -> unit

val get_root_info : session -> Sundials.Roots.t -> unit

val normal : session -> float -> val_array -> float * solver_result
val one_step : session -> float -> val_array -> float * solver_result

val get_dky : session -> float -> int -> nvec -> unit

val get_work_space          : session -> int * int
val get_num_steps           : session -> int
val get_num_rhs_evals       : session -> int
val get_num_lin_solv_setups : session -> int
val get_num_err_test_fails  : session -> int
val get_last_order          : session -> int
val get_current_order       : session -> int

val get_actual_init_step    : session -> float
val get_last_step           : session -> float
val get_current_step        : session -> float
val get_current_time        : session -> float

val get_integrator_stats    : session -> integrator_stats
val print_integrator_stats  : session -> unit
val last_step_size          : session -> float
val next_step_size          : session -> float

(* optional input functions *)

(* path to an error file, whether to truncate (true) or append (false) *)
val set_error_file : session -> string -> bool -> unit

(* Error handler function *)
val set_err_handler_fn : session -> (error_details -> unit) -> unit
val clear_err_handler_fn : session -> unit

(* Maximum order for BDF or Adams method *)
val set_max_ord : session -> int -> unit

(* Maximum no. of internal steps before tout *)
val set_max_num_steps : session -> int -> unit

(* Maximum no. of warnings for tn + h = tn *)
val set_max_hnil_warns : session -> int -> unit

(* Flag to activate stability limit detection *)
val set_stab_lim_det : session -> bool -> unit

(* Initial step size *)
val set_init_step : session -> float -> unit

(* Minimum absolute step size *)
val set_min_step : session -> float -> unit

(* Maximum absolute step size *)
val set_max_step : session -> float -> unit

(* Value of tstop *)
val set_stop_time : session -> float -> unit

(* Maximum no. of error test failures *)
val set_max_err_test_fails : session -> int -> unit

(* Maximum no. of nonlinear iterations *)
val set_max_nonlin_iters : session -> int -> unit

(* Maximum no. of convergence failures *)
val set_max_conv_fails : session -> int -> unit

(* Coefficient in the nonlinear convergence test *)
val set_nonlin_conv_coef : session -> float -> unit

(* Nonlinear iteration type *)
val set_iter_type : session -> iter -> unit

(* Direction of zero-crossing *)
val set_root_direction : session -> root_direction array -> unit
val set_all_root_directions : session -> root_direction -> unit

(* Disable rootfinding warnings *)
val set_no_inactive_root_warn : session -> unit

(* Optional output functions *)

(* No. of order reductions due to stability limit detection *)
val get_num_stab_lim_order_reds : session -> int

(* Suggested factor for tolerance scaling *)
val get_tol_scale_factor : session -> float

(* Error weight vector for state variables *)
val get_err_weights : session -> nvec -> unit

(* Estimated local error vector *)
val get_est_local_errors : session -> nvec -> unit

(* No. of nonlinear solver iterations  *)
val get_num_nonlin_solv_iters : session -> int

(* No. of nonlinear convergence failures  *)
val get_num_nonlin_solv_conv_fails : session -> int

(* No. of calls to user root function  *)
val get_num_g_evals : session -> int

(* direct linear solvers functions *)

type 't jacobian_arg =
  {
    jac_t   : float;
    jac_y   : val_array;
    jac_fy  : val_array;
    jac_tmp : 't
  }

type triple_tmp = val_array * val_array * val_array

module Dls :
  sig
    val set_dense_jac_fn :
         session
      -> (triple_tmp jacobian_arg -> Densematrix.t -> unit)
      -> unit

    val clear_dense_jac_fn : session -> unit

    val set_band_jac_fn :
         session
      -> (triple_tmp jacobian_arg -> int -> int -> Bandmatrix.t -> unit)
      -> unit

    val clear_band_jac_fn : session -> unit

    val get_work_space : session -> int * int

    (* No. of Jacobian evaluations *)
    val get_num_jac_evals : session -> int

    (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
    val get_num_rhs_evals : session -> int
  end

module Diag :
  sig
    val get_work_space : session -> int * int

    (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
    val get_num_rhs_evals : session -> int
  end

module BandPrec :
  sig
    val get_work_space : session -> int * int

    (* No. of r.h.s. calls for finite diff. banded Jacobian evals. *)
    val get_num_rhs_evals : session -> int
  end

(* iterative linear solvers *)
module Spils :
  sig
    type solve_arg =
      {
        rhs   : val_array;
        gamma : float;
        delta : float;
        left  : bool; (* true: left, false: right *)
      }

    type single_tmp = val_array

    type gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    val set_preconditioner :
      session
      -> (triple_tmp jacobian_arg -> bool -> float -> bool)
      -> (single_tmp jacobian_arg -> solve_arg -> val_array -> unit)
      -> unit

    val set_jac_times_vec_fn :
      session
      -> (single_tmp jacobian_arg
          -> val_array (* v *)
          -> val_array (* Jv *)
          -> unit)
      -> unit
    val clear_jac_times_vec_fn : session -> unit

    val set_prec_type : session -> preconditioning_type -> unit

    val set_gs_type :
      session -> gramschmidt_type -> unit

    val set_eps_lin : session -> float -> unit

    val set_maxl : session -> int -> unit

    val get_work_space       : session -> int * int
    val get_num_prec_evals   : session -> int
    val get_num_prec_solves  : session -> int
    val get_num_lin_iters    : session -> int
    val get_num_conv_fails   : session -> int
    val get_num_jtimes_evals : session -> int
    val get_num_rhs_evals    : session -> int
  end

