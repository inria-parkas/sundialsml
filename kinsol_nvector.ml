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

type 'a session

type 'a nvector = 'a Nvector.nvector

type 'a single_tmp = 'a
type 'a double_tmp = 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_u   : 'a nvector;
    jac_fu  : 'a nvector;
    jac_tmp : 't
  }

type 'a linear_solver =
  (* TODO: remove from the nvector version. *)
  | Dense of ('a dense_jac_fn) option
  (* TODO: remove from the nvector version. *)
  | LapackDense of ('a dense_jac_fn) option
  (* TODO: remove from the nvector version. *)
  | Band of bandrange * ('a band_jac_fn) option
  (* TODO: remove from the nvector version. *)
  | LapackBand of bandrange * ('a band_jac_fn) option
  | Spgmr of int option * 'a spils_callbacks
  | Spbcg of 'a spils_callbacks
  | Sptfqmr of 'a spils_callbacks

(* TODO: remove from the nvector version. *)
and 'a dense_jac_fn = ('a double_tmp, 'a) jacobian_arg
                      -> Dls.DenseMatrix.t -> unit

(* TODO: remove from the nvector version. *)
and 'a band_jac_fn = ('a double_tmp, 'a) jacobian_arg -> int
                      -> int -> Dls.BandMatrix.t -> unit

(* TODO: remove from the nvector version. *)
and bandrange = { mupper : int; mlower : int; }

and 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                      -> 'a nvector -> unit) option;

    prec_setup_fn : (('a double_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                     -> bool) option;

    jac_times_vec_fn : ('a -> 'a -> 'a -> bool -> bool) option;
  }

and 'a prec_solve_arg = { uscale : 'a nvector; fscale : 'a nvector; }

let spils_no_precond = {
  prec_solve_fn = None;
  prec_setup_fn = None;
  jac_times_vec_fn = None;
}

type print_level =
  | NoInformation     (* 0 *)
  | ShowScaledNorms   (* 1 *)
  | ShowScaledDFNorm  (* 2 *)
  | ShowGlobalValues  (* 3 *)

type eta_params = {
  egamma : float option;
  ealpha : float option;
}

type eta_choice =
  | EtaChoice1                  (* KIN_ETACHOICE1 *)
  | EtaChoice2 of eta_params    (* KIN_ETACHOICE2 *)
  | EtaConstant of float option (* KIN_ETACONSTANT *)

(* TODO: remove from the nvector version. *)
module Dls =
  struct
    (* TODO: Eliminate this function and just provide reinit? *)
    external set_dense_jac_fn : 'a session -> 'a dense_jac_fn -> unit
      = "c_kinsol_dls_set_dense_jac_fn"

    (* TODO: Eliminate this function and just provide reinit? *)
    external clear_dense_jac_fn : 'a session -> unit
      = "c_kinsol_dls_clear_dense_jac_fn"

    (* TODO: Eliminate this function and just provide reinit? *)
    external set_band_jac_fn : 'a session -> 'a band_jac_fn -> unit
      = "c_kinsol_dls_set_band_jac_fn"

    (* TODO: Eliminate this function and just provide reinit? *)
    external clear_band_jac_fn : 'a session -> unit
      = "c_kinsol_dls_clear_band_jac_fn"

    external get_work_space : 'a session -> int * int
      = "c_kinsol_dls_get_work_space"

    external get_num_jac_evals : 'a session -> int
      = "c_kinsol_dls_get_num_jac_evals"

    external get_num_func_evals : 'a session -> int
      = "c_kinsol_dls_get_num_func_evals"
  end

module Spils =
  struct
    type 'a solve_arg = 'a prec_solve_arg =
      {
        uscale : 'a nvector;
        fscale : 'a nvector;
      }

    (* TODO: CVODE/IDA versions to allow optional setup function. *)
    (* TODO: Eliminate this function and just provide reinit? *)
    external set_preconditioner :
            'a session
            -> (('a double_tmp, 'a) jacobian_arg -> 'a prec_solve_arg -> bool) option
            -> (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                  -> 'a nvector -> unit)
            -> unit
        = "c_nvec_kinsol_spils_set_preconditioner"

    (* TODO: add a clear_preconditioner function? *)

    (* TODO: Eliminate this function and just provide reinit? *)
    external set_jac_times_vec_fn :
      'a session -> ('a -> 'a -> 'a -> bool -> bool) -> unit
        = "c_nvec_kinsol_spils_set_jac_times_vec_fn"

    (* TODO: add a clear_jac_times_vec_fn function? *)

    external get_work_space       : 'a session -> int * int
        = "c_kinsol_spils_get_work_space"

    external get_num_lin_iters    : 'a session -> int
        = "c_kinsol_spils_get_num_lin_iters"

    external get_num_conv_fails   : 'a session -> int
        = "c_kinsol_spils_get_num_conv_fails"

    external get_num_prec_evals   : 'a session -> int
        = "c_kinsol_spils_get_num_prec_evals"

    external get_num_prec_solves  : 'a session -> int
        = "c_kinsol_spils_get_num_prec_solves"

    external get_num_jtimes_evals : 'a session -> int
        = "c_kinsol_spils_get_num_jtimes_evals"

    external get_num_func_evals    : 'a session -> int
        = "c_kinsol_spils_get_num_func_evals"
  end

external init :
    'a linear_solver
    -> ('a -> 'a -> unit)
    -> 'a nvector
    -> 'a session
  = "c_nvec_kinsol_init"

type result =
  | Success           (** KIN_SUCCESS *)
  | InitialGuessOK    (** KIN_INITIAL_GUESS_OK *)
  | StoppedOnStepTol  (** KIN_STEP_LT_STPTOL *)

external solve :
    'a session
    -> 'a nvector
    -> bool
    -> 'a nvector
    -> 'a nvector
    -> result
  = "c_nvec_kinsol_solve"

external set_error_file : 'a session -> string -> bool -> unit
  = "c_kinsol_set_error_file"

external set_err_handler_fn
  : 'a session -> (Sundials.error_details -> unit) -> unit
  = "c_kinsol_set_err_handler_fn"

external clear_err_handler_fn : 'a session -> unit
  = "c_kinsol_clear_err_handler_fn"

external set_info_file : 'a session -> string -> bool -> unit
  = "c_kinsol_set_info_file"

external set_info_handler_fn
  : 'a session -> (Sundials.error_details -> unit) -> unit
  = "c_kinsol_set_info_handler_fn"

external clear_info_handler_fn : 'a session -> unit
  = "c_kinsol_clear_info_handler_fn"

external set_print_level : 'a session -> print_level -> unit
  = "c_kinsol_set_print_level"

external set_num_max_iters : 'a session -> int -> unit
  = "c_kinsol_set_num_max_iters"

external set_no_init_setup : 'a session -> bool -> unit
  = "c_kinsol_set_no_init_setup"

(* TODO: remove from the nvector version. *)
external set_no_res_mon : 'a session -> bool -> unit
  = "c_kinsol_set_no_res_mon"

external set_max_setup_calls : 'a session -> int option -> unit
  = "c_kinsol_set_max_setup_calls"

(* TODO: remove from the nvector version. *)
external set_max_sub_setup_calls : 'a session -> int option -> unit
  = "c_kinsol_set_max_sub_setup_calls"

external set_eta_choice : 'a session -> eta_choice -> unit
  = "c_kinsol_set_eta_choice"

external set_res_mon_const_value : 'a session -> float option -> unit
  = "c_kinsol_set_res_mon_const_value"

external set_res_mon_params : 'a session -> float option -> float option -> unit
  = "c_kinsol_set_res_mon_params"

external set_no_min_eps : 'a session -> bool -> unit
  = "c_kinsol_set_no_min_eps"

external set_max_newton_step : 'a session -> float option -> unit
  = "c_kinsol_set_max_newton_step"

external set_max_beta_fails : 'a session -> float option -> unit
  = "c_kinsol_set_max_beta_fails"

external set_rel_err_func : 'a session -> float option -> unit
  = "c_kinsol_set_rel_err_func"

external set_func_norm_tol : 'a session -> float option -> unit
  = "c_kinsol_set_func_norm_tol"

external set_scaled_step_tol : 'a session -> float option -> unit
  = "c_kinsol_set_scaled_step_tol"

external set_constraints : 'a session -> 'a -> unit
  = "c_kinsol_set_constraints"

external set_sys_func : 'a session -> 'a -> unit
  = "c_kinsol_set_sys_func"

external get_work_space : 'a session -> int * int
  = "c_kinsol_get_work_space"

external get_num_func_evals : 'a session -> int
  = "c_kinsol_get_num_func_evals"

external get_num_nonlin_solv_iters : 'a session -> int
  = "c_kinsol_get_num_nonlin_solv_iters"

external get_num_beta_cond_fails : 'a session -> int
  = "c_kinsol_get_num_beta_cond_fails"

external get_num_backtrack_ops : 'a session -> int
  = "c_kinsol_get_num_backtrack_ops"

external get_func_norm : 'a session -> float
  = "c_kinsol_get_func_norm"

external get_step_length : 'a session -> float
  = "c_kinsol_get_step_length"

