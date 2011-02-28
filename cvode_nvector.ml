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

include Cvode

type 'a nvector = 'a Nvector.nvector

type root_array = Sundials.Roots.t
type root_val_array = Sundials.Roots.val_array

type 'a session

(* interface *)

type handler =
  | RhsFn
  | RootsFn
  | ErrorHandler
  | ErrorWeight
  | JacFn
  | BandJacFn
  | PreSetupFn
  | PreSolveFn
  | JacTimesFn

let handler_name h = match h with
  | RhsFn        -> "cvode_serial_callback_rhsfn"
  | RootsFn      -> "cvode_serial_callback_rootsfn"
  | ErrorHandler -> "cvode_serial_callback_errorhandler"
  | ErrorWeight  -> "cvode_serial_callback_errorweight"
  | JacFn        -> "cvode_serial_callback_jacfn"
  | BandJacFn    -> "cvode_serial_callback_bandjacfn"
  | PreSetupFn   -> "cvode_serial_callback_presetupfn"
  | PreSolveFn   -> "cvode_serial_callback_presolvefn"
  | JacTimesFn   -> "cvode_serial_callback_jactimesfn"

external external_register_handler : 'a session -> handler -> unit
    = "c_register_handler"

let register_handler s h f =
  Callback.register (handler_name h) f;
  external_register_handler s h

external external_init
    : lmm -> iter -> 'a nvector -> int -> float -> 'a session
    = "c_nvec_init"

external nroots : 'a session -> int
    = "c_nroots"

external neqs   : 'a session -> int
    = "c_neqs"

external reinit
    : 'a session -> float -> 'a nvector -> unit
    = "c_nvec_reinit"

external sv_tolerances  : 'a session -> float -> 'a nvector -> unit
    = "c_nvec_sv_tolerances"
external ss_tolerances  : 'a session -> float -> float -> unit
    = "c_ss_tolerances"
external wf_tolerances  : 'a session -> ('a -> 'a -> unit)
                         -> unit
    = "c_nvec_wf_tolerances"

let wf_tolerances s efun =
  register_handler s ErrorWeight efun;
  wf_tolerances s efun

external get_root_info  : 'a session -> Roots.t -> unit
    = "c_get_root_info"

external normal
    : 'a session -> float -> 'a nvector -> float * solver_result
    = "c_nvec_normal"

external one_step
    : 'a session -> float -> 'a nvector -> float * solver_result
    = "c_nvec_one_step"

external get_dky
    : 'a session -> float -> int -> 'a nvector -> unit
    = "c_nvec_get_dky"

let init' lmm iter f (num_roots, roots) y0 t0 =
  let s = external_init lmm iter y0 num_roots t0 in
  register_handler s RhsFn f;
  register_handler s RootsFn roots;
  s

let init lmm iter f roots y0 = init' lmm iter f roots y0 0.0

external get_integrator_stats : 'a session -> integrator_stats
    = "c_get_integrator_stats"

external get_work_space         : 'a session -> int * int
    = "c_get_work_space"

external get_num_steps          : 'a session -> int
    = "c_get_num_steps"

external get_num_rhs_evals      : 'a session -> int
    = "c_get_num_rhs_evals"

external get_num_lin_solv_setups : 'a session -> int
    = "c_get_num_lin_solv_setups"

external get_num_err_test_fails : 'a session -> int
    = "c_get_num_err_test_fails"

external get_last_order         : 'a session -> int
    = "c_get_last_order"

external get_current_order      : 'a session -> int
    = "c_get_current_order"

external get_actual_init_step   : 'a session -> float
    = "c_get_actual_init_step"

external get_last_step          : 'a session -> float
    = "c_get_last_step"

external get_current_step       : 'a session -> float
    = "c_get_current_step"

external get_current_time       : 'a session -> float
    = "c_get_current_time"

let print_integrator_stats s =
  let stats = get_integrator_stats s
  in
    Printf.printf "num_steps = %d\n"           stats.num_steps;
    Printf.printf "num_rhs_evals = %d\n"       stats.num_rhs_evals;
    Printf.printf "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.printf "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.printf "last_order = %d\n"          stats.last_order;
    Printf.printf "current_order = %d\n"       stats.current_order;
    Printf.printf "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.printf "last_step = %e\n"           stats.last_step;
    Printf.printf "current_step = %e\n"        stats.current_step;
    Printf.printf "current_time = %e\n"        stats.current_time;

external set_error_file : 'a session -> string -> bool -> unit
    = "c_set_error_file"

external enable_err_handler_fn : 'a session -> unit
    = "c_enable_err_handler_fn"

external clear_err_handler_fn  : 'a session -> unit
    = "c_disable_err_handler_fn"

let set_err_handler_fn s errh =
  register_handler s ErrorHandler errh;
  enable_err_handler_fn s

external set_max_ord            : 'a session -> int -> unit
    = "c_set_max_ord"
external set_max_num_steps      : 'a session -> int -> unit
    = "c_set_max_num_steps"
external set_max_hnil_warns     : 'a session -> int -> unit
    = "c_set_max_hnil_warns"
external set_stab_lim_det       : 'a session -> bool -> unit
    = "c_set_stab_lim_det"
external set_init_step          : 'a session -> float -> unit
    = "c_set_init_step"
external set_min_step           : 'a session -> float -> unit
    = "c_set_min_step"
external set_max_step           : 'a session -> float -> unit
    = "c_set_max_step"
external set_stop_time          : 'a session -> float -> unit
    = "c_set_stop_time"
external set_max_err_test_fails : 'a session -> int -> unit
    = "c_set_max_err_test_fails"
external set_max_nonlin_iters   : 'a session -> int -> unit
    = "c_set_max_nonlin_iters"
external set_max_conv_fails     : 'a session -> int -> unit
    = "c_set_max_conv_fails"
external set_nonlin_conv_coef   : 'a session -> float -> unit
    = "c_set_nonlin_conv_coef"
external set_iter_type          : 'a session -> iter -> unit
    = "c_set_iter_type"

external set_root_direction' : 'a session -> int_array -> unit
    = "c_set_root_direction"

let int_of_root_direction x =
  match x with
  | Increasing -> 1l
  | Decreasing -> -1l
  | IncreasingOrDecreasing -> 0l

let set_root_direction s rda =
  let n = nroots s in
  let rdirs = make_int_array n in
  if (n > Array.length rda)
    then Bigarray.Array1.fill rdirs
            (int_of_root_direction IncreasingOrDecreasing);
  Array.iteri (fun i v -> rdirs.{i} <- int_of_root_direction v) rda;
  set_root_direction' s rdirs

let set_all_root_directions s rd =
  let rdirs = make_int_array (nroots s) in
  Bigarray.Array1.fill rdirs (int_of_root_direction rd);
  set_root_direction' s rdirs

external set_no_inactive_root_warn      : 'a session -> unit
    = "c_set_no_inactive_root_warn"

external get_num_stab_lim_order_reds    : 'a session -> int
    = "c_get_num_stab_lim_order_reds"

external get_tol_scale_factor           : 'a session -> float
    = "c_get_tol_scale_factor"

external get_err_weights                : 'a session -> 'a nvector -> unit
    = "c_nvec_get_err_weights"

external get_est_local_errors           : 'a session -> 'a nvector -> unit
    = "c_nvec_get_est_local_errors"

external get_num_nonlin_solv_iters      : 'a session -> int
    = "c_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : 'a session -> int
    = "c_get_num_nonlin_solv_conv_fails"

external get_num_g_evals                : 'a session -> int
    = "c_get_num_g_evals"

type 'a single_tmp = 'a
type 'a triple_tmp = 'a * 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;
    jac_y   : 'a;
    jac_fy  : 'a;
    jac_tmp : 't
  }

module Dls =
  struct
    external enable_dense_jac_fn    : 'a session -> unit
        = "c_nvec_dls_enable_dense_jac_fn"

    external clear_dense_jac_fn : 'a session -> unit
        = "c_nvec_dls_disable_dense_jac_fn"

    let set_dense_jac_fn s f =
        register_handler s JacFn f;
        enable_dense_jac_fn s

    external enable_band_jac_fn     : 'a session -> unit
        = "c_nvec_dls_enable_band_jac_fn"

    external clear_band_jac_fn : 'a session -> unit
        = "c_nvec_dls_disable_band_jac_fn"

    let set_band_jac_fn s f =
        register_handler s BandJacFn f;
        enable_band_jac_fn s

    external get_work_space : 'a session -> int * int
        = "c_dls_get_work_space"

    external get_num_jac_evals      : 'a session -> int
        = "c_dls_get_num_jac_evals"

    external get_num_rhs_evals      : 'a session -> int
        = "c_dls_get_num_rhs_evals"
  end

module Diag =
  struct
    external get_work_space         : 'a session -> int * int
        = "c_diag_get_work_space"

    external get_num_rhs_evals      : 'a session -> int
        = "c_diag_get_num_rhs_evals"
  end

module BandPrec =
  struct
    external get_work_space : 'a session -> int * int
        = "c_bandprec_get_work_space"

    external get_num_rhs_evals      : 'a session -> int
        = "c_bandprec_get_num_rhs_evals"
  end

module Spils =
  struct
    type 'a solve_arg =
      {
        rhs   : 'a;
        gamma : float;
        delta : float;
        left  : bool;
      }

    type gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    external enable_preconditioner      : 'a session -> unit
        = "c_nvec_enable_preconditioner"

    let set_preconditioner s fsetup fsolve =
        register_handler s PreSetupFn fsetup;
        register_handler s PreSolveFn fsolve;
        enable_preconditioner s

    external clear_jac_times_vec_fn : 'a session -> unit
        = "c_nvec_disable_jac_times_vec_fn"

    external enable_jac_times_vec_fn    : 'a session -> unit
        = "c_nvec_enable_jac_times_vec_fn"

    let set_jac_times_vec_fn s f =
        register_handler s JacTimesFn f;
        enable_jac_times_vec_fn s

    external set_prec_type : 'a session -> preconditioning_type -> unit
        = "c_set_prec_type"

    external set_gs_type : 'a session -> gramschmidt_type -> unit
        = "c_set_gs_type"

    external set_eps_lin            : 'a session -> float -> unit
        = "c_set_eps_lin"

    external set_maxl               : 'a session -> int -> unit
        = "c_set_maxl"

    external get_num_lin_iters      : 'a session -> int
        = "c_spils_get_num_lin_iters"

    external get_num_conv_fails     : 'a session -> int
        = "c_spils_get_num_conv_fails"

    external get_work_space         : 'a session -> int * int
        = "c_spils_get_work_space"

    external get_num_prec_evals     : 'a session -> int
        = "c_spils_get_num_prec_evals"

    external get_num_prec_solves    : 'a session -> int
        = "c_spils_get_num_prec_solves"

    external get_num_jtimes_evals   : 'a session -> int
        = "c_spils_get_num_jtimes_evals"

    external get_num_rhs_evals      : 'a session -> int
        = "c_spils_get_num_rhs_evals"

  end

