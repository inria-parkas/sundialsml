(* Aug 2010, Timothy Bourke (INRIA) *)

include Cvode

type nvec = Sundials.Carray.t
type val_array = Sundials.Carray.t
type der_array = Sundials.Carray.t

type session

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

external external_register_handler : session -> handler -> unit
    = "c_register_handler"

let register_handler s h f =
  Callback.register (handler_name h) f;
  external_register_handler s h

external external_init
    : lmm -> iter -> val_array -> int -> float -> session
    = "c_ba_init"

external nroots         : session -> int
    = "c_nroots"

external neqs           : session -> int
    = "c_neqs"

external reinit
    : session -> float -> val_array -> unit
    = "c_ba_reinit"

external sv_tolerances  : session -> float -> nvec -> unit
    = "c_ba_sv_tolerances"
external ss_tolerances  : session -> float -> float -> unit
    = "c_ss_tolerances"
external wf_tolerances  : session -> (val_array -> nvec -> unit) -> unit
    = "c_ba_wf_tolerances"

let wf_tolerances s efun =
  register_handler s ErrorWeight efun;
  wf_tolerances s efun

external get_root_info  : session -> Roots.t -> unit
    = "c_get_root_info"

external normal
    : session -> float -> val_array -> float * solver_result
    = "c_ba_normal"

external one_step
    : session -> float -> val_array -> float * solver_result
    = "c_ba_one_step"

external get_dky
    : session -> float -> int -> nvec -> unit
    = "c_ba_get_dky"

let init' lmm iter f (num_roots, roots) y0 t0 =
  let s = external_init lmm iter y0 num_roots t0 in
  register_handler s RhsFn f;
  register_handler s RootsFn roots;
  s

let init lmm iter f roots y0 = init' lmm iter f roots y0 0.0

external get_integrator_stats   : session -> integrator_stats
    = "c_get_integrator_stats"

external last_step_size         : session -> float
    = "c_last_step_size"

external next_step_size         : session -> float
    = "c_next_step_size"

external get_work_space         : session -> int * int
    = "c_get_work_space"

external get_num_steps          : session -> int
    = "c_get_num_steps"

external get_num_rhs_evals      : session -> int
    = "c_get_num_rhs_evals"

external get_num_lin_solv_setups : session -> int
    = "c_get_num_lin_solv_setups"

external get_num_err_test_fails : session -> int
    = "c_get_num_err_test_fails"

external get_last_order         : session -> int
    = "c_get_last_order"

external get_current_order      : session -> int
    = "c_get_current_order"

external get_actual_init_step   : session -> float
    = "c_get_actual_init_step"

external get_last_step          : session -> float
    = "c_get_last_step"

external get_current_step       : session -> float
    = "c_get_current_step"

external get_current_time       : session -> float
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

external set_error_file         : session -> string -> bool -> unit
    = "c_set_error_file"

external enable_err_handler_fn  : session -> unit
    = "c_enable_err_handler_fn"

external clear_err_handler_fn  : session -> unit
    = "c_disable_err_handler_fn"

let set_err_handler_fn s errh =
  register_handler s ErrorHandler errh;
  enable_err_handler_fn s

external set_max_ord            : session -> int -> unit
    = "c_set_max_ord"
external set_max_num_steps      : session -> int -> unit
    = "c_set_max_num_steps"
external set_max_hnil_warns     : session -> int -> unit
    = "c_set_max_hnil_warns"
external set_stab_lim_det       : session -> bool -> unit
    = "c_set_stab_lim_det"
external set_init_step          : session -> float -> unit
    = "c_set_init_step"
external set_min_step           : session -> float -> unit
    = "c_set_min_step"
external set_max_step           : session -> float -> unit
    = "c_set_max_step"
external set_stop_time          : session -> float -> unit
    = "c_set_stop_time"
external set_max_err_test_fails : session -> int -> unit
    = "c_set_max_err_test_fails"
external set_max_nonlin_iters   : session -> int -> unit
    = "c_set_max_nonlin_iters"
external set_max_conv_fails     : session -> int -> unit
    = "c_set_max_conv_fails"
external set_nonlin_conv_coef   : session -> float -> unit
    = "c_set_nonlin_conv_coef"
external set_iter_type          : session -> iter -> unit
    = "c_set_iter_type"

external set_root_direction'    : session -> int_array -> unit
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
  let n = nroots s in
  if (n > 0) then begin
    let rdirs = make_int_array n in
    Bigarray.Array1.fill rdirs (int_of_root_direction rd);
    set_root_direction' s rdirs
  end; ()

external set_no_inactive_root_warn      : session -> unit
    = "c_set_no_inactive_root_warn"

external get_num_stab_lim_order_reds    : session -> int
    = "c_get_num_stab_lim_order_reds"

external get_tol_scale_factor           : session -> float
    = "c_get_tol_scale_factor"

external get_err_weights                : session -> nvec -> unit
    = "c_ba_get_err_weights"

external get_est_local_errors           : session -> nvec -> unit
    = "c_ba_get_est_local_errors"

external get_num_nonlin_solv_iters      : session -> int
    = "c_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : session -> int
    = "c_get_num_nonlin_solv_conv_fails"

external get_num_g_evals                : session -> int
    = "c_get_num_g_evals"

type 't jacobian_arg =
  {
    jac_t   : float;
    jac_y   : val_array;
    jac_fy  : val_array;
    jac_tmp : 't
  }

type triple_tmp = val_array * val_array * val_array

module Dls =
  struct
    external enable_dense_jac_fn  : session -> unit
        = "c_ba_dls_enable_dense_jac_fn"

    external clear_dense_jac_fn : session -> unit
        = "c_ba_dls_disable_dense_jac_fn"

    let set_dense_jac_fn s f =
        register_handler s JacFn f;
        enable_dense_jac_fn s

    external enable_band_jac_fn   : session -> unit
        = "c_ba_dls_enable_band_jac_fn"

    external clear_band_jac_fn : session -> unit
        = "c_ba_dls_disable_band_jac_fn"

    let set_band_jac_fn s f =
        register_handler s BandJacFn f;
        enable_band_jac_fn s

    external get_work_space : session -> int * int
        = "c_dls_get_work_space"

    external get_num_jac_evals    : session -> int
        = "c_dls_get_num_jac_evals"

    external get_num_rhs_evals    : session -> int
        = "c_dls_get_num_rhs_evals"
  end

module Diag =
  struct
    external get_work_space       : session -> int * int
        = "c_diag_get_work_space"

    external get_num_rhs_evals    : session -> int
        = "c_diag_get_num_rhs_evals"
  end

module BandPrec =
  struct
    external get_work_space : session -> int * int
        = "c_bandprec_get_work_space"

    external get_num_rhs_evals    : session -> int
        = "c_bandprec_get_num_rhs_evals"
  end

module Spils =
  struct
    type solve_arg =
      {
        rhs   : val_array;
        gamma : float;
        delta : float;
        left  : bool;
      }

    type single_tmp = val_array

    type gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    external enable_preconditioner  : session -> unit
        = "c_ba_enable_preconditioner"

    let set_preconditioner s fsetup fsolve =
        register_handler s PreSetupFn fsetup;
        register_handler s PreSolveFn fsolve;
        enable_preconditioner s

    external enable_jac_times_vec_fn : session -> unit
        = "c_ba_enable_jac_times_vec_fn"

    external clear_jac_times_vec_fn : session -> unit
        = "c_ba_disable_jac_times_vec_fn"

    let set_jac_times_vec_fn s f =
        register_handler s JacTimesFn f;
        enable_jac_times_vec_fn s

    external set_prec_type
        : session -> preconditioning_type -> unit
        = "c_set_prec_type"

    external set_gs_type : session -> gramschmidt_type -> unit
        = "c_set_gs_type"

    external set_eps_lin : session -> float -> unit
        = "c_set_eps_lin"

    external set_maxl   : session -> int -> unit
        = "c_set_maxl"

    external get_num_lin_iters      : session -> int
        = "c_spils_get_num_lin_iters"

    external get_num_conv_fails     : session -> int
        = "c_spils_get_num_conv_fails"

    external get_work_space         : session -> int * int
        = "c_spils_get_work_space"

    external get_num_prec_evals     : session -> int
        = "c_spils_get_num_prec_evals"

    external get_num_prec_solves    : session -> int
        = "c_spils_get_num_prec_solves"

    external get_num_jtimes_evals   : session -> int
        = "c_spils_get_num_jtimes_evals"

    external get_num_rhs_evals      : session -> int
        = "c_spils_get_num_rhs_evals"

  end

