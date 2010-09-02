(* Aug 2010, Timothy Bourke (INRIA) *)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_serial.h (and code in cvode_serial.c) must also be updated.
 *)

module Carray =
  struct
    type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

    let kind = Bigarray.float64
    let layout = Bigarray.c_layout
    type c_array =
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
    let empty = Bigarray.Array1.create kind layout 0

    let create = Bigarray.Array1.create kind layout
    let of_array = Bigarray.Array1.of_array kind layout

    let fill = Bigarray.Array1.fill

    let length = Bigarray.Array1.dim

    let print_with_time t v =
      Printf.printf "% .8f" t;
      for i = 0 to (length v - 1) do
        Printf.printf "\t% .8f" v.{i}
      done;
      print_newline ()
  end

type val_array = Carray.t
type der_array = Carray.t
type rootval_array = Carray.t

(* root arrays *)

type int_array = (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t
let create_int_array = Bigarray.Array1.create Bigarray.int Carray.layout

module Roots =
  struct
    type t = int_array

    let create = create_int_array
    let empty = create 0

    let get roots i = roots.{i} <> 0
    let get' roots i = roots.{i}

    let set a i v = Bigarray.Array1.set a i (if v then 1 else 0)

    let print v =
      let isroot = get v in
      for i = 0 to (Carray.length v - 1) do
        if i > 0 then print_string " "; 
        print_string (if (isroot i) then "1" else "0")
      done;
      print_newline ()

    let length = Bigarray.Array1.dim
    let reset v = Bigarray.Array1.fill v 0

    let fold_left f a v =
      let rec check (i, a) =
        if i < 0 then a
        else check (i - 1, f a v.{i})
      in
      check (Bigarray.Array1.dim v - 1, a)

    let exists = fold_left (fun a x -> a || x <> 0) false
  end

let no_roots = (0, (fun _ _ _ -> ()))

(* TODO: get these types working in the initialization function *)
type lmm =
| Adams
| BDF

type preconditioning_type =
| PrecNone
| PrecLeft
| PrecRight
| PrecBoth

type bandrange = {mupper : int; mlower : int}
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

exception BadK
exception BadT
exception BadDky

type session
exception RecoverableFailure
exception StopTimeReached
let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("cvode_RecoverableFailure",      RecoverableFailure);

    ("cvode_StopTimeReached",         StopTimeReached);
    ("cvode_IllInput",                IllInput);
    ("cvode_TooClose",                TooClose);
    ("cvode_TooMuchWork",             TooMuchWork);
    ("cvode_TooMuchAccuracy",         TooMuchAccuracy);
    ("cvode_ErrFailure",              ErrFailure);
    ("cvode_ConvergenceFailure",      ConvergenceFailure);
    ("cvode_LinearInitFailure",       LinearInitFailure);
    ("cvode_LinearSetupFailure",      LinearSetupFailure);
    ("cvode_LinearSolveFailure",      LinearSolveFailure);
    ("cvode_RhsFuncErr",              RhsFuncErr);
    ("cvode_FirstRhsFuncFailure",     FirstRhsFuncFailure);
    ("cvode_RepeatedRhsFuncErr",      RepeatedRhsFuncErr);
    ("cvode_UnrecoverableRhsFuncErr", UnrecoverableRhsFuncErr);
    ("cvode_RootFuncFailure",         RootFuncFailure);

    ("cvode_BadK",                    BadK);
    ("cvode_BadT",                    BadT);
    ("cvode_BadDky",                  BadDky);
  ]

(* passing callbacks to c *)

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

(* interface *)

external external_init
    : lmm -> iter -> val_array -> int -> float -> session
    = "c_init"

external nroots : session -> int
    = "c_nroots"

external neqs : session -> int
    = "c_neqs"

external reinit : session -> float -> val_array -> unit
    = "c_re_init"

external sv_tolerances : session -> float -> Carray.t -> unit
    = "c_sv_tolerances"
external ss_tolerances : session -> float -> float -> unit
    = "c_ss_tolerances"

external wf_tolerances : session -> (val_array -> Carray.t -> unit) -> unit
    = "c_wf_tolerances"

let set_wf_tolerances s efun =
  register_handler s ErrorWeight efun;
  wf_tolerances s

external get_root_info : session -> Roots.t -> unit
    = "c_get_root_info"

external free : session -> unit
    = "c_free"

external normal : session -> float -> val_array -> float * solver_result
    = "c_normal"

external one_step : session -> float -> val_array -> float * solver_result
    = "c_one_step"

external get_dky : session -> float -> int -> Carray.t -> unit
    = "c_get_dky"

let init' lmm iter f (num_roots, roots) y0 t0 =
  let s = external_init lmm iter y0 num_roots t0 in
  register_handler s RhsFn f;
  register_handler s RootsFn roots;
  s

let init lmm iter f roots y0 = init' lmm iter f roots y0 0.0

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

external get_integrator_stats : session -> integrator_stats
    = "c_get_integrator_stats"

external last_step_size : session -> float
    = "c_last_step_size"

external next_step_size : session -> float
    = "c_next_step_size"

external get_num_steps : session -> int
    = "c_get_num_steps"

external get_num_rhs_evals : session -> int
    = "c_get_num_rhs_evals"

external get_num_lin_solv_setups : session -> int
    = "c_get_num_lin_solv_setups"

external get_num_err_test_fails : session -> int
    = "c_get_num_err_test_fails"

external get_last_order : session -> int
    = "c_get_last_order"

external get_current_order : session -> int
    = "c_get_current_order"

external get_actual_init_step : session -> float
    = "c_get_actual_init_step"

external get_last_step : session -> float
    = "c_get_last_step"

external get_current_step : session -> float
    = "c_get_current_step"

external get_current_time : session -> float
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

external set_error_file : session -> string -> bool -> unit 
    = "c_set_error_file"

external enable_err_handler_fn : session -> unit 
    = "c_enable_err_handler_fn"

let set_err_handler_fn s errh =
  register_handler s ErrorHandler errh;
  enable_err_handler_fn s

external set_max_ord : session -> int -> unit 
    = "c_set_max_ord"
external set_max_num_steps : session -> int -> unit 
    = "c_set_max_num_steps"
external set_max_hnil_warns : session -> int -> unit 
    = "c_set_max_hnil_warns"
external set_stab_lim_det : session -> bool -> unit 
    = "c_set_stab_lim_det"
external set_init_step: session -> float -> unit 
    = "c_set_init_step"
external set_min_step : session -> float -> unit 
    = "c_set_min_step"
external set_max_step : session -> float -> unit 
    = "c_set_max_step"
external set_stop_time : session -> float -> unit 
    = "c_set_stop_time"
external set_max_err_test_fails : session -> int -> unit 
    = "c_set_max_err_test_fails"
external set_max_nonlin_iters : session -> int -> unit 
    = "c_set_max_nonlin_iters"
external set_max_conv_fails : session -> int -> unit 
    = "c_set_max_conv_fails"
external set_nonlin_conv_coef : session -> float -> unit 
    = "c_set_nonlin_conv_coef "
external set_iter_type : session -> iter -> unit 
    = "c_set_iter_type"

external set_root_direction' : session -> int_array -> unit 
    = "c_set_root_direction"

let int_of_root_direction x =
  match x with
  | Increasing -> 1
  | Decreasing -> -1
  | IncreasingOrDecreasing -> 0
    
let set_root_direction s rda =
  let n = nroots s in
  let rdirs = create_int_array n in
  if (n > Array.length rda)
    then Bigarray.Array1.fill rdirs
            (int_of_root_direction IncreasingOrDecreasing);
  Array.iteri (fun i v -> rdirs.{i} <- int_of_root_direction v) rda;
  set_root_direction' s rdirs

let set_all_root_directions s rd =
  let rdirs = create_int_array (nroots s) in
  Bigarray.Array1.fill rdirs (int_of_root_direction rd);
  set_root_direction' s rdirs

external set_no_inactive_root_warn : session -> unit 
    = "c_set_no_inactive_root_warn"

external get_num_stab_lim_order_reds : session -> int
    = "c_get_num_stab_lim_order_reds"

external get_tol_scale_factor : session -> float
    = "c_get_tol_scale_factor"

external get_err_weights : session -> Carray.t -> unit
    = "c_get_err_weights"

external get_est_local_errors : session -> Carray.t -> unit
    = "c_get_est_local_errors"

external get_num_nonlin_solv_iters : session -> int
    = "c_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : session -> int
    = "c_get_num_nonlin_solv_conv_fails"

external get_num_g_evals : session -> int
    = "c_get_num_g_evals"

(* note: uses DENSE_ELEM rather than the more efficient DENSE_COL. *)
module Densematrix =
  struct
    type t

    external get : t -> (int * int) -> float
        = "c_densematrix_get"

    external set : t -> (int * int) -> float -> unit
        = "c_densematrix_set"
  end

(* note: uses BAND_ELEM rather than the more efficient BAND_COL/BAND_COL_ELEM *)
module Bandmatrix =
  struct
    type t

    external get : t -> (int * int) -> float
        = "c_bandmatrix_get"

    external set : t -> (int * int) -> float -> unit
        = "c_bandmatrix_set"
  end

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
    external enable_dense_jac_fn : session -> unit
        = "c_dls_enable_dense_jac_fn"

    let set_dense_jac_fn s f =
        register_handler s JacFn f;
        enable_dense_jac_fn s

    external enable_band_jac_fn : session -> unit
        = "c_dls_enable_band_jac_fn"

    let set_band_jac_fn s f =
        register_handler s BandJacFn f;
        enable_band_jac_fn s

    external get_num_jac_evals : session -> int
        = "c_dls_get_num_jac_evals"

    external get_num_rhs_evals : session -> int
        = "c_dls_get_num_rhs_evals"
  end

module Diag =
  struct
    external get_num_rhs_evals : session -> int
        = "c_diag_get_num_rhs_evals"
  end

module BandPrec =
  struct
    external get_num_rhs_evals : session -> int
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

    external enable_preconditioner : session -> unit
        = "c_enable_preconditioner"

    let set_preconditioner s fsetup fsolve =
        register_handler s PreSetupFn fsetup;
        register_handler s PreSolveFn fsolve;
        enable_preconditioner s

    external enable_jac_times_vec_fn : session -> unit
        = "c_enable_jac_times_vec_fn"

    let set_jac_times_vec_fn s f =
        register_handler s JacTimesFn f;
        enable_jac_times_vec_fn s

    external set_prec_type : session -> preconditioning_type -> unit
        = "c_set_prec_type"

    external set_gs_type :
        session -> gramschmidt_type -> unit
        = "c_set_gs_type"

    external set_eps_lin : session -> float -> unit
        = "c_set_eps_lin"

    external set_maxl : session -> int -> unit
        = "c_set_maxl"

    external get_num_lin_iters : session -> int
        = "c_spils_get_num_lin_iters"

    external get_num_conv_fails : session -> int
        = "c_spils_get_num_conv_fails"

    external get_num_prec_evals : session -> int
        = "c_spils_get_num_prec_evals"

    external get_num_prec_solves : session -> int
        = "c_spils_get_num_prec_solves"

    external get_num_jtimes_evals : session -> int
        = "c_spils_get_num_jtimes_evals"

    external get_num_rhs_evals : session -> int
        = "c_spils_get_num_rhs_evals"

  end

