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
      Printf.printf "%.8f" t;
      for i = 0 to (length v - 1) do
        Printf.printf "\t%f" v.{i}
      done;
      print_newline ()
  end

type val_array = Carray.t
type der_array = Carray.t
type rootval_array = Carray.t

(* root arrays *)

type int_array = (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
let create_int_array = Bigarray.Array1.create Bigarray.int32 Carray.layout

module Roots =
  struct
    type t = int_array

    let create = create_int_array
    let empty = create 0

    let get roots i = roots.{i} <> 0l

    let set a i v = Bigarray.Array1.set a i (if v then 1l else 0l)

    let print v =
      let isroot = get v in
      for i = 0 to (Carray.length v - 1) do
        if i > 0 then print_string " "; 
        print_string (if (isroot i) then "1" else "0")
      done;
      print_newline ()

    let length = Bigarray.Array1.dim
    let reset v = Bigarray.Array1.fill v 0l

    let fold_left f a v =
      let rec check (i, a) =
        if i < 0 then a
        else check (i - 1, f a v.{i})
      in
      check (Bigarray.Array1.dim v - 1, a)

    let exists = fold_left (fun a x -> a || x <> 0l) false
  end

let no_roots = (0, (fun _ _ _ -> 0))

(* TODO: get these types working in the initialization function *)
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
| JacFn
| BandJacFn
| PreSetupFn
| PreSolveFn
| JacTimesFn

let handler_name h = match h with
  | RhsFn        -> "cvode_serial_callback_rhsfn"
  | RootsFn      -> "cvode_serial_callback_rootsfn"
  | ErrorHandler -> "cvode_serial_callback_errorhandler"
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
    = "c_reinit"

external set_tolerances : session -> float -> Carray.t -> unit
    = "c_set_tolerances"

external get_roots : session -> Roots.t -> unit
    = "c_get_roots"

external free : session -> unit
    = "c_free"

external advance : session -> float -> val_array -> float * solver_result
    = "c_advance"

external step : session -> float -> val_array -> float * solver_result
    = "c_step"

external get_dky : session -> float -> int -> Carray.t -> unit
    = "c_get_dky"

let init' lmm iter f (num_roots, roots) y0 t0 =
  let s = external_init lmm iter y0 num_roots t0 in
  register_handler s RhsFn f;
  register_handler s RootsFn roots;
  s

let init lmm iter f roots y0 = init' lmm iter f roots y0 0.0

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

external integrator_stats : session -> integrator_stats
    = "c_integrator_stats"

external last_step_size : session -> float
    = "c_last_step_size"

external next_step_size : session -> float
    = "c_next_step_size"

let print_integrator_stats s =
  let stats = integrator_stats s
  in
    Printf.printf "steps = %d\n"                stats.steps;
    Printf.printf "rhs_evals = %d\n"            stats.rhs_evals;
    Printf.printf "linear_solver_setups = %d\n" stats.linear_solver_setups;
    Printf.printf "error_test_failures = %d\n"  stats.error_test_failures;
    Printf.printf "last_internal_order = %d\n"  stats.last_internal_order;
    Printf.printf "next_internal_order = %d\n"  stats.next_internal_order;
    Printf.printf "initial_step_size = %e\n"    stats.initial_step_size;
    Printf.printf "last_step_size = %e\n"       stats.last_step_size;
    Printf.printf "next_step_size = %e\n"       stats.next_step_size;
    Printf.printf "internal_time = %e\n"        stats.internal_time;

external set_error_file : session -> string -> bool -> unit 
    = "c_set_error_file"

external enable_error_handler : session -> unit 
    = "c_enable_error_handler"

let set_error_handler s errh =
  register_handler s ErrorHandler errh;
  enable_error_handler s

external set_max_ord : session -> int -> unit 
    = "c_set_max_ord"
external set_max_num_steps : session -> int -> unit 
    = "c_set_max_num_steps"
external set_max_hnil_warns : session -> int -> unit 
    = "c_set_max_hnil_warns"
external set_stability_limit_detection : session -> bool -> unit 
    = "c_set_stability_limit_detection"
external set_initial_step_size : session -> float -> unit 
    = "c_set_initial_step_size"
external set_min_abs_step_size : session -> float -> unit 
    = "c_set_min_abs_step_size"
external set_max_abs_step_size : session -> float -> unit 
    = "c_set_max_abs_step_size"
external set_stop_time : session -> float -> unit 
    = "c_set_stop_time"
external set_max_error_test_failures : session -> int -> unit 
    = "c_set_max_error_test_failures"
external set_max_nonlinear_iterations : session -> int -> unit 
    = "c_set_max_nonlinear_iterations"
external set_max_convergence_failures : session -> int -> unit 
    = "c_set_max_convergence_failures"
external set_nonlinear_convergence_coeffficient : session -> float -> unit 
    = "c_set_nonlinear_convergence_coeffficient"
external set_nonlinear_iteration_type : session -> iter -> unit 
    = "c_set_nonlinear_iteration_type"

external set_root_direction' : session -> int_array -> unit 
    = "c_set_root_direction"

let int32_of_root_direction x =
  match x with
  | Increasing -> 1l
  | Decreasing -> -1l
  | IncreasingOrDecreasing -> 0l
    
let set_root_direction s rda =
  let n = nroots s in
  let rdirs = create_int_array n in
  if (n > Array.length rda)
    then Bigarray.Array1.fill rdirs
            (int32_of_root_direction IncreasingOrDecreasing);
  Array.iteri (fun i v -> rdirs.{i} <- int32_of_root_direction v) rda;
  set_root_direction' s rdirs

let set_all_root_directions s rd =
  let rdirs = create_int_array (nroots s) in
  Bigarray.Array1.fill rdirs (int32_of_root_direction rd);
  set_root_direction' s rdirs

external disable_inactive_root_warnings : session -> unit 
    = "c_disable_inactive_root_warnings"

(*
let print_stats s =
  let in_stats = integrator_stats s
  (* and ls_stats = TODO *)
  (* and jac_evals = TODO CVDlsGetNumJacEvals *)
  (* and root_evals = TODO CVodeGetNumGEvals *)
  and printf = Printf.printf
  in
    printf("\nFinal Statistics:\n");
    printf "nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n"
      in_stats.steps
      in_stats.rhs_evals
      in_stats.linear_solver_setups
      ls_stats.rhs_evals
      jac_evals;
    printf "nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n\n"
      ls_stats.iterations
      ls_stats.convergence_failures
      in_stats.error_test_failures
      root_evals
*) 

(* direct linear solvers optional input functions *)

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

external enable_dense_jacobian_fn : session -> unit
    = "c_enable_dense_jacobian_fn"

let set_dense_jacobian_fn s f =
    register_handler s JacFn f;
    enable_dense_jacobian_fn s

external enable_band_jacobian_fn : session -> unit
    = "c_enable_band_jacobian_fn"

let set_band_jacobian_fn s f =
    register_handler s BandJacFn f;
    enable_band_jacobian_fn s

(* iterative linear solvers optional input functions *)
module Spils =
  struct
    type solve_arg =
      {
        rhs   : val_array;
        gamma : float;
        delta : float;
        left  : bool; (* true: left, false: right *)
      }

    type single_tmp = val_array

    type preconditioning_type =
    | PrecNone
    | PrecLeft
    | PrecRight
    | PrecBoth

    type gramschmidt_type =
    | ModifiedGS
    | ClassicalGS

    external enable_preconditioner_fns : session -> unit
        = "c_enable_preconditioner_fns"

    let set_preconditioner_fns s fsetup fsolve =
        register_handler s PreSetupFn fsetup;
        register_handler s PreSolveFn fsolve;
        enable_preconditioner_fns s

    external enable_jacobian_times_vector_fn : session -> unit
        = "c_enable_jacobian_times_vector_fn"

    let set_jacobian_times_vector_fn s f =
        register_handler s JacTimesFn f;
        enable_jacobian_times_vector_fn s

    external set_preconditioning_type : session -> preconditioning_type -> unit
        = "c_set_preconditioning_type"

    external set_gramschmidt_orthogonalization :
        session -> gramschmidt_type -> unit
        = "c_set_gramschmidt_orthogonalization"

    external set_eps_linear_convergence_factor : session -> float -> unit
        = "c_set_eps_linear_convergence_factor"

    external set_max_subspace_dimension : session -> int -> unit
        = "c_set_max_subspace_dimension"

  end

