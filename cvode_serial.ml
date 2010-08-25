(* Aug 2010, Timothy Bourke (INRIA) *)

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
      let found = ref false in
      for i = 0 to (Carray.length v - 1) do
        if (isroot i) then (Printf.printf " root-%03i" i; found := true)
      done;
      if (!found) then print_newline ()

    let length = Bigarray.Array1.dim

    let reset v = Bigarray.Array1.fill v 0l
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

type session
exception RecoverableFailure
let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("cvode_RecoverableFailure",      RecoverableFailure);
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
  ]

external init' : lmm -> iter -> val_array -> int -> session
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

external advance : session -> float -> val_array -> float * bool
    = "c_advance"

external step : session -> float -> val_array -> float * bool
    = "c_step"

let init lmm iter f (num_roots, roots) y0 =
  Callback.register "cvode_serial_callback_f" f;
  Callback.register "cvode_serial_callback_roots" roots;
  init' lmm iter y0 num_roots

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
  in let _ = print_endline "--"
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

external set_error_handler' : session -> unit 
    = "c_set_error_handler"
let set_error_handler s errh =
  Callback.register "cvode_serial_callback_errh" errh;
  set_error_handler' s

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

