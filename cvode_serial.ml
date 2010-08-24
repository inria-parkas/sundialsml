(* Aug 2010, Timothy Bourke (INRIA) *)

let kind = Bigarray.float64
let layout = Bigarray.c_layout
type c_array =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
type int_array =
  (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t

type val_array = c_array
type der_array = c_array
type root_array = c_array

let create = Bigarray.Array1.create kind layout
let of_array = Bigarray.Array1.of_array kind layout
let int_array = Bigarray.Array1.create Bigarray.int32 layout

let length = Bigarray.Array1.dim

let print_results t v =
  Printf.printf "%.8f" t;
  for i = 0 to (length v - 1) do
    Printf.printf "\t%f" v.{i}
  done;
  print_newline ()

let print_roots v =
  let found = ref false in
  for i = 0 to (length v - 1) do
    if (v.{i} <> 0l) then (Printf.printf " root-%03i" i; found := true)
  done;
  if (!found) then print_newline ()

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

external reinit : session -> float -> val_array -> unit
    = "c_reinit"

external set_tolerances : session -> float -> c_array -> unit
    = "c_set_tolerances"

external get_roots : session -> int_array -> unit
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

