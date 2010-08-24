(* Aug 2010, Timothy Bourke (INRIA) *)

let kind = Bigarray.float64
let layout = Bigarray.c_layout
type c_array =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
let empty = Bigarray.Array1.create kind layout 0

type val_array = c_array
type der_array = c_array
type rootval_array = c_array

let create = Bigarray.Array1.create kind layout
let of_array = Bigarray.Array1.of_array kind layout

let fill = Bigarray.Array1.fill

let length = Bigarray.Array1.dim

let print_results t v =
  Printf.printf "%.8f" t;
  for i = 0 to (length v - 1) do
    Printf.printf "\t%f" v.{i}
  done;
  print_newline ()

(* root arrays *)

module Roots =
  struct
    type t = (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t

    let create = Bigarray.Array1.create Bigarray.int32 layout
    let empty = create 0

    let get roots i = roots.{i} <> 0l

    let set a i v = Bigarray.Array1.set a i (if v then 1l else 0l)

    let print v =
      let isroot = get v in
      let found = ref false in
      for i = 0 to (length v - 1) do
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

