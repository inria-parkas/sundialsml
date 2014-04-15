(* Compile with:
    ocamlc -o ramp.byte -I +sundials -dllpath +sundials \
           bigarray.cma unix.cma sundials.cma ramp.ml
 *)

module Cvode = Cvode_serial

type mode = RampingUp | Flat | RampingDown

let string_of_mode m =
  match m with
  | RampingUp -> "RampingUp"
  | Flat -> "Flat"
  | RampingDown -> "RampingDown"

let debug = ref false
let printf = Printf.printf

let t_i = 0.0
let y_i = 0.0

let idx_t = 0
let idx_y = 1

let thresh_1 = 3.0
let thresh_2 = 7.0

let max_sim_t = 10.0

let disc_state = ref RampingUp
let f1 t_s y yd =
  if !debug then
    printf "> f called at %.15e (%s)\n" t_s (string_of_mode !disc_state);
  yd.{idx_t} <- 1.0;
  yd.{idx_y} <-
    match !disc_state with
    | RampingUp -> 1.0
    | Flat -> 0.0
    | RampingDown -> -1.0

let f2 t_s y yd =
  if !debug then
    printf "> f called at %.15e (%s)\n" t_s (string_of_mode !disc_state);
  yd.{idx_t} <- 1.0;
  yd.{idx_y} <-
    if t_s <= 3.0 then 1.0
    else if t_s <= 7.0 then 0.0
    else -1.0

let g t_s y gout =
  gout.{0} <- y.{idx_t} -. thresh_1; (* up(t -. thresh_1) *)
  gout.{1} <- y.{idx_t} -. thresh_2  (* up(t -. thresh_2) *)

let handle_roots r =
  if !debug then printf "handle_roots: ";
  if Cvode.Roots.detected r 0
  then (if !debug then printf "up0"; disc_state := Flat);
  if Cvode.Roots.detected r 1
  then (if !debug then printf "up1"; disc_state := RampingDown);
  if !debug then printf "\n"

(* simulation *)

let rootdata = Cvode.Roots.create 2
exception Done

let run_experiment with_zeros with_reinit =
  disc_state := RampingUp;
  let y = Cvode.Carray.of_array [| t_i; y_i |] in

  let gg = if with_zeros then (2, g) else Cvode.no_roots in
  let ff = if with_zeros then f1 else f2 in

  let s = Cvode.init Cvode.Adams Cvode.Functional ff ~roots:gg y in
  let _ = Cvode.set_stop_time s max_sim_t in

  printf "t_sim\t\t\tt\t\t\ty\t\t\ty (ideal)\t\tattempted step\t\tactual step\t\terr test fails\n";
  printf "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%d\n" 0.0 y.{idx_t} y.{idx_y} 0.0 0.0 0.0 0;
  try
    let pre_err_test_fails = ref 0 in
    while true do
      let curr_step_size = Cvode.get_current_step s in
      let (t', result) = Cvode.solve_one_step s max_sim_t y in
      let last_step_size = Cvode.get_last_step s in

      let y_ideal =
        if t' <= thresh_1 then t'
        else if t' <= thresh_2 then 3.0
        else 10.0 -. t' in

      let err_test_fails = Cvode.get_num_err_test_fails s in
      printf "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%d\n"
        t' y.{idx_t} y.{idx_y} y_ideal curr_step_size last_step_size
        (err_test_fails - !pre_err_test_fails);
      pre_err_test_fails := err_test_fails;

      match result with
      | Cvode.RootsFound ->
          Cvode.get_root_info s rootdata;
          handle_roots rootdata;
          if with_reinit then Cvode.reinit s t' y
      | Cvode.StopTimeReached -> raise Done
      | Cvode.Continue -> ()
    done
  with Done -> ()

let opt_with_zeros = ref false
let opt_with_reinit = ref true

let args = [
    ("-zeros",
     Arg.Set opt_with_zeros,
     "Use a discrete state and zero-crossings.");
    ("-no-reinit",
     Arg.Clear opt_with_reinit,
     "Do not reinitialize the solver after zero-crossings.");
    ("-debug",
     Arg.Set debug,
     "Print internal debugging messages.");
  ]

let _ =
  Arg.parse args (fun _ -> ()) "ramp: a simple test of solver behaviour.";
  run_experiment !opt_with_zeros !opt_with_reinit

