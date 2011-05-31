
module Cvode = Cvode_serial

let root_dir  = ref Cvode.Increasing

let args = [
    ("-e", Arg.Unit (fun () -> root_dir := Cvode.IncreasingOrDecreasing),
     "Use either instead of up()");
]

let _ = Arg.parse args
    (fun _ -> ()) "sinezero: different types of root detection on a sine wave"

(* indices *)
let x = 0

let x_i = 0.0

let pi = 4.0 *. atan 1.0

let max_sim_t = 4.0 *. pi +. 1.0

let f t_s y yd =
  yd.{x} <- cos(t_s);
  Printf.printf "f(% e, [% e]) = [% e]\n" t_s y.{x} yd.{x}

let g t_s y gout =
  gout.{0} <- y.{x};
  Printf.printf "g(% e, [% e]) = [% e]\n" t_s y.{x} gout.{0}

(* simulation *)

let y = Cvode.Carray.of_array [| x_i |]

let s = Cvode.init Cvode.Adams Cvode.Functional f (1, g) y
let rootdata = Cvode.Roots.create 1

let _ = Cvode.set_all_root_directions s !root_dir
let _ = Cvode.set_stop_time s max_sim_t
let _ = Cvode.set_max_step s 5.0

exception Done

let _ =
  Printf.printf "time\t\t t\n";
  Printf.printf "------------------------------------\n";
  Cvode.Carray.print_with_time 0.0 y;
  try
    let i = ref 0 in
    while true do
      let (t', result) = Cvode.one_step s max_sim_t y in

      Printf.printf "\nstep %3d.\n" !i;
      incr i;

      Cvode.Carray.print_with_time t' y;
      Printf.printf "\t\t(step size = %e)\n" (Cvode.get_last_step s);
        
      match result with
      | Cvode.RootsFound -> print_endline "** root found"
      | Cvode.StopTimeReached -> raise Done
      | Cvode.Continue -> ()
    done
  with Done -> ()

