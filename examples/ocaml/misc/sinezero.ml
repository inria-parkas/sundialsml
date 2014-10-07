
let root_dir  = ref Sundials.RootDirs.Increasing

let args = [
    ("-e", Arg.Unit (fun () ->
             root_dir := Sundials.RootDirs.IncreasingOrDecreasing),
     "Use either instead of up()");
]

let _ = Arg.parse args
    (fun _ -> ()) "sinezero: different types of root detection on a sine wave"

(* indices *)
let x = 0

let x_i = 0.0

let pi = 4.0 *. atan 1.0

let max_sim_t = 4.0 *. pi +. 1.0

let print_with_time t v =
  Printf.printf "%e" t;
  Sundials.RealArray.iter (Printf.printf "\t% e") v;
  print_newline ()

let f t_s y yd =
  yd.{x} <- cos(t_s);
  Printf.printf "f(% e, [% e]) = [% e]\n" t_s y.{x} yd.{x}

let g t_s y gout =
  gout.{0} <- y.{x};
  Printf.printf "g(% e, [% e]) = [% e]\n" t_s y.{x} gout.{0}

(* simulation *)

let y = Sundials.RealArray.of_array [| x_i |]
let y_nvec = Nvector_serial.wrap y

let s = Cvode.init Cvode.Adams Cvode.Functional Cvode.default_tolerances
                   f ~roots:(1, g) 0. y_nvec
let rootdata = Sundials.Roots.create 1

let _ = Cvode.set_all_root_directions s !root_dir
let _ = Cvode.set_stop_time s max_sim_t
let _ = Cvode.set_max_step s 5.0

exception Done

let _ =
  Printf.printf "time\t\t t\n";
  Printf.printf "------------------------------------\n";
  print_with_time 0.0 y;
  try
    let i = ref 0 in
    while true do
      let (t', result) = Cvode.solve_one_step s max_sim_t y_nvec in

      Printf.printf "\nstep %3d.\n" !i;
      incr i;

      print_with_time t' y;
      Printf.printf "\t\t(step size = %e)\n" (Cvode.get_last_step s);
        
      match result with
      | Sundials.RootsFound -> print_endline "** root found"
      | Sundials.StopTimeReached -> raise Done
      | Sundials.Continue -> ()
    done
  with Done -> ()

