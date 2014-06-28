(* on essaie d'integrer de t+h a t *)

(* indices *)
let x = 0

let x_i = -1.0

let max_sim_t = 5.0

let f t_s y yd =
  yd.{x} <- 1.0;
  Printf.printf "f(% e, [% e]) = [% e]\n" t_s y.{x} yd.{x}

let g t_s y gout =
  gout.{0} <- 1.0 /. 0.0;
  if y.{x} <> 2.0 then Printf.printf "toti %.15e " y.{x}; gout.{0} <- 1.0 /. (2.0 -. y.{x});
  Printf.printf "g(% e, [% e]) = [% e]\n" t_s y.{x} gout.{0}

(* simulation *)

let y = Sundials.RealArray.of_array [| x_i |]
let y_nvec= Nvector_serial.wrap y

let s = Cvode.init Cvode.Adams Cvode.Functional Cvode.default_tolerances
                   f ~roots:(1, g) y_nvec
let rootdata = Sundials.Roots.create 1

let _ = Cvode.set_stop_time s max_sim_t

exception Done

let _ =
  Printf.printf "time\t\t t\n";
  Printf.printf "------------------------------------\n";
  Sundials.RealArray.print_with_time 0.0 y;
  try
    let i = ref 0 in
    while true do
      let (t', result) = Cvode.solve_one_step s max_sim_t y_nvec in

      Printf.printf "\nstep %3d.\n" !i;
      incr i;

      Sundials.RealArray.print_with_time t' y;
      Printf.printf "\t\t(step size = %e)\n" (Cvode.get_last_step s);
        
      match result with
      | Sundials.RootsFound -> print_endline "** root found"
      | Sundials.StopTimeReached -> raise Done
      | Sundials.Continue -> ()
    done
  with Done -> ()

