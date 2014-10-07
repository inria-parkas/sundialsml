
let t = 0
let x = 1

let thresh1 = 0.987
let thresh2 = 1.987

let t_i = 0.0
let x_i = 0.0

let max_sim_t = 5.0

let print_with_time t v =
  Printf.printf "%e" t;
  Sundials.RealArray.iter (Printf.printf "\t% e") v;
  print_newline ()

(*
       der t = 1.0 init t_i
   and r = if t < thresh1 then 1.0
           else if t < thresh2 then -. 1.0
           else 2.0
   and der x = r init x_i

 *)

let f t_s y yd =
  let r =
    if y.{t} < thresh1 then 1.0
    else if y.{t} < thresh2 then -1.0
    else 2.0
  in
  yd.{t} <- 1.0;
  yd.{x} <- r

(* z = up(t -. thresh2) *)

let g t_s y gout =
  gout.{0} <- y.{t} -. thresh2

(* simulation *)

let y = Sundials.RealArray.of_array [| t_i; x_i |]
let y_nvec= Nvector_serial.wrap y

let s = Cvode.init Cvode.Adams Cvode.Functional Cvode.default_tolerances
                   f ~roots:(1, g) 0. y_nvec
let rootdata = Sundials.Roots.create 1

let _ = Cvode.set_stop_time s max_sim_t

exception Done

let _ =
  Printf.printf "time\t\t t\t\t x\n";
  Printf.printf "----------------------------------------\n";
  print_with_time 0.0 y;
  try
    let i = ref 0 in
    while true do
      let (t', result) = Cvode.solve_one_step s max_sim_t y_nvec in

      Printf.printf "\nstep %3d.\n" !i;
      incr i;

      print_with_time t' y;
      Printf.printf "%e\tstep size = %e\n" t' (Cvode.get_last_step s);
        
      match result with
      | Sundials.RootsFound -> print_endline "** root found"
      | Sundials.StopTimeReached -> raise Done
      | Sundials.Continue -> ()
    done
  with Done -> ()

