
module Cvode = Cvode.Serial

let t = 0
let x = 1

let thresh1 = 0.987
let thresh2 = 1.987

let t_i = 0.0
let x_i = 0.0

let max_sim_t = 5.0

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

let y = Cvode.Carray.of_array [| t_i; x_i |]

let s = Cvode.init Cvode.Adams Cvode.Functional f (1, g) y
let rootdata = Cvode.Roots.create 1

let _ = Cvode.set_stop_time s max_sim_t

exception Done

let _ =
  Printf.printf "time\t\t t\t\t x\n";
  Printf.printf "----------------------------------------\n";
  Cvode.Carray.print_with_time 0.0 y;
  try
    let i = ref 0 in
    while true do
      let (t', result) = Cvode.one_step s max_sim_t y in

      Printf.printf "\nstep %3d.\n" !i;
      incr i;

      Cvode.Carray.print_with_time t' y;
      Printf.printf "%e\tstep size = %e\n" t' (Cvode.get_last_step s);
        
      match result with
      | Cvode.RootsFound -> print_endline "** root found"
      | Cvode.StopTimeReached -> raise Done
      | Cvode.Continue -> ()
    done
  with Done -> ()

let _ = Cvode.free s

