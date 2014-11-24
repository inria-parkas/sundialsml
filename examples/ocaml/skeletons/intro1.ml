(* Changes to this file must be duplicated in intro.doc. *)
#use "topfind";;
#require "sundialsml";;

let f t y yd = yd.{0} <- 1.;;
let g t y gout = gout.{0} <- y.{0};;
let y = Sundials.RealArray.of_array [| -1.0 |];;
let yvec = Nvector_serial.wrap y;;
let s = Cvode.init Cvode.Adams Cvode.Functional Cvode.default_tolerances
                   f ~roots:(1, g) 0. yvec;;
Cvode.set_stop_time s 2.;;

(* repeat the commands below to advance the simulation until t = 2.0 *)
let (t', result) = Cvode.solve_normal s 2. yvec;;
Printf.printf "%e: " t';
Sundials.RealArray.iter (Printf.printf "\t%e") y;
Printf.printf "\n";;
