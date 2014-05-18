
module Cvode = Cvode_serial

let f t y yd =
  yd.{0} <- if y.{0} >= 0.0 then -1.0 else 1.0

let y = Cvode.Carray.of_array [| 1.0 |]

let s = Cvode.init Cvode.Adams Cvode.Functional
                   Cvode.default_tolerances f y

(* let _ = Cvode.set_stop_time s 10.0 *)

let _ =
  Cvode.Carray.print_with_time 0.0 y;
  (* for i = 1 to 200 do *)
  let t = ref 0.1 in
  let keep_going = ref true in
  while !keep_going do
    let (t', result) = Cvode.solve_normal s !t y in
        Cvode.Carray.print_with_time t' y;
        t := t' +. 0.1;
        match result with
        | Cvode.RootsFound -> failwith "There are no roots!"
        | Cvode.StopTimeReached -> keep_going := false
        | Cvode.Continue -> ();
  done

