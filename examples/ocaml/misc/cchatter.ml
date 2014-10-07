
let print_with_time t v =
  Printf.printf "%e" t;
  Sundials.RealArray.iter (Printf.printf "\t% e") v;
  print_newline ()

let f t y yd =
  yd.{0} <- if y.{0} >= 0.0 then -1.0 else 1.0

let y = Sundials.RealArray.of_array [| 1.0 |]
let y_nvec = Nvector_serial.wrap y

let s = Cvode.init Cvode.Adams Cvode.Functional
                   Cvode.default_tolerances f 0. y_nvec

(* let _ = Cvode.set_stop_time s 10.0 *)

let _ =
  print_with_time 0.0 y;
  (* for i = 1 to 200 do *)
  let t = ref 0.1 in
  let keep_going = ref true in
  while !keep_going do
    let (t', result) = Cvode.solve_normal s !t y_nvec in
        print_with_time t' y;
        t := t' +. 0.1;
        match result with
        | Sundials.RootsFound -> failwith "There are no roots!"
        | Sundials.StopTimeReached -> keep_going := false
        | Sundials.Continue -> ();
  done

