
module Cvode = Cvode.Serial

let f t y yd =
  yd.{0} <- 1.0

let g t y gout =
  gout.{0} <- if t < 2.000012345 then -1.0 else 1.0

let y = Cvode.Carray.of_array [| 0.0 |]

let s = Cvode.init Cvode.Adams Cvode.Functional f (1, g) y
let rootdata = Cvode.Roots.create 1

let _ = Cvode.set_stop_time s 5.0

let _ =
  Cvode.Carray.print_with_time 0.0 y;
  let t = ref 0.1 in
  let keep_going = ref true in
  while !keep_going do
    let (t', result) = Cvode.normal s !t y in
        Cvode.Carray.print_with_time t' y;
        t := t' +. 0.1;
        match result with
        | Cvode.RootsFound -> begin
              Cvode.get_root_info s rootdata;
              Cvode.Roots.print rootdata
            end
        | Cvode.StopTimeReached -> keep_going := false
        | Cvode.Continue -> ();
  done

let _ = Cvode.free s

