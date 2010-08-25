
module Cvode = Cvode_serial

let f t y yd =
  yd.{0} <- cos(y.{1});
  yd.{1} <- 1.0;
  yd.{2} <- 0.0

let g t y gout =
  gout.{0} <- y.{0};
  gout.{1} <- y.{1}

let y = Cvode.Carray.of_array [| 0.0; 0.0; 0.0 |]

let s = Cvode.init Cvode.Adams Cvode.Functional f (2, g) y
let rootdata = Cvode.Roots.create 2

let _ =
  Cvode.Carray.print_with_time 0.0 y;
  (* for i = 1 to 200 do *)
  let t = ref 0.1 in
  while true do
    let (t', roots) = Cvode.advance s !t y in
        Cvode.Carray.print_with_time t' y;
        t := t' +. 0.1;
        if (roots) then begin
          Cvode.get_roots s rootdata;
          Cvode.Roots.print rootdata
        end
  done

let _ = Cvode.free s

