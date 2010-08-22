
let f t y yd =
  yd.{0} <- cos(y.{1});
  yd.{1} <- 1.0;
  yd.{2} <- 0.0

let g t y gout =
  gout.{0} <- y.{0};
  gout.{1} <- y.{1}

let y = Cvode_serial.of_array [| 0.0; 0.0; 0.0 |]

let s = Cvode_serial.init f (2, g) y
let rootdata = Cvode_serial.int_array 2

let _ =
  Cvode_serial.print_results 0.0 y;
  for i = 1 to 200 do
    let (t', roots) = Cvode_serial.advance s (0.1 *. float(i)) y in
        Cvode_serial.print_results t' y;
        if (roots) then begin
          Cvode_serial.get_roots s rootdata;
          Cvode_serial.print_roots rootdata
        end
  done

let _ = Cvode_serial.free s

