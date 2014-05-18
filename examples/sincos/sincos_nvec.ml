
module Cvode = Cvode_nvector

let f t y yd =
  yd.(0) <- cos(y.(1));
  yd.(1) <- 1.0;
  yd.(2) <- 0.0

let g t y gout =
  gout.{0} <- y.(0);
  gout.{1} <- y.(1)

let y = [| 0.0; 0.0; 0.0 |]
let y_nv = Nvector_array.wrap y

let s = Cvode.init Cvode.Adams Cvode.Functional Cvode.default_tolerances
                   f ~roots:(2, g) (2, y_nv)
let rootdata = Cvode.Roots.create 2

(* let _ = Cvode.set_stop_time s 10.0 *)

let print_with_time t v =
  Cvode.print_time ("", "") t;
  Array.iter (Printf.printf "\t% .8f") v;
  print_newline ()

let _ =
  print_with_time 0.0 y;
  (* for i = 1 to 200 do *)
  let t = ref 0.1 in
  let keep_going = ref true in
  while !keep_going do
    let (t', result) = Cvode.solve_normal s !t y_nv in
        print_with_time t' y;
        t := t' +. 0.1;
        match result with
        | Cvode.RootsFound -> begin
              Cvode.get_root_info s rootdata;
              Cvode.print_time ("R: ", "") t';
              Cvode.Roots.print rootdata
            end
        | Cvode.StopTimeReached -> keep_going := false
        | Cvode.Continue -> ();
  done

