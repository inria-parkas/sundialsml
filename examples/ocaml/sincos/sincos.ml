
let print_with_time t v =
  Printf.printf "%e" t;
  Sundials.RealArray.iter (Printf.printf "\t% e") v;
  print_newline ()

let string_of_root x =
  match x with
  | Sundials.Roots.Rising  -> "rising"
  | Sundials.Roots.Falling -> "falling"
  | Sundials.Roots.NoRoot  -> "noroot"

let print_roots vs =
  Sundials.Roots.iter (fun x ->
    Printf.printf "\t%s" (string_of_root x)) vs;
  print_newline ()

let f t y yd =
  yd.{0} <- cos(y.{1});
  yd.{1} <- 1.0;
  yd.{2} <- 0.0

let g t y gout =
  gout.{0} <- y.{0};
  gout.{1} <- y.{1}

(* let f t y yd =
  yd.{0} <- if t <= 10.0 then -1.0 else 1.0 *)

let y = Sundials.RealArray.of_array [| 0.0; 0.0; 0.0 |]
let y_nvec= Nvector_serial.wrap y

let s = Cvode.init Cvode.Adams Cvode.Functional Cvode.default_tolerances
                   f 0.0 ~roots:(2, g) y_nvec
let rootdata = Sundials.Roots.create 2

(* let _ = Cvode.set_stop_time s 20.0 *)

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
        | Cvode.RootsFound -> begin
              Cvode.get_root_info s rootdata;
              Printf.printf "R: %e" t';
              print_roots rootdata
            end
        | Cvode.StopTimeReached -> keep_going := false
        | Cvode.Success -> ();
  done

