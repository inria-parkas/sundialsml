
let p = 1.0 (* 1.0 or 0.2 *)
let step_inc = 0.1 (* not -1 zero-crossing if set to 100.0 *)

(* indices *)
let x = 0
and t = 1

let x_i = -1.0
and t_i = -. p

let yc = ref 1.0

let print_with_time t v =
  Printf.printf "%.15e" t;
  Sundials.RealArray.iter (Printf.printf "\t% .15e") v;
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

let f _ _ yd =
  yd.{x} <- !yc;
  yd.{t} <- 1.0

let num_roots = 3

let z_x = 0
and z_nx = 1
and z_t  = 2

let g _ y gout =
  gout.{z_x}  <- y.{x};
  gout.{z_nx} <- -. y.{x};
  gout.{z_t}  <- y.{t}
  (* ; Printf.printf "->%.15f (z_x=%.15f z_nx=%.15f z_t=%.15f)\n" t_s gout.{z_x}
     gout.{z_nx} gout.{z_t} *)

let d _ r ys =
  let root = Sundials.Roots.detected r in

  if root z_x then yc := -1.0
  else if root z_nx then yc := 1.0
  else ();

  if root z_t then begin
    ys.{t} <- -. p;
    print_endline "** timer **"
  end

(* simulation *)

let y = Sundials.RealArray.of_array [| x_i; t_i |]
let y_nvec = Nvector_serial.wrap y

let s = Cvode.(init Adams default_tolerances f ~roots:(num_roots, g) 0. y_nvec)
let rootdata = Sundials.Roots.create num_roots
let _ = Cvode.set_all_root_directions s Sundials.RootDirs.Increasing

exception Done

let _ =
  Printf.printf "---period=%f\n" p;
  Printf.printf "time\t\t t\n";
  Printf.printf "------------------------------------\n";
  print_with_time 0.0 y;
  try
    let t = ref 0.0 in
    while true do
      let (t', result) = Cvode.solve_normal s (!t +. step_inc) y_nvec in
      (* let (t', result) = Cvode.one_step s 100.0 y in *)
      t := t';

      print_with_time t' y;
        
      match result with
      | Cvode.RootsFound -> begin
            Cvode.get_root_info s rootdata;
            print_roots rootdata;
            d t' rootdata y;
            print_endline "after discrete step:";
            print_with_time t' y;
            Cvode.reinit s t' y_nvec
          end
      | Cvode.StopTimeReached -> raise Done
      | Cvode.Success -> ()
    done
  with Done -> ()

