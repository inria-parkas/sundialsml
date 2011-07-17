
module Cvode = Cvode_serial

let p = 1.0 (* 1.0 or 0.2 *)
let step_inc = 0.1 (* not -1 zero-crossing if set to 100.0 *)

(* indices *)
let x = 0
and t = 1

let x_i = -1.0
and t_i = -. p

let yc = ref 1.0

let f t_s y yd =
  yd.{x} <- !yc;
  yd.{t} <- 1.0

let num_roots = 3

let z_x = 0
and z_nx = 1
and z_t  = 2

let g t_s y gout =
  gout.{z_x}  <- y.{x};
  gout.{z_nx} <- -. y.{x};
  gout.{z_t}  <- y.{t}
  (* ; Printf.printf "->%.15f (z_x=%.15f z_nx=%.15f z_t=%.15f)\n" t_s gout.{z_x}
     gout.{z_nx} gout.{z_t} *)

let d t_s r ys =
  let root = Cvode.Roots.get r in

  if root z_x then yc := -1.0
  else if root z_nx then yc := 1.0
  else ();

  if root z_t then begin
    ys.{t} <- -. p;
    print_endline "** timer **"
  end

(* simulation *)

let y = Cvode.Carray.of_array [| x_i; t_i |]

let s = Cvode.init Cvode.Adams Cvode.Functional f (num_roots, g) y
let rootdata = Cvode.Roots.create num_roots
let _ = Cvode.set_all_root_directions s Cvode.Increasing
let _ = Cvode.extra_time_precision := true

exception Done

let _ =
  Printf.printf "---period=%f\n" p;
  Printf.printf "time\t\t t\n";
  Printf.printf "------------------------------------\n";
  Cvode.Carray.print_with_time'' 0.0 y;
  try
    let t = ref 0.0 in
    while true do
      let (t', result) = Cvode.normal s (!t +. step_inc) y in
      (* let (t', result) = Cvode.one_step s 100.0 y in *)
      t := t';

      Cvode.Carray.print_with_time'' t' y;
        
      match result with
      | Cvode.RootsFound -> begin
            Cvode.get_root_info s rootdata;
            Cvode.Roots.print rootdata;
            d t' rootdata y;
            print_endline "after discrete step:";
            Cvode.Carray.print_with_time'' t' y;
            Cvode.reinit s t' y
          end
      | Cvode.StopTimeReached -> raise Done
      | Cvode.Continue -> ()
    done
  with Done -> ()

