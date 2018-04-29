(* A pendulum attached to a rigid rod banging against a wall, modeled by a DAE
   in Cartesian coordinates.  The wall extends 30 degrees to the left of the
   plumb line extending from the pivot.

      .
     //\          // is the wall, \ is the rod, . is the pivot, o is a mass.
    //  \
   //    o

 *)

module RealArray = Sundials.RealArray
module Nvector = Nvector_serial

let (+) = (+.)
let (-) = (-.)
let ( * ) = ( *. )
let (/) = (/.)
let neg x = -. x

let dt = (1. / 60.)
let t_end = 10.0

let pi = 4. * atan (1.)

let r = 1.0    (* length of rod [m] *)
let g = 9.8    (* gravitational acceleration [m/s^2] *)
let k = -0.5   (* elasticity of collision with wall *)

(* direction of the wall relative to the pivot *)
let wall_angle = -. pi / 6.
let wall = (sin wall_angle, -. cos wall_angle)

(* Position of pivot on the screen, scaled so that the screen is 1.0x1.0.  The
   solver will take the pivot as the origin. *)
let pivot = (0.5, 0.8)

(* There are 2 variables.  *)
let theta, theta' = 0, 1

let rhs t y yd =
  yd.{theta}  <- y.{theta'};
  yd.{theta'} <- -. g * sin y.{theta}

let roots t y r =
  r.{0} <- wall_angle - y.{theta}

let main () =
  let y = RealArray.of_list [ pi/2. ; 0. ] in
  let nv_y = Nvector.(wrap y) in

  let s = Cvode.(init Adams Functional default_tolerances
                      rhs ~roots:(1, roots) 0.0 nv_y)
  in
  Cvode.set_stop_time s 10.0;
  Cvode.set_all_root_directions s Sundials.RootDirs.Increasing;

  let rec stepto tnext t =
    if t >= tnext then t else
    match Cvode.solve_normal s tnext nv_y with
    | (tret, Cvode.RootsFound) ->
        y.{theta'} <- k * y.{theta'};
        Cvode.reinit s tret nv_y;
        stepto tnext tret
    | (tret, _) -> tret
  in
  let rec showloop t = if t < t_end then begin
    Showpendulum.show (r * sin y.{theta}, -. r * cos y.{theta});
    showloop (stepto (t + dt) t)
  end in
  showloop 0.0

let _ =
  try
    Showpendulum.start (1.5 *. r) dt false pivot wall;
    main ();
    Showpendulum.stop ()
  with Graphics.Graphic_failure _ -> ()

