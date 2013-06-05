(* A pendulum attached to a rigid rod banging against a wall, modeled by a DAE
   in Cartesian coordinates.  The wall extends 30 degrees to the left of the
   plumb line extending from the pivot.

      .
     //\          // is the wall, \ is the rod, . is the pivot, o is a mass.
    //  \
   //    o

   The Cartesian-coordinate DAE governing this system is

     x'' = p x
     y'' = p y - g
     x^2 + y^2 = r^2

   where g is the gravitational acceleration, r is the length of the rod, x,y
   is the position of the mass, and p is the mass' (signed) acceleration due to
   the rod pulling on the mass.  The rod's pull is parallel to the direction of
   the rod, so it is a scalar multiple of (x,y), i.e. (p x, p y) where p is a
   scalar variable.

   IDA needs the system to be put into the form F(t,X,X') = 0, so we introduce
   auxiliary variables vx and vy, representing the mass' velocity.

     vx - x'         = 0
     vy - y'         = 0
     vx' - p x       = 0
     vy' - p y + g   = 0
     x^2 + y^2 - r^2 = 0

   This is a system with five unknowns (x,y,vx,vy,p) and five equations, so it
   is properly constrained (i.e. neither over- nor under-constrained).

   p is an algebraic variable (its derivative p' does not appear), while all
   other variables are differential (their derivatives appear in the
   equations).

 *)

(* Setup & auxiliary functions *)
module Ida = Ida_serial
let pi = 4. *. atan (1.)
let degree_to_radian x = x *. pi /. 180.

(* Problem constants *)
let r = 1.0                            (* length of rod [m] *)
let g = 9.8                            (* gravitational acceleration [m/s^2] *)
let k = 0.5                            (* elasticity of collision with wall *)
(* direction of the wall relative to the pivot *)
let wall_angle = degree_to_radian 30.
let wall = (-. sin wall_angle, -. cos wall_angle)

(* Initial conditions; because p is an algebraic variable, its initial values
   can be computed from the rest of the initial values using calc_ic_ya_yd'.
   (But note that the system is index 3, so the algorithm is not guaranteed
    to converge.  Indeed, poorly chosen initial guesses can result in
    convergence failure.)  *)
let theta0 = degree_to_radian 80.       (* angle with vertical axis *)
let (x0, y0) = (r *. sin theta0, -. r *. cos theta0)
let (x'0, y'0) = (0., 0.)

(* There are 5 variables.  *)
let x_i   = 0
let y_i   = 1
let vx_i  = 2
let vy_i  = 3
let p_i   = 4
let nvars = 5
(* There are 5 residues (i.e. 5 equations).  *)
let vx_x  = 0                           (* Relationship between vx and x. *)
let vy_y  = 1                           (* Relationship between vy and y. *)
let acc_x  = 2                          (* Acceleration in the x direction. *)
let acc_y  = 3                          (* Acceleration in the y direction. *)
let constr = 4                          (* Constraint. *)
let neqs   = 5
(* Number of residues should match the number of variables. *)
let _ = assert (neqs = nvars)

(* The residual function, F.  *)
let res_pendulum t vars vars' res =
  let x  = vars.{x_i}  and x'  = vars'.{x_i}
  and y  = vars.{y_i}  and y'  = vars'.{y_i}
  and vx = vars.{vx_i} and vx' = vars'.{vx_i}
  and vy = vars.{vy_i} and vy' = vars'.{vy_i}
  and p  = vars.{p_i} in
  res.{vx_x}   <- vx -. x';
  res.{vy_y}   <- vy -. y';
  res.{acc_x}  <- vx' -. p *. x;
  res.{acc_y}  <- vy' -. p *. y +. g;
  res.{constr} <- vx*.x +. vy*.y
    (* vx'*.x +. vx*.vx +. vy'*.y +. vy*.vy *) (* <-- works, index 1 *)
    (* vx*.x +. vy*.y *)                       (* <-- works, index 2 *)
    (* x*.x +. y*.y -. r*.r *)                 (* <-- convergence failures *)

(* The root function -- computes distance from the wall along x axis.  The wall
   is a line that connects pivot = (0,0) and wall = (wx,wy) where wall is a
   constant defined above.  *)
let distance_from_wall vars =
  let x  = vars.{x_i}
  and y  = vars.{y_i}
  and (wx,wy) = wall in
  y *. (wx /. wy) -. x
let roots t vars vars' r =
  r.{0} <- distance_from_wall vars

(* The initial value of p' can be any garbage -- it's not used.  The initial
   value of p should be close enough to its correct value so that an iterative
   method converges.  *)
let vars = Ida.Carray.of_array [|x0; y0; x'0; y'0; g*.y0 /. (r*.r)|]
let vars' = Ida.Carray.of_array [|x'0; y'0; g*.y0*.x0 /. (r*.r); g*.y0*.y0/.(r*.r) -.g; g*.y'0 /. (r*.r)|]
let var_types =
  let d = Ida.VarTypes.Differential and a = Ida.VarTypes.Algebraic in
  Ida.VarTypes.of_array [|d; d; d; d; a|]

(* Rendering parameters  *)

(* Position of pivot on the screen, scaled so that the screen is 1.0x1.0.  The
   solver will take the pivot as the origin. *)
let pivot = (0.5, 0.8)

(* Options *)
let trace = ref false
let log = ref false
let show = ref true
let delay = ref true
let dt = ref (1./.60.)                  (* time step [s] *)
let t_end = ref 10.0                    (* simulation end time [s] *)

let args = [
    ("-trace", Arg.Set trace, "Show a trace of mass positions.");
    ("-d", Arg.Set_float dt, "Set the default time step.");
    ("-noshow", Arg.Clear show, "Disable the graphical display.");
    ("-nodelay", Arg.Clear delay, "No delays between frames.");
    ("-log",
     Arg.Unit (fun () -> log := true; show := false; delay := false),
     "Log state variables to stdout (implies -noshow and -nodelay).");
]

let show_nvector a =
  let n = Bigarray.Array1.dim a in
  let str = ref "[" in
  for i = 0 to n-2 do
    str := !str ^ Printf.sprintf "%g," a.{i}
  done;
  if 0 <= n-1 then str := !str ^ Printf.sprintf "%g" a.{n-1};
  str := !str ^ "]";
  !str
let print_nvector a = Printf.printf "%s" (show_nvector a)

let _ =
  Arg.parse args (fun _ -> ()) "pendulum: simulate a pendulum hitting against an oblique wall";
  if !log then Ida.Carray.print_with_time 0.0 vars;
  let t_delay = if !delay then !dt else 0. in
  if !show then Showpendulum.start (1.5 *. r) t_delay !trace pivot wall;
  let frames = int_of_float (!t_end /. !dt) in

  let ida = Ida.init Ida.Dense res_pendulum (1, roots) vars vars' in
  Ida.calc_ic_ya_yd' ~y:vars ~y':vars' ida var_types !dt;
  Ida.set_suppress_alg ida true;

  let t = ref 0.
  and tnext = ref !dt
  in
  for i = 1 to frames do
    tnext := !t +. !dt;
    while !t < !tnext do
      let (tret, flag) = Ida.solve_normal ida !tnext vars vars' in
      t := tret;
      if !show then Showpendulum.show (vars.{x_i}, vars.{y_i});
      if !log then Ida.Carray.print_with_time tret vars;

      (* IDA reports spurious zero-crossing right after a reset (probably due
         to floating point error), so we have to guard the reset with a
         heuristic to judge if this is a genuine zero-crossing.  *)
      if flag = Ida.RootsFound && vars.{vx_i} < 0. then
        (Printf.printf "Bang!  Hit against the wall.\n"; flush stdout;
         vars.{vx_i} <- -. k *. vars.{vx_i};
         vars.{vy_i} <- -. k *. vars.{vy_i};
         vars'.{x_i} <- -. k *. vars'.{x_i};
         vars'.{y_i} <- -. k *. vars'.{y_i};
         Ida.reinit ida !t vars vars'
        )
    done
  done;
  Showpendulum.stop ()
