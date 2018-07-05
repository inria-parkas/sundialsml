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

   The system as given above has index 3.  By differentiating the last
   constraint we can lower the index.

     index 2:
         vx - x'         = 0
         vy - y'         = 0
         vx' - p x       = 0
         vy' - p y + g   = 0
         x vx + y vy     = 0

     index 1:
         vx - x'         = 0
         vy - y'         = 0
         vx' - p x       = 0
         vy' - p y + g   = 0
         x vx' + y vy' + vx^2 + vy^2 = 0

   Further variations can be obtained by replacing some occurrences of vx in
   the last constraint equation with x' (but not occurrences of vx' with x'',
   because that would make the system second-order again) and replacing vy by
   y'.  Conceptually, one would expect that this has no effect on the solution,
   but it seems to affect the precision of the solution and some formulations
   can cause convergence failures.

   There are total of 14 different formulations for the last constraint, which
   can be switched from the command line.  A constructor of the constraint_form
   type showing the desired formula should be given varbatim as the argument to
   the -constraint commandline option.

   This file supplies analytically derived Jacobian function and initial
   condition calculator.  The analytical Jacobian function is more precise than
   the fallback implementation provided by IDA, while the initial condition
   calculator is precise and never fails whereas IDA's initial condition
   calculator fails depending on the format of the constraint formula.  These
   functions are disabled by default and are activated through -analytical-jac
   and -analytical-ic, respectively.

 *)

open Sundials
module Nvector = Nvector_serial

let (+) = (+.)
let (-) = (-.)
let ( * ) = ( *. )
let (/) = (/.)
let neg x = -. x

let dt = (1. / 60.) (* time step [s] *)
let t_end = 10.0    (* simulation end time [s] *)

(* Setup & auxiliary functions *)

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

(* Initial conditions.  Because p is the only algebraic variable, it can be
   computed from the rest of the initial values using calc_ic_ya_yd'; however,
   see the problem with higher-index formulation noted below, after the
   residual function.  *)
let theta0 = pi / 2.       (* angle with vertical axis *)
let (x0, y0) = (r * sin theta0, -. r * cos theta0)

(* There are 5 variables.  *)
let x, y, vx, vy, p = 0, 1, 2, 3, 4

(* There are 5 residues (i.e. 5 equations).  *)
let vx_x   = 0      (* Relationship between vx and x. *)
let vy_y   = 1      (* Relationship between vy and y. *)
let acc_x  = 2      (* Acceleration in the x direction. *)
let acc_y  = 3      (* Acceleration in the y direction. *)
let constr = 4      (* Constraint. *)
let neqs   = 5

let residual t vars vars' res =
  res.{vx_x}   <-  vars.{vx}  -  vars'.{x};
  res.{vy_y}   <-  vars.{vy}  -  vars'.{y};
  res.{acc_x}  <- vars'.{vx}  -  vars.{p} * vars.{x};
  res.{acc_y}  <- vars'.{vy}  -  vars.{p} * vars.{y}  +  g;
  res.{constr} <- vars.{x} * vars'.{vx} + vars.{y} * vars'.{vy}
                  + vars.{vx} * vars'.{x} + vars.{vy} * vars'.{y}

let jac Ida.({ jac_y = vars ; jac_y' = vars'; jac_coef = c }) out =
  let out = Matrix.Dense.unwrap out in
  out.{x,  vx_x}   <- -.c;
  out.{vx, vx_x}   <- 1.;
  out.{y,  vy_y}   <- -.c;
  out.{vy, vy_y}   <- 1.;
  out.{x,  acc_x}  <- -. vars.{p};
  out.{vx, acc_x}  <- c;
  out.{p,  acc_x}  <- -.vars.{x};
  out.{y,  acc_y}  <- -.vars.{p};
  out.{vy, acc_y}  <- c;
  out.{p,  acc_y}  <- -.vars.{y};
  out.{x,  constr} <- c * vars.{vx}  +  vars'.{vx};
  out.{y,  constr} <- c * vars.{vy}  +  vars'.{vy};
  out.{vx, constr} <- c * vars.{x}   +  vars'.{x};
  out.{vy, constr} <- c * vars.{y}   +  vars'.{y}

let roots t vars vars' r =
  let (wx,wy) = wall in
  r.{0} <- vars.{x} - vars.{y} * (wx / wy)

let d = Ida.VarId.differential
let a = Ida.VarId.algebraic

let main () =
  let vars  = RealArray.of_list [x0; y0; 0.; 0.; 0.] in
  let vars' = RealArray.make neqs 0. in
  let var_types = Nvector.wrap (RealArray.of_list [ d; d; d; d; a ]) in
  let nv_vars, nv_vars' = Nvector.(wrap vars, wrap vars') in

  let s = Ida.(init Dls.(solver ~jac (dense nv_vars (Matrix.dense 5)))
                    (SStolerances (1e-9, 1e-9))
                    residual ~roots:(1, roots) 0. nv_vars nv_vars')
  in
  Ida.set_id s var_types;
  Ida.set_suppress_alg s true;
  Ida.calc_ic_ya_yd' s ~y:nv_vars ~y':nv_vars' ~varid:var_types dt;

  let rec stepto tnext t =
    if t >= tnext then t else
    match Ida.solve_normal s tnext nv_vars nv_vars' with
    | (tret, Ida.RootsFound) ->
        vars.{vx} <- k * vars.{vx};
        vars.{vy} <- k * vars.{vy};
        Ida.reinit s tret nv_vars nv_vars';
        Ida.calc_ic_ya_yd' s ~y:nv_vars ~y':nv_vars' ~varid:var_types (t + dt);
        stepto tnext tret
    | (tret, _) -> tret
  in
  let rec showloop t = if t < t_end then begin
    Showpendulum.show (vars.{x}, vars.{y});
    showloop (stepto (t + dt) t)
  end in
  showloop 0.0

let _ =
  try
    Showpendulum.start (1.5 * r) dt false pivot wall;
    main ();
    Showpendulum.stop ()
  with Graphics.Graphic_failure _ -> ()

