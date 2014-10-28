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

module RealArray = Sundials.RealArray
module Nvector = Nvector_serial

(* Commandline Options *)
type constraint_form =
| C_xx_yy_rr
| C_xx'_yy'
| C_vxx_yy'
| C_xx'_vyy
| C_vxx_vyy
| C_vx'x_vy'y_x'x'_y'y'
| C_vx'x_vy'y_x'x'_vyy'
| C_vx'x_vy'y_x'x'_vyvy
| C_vx'x_vy'y_vxx'_y'y'
| C_vx'x_vy'y_vxx'_vyy'
| C_vx'x_vy'y_vxx'_vyvy
| C_vx'x_vy'y_vxvx_y'y'
| C_vx'x_vy'y_vxvx_vyy'
| C_vx'x_vy'y_vxvx_vyvy
let constraint_form           = ref C_vx'x_vy'y_vxx'_vyy'
let use_analytical_correction = ref false
let check_consistency         = ref false
let use_analytical_jac        = ref false
let trace = ref false
let log = ref false
let show = ref true
let delay = ref true
let dt = ref (1./.60.)                  (* time step [s] *)
let t_end = ref 10.0                    (* simulation end time [s] *)

let args = [
    ("-t", Arg.Set_float t_end, "Duration in seconds.");
    ("-trace", Arg.Set trace, "Show a trace of mass positions.");
    ("-d", Arg.Set_float dt, "Set the default time step.");
    ("-noshow", Arg.Clear show, "Disable the graphical display.");
    ("-nodelay", Arg.Clear delay, "No delays between frames.");
    ("-constraint", Arg.String (fun s ->
      constraint_form :=
        match s with
        | "C_xx_yy_rr" -> C_xx_yy_rr
        | "C_xx'_yy'" -> C_xx'_yy'
        | "C_vxx_yy'" -> C_vxx_yy'
        | "C_xx'_vyy" -> C_xx'_vyy
        | "C_vxx_vyy" -> C_vxx_vyy
        | "C_vx'x_vy'y_x'x'_y'y'" -> C_vx'x_vy'y_x'x'_y'y'
        | "C_vx'x_vy'y_x'x'_vyy'" -> C_vx'x_vy'y_x'x'_vyy'
        | "C_vx'x_vy'y_x'x'_vyvy" -> C_vx'x_vy'y_x'x'_vyvy
        | "C_vx'x_vy'y_vxx'_y'y'" -> C_vx'x_vy'y_vxx'_y'y'
        | "C_vx'x_vy'y_vxx'_vyy'" -> C_vx'x_vy'y_vxx'_vyy'
        | "C_vx'x_vy'y_vxx'_vyvy" -> C_vx'x_vy'y_vxx'_vyvy
        | "C_vx'x_vy'y_vxvx_y'y'" -> C_vx'x_vy'y_vxvx_y'y'
        | "C_vx'x_vy'y_vxvx_vyy'" -> C_vx'x_vy'y_vxvx_vyy'
        | "C_vx'x_vy'y_vxvx_vyvy" -> C_vx'x_vy'y_vxvx_vyvy
        | _ -> raise (Invalid_argument
                        ("unrecognized constraint format: " ^ s))),
     "Formula for the constraint (see constraint_form type in the source)");
    ("-analytical-ic", Arg.Set use_analytical_correction,
     "Use analytically derived initial condition calculator");
    ("-check", Arg.Set check_consistency,
     "Check consistency upon every reset");
    ("-analytical-jac", Arg.Set use_analytical_jac,
     "Use analytically derived Jacobian instead of IDA's fallback");
    ("-log",
     Arg.Unit (fun () -> log := true; show := false; delay := false),
     "Log state variables to stdout (implies -noshow and -nodelay).");
]


(* Setup & auxiliary functions *)
module Matrix = Dls.DenseMatrix

let pi = 4. *. atan (1.)

let degree_to_radian x = x *. pi /. 180.

let show_nvector (a : RealArray.t) =
  let n = Bigarray.Array1.dim a in
  let str = ref "[" in
  for i = 0 to n-2 do
    str := !str ^ Printf.sprintf "%g," a.{i}
  done;
  if 0 <= n-1 then str := !str ^ Printf.sprintf "%g" a.{n-1};
  str := !str ^ "]";
  !str

let print_nvector a = Printf.printf "%s" (show_nvector a)

let print_with_time t v =
  Printf.printf "%.15e" t;
  Sundials.RealArray.iter (Printf.printf "\t% .15e") v;
  print_newline ()

(* Problem constants *)
let r = 1.0                            (* length of rod [m] *)
let g = 9.8                            (* gravitational acceleration [m/s^2] *)
let k = 0.5                            (* elasticity of collision with wall *)
(* direction of the wall relative to the pivot *)
let wall_angle = degree_to_radian (30.)
let wall = (-. sin wall_angle, -. cos wall_angle)

(* Initial conditions.  Because p is the only algebraic variable, it can be
   computed from the rest of the initial values using calc_ic_ya_yd'; however,
   see the problem with higher-index formulation noted below, after the
   residual function.  *)
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

let vars = RealArray.create neqs
let vars' = RealArray.create neqs
let var_types =
  let d = Ida.VarId.Differential and a = Ida.VarId.Algebraic in
  Nvector.wrap (RealArray.of_list
    (List.map Ida.VarId.to_float [d; d; d; d; a]))

(* The residual function F.  *)
let residual t vars vars' res =
  let x  = vars.{x_i}  and x'  = vars'.{x_i}
  and y  = vars.{y_i}  and y'  = vars'.{y_i}
  and vx = vars.{vx_i} and vx' = vars'.{vx_i}
  and vy = vars.{vy_i} and vy' = vars'.{vy_i}
  and p  = vars.{p_i} in
  res.{vx_x}   <- vx -. x';
  res.{vy_y}   <- vy -. y';
  res.{acc_x}  <- vx' -. p *. x;
  res.{acc_y}  <- vy' -. p *. y +. g;
  res.{constr} <-
    match !constraint_form with
    | C_xx_yy_rr -> x*.x +. y*.y -. r*.r
    | C_xx'_yy' -> x*.x' +. y*.y'
    | C_vxx_yy' -> vx*.x +. y*.y'
    | C_xx'_vyy -> vy*.y +. x*.x'
    | C_vxx_vyy -> vx*.x +. vy*.y
    | C_vx'x_vy'y_x'x'_y'y' -> vx'*.x +. vy'*.y +. x'*.x' +. y'*.y'
    | C_vx'x_vy'y_x'x'_vyy' -> vx'*.x +. vy*.y' +. vy'*.y +. x'*.x'
    | C_vx'x_vy'y_x'x'_vyvy -> vx'*.x +. vy*.vy +. vy'*.y +. x'*.x'
    | C_vx'x_vy'y_vxx'_y'y' -> vx*.x' +. vx'*.x +. vy'*.y +. y'*.y'
    | C_vx'x_vy'y_vxx'_vyy' -> vx*.x' +. vx'*.x +. vy*.y' +. vy'*.y
    | C_vx'x_vy'y_vxx'_vyvy -> vx*.x' +. vx'*.x +. vy*.vy +. vy'*.y
    | C_vx'x_vy'y_vxvx_y'y' -> vx*.vx +. vx'*.x +. vy'*.y +. y'*.y'
    | C_vx'x_vy'y_vxvx_vyy' -> vx*.vx +. vx'*.x +. vy*.y' +. vy'*.y
    | C_vx'x_vy'y_vxvx_vyvy -> vx*.vx +. vx'*.x +. vy*.vy +. vy'*.y

let jac params out =
  let vars  = params.Ida.jac_y
  and vars' = params.Ida.jac_y'
  and c = params.Ida.jac_coef in
  let x  = vars.{x_i}  and x'  = vars'.{x_i}
  and y  = vars.{y_i}  and y'  = vars'.{y_i}
  and vx = vars.{vx_i} and vx' = vars'.{vx_i}
  and vy = vars.{vy_i} and vy' = vars'.{vy_i}
  and p  = vars.{p_i} in
  Matrix.set out vx_x  x_i  (-.c);
  Matrix.set out vx_x  y_i  (0.);
  Matrix.set out vx_x  vx_i (1.);
  Matrix.set out vx_x  vy_i (0.);
  Matrix.set out vx_x  p_i  (0.);
  Matrix.set out vy_y  x_i  (0.);
  Matrix.set out vy_y  y_i  (-.c);
  Matrix.set out vy_y  vx_i (0.);
  Matrix.set out vy_y  vy_i (1.);
  Matrix.set out vy_y  p_i  (0.);
  Matrix.set out acc_x x_i  (-.p);
  Matrix.set out acc_x y_i  (0.);
  Matrix.set out acc_x vx_i (c);
  Matrix.set out acc_x vy_i (0.);
  Matrix.set out acc_x p_i  (-.x);
  Matrix.set out acc_y x_i  (0.);
  Matrix.set out acc_y y_i  (-.p);
  Matrix.set out acc_y vx_i (0.);
  Matrix.set out acc_y vy_i (c);
  Matrix.set out acc_y p_i  (-.y);

  match !constraint_form with
  | C_xx_yy_rr ->
    (Matrix.set out constr  x_i (2.*.x);
     Matrix.set out constr  y_i (2.*.y);
     Matrix.set out constr  vx_i (0.);
     Matrix.set out constr  vy_i (0.);
     Matrix.set out constr  p_i (0.))
  | C_xx'_yy' ->
    (Matrix.set out constr  x_i (c*.x +. x');
     Matrix.set out constr  y_i (c*.y +. y');
     Matrix.set out constr  vx_i (0.);
     Matrix.set out constr  vy_i (0.);
     Matrix.set out constr  p_i (0.))
  | C_vxx_yy' ->
    (Matrix.set out constr  x_i (vx);
     Matrix.set out constr  y_i (c*.y +. y');
     Matrix.set out constr  vx_i (x);
     Matrix.set out constr  vy_i (0.);
     Matrix.set out constr  p_i (0.))
  | C_xx'_vyy ->
    (Matrix.set out constr  x_i (c*.x +. x');
     Matrix.set out constr  y_i (vy);
     Matrix.set out constr  vx_i (0.);
     Matrix.set out constr  vy_i (y);
     Matrix.set out constr  p_i (0.))
  | C_vxx_vyy ->
    (Matrix.set out constr  x_i (vx);
     Matrix.set out constr  y_i (vy);
     Matrix.set out constr  vx_i (x);
     Matrix.set out constr  vy_i (y);
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_x'x'_y'y' ->
    (Matrix.set out constr  x_i (2.*.c*.x' +. vx');
     Matrix.set out constr  y_i (2.*.c*.y' +. vy');
     Matrix.set out constr  vx_i (c*.x);
     Matrix.set out constr  vy_i (c*.y);
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_x'x'_vyy' ->
    (Matrix.set out constr  x_i (2.*.c*.x' +. vx');
     Matrix.set out constr  y_i (c*.vy +. vy');
     Matrix.set out constr  vx_i (c*.x);
     Matrix.set out constr  vy_i (c*.y +. y');
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_x'x'_vyvy ->
    (Matrix.set out constr  x_i (2.*.c*.x' +. vx');
     Matrix.set out constr  y_i (vy');
     Matrix.set out constr  vx_i (c*.x);
     Matrix.set out constr  vy_i (c*.y +. 2.*.vy);
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_vxx'_y'y' ->
    (Matrix.set out constr  x_i (c*.vx +. vx');
     Matrix.set out constr  y_i (2.*.c*.y' +. vy');
     Matrix.set out constr  vx_i (c*.x +. x');
     Matrix.set out constr  vy_i (c*.y);
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_vxx'_vyy' ->
    (Matrix.set out constr  x_i (c*.vx +. vx');
     Matrix.set out constr  y_i (c*.vy +. vy');
     Matrix.set out constr  vx_i (c*.x +. x');
     Matrix.set out constr  vy_i (c*.y +. y');
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_vxx'_vyvy ->
    (Matrix.set out constr  x_i (c*.vx +. vx');
     Matrix.set out constr  y_i (vy');
     Matrix.set out constr  vx_i (c*.x +. x');
     Matrix.set out constr  vy_i (c*.y +. 2.*.vy);
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_vxvx_y'y' ->
    (Matrix.set out constr  x_i (vx');
     Matrix.set out constr  y_i (2.*.c*.y' +. vy');
     Matrix.set out constr  vx_i (c*.x +. 2.*.vx);
     Matrix.set out constr  vy_i (c*.y);
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_vxvx_vyy' ->
    (Matrix.set out constr  x_i (vx');
     Matrix.set out constr  y_i (c*.vy +. vy');
     Matrix.set out constr  vx_i (c*.x +. 2.*.vx);
     Matrix.set out constr  vy_i (c*.y +. y');
     Matrix.set out constr  p_i (0.))
  | C_vx'x_vy'y_vxvx_vyvy ->
    (Matrix.set out constr  x_i (vx');
     Matrix.set out constr  y_i (vy');
     Matrix.set out constr  vx_i (c*.x +. 2.*.vx);
     Matrix.set out constr  vy_i (c*.y +. 2.*.vy);
     Matrix.set out constr  p_i (0.))

(* The root function -- computes distance from the wall along x axis.  The wall
   is a line that connects pivot = (0,0) and wall = (wx,wy) where wall is a
   constant defined above.  *)
let distance_from_wall vars =
  let x       = vars.{x_i}
  and y       = vars.{y_i}
  and (wx,wy) = wall in
  x -. y *. (wx /. wy)
let roots t vars vars' r =
  r.{0} <- distance_from_wall vars

(* A hand-written function to compute initial conditions from given values for
   x,y,vx,vy.  The vars vector's x,y,vx,vy fields must be initialized with
   those values, and this function computes the rest.  *)
let init_from_xy_vxvy vars vars' =
  let x  = vars.{x_i}  and y  = vars.{y_i}
  and vx = vars.{vx_i} and vy = vars.{vy_i} in
  if !check_consistency && abs_float (x*.vx +. y*.vy) > 1e-6 then
    (Printf.printf "vars  = %s\nvars' = %s\n"
       (show_nvector vars) (show_nvector vars');
     raise (Invalid_argument ("(vx,vy) dot (x,y) is too large: "
                              ^ string_of_float (x*.vx +. y*.vy))));
  if !check_consistency && abs_float (x*.x +. y*.y -. r*.r) > 1e-6 then
    raise (Invalid_argument "x^2 + y^2 - r^2 is too large");
  let p  = (g*.y -. vx*.vx -. vy*.vy) /. (r*.r) in
  vars.{p_i} <- p;
  vars'.{x_i} <- vx;
  vars'.{y_i} <- vy;
  vars'.{vx_i} <- p*.x;
  vars'.{vy_i} <- p*.y -. g;
  let vx' = vars'.{vx_i} and vy' = vars'.{vy_i} in
  vars'.{p_i} <- (g*.vy -. 2.*.vx'*.vx -. 2.*.vy'*.vy) /. (r*.r)

(* Rendering parameters  *)

(* Position of pivot on the screen, scaled so that the screen is 1.0x1.0.  The
   solver will take the pivot as the origin. *)
let pivot = (0.5, 0.8)

let check_satisfaction =
  let res = RealArray.create neqs in
  fun t vars vars' ->
    residual t vars vars' res;
    if !check_consistency then
      begin
        let norm = Nvector_serial.DataOps.n_vmaxnorm res in
        if norm > 1e-5 then
          raise (Failure 
                   (Printf.sprintf "initial residue too large: t = %g, ||res||=%g\nvars  = %s\nvars' = %s\nres   = %s\n"
                      t norm (show_nvector vars)
                      (show_nvector vars') (show_nvector res)))
      end

let print_with_time t v =
  Printf.printf "%e" t;
  Sundials.RealArray.iter (Printf.printf "\t% e") v;
  print_newline ()

let main () =
  Arg.parse args (fun _ -> ()) "pendulum: simulate a pendulum hitting against an oblique wall";
  if !log then print_with_time 0.0 vars;
  let t_delay = if !delay then !dt else 0. in
  if !show then Showpendulum.start (1.5 *. r) t_delay !trace pivot wall;
  let frames = int_of_float (!t_end /. !dt) in

  vars.{x_i} <- x0;
  vars.{y_i} <- y0;
  vars.{vx_i} <- x'0;
  vars.{vy_i} <- y'0;

  let nv_vars  = Nvector.wrap vars in
  let nv_vars' = Nvector.wrap vars' in

  if !use_analytical_correction then init_from_xy_vxvy vars vars';

  let solver =
    if !use_analytical_jac then Ida.Dls.dense ~jac:jac ()
    else Ida.Dls.dense ()
  in
  let ida = Ida.init solver (Ida.SStolerances (1e-9, 1e-9)) residual
                     ~roots:(1, roots) 0. nv_vars nv_vars'
  in
  Ida.set_all_root_directions ida Sundials.RootDirs.Decreasing;
  if !use_analytical_jac then Ida.Dls.set_dense_jac_fn ida jac;

  Ida.set_id ida var_types;
  Ida.set_suppress_alg ida true;
  if not !use_analytical_correction
  then Ida.calc_ic_ya_yd' ~y:nv_vars ~y':nv_vars' ida ~varid:var_types !dt;

  check_satisfaction 0. vars vars';

  let t = ref 0.
  and tnext = ref !dt
  in
  for i = 1 to frames do
    tnext := !t +. !dt;
    while !t < !tnext do
      let (tret, flag) = Ida.solve_normal ida !tnext nv_vars nv_vars' in
      t := tret;
      if !show then Showpendulum.show (vars.{x_i}, vars.{y_i});
      if !log then print_with_time tret vars;

      if flag = Ida.RootsFound then
        (Printf.printf "Bang!  Hit against the wall.\n";
         Printf.printf "vars  = %s\nvars' = %s\n"
           (show_nvector vars) (show_nvector vars');
         vars.{vx_i} <- -. k *. vars.{vx_i};
         vars.{vy_i} <- -. k *. vars.{vy_i};

         if !use_analytical_correction
         then (init_from_xy_vxvy vars vars';
               Ida.reinit ida !t nv_vars nv_vars')
         else (Ida.reinit ida !t nv_vars nv_vars';
               Ida.calc_ic_ya_yd' ~y:nv_vars ~y':nv_vars' ida
               ~varid:var_types (!t +. !dt));

         check_satisfaction !t vars vars';
        )
    done
  done;
  Printf.printf "error test failures = %d\n"
    (Ida.get_num_err_test_fails ida);
  Showpendulum.stop ()

let _ =
  try main ()
  with
    (* This failure happens when you close the window before the simulation
       ends.  *)
    Graphics.Graphic_failure _ -> ()
