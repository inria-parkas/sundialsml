(* ----------------------------------------------------------------------------- {{{
 * Programmer(s): Radu Serban and David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Oct 2021.
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves a simple pendulum equation in Cartesian coordinates where
 * the pendulum bob has mass 1 and is suspended from the origin with a rod of
 * length 1. The governing equations are
 *
 * x'  = vx
 * y'  = vy
 * vx' = -x * T
 * vy' = -y * T - g
 *
 * with the constraints
 *
 * x^2 + y^2 - 1 = 0
 * x * vx + y * vy = 0
 *
 * where x and y are the pendulum bob position, vx and vy are the bob velocity
 * in the x and y directions respectively, T is the tension in the rod, and
 * g is acceleration due to gravity chosen such that the pendulum has period 2.
 * The initial condition at t = 0 is x = 1, y = 0, vx = 0, and vy = 0.
 *
 * A reference solution is computed using the pendulum equation in terms of the
 * angle between the x-axis and the pendulum rod i.e., theta in [0, -pi]. The
 * governing equations are
 *
 * theta'  = vtheta
 * vtheta' = -g * cos(theta)
 *
 * where theta is the angle from the x-axis, vtheta is the angular velocity, and
 * g the same acceleration due to gravity from above. The initial condition at
 * t = 0 is theta = 0 and vtheta = 0.
 *
 * The Cartesian formulation is run to a final time tf (default 30) with and
 * without projection for various integration tolerances. The error in the
 * position and velocity at tf compared to the reference solution, the error in
 * the position constraint equation, and various integrator statistics are
 * printed to the screen for each run.
 *
 * When projection is enabled a user-supplied function is used to project the
 * position, velocity, and error to the constraint manifold.
 *
 * Optional command line inputs may be used to change the final simulation time
 * (default 30), the initial tolerance (default 1e-5), the number of outputs
 * (default 1), or disable error projection. Use the option --help for a list
 * of the command line flags.
 * --------------------------------------------------------------------------- }}}*)

open Sundials
let printf = Printf.printf
let fprintf = Printf.fprintf

(* Problem Constants *)
let grav = 13.750371636040745654980191559621114395801712

(* -----------------------------------------------------------------------------
 * Functions provided to CVODE
 * ---------------------------------------------------------------------------*)

(* ODE RHS function for the reference system *)
let fref t yydata fydata =
  fydata.{0} <- yydata.{1};                (* theta'          *)
  fydata.{1} <- -. grav *. cos yydata.{0}  (* -g * cos(theta) *)

(* ODE RHS function for the Cartesian system *)
let f t (yydata : RealArray.t) (fydata : RealArray.t) =
  (* Get vector components *)
  let x  = yydata.{0} in
  let y  = yydata.{1} in
  let xd = yydata.{2} in
  let yd = yydata.{3} in

  (* Compute tension *)
  let tmp = xd *. xd +. yd *. yd -. grav *. y in

  (* Compute RHS *)
  fydata.{0} <- xd;
  fydata.{1} <- yd;
  fydata.{2} <- -. x *. tmp;
  fydata.{3} <- -. y *. tmp -. grav

(* Projection function *)
let proj t (yydata : RealArray.t) (cdata : RealArray.t)
           epsProj (edata : RealArray.t) =
  (* Extract current solution *)
  let x  = yydata.{0} in
  let y  = yydata.{1} in
  let xd = yydata.{2} in
  let yd = yydata.{3} in

  (* Project positions *)

  let r = sqrt (x *. x +. y *. y) in

  let x_new = x /. r in
  let y_new = y /. r in

  (* Project velocities
   *
   *        +-            -+  +-    -+
   *        |  y*y    -x*y |  |  xd  |
   *  P v = |              |  |      |
   *        | -x*y     x*x |  |  yd  |
   *        +-            -+  +-    -+
   *)

  let xd_new =    xd *. y_new *. y_new -. yd *. x_new *. y_new in
  let yd_new = -. xd *. x_new *. y_new +. yd *. x_new *. x_new in

  (* Return position and velocity corrections *)

  cdata.{0} <- x_new  -. x;
  cdata.{1} <- y_new  -. y;
  cdata.{2} <- xd_new -. xd;
  cdata.{3} <- yd_new -. yd;

  (* Project error P * err *)
  let e1 = edata.{0} in
  let e2 = edata.{1} in
  let e3 = edata.{2} in
  let e4 = edata.{3} in

  let e1_new =    y_new *. y_new *. e1 -. x_new *. y_new *. e2 in
  let e2_new = -. x_new *. y_new *. e1 +. x_new *. x_new *. e2 in

  let e3_new =    y_new *. y_new *. e3 -. x_new *. y_new *. e4 in
  let e4_new = -. x_new *. y_new *. e3 +. x_new *. x_new *. e4 in

  edata.{0} <- e1_new;
  edata.{1} <- e2_new;
  edata.{2} <- e3_new;
  edata.{3} <- e4_new

(* -----------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*)

(* Print command line options *)
let input_help =
  "\nCommand line options:\n\
     --tol <tol>      : relative and absolute tolerance\n\
     --tf <time>      : final simulation time\n\
     --nout <outputs> : number of outputs\n\
     --noerrproj      : disable error projection\n"

type inputs = {
  tol : float;
  tf  : float;
  nout : int;
  projerr : bool;
}

(* Read command line unputs *)
let read_inputs () =
  let inputs = ref {
    tol     = 1.0e-5; (* integration tolerance   *)
    tf      = 30.0;   (* final integration time  *)
    nout    = 1;      (* number of outputs       *)
    projerr = true;   (* enable error projection *)
  } in
  let args = Arg.[
    "--tol",  Float (fun tol  -> inputs := { !inputs with tol })  , "";
    "--tf",   Float (fun tf   -> inputs := { !inputs with tf })   , "";
    "--nout", Int   (fun nout -> inputs := { !inputs with nout }) , "";
    "--noerrproj", Unit (fun () -> inputs := { !inputs with projerr = false }) , "";
  ] in
  let anon_fn s = fprintf stderr "ERROR: Invalid input %s\n%s" s input_help in
  Arg.parse args anon_fn input_help;
  !inputs

(* -----------------------------------------------------------------------------
 * Functions to integrate the Cartesian and reference systems
 * ---------------------------------------------------------------------------*)

(* Compute the Cartesian system solution *)
let get_sol cvode_mem yy0 tol tf nout proj projerr yref =
  (* Enable or disable projection *)
  if proj then begin
    printf("  YES   ");
    Cvode.set_proj_frequency cvode_mem 1;
    (* Enable or disable error projection *)
      Cvode.set_proj_err_est cvode_mem projerr
  end else begin
    Cvode.set_proj_frequency cvode_mem 0;
    printf "  NO    ";
  end;

  (* Create vector to store the solution *)
  let yy = Nvector_serial.make 4 0.0 in

  (* Copy initial condition into solution vector *)
  Nvector_serial.Ops.scale 1.0 yy0 yy;

  (* Get pointer to vector data *)
  let yydata = Nvector.unwrap yy in

  (* Reinitialize CVODE for this run *)
  Cvode.reinit cvode_mem 0.0 yy0;

  (* Set integration tolerances for this run *)
  Cvode.(set_tolerances cvode_mem (SStolerances (tol, tol)));

  (* Open output file *)
  let outname =
    if proj then Printf.sprintf "cvPendulum_dns_tol_%03.2e_proj.txt" tol
    else Printf.sprintf "cvPendulum_dns_tol_%03.2e.txt" tol
  in
  let fid = open_out outname in

  (* Output initial condition *)
  fprintf fid "%0.4e %14.6e %14.6e %14.6e %14.6e\n"
          0.0 yydata.{0} yydata.{1} yydata.{2} yydata.{3};

  (* Integrate to tf and peridoically output the solution *)
  let dtout = tf /. float nout in
  let tout = ref dtout in

  for out = 0 to nout - 1 do
    (* Set stop time (do not interpolate output) *)
    Cvode.set_stop_time cvode_mem !tout;

    (* Integrate to tout *)
    let t, _ = Cvode.solve_normal cvode_mem !tout yy in

    (* Write output *)
    fprintf fid "%0.4e %14.6e %14.6e %14.6e %14.6e\n"
            t yydata.{0} yydata.{1} yydata.{2} yydata.{3};

    (* Update output time *)
    tout := if out < nout - 1 then !tout +. dtout else tf
  done;

  (* Close output file *)
  close_out fid;

  (* Compute the constraint violation *)
  let x = yydata.{0} in
  let y = yydata.{1} in
  let g = abs_float (x *. x +. y *. y -. 1.0) in

  (* Compute the absolute error compared to the reference solution *)
  Nvector_serial.Ops.linearsum 1.0 yy (-1.0) yref yy;
  Nvector_serial.Ops.abs yy yy;

  let x  = yydata.{0} in
  let y  = yydata.{1} in
  let xd = yydata.{2} in
  let yd = yydata.{3} in

  (* Output errors *)
  printf "%8.2e  %8.2e  %8.2e  %8.2e  |  %8.2e  |" x y xd yd g;

  (* Get integrator stats *)
  let open Cvode in
  let nst     = get_num_steps cvode_mem
  and nfe     = get_num_rhs_evals cvode_mem
  and nsetups = get_num_lin_solv_setups cvode_mem
  and netf    = get_num_err_test_fails cvode_mem
  and ncfn    = get_num_nonlin_solv_conv_fails cvode_mem
  and nje     = Dls.get_num_jac_evals cvode_mem
  and nfeLS   = Dls.get_num_lin_rhs_evals cvode_mem
  in
  (* Output stats *)
  printf " %6d   %6d+%-4d     %4d (%3d)     |  %3d  %3d\n"
    nst nfe nfeLS nsetups nje ncfn netf

(* Compute the reference system solution *)
let ref_sol tf yref nout =
  let tol = 1.0e-14 in
  (* Set the initial condition *)
  let yydata = RealArray.of_array [|
                  0.0; (* theta *)
                  0.0  (* theta' *)
  |] in
  (* Create the solution vector *)
  let yy = Nvector_serial.wrap yydata in

  (* Create CVODE memory *)
  (* Initialize CVODE *)
  (* Set integration tolerances *)
  (* Create dense SUNMatrix for use in linear solves *)
  (* Create dense SUNLinearSolver object *)
  (* Attach the matrix and linear solver to CVODE *)
  let m = Matrix.dense 2 in
  let cvode_mem = Cvode.(init BDF
                           ~lsolver:Dls.(solver (dense yy m))
                           (SStolerances (tol, tol))
                           fref 0.0 yy)
  in
  (* Set CVODE optional inputs *)
  Cvode.set_max_num_steps cvode_mem 100000;
  Cvode.set_stop_time cvode_mem tf;

  (* Open output file *)
  let fid = open_out "cvPendulum_dns_ref.txt" in

  (* Output initial condition *)
  let th  = yydata.{0} in
  let thd = yydata.{1} in
  fprintf fid "%0.4e %14.6e %14.6e %14.6e %14.6e\n"
          0.0 (cos th) (sin th) (-. thd *. sin th) (thd *. cos th);

  (* Integrate to tf and periodically output the solution *)
  let dtout = tf /. float nout in
  let tout = ref dtout in

  for out = 0 to nout - 1 do
    (* Set stop time (do not interpolate output) *)
    Cvode.set_stop_time cvode_mem !tout;

    (* Integrate to tout *)
    let t, _ = Cvode.solve_normal cvode_mem tf yy in

    (* Write output *)
    let th  = yydata.{0} in
    let thd = yydata.{1} in
    fprintf fid "%0.4e %14.6e %14.6e %14.6e %14.6e\n"
            t (cos th) (sin th) (-.thd *. sin th) (thd *. cos th);

    (* Update output time *)
    tout := if out < nout - 1 then !tout +. dtout else tf
  done;

  (* Close output file *)
  close_out fid;

  (* Get solution components *)
  let th  = yydata.{0} in
  let thd = yydata.{1} in

  (* Convert to Cartesian reference solution *)
  let yydata = Nvector.unwrap yref in
  yydata.{0} <- cos th;
  yydata.{1} <- sin th;
  yydata.{2} <- -. thd *. sin th;
  yydata.{3} <-    thd *. cos th

(* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*)

let main () =
  (* Read command line inputs *)
  let { tol; tf; nout; projerr } as inputs = read_inputs () in

  (* Compute reference solution *)
  let yref = Nvector_serial.make 4 0.0 in

  ref_sol tf yref nout;

  (* Create serial vector to store the initial condition *)
  (* Set the initial condition values *)
  let yy0 = Nvector_serial.wrap
      (RealArray.of_array [| 1.0; (* x *)
                             0.0; (* y *)
                             0.0; (* xd *)
                             0.0; (* yd *) |])
  in

  (* Create CVODE memory *)
  (* Initialize CVODE *)
  (* Set integration tolerances *)
  (* Create dense SUNMatrix for use in linear solves *)
  (* Create dense SUNLinearSolver object *)
  (* Attach the matrix and linear solver to CVODE *)
  (* Set a user-supplied projection function *)
  let a = Matrix.dense 4 in
  let cvode_mem =
    Cvode.(init BDF ~lsolver:Dls.(solver (dense yy0 a))
                    (SStolerances (tol, tol))
                    f ~projfn:proj 0.0 yy0)
  in
  (* Set maximum number of steps between outputs *)
  Cvode.set_max_num_steps cvode_mem 50000;

  (* Compute the solution with various tolerances *)
  let tol = ref tol in
  for i = 0 to 4 do

    (* Output tolerance and output header for this run *)
    printf "\n\nTol = %8.2e\n" !tol;
    printf "Project    x         y";
    printf "         x'        y'     |     g      |    ";
    printf "nst     rhs eval    setups (J eval)  |   cf   ef\n";

    (* Compute solution with projection *)
    get_sol cvode_mem yy0 !tol tf nout true projerr yref;

    (* Compute solution without projection *)
    get_sol cvode_mem yy0 !tol tf nout false false yref;

    (* Reduce tolerance for next run *)
    tol := !tol /. 10.0
  done

(* Check environment variables for extra arguments.  *)
let reps =
  try int_of_string (Unix.getenv "NUM_REPS")
  with Not_found | Failure _ -> 1
let gc_at_end =
  try int_of_string (Unix.getenv "GC_AT_END") <> 0
  with Not_found | Failure _ -> false
let gc_each_rep =
  try int_of_string (Unix.getenv "GC_EACH_REP") <> 0
  with Not_found | Failure _ -> false

(* Entry point *)
let _ =
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()

