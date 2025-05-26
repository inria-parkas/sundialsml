(* clang-format off *)
(* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *-----------------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2025.
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * In this example we consider the simple harmonic oscillator
 *    x''(t) + omega^2*x(t) = 0.
 * We rewrite the second order ODE as the first order ODE model
 *    x'(t) = v(t)
 *    v'(t) = -omega^2*x(t).
 * With the initial conditions x(0) = x0 and v(0) = v0,
 * the analytical solution is
 *    x(t) = A*cos(t*omega + phi),
 *    v(t) = -A*omega*sin(t*omega + phi)
 * where A = sqrt(x0^2 + v0^2/omega) and tan(phi) = v0/(omega*x0).
 * The total energy (potential + kinetic) in this system is
 *    E = (v^2 + omega^2*x^2) / 2
 * E is conserved and is the system Hamiltonian.
 * We simulate the problem on t = [0, 2pi] using the symplectic methods
 * in SPRKStep. Symplectic methods will approximately conserve E.
 *
 * The example has the following command line arguments:
 *   --order <int>               the order of the method to use (default 4)
 *   --dt <Real>                 the fixed-time step size to use (default 0.01)
 *   --nout <int>                the number of output times (default 100)
 *   --use-compensated-sums      turns on compensated summation in ARKODE where
 *                               applicable
 *   --disable-tstop             turns off tstop mode
 * --------------------------------------------------------------------------*)

open Sundials
module SPRKStep = Arkode.SPRKStep

let printf = Printf.printf
let fprintf = Printf.fprintf
let eprintf = Printf.eprintf

let pi = 3.14159265358979323846264338327950

let order            = ref 4
let num_output_times = ref 8
let use_compsums     = ref false
let use_tstop        = ref true
let dt               = ref 1e-3
let tf               = ref (2.0 *. pi)

let t0 = 0.0

let parse_args () =
  Arg.parse [
      ("--order", Arg.Set_int order,
       "the order of the method to use (default 4)");
      ("--dt", Arg.Set_float dt,
       "the fixed-time step size to use (default 0.01)\n");
      ("--nout", Arg.Set_int num_output_times,
       "the number of output times (default 100)");
      ("--use-compensated-sums", Arg.Set use_compsums,
       "turns on compensated summation in ARKODE where applicable");
      ("--disable-tstop", Arg.Clear use_tstop,
       "turns off tstop mode");
    ] (fun _ -> raise Arg.(Bad "unrecognized option"))
      ( "ark_harmonic_symplectic: an ARKODE example demonstrating "
        ^ "the SPRKStep time-stepping module solving a simple harmonic "
        ^ "oscillator\n")

type user_data = {
    a     : float;
    phi   : float;
    omega : float;
  }

let solution { a; phi; omega } t _y (sol : RealArray.t) =
  (* compute solution *)
  sol.{0} <- a *. cos(omega *. t +. phi);
  sol.{1} <- -. a *. omega *. sin(omega *. t +. phi)

let energy { omega; _ } (y : RealArray.t) _dt =
  let x = y.{0} in
  let v = y.{1} in
  let omega2 = omega *. omega in
  (v *. v +. omega2 *. x *. x) /. 2.0

let xdot _ud _t (y : RealArray.t) (ydot : RealArray.t) =
  let v = y.{1} in
  ydot.{0} <- v

let vdot { omega; _ } _t (y : RealArray.t) (ydot : RealArray.t) =
  let x = y.{0} in
  let omega2 = omega *. omega in
  ydot.{1} <- -. omega2 *. x

let main () =
  (* Parse the command line arguments *)
  parse_args ();

  let a = 10.0 in
  let phi = 0.0 in
  let omega = 1.0 in

  (* Default integrator options and problem parameters *)
  let order            = !order in
  let use_compsums     = !use_compsums in
  let num_output_times = !num_output_times in
  let tf    = !tf in
  let dt    = !dt in
  let dTout = (tf -. t0) /. float_of_int num_output_times in

  (* Create the SUNDIALS context object for this simulation *)
  let context = Context.make () in

  printf "\n   Begin simple harmonic oscillator problem\n\n";
  let user_data = { a; phi; omega } in

  (* Allocate our state vector *)
  let y = Nvector_serial.make ~context 2 0.0 in
  let sol = Nvector.clone y in
  let soldata = Nvector.unwrap sol in

  (* Fill the initial conditions (x0 then v0) *)
  let ydata = Nvector.unwrap y in
  ydata.{0} <- a *. cos(phi);
  ydata.{1} <- -. a *. omega *. sin(phi);

  (* Create SPRKStep integrator *)
  let arkode_mem =
    SPRKStep.init ~context ~step:dt ~order
                  ~f1:(xdot user_data) ~f2:(vdot user_data) t0 y
  in
  SPRKStep.set_use_compensated_sums arkode_mem use_compsums;
  SPRKStep.set_max_num_steps arkode_mem (int_of_float (ceil(tf /. dt)) + 2);

  (* Print out starting energy, momentum before integrating *)
  let tret = t0 in
  let tout = ref (t0 +. dTout) in
  (* Output current integration status *)
  printf "t = %.6f, x(t) = %.6f, E = %.6f, sol. err = %.6f\n"
         tret ydata.{0} (energy user_data ydata dt) 0.0;

  (* Do integration *)
  for _iout = 0 to num_output_times - 1 do
    if !use_tstop then SPRKStep.set_stop_time arkode_mem !tout;
    let tret, _ =SPRKStep.evolve_normal arkode_mem !tout y in

    (* Compute the anaytical solution *)
    solution user_data tret y soldata;

    (* Compute L2 error *)
    Nvector_serial.Ops.linearsum 1.0 y (-1.0) sol sol;
    let err = sqrt(Nvector_serial.Ops.dotprod sol sol) in

    (* Output current integration status *)
    printf "t = %.6f, x(t) = %.6f, E = %.6f, sol. err = %.16e\n"
           tret ydata.{0} (energy user_data ydata dt) err;

    (* Check that solution error is within tolerance *)
    if err > max (dt /. Float.pow 10.0 (float_of_int (order - 2)))
                 (1000.0 *. Config.unit_roundoff)
    then (eprintf "FAILURE: solution error is too high\n"; exit 1);

    (* Check if the solve was successful, if so, update the time and continue *)
    tout := min tf (!tout +. dTout)
  done;

  printf "\n%!";
  SPRKStep.print_all_stats arkode_mem OutputTable

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
  for _ = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()

