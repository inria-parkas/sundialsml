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
 * In this example we consider the time-dependent damped harmonic oscillator
 *    q'(t) = p(t) exp(-F(t))
 *    p'(t) = -(F(t) * p + omega^2(t) * q)
 * With the initial conditions q(0) = 1, p(0) = 0.
 * The Hamiltonian for the system is
 *    H(p,q,t) = (p^2 * exp(-F(t)))/2 + (omega^2(t) * q^2 * exp(F(t)))/2
 * where omega(t) = cos(t/2), F(t) = 0.018*sin(t/pi).
 * We simulate the problem on t = [0, 30] using the symplectic methods in
 * SPRKStep.
 *
 * This is example 7.2 from:
 * Struckmeier, J., & Riedel, C. (2002). Canonical transformations and exact
 * invariants for time‐dependent Hamiltonian systems. Annalen der Physik, 11(1),
 * 15-38.
 *
 * The example has the following command line arguments:
 *   --order <int>               the order of the method to use (default 4)
 *   --dt <Real>                 the fixed-time step size to use (default 0.01)
 *   --nout <int>                the number of output times (default 100)
 *   --disable-tstop             turns off tstop mode
 *   --use-compensated-sums      turns on compensated summation in ARKODE where
 *                               applicable
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
let tf               = ref (10.0 *. pi)
let dt               = ref 1e-3

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
      ( "ark_damped_harmonic_symplectic: an ARKODE example demonstrating "
        ^ "the SPRKStep time-stepping module solving a time-dependent "
        ^ "damped harmonic oscillator\n")

let omega t = cos(t /. 2.0)

let f t = 0.018 *. sin(t /. pi)

let hamiltonian (yvec : Nvector_serial.t) t =
  let y = Nvector.unwrap yvec in
  let p = y.{0} in
  let q = y.{1} in
  (p *. p *. exp(-. f(t))) /. 2.0 +.
  (omega(t) *. omega(t) *. q *. q *. exp(f(t))) /. 2.0

let qdot t (y : RealArray.t) (ydot : RealArray.t) =
  let p = y.{0} in
  ydot.{1} <- p *. exp(-. f(t))

let pdot t (y : RealArray.t) (ydot : RealArray.t) =
  let p = y.{0} in
  let q = y.{1} in
  ydot.{0} <- -.(f(t) *. p +. omega(t) *. omega(t) *. q)

let main () =
  (* Parse the command line arguments *)
  parse_args ();

  (* Default integrator options *)
  let order            = !order in
  let use_compsums     = !use_compsums in
  let num_output_times = !num_output_times in

  (* Default problem parameters *)
  let tf    = !tf in
  let dt    = !dt in
  let dTout = (tf -. t0) /. float_of_int num_output_times in

  (* Create the SUNDIALS context object for this simulation *)
  let context = Context.make () in

  printf "\n   Begin time-dependent damped harmonic oscillator problem\n\n";

  (* Allocate our state vector *)
  let y = Nvector_serial.make ~context 2 0.0 in

  (* Fill the initial conditions *)
  let ydata = Nvector.unwrap y in
  ydata.{0} <- 0.0; (* \dot{q} = p *)
  ydata.{1} <- 1.0; (* \ddot{q} = \dot{p} *)

  (* Create SPRKStep integrator *)
  let arkode_mem =
    SPRKStep.init ~context ~step:dt ~order ~f1:qdot ~f2:pdot t0 y
  in
  SPRKStep.set_use_compensated_sums arkode_mem use_compsums;
  SPRKStep.set_max_num_steps arkode_mem (int_of_float (ceil(tf /. dt)) + 2);

  (* Print out starting Hamiltonian before integrating *)
  let tret = t0 in
  let tout = ref (t0 +. dTout) in
  (* Output current integration status *)
  printf "t = %.6f, q(t) = %.6f, H = %.6f\n" tret ydata.{1} (hamiltonian y tret);

  (* Do integration *)
  for _iout = 0 to num_output_times - 1 do
    if !use_tstop then SPRKStep.set_stop_time arkode_mem !tout;
    let tret, _ =SPRKStep.evolve_normal arkode_mem !tout y in

    (* Output current integration status *)
    printf "t = %.6f, q(t) = %.6f, H = %.6f\n" tret ydata.{1} (hamiltonian y tret);

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

