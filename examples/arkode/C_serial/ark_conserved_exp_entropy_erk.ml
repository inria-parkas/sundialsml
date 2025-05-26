(* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *---------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2025.
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example problem is adapted from:
 *
 * H. Ranocha, M. Sayyari, L. Dalcin, M. Parsani, and D.I. Ketcheson,
 * "Relaxation Runge-Kutta Methods: Fully-Discrete Explicit Entropy-Stable
 * Schemes for the Compressible Euler and Navier-Stokes Equations," SIAM Journal
 * on Scientific Computing, 42(2), 2020, https://doi.org/10.1137/19M1263480.
 * -----------------------------------------------------------------------------
 * This example evolves system
 *
 *   du/dt = -exp(v)
 *   dv/dt =  exp(u)
 *
 * for t in the interval [0, 5] with the initial condition
 *
 *   u(0) = 1.0
 *   v(0) = 0.5
 *
 * The system has the analytic solution
 *
 *   u = log(e + e^(3/2)) - log(b)
 *   v = log(a * e^(a * t)) - log(b)
 *
 * where log is the natural logarithm, a = sqrt(e) + e, and
 * b = sqrt(e) + e^(a * t).
 *
 * The conserved exponential entropy for the system is given by
 * ent(u,v) = exp(u) + exp(v) with the Jacobian
 * ent'(u,v) = [ de/du de/dv ]^T = [ exp(u) exp(v) ]^T.
 *
 * The problem is advanced in time with an explicit or implicit relaxed
 * Runge-Kutta method from ERKStep to ensure conservation of the entropy.
 * ---------------------------------------------------------------------------*)

open Sundials
module ERKStep = Arkode.ERKStep
module DM = Matrix.Dense

let printf = Printf.printf
let fprintf = Printf.fprintf

(* Value of the natural number e *)
let eval = 2.718281828459045235360287471352662497757247093699959574966

(* ----------------------- *
 * User-supplied functions *
 * ----------------------- *)

(* ODE RHS function f(t,y). *)
let f _t (ydata : RealArray.t) (fdata : RealArray.t) =
  fdata.{0} <- -. exp(ydata.{1});
  fdata.{1} <- exp(ydata.{0})

(* Entropy function e(y) *)
let ent (ydata : RealArray.t) =
  exp ydata.{0} +. exp ydata.{1}

(* Entropy function Jacobian Je(y) = de/dy *)
let jac_ent (ydata : RealArray.t) (jdata : RealArray.t) =
  jdata.{0} <- exp(ydata.{0});
  jdata.{1} <- exp(ydata.{1})

(* ----------------- *
 * Utility functions *
 * ----------------- *)

(* Analytic solution *)
let ans t (ydata : RealArray.t) =
  let a = sqrt eval +. eval in
  let b = sqrt eval +. exp(a *. t) in
  ydata.{0} <- log(eval +. exp(1.5)) -. log(b);
  ydata.{1} <- log(a *. exp(a *. t)) -. log(b)

(* ------------ *
 * Main Program *
 * ------------ *)

let main () =
  (* Initial and final times *)
  let t0 = 0.0 in
  let tf = 5.0 in

  (* Relative and absolute tolerances *)
  let reltol = 1.0e-6 in
  let abstol = 1.0e-10 in

  let argc = Array.length Sys.argv in

  (* -------------------- *
   * Output Problem Setup *
   * -------------------- *)

  (* enable relaxation *)
  let relax  = if argc > 1 then int_of_string Sys.argv.(1) <> 0 else true in
  (* adaptive stepping *)
  let fixed_h = if argc > 2 then float_of_string Sys.argv.(2) else 0.0 in

  printf "\nConserved Exponential Entropy problem:\n";
  printf "   method     = ERK\n";
  printf "   reltol     = %.1e\n" reltol;
  printf "   abstol     = %.1e\n" abstol;
  if fixed_h > 0.0 then printf "   fixed h    = %.1e\n" fixed_h;
  if relax then printf "   relaxation = ON\n"
  else printf "   relaxation = OFF\n";
  printf "\n";

  (* ------------ *
   * Setup ARKODE *
   * ------------ *)

  (* Create the SUNDIALS context object for this simulation *)
  let context = Context.make () in

  (* Create serial vector and set the initial condition values *)
  let y = Nvector_serial.make ~context 2 0.0 in
  let ydata = Nvector.unwrap y in
  ydata.{0} <- 1.0;
  ydata.{1} <- 0.5;

  let ytrue = Nvector.clone y in
  let ytdata = Nvector.unwrap ytrue in

  (* Initialize ERKStep *)
  let arkode_mem = ERKStep.(init ~context (SStolerances (reltol, abstol)) f t0 y) in

  (* Enable relaxation methods *)
  if relax then ERKStep.Relax.enable arkode_mem ent jac_ent;

  if fixed_h > 0.0 then ERKStep.set_fixed_step arkode_mem (Some fixed_h);

  (* Open output stream for results, output comment line *)
  let ufid = open_out "ark_conserved_exp_entropy_erk.txt" in
  fprintf ufid "# vars: t u v entropy u_err v_err entropy_error\n";

  (* --------------- *
   * Advance in Time *
   * --------------- *)

  (* Output the initial condition and entropy *)
  let ent0 = ent ydata in

  fprintf ufid
    "%23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n"
    t0 ydata.{0} ydata.{1} ent0 0.0 0.0 0.0;

  printf " step   t              u              v              e              delta e\n";
  printf " -------------------------------------------------------------------------------\n";
  printf "%5d %14.6e %14.6e %14.6e %14.6e %14.6e\n" 0 t0 ydata.{0} ydata.{1} ent0 0.0;

  let rec loop t =
    if t >= tf then ()
    else
      (* Evolve in time *)
      let t, _ = ERKStep.evolve_one_step arkode_mem tf y in

      (* Output solution and errors *)
      let entn = ent ydata in
      ans t ytdata;

      let ent_err = entn -. ent0 in
      let u_err = ydata.{0} -. ytdata.{0} in
      let v_err = ydata.{1} -. ytdata.{1} in

      (* Output to the screen periodically *)
      let nst = ERKStep.get_num_steps arkode_mem in

      if nst mod 40 = 0 then
        printf "%5d %14.6e %14.6e %14.6e %14.6e %14.6e\n"
               nst t ydata.{0} ydata.{1} entn ent_err;

      (* Write all steps to file *)
      fprintf ufid
              "%23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n"
              t ydata.{0} ydata.{1} entn u_err v_err ent_err;
      loop t
  in
  loop t0;
  printf " -------------------------------------------------------------------------------\n";

  (* ------------ *
   * Output Stats *
   * ------------ *)

  (* Get final statistics on how the solve progressed *)
  let nst = ERKStep.get_num_steps arkode_mem
  and nst_a = ERKStep.get_num_step_attempts arkode_mem
  and netf = ERKStep.get_num_err_test_fails arkode_mem
  and nfe = ERKStep.get_num_rhs_evals arkode_mem
  in
  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total number of error test failures = %d\n" netf;
  printf "   Total RHS evals = %d\n" nfe;

  if relax then begin
    let nre    = ERKStep.Relax.get_num_fn_evals arkode_mem
    and nrje   = ERKStep.Relax.get_num_jac_evals arkode_mem
    and nrf    = ERKStep.Relax.get_num_fails arkode_mem
    and nrbf   = ERKStep.Relax.get_num_bound_fails arkode_mem
    and nrnlsf = ERKStep.Relax.get_num_solve_fails arkode_mem
    and nrnlsi = ERKStep.Relax.get_num_solve_iters arkode_mem
    in
    printf "   Total Relaxation Fn evals    = %d\n" nre;
    printf "   Total Relaxation Jac evals   = %d\n" nrje;
    printf "   Total Relaxation fails       = %d\n" nrf;
    printf "   Total Relaxation bound fails = %d\n" nrbf;
    printf "   Total Relaxation NLS fails   = %d\n" nrnlsf;
    printf "   Total Relaxation NLS iters   = %d\n" nrnlsi;
  end;
  printf "\n"

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

