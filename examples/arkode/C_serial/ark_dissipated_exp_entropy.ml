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
 * This example evolves the equation du/dt = -exp(u) for t in the interval
 * [0, 5] with the initial condition u(0) = 0.5. The equation has the analytic
 * solution u(t) = -log(e^{-0.5} + t) and the dissipated exponential entropy is
 * given by ent(u) = exp(u) with Jacobian ent'(u) = de/du = exp(u).
 *
 * The problem is advanced in time with an explicit or implicit relaxed
 * Runge-Kutta method to ensure dissipation of the entropy.
 * ---------------------------------------------------------------------------*)

open Sundials
module ARKStep = Arkode.ARKStep
module DM = Matrix.Dense

let ge670 = match Sundials.Config.version with
            | 6, m, _, _ -> m >= 7
            | m, _, _, _ -> m > 6

let printf = Printf.printf
let fprintf = Printf.fprintf

(* ----------------------- *
 * User-supplied functions *
 * ----------------------- *)

(* ODE RHS function f(t,y). *)
let f _t (ydata : RealArray.t) (fdata : RealArray.t) =
  fdata.{0} <- -. exp(ydata.{0})

(* ODE RHS Jacobian function J(t,y) = df/dy. *)
let jac ARKStep.{ jac_y = (ydata : RealArray.t); _ } jdata =
  DM.set jdata 0 0 (-. exp(ydata.{0}))

(* Entropy function e(y) *)
let ent (ydata : RealArray.t) = exp ydata.{0}

(* Entropy function Jacobian Je(y) = de/dy *)
let jac_ent (ydata : RealArray.t) (jdata : RealArray.t) =
  jdata.{0} <- exp(ydata.{0})

(* ----------------- *
 * Utility functions *
 * ----------------- *)

(* Analytic solution *)
let ans t (ydata : RealArray.t) =
  ydata.{0} <- log(exp(-0.5) +. t)

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
  (* implicit          *)
  let implicit = if argc > 2 then int_of_string Sys.argv.(2) <> 0 else true in
  (* adaptive stepping *)
  let fixed_h = if argc > 3 then float_of_string Sys.argv.(3) else 0.0 in

  printf "\nDissipated Exponential Entropy problem:\n";
  if implicit then printf "   method     = DIRK\n"
  else printf "   method     = ERK\n";
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

  (* Initialize ARKStep *)
  let problem =
    if implicit then
      let a = Matrix.dense ~context 2 in
      let lsolver = ARKStep.Dls.solver ~jac
                      (LinearSolver.Direct.dense ~context y a)
      in
      ARKStep.implicit ~lsolver f
    else
      ARKStep.explicit f
  in
  let arkode_mem = ARKStep.(init ~context problem
                                 (SStolerances (reltol, abstol)) t0 y)
  in
  (* Select a Butcher table with non-negative b values *)
  if implicit then begin
    ARKStep.(set_table_name arkode_mem
                            ~itable:"ARKODE_ARK2_DIRK_3_1_2"
                            ~etable:"ARKODE_ERK_NONE"
                            ());
    if ge670 then ARKStep.set_nonlin_conv_coef arkode_mem 0.01
  end;

  (* Enable relaxation methods *)
  if relax then ARKStep.Relax.enable arkode_mem ent jac_ent;

  if fixed_h > 0.0 then ARKStep.set_fixed_step arkode_mem (Some fixed_h);

  (* Open output stream for results, output comment line *)
  let ufid = open_out "ark_dissipated_exp_entropy.txt" in
  fprintf ufid "# vars: t u entropy u_err delta_entropy\n";

  (* --------------- *
   * Advance in Time *
   * --------------- *)

  (* Output the initial condition and entropy *)
  let ent0 = ent ydata in

  fprintf ufid
    "%23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n"
    t0 ydata.{0} ydata.{1} ent0 0.0 0.0 0.0;

  fprintf ufid "%23.16e %23.16e %23.16e %23.16e %23.16e\n" t0 ydata.{0} ent0 0.0 0.0;

  printf " step   t              u              e              u_err          delta e\n";
  printf " -------------------------------------------------------------------------------\n";
  printf "%5d %14.6e %14.6e %14.6e %14.6e %14.6e\n" 0 t0 ydata.{0} ent0 0.0 0.0;

  let rec loop t =
    if t >= tf then ()
    else
      (* Evolve in time *)
      let t, _ = ARKStep.evolve_one_step arkode_mem tf y in

      (* Output solution and errors *)
      let entn = ent ydata in
      ans t ytdata;

      let delta_ent = entn -. ent0 in
      let u_err = ydata.{0} -. ytdata.{0} in

      (* Output to the screen periodically *)
      let nst = ARKStep.get_num_steps arkode_mem in

      if nst mod 40 = 0 then
        printf "%5d %14.6e %14.6e %14.6e %14.6e %14.6e\n"
               nst t ydata.{0} entn u_err delta_ent;

      (* Write all steps to file *)
      fprintf ufid "%23.16e %23.16e %23.16e %23.16e %23.16e\n"
                   t ydata.{0} entn u_err delta_ent;
      loop t
  in
  loop t0;
  printf " -------------------------------------------------------------------------------\n";

  (* ------------ *
   * Output Stats *
   * ------------ *)

  (* Get final statistics on how the solve progressed *)
  let nst = ARKStep.get_num_steps arkode_mem
  and nst_a = ARKStep.get_num_step_attempts arkode_mem
  and netf = ARKStep.get_num_err_test_fails arkode_mem
  and nfe, nfi = ARKStep.get_num_rhs_evals arkode_mem
  in
  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total number of error test failures = %d\n" netf;
  printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe nfi;

  if implicit then begin
    let nni     = ARKStep.get_num_nonlin_solv_iters arkode_mem
    and ncfn    = ARKStep.get_num_nonlin_solv_conv_fails arkode_mem
    and nsetups = ARKStep.get_num_lin_solv_setups arkode_mem
    and nje     = ARKStep.Dls.get_num_jac_evals arkode_mem
    and nfeLS   = ARKStep.Dls.get_num_lin_rhs_evals arkode_mem
    in
    printf "   Total number of Newton iterations = %d\n" nni;
    printf "   Total number of linear solver convergence failures = %d\n" ncfn;
    printf "   Total linear solver setups = %d\n" nsetups;
    printf "   Total number of Jacobian evaluations = %d\n" nje;
    printf "   Total RHS evals for setting up the linear system = %d\n" nfeLS
  end;

  if relax then begin
    let nre    = ARKStep.Relax.get_num_fn_evals arkode_mem
    and nrje   = ARKStep.Relax.get_num_jac_evals arkode_mem
    and nrf    = ARKStep.Relax.get_num_fails arkode_mem
    and nrbf   = ARKStep.Relax.get_num_bound_fails arkode_mem
    and nrnlsf = ARKStep.Relax.get_num_solve_fails arkode_mem
    and nrnlsi = ARKStep.Relax.get_num_solve_iters arkode_mem
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

