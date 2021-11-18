(*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *-----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Nov 2021.
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with analytical
 * solution,
 *    dy/dt = lamda*y + 1/(1+t^2) - lamda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lamda".  The value of lamda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the DIRK method and a
 * custom 'matrix-embedded' SUNLinearSolver. Output is printed
 * every 1.0 units of time (10 total).  Run statistics (optional
 * outputs) are printed at the end.
 *-----------------------------------------------------------------*)

open Sundials
module ARKStep = Arkode.ARKStep

let printf = Printf.printf

(*-------------------------------
 * Functions called by the solver
 *-------------------------------*)

(* f routine to compute the ODE RHS function f(t,y). *)
let f rdata t y ydot =
  let lambda = rdata.{0} in (* set shortcut for stiffness parameter *)
  let u = y.{0} in          (* access current solution value *)

  (* fill in the RHS function: "NV_Ith_S" accesses the 0th entry of ydot *)
  ydot.{0} <- lambda *. u +. 1.0 /. (1.0 +. t *. t) -. lambda *. atan t

(*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*)

type matrix_embedded_ls_content = {
  rdata : RealArray.t;
  mutable arkode_mem : (RealArray.t, Nvector_serial.kind) ARKStep.session option;
}

(* linear solve routine *)
let matrix_embedded_ls_solve content () x b tol =
  let { rdata; arkode_mem } = content in
  match arkode_mem with
  | None ->
      Printf.eprintf "internal error: linear solver not properly configured\n";
      exit(-1)
  | Some arkode_mem ->
      (* retrieve implicit system data from ARKStep *)
      let { ARKStep.gamma; _ } = ARKStep.get_nonlin_system_data arkode_mem in
      (* extract stiffness parameter from user_data *)
      let lambda = rdata.{0} in
      (* perform linear solve: (1-gamma*lamda)*x = b *)
      x.{0} <- b.{0} /. (1.0 -. gamma *. lambda)

(* constructor *)
let matrix_embedded_ls rdata =
  let content = { rdata; arkode_mem = None } in
  (fun session -> content.arkode_mem <- Some session),
  LinearSolver.Custom.(
    make_without_matrix (make_ops ~solver_type:MatrixEmbedded
                                  ~solve:matrix_embedded_ls_solve ())
                        content)

(* check the computed solution *)
let check_ans y t rtol atol =
  (* compute solution error *)
  let ans = atan(t) in
  let ewt = 1.0 /. (rtol *. abs_float ans +. atol) in
  let err = ewt *. abs_float(y.{0} -. ans) in

  (* is the solution within the tolerances? *)
  if err >= 1.0 then begin
    printf "\nSUNDIALS_WARNING: check_ans error=%g\n\n" err;
    false
  end else true

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in         (* initial time *)
  let tf = 10.0 in        (* final time *)
  let dTout = 1.0 in      (* time between outputs *)
  let neq = 1 in          (* number of dependent vars. *)
  let reltol = 1.0e-6 in  (* tolerances *)
  let abstol = 1.0e-10 in
  let lambda  = -100.0 in (* stiffness parameter *)

  (* Initial diagnostics output *)
  printf "\nAnalytical ODE test problem:\n";
  printf "    lamda = %g\n"     lambda;
  printf "   reltol = %.1e\n"   reltol;
  printf "   abstol = %.1e\n\n" abstol;

  (* Initialize data structures *)
  let y = Nvector_serial.make neq 0.0 in (* Create serial vector for solution *)
                                         (* Specify initial condition *)

  (* Pass lamda to user functions *)
  let rdata = RealArray.of_list [ lambda ] in
  (* Initialize custom matrix-embedded linear solver *)
  let set_ls_session, ls = matrix_embedded_ls rdata in
  (* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  (* Specify tolerances *)
  (* Attach linear solver *)
  (* Specify linearly implicit RHS, with non-time-dependent Jacobian *)
  let arkode_mem =
    ARKStep.(init (implicit (f rdata)
                            ~lsolver:(matrix_embedded_solver ls)
                            ~linearity:(Linear false))
                  (SStolerances (reltol, abstol))
                  t0 y)
  in
  set_ls_session arkode_mem;

  (* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached. *)
  printf "        t           u\n";
  printf "   ---------------------\n";
  let ydata = Nvector.unwrap y in
  let rec loop t tout =
    if tf -. t <= 1.0e-15 then t
    else begin
      let t, _ = ARKStep.evolve_normal arkode_mem tout y in (* call integrator *)
      printf "  %10.6f  %10.6f\n" t ydata.{0};      (* access/print solution *)
      loop t (min tf (tout +. dTout))
    end
  in
  let t =
    try loop t0 (t0 +. dTout);
    with _ -> (Format.eprintf "Solver failure, stopping integration\n"; 0.0)
  in
  printf "   ---------------------\n";

  (* Get/print some final statistics on how the solve progressed *)
  let nst = ARKStep.get_num_steps arkode_mem in
  let nst_a = ARKStep.get_num_step_attempts arkode_mem in
  let nfe, nfi = ARKStep.get_num_rhs_evals arkode_mem in
  let nsetups = ARKStep.get_num_lin_solv_setups arkode_mem in
  let netf = ARKStep.get_num_err_test_fails arkode_mem in
  let nni = ARKStep.get_num_nonlin_solv_iters arkode_mem in
  let ncfn = ARKStep.get_num_nonlin_solv_conv_fails arkode_mem in
  let nje = ARKStep.Dls.get_num_jac_evals arkode_mem in
  let nfeLS = ARKStep.Dls.get_num_lin_rhs_evals arkode_mem in

  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe nfi;
  printf "   Total linear solver setups = %d\n" nsetups;
  printf "   Total RHS evals for setting up the linear system = %d\n" nfeLS;
  printf "   Total number of Jacobian evaluations = %d\n" nje;
  printf "   Total number of Newton iterations = %d\n" nni;
  printf "   Total number of linear solver convergence failures = %d\n" ncfn;
  printf "   Total number of error test failures = %d\n\n" netf;

  (* check the solution error *)
  ignore (check_ans ydata t reltol abstol)

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

