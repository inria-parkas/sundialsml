(*--------------------------------------------------------------- {{{
 * Programmer(s): Daniel R. Reynolds @ SMU
 * --------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Oct 2021.
 *---------------------------------------------------------------
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
 * The following test simulates the Robertson problem,
 * corresponding to the kinetics of an autocatalytic reaction.
 * This is an ODE system with 3 components, Y = [u,v,w], satisfying
 * the equations,
 *    du/dt = -0.04*u + 1e4*v*w
 *    dv/dt = 0.04*u - 1e4*v*w - 3e7*v^2
 *    dw/dt = 3e7*v^2
 * for t in the interval [0.0, 1e11], with initial conditions
 * Y0 = [1,0,0].
 *
 * This program solves the problem with one of the solvers, ERK,
 * DIRK or ARK.  For DIRK and ARK, implicit subsystems are solved
 * using a Newton iteration with the dense SUNLinearSolver, and a
 * user-supplied Jacobian routine. The constraint y_i >= 0 is
 * posed for all components.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *--------------------------------------------------------------- }}} *)

open Sundials
module ARKStep = Arkode.ARKStep

let printf = Printf.printf
let fprintf = Printf.fprintf
let eprintf = Printf.eprintf

(*-------------------------------
 * Functions called by the solver
 *-------------------------------*)

(* f routine to compute the ODE RHS function f(t,y). *)
let f _ (y : RealArray.t) (ydot : RealArray.t) =
  let u = y.{0} in   (* access current solution *)
  let v = y.{1} in
  let w = y.{2} in

  (* Fill in ODE RHS function *)
  ydot.{0} <- -0.04 *. u +. 1.e4 *. v *. w;
  ydot.{1} <- 0.04 *. u -. 1.e4 *. v *. w -. 3.e7 *. v *. v;
  ydot.{2} <- 3.e7 *. v *. v

(* Jacobian routine to compute J(t,y) = df/dy. *)
let jac { ARKStep.jac_y = (y : RealArray.t); _ } jac =
  let v = y.{1} in   (* access current solution *)
  let w = y.{2} in
  Matrix.Dense.set_to_zero jac; (* initialize Jacobian to zero *)

  (* Fill in the Jacobian of the ODE RHS function *)
  let open Matrix.Dense in
  set jac 0 0 (-0.04);
  set jac 0 1 (1.e4 *. w);
  set jac 0 2 (1.e4 *. v);

  set jac 1 0 0.04;
  set jac 1 1 (-1.e4 *. w -. 6.e7 *. v);
  set jac 1 2 (-1.e4 *. v);

  set jac 2 1 (6.e7 *. v)

(*-------------------------------
 * Private helper functions
 *-------------------------------*)

(* compare the solution at the final time 1e11s to a reference solution computed
   using a relative tolerance of 1e-8 and absoltue tolerance of 1e-14 *)
let check_ans (y : Nvector_serial.t) _ rtol atol =

  (* create reference solution and error weight vectors *)
  let refv = Nvector.clone y in
  let refd = Nvector.unwrap refv in
  let ewtv = Nvector.clone y in

  (* set the reference solution data *)
  refd.{0} <- 2.0833403356917897e-08;
  refd.{1} <- 8.1470714598028223e-14;
  refd.{2} <- 9.9999997916651040e-01;

  (* compute the error weight vector *)
  let open Nvector_serial in
  Ops.abs refv ewtv;
  Ops.scale rtol ewtv ewtv;
  Ops.addconst ewtv atol ewtv;
  if Ops.min ewtv <= 0.0 then
    failwith "SUNDIALS_ERROR: check_ans failed - ewt <= 0";
  Ops.inv ewtv ewtv;

  (* compute the solution error *)
  Ops.linearsum 1.0 y (-1.0) refv refv;
  let err = Ops.wrmsnorm refv ewtv in

  (* is the solution within the tolerances? *)
  if not (err < 1.0) then
    eprintf "\nSUNDIALS_WARNING: check_ans error=%g\n\n" err

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in                               (* initial time *)
  let tf = 1.e11 in                             (* final time *)
  let dTout = (tf -. t0) /. 100.0 in            (* time between outputs *)
  let nt = int_of_float (ceil (tf /. dTout)) in (* number of output times *)
  let neq = 3 in                                (* number of dependent vars. *)

  (* set up the initial conditions, tolerances, initial time step size *)
  let u0 = 1.0 in
  let v0 = 0.0 in
  let w0 = 0.0 in
  let reltol = 1.e-3 in
  let abstol = 1.e-7 in
  let h0 = 1.e-4 *. reltol in

  (* Initial problem output *)
  printf "\nRobertson ODE test problem:\n";
  printf "    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n" u0 v0 w0;

  (* Create serial vector for solution *)
  (* Set initial conditions into y *)
  let yd = RealArray.of_array [| u0; v0; w0 |] in
  let y  = Nvector_serial.wrap yd in

  (* Set constraints to all 1's for nonnegative solution values. *)
  let constraints = Nvector_serial.make neq 1.0 in

  (* Initialize dense matrix data structure and solver *)
  let a = Matrix.dense neq in
  let ls = LinearSolver.Direct.dense y a in

  (* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  (* Attach matrix and linear solver *)
  (* Set the Jacobian routine *)
 (* Specify tolerances *)
  let arkode_mem = ARKStep.(init
                             (implicit ~lsolver:(Dls.solver ~jac ls) f)
                             (SStolerances (reltol, abstol))
                             t0 y)
  in
  ARKStep.set_init_step arkode_mem h0;            (* Set custom initial step *)
  ARKStep.set_max_err_test_fails arkode_mem 20;   (* Increase max error test fails *)
  ARKStep.set_max_nonlin_iters arkode_mem 8;      (* Increase max nonlin iters  *)
  ARKStep.set_nonlin_conv_coef arkode_mem 1.e-7;  (* Set nonlinear convergence coeff. *)
  ARKStep.set_max_num_steps arkode_mem 100000;    (* Increase max num steps *)
  ARKStep.(set_predictor_method arkode_mem MaximumOrderPredictor);
                                                  (* Specify maximum-order predictor *)
  ARKStep.set_constraints arkode_mem constraints; (* Set constraints *)

  (* Open output stream for results, output comment line *)
  let ufid = open_out "solution.txt" in
  fprintf ufid "# t u v w\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16e %.16e %.16e %.16e\n" t0 yd.{0} yd.{1} yd.{2};

  (* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let t = t0 in
  let tout = ref (t0 +. dTout) in
  printf "        t           u           v           w\n";
  printf "   --------------------------------------------------\n";
  printf "  %10.3e  %12.5e  %12.5e  %12.5e\n" t yd.{0} yd.{1} yd.{2};
  for _ = 0 to nt - 1 do
    let t, flag = ARKStep.evolve_normal arkode_mem !tout y in     (* call integrator *)

    printf "  %10.3e  %12.5e  %12.5e  %12.5e\n" t yd.{0} yd.{1} yd.{2}; (* access/print solution *)
    fprintf ufid " %.16e %.16e %.16e %.16e\n" t yd.{0} yd.{1} yd.{2};
    (match flag with
     | ARKStep.Success ->           (* successful solve: update time *)
        tout := min tf (!tout +. dTout);
     | _ ->                         (* unsuccessful solve: break *)
      eprintf "Solver failure, stopping integration\n")
  done;
  printf "   --------------------------------------------------\n";
  close_out ufid;

  (* Print some final statistics *)
  let open ARKStep in
  let nst      = get_num_steps arkode_mem in
  let nst_a    = get_num_step_attempts arkode_mem in
  let nfe, nfi = get_num_rhs_evals arkode_mem in
  let nsetups  = get_num_lin_solv_setups arkode_mem in
  let netf     = get_num_err_test_fails arkode_mem in
  let nni      = get_num_nonlin_solv_iters arkode_mem in
  let ncfn     = get_num_nonlin_solv_conv_fails arkode_mem in
  let nje      = Dls.get_num_jac_evals arkode_mem in
  let nfeLS    = Dls.get_num_lin_rhs_evals arkode_mem in
  let nctf     = get_num_constr_fails arkode_mem in

  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe nfi;
  printf "   Total linear solver setups = %d\n" nsetups;
  printf "   Total RHS evals for setting up the linear system = %d\n" nfeLS;
  printf "   Total number of Jacobian evaluations = %d\n" nje;
  printf "   Total number of Newton iterations = %d\n" nni;
  printf "   Total number of nonlinear solver convergence failures = %d\n" ncfn;
  printf "   Total number of error test failures = %d\n" netf;
  printf "   Total number of constraint test failures = %d\n" nctf;

  (* check the solution error *)
  check_ans y t reltol abstol

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

