(*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jan 2016.
 *---------------------------------------------------------------
 * Copyright (c) 2015, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is an ODE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions
 * Y0 = [u0,v0,w0].
 *
 * We have 3 different testing scenarios:
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 *    Here, all three components exhibit a rapid transient change
 *    during the first 0.2 time units, followed by a slow and
 *    smooth evolution.
 *
 * Test 2:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
 *    Here, w experiences a fast initial transient, jumping 0.5
 *    within a few steps.  All values proceed smoothly until
 *    around t=6.5, when both u and v undergo a sharp transition,
 *    with u increaseing from around 0.5 to 5 and v decreasing
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * Test 3:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3,  ep=5.0e-4
 *    Here, all components undergo very rapid initial transients
 *    during the first 0.3 time units, and all then proceed very
 *    smoothly for the remainder of the simulation.
 *
 * This file is hard-coded to use test 2.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the ARKDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *-----------------------------------------------------------------*)

open Sundials
module ARKStep = Arkode.ARKStep

let compat2_3 =
  match Config.sundials_version with
  | 2,_,_ -> true
  | 3,_,_ -> true
  | _ -> false

let sungte500 =
  let n, _, _ = Config.sundials_version in
  n >= 5

let printf = Printf.printf
let fprintf = Printf.fprintf

(* f routine to compute the ODE RHS function f(t,y). *)
let f rdata t (y : RealArray.t) (ydot : RealArray.t) =
  let a  = rdata.(0) in     (* access data entries *)
  let b  = rdata.(1) in
  let ep = rdata.(2) in
  let u = y.{0} in          (* access solution values *)
  let v = y.{1} in
  let w = y.{2} in

  (* fill in the RHS function *)
  ydot.{0} <- a -. (w+.1.0)*.u +. v*.u*.u;
  ydot.{1} <- w*.u -. v*.u*.u;
  ydot.{2} <- (b-.w)/.ep -. w*.u

(* Jacobian routine to compute J(t,y) = df/dy. *)
let jac rdata { ARKStep.jac_y = y } j =
  let ep = rdata.(2) in   (* access data entries *)
  let u = y.{0} in        (* access solution values *)
  let v = y.{1} in
  let w = y.{2} in

  (* fill in the Jacobian *)
  let open Matrix.Dense in
  set j 0 0 (-.(w+.1.0) +. 2.0*.u*.v);
  set j 0 1 (u*.u);
  set j 0 2 (-.u);

  set j 1 0 (w -. 2.0*.u*.v);
  set j 1 1 (-.u*.u);
  set j 1 2 (u);

  set j 2 0 (-.w);
  set j 2 1 (0.0);
  set j 2 2 (-.1.0/.ep -. u)

type test = {
    u0 : float;
    v0 : float;
    w0 : float;
    a  : float;
    b  : float;
    ep : float;
  }

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in                            (* initial time *)
  let tf = 10.0 in                           (* final time *)
  let dTout = 1.0 in                         (* time between outputs *)
  let nt = int_of_float (ceil(tf/.dTout)) in (* number of output times *)
  let test = 2 in                            (* test problem to run *)
  let reltol = 1.0e-6 in                     (* tolerances *)
  let abstol = 1.0e-10 in

  (* general problem variables *)
  (* set up the test problem according to the desired test *)
  let testdata =
    match test with
    | 1 -> { u0 = 3.9;
             v0 = 1.1;
             w0 = 2.8;
             a  = 1.2;
             b  = 2.5;
             ep = 1.0e-5 }
    | 3 -> { u0 = 3.0;
             v0 = 3.0;
             w0 = 3.5;
             a  = 0.5;
             b  = 3.0;
             ep = 5.0e-4; }
    | _ -> { u0 = 1.2;
             v0 = 3.1;
             w0 = 3.0;
             a  = 1.0;
             b  = 3.5;
             ep = 5.0e-6; }
  in

  (* Initial problem output *)
  printf "\nBrusselator ODE test problem:\n";
  printf "    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n"
                                          testdata.u0 testdata.v0 testdata.w0;
  printf "    problem parameters:  a = %g,  b = %g,  ep = %g\n"
                                          testdata.a  testdata.b  testdata.ep;
  printf "    reltol = %.1e,  abstol = %.1e\n\n" reltol abstol;

  (* Initialize data structures *)
  let rdata = [| testdata.a; testdata.b; testdata.ep |] in
  let data = RealArray.of_array [|    (* Set initial conditions *)
                testdata.u0;
                testdata.v0;
                testdata.w0;
              |] in
  let y = Nvector_serial.wrap data in (* Create serial vector for solution *)

  (* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time t0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  let m = Matrix.dense 3 in
  let arkode_mem = ARKStep.(
    init (implicit (f rdata)
                   ~lsolver:Dls.(solver ~jac:(jac rdata) (dense y m)))
         (SStolerances (reltol, abstol))
         t0
         y
  ) in

  (* Specify stiff interpolant *)
  if sungte500 then ARKStep.(set_interpolant_type arkode_mem Lagrange);

  (* Open output stream for results, output comment line *)
  let ufid = open_out "solution.txt" in
  fprintf ufid "# t u v w\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16e %.16e %.16e %.16e\n" t0 data.{0} data.{1} data.{2};

  (* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let tout =ref (t0 +. dTout) in
  printf "        t           u           v           w\n";
  printf "   -------------------------------------------\n";
  (try
     if not compat2_3 then begin
       printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" 0. data.{0} data.{1} data.{2};
       fprintf ufid " %.16e %.16e %.16e %.16e\n" 0. data.{0} data.{1} data.{2}
     end;
     for iout=0 to nt-1 do
       (* call integrator *)
       let t, _ = ARKStep.evolve_normal arkode_mem !tout y in
       (* access/print solution *)
       printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" t data.{0} data.{1} data.{2};
       fprintf ufid " %.16e %.16e %.16e %.16e\n" t data.{0} data.{1} data.{2};
       (* successful solve: update time *)
       tout := min (!tout +. dTout) tf
     done
   with _ -> (* unsuccessful solve: break *)
             fprintf stderr "Solver failure, stopping integration\n");
  printf "   -------------------------------------------\n";
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

  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe nfi;
  printf "   Total linear solver setups = %d\n" nsetups;
  printf "   Total RHS evals for setting up the linear system = %d\n" nfeLS;
  printf "   Total number of Jacobian evaluations = %d\n" nje;
  printf "   Total number of Newton iterations = %d\n" nni;
  printf "   Total number of linear solver convergence failures = %d\n" ncfn;
  printf "   Total number of error test failures = %d\n\n" netf

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
