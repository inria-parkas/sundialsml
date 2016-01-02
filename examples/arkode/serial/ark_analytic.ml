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
 * This program solves the problem with the DIRK method,
 * Newton iteration with the ARKDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*)

module RealArray = Sundials.RealArray
let printf = Printf.printf
let fprintf = Printf.fprintf

(* Functions called by the solver *)

(* f routine to compute the ODE RHS function f(t,y). *)
let f lamda t y ydot =
  let u = y.{0} in (* access current solution value *)
  (* fill in the RHS function *)
  ydot.{0} <- lamda*.u +. 1.0/.(1.0+.t*.t) -. lamda*.atan(t)

(* Jacobian routine to compute J(t,y) = df/dy. *)
let jac lamda _ j =
  (* Fill in Jacobian of f: "DENSE_ELEM" accesses the (0,0) entry of J *)
  Dls.DenseMatrix.set j 0 0 lamda

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0     = 0.0 in      (* initial time *)
  let tf     = 10.0 in     (* final time *)
  let dTout  = 1.0 in      (* time between outputs *)
  let neq    = 1 in        (* number of dependent vars. *)
  let reltol = 1.0e-6 in   (* tolerances *)
  let abstol = 1.0e-10 in
  let lamda  = -100.0 in   (* stiffness parameter *)

  (* Initial diagnostics output *)
  printf "\nAnalytical ODE test problem:\n";
  printf "    lamda = %g\n"     lamda;
  printf "   reltol = %.1e\n"   reltol;
  printf "   abstol = %.1e\n\n" abstol;

  (* Initialize data structures *)
  let data = RealArray.make neq 0.0 in (* Set initial conditions *)
  let y = Nvector_serial.wrap data in  (* Create serial vector for solution *)

  (* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time t0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  let arkode_mem =
    Arkode.init
      (Arkode.Implicit
        (f lamda,
         Arkode.Newton (Arkode.Dls.dense ~jac:(jac lamda) ()),
         Arkode.Linear false))
      (Arkode.SStolerances (reltol, abstol))
      t0
      y
  in
  (* Open output stream for results, output comment line *)
  let ufid = open_out "solution.txt" in
  fprintf ufid "# t u\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16e %.16e\n" t0 data.{0};  

  (* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let t = ref t0 in
  let tout = ref (t0 +. dTout) in
  printf "        t           u\n";
  printf "   ---------------------\n";
  try
    while (tf -. !t > 1.0e-15) do
      (* call integrator *)
      let t', flag = Arkode.solve_normal arkode_mem !tout y in
      t := t';
      printf "  %10.6f  %10.6f\n" t' data.{0};      (* access/print solution *)
      fprintf ufid " %.16e %.16e\n" t' data.{0};  
      (* successful solve: update time *)
      tout := !tout +. dTout;
      if !tout > tf then tout := tf;
    done
  with _ -> (* unsuccessful solve: break *)
            fprintf stderr "Solver failure, stopping integration\n";
  printf "   ---------------------\n";
  close_out ufid;

  (* Get/print some final statistics on how the solve progressed *)
  let nst      = Arkode.get_num_steps arkode_mem in
  let nst_a    = Arkode.get_num_step_attempts arkode_mem in
  let nfe, nfi = Arkode.get_num_rhs_evals arkode_mem in
  let nsetups  = Arkode.get_num_lin_solv_setups arkode_mem in
  let netf     = Arkode.get_num_err_test_fails arkode_mem in
  let nni      = Arkode.get_num_nonlin_solv_iters arkode_mem in
  let ncfn     = Arkode.get_num_nonlin_solv_conv_fails arkode_mem in
  let nje      = Arkode.Dls.get_num_jac_evals arkode_mem in
  let nfeLS    = Arkode.Dls.get_num_rhs_evals arkode_mem in

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
  with Not_found | Failure "int_of_string" -> 1
let gc_at_end =
  try int_of_string (Unix.getenv "GC_AT_END") <> 0
  with Not_found | Failure "int_of_string" -> false
let gc_each_rep =
  try int_of_string (Unix.getenv "GC_EACH_REP") <> 0
  with Not_found | Failure "int_of_string" -> false

(* Entry point *)
let _ =
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
