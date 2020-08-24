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
 *     dy/dt = (t+1)*exp(-y)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 * This has analytical solution
 *      y(t) = log(0.5*t^2 + t + 1)
 *
 * This program solves the problem with the ERK method.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*)

open Sundials
module ERKStep = Arkode.ERKStep

let printf = Printf.printf
let fprintf = Printf.fprintf

(* Functions called by the solver *)

(* f routine to compute the ODE RHS function f(t,y). *)
let f t (y : RealArray.t) (ydot : RealArray.t) =
  ydot.{0} <- (t+.1.0)*.exp(-.y.{0})

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0     = 0.0 in      (* initial time *)
  let tf     = 10.0 in     (* final time *)
  let dTout  = 1.0 in      (* time between outputs *)
  let neq    = 1 in        (* number of dependent vars. *)
  let reltol = 1.0e-6 in   (* tolerances *)
  let abstol = 1.0e-10 in

  (* Initial problem output *)
  printf "\nAnalytical ODE test problem:\n";
  printf "   reltol = %.1e\n"   reltol;
  printf "   abstol = %.1e\n\n" abstol;

  (* Initialize data structures *)
  let data = RealArray.make neq 0.0 in (* Set initial conditions *)
  let y = Nvector_serial.wrap data in  (* Create serial vector for solution *)

  (* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully explicit, we set f_U to NULL and f_E to f. *)
  let arkode_mem =
    ERKStep.(init (SStolerances (reltol, abstol)) f t0 y)
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
  (try
     while (tf -. !t > 1.0e-15) do
       (* call integrator *)
       let t', _ = ERKStep.solve_normal arkode_mem !tout y in
       t := t';
       printf "  %10.6f  %10.6f\n" t' data.{0};      (* access/print solution *)
       fprintf ufid " %.16e %.16e\n" t' data.{0};
       (* successful solve: update time *)
       tout := min (!tout +. dTout) tf
     done
   with _ -> (* unsuccessful solve: break *)
             fprintf stderr "Solver failure, stopping integration\n");
  printf "   ---------------------\n";
  close_out ufid;

  (* Get/print some final statistics *)
  let open ERKStep in
  let nst      = get_num_steps arkode_mem in
  let nst_a    = get_num_step_attempts arkode_mem in
  let nfe      = get_num_rhs_evals arkode_mem in
  let netf     = get_num_err_test_fails arkode_mem in

  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total RHS evals = %d\n" nfe;
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
