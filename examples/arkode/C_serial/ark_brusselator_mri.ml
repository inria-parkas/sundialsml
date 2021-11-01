(* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2020.
 *---------------------------------------------------------------
 * Based on ark_brusselator.c by Daniel R. Reynolds @ SMU
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics. This is an ODE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 *
 * for t in the interval [0.0, 2.0], with parameter values a=1,
 * b=3.5, and ep=1.0e-2. The initial conditions Y0 = [u0,v0,w0] are
 * u0=1.2, v0=3.1, and w0=3.
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 * ----------------------------------------------------------------*)

open Sundials
module ARKStep = Arkode.ARKStep
module MRIStep = Arkode.MRIStep

let printf = Printf.printf
let fprintf = Printf.fprintf

(* ff routine to compute the fast portion of the ODE RHS. *)
let ff rdata t (y : RealArray.t) (ydot : RealArray.t) =
  let b  = rdata.(1) in
  let ep = rdata.(2) in
  let w = y.{2} in
  (* fill in the RHS function *)
  ydot.{0} <- 0.0;
  ydot.{1} <- 0.0;
  ydot.{2} <- (b-.w)/.ep

(* fs routine to compute the slow portion of the ODE RHS. *)
let fs rdata t (y : RealArray.t) (ydot : RealArray.t) =
  let a  = rdata.(0) in
  let u = y.{0} in
  let v = y.{1} in
  let w = y.{2} in
  (* fill in the RHS function *)
  ydot.{0} <- a -. (w+.1.0)*.u +. v*.u*.u;
  ydot.{1} <- w*.u -. v*.u*.u;
  ydot.{2} <- -. w*.u

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in                            (* initial time *)
  let tf = 2.0 in                            (* final time *)
  let dTout = 0.1 in                         (* time between outputs *)
  let nt = int_of_float (ceil(tf/.dTout)) in (* number of output times *)
  let hs = 0.025 in                          (* slow step size *)
  let hf = 0.001 in                          (* fast step size *)

  (* general problem variables *)
  let a  = 1.0
  and b  = 3.5
  and ep = 1.0e-2
  and u0 = 1.2
  and v0 = 3.1
  and w0 = 3.0
  in

  (* Initial problem output *)
  printf "\nBrusselator ODE test problem:\n";
  printf "    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n" u0 v0 w0;
  printf "    problem parameters:  a = %g,  b = %g,  ep = %g\n" a  b  ep;
  printf "    hs = %g,  hf = %g\n\n" hs hf;

  (* Initialize data structures *)
  let rdata = [| a; b; ep |] in
  let data = RealArray.of_array [| u0; v0; w0; |] in
  let y = Nvector_serial.wrap data in (* Create serial vector for solution *)

  (* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, and the
     initial dependent variable vector y. *)
  let inner_arkode_mem = ARKStep.(init (explicit (ff rdata))
                                   default_tolerances t0 y)
  in
  ARKStep.set_erk_table_num inner_arkode_mem
    Arkode.ButcherTable.Knoth_Wolke_3_3;
  ARKStep.set_fixed_step inner_arkode_mem (Some hf);

  (* Call MRIStepCreate to initialize the MRI timestepper module and
     specify the right-hand side functions in y'=fs(t,y)+ff(t,y),
     the inital time T0, and the initial dependent variable vector y. *)
  let arkode_mem = MRIStep.(init
                              InnerStepper.(from_arkstep inner_arkode_mem)
                              default_tolerances
                              (fs rdata) ~slowstep:hs t0 y) in

  (* Open output stream for results, output comment line *)
  let ufid = open_out "ark_brusselator_mri_solution.txt" in
  fprintf ufid "# t u v w\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16e %.16e %.16e %.16e\n" t0 data.{0} data.{1} data.{2};

  (* Main time-stepping loop: calls MRIStepEvolve to perform the
     integration, then prints results. Stops when the final time
     has been reached *)
  let tout =ref (t0 +. dTout) in
  printf "        t           u           v           w\n";
  printf "   ----------------------------------------------\n";
  (try
     printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" 0. data.{0} data.{1} data.{2};
     fprintf ufid " %.16e %.16e %.16e %.16e\n" 0. data.{0} data.{1} data.{2};
     for iout=0 to nt-1 do
       (* call integrator *)
       let t, _ = MRIStep.solve_normal arkode_mem !tout y in
       (* access/print solution *)
       printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" t data.{0} data.{1} data.{2};
       fprintf ufid " %.16e %.16e %.16e %.16e\n" t data.{0} data.{1} data.{2};
       (* successful solve: update time *)
       tout := min (!tout +. dTout) tf
     done
   with _ -> (* unsuccessful solve: break *)
             fprintf stderr "Solver failure, stopping integration\n");
  printf "   ----------------------------------------------\n";
  close_out ufid;

  (* Print some final statistics *)
  let nsts = MRIStep.get_num_steps arkode_mem in
  let nfs = MRIStep.get_num_rhs_evals arkode_mem in
  let nstf = ARKStep.get_num_steps inner_arkode_mem in
  let nff, _ = ARKStep.get_num_rhs_evals inner_arkode_mem in

  printf "\nFinal Solver Statistics:\n";
  printf "   Steps: nsts = %d, nstf = %d\n" nsts nstf;
  printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfs nff

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
