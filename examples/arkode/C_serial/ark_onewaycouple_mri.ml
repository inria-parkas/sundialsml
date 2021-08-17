(* ------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2020.
 *---------------------------------------------------------------
 * Based a linear test problem from Estep, Ginting, and Tavener,
 * "A Posteriori analysis of a multirate numerical method for
 * ordinary differential equations," 2012 and an example program by
 * Rujeko Chinomona @ SMU.
 * ------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ------------------------------------------------------------------
 * Example problem:
 *
 * This example simulates an ODE system with 3 components,
 * Y = [u,v,w], given by the equations,
 *
 *   du/dt = -50v
 *   dv/dt =  50u
 *   dw/dt = -w+u+v
 *
 * for t in the interval [0.0, 1.0] with intial conditions u(0)=1.0,
 * v(0)=0.0, and w(0)=2.0. In this problem the slow time scale (w)
 * depends on the fast components (u and v), but the fast components
 * are independent of the slow component. This system has the
 * analytic soltuion,
 *
 *   u(t) = cos(50t)
 *   v(t) = sin(50t)
 *   w(t) = 5051/2501*exp(-t) - 49/2501*cos(50t) + 51/2501*sin(50t)
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

(* ------------------------------
 * Functions called by the solver
 * ------------------------------*)

(* ff routine to compute the fast portion of the ODE RHS. *)
let ff t (y : RealArray.t) (ydot : RealArray.t) =
  let w = y.{2} in
  ydot.{0} <- 0.0;
  ydot.{1} <- 0.0;
  ydot.{2} <- -.w

(* fs routine to compute the slow portion of the ODE RHS. *)
let fs t (y : RealArray.t) (ydot : RealArray.t) =
  let c1 = 50.0 in
  let u = y.{0} in
  let v = y.{1} in
  ydot.{0} <- -.c1*.v;
  ydot.{1} <- c1*.u;
  ydot.{2} <- u+.v

(* ------------------------------------
 * Private solution and error functions
 * ------------------------------------*)

(* function to compute the analytic solution of the ODE *)
let ans t (ytrue : Nvector_serial.t) =
  let ytrue = Nvector_serial.unwrap ytrue in
  let c1 = 50.0 in
  let c2 = 5051.0/.2501.0 in
  let c3 = 49.0/.2501.0 in
  let c4 = 51.0/.2501.0 in
  ytrue.{0} <- cos(c1*.t);
  ytrue.{1} <- sin(c1*.t);
  ytrue.{2} <- c2*.exp(-.t) -. c3*.cos(c1*.t) +. c4*.sin(c1*.t)

(* function to compute the max error in the solution *)
let err (y : Nvector_serial.t) (ytrue : Nvector_serial.t) =
  (* compute the error and store it in ytrue *)
  Nvector_serial.Ops.linearsum 1.0 y (-.1.0) ytrue ytrue;
  (* compute the max norm of the error *)
  Nvector_serial.Ops.maxnorm ytrue

let main () =
  (* general problem parameters *)
  let t0 = 0.0 in     (* initial time *)
  let tf = 1.0 in     (* final time *)
  let dTout = 0.1 in  (* time between outputs *)
  let nt = Int.of_float (ceil(tf/.dTout)) in (* number of output times *)
  let hs = 0.001 in   (* slow step size *)
  let hf = 0.0001 in  (* fast step size *)

  (*
   * Initialization
   *)

  (* Set the initial contions *)
  let u0 = 1.0
  and v0 = 0.0
  and w0 = 2.0
  in

  (* Initial problem output *)
  printf "\nOne way coupling ODE test problem:\n";
  printf "    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n" u0 v0 w0;
  printf "    hs = %g,  hf = %g\n\n" hs hf;

  (* Create and initialize serial vector for the solution *)
  let ydata = RealArray.of_array [| u0; v0; w0; |] in
  let y = Nvector_serial.wrap ydata in (* Create serial vector for solution *)

  (* Create serial vector for the analytic solution *)
  let ytrue = Nvector_serial.Ops.clone y in

  (* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, and the
     initial dependent variable vector y. *)
  let inner_arkode_mem = ARKStep.(init
                                   (explicit ff)
                                   Arkode.default_tolerances
                                   t0
                                   y)
  in
  ARKStep.set_erk_table_num inner_arkode_mem
    Arkode.ButcherTable.Knoth_Wolke_3_3;
  ARKStep.set_fixed_step inner_arkode_mem (Some hf);

  (* Call MRIStepCreate to initialize the MRI timestepper module and
     specify the right-hand side functions in y'=fs(t,y)+ff(t,y),
     the inital time T0, and the initial dependent variable vector y. *)
  (* Specify slow and fast step sizes *)
  let arkode_mem = MRIStep.(init inner_arkode_mem fs hs t0 y) in

  (*
   * Integrate ODE
   *)

  (* Open output stream for results, output comment line *)
  let ufid = open_out "ark_onewaycouple_mri_solution.txt" in
  fprintf ufid "# t u v w maxerr\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16f %.16f %.16f %.16f %.16f\n"
          t0 ydata.{0} ydata.{1} ydata.{2} 0.;

  (* Main time-stepping loop: calls MRIStepEvolve to perform the
     integration, then prints results. Stops when the final time
     has been reached *)
  printf "        t           u           v           w       max err\n";
  printf "   ----------------------------------------------------------\n";
  printf "  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n"
            t0 ydata.{0} ydata.{1} ydata.{2} 0.;

  let rec loop iout tout =
    if iout = nt then ()
    else begin
      (* call integrator *)
      let t, _ = MRIStep.solve_normal arkode_mem tout y in

      (* compute the analytic solution *)
      ans t ytrue;

      (* compute the error compared to the analytic solution *)
      let error = err y ytrue in

      (* access/print solution and error *)
      printf "  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n"
             t ydata.{0} ydata.{1} ydata.{2} error;
      fprintf ufid " %.16f %.16f %.16f %.16f %.16f\n"
             t ydata.{0} ydata.{1} ydata.{2} error;

      (* successful solve: update time *)
      loop (iout+1) (min (tout+.dTout) tf)
    end
  in
  loop 0 (t0+.dTout);
  printf "   ----------------------------------------------------------\n";

  (*
   * Finalize
   *)

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
