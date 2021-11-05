(* ------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2020.
 *---------------------------------------------------------------
 * Based a linear example program by Rujeko Chinomona @ SMU.
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
 * This example simulates an ODE system with 3 components,
 * Y = [u,v,w], given by the equations,
 *
 *   du/dt =  100v+w
 *   dv/dt = -100u
 *   dw/dt = -w+u
 *
 * for t in the interval [0.0, 2.0] with intial conditions
 * u(0)=9001/10001, v(0)=-1e-5/10001, and w(0)=1000. In this problem
 * the slow (w) and fast (u and v) components depend on one another.
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
  let c1 = 100.0 in
  let u = y.{0} in
  let v = y.{1} in
  ydot.{0} <- c1*.v;
  ydot.{1} <- -.c1*.u;
  ydot.{2} <- u

(* fs routine to compute the slow portion of the ODE RHS. *)
let fs t (y : RealArray.t) (ydot : RealArray.t) =
  let w = y.{2} in
  ydot.{0} <- w;
  ydot.{1} <- 0.0;
  ydot.{2} <- -.w

let main () =
  (* general problem parameters *)
  let t0 = 0.0 in      (* initial time *)
  let tf = 2.0 in      (* final time *)
  let dTout = 0.1 in   (* time between outputs *)
  let nt = Int.of_float (ceil(tf/.dTout)) in (* number of output times *)
  let hs = 0.001 in    (* slow step size *)
  let hf = 0.00002 in  (* fast step size *)

  (*
   * Initialization
   *)

  (* Set the initial contions *)
  let u0 = 9001.0/.10001.0
  and v0 = -.1.0e5/.10001.0
  and w0 = 1000.0
  in

  (* Initial problem output *)
  printf "\nTwo way coupling ODE test problem:\n";
  printf "    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n" u0 v0 w0;
  printf "    hs = %g,  hf = %g\n\n" hs hf;

  (* Create and initialize serial vector for the solution *)
  let ydata = RealArray.of_array [| u0; v0; w0; |] in
  let y = Nvector_serial.wrap ydata in (* Create serial vector for solution *)

  (* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, and the
     initial dependent variable vector y. *)
  let inner_arkode_mem = ARKStep.(init (explicit ff) default_tolerances t0 y) in
  ARKStep.set_erk_table_num inner_arkode_mem
    Arkode.ButcherTable.Knoth_Wolke_3_3;
  ARKStep.set_fixed_step inner_arkode_mem (Some hf);

  (* Call MRIStepCreate to initialize the MRI timestepper module and
     specify the right-hand side functions in y'=fs(t,y)+ff(t,y),
     the inital time T0, and the initial dependent variable vector y. *)
  (* Specify slow and fast step sizes *)
  let arkode_mem = MRIStep.(init
                              InnerStepper.(from_arkstep inner_arkode_mem)
                              default_tolerances
                              fs ~slowstep:hs t0 y) in

  (*
   * Integrate ODE
   *)

  (* Open output stream for results, output comment line *)
  let ufid = open_out "ark_twowaycouple_mri_solution.txt" in
  fprintf ufid "# t u v w maxerr\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16f %.16f %.16f %.16f\n"
          t0 ydata.{0} ydata.{1} ydata.{2};

  (* Main time-stepping loop: calls MRIStepEvolve to perform the
     integration, then prints results. Stops when the final time
     has been reached *)
  printf "        t           u           v           w\n";
  printf "   -----------------------------------------------\n";
  printf "  %10.6f  %10.6f  %10.6f  %10.6f\n"
            t0 ydata.{0} ydata.{1} ydata.{2};

  let rec loop iout tout =
    if iout = nt then ()
    else begin
      (* call integrator *)
      let t, _ = MRIStep.evolve_normal arkode_mem tout y in

      (* access/print solution and error *)
      printf "  %10.6f  %10.6f  %10.6f  %10.6f\n"
             t ydata.{0} ydata.{1} ydata.{2};
      fprintf ufid " %.16f %.16f %.16f %.16f\n"
             t ydata.{0} ydata.{1} ydata.{2};

      (* successful solve: update time *)
      loop (iout+1) (min (tout+.dTout) tf)
    end
  in
  loop 0 (t0+.dTout);
  printf "   -----------------------------------------------\n";

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
