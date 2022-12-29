(* ------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2020.
 * ------------------------------------------------------------------
 * Based an example program by Rujeko Chinomona @ SMU.
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
 * The following test simulates a simple 1D reaction-diffusion
 * equation,
 *
 *   y_t = k * y_xx + y^2 * (1-y)
 *
 * for t in [0, 3], x in [0, L] with boundary conditions,
 *
 *   y_x(0,t) = y_x(L,t) = 0
 *
 * and initial condition,
 *
 *   y(x,0) = (1 + exp(lambda*(x-1))^(-1),
 *
 * with parameter k = 1e-4/ep, lambda = 0.5*sqrt(2*ep*1e4),
 * ep = 1e-2, and L = 5.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
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

type userdata = {
  n   : int;   (* number of intervals   *)
  k   : float; (* diffusion coefficient *)
  dx  : float; (* mesh spacing          *)
  lam : float;
}

(* ff routine to compute the fast portion of the ODE RHS. *)
let ff { n } _ (y : RealArray.t) (ydot : RealArray.t) =
  (* iterate over domain, computing reaction term *)
  for i = 0 to n - 1 do
    ydot.{i} <- y.{i} *. y.{i} *. (1.0 -. y.{i})
  done

(* fs routine to compute the slow portion of the ODE RHS. *)
let fs { n; k; dx } _ (y : RealArray.t) (ydot : RealArray.t) =
  (* iterate over domain, computing diffusion term *)
  let c1 = k/.dx/.dx in
  let c2 = 2.0*.k/.dx/.dx in

  (* left boundary condition *)
  ydot.{0} <- c2*.(y.{1} -. y.{0});

  (* interior points *)
  for i=1 to n-2 do
    ydot.{i} <- c1*.y.{i-1} -. c2*.y.{i} +. c1*.y.{i+1}
  done;

  (* right boundary condition *)
  ydot.{n-1} <- c2*.(y.{n-2} -. y.{n-1})

(* -----------------------------------------
 * Private function to set initial condition
 * -----------------------------------------*)

let set_initial_condition { n; dx; lam } y =
  let y = Nvector_serial.unwrap y in
  (* set initial condition *)
  for i = 0 to n-1 do
    y.{i} <- 1.0/.(1. +. exp(lam*.((float i)*.dx-.1.0)))
  done

(* Main Program *)
let main () =

  (* general problem parameters *)
  let t0 = 0.0 in     (* initial time *)
  let tf = 3.0 in     (* final time *)
  let dTout = 0.1 in  (* time between outputs *)
  let nt = int_of_float(ceil(tf/.dTout)) in (* number of output times *)
  let hs = 0.001 in   (* slow step size *)
  let hf = 0.00002 in (* fast step size *)

  let l = 5.0 in      (* domain length *)
  let n = 1001 in     (* number of mesh points *)
  let ep = 1e-2 in

  (*
   * Initialization
   *)

  (* allocate and fill user data structure *)
  let udata = {
    n;
    dx = l /. (float n -. 1.0);
    k = 1e-4/.ep;
    lam = 0.5*.sqrt(2.0 *. ep *. 1e4)
  } in

  (* Initial problem output *)
  printf "\n1D reaction-diffusion PDE test problem:\n";
  printf "  N = %d\n" udata.n;
  printf "  diffusion coefficient:  k = %g\n" udata.k;

  (* Create and initialize serial vector for the solution *)
  let y = Nvector_serial.make n 0. in (* Create serial vector for solution *)
  set_initial_condition udata y;

  (* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, and the
     initial dependent variable vector y. *)
  let inner_arkode_mem = ARKStep.(init (explicit (ff udata)) default_tolerances t0 y) in
  ARKStep.set_erk_table_num inner_arkode_mem
    Arkode.ButcherTable.Knoth_Wolke_3_3;
  ARKStep.set_fixed_step inner_arkode_mem (Some hf);

  (* Call MRIStepCreate to initialize the MRI timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  (* Pass udata to user functions *)
  (* Specify slow and fast step sizes *)
  let arkode_mem = MRIStep.(init
                              (explicit (fs udata))
                              default_tolerances
                              InnerStepper.(from_arkstep inner_arkode_mem)
                              ~slowstep:hs t0 y)
  in
  (* Increase max num steps  *)
  MRIStep.set_max_num_steps arkode_mem 10000;

  (*
   * Integrate ODE
   *)

  (* output mesh to disk *)
  let fid = open_out "heat_mesh.txt" in
  for i=0 to n-1 do
    fprintf fid "  %.16e\n" (udata.dx*.float i)
  done;

  (* Open output stream for results, access data array *)
  let ufid = open_out "heat1D.txt" in
  let data = Nvector_serial.unwrap y in

  (* output initial condition to disk *)
  for i=0 to n-1 do
    fprintf ufid " %.16e" data.{i}
  done;
  fprintf ufid "\n";

  (* Main time-stepping loop: calls MRIStepEvolve to perform the integration, then
     prints results. Stops when the final time has been reached *)
  printf "        t      ||u||_rms\n";
  printf "   -------------------------\n";
  printf "  %10.6f  %10.6f\n" t0
    (sqrt((Nvector_serial.Ops.dotprod y y)/.float n));

  let dTout = (tf-.t0)/.(float nt) in
  let rec loop iout tout =
    if iout = nt then ()
    else begin
      (* call integrator *)
      let t, _ = MRIStep.evolve_normal arkode_mem tout y in

      (* print solution stats and output results to disk *)
      printf "  %10.6f  %10.6f\n" t
        (sqrt((Nvector_serial.Ops.dotprod y y)/.float n));

      for i=0 to n-1 do
        fprintf ufid " %.16e" data.{i}
      done;
      fprintf ufid "\n";

      (* successful solve: update output time *)
      loop (iout+1) (min (tout+.dTout) tf)
    end
  in
  loop 0 (t0+.dTout);
  printf "   -------------------------\n";

  (* Print some final statistics *)
  let nsts = MRIStep.get_num_steps arkode_mem in
  let nfse, _ = MRIStep.get_num_rhs_evals arkode_mem in
  let nstf = ARKStep.get_num_steps inner_arkode_mem in
  let nff, _ = ARKStep.get_num_rhs_evals inner_arkode_mem in
  if Sundials_impl.Version.lt620 then begin
    printf "\nFinal Solver Statistics:\n";
    printf "   Steps: nsts = %d, nstf = %d\n" nsts nstf;
    printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfse nff
  end else begin
    printf "\nFinal Slow Statistics:\n";
    flush stdout;
    MRIStep.print_all_stats arkode_mem Logfile.stdout Sundials.OutputTable;
    Logfile.flush Logfile.stdout;
    printf "\nFinal Fast Statistics:\n";
    flush stdout;
    ARKStep.print_all_stats inner_arkode_mem Logfile.stdout Sundials.OutputTable;
    Logfile.flush Logfile.stdout;

    let fid = Logfile.openfile "ark_reaction_diffusion_mri_slow_stats.csv" in
    MRIStep.print_all_stats arkode_mem fid Sundials.OutputCSV;
    Logfile.close fid;

    let fid = Logfile.openfile "ark_reaction_diffusion_mri_fast_stats.csv" in
    ARKStep.print_all_stats inner_arkode_mem fid Sundials.OutputCSV;
    Logfile.close fid
  end

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
