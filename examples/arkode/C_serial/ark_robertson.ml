(*---------------------------------------------------------------
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
 * using a Newton iteration with the ARKDENSE dense linear solver,
 * and a user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *---------------------------------------------------------------*)

open Sundials
module ARKStep = Arkode.ARKStep

let printf = Printf.printf
let fprintf = Printf.fprintf

let sundials_270_or_later =
  match Config.sundials_version with
  | 2,5,_ | 2,6,_ -> false
  | _ -> true

(* f routine to compute the ODE RHS function f(t,y). *)
let f _ (y : RealArray.t) (ydot : RealArray.t) =
  let u = y.{0} in   (* access current solution *)
  let v = y.{1} in
  let w = y.{2} in

  (* Fill in ODE RHS function *)
  ydot.{0} <- -0.04*.u +. 1.e4*.v*.w;
  ydot.{1} <- 0.04*.u -. 1.e4*.v*.w -. 3.e7*.v*.v;
  ydot.{2} <- 3.e7*.v*.v

(* Jacobian routine to compute J(t,y) = df/dy. *)
let jac { ARKStep.jac_y = (y : RealArray.t) } j =
  let v = y.{1} in   (* access current solution *)
  let w = y.{2} in

  (* Fill in the Jacobian of the ODE RHS function *)
  let open Matrix.Dense in
  set j 0 0 (-0.04);
  set j 0 1 (1.e4*.w);
  set j 0 2 (1.e4*.v);

  set j 1 0 (0.04);
  set j 1 1 (-1.e4*.w -. 6.e7*.v);
  set j 1 2 (-1.e4*.v);

  set j 2 1 (6.e7*.v)

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in                             (* initial time *)
  let tf = 1.e11 in                           (* final time *)
  let dTout = (tf-.t0)/.100.0 in              (* time between outputs *)
  let nt = int_of_float (ceil(tf/.dTout)) in  (* number of output times *)

  (* set up the initial conditions, tolerances, initial time step size *)
  let u0 = 1.0 in
  let v0 = 0.0 in
  let w0 = 0.0 in
  let reltol = 1.e-4 in
  let abstol = 1.e-11 in
  let h0 = 1.e-4 *. reltol in

  (* Initial problem output *)
  printf "\nRobertson ODE test problem:\n";
  printf "    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n" u0 v0 w0;

  (* Create serial vector for solution *)
  (* Set initial conditions into y *)
  let y = RealArray.of_array [| u0; v0; w0 |] in
  let y_nv = Nvector_serial.wrap y in

  (* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time t0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  let m = Matrix.dense 3 in
  let arkode_mem = ARKStep.(
    init
      (implicit ~lsolver:Dls.(solver ~jac:jac (dense y_nv m)) f)
      (SStolerances (reltol, abstol))
      t0
      y_nv
  ) in
  ARKStep.set_init_step arkode_mem h0;         (* Set custom initial step *)
  ARKStep.set_max_err_test_fails arkode_mem 20;(* Increase max error test fails*)
  ARKStep.set_max_nonlin_iters arkode_mem 8;   (* Increase max nonlin iters  *)
  ARKStep.set_nonlin_conv_coef arkode_mem 1.e-7;(* Nonlinear convergence coeff.*)
  ARKStep.set_max_num_steps arkode_mem 100000; (* Increase max num steps *)

  if sundials_270_or_later then
    ARKStep.(set_predictor_method arkode_mem MaximumOrderPredictor);

  (* Open output stream for results, output comment line *)
  let ufid = open_out "solution.txt" in
  fprintf ufid "# t u v w\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16e %.16e %.16e %.16e\n" t0 y.{0} y.{1} y.{2};

  (* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let tout = ref (t0 +. dTout) in
  printf "        t           u           v           w\n";
  printf "   --------------------------------------------------\n";
  printf "  %10.3e  %12.5e  %12.5e  %12.5e\n" t0 y.{0} y.{1} y.{2};
  (try
     for _ = 0 to nt-1 do
       (* call integrator *)
       let t, _ = ARKStep.evolve_normal arkode_mem !tout y_nv in

       (* access/print solution *)
       printf "  %10.3e  %12.5e  %12.5e  %12.5e\n" t y.{0} y.{1} y.{2};
       fprintf ufid " %.16e %.16e %.16e %.16e\n" t y.{0} y.{1} y.{2};
       (* successful solve: update time *)
       tout := min (!tout +. dTout) tf
     done
   with _ -> (* unsuccessful solve: break *)
             fprintf stderr "Solver failure, stopping integration\n");
  printf("   --------------------------------------------------\n");
  close_out ufid;

  (* Print some final statistics *)
  if Sundials_impl.Version.lt620 then begin
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
    printf "   Total number of nonlinear solver convergence failures = %d\n" ncfn;
    printf "   Total number of error test failures = %d\n" netf
  end else begin
    printf "\nFinal Statistics:\n";
    ARKStep.print_all_stats arkode_mem Logfile.stdout Sundials.OutputTable;
    let fid = Logfile.openfile "ark_robertson_stats.csv" in
    ARKStep.print_all_stats arkode_mem fid Sundials.OutputCSV;
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
