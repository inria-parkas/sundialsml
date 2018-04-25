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
 * The following test simulates a simple 1D heat equation,
 *    u_t = k*u_xx + f
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  0
 * Dirichlet boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0,
 * and a point-source heating term,
 *    f = 1 for x=0.5.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * This program solves the problem with either an ERK or DIRK
 * method.  For the DIRK method, we use a Newton iteration with
 * the PCG linear solver, and a user-supplied Jacobian-vector
 * product routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *---------------------------------------------------------------*)

module RealArray = Sundials.RealArray
let printf = Printf.printf
let fprintf = Printf.fprintf
let n_vdotprod = Nvector_serial.Ops.n_vdotprod

let sundials_270_or_later =
  match Sundials.sundials_version with
  | 2,5,_ | 2,6,_ -> false
  | _ -> true

(* user data structure *)
type user_data = {
    n   : int;    (* number of intervals   *)
    dx  : float;  (* mesh spacing          *)
    k   : float;  (* diffusion coefficient *)
  }

(* Functions called by the solver *)

(* f routine to compute the ODE RHS function f(t,y). *)
let f { n; dx; k } t (y : RealArray.t) (ydot : RealArray.t) =
  RealArray.fill ydot 0.0;    (* Initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  let c1 = k/.dx/.dx in
  let c2 = -2.0*.k/.dx/.dx in
  let isource = n/2 in
  ydot.{0} <- 0.0;               (* left boundary condition *)
  for i=1 to n-1-1 do
    ydot.{i} <- c1*.y.{i-1} +. c2*.y.{i} +. c1*.y.{i+1}
  done;
  ydot.{n-1} <- 0.0;             (* right boundary condition *)
  ydot.{isource} <- ydot.{isource} +. 0.01/.dx (* source term *)

(* Jacobian routine to compute J(t,y) = df/dy. *)
let jac { n; dx; k } _ (v : RealArray.t) (jv : RealArray.t) =
  (* iterate over domain, computing all Jacobian-vector products *)
  let c1 = k/.dx/.dx in
  let c2 = -2.0*.k/.dx/.dx in
  jv.{0} <- 0.0;
  for i=1 to n-1-1 do
    jv.{i} <- c1*.v.{i-1} +. c2*.v.{i} +. c1*.v.{i+1}
  done;
  jv.{n-1} <- 0.0

(* Main Program *)
let main () =

  (* general problem parameters *)
  let t0 = 0.0 in       (* initial time *)
  let tf = 1.0 in       (* final time *)
  let nt = 10 in        (* total number of output times *)
  let rtol = 1.e-6 in   (* relative tolerance *)
  let atol = 1.e-10 in  (* absolute tolerance *)
  let mesh_n = 201 in   (* spatial mesh size *)
  let heat_k = 0.5 in   (* heat conductivity *)

  (* general problem variables *)

  (* allocate and fill udata structure *)
  let udata = {
    n  = mesh_n;
    dx = 1.0/.(1.0*.float(mesh_n)-.1.0);     (* mesh spacing *)
    k  = heat_k;
  } in

  (* Initial problem output *)
  printf "\n1D Heat PDE test problem:\n";
  printf "  N = %d\n" udata.n;
  printf "  diffusion coefficient:  k = %g\n" udata.k;

  (* Initialize data structures *)
  let data = RealArray.make mesh_n 0.0 in(* Set initial conditions *)
  let y = Nvector_serial.wrap data in (* Create serial vector for solution *)

  (* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time t0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  let arkode_mem = Arkode.(
    init
      (Implicit
        (f udata,
         Newton Spils.(solver (pcg ~maxl:mesh_n y)
                              ~jac_times_vec:(None, jac udata)
                              prec_none),
         Linear true))
      (SStolerances (rtol, atol))
      t0
      y
  ) in
  (* Set routines *)
  Arkode.set_max_num_steps arkode_mem 10000;   (* Increase max num steps  *)

  if sundials_270_or_later then
    Arkode.(set_predictor_method arkode_mem MaximumOrderPredictor);

  (* output mesh to disk *)
  let fid = open_out "heat_mesh.txt" in
  for i=0 to mesh_n-1 do
    fprintf fid "  %.16e\n" (udata.dx*.float i)
  done;
  close_out fid;

  (* Open output stream for results, access data array *)
  let ufid = open_out "heat1D.txt" in

  (* output initial condition to disk *)
  for i=0 to mesh_n-1 do
    fprintf ufid " %.16e" data.{i}
  done;
  fprintf ufid "\n";

  (* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let dTout = (tf-.t0)/.float nt in
  let tout  = ref (t0+.dTout) in
  printf "        t      ||u||_rms\n";
  printf "   -------------------------\n";
  printf "  %10.6f  %10.6f\n" t0 (sqrt((n_vdotprod y y)/.float mesh_n));
  (try
     for iout=0 to nt-1 do

       (* call integrator *)
       let t, _ = Arkode.solve_normal arkode_mem !tout y in
       (* print solution stats *)
       printf "  %10.6f  %10.6f\n" t (sqrt((n_vdotprod y y)/.float mesh_n));
       (* successful solve: update output time *)
       tout := !tout +. dTout;
       if !tout > tf then tout := tf;

       (* output results to disk *)
       for i=0 to mesh_n-1 do
         fprintf ufid " %.16e" data.{i}
       done;
       fprintf ufid "\n"
     done
   with _ ->
     (* unsuccessful solve: break *)
     fprintf stderr "Solver failure, stopping integration\n");
  printf "   -------------------------\n";
  close_out ufid;

  (* Print some final statistics *)
  let open Arkode in
  let nst      = get_num_steps arkode_mem in
  let nst_a    = get_num_step_attempts arkode_mem in
  let nfe, nfi = get_num_rhs_evals arkode_mem in
  let nsetups  = get_num_lin_solv_setups arkode_mem in
  let netf     = get_num_err_test_fails arkode_mem in
  let nni      = get_num_nonlin_solv_iters arkode_mem in
  let ncfn     = get_num_nonlin_solv_conv_fails arkode_mem in
  let nli      = Spils.get_num_lin_iters arkode_mem in
  let nJv      = Spils.get_num_jtimes_evals arkode_mem in
  let nlcf     = Spils.get_num_conv_fails arkode_mem in

  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe nfi;
  printf "   Total linear solver setups = %d\n" nsetups;
  printf "   Total linear iterations = %d\n" nli;
  printf "   Total number of Jacobian-vector products = %d\n" nJv;
  printf "   Total number of linear solver convergence failures = %d\n" nlcf;
  printf "   Total number of Newton iterations = %d\n" nni;
  printf "   Total number of nonlinear solver convergence failures = %d\n" ncfn;
  printf "   Total number of error test failures = %d\n" netf

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
