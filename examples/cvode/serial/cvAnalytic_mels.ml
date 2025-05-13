(* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Nov 2021.
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with analytical solution,
 * solution,
 *    dy/dt = lamda*y + 1/(1+t^2) - lamda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lamda".  The value of lamda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the BDF method, Newton
 * iteration, and a custom 'matrix-embedded' SUNLinearSolver. Output
 * is printed every 1.0 units of time (10 total).  Run statistics
 * (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*)

open Sundials

let printf = Printf.printf

(*-------------------------------
 * Functions called by the solver
 *-------------------------------*)

(* f routine to compute the ODE RHS function f(t,y). *)
let f lambda t (y : RealArray.t) (ydot : RealArray.t) =
  let u = y.{0} in          (* access current solution value *)
  (* fill in the RHS function: "NV_Ith_S" accesses the 0th entry of ydot *)
  ydot.{0} <- lambda *. u +. 1.0 /. (1.0 +. t *. t) -. lambda *. atan t

(*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*)

type matrix_embedded_ls_content = {
  lambda : float;
  mutable cvode_mem : (RealArray.t, Nvector_serial.kind) Cvode.session option;
}

(* linear solve routine *)
let matrix_embedded_ls_solve content _m (x : RealArray.t) (b : RealArray.t) _ =
  let { lambda; cvode_mem } = content in
  match cvode_mem with
  | None -> failwith "linear solver not properly configure"
  | Some cvode_mem ->
      (* retrieve implicit system data from ARKStep *)
      let gamma = Cvode.get_current_gamma cvode_mem in
      (* perform linear solve: (1-gamma*lamda)*x = b *)
      x.{0} <- b.{0} /. (1.0 -. gamma *. lambda)

(* constructor *)
let matrix_embedded_ls lambda =
  let content = { lambda; cvode_mem = None } in
  (fun session -> content.cvode_mem <- Some session),
  LinearSolver.Custom.(
    make_without_matrix (make_ops ~solver_type:LinearSolver.MatrixEmbedded
                                  ~solve:matrix_embedded_ls_solve ())
                        content)

(*-------------------------------
 * Private helper functions
 *-------------------------------*)

(* check the computed solution *)
let check_ans (y : RealArray.t) t rtol atol =
  (* compute solution error *)
  let ans = atan(t) in
  let ewt = 1.0 /. (rtol *. abs_float ans +. atol) in
  let err = ewt *. abs_float(y.{0} -. ans) in

  (* is the solution within the tolerances? *)
  if err >= 1.0 then begin
    printf "\nSUNDIALS_WARNING: check_ans error=%g\n\n" err;
    false
  end else true

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in         (* initial time *)
  let tf = 10.0 in        (* final time *)
  let dTout = 1.0 in      (* time between outputs *)
  let neq = 1 in          (* number of dependent vars. *)
  let reltol = 1.0e-6 in  (* tolerances *)
  let abstol = 1.0e-10 in
  let lambda  = -100.0 in (* stiffness parameter *)

  (* Initial diagnostics output *)
  print_string "\nAnalytical ODE test problem:\n";
  printf "    lamda = %g\n\
         \   reltol = %.1e\n\
         \   abstol = %.1e\n\n" lambda reltol abstol;

  (* Initialize data structures *)
  let y = Nvector_serial.make neq 0.0 in (* Create serial vector for solution *)
                                         (* Specify initial condition *)

  (* Create custom matrix-embedded linear solver *)
  let set_ls_session, ls = matrix_embedded_ls lambda in
  (* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula *)
  (* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. *)
  (* Call CVodeSStolerances to specify the scalar relative and absolute tolerances *)
  (* Call CVodeSetLinearSolver to attach the linear solver to CVode *)
  let cvode_mem = Cvode.(init
                           BDF
                           (SStolerances (reltol, abstol))
                           ~lsolver:(matrix_embedded_solver ls)
                           (f lambda)
                           t0 y)
  in
  set_ls_session cvode_mem;

  (* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  *)
  print_string "        t           u\n";
  print_string "   ---------------------\n";
  let ydata = Nvector.unwrap y in
  let rec loop t tout =
    if tf -. t <= 1.0e-15 then t
    else begin
      let t, _ = Cvode.solve_normal cvode_mem tout y in (* call integrator *)
      printf "  %10.6f  %10.6f\n" t ydata.{0};      (* access/print solution *)
      loop t (min tf (tout +. dTout))
    end
  in
  let t =
    try loop t0 (t0 +. dTout);
    with _ -> (Format.eprintf "Solver failure, stopping integration\n"; 0.0)
  in
  print_string "   ---------------------\n";

  (* Get/print some final statistics on how the solve progressed *)
  let nst = Cvode.get_num_steps cvode_mem in
  let nfe = Cvode.get_num_rhs_evals cvode_mem in
  let nsetups = Cvode.get_num_lin_solv_setups cvode_mem in
  let netf = Cvode.get_num_err_test_fails cvode_mem in
  let nni = Cvode.get_num_nonlin_solv_iters cvode_mem in
  let ncfn = Cvode.get_num_nonlin_solv_conv_fails cvode_mem in
  let nje = Cvode.Dls.get_num_jac_evals cvode_mem in
  let nfeLS = Cvode.Dls.get_num_lin_rhs_evals cvode_mem in

  print_string "\nFinal Solver Statistics:\n";
  print_string "   Internal solver steps = "; print_int nst;
  print_string "\n   Total RHS evals = "; print_int nfe;
  print_string "\n   Total linear solver setups = "; print_int nsetups;
  print_string "\n   Total RHS evals for setting up the linear system = "; print_int nfeLS;
  print_string "\n   Total number of Jacobian evaluations = "; print_int nje;
  print_string "\n   Total number of Newton iterations = "; print_int nni;
  print_string "\n   Total number of linear solver convergence failures = "; print_int ncfn;
  print_string "\n   Total number of error test failures = "; print_int netf;
  print_string "\n\n";

  (* check the solution error *)
  ignore (check_ans ydata t reltol abstol)

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

