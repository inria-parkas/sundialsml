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
 * The following is a simple example problem with analytical
 * solution adapted from example 10.2 of Ascher & Petzold, "Computer
 * Methods for Ordinary Differential Equations and
 * Differential-Algebraic Equations," SIAM, 1998, page 267:
 *    x1'(t) = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t)
 *         0 = (t+2)*x1 - (t+2)*exp(t)
 * for t in the interval [0.0, 1.0], with initial condition:
 *    x1(0) = 1   and   x2(0) = -1/2.
 * The problem has true solution
 *    x1(t) = exp(t)  and  x2(t) = exp(t)/(t-2)
 *
 * This program solves the problem with IDA using a custom
 * 'matrix-embedded' SUNLinearSolver. Output is printed
 * every 0.1 units of time (10 total).  Run statistics (optional
 * outputs) are printed at the end.
 * -----------------------------------------------------------------*)

open Sundials

let printf = Printf.printf

(*-------------------------------
 * Functions called by the solver
 *-------------------------------*)

(* System residual function:
      0 = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t) - x1'(t)
      0 = (t+2)*x1 - (t+2)*exp(t)
*)

let fres alpha t (yy : RealArray.t) (yp : RealArray.t) (rr : RealArray.t) =
  let x1 = yy.{0} in               (* access current solution values *)
  let x2 = yy.{1} in
  let x1p = yp.{0} in              (* access current derivative values *)
  rr.{0} <- (1.0 -. alpha) /. (t -. 2.0) *. x1
            -. x1
            +. (alpha -. 1.0) *. x2
            +. 2.0 *. exp t
            -. x1p;
  rr.{1} <- (t +. 2.0) *. x1 -. (t +. 2.0) *. exp t

(*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*)

type matrix_embedded_ls_content = {
  alpha : float;
  mutable ida_mem : (RealArray.t, Nvector_serial.kind) Ida.session option;
}

(* linear solve routine *)
let matrix_embedded_ls_solve content () (x : RealArray.t) (b : RealArray.t) tol =
  let { alpha; ida_mem } = content in
  match ida_mem with
  | None -> failwith "linear solver not properly configure"
  | Some ida_mem ->
      (* retrieve implicit system data from IDA *)
      let { Ida.cj; Ida.tn = tcur; _ } = Ida.get_nonlin_system_data ida_mem in
      (* perform linear solve: A*x=b
             A = df/dy + cj*df/dyp
          =>
             A = [ - cj - (alpha - 1)/(t - 2) - 1, alpha - 1]
                 [                          t + 2,         0]

       *)
      let a11 = -. cj -. (alpha -. 1.0) /. (tcur -. 2.0) -. 1.0 in
      let a12 = alpha -. 1.0 in
      let a21 = tcur +. 2.0 in
      let b1 = b.{0} in
      let b2 = b.{1} in
      x.{0} <- b2 /. a21;
      x.{1} <- -. (a11 *. b2 -. a21 *. b1) /. (a12 *. a21)

(* constructor *)
let matrix_embedded_ls alpha =
  let content = { alpha; ida_mem = None } in
  (fun session -> content.ida_mem <- Some session),
  LinearSolver.Custom.(
    make_without_matrix (make_ops ~solver_type:MatrixEmbedded
                                  ~solve:matrix_embedded_ls_solve ())
                        content)

(*-------------------------------
 * Private helper functions
 *-------------------------------*)

(* routine to fill analytical solution and its derivative *)
let analytical_solution t (y : Nvector_serial.t) (yp : Nvector_serial.t) =
  let y = Nvector.unwrap y in
  let yp = Nvector.unwrap yp in
  y.{0} <- exp t;
  y.{1} <- exp t /. (t -. 2.0);
  yp.{0} <- exp t;
  yp.{1} <- exp t /. (t -. 2.0) -. exp t /. (t -. 2.0) /. (t -. 2.0)

(* check the computed solution *)
let check_ans (y : Nvector_serial.t) t rtol atol =
  (* create solution and error weight vectors *)
  let ytrue = Nvector.clone y in
  let ewt = Nvector.clone y in
  let abstol = Nvector.clone y in

  (* set the solution data *)
  analytical_solution t ytrue abstol;

  (* compute the error weight vector, loosen atol *)
  Nvector_serial.Ops.const atol abstol;
  Nvector_serial.Ops.abs ytrue ewt;
  Nvector_serial.Ops.linearsum rtol ewt 10.0 abstol ewt;

  if Nvector_serial.Ops.min ewt <= 0.0 then
    (Format.eprintf "\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n\n"; false)
  else begin
    Nvector_serial.Ops.inv ewt ewt;

    (* compute the solution error *)
    Nvector_serial.Ops.linearsum 1.0 y (-1.0) ytrue ytrue;
    let err = Nvector_serial.Ops.wrmsnorm ytrue ewt in

    (* is the solution within the tolerances? *)
    if err >= 1.0 then
      (Format.eprintf "\nSUNDIALS_WARNING: check_ans error=%g\n\n" err; false)
    else true
  end

(*---- end of file ----*)

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in         (* initial time *)
  let tf = 1.0 in         (* final time *)
  let dTout = 0.1 in      (* time between outputs *)
  let neq = 2 in          (* number of dependent vars. *)
  let reltol = 1.0e-4 in  (* tolerances *)
  let abstol = 1.0e-9 in
  let alpha  = 10.0 in    (* stiffness parameter *)

  (* Initial diagnostics output *)
  printf "\nAnalytical DAE test problem:\n";
  printf "    alpha = %g\n"     alpha;
  printf "   reltol = %.1e\n"   reltol;
  printf "   abstol = %.1e\n\n" abstol;

  (* Initialize data structures *)
  let yy = Nvector_serial.make neq 0.0 in (* Create serial vector for solution *)
  let yp = Nvector.clone yy in
  analytical_solution t0 yy yp;           (* Specify initial conditions *)

  (* Create custom matrix-embedded linear solver *)
  let set_ls_session, ls = matrix_embedded_ls alpha in
  (* Call IDACreate and IDAInit to initialize IDA memory *)
  (* Attach the linear solver *)
  let ida_mem = Ida.(init (SStolerances (reltol, abstol))
                          ~lsolver:(matrix_embedded_solver ls)
                          (fres alpha) t0 yy yp)
  in
  set_ls_session ida_mem;

  (* In loop, call IDASolve, print results, and test for error.
     Stops when the final time has been reached. *)
  printf "        t          x1         x2\n";
  printf "   ----------------------------------\n";
  let yydata = Nvector.unwrap yy in
  let rec loop t tout =
    if tf -. t <= 1.0e-15 then t
    else begin
      let t, _ = Ida.solve_normal ida_mem tout yy yp in   (* call integrator *)
      printf "  %10.6f  %10.6f  %10.6f\n" t yydata.{0} yydata.{1}; (* access/print solution *)
      loop t (min tf (tout +. dTout))
    end
  in
  let t =
    try loop t0 (t0 +. dTout);
    with _ -> (Format.eprintf "Solver failure, stopping integration\n"; 0.0)
  in
  printf "   ----------------------------------\n";

  (* Get/print some final statistics on how the solve progressed *)
  let nst = Ida.get_num_steps ida_mem in
  let nre = Ida.get_num_res_evals ida_mem in
  let nni = Ida.get_num_nonlin_solv_iters ida_mem in
  let netf = Ida.get_num_err_test_fails ida_mem in
  let ncfn = Ida.get_num_nonlin_solv_conv_fails ida_mem in
  let nreLS = Ida.Dls.get_num_lin_res_evals ida_mem in

  print_string "\nFinal Solver Statistics: \n";
  print_string "\nNumber of steps                    = "; print_int nst;
  print_string "\nNumber of residual evaluations     = "; print_int (nre + nreLS);
  print_string "\nNumber of nonlinear iterations     = "; print_int nni;
  print_string "\nNumber of error test failures      = "; print_int netf;
  print_string "\nNumber of nonlinear conv. failures = "; print_int ncfn;
  print_string "\n";

  (* check the solution error *)
  ignore (check_ans yy t reltol abstol)

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

