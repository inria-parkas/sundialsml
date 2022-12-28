(* ----------------------------------------------------------------- {{{
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Oct 2021.
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
 * The following is a simple example problem based off of
 * cvRoberts_klu.c. We simulate a scenario where a set of independent
 * ODEs are grouped together to form a larger system. For simplicity,
 * each set of ODEs is the same problem. The problem is from chemical
 * kinetics, and consists of the following three rate equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration, the KLU sparse direct linear solver, and a user-supplied
 * Jacobian routine. It uses a scalar relative tolerance and a vector
 * absolute tolerance. Output is printed in decades from t = .4 to t =
 * 4.e10. Run statistics (optional outputs) are printed at the end.
 *
 * The program takes one optional argument, the number of groups
 * of independent ODE systems:
 *
 *    ./cvRoberts_block_klu [number of groups]
 *
 * The problem is comparable to the CUDA version -
 * cvRoberts_block_cusolversp_batchqr.cu. It was based off of the
 * cvRoberts_klu.c example.
 * ----------------------------------------------------------------- }}} *)

open Sundials

let printf = Printf.printf

(* User-defined vector and matrix accessor macro: Ith *)

(* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..neq] and neq is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0. *)

let ith (v : RealArray.t) i = v.{i - 1}
let set_ith (v : RealArray.t) i e = v.{i - 1} <- e

(* Problem Constants *)

let groupsize = 3    (* number of equations per group *)
let y1 =    1.0      (* initial y components *)
let y2 =    0.0
let y3 =    0.0
let rtol =  1.0e-4   (* scalar relative tolerance            *)
let atol1 = 1.0e-8   (* vector absolute tolerance components *)
let atol2 = 1.0e-14
let atol3 = 1.0e-6
let t0 =    0.0      (* initial time           *)
let t1 =    0.4      (* first output time      *)
let tmult = 10.0     (* output time factor     *)
let nout =  12       (* number of output times *)

(* user data structure *)
type user_data = {
  ngroups : int;
  neq : int;
}

(*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 *)
let f { neq; _ } _ (y : RealArray.t) (ydot : RealArray.t) =
  let rec dogroup groupj =
    if groupj = neq then ()
    else begin
      let y1 = ith y (1+groupj) in
      let y2 = ith y (2+groupj) in
      let y3 = ith y (3+groupj) in
      let yd1 = -0.04 *. y1 +. 1.0e4 *. y2 *. y3 in
      let yd3 = 3.0e7 *. y2 *. y2 in

      set_ith ydot (1+groupj) yd1;
      set_ith ydot (3+groupj) yd3;
      set_ith ydot (2+groupj) (-. yd1 -. yd3);

      dogroup (groupj + groupsize)
    end
  in dogroup 0

(*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 *)

let jac { ngroups; _ } {Cvode.jac_y = (ydata : RealArray.t)} smat =
  let nnzper = groupsize * groupsize in
  Matrix.Sparse.set_to_zero smat;
  let colvals, rowptrs, data = Matrix.Sparse.unwrap smat in
  let idx = Index.of_int in

  rowptrs.{0} <- idx 0;
  for groupj = 0 to ngroups - 1 do
    (* get y values *)
    let y2 = ydata.{groupsize * groupj + 1} in
    let y3 = ydata.{groupsize * groupj + 2} in

    (* there are 3 entries per row *)
    rowptrs.{1 + groupsize * groupj}     <- idx (3 + nnzper * groupj);
    rowptrs.{1 + groupsize * groupj + 1} <- idx (6 + nnzper * groupj);
    rowptrs.{1 + groupsize * groupj + 2} <- idx (9 + nnzper * groupj);

    (* first row of block *)
    data.{nnzper * groupj}     <- -0.04;
    data.{nnzper * groupj + 1} <- 1.0e4 *. y3;
    data.{nnzper * groupj + 2} <- 1.0e4 *. y2;
    colvals.{nnzper * groupj}     <- idx (groupsize * groupj);
    colvals.{nnzper * groupj + 1} <- idx (groupsize * groupj + 1);
    colvals.{nnzper * groupj + 2} <- idx (groupsize * groupj + 2);

    (* second row of block *)
    data.{nnzper * groupj + 3} <- 0.04;
    data.{nnzper * groupj + 4} <- (-1.0e4 *. y3) -. (6.0e7 *. y2);
    data.{nnzper * groupj + 5} <- -1.0e4 *. y2;
    colvals.{nnzper * groupj + 3} <- idx (groupsize * groupj);
    colvals.{nnzper * groupj + 4} <- idx (groupsize * groupj + 1);
    colvals.{nnzper * groupj + 5} <- idx (groupsize * groupj + 2);

    (* third row of block *)
    data.{nnzper * groupj + 6} <- 0.0;
    data.{nnzper * groupj + 7} <- 6.0e7 *. y2;
    data.{nnzper * groupj + 8} <- 0.0;
    colvals.{nnzper * groupj + 6} <- idx (groupsize * groupj);
    colvals.{nnzper * groupj + 7} <- idx (groupsize * groupj + 1);
    colvals.{nnzper * groupj + 8} <- idx (groupsize * groupj + 2)
  done

(*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 *)

let print_output t y1 y2 y3 =
  printf "At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n" t y1 y2 y3

(*
 * Get and print some final statistics
 *)

let print_final_stats cvode_mem =
  let open Cvode in
  let nst      = get_num_steps cvode_mem in
  let nfe      = get_num_rhs_evals cvode_mem in
  let nsetups  = get_num_lin_solv_setups cvode_mem in
  let netf     = get_num_err_test_fails cvode_mem in
  let nni      = get_num_nonlin_solv_iters cvode_mem in
  let nnf      = get_num_nonlin_solv_conv_fails cvode_mem in
  let nje      = Dls.get_num_jac_evals cvode_mem in
  let nge      = get_num_g_evals cvode_mem in
  printf "\nFinal Statistics:\n";
  printf "nst = %-6d nfe  = %-6d nsetups = %-6d nje = %d\n" nst nfe nsetups nje;
  if Sundials_impl.Version.lt620
  then printf "nni = %-6d ncfn = %-6d netf = %-6d    nge = %d\n \n"
       nni nnf netf nge
  else let ncfn   = get_num_step_solve_fails cvode_mem in
       printf "nni = %-6d nnf = %-6d netf = %-6d    ncfn = %-6d nge = %d\n\n"
              nni nnf netf ncfn nge

(*
 *-------------------------------
 * Main Program
 *-------------------------------
 *)

let main () =
  (* Parse command line arguments *)
  let ngroups =
    if Array.length Sys.argv > 1 then int_of_string Sys.argv.(1)
    else 1000
  in
  let neq = ngroups * groupsize in
  let udata = { neq; ngroups } in

  (* Create serial vector of length neq for I.C. and abstol *)
  let y = Nvector_serial.make neq 0.0 in
  let ydata = Nvector.unwrap y in
  let abstol = Nvector_serial.make neq 0.0 in
  let abstoldata = Nvector.unwrap abstol in

  (* Initialize y *)
  let rec dogroup groupj =
    if groupj < neq then begin
      set_ith ydata (1+groupj) y1;
      set_ith ydata (2+groupj) y2;
      set_ith ydata (3+groupj) y3;
      dogroup (groupj + groupsize)
    end else ()
  in
  dogroup 0;

  (* Set the scalar relative tolerance *)
  let reltol = rtol in

  (* Set the vector absolute tolerance *)
  let rec dogroup groupj =
    if groupj < neq then begin
      set_ith abstoldata (1+groupj) atol1;
      set_ith abstoldata (2+groupj) atol2;
      set_ith abstoldata (3+groupj) atol3;
      dogroup (groupj + groupsize)
    end else ()
  in
  dogroup 0;

  (* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula *)
  (* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. *)
  (* Call CVodeSetUserData to attach the user data structure *)
  (* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances *)
  (* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode *)
  (* Create sparse SUNMatrix for use in linear solves *)
  (* Create KLU solver object for use by CVode *)
  (* Set the user-supplied Jacobian routine Jac *)
  let nnz = groupsize * groupsize * ngroups in
  let a = Matrix.sparse_csr ~nnz neq in
  let cvode_mem =
    Cvode.(init BDF
             ~lsolver:Dls.(solver ~jac:(jac udata) (klu y a))
             (SVtolerances (reltol, abstol)) (f udata) t0 y)
  in

  (* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  *)
  printf " \nGroup of independent 3-species kinetics problems\n\n";
  printf "number of groups = %d\n\n" ngroups;

  let iout = ref 0 in
  let tout = ref t1 in
  (try
    while (!iout <> nout) do
      let t, retval = Cvode.solve_normal cvode_mem !tout y in

      printf "group %d: " 0;
      print_output t (ith ydata 1) (ith ydata 2) (ith ydata 3);

      match retval with
      | Cvode.Success ->
          incr iout;
          tout := !tout *. tmult
      | _ -> raise Exit
    done
  with Exit -> ());

  (* Print some final statistics *)
  print_final_stats cvode_mem

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

