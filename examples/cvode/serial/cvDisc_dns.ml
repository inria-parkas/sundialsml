(*
 * -----------------------------------------------------------------
 * Programmers: Radu Serban, Alan Hindmarsh, and Cody Balos @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Feb 2019.
 * -----------------------------------------------------------------
 * Simple 1D example to illustrate integrating over discontinuities:
 *
 * A) Discontinuity in solution
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       y' = -y   ; y(1) = 1    ; t = [1,2]
 *
 * B) Discontinuity in RHS (y')
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       z' = -5*z ; z(1) = y(1) ; t = [1,2]
 *    This case is solved twice, first by explicitly treating the
 *    discontinuity point and secondly by letting the integrator
 *    deal with the discontinuity.
 * -----------------------------------------------------------------
 *)

open Sundials

let printf = Printf.printf

(* XXX
#include <cvode/cvode.h>               /* prototypes for CVODE functions and const */
#include <cvode/cvode_direct.h>        /* access to CVDls interface                */
#include <nvector/nvector_serial.h>    /* access to serial NVector                 */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver          */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype          */
*)

(* Problem Constants *)
let neq = 1 (* number of equations *)
type flag = RHS1 | RHS2

(*
 * RHS function
 * The form of the RHS function is controlled by the flag passed as f_data:
 *   flag = RHS1 -> y' = -y
 *   flag = RHS2 -> y' = -5*y
 *)
let f flag t (y : RealArray.t) (ydot : RealArray.t) =
  match !flag with
  | RHS1 -> ydot.{0} <- -. y.{0}
  | RHS2 -> ydot.{0} <- -.5.0 *. y.{0}

let solve cvode_mem ynv y t_stop =
  let rec go t =
    if t < t_stop then begin
      (* advance solver just one internal step *)
      let (t', _) = Cvode.solve_one_step cvode_mem t_stop ynv in
      printf "%12.8e  %12.8e\n" t' y.{0};
      go t'
    end
  in
  go

let main () =
(* XXX
  void *cvode_mem;
  SUNMatrix A;
  SUNLinearSolver LS;
  N_Vector y;
  int flag, ret;
  realtype reltol, abstol, t0, t1, t2, t;
  long int nst1, nst2, nst;
*)
  let reltol = 1.0e-3 in
  let abstol = 1.0e-4 in
  let t0 = 0.0 in
  let t1 = 1.0 in
  let t2 = 2.0 in

  (* Allocate the vector of initial conditions *)
  (* Set initial condition *)
  let ynv = Nvector_serial.make neq 1.0 in
  let y = Nvector_serial.unwrap ynv in

  (*
   * ------------------------------------------------------------
   *  Shared initialization and setup
   * ------------------------------------------------------------
   *)

  (* Provide RHS flag as user data which can be access in user provided routines *)
  let flag = ref RHS1 in
  (* Create dense SUNMatrix for use in linear solver *)
  let a = Matrix.dense neq in
  (* Call CVodeCreate to create CVODE memory block and specify the
   * Backward Differentiaion Formula and the use of a Newton Iteration *)
  (* Call CVodeInit to initialize integrator memory and specify the
   * user's right hand side function y'=f(t,y), the initial time T0
   * and the initial condiition vector y. *)
  (* Call CVodeSStolerances to specify integration tolereances,
   * specifically the scalar relative and absolute tolerance. *)
  (* Create dense linear solver for use by CVode *)
  (* Attach the linear solver and matrix to CVode by calling CVDlsSetLinearSolver *)
  let cvode_mem = Cvode.(init
    BDF (Newton Dls.(solver (dense ynv a)))
    (SStolerances (reltol, abstol))
    (f flag) t0 ynv)
  in

  (*
   * ---------------------------------------------------------------
   * Discontinuity in the solution
   *
   * 1) Integrate to the discontinuity
   * 2) Integrate from the discontinuity
   * ---------------------------------------------------------------
   *)

  (* ---- Integrate to the discontinuity *)
  printf "\nDiscontinuity in solution\n\n";

  (* set TSTOP (max time solution proceeds to) - this is not required *)
  Cvode.set_stop_time cvode_mem t1;

  flag := RHS1; (* use -y for RHS *)

  printf "%12.8e  %12.8e\n" t0 y.{0};
  solve cvode_mem ynv y t1 t0;
  (* Get the number of steps the solver took to get to the discont. *)
  let nst1 = Cvode.get_num_steps cvode_mem in

  (* ---- Integrate from the discontinuity *)

  (* Include discontinuity *)
  y.{0} <- 1.0;

  (* Reinitialize the solver *)
  Cvode.reinit cvode_mem t1 ynv;

  (* set TSTOP (max time solution proceeds to) - this is not required *)
  Cvode.set_stop_time cvode_mem t2;

  flag := RHS1; (* use -y for RHS *)

  printf "%12.8e  %12.8e\n" t1 y.{0};
  solve cvode_mem ynv y t2 t1;

  (* Get the number of steps the solver took after the discont. *)
  let nst2 = Cvode.get_num_steps cvode_mem in

  (* Print statistics *)
  let nst = nst1 + nst2 in
  printf "\nNumber of steps: %d + %d = %d\n" nst1 nst2 nst;

  (*
   * ---------------------------------------------------------------
   * Discontinuity in RHS: Case 1 - explicit treatment
   * Note that it is not required to set TSTOP, but without it
   * we would have to find y(t1) to reinitialize the solver.
   * ---------------------------------------------------------------
   *)

  printf "\nDiscontinuity in RHS: Case 1 - explicit treatment\n\n";

  (* Set initial condition *)
  y.{0} <- 1.0;

  (* Reinitialize the solver. CVodeReInit does not reallocate memory
   * so it can only be used when the new problem size is the same as
   * the problem size when CVodeCreate was called. *)
  Cvode.reinit cvode_mem t0 ynv;

  (* ---- Integrate to the discontinuity *)

  (* Set TSTOP (max time solution proceeds to) to location of discont. *)
  Cvode.set_stop_time cvode_mem t1;

  flag := RHS1; (* use -y for RHS *)

  printf "%12.8e  %12.8e\n" t0 y.{0};
  solve cvode_mem ynv y t1 t0;

  (* Get the number of steps the solver took to get to the discont. *)
  let nst1 = Cvode.get_num_steps cvode_mem in

  (* If TSTOP was not set, we'd need to find y(t1): *)
  (* CVodeGetDky(cvode_mem, t1, 0, y); *)

  (* ---- Integrate from the discontinuity *)

  (* Reinitialize solver *)
  Cvode.reinit cvode_mem t1 ynv;

  (* set TSTOP (max time solution proceeds to) - this is not required *)
  Cvode.set_stop_time cvode_mem t2;

  flag := RHS2; (* use -5y for RHS *)

  printf "%12.8e  %12.8e\n" t1 y.{0};
  solve cvode_mem ynv y t2 t1;

  (* Get the number of steps the solver took after the discont. *)
  let nst2 = Cvode.get_num_steps cvode_mem in

  (* Print statistics *)
  let nst = nst1 + nst2 in
  printf "\nNumber of steps: %d + %d = %d\n" nst1 nst2 nst;

  (*
   * ---------------------------------------------------------------
   * Discontinuity in RHS: Case 2 - let CVODE deal with it
   * Note that here we MUST set TSTOP to ensure that the
   * change in the RHS happens at the appropriate time
   * ---------------------------------------------------------------
   *)
  printf "\nDiscontinuity in RHS: Case 2 - let CVODE deal with it\n\n";

  (* Set initial condition *)
  y.{0} <- 1.0;

  (* Reinitialize the solver. CVodeReInit does not reallocate memory
   * so it can only be used when the new problem size is the same as
   * the problem size when CVodeCreate was called. *)
  Cvode.reinit cvode_mem t0 ynv;

  (* ---- Integrate to the discontinuity *)

  (* Set TSTOP (max time solution proceeds to) to location of discont. *)
  Cvode.set_stop_time cvode_mem t1;

  flag := RHS1; (* use -y for RHS *)

  printf "%12.8e  %12.8e\n" t0 y.{0};
  solve cvode_mem ynv y t1 t0;

  (* Get the number of steps the solver took to get to the discont. *)
  let nst1 = Cvode.get_num_steps cvode_mem in

  (* ---- Integrate from the discontinuity *)

  (* set TSTOP (max time solution proceeds to) - this is not required *)
  Cvode.set_stop_time cvode_mem t2;

  flag := RHS2; (* use -5y for RHS *)

  printf "%12.8e  %12.8e\n" t1 y.{0};
  solve cvode_mem ynv y t2 t1;

  (* Get the number of steps the solver took after the discont. *)
  let nst = Cvode.get_num_steps cvode_mem in

  (* Print statistics *)
  let nst2 = nst - nst1 in
  printf "\nNumber of steps: %d + %d = %d\n" nst1 nst2 nst

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

