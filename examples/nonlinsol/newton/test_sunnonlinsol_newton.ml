(* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2020.
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the testing routine to check the SUNNonlinearSolver Newton module
 * ---------------------------------------------------------------------------*)

open Sundials
module LS = LinearSolver
module DLS = LinearSolver.Direct
module NLS = NonlinearSolver
module NVOps = Nvector_serial.DataOps

let printf = Printf.printf

let neq = 3        (* number of equations        *)
let tol = 1.0e-2   (* nonlinear solver tolerance *)
let maxit = 10     (* max nonlinear iterations   *)

let zero = 0.0
let half = 0.5
let one = 1.0
let two = 2.0
let three = 3.0
let four = 4.0
let six = 6.0

(* approximate solution *)
let y1 = 0.785196933062355226
let y2 = 0.496611392944656396
let y3 = 0.369922830745872357

(*
 * Proxies for integrator memory struct and functions
 *)

(* -----------------------------------------------------------------------------
 * Jacobian of the nonlinear residual function
 *
 *            ( 2x  2y  2z )
 * J(x,y,z) = ( 4x  2y  -4 )
 *            ( 6x  -4  2z )
 *
 * ---------------------------------------------------------------------------*)
let jac t y j =
  j.{0,0} <- two*.y.{0};
  j.{1,0} <- two*.y.{1};
  j.{2,0} <- two*.y.{2};

  j.{0,1} <- four*.y.{0};
  j.{1,1} <- two*.y.{1};
  j.{2,1} <- -.four;

  j.{0,2} <- six*.y.{0};
  j.{1,2} <- -.four;
  j.{2,2} <- two*.y.{2}

(* Integrator memory structure *)
type 'k integrator_mem_rec = {
  y0   : Nvector_serial.t;
  ycur : Nvector_serial.t;
  ycor : Nvector_serial.t;
  w    : Nvector_serial.t;

  x    : Nvector_serial.t;
  a    : 'k Matrix.dense;
  ls   : (Matrix.Dense.t, 'k, [`Dls]) LS.serial_t;
}

(* Proxy for integrator lsetup function *)
let lsetup imem jbad _ =
  (* compute the Jacobian *)
  jac zero (Nvector.unwrap imem.ycur) (Matrix.(Dense.unwrap (unwrap imem.a)));

  (* setup the linear solver *)
  LS.setup imem.ls imem.a;

  (* update Jacobian status *)
  true

(* Proxy for integrator lsolve function *)
let lsolve imem b _ =
  LS.solve imem.ls imem.a imem.x (Nvector_serial.wrap b) zero;
  NVOps.scale one (Nvector_serial.unwrap imem.x) b

(* Proxy for integrator convergence test function *)
let conv_test imem y del tol ewt _ =
  (* compute the norm of the correction *)
  let delnrm = NVOps.wrmsnorm del ewt in
  if delnrm <= tol then NLS.Success else NLS.Continue

(* -----------------------------------------------------------------------------
 * Nonlinear residual function
 *
 * f1(x,y,z) = x^2 + y^2 + z^2 - 1 = 0
 * f2(x,y,z) = 2x^2 + y^2 - 4z     = 0
 * f3(x,y,z) = 3x^2 - 4y + z^2     = 0
 *
 * ---------------------------------------------------------------------------*)
let res imem y f _ =
  (* update state based on current correction *)
  Nvector_serial.Ops.linearsum one imem.y0 one imem.ycor imem.ycur;
  (* compute the residual function *)
  let ycur = Nvector.unwrap imem.ycur in
  let y1, y2, y3 = ycur.{0}, ycur.{1}, ycur.{2} in
  f.{0} <- y1*.y1 +. y2*.y2 +. y3*.y3 -. one;
  f.{1} <- two *. y1*.y1 +. y2*.y2 -. four *. y3;
  f.{2} <- three *. (y1*.y1) -. four *. y2 +. y3*.y3

(* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*)
let main () =
  (* create vector *)
  let x = Nvector_serial.make neq 0.0 in
  let y0 = Nvector_serial.wrap (RealArray.of_array [| half; half; half |]) in
  let ycur = Nvector.clone y0 in
  let ydata = Nvector.unwrap ycur in
  (* set initial guess for the state *)
  let ycor = Nvector_serial.make neq zero in
  (* set weights *)
  let w = Nvector_serial.make neq one in

  (* create dense matrix *)
  let a = Matrix.dense neq in

  (* create dense linear solver *)
  (* initialize the linear solver *)
  let ls = DLS.dense y0 a in
  LS.init ls;

  (* set integrator memory *)
  let imem = { y0; ycur; ycor; w; x; a; ls } in

  (* create nonlinear solver *)
  let nls = NLS.Newton.make y0 in

  (* set the nonlinear residual function *)
  NLS.set_sys_fn nls (res imem);

  (* set the wrapper functions to linear solver setup and solve functions *)
  NLS.set_lsetup_fn nls (lsetup imem);
  NLS.set_lsolve_fn nls (lsolve imem);
  NLS.set_convtest_fn nls (conv_test imem);

  (* set the maximum number of nonlinear iterations *)
  NLS.set_max_iters nls maxit;

  (* solve the nonlinear system *)
  NLS.solve nls ~y0 ~ycor ~w tol true;

  (* update the initial guess with the final correction *)
  Nvector_serial.Ops.linearsum one y0 one ycor ycur;

  (* print the solution *)
  printf "Solution:\n";
  printf "y1 = %g\n" ydata.{0};
  printf "y2 = %g\n" ydata.{1};
  printf "y3 = %g\n" ydata.{2};

  (* print the solution error *)
  printf("Solution Error:\n");
  printf "e1 = %g\n" (ydata.{0} -. y1);
  printf "e2 = %g\n" (ydata.{1} -. y2);
  printf "e3 = %g\n" (ydata.{2} -. y3);

  (* get the number of linear iterations *)
  let niters = NLS.get_num_iters nls in
  printf "Number of nonlinear iterations: %d\n" niters;
  printf("SUCCESS\n")

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

