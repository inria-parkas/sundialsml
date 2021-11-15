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
 * This is the testing routine to check the SUNNonlinearSolver fixed point
 * module. This test solves the nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * where the fixed point function is
 *
 * g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * This system has the analytic solution x = 1/2, y = 1, z = -pi/6.
 * ---------------------------------------------------------------------------*)

open Sundials
module NLS = NonlinearSolver

let sungte500 =
  let n, _, _ = Sundials_configuration.sundials_version in
  n >= 5

let printf = Printf.printf

let neq = 3              (* number of equations        *)
let tol = if sungte500 then 100.0 *. sqrt(Config.unit_roundoff) else 1.0e-4
let mxiter = if sungte500 then 20 else 10

let zero = 0.0
let ptone = 0.1
let half = 0.5
let ptnine = 0.9
let one = 1.0
let oneptzerosix = 1.06
let three = 3.0
let six = 6.0
let nine = 9.0
let ten = 10.0
let twenty = 20.0
let sixty = 60.0
let pi = 3.1415926535898

(* analytic (>=5.0.0) or approximate (<5.0.0) solution *)
let xtrue = half
let ytrue = if sungte500 then one else zero
let ztrue = -.pi/.six

(* Proxy for integrator convergence test function *)
let conv_test (y : RealArray.t) (del : RealArray.t) tol (ewt : RealArray.t) _ =
  (* compute the norm of the correction *)
  if Nvector_serial.DataOps.maxnorm del <= tol
  then NLS.Success else NLS.Continue

(* -----------------------------------------------------------------------------
 * Nonlinear system
 *
 * 3x - cos(yz) - 1/2 = 0
 * x^2 - 81(y+0.1)^2 + sin(z) + 1.06 = 0
 * exp(-xy) + 20z + (10 pi - 3)/3 = 0
 *
 * Nonlinear fixed point function
 *
 * g1(x,y,z) = 1/3 cos(yz) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) - 0.1
 * g3(x,y,z) = -1/20 exp(-xy) - (10 pi - 3) / 60
 *
 * ---------------------------------------------------------------------------*)
let fp_function_lt500 y f _ =
  f.{0} <- (one/.three) *. cos(y.{1}*.y.{2}) +. (one/.six);
  f.{1} <- (one/.nine) *. sqrt(y.{0}*.y.{0} +. sin(y.{2}) +. oneptzerosix)
              -. ptone;
  f.{2} <- -.(one/.twenty) *. exp(-.y.{0}*.y.{1})
                -. (ten *. pi -. three) /. sixty

(* -----------------------------------------------------------------------------
 * Nonlinear system F(x,y,z):
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * Nonlinear fixed point function G(x,y,z):
 *
 * G1(x,y,z) = 1/3 cos((y-1)yz) + 1/6
 * G2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * G3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * Corrector form g(x,y,z):
 *
 * g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6 - x0
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9 - y0
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60 - z0
 *
 * ---------------------------------------------------------------------------*)
let fp_function (y0, ycur) ycor gvec _ =
  (* update state based on current correction *)
  Nvector_serial.DataOps.linearsum one y0 one ycor ycur;
  let x = ycur.{0} in
  let y = ycur.{1} in
  let z = ycur.{2} in
  gvec.{0} <- (one/.three) *. cos((y -. one) *. z) +. (one/.six);
  gvec.{1} <- (one/.nine) *. sqrt(x *. x +. sin(z) +. oneptzerosix) +. ptnine;
  gvec.{2} <- -.(one/.twenty) *. exp(-. x *. (y -. one))
              -. (ten *. pi -. three) /. sixty;
  Nvector_serial.DataOps.linearsum one gvec (-.one) y0 gvec

let check_ans data tol =
  (* print the solution *)
  printf "Computed solution:\n";
  printf "    y1 = %g\n" data.{0};
  printf "    y2 = %g\n" data.{1};
  printf "    y3 = %g\n" data.{2};

  (* solution error *)
  let ex = Float.abs(data.{0} -. xtrue) in
  let ey = Float.abs(data.{1} -. ytrue) in
  let ez = Float.abs(data.{2} -. ztrue) in

  (* print the solution error *)
  printf "Solution error:\n";
  printf "    ex = %g\n" ex;
  printf "    ey = %g\n" ey;
  printf "    ez = %g\n" ez;

  let tol = tol *. ten in
  if ex > tol || ey > tol || ez > tol then printf "FAIL\n"
  else printf "PASS\n"

(* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*)
let main () =
  let argc = Array.length Sys.argv in
  let maa = if argc > 1 then int_of_string Sys.argv.(1) else 0 in
  let damping = if argc > 2 then float_of_string Sys.argv.(2) else 1.0 in

  if sungte500 then begin
    (* Print problem description *)
    printf "Solve the nonlinear system:\n";
    printf "    3x - cos((y-1)z) - 1/2 = 0\n";
    printf "    x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0\n";
    printf "    exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0\n";
    printf "Analytic solution:\n";
    printf "    x = %g\n" xtrue;
    printf "    y = %g\n" ytrue;
    printf "    z = %g\n" ztrue;
    printf "Solution method: Anderson accelerated fixed point iteration.\n";
    printf "    tolerance = %g\n" tol;
    printf "    max iters = %d\n" mxiter;
    printf "    accel vec = %d\n" maa;
    printf "    damping   = %g\n" damping
  end;

  (* set initial guess *)
  let y0data = RealArray.of_array [| ptone; ptone; -.ptone |] in
  let y0 = Nvector_serial.wrap y0data in
  let ycor = Nvector_serial.make neq 0.0 in
  let ycordata = Nvector_serial.unwrap ycor in
  let ycur = Nvector.clone y0 in
  let ycurdata = Nvector.unwrap ycur in
  (* set weights *)
  let w = Nvector_serial.make neq one in

  (* create nonlinear solver *)
  let nls = NLS.FixedPoint.make ~acceleration_vectors:maa y0 in

  (* set the nonlinear residual function *)
  NLS.set_sys_fn nls
    (if sungte500 then (fp_function (y0data, ycordata))
     else fp_function_lt500);

  (* set the convergence test function *)
  NLS.(set_convtest_fn nls (OConvTest conv_test));

  (* set the maximum number of nonlinear iterations *)
  NLS.set_max_iters nls mxiter;

  if sungte500 then NLS.FixedPoint.set_damping nls damping;

  (* solve the nonlinear system *)
  NLS.solve nls ~y0 ~ycor ~w tol true ();

  if sungte500 then
    (* update the initial guess with the final correction *)
    Nvector_serial.DataOps.linearsum one y0data one ycordata ycurdata
  else begin
    (* print the solution *)
    printf "Solution:\n";
    printf "y1 = %g\n" ycordata.{0};
    printf "y2 = %g\n" ycordata.{1};
    printf "y3 = %g\n" ycordata.{2};

    (* print the solution error *)
    printf "Solution Error:\n";
    printf "e1 = %g\n" (ycordata.{0} -. xtrue);
    printf "e2 = %g\n" (ycordata.{1} -. ytrue);
    printf "e3 = %g\n" (ycordata.{2} -. ztrue)
  end;

  (* get the number of linear iterations *)
  let niters = NLS.get_num_iters nls in

  printf "Number of nonlinear iterations: %d\n" niters;

  if sungte500 then check_ans ycurdata tol
  else printf "SUCCESS\n"

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

