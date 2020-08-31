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
 * module
 * ---------------------------------------------------------------------------*)

open Sundials
module NLS = NonlinearSolver

let printf = Printf.printf

let neq = 3              (* number of equations        *)
let tol = 1.0e-4         (* nonlinear solver tolerance *)
let maxit = 10           (* max nonlinear iterations   *)

let zero = 0.0
let ptone = 0.1
let half = 0.5
let one = 1.0
let oneptzerosix = 1.06
let three = 3.0
let four = 4.0
let six = 6.0
let nine = 9.0
let ten = 10.0
let twenty = 20.0
let sixty = 60.0
let eightyone = 81.0
let pi = 3.1415926535898

(* approximate solution *)
let y1 = half
let y2 = zero
let y3 = -.pi/.six

(* Proxy for integrator convergence test function *)
let conv_test (y : RealArray.t) (del : RealArray.t) tol (ewt : RealArray.t) _ =
  (* compute the norm of the correction *)
  if Nvector_serial.DataOps.n_vmaxnorm del <= tol
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
let fp_function y f _ =
  f.{0} <- (one/.three) *. cos(y.{1}*.y.{2}) +. (one/.six);
  f.{1} <- (one/.nine) *. sqrt(y.{0}*.y.{0} +. sin(y.{2}) +. oneptzerosix)
              -. ptone;
  f.{2} <- -.(one/.twenty) *. exp(-.y.{0}*.y.{1})
                -. (ten *. pi -. three) /. sixty

(* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*)
let main () =
  (* create vectors *)
  let x = Nvector_serial.make neq 0.0 in
  (* set initial guess *)
  let y0 = Nvector_serial.wrap (RealArray.of_array [| ptone; ptone; -.ptone |]) in
  let y = Nvector_serial.Ops.n_vclone x in
  let ydata = Nvector_serial.unwrap y in
  (* set weights *)
  let w = Nvector_serial.make neq one in

  (* create nonlinear solver *)
  let nls = NLS.FixedPoint.make y in

  (* set the nonlinear residual function *)
  NLS.set_sys_fn nls fp_function;

  (* set the convergence test function *)
  NLS.set_convtest_fn nls conv_test;

  (* set the maximum number of nonlinear iterations *)
  NLS.set_max_iters nls maxit;

  (* solve the nonlinear system *)
  NLS.solve nls ~y0 ~y ~w tol true;

  (* print the solution *)
  printf "Solution:\n";
  printf "y1 = %g\n" ydata.{0};
  printf "y2 = %g\n" ydata.{1};
  printf "y3 = %g\n" ydata.{2};

  (* print the solution error *)
  printf "Solution Error:\n";
  printf "e1 = %g\n" (ydata.{0} -. y1);
  printf "e2 = %g\n" (ydata.{1} -. y2);
  printf "e3 = %g\n" (ydata.{2} -. y3);

  (* get the number of linear iterations *)
  let niters = NLS.get_num_iters nls in

  printf "Number of nonlinear iterations: %d\n" niters;

  (* Print result *)
  printf "SUCCESS\n"

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

