(* ----------------------------------------------------------------------------- {{{
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Oct 2021.
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves the nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * using the accelerated fixed pointer solver in KINSOL. The nonlinear fixed
 * point function is
 *
 * g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * This system has the analytic solution x = 1/2, y = 1, z = -pi/6.
 * --------------------------------------------------------------------------- }}} *)

open Sundials
let printf = Printf.printf

let ge580 =
  let n, m, _ = Config.sundials_version in
  n > 5 || (n = 5 && m >= 8)

let ge600 =
  let n, _, _ = Config.sundials_version in
  n >= 6

(* problem constants *)
let neq          = 3 (* number of equations *)

let zero         = 0.0             (* real 0.0  *)
let ptone        = 0.1             (* real 0.1  *)
let half         = 0.5             (* real 0.5  *)
let ptnine       = 0.9             (* real 0.9  *)
let one          = 1.0             (* real 1.0  *)
let oneptzerosix = 1.06            (* real 1.06 *)
let oneptone     = 1.1             (* real 1.1  *)
let three        = 3.0             (* real 3.0  *)
let four         = 4.0             (* real 4.0  *)
let six          = 6.0             (* real 6.0  *)
let nine         = 9.0             (* real 9.0  *)
let ten          = 10.0            (* real 10.0 *)
let twenty       = 20.0            (* real 20.0 *)
let sixty        = 60.0            (* real 60.0 *)
let eightyone    = 81.0            (* real 81.0 *)
let pi           = 3.1415926535898 (* real pi   *)

(* analytic solution *)
let xtrue = half
let ytrue = one
let ztrue = -. pi /. six

let int_of_orthaa = function
  | Kinsol.MGS -> 0
  | Kinsol.ICWY -> 1
  | Kinsol.CGS2 -> 2
  | Kinsol.DCGS2 -> 3

let orthaa_of_int = function
  | 0 -> Kinsol.MGS
  | 1 -> Kinsol.ICWY
  | 2 -> Kinsol.CGS2
  | 3 -> Kinsol.DCGS2
  | _ -> failwith "invalid orthaa type"

(* -----------------------------------------------------------------------------
 * Nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * Nonlinear fixed point function
 *
 * g1(x,y,z) = 1/3 cos((y-1)z) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * ---------------------------------------------------------------------------*)
let fpfunction (udata : RealArray.t) (gdata: RealArray.t) =
  let x = udata.{0} in
  let y = udata.{1} in
  let z = udata.{2} in

  gdata.{0} <- (one /. three) *. cos((y -. one) *. z) +. (one /. six);
  gdata.{1} <- (one /. nine) *. sqrt(x *. x +. sin(z) +. oneptzerosix) +. ptnine;
  gdata.{2} <- -.(one /. twenty) *. exp(-. x *. (y -. one))
               -. (ten *. pi -. three) /. sixty

(* -----------------------------------------------------------------------------
 * Check the solution of the nonlinear system and return PASS or FAIL
 * ---------------------------------------------------------------------------*)
let check_ans (data : RealArray.t) tol =
  (* print the solution *)
  printf "Computed solution:\n";
  printf "    x = %g\n" data.{0};
  printf "    y = %g\n" data.{1};
  printf "    z = %g\n" data.{2};

  (* solution error *)
  let ex = abs_float (data.{0} -. xtrue) in
  let ey = abs_float (data.{1} -. ytrue) in
  let ez = abs_float (data.{2} -. ztrue) in

  (* print the solution error *)
  printf "Solution error:\n";
  printf "    ex = %g\n" ex;
  printf "    ey = %g\n" ey;
  printf "    ez = %g\n" ez;

  let tol = tol *. ten in
  if ex > tol || ey > tol || ez > tol then failwith "FAIL";

  printf "PASS\n"

(* -----------------------------------------------------------------------------
 * Main program
 * ---------------------------------------------------------------------------*)
let main () =
  let a_tol = ref (100.0 *. sqrt Config.unit_roundoff) in
  let a_mxiter = ref (if Sundials_impl.Version.lt580 then 10 else 30) in
  let a_maa = ref 0 in          (* no acceleration *)
  let a_delay_aa = ref 0 in     (* no delay *)
  let a_damping_fp = ref 1.0 in (* no FP damping *)
  let a_damping_aa = ref 1.0 in (* no damping *)
  let a_orth_aa = ref 0 in      (* MGS *)
  Arg.(parse
    [
      "--tol", Set_float a_tol,
      "nonlinear solver tolerance\n";

      "--maxiter", Set_int a_mxiter,
      "max number of nonlinear iterations\n";

      "--m_aa", Set_int a_maa,
      "number of Anderson acceleration vectors\n";

      "--delay_aa", Set_int a_delay_aa,
      "Anderson acceleration delay\n";

      "--damping_fp", Set_float a_damping_fp,
      "fixed point damping parameter\n";

      "--damping_aa", Set_float a_damping_aa,
      "Anderson acceleration damping parameter\n";

      "---orth_aa", Set_int a_orth_aa,
      "Anderson acceleration orthogonalization method\n"
    ] (fun _ -> ()) "\n Command line options:\n");

  let tol = !a_tol in
  let mxiter = !a_mxiter in
  let maa = !a_maa in
  let delay_aa = !a_delay_aa in
  let damping_fp = !a_damping_fp in
  let damping_aa = !a_damping_aa in
  let orthaa = orthaa_of_int !a_orth_aa in

  (* -------------------------
   * Print problem description
   * ------------------------- *)

  printf "Solve the nonlinear system:\n";
  printf "    3x - cos((y-1)z) - 1/2 = 0\n";
  printf "    x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0\n";
  printf "    exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0\n";
  printf "Analytic solution:\n";
  printf "    x = %g\n" xtrue;
  printf "    y = %g\n" ytrue;
  printf "    z = %g\n" ztrue;
  printf "Solution method: Anderson accelerated fixed point iteration.\n";
  if ge600 then begin
    printf "    tolerance    = %g\n" tol;
    printf "    max iters    = %d\n" mxiter;
    printf "    m_aa         = %d\n" maa;
    printf "    delay_aa     = %d\n" delay_aa;
    printf "    damping_aa   = %g\n" damping_aa;
    printf "    damping_fp   = %g\n" damping_fp;
    printf "    orth routine = %d\n" (int_of_orthaa orthaa)
  end else if ge580 then begin
    printf "    tolerance  = %g\n" tol;
    printf "    max iters  = %d\n" mxiter;
    printf "    m_aa       = %d\n" maa;
    printf "    delay_aa   = %d\n" delay_aa;
    printf "    damping_aa = %g\n" damping_aa;
    printf "    damping_fp = %g\n" damping_fp
  end else begin
    printf "    tolerance = %g\n" tol;
    printf "    max iters = %d\n" mxiter;
    printf "    accel vec = %d\n" maa;
    printf "    damping   = %g\n" damping_aa
  end;

  (* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- *)

  let u = Nvector_serial.make neq 0.0 in
  let scale = Nvector.clone u in

  (* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * ----------------------------------------- *)

  (* Set number of prior residuals used in Anderson acceleration *)
  (* Set maximum number of iterations *)
  let kmem = Kinsol.init ~max_iters:mxiter ~orthaa ~maa fpfunction u in

  (* -------------------
   * Set optional inputs
   * ------------------- *)

  (* Specify stopping tolerance based on residual *)
  Kinsol.set_func_norm_tol kmem tol;

  (* Set Anderson acceleration damping parameter *)
  Kinsol.set_damping_aa kmem damping_aa;

  (* -------------
   * Initial guess
   * ------------- *)

  (* Get vector data array *)
  let data = Nvector.unwrap u in
  data.{0} <-  ptone;
  data.{1} <-  ptone;
  data.{2} <- -. ptone;

  (* ----------------------------
   * Call KINSol to solve problem
   * ---------------------------- *)

  (* No scaling used *)
  Nvector_serial.Ops.const one scale;

  (* Call main solver *)
  ignore (Kinsol.(solve
                   kmem         (* KINSol memory block *)
                   u            (* initial guess on input; solution vector *)
                   FixedPoint   (* global strategy choice *)
                   scale        (* scaling vector, for the variable cc *)
                   scale));     (* scaling vector for function values fval *)

  (* ------------------------------------
   * Get solver statistics
   * ------------------------------------ *)

  (* get solver stats *)
  let nni = Kinsol.get_num_nonlin_solv_iters kmem in
  let nfe = Kinsol.get_num_func_evals kmem in

  printf "\nFinal Statistics:\n";
  printf "Number of nonlinear iterations: %6d\n" nni;
  printf "Number of function evaluations: %6d\n" nfe;

  (* ------------------------------------
   * Print solution and check error
   * ------------------------------------ *)

  (* check solution *)
  check_ans (Nvector.unwrap u) tol

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

