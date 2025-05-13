
(*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2025.
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file contains the prototypes for functions to
 * test SUNLinearSolver module implementations.
 * -----------------------------------------------------------------
 *)

open Sundials
module LS = LinearSolver

module type LINSOL_TESTS = sig
  type nd
  type nk

  type nvec = (nd, nk) Nvector.t
  val check_vector : nvec -> nvec -> float -> bool
end

module Test (LT : LINSOL_TESTS) (NV : Nvector.NVECTOR_OPS with type t = LT.nvec) =
struct

(* private functions *)

let printf = Format.printf

let print_time_flag = ref false

let print_time format time =
  if !print_time_flag then
    printf format time

let get_time = Sundials.Util.get_monotonic_time

let set_timing b =
  print_time_flag := b

(* ----------------------------------------------------------------------
 * SUNLinSolGetType Test
 * --------------------------------------------------------------------*)
let sunlinsolgettype lsolv suntype myid =
  let start_time = get_time () in
  let mysuntype = LS.get_type lsolv in
  let stop_time = get_time () in
  if suntype <> mysuntype then begin
    printf ">>> FAILED test -- SUNLinSolGetType, Proc %d @\n" myid;
    print_time "    SUNLinSolGetType Time: %22.15e @\n @\n" (stop_time -. start_time)
  end else if myid = 0 then begin
    printf "    PASSED test -- SUNLinSolGetType @\n";
    print_time "    SUNLinSolGetType Time: %22.15e @\n @\n" (stop_time -. start_time)
  end;
  if suntype = mysuntype then 0 else 1

(* ----------------------------------------------------------------------
 * SUNLinSolGetID Test
 * --------------------------------------------------------------------*)
let sunlinsolgetid lsolv sunid myid =
  let start_time = get_time () in
  let mysunid = LS.get_id lsolv in
  let stop_time = get_time () in
  if sunid <> mysunid then begin
    printf ">>> FAILED test -- SUNLinSolGetID, Proc %d @\n" myid;
    print_time "    SUNLinSolGetID Time: %22.15e @\n @\n" (stop_time -. start_time)
  end else if myid = 0 then begin
    printf "    PASSED test -- SUNLinSolGetID @\n";
    print_time "    SUNLinSolGetID Time: %22.15e @\n @\n" (stop_time -. start_time)
  end;
  if sunid = mysunid then 0 else 1

(* ----------------------------------------------------------------------
 * SUNLinSolLastFlag Test
 * --------------------------------------------------------------------*)
let sunlinsollastflag lsolv myid =
  let start_time = get_time () in
  let lastflag = LS.get_last_flag lsolv in
  let stop_time = get_time () in
  if myid = 0 then begin
    printf "    PASSED test -- SUNLinSolLastFlag (%d) @\n" lastflag;
    print_time "    SUNLinSolLastFlag Time: %22.15e @\n @\n" (stop_time -. start_time)
  end;
  0

(* ----------------------------------------------------------------------
 * SUNLinSolSpace Test
 * --------------------------------------------------------------------*)
let sunlinsolspace lsolv myid =
  try
    (* call SUNLinSolSpace (failure based on output flag) *)
    let start_time = get_time () in
    let lenrw, leniw = try LS.get_work_space lsolv with _ -> begin
        printf ">>> FAILED test -- SUNLinSolSpace, Proc %d @\n" myid;
        let stop_time = get_time () in
        print_time "    SUNLinSolSpace Time: %22.15e @\n @\n" (stop_time -. start_time);
        raise Exit
      end
    in
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSpace, lenrw = %d, leniw = %d@\n" lenrw leniw;
      print_time "    SUNLinSolSpace Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolNumIters Test
 * --------------------------------------------------------------------*)
let sunlinsolnumiters lsolv myid =
  (* the only way to fail this test is if the function is NULL,
     which will cause a seg-fault *)
  let start_time = get_time () in
  let numiters = LS.get_num_iters lsolv in
  let stop_time = get_time () in
  if myid = 0 then begin
    printf "    PASSED test -- SUNLinSolNumIters (%d) @\n" numiters;
    print_time "    SUNLinSolNumIters Time: %22.15e @\n @\n" (stop_time -. start_time)
  end;
  0

(* ----------------------------------------------------------------------
 * SUNLinSolResNorm Test
 * --------------------------------------------------------------------*)
let sunlinsolresnorm lsolv myid =
  (* this test can fail if the function is NULL, which will cause a seg-fault *)
  let start_time = get_time () in
  let resnorm = LS.get_res_norm lsolv in
  let stop_time = get_time () in
  (* this test can also fail if the return value is negative *)
  if resnorm < 0.0 then
    printf ">>> FAILED test -- SUNLinSolResNorm returned %g on Proc %d @\n" resnorm myid
  else if myid = 0 then begin
    printf "    PASSED test -- SUNLinSolResNorm@\n";
    print_time "    SUNLinSolResNorm Time: %22.15e @\n @\n" (stop_time -. start_time)
  end;
  if resnorm >= 0.0 then 0 else 1

(* ----------------------------------------------------------------------
 * SUNLinSolResid Test
 * --------------------------------------------------------------------*)
let sunlinsolresid lsolv myid =
  try
    (* this test can fail if the function returns NULL *)
    let start_time = get_time () in
    let _resid = try LS.get_res_id lsolv with _ -> begin
        printf ">>> FAILED test -- SUNLinSolResid returned NULL N_Vector on Proc %d @\n" myid;
        raise Exit
      end
    in
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolResid@\n";
      print_time "    SUNLinSolResid Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolSetATimes Test
 * --------------------------------------------------------------------*)
let sunlinsolsetatimes lsolv atimes myid =
  try
    (* try calling SetATimes routine: should pass/fail based on expected input *)
    let start_time = get_time () in
    (try
       LS.set_atimes lsolv atimes
     with e -> begin
       printf ">>> FAILED test -- SUNLinSolSetATimes failed with %s on Proc %d @\n"
         (Printexc.to_string e) myid;
       raise Exit
     end);
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSetATimes @\n";
      print_time "    SUNLinSolSetATimes Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolSetPreconditioner
 * --------------------------------------------------------------------*)
let sunlinsolsetpreconditioner lsolv psetup psolve myid =
  try
    (* try calling SetPreconditioner routine: should pass/fail based on expected input *)
    let start_time = get_time () in
    (try
       LS.set_preconditioner lsolv psetup psolve
     with e -> begin
       printf ">>> FAILED test -- SUNLinSolSetPreconditioner failed with %s on Proc %d @\n"
         (Printexc.to_string e) myid;
       raise Exit
     end);
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSetPreconditioner @\n";
      print_time "    SUNLinSolSetPreconditioner Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolSetScalingVectors
 * --------------------------------------------------------------------*)
let sunlinsolsetscalingvectors lsolv s1 s2 myid =
  try
    (* try calling SetScalingVectors routine: should pass/fail based on expected input *)
    let start_time = get_time () in
    (try
       LS.set_scaling_vectors lsolv s1 s2
     with e -> begin
       printf ">>> FAILED test -- SUNLinSolSetScalingVectors failed with %s on Proc %d @\n"
         (Printexc.to_string e) myid;
       raise Exit
     end);
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSetScalingVectors @\n";
      print_time "    SUNLinSolSetScalingVectors Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolSetZeroGuess
 * --------------------------------------------------------------------*)
let sunlinsolsetzeroguess lsolv myid =
  try
    (* try calling SetZeroGuess routine: should pass/fail based on expected input *)
    let start_time = get_time () in
    (try
       LS.set_zero_guess lsolv true
     with e -> begin
       printf ">>> FAILED test -- SUNLinSolSetZeroGuess failed with %s on Proc %d @\n"
         (Printexc.to_string e) myid;
       raise Exit
     end);
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSetZeroGuess @\n";
      print_time "    SUNLinSolSetZeroGuess Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;

    (* try calling SetZeroGuess routine: should pass/fail based on expected input *)
    let start_time = get_time () in
    (try
       LS.set_zero_guess lsolv false
     with e -> begin
       printf ">>> FAILED test -- SUNLinSolSetZeroGuess failed with %s on Proc %d @\n"
         (Printexc.to_string e) myid;
       raise Exit
     end);
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSetZeroGuess @\n";
      print_time "    SUNLinSolSetZeroGuess Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolInitialize Test
 * --------------------------------------------------------------------*)
let sunlinsolinitialize lsolv myid =
  try
    let start_time = get_time () in
    (try
       LS.init lsolv
     with _ -> begin
       printf ">>> FAILED test -- SUNLinSolInitialize check, Proc %d @\n" myid;
       let stop_time = get_time () in
       print_time "    SUNLinSolInitialize Time: %22.15e @\n @\n" (stop_time -. start_time);
       raise Exit
     end);
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolInitialize @\n";
      print_time "    SUNLinSolInitialize Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolSetup Test
 *
 * This test must follow SUNLinSolInitialize
 * --------------------------------------------------------------------*)
let sunlinsolsetup lsolv oa myid =
  try
    let start_time = get_time () in
    (try
       LS.setup lsolv oa
     with _ -> begin
       printf ">>> FAILED test -- SUNLinSolSetup check, Proc %d @\n" myid;
       let stop_time = get_time () in
       print_time "    SUNLinSolSetup Time: %22.15e @\n @\n" (stop_time -. start_time);
       raise Exit
     end);
    let stop_time = get_time () in
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSetup @\n";
      print_time "    SUNLinSolSetup Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNLinSolSolve Test
 *
 * This test must follow SUNLinSolSetup.  Also, x must be the
 * solution to the linear system A*x = b (for the original A matrix);
 * while the 'A' that is supplied to this function should have been
 * 'setup' by the SUNLinSolSetup() function prior to this call.
 * --------------------------------------------------------------------*)
let sunlinsolsolve lsolv oa x b tol zeroguess myid =
  try
    (* clone to create solution vector *)
    let y = NV.clone x in

    (* set initial guess for the linear system *)
    if zeroguess then NV.const 0.0 y
    else NV.addconst x (sqrt Sundials.Config.unit_roundoff) y;

    (* signal if the initial guess is zero or non-zero *)
    (try
       LS.set_zero_guess lsolv zeroguess
     with e -> begin
       printf ">>> FAILED test -- SUNLinSolSetZeroGuess failed with %s on Proc %d @\n"
         (Printexc.to_string e) myid;
       raise Exit
     end);

    (* perform solve *)
    let start_time = get_time () in
    (try
      LinearSolver.solve lsolv oa y b tol
     with e -> begin
       printf ">>> FAILED test -- SUNLinSolSolve failed with %s on Proc %d @\n"
         (Printexc.to_string e) myid;
       raise Exit
     end);
    let stop_time = get_time () in

    (* Check solution *)
    let failure = LT.check_vector x y (10.0 *. tol) in
    if failure then begin
      printf ">>> FAILED test -- SUNLinSolSolve check, Proc %d @\n" myid;
      print_time "    SUNLinSolSolve Time: %22.15e @\n @\n" (stop_time -. start_time);
      raise Exit
    end;
    if myid = 0 then begin
      printf "    PASSED test -- SUNLinSolSolve @\n";
      print_time "    SUNLinSolSolve Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

end

