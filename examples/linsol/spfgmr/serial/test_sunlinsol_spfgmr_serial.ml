(*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This is the testing routine to check the SUNLinSol SPFGMR module
 * implementation.
 * -----------------------------------------------------------------
 *)

module LS = Sundials.LinearSolver
module Matrix = Sundials.Matrix
module RealArray = Sundials.RealArray

let printf = Format.printf
let (+=) r x = r := !r + x
let fneq a b tol = if Sundials.Util.compare_float ~tol a b then 0 else 1

(* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*)
module SPFGMR_tests : (Test_linsol.LINSOL_TESTS
                         with type nd = Nvector_serial.data
                          and type nk = Nvector_serial.kind) = struct
  type nd = Nvector_serial.data
  type nk = Nvector_serial.kind
  type nvec = (nd, nk) Nvector.t

  let check_vector x y tol =
    let xdata, ydata = Nvector.(unwrap x, unwrap y) in
    let nx, ny = Sundials.RealArray.(length xdata, length ydata) in
    if nx <> ny then
      (printf ">>> ERROR: check_vector: Different data array lengths @\n"; true)
    else begin
      let failure = ref 0 in
      for i = 0 to nx - 1 do
        failure += fneq xdata.{i} ydata.{i} tol
      done;
      if !failure > 0 then begin
        printf "Check_vector failures:@\n";
        for i = 0 to nx - 1 do
          if fneq xdata.{i} ydata.{i} tol <> 0 then
            printf "  xdata.{%d} = %g != %g (err = %g)@\n" i
                   xdata.{i} ydata.{i} (abs_float (xdata.{i} -. ydata.{i}))
        done;
        true
      end
      else false
    end

end

module Test = Test_linsol.Test (SPFGMR_tests) (Nvector_serial.Ops)

let sync_device () = ()

(* user data structure *)
type user_data = {
  n : int;          (* problem size *)
  d : RealArray.t;  (* matrix diagonal *)
  s1 : RealArray.t; (* scaling vectors supplied to SPFGMR *)
  s2 : RealArray.t;
}

(* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*)

(* matrix-vector product  *)
let atimes { n; s1; s2; _ } v z =
  (* perform product at the left domain boundary (note: v is zero at the boundary)*)
  z.{0} <- (5.0 *. v.{0} *. s2.{0} -. v.{1} *. s2.{1}) /. s1.{0};

  (* iterate through interior of local domain, performing product *)
  for i=1 to n-2 do
    z.{i} <- ((-.v.{i-1}) *. s2.{i-1} +. 5.0 *. v.{i} *. s2.{i} -. v.{i+1} *. s2.{i+1}) /. s1.{i}
  done;

  (* perform product at right boundary (note: v is zero at the boundary)*)
  z.{n-1} <- ((-.v.{n-2}) *. s2.{n-2} +. 5.0 *. v.{n-1} *. s2.{n-1}) /. s1.{n-1}

(* preconditioner setup -- nothing to do here since everything is already stored *)
let psetup _ () = ()

(* preconditioner solve *)
let psolve { d; n } r z _tol _lr =
  (* iterate through domain, performing Jacobi solve *)
  for i=0 to n-1 do
    z.{i} <- r.{i} /. d.{i}
  done

(* uniform random number generator in .{0,1} *)
let urand () =
  float_of_int (Sundials.Util.rand ())
  /. float_of_int (Sundials.Util.rand_max)

(* ----------------------------------------------------------------------
 * SUNLinSol_SPFGMR Linear Solver Testing Routine
 *
 * We run multiple tests to exercise this solver:
 * 1. simple tridiagonal system (no preconditioning)
 * 2. simple tridiagonal system (Jacobi preconditioning)
 * 3. tridiagonal system w/ scale vector s1 (no preconditioning)
 * 4. tridiagonal system w/ scale vector s1 (Jacobi preconditioning)
 * 5. tridiagonal system w/ scale vector s2 (no preconditioning)
 * 6. tridiagonal system w/ scale vector s2 (Jacobi preconditioning)
 *
 * Note: We construct a tridiagonal matrix Ahat, a random solution xhat,
 *       and a corresponding rhs vector bhat = Ahat*xhat, such that each
 *       of these is unit-less.  To test row/column scaling, we use the
 *       matrix A = S1-inverse Ahat S2, rhs vector b = S1-inverse bhat,
 *       and solution vector x = (S2-inverse) xhat; hence the linear
 *       system has rows scaled by S1-inverse and columns scaled by S2,
 *       where S1 and S2 are the diagonal matrices with entries from the
 *       vectors s1 and s2, the 'scaling' vectors supplied to SPFGMR
 *       having strictly positive entries.  When this is combined with
 *       preconditioning, assume that Phat is the desired preconditioner
 *       for Ahat, then our preconditioning matrix P \approx A should be
 *         left prec:  P-inverse \approx S1-inverse Ahat-inverse S1
 *         right prec:  P-inverse \approx S2-inverse Ahat-inverse S2.
 *       Here we use a diagonal preconditioner D, so the S*-inverse
 *       and S* in the product cancel one another.
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in
  let passfail = ref 0 in

  (* check inputs: local problem size, timing flag *)
  if Array.length Sys.argv < 6 then
    (printf "ERROR: FIVE (5) Inputs required:@\n";
     printf "  Problem size should be >0@\n";
     printf "  Gram-Schmidt orthogonalization type should be 1 or 2@\n";
     printf "  Maximum Krylov subspace dimension should be >0@\n";
     printf "  Solver tolerance should be >0@\n";
     printf "  timing output flag should be 0 or 1 @\n";
     exit (-1));

  let problem_size = int_of_string Sys.argv.(1) in
  if problem_size <= 0 then begin
    printf "ERROR: Problem size must be a positive integer@\n";
    exit 1
  end;

  let gstype = match int_of_string Sys.argv.(2) with
    | 1 -> LS.Iterative.ModifiedGS
    | 2 -> LS.Iterative.ClassicalGS
    | _ -> printf "ERROR: Gram-Schmidt orthogonalization type must be either 1 or 2@\n"; exit 1
  in

  let maxl = int_of_string Sys.argv.(3) in
  if maxl <= 0 then begin
    printf "ERROR: Maximum Krylov subspace dimension must be a positive integer@\n";
    exit 1
  end;

  let tol = float_of_string Sys.argv.(4) in
  if tol <= 0.0 then begin
    printf "ERROR: Solver tolerance must be a positive real number@\n";
    exit 1
  end;

  let print_timing = int_of_string Sys.argv.(5) in
  let _ = Test.set_timing (print_timing <> 0) in

  printf "@\nSPFGMR linear solver test:@\n";
  printf "  problem size = %d@\n" problem_size;
  printf "  Gram-Schmidt orthogonalization type = %d@\n"
    LS.Iterative.(match gstype with ModifiedGS -> 1 | ClassicalGS -> 2);
  printf "  Maximum Krylov subspace dimension = %d@\n" maxl;
  printf "  Solver Tolerance = %g@\n" tol;
  printf "  timing output flag = %d@\n@\n" print_timing;

  (* Create vectors *)
  let xvec = Nvector_serial.make problem_size 0.0
  and xhat = Nvector_serial.wrap
               (* Fill xhat vector with uniform random data in .{1,2} *)
               (RealArray.init problem_size (fun _ -> 1.0 +. urand ()))
  and bvec = Nvector_serial.make problem_size 0.0
  in
  let prob_data = {
    d  = (* Fill Jacobi vector with matrix diagonal *)
         RealArray.make problem_size 5.0;
    s1 = RealArray.make problem_size 0.0;
    s2 = RealArray.make problem_size 0.0;
    n  = problem_size }
  in
  let x = Nvector_serial.unwrap xvec in
  let b = Nvector_serial.unwrap bvec in
  let s1vec = Nvector_serial.wrap prob_data.s1 in
  let s2vec = Nvector_serial.wrap prob_data.s2 in

  (* Create SPFGMR linear solver *)
  let ls = LS.Iterative.spfgmr ~maxl xvec in
  LS.Iterative.(set_prec_type ls PrecRight);

  (* Run Tests *)
  fails += Test.sunlinsolgettype ls LS.Iterative 0;
  fails += Test.sunlinsolgetid ls LS.Spfgmr 0;
  fails += Test.sunlinsolsetatimes ls (atimes prob_data) 0;
  fails += Test.sunlinsolsetpreconditioner ls
             (psetup prob_data) (psolve prob_data) 0;
  fails += Test.sunlinsolsetscalingvectors ls s1vec s2vec 0;
  fails += Test.sunlinsolsetzeroguess ls 0;
  fails += Test.sunlinsolinitialize ls 0;
  fails += Test.sunlinsolspace ls 0;
  (try LS.Iterative.set_gs_type ls gstype with _ -> fails += 1);
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_SPFGMR module failed %d initialization tests@\n@\n" !fails;
    exit 1
  end else
    printf "SUCCESS: SUNLinSol_SPFGMR module passed all initialization tests@\n@\n";

  (*** Test 1: simple Poisson-like solve (no preconditioning) ***)

  (* set scaling vectors *)
  Nvector_serial.Ops.const 1.0 s1vec;
  Nvector_serial.Ops.const 1.0 s2vec;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.div xhat s2vec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecNone) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolnumiters ls 0;
  fails += Test.sunlinsolresnorm ls 0;
  fails += Test.sunlinsolresid ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_SPFGMR module, problem 1, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_SPFGMR module, problem 1, passed all tests@\n@\n";

  (*** Test 2: simple Poisson-like solve (Jacobi preconditioning) ***)

  (* set scaling vectors *)
  Nvector_serial.Ops.const 1.0 s1vec;
  Nvector_serial.Ops.const 1.0 s2vec;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.div xhat s2vec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecRight) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolnumiters ls 0;
  fails += Test.sunlinsolresnorm ls 0;
  fails += Test.sunlinsolresid ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_SPFGMR module, problem 2, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_SPFGMR module, problem 2, passed all tests@\n@\n";


  (*** Test 3: Poisson-like solve w/ scaled rows (no preconditioning) ***)

  (* set scaling vectors *)
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s1;
  Nvector_serial.Ops.const 1.0 s2vec;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.div xhat s2vec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecNone) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolnumiters ls 0;
  fails += Test.sunlinsolresnorm ls 0;
  fails += Test.sunlinsolresid ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_SPFGMR module, problem 3, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_SPFGMR module, problem 3, passed all tests@\n@\n";


  (*** Test 4: Poisson-like solve w/ scaled rows (Jacobi preconditioning) ***)

  (* set scaling vectors *)
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s1;
  Nvector_serial.Ops.const 1.0 s2vec;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.div xhat s2vec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecRight) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolnumiters ls 0;
  fails += Test.sunlinsolresnorm ls 0;
  fails += Test.sunlinsolresid ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_SPFGMR module, problem 4, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_SPFGMR module, problem 4, passed all tests@\n@\n";

  (*** Test 5: Poisson-like solve w/ scaled columns (no preconditioning) ***)

  (* set scaling vectors *)
  Nvector_serial.Ops.const 1.0 s1vec;
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s2;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.div xhat s2vec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecNone) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolnumiters ls 0;
  fails += Test.sunlinsolresnorm ls 0;
  fails += Test.sunlinsolresid ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_SPFGMR module, problem 5, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_SPFGMR module, problem 5, passed all tests@\n@\n";

  (*** Test 6: Poisson-like solve w/ scaled columns (Jacobi preconditioning) ***)

  (* set scaling vector, Jacobi solver vector *)
  Nvector_serial.Ops.const 1.0 s1vec;
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s2;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.div xhat s2vec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecRight) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true 0;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolnumiters ls 0;
  fails += Test.sunlinsolresnorm ls 0;
  fails += Test.sunlinsolresid ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_SPFGMR module, problem 6, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_SPFGMR module, problem 6, passed all tests@\n@\n";

  !passfail

let _ = main ()

