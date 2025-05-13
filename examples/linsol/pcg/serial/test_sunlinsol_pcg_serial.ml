(*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This is the testing routine to check the SUNLinSol PCG module
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
module PCG_tests : (Test_linsol.LINSOL_TESTS
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

module Test = Test_linsol.Test (PCG_tests) (Nvector_serial.Ops)

let sync_device () = ()

(* user data structure *)
type user_data = {
  n : int;         (* problem size *)
  d : RealArray.t; (* matrix diagonal *)
  s : RealArray.t; (* scaling vector supplied to PCG *)
}

(* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*)

(* matrix-vector product  *)
let atimes { n; s; _ } v z =
  (* perform product at left boundary (note: v is zero at the boundary)*)
  z.{0} <- (5.0 *. v.{0} /. s.{0} -. v.{1} /. s.{1}) /. s.{0};

  (* iterate through interior of domain, performing product *)
  for i=1 to n-2 do
    z.{i} <- ((-.v.{i-1}) /. s.{i-1} +. 5.0 *. v.{i} /. s.{i} -. v.{i+1} /. s.{i+1}) /. s.{i}
  done;

  (* perform product at right boundary (note: v is zero at the boundary)*)
  z.{n-1} <- ((-.v.{n-2}) /. s.{n-2} +. 5.0 *. v.{n-1} /. s.{n-1}) /. s.{n-1}

(* preconditioner setup -- nothing to do here since everything is already stored *)
let psetup _ () = ()

(* preconditioner solve *)
let psolve { d; s; n } r z _tol _lr =
  (* iterate through domain, performing Jacobi solve *)
  for i=0 to n-1 do
    z.{i} <- s.{i} *. s.{i} *. r.{i} /. d.{i}
  done

(* uniform random number generator in .{0,1} *)
let urand () =
  float_of_int (Sundials.Util.rand ())
  /. float_of_int (Sundials.Util.rand_max)

(* ----------------------------------------------------------------------
 * SUNOCG Linear Solver Testing Routine
 *
 * We run multiple tests to exercise this solver:
 * 1. simple tridiagonal system (no preconditioning)
 * 2. simple tridiagonal system (Jacobi preconditioning)
 * 3. tridiagonal system w/ scale vector s (no preconditioning)
 * 4. tridiagonal system w/ scale vector s (Jacobi preconditioning)
 *
 * Note: We construct a tridiagonal matrix Ahat, a random solution
 *       xhat, and a corresponding rhs vector bhat = Ahat*xhat, such
 *       that each of these is unit-less.  To test scaling, we use
 *       the matrix
 *             A = (S-inverse) Ahat (S-inverse),
 *       solution vector
 *             x = S xhat;
 *       and construct b = A*x.  Hence the linear system has both rows
 *       and columns scaled by (S-inverse), where S is the diagonal
 *       matrix with entries from the vector s, the 'scaling' vector
 *       supplied to PCG having strictly positive entries.
 *
 *       When this is combined with preconditioning, we construct
 *       P \approx (A-inverse) by taking a unit-less preconditioner
 *       Phat \approx (Ahat-inverse), and constructing the operator
 *       P via
 *             P = S Phat S \approx S (Ahat-inverse) S = A-inverse
 *       We apply this via the steps:
 *             z = Pr = S Phat S r
 *       Since both S and Phat are diagonal matrices, this is
 *       equivalent to
 *             z(i) = s(i)^2 Phat(i) r(i)
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in
  let passfail = ref 0 in

  (* check inputs: local problem size, timing flag *)
  if Array.length Sys.argv < 5 then
    (printf "ERROR: FOUR (4) Inputs required:@\n";
     printf "  Problem size should be >0@\n";
     printf "  Maximum Krylov subspace dimension should be >0@\n";
     printf "  Solver tolerance should be >0@\n";
     printf "  timing output flag should be 0 or 1 @\n";
     exit (-1));

  let problem_size = int_of_string Sys.argv.(1) in
  if problem_size <= 0 then begin
    printf "ERROR: Problem size must be a positive integer@\n";
    exit 1
  end;

  let maxl = int_of_string Sys.argv.(2) in
  if maxl <= 0 then begin
    printf "ERROR: Maximum Krylov subspace dimension must be a positive integer@\n";
    exit 1
  end;

  let tol = float_of_string Sys.argv.(3) in
  if tol <= 0.0 then begin
    printf "ERROR: Solver tolerance must be a positive real number@\n";
    exit 1
  end;

  let print_timing = int_of_string Sys.argv.(4) in
  let _ = Test.set_timing (print_timing <> 0) in

  printf "@\nPCG linear solver test:@\n";
  printf "  Problem size = %d@\n" problem_size;
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
    d = (* Fill Jacobi vector with matrix diagonal *)
        RealArray.make problem_size 5.0;
    s = RealArray.make problem_size 0.0;
    n = problem_size }
  in
  let x = Nvector_serial.unwrap xvec in
  let b = Nvector_serial.unwrap bvec in
  let svec = Nvector_serial.wrap prob_data.s in
  let null = Nvector_serial.make 0 0.0 in

  (* Create PCG linear solver *)
  let ls = LS.Iterative.pcg ~maxl xvec in
  LS.Iterative.(set_prec_type ls PrecRight);

  (* Run Tests *)
  fails += Test.sunlinsolgettype ls LS.Iterative 0;
  fails += Test.sunlinsolgetid ls LS.Pcg 0;
  fails += Test.sunlinsolsetatimes ls (atimes prob_data) 0;
  fails += Test.sunlinsolsetpreconditioner ls
             (psetup prob_data) (psolve prob_data) 0;
  fails += Test.sunlinsolsetscalingvectors ls svec null 0;
  fails += Test.sunlinsolsetzeroguess ls 0;
  fails += Test.sunlinsolinitialize ls 0;
  fails += Test.sunlinsolspace ls 0;
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_PCG module failed %d initialization tests@\n@\n" !fails;
    exit 1
  end else
    printf "SUCCESS: SUNLinSol_PCG module passed all initialization tests@\n@\n";

  (*** Test 1: simple Poisson-like solve (no preconditioning) ***)

  (* set scaling vector *)
  Nvector_serial.Ops.const 1.0 svec;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.prod xhat svec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run test with this setup *)
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
    printf "FAIL: SUNLinSol_PCG module, problem 1, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_PCG module, problem 1, passed all tests@\n@\n";

  (*** Test 2: simple Poisson-like solve (Jacobi preconditioning) ***)

  (* set scaling vector *)
  Nvector_serial.Ops.const 1.0 svec;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.prod xhat svec xvec;

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
    printf "FAIL: SUNLinSol_PCG module, problem 2, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_PCG module, problem 2, passed all tests@\n@\n";


  (*** Test 3: Poisson-like solve w/ scaling (no preconditioning) ***)

  (* set scaling vector *)
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.prod xhat svec xvec;

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
    printf "FAIL: SUNLinSol_PCG module, problem 3, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_PCG module, problem 3, passed all tests@\n@\n";


  (*** Test 4: Poisson-like solve w/ scaling (Jacobi preconditioning) ***)

  (* set scaling vectors *)
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s;

  (* Fill x vector with scaled version *)
  Nvector_serial.Ops.prod xhat svec xvec;

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
    printf "FAIL: SUNLinSol_PCG module, problem 4, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else
    printf "SUCCESS: SUNLinSol_PCG module, problem 4, passed all tests@\n@\n";

  !passfail

let _ = main ()

