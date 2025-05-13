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
                      with type nd = Nvector_parallel.data
                       and type nk = Nvector_parallel.kind) = struct
  type nd = Nvector_parallel.data
  type nk = Nvector_parallel.kind
  type nvec = (nd, nk) Nvector.t

  let check_vector x y tol =
    let (xdata, _, _), (ydata, _, _) = Nvector.(unwrap x, unwrap y) in
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

module Test = Test_linsol.Test (PCG_tests) (Nvector_parallel.Ops)

let sync_device () = ()

(* user data structure *)
type user_data = {
  nloc   : int;              (* local problem size *)
  d      : RealArray.t;      (* matrix diagonal *)
  s      : RealArray.t;      (* scaling vector supplied to PCG *)
  comm   : Mpi.communicator; (* communicator object *)
  myid   : int;              (* MPI process ID *)
  nprocs : int;              (* total number of MPI processes *)
}

(* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*)

let float_nbytes = Bytes.length (Marshal.to_bytes 0.0 [])

(* matrix-vector product  *)
let atimes { nloc = n; s; comm; myid; nprocs; _ } (v, _, _) (z, _, _) =

  (* send/recv boundary data with neighbors *)
  let vsL = v.{0} /. s.{0} in
  let vsR = v.{n-1} /. s.{n-1} in

  let recv_req_l, _send_req_l =
    if myid > 0 then (* left neighbor exists *)
      Mpi.ireceive float_nbytes (myid - 1) Mpi.any_tag comm,
      Mpi.isend vsL (myid - 1) 0 comm
    else Mpi.null_request, Mpi.null_request
  in
  let recv_req_r, _send_req_r =
    if myid < nprocs-1 then (* right neighbor exists *)
      Mpi.ireceive float_nbytes (myid + 1) Mpi.any_tag comm,
      Mpi.isend vsR (myid + 1) 1 comm
    else Mpi.null_request, Mpi.null_request
  in

  (* iterate through interior of local domain, performing product *)
  for i=1 to n-2 do
    z.{i} <- ((-.v.{i-1}) /. s.{i-1} +. 5.0 *. v.{i} /. s.{i} -. v.{i+1} /. s.{i+1}) /. s.{i}
  done;

  (* wait on neighbor data to arrive *)
  let vL = if myid > 0 then Mpi.wait_receive recv_req_l (* left neighbor exists *) else 0.0
  in
  let vR = if myid < nprocs-1 then Mpi.wait_receive recv_req_r (* right neighbor exists *) else 0.0
  in

  (* perform product at subdomain boundaries (note: vL/vR are zero at boundary)*)
  z.{0} <- ((-. vL) +. 5.0 *. v.{0} /. s.{0} -. v.{1} /. s.{1}) /. s.{0};
  z.{n-1} <- ((-.v.{n-2}) /. s.{n-2} +. 5.0 *. v.{n-1} /. s.{n-1} -. vR) /. s.{n-1}

(* preconditioner setup -- nothing to do here since everything is already stored *)
let psetup _ () = ()

(* preconditioner solve *)
let psolve { d; s; nloc = n } (r, _, _) (z, _, _) _tol _lr =
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

  (* Set up MPI environment *)
  let comm   = Mpi.comm_world in
  let nprocs = Mpi.comm_size comm in
  let myid   = Mpi.comm_rank comm in

  (* check inputs: local problem size, timing flag *)
  if Array.length Sys.argv < 5 then
    (printf "ERROR: FOUR (4) Inputs required:@\n";
     printf "  Local problem size should be >0@\n";
     printf "  Maximum Krylov subspace dimension should be >0@\n";
     printf "  Solver tolerance should be >0@\n";
     printf "  timing output flag should be 0 or 1 @\n";
     exit (-1));

  let nloc = int_of_string Sys.argv.(1) in
  if nloc <= 0 then begin
    printf "ERROR: Local problem size must be a positive integer@\n";
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

  if myid = 0 then begin
    printf "@\nPCG linear solver test:@\n";
    printf "  nprocs = %d@\n" nprocs;
    printf "  local/global problem sizes = %d/%d@\n" nloc (nprocs * nloc);
    printf "  Maximum Krylov subspace dimension = %d@\n" maxl;
    printf "  Solver Tolerance = %g@\n" tol;
    printf "  timing output flag = %d@\n@\n" print_timing;
  end;

  (* Create vectors *)
  let xvec = Nvector_parallel.make nloc (nprocs * nloc) comm 0.0
  and xhat = Nvector_parallel.wrap
               (* Fill xhat vector with uniform random data in .{1,2} *)
               (RealArray.init nloc (fun _ -> 1.0 +. urand ()),
                nprocs * nloc, comm)
  and bvec = Nvector_parallel.make nloc (nprocs * nloc) comm 0.0
  in
  let prob_data = {
    d    = (* Fill Jacobi vector with matrix diagonal *)
           RealArray.make nloc 5.0;
    s    = RealArray.make nloc 0.0;
    nloc; comm; nprocs; myid }
  in
  let x = Nvector_parallel.unwrap xvec in
  let b = Nvector_parallel.unwrap bvec in
  let svec = Nvector_parallel.wrap (prob_data.s, nprocs * nloc, comm) in
  let null = Nvector_parallel.make 0 0 comm 0.0 in

  (* Create PCG linear solver *)
  let ls = LS.Iterative.pcg ~maxl xvec in
  LS.Iterative.(set_prec_type ls PrecRight);

  (* Run Tests *)
  fails += Test.sunlinsolgettype ls LS.Iterative myid;
  fails += Test.sunlinsolgetid ls LS.Pcg myid;
  fails += Test.sunlinsolsetatimes ls (atimes prob_data) myid;
  fails += Test.sunlinsolsetpreconditioner ls
             (psetup prob_data) (psolve prob_data) myid;
  fails += Test.sunlinsolsetscalingvectors ls svec null myid;
  fails += Test.sunlinsolsetzeroguess ls myid;
  fails += Test.sunlinsolinitialize ls myid;
  fails += Test.sunlinsolspace ls myid;
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_PCG module failed %d initialization tests@\n@\n" !fails;
    exit 1
  end else if myid = 0 then
    printf "SUCCESS: SUNLinSol_PCG module passed all initialization tests@\n@\n";

  (*** Test 1: simple Poisson-like solve (no preconditioning) ***)

  (* set scaling vector *)
  Nvector_parallel.Ops.const 1.0 svec;

  (* Fill x vector with scaled version *)
  Nvector_parallel.Ops.prod xhat svec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run test with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecNone) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false myid;
  fails += Test.sunlinsollastflag ls myid;
  fails += Test.sunlinsolnumiters ls myid;
  fails += Test.sunlinsolresnorm ls myid;
  fails += Test.sunlinsolresid ls myid;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_PCG module, problem 1, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else if myid = 0 then
    printf "SUCCESS: SUNLinSol_PCG module, problem 1, passed all tests@\n@\n";

  (*** Test 2: simple Poisson-like solve (Jacobi preconditioning) ***)

  (* set scaling vector *)
  Nvector_parallel.Ops.const 1.0 svec;

  (* Fill x vector with scaled version *)
  Nvector_parallel.Ops.prod xhat svec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecRight) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false myid;
  fails += Test.sunlinsollastflag ls myid;
  fails += Test.sunlinsolnumiters ls myid;
  fails += Test.sunlinsolresnorm ls myid;
  fails += Test.sunlinsolresid ls myid;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_PCG module, problem 2, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else if myid = 0 then
    printf "SUCCESS: SUNLinSol_PCG module, problem 2, passed all tests@\n@\n";


  (*** Test 3: Poisson-like solve w/ scaling (no preconditioning) ***)

  (* set scaling vector *)
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s;

  (* Fill x vector with scaled version *)
  Nvector_parallel.Ops.prod xhat svec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecNone) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false myid;
  fails += Test.sunlinsollastflag ls myid;
  fails += Test.sunlinsolnumiters ls myid;
  fails += Test.sunlinsolresnorm ls myid;
  fails += Test.sunlinsolresid ls myid;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_PCG module, problem 3, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else if myid = 0 then
    printf "SUCCESS: SUNLinSol_PCG module, problem 3, passed all tests@\n@\n";


  (*** Test 4: Poisson-like solve w/ scaling (Jacobi preconditioning) ***)

  (* set scaling vectors *)
  RealArray.map (fun _ -> 1.0 +. 1000.0 *. urand ()) prob_data.s;

  (* Fill x vector with scaled version *)
  Nvector_parallel.Ops.prod xhat svec xvec;

  (* Fill b vector with result of matrix-vector product *)
  fails := 0;
  atimes prob_data x b;

  (* Run tests with this setup *)
  (try LS.Iterative.(set_prec_type ls PrecRight) with _ -> fails += 1);
  fails += Test.sunlinsolsetup ls None myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol true myid;
  fails += Test.sunlinsolsolve ls None xvec bvec tol false myid;
  fails += Test.sunlinsollastflag ls myid;
  fails += Test.sunlinsolnumiters ls myid;
  fails += Test.sunlinsolresnorm ls myid;
  fails += Test.sunlinsolresid ls myid;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol_PCG module, problem 4, failed %d tests@\n@\n" !fails;
    passfail += 1
  end else if myid = 0 then
    printf "SUCCESS: SUNLinSol_PCG module, problem 4, passed all tests@\n@\n";

  (* check if any other process failed *)
  Mpi.(allreduce_int 1 Max comm)

let _ = main ()

