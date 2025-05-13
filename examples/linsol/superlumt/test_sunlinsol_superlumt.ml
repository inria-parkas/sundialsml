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
 * This is the testing routine to check the SUNLinSol SuperLUMT
 * module implementation.
 * -----------------------------------------------------------------
 *)

module LS = Sundials.LinearSolver
module Matrix = Sundials.Matrix

let f_rand () = float_of_int (Sundials.Util.rand ())
let f_rand_max = float_of_int (Sundials.Util.rand_max)
let unit_roundoff = Sundials.Config.unit_roundoff

let printf = Format.printf
let (+=) r x = r := !r + x
let fneq a b tol = if Sundials.Util.compare_float ~tol a b then 0 else 1

(* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*)
module Superlumt_tests : (Test_linsol.LINSOL_TESTS
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
            printf "  xdata[%d] = %g != %g (err = %g)@\n" i
                   xdata.{i} ydata.{i} (abs_float (xdata.{i} -. ydata.{i}))
        done;
        true
      end
      else false
    end

end

module Test = Test_linsol.Test (Superlumt_tests) (Nvector_serial.Ops)

let sync_device () = ()

let int_of_mattype (type s) : s Matrix.Sparse.sformat -> int = function
  | Matrix.Sparse.CSC -> 0
  | Matrix.Sparse.CSR -> 1

(* ----------------------------------------------------------------------
 * SUNLinSol_SuperLUMT Testing Routine
 * --------------------------------------------------------------------*)
let rec main () =
  (* check input and set matrix dimensions *)
  if Array.length Sys.argv < 5 then
    (printf "ERROR: FOUR (4) Inputs required: matrix size, matrix type (0 = CSC / 1 = CSR), num_threads, print timing @\n";
     exit (-1));

  let n = int_of_string Sys.argv.(1) in
  if n <= 0 then
    (printf "ERROR: matrix size must be a positive integer @\n";
     exit (-1));

  let mattype = int_of_string Sys.argv.(2) in
  if (mattype != 0) && (mattype != 1) then begin
    printf "ERROR: matrix type must be 0 or 1 @\n";
    exit (-1);
  end;

  let num_threads = int_of_string Sys.argv.(3) in
  if num_threads <= 0 then begin
    printf "ERROR: number_threads must be a positive integer @\n";
    exit (-1)
  end;

  let print_timing = int_of_string Sys.argv.(4) in
  let _ = Test.set_timing (print_timing <> 0) in

  if mattype = 0
  then main_with_type n num_threads Matrix.Sparse.CSC
  else main_with_type n num_threads Matrix.Sparse.CSR

and main_with_type : type s. int -> int -> s Matrix.Sparse.sformat -> int
  = fun n nthreads mattype ->
  let fails = ref 0 in

  printf "@\nSuperLUMT linear solver test: size %d, type %d, num_threads %d@\n@\n"
         n (int_of_mattype mattype) nthreads;

  (* Create matrices and vectors *)
  let matb = Matrix.dense n
  and x = Nvector_serial.make n 0.0
  and y = Nvector_serial.make n 0.0
  and b = Nvector_serial.make n 0.0
  in

  (* Fill matrix with uniform random data in [0,1/N] *)
  let matdata = Matrix.(Dense.unwrap (unwrap matb)) in
  for _k = 0 to (5*n - 1) do
    let i = Sundials.Util.rand () mod n in
    let j = Sundials.Util.rand () mod n in
    matdata.{j, i} <- f_rand () /. f_rand_max /. float_of_int n
  done;

  (* Add identity to matrix *)
  (try
     Matrix.scale_addi 1.0 matb
   with _ -> begin
     printf "FAIL: SUNLinSol SUNMatScaleAddI failure@\n";
     exit 1;
   end);

  (* Fill x vector with uniform random data in [0,1] *)
  let xdata = Nvector.unwrap x in
  Sundials.RealArray.map (fun _ -> f_rand() /. f_rand_max) xdata;

  (* Create sparse matrix from dense, and destroy B *)
  let mata = Matrix.(wrap_sparse (Sparse.from_dense mattype 0.0 (unwrap matb))) in

  (* copy x into y to print in case of solver failure *)
  Nvector_serial.Ops.scale 1.0 x y;

  (* create right-hand side vector for linear solve *)
  (try
     Matrix.matvec mata x b
   with _ -> begin
     printf "FAIL: SUNLinSol SUNMatMatvec failure@\n";
     exit 1
   end);

  (* Create dense linear solver *)
  let ls = LS.Direct.superlumt ~nthreads x mata in

  (* Run Tests *)
  fails += Test.sunlinsolinitialize ls 0;
  fails += Test.sunlinsolsetup ls (Some mata) 0;
  fails += Test.sunlinsolsolve ls (Some mata) x b (100.0 *. unit_roundoff) true 0;

  fails += Test.sunlinsolgettype ls LS.Direct 0;
  fails += Test.sunlinsolgetid ls LS.Superlumt 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolspace ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol module failed %d tests @\n @\n" !fails;
    printf "@\nA =@\n%a@\n" Matrix.Sparse.pp (Matrix.unwrap mata);
    printf "x (original) =@\n%a@\n" Nvector_serial.pp y;
    printf "b =@\n%a@\n" Nvector_serial.pp b;
    printf "x (computed) =@\n%a@" Nvector_serial.pp x
  end else
    printf "SUCCESS: SUNLinSol module passed all tests @\n @\n";

  !fails

let _ = main ()

