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
 * This is the testing routine to check the SUNLinSol LapackDense
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
module LapackDense_tests : (Test_linsol.LINSOL_TESTS
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

module Test = Test_linsol.Test (LapackDense_tests) (Nvector_serial.Ops)

let sync_device () = ()

(* ----------------------------------------------------------------------
 * SUNLinSol_Dense Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in

  (* check input and set matrix dimensions *)
  if Array.length Sys.argv < 3 then
    (printf "ERROR: THREE (3) Inputs required: matrix cols, print matrix on fail, print timing \n";
     exit (-1));

  let cols = int_of_string Sys.argv.(1) in
  if cols <= 0 then
    (printf "ERROR: number of matrix columns must be a positive integer @\n";
     exit (-1));

  let rows = cols in

  let print_matrix_on_fail = int_of_string Sys.argv.(2) <> 0 in

  let print_timing = int_of_string Sys.argv.(3) in
  let _ = Test.set_timing (print_timing <> 0) in

  printf "@\nLapackDense linear solver test: size %d@\n@\n" cols;

  (* Create matrices and vectors *)
  let mata = Matrix.dense ~m:rows cols
  and matb = Matrix.dense ~m:rows cols
  and mati = Matrix.dense ~m:rows cols
  and x = Nvector_serial.make cols 0.0
  and y = Nvector_serial.make cols 0.0
  and b = Nvector_serial.make cols 0.0
  in

  (* Fill A matrix with uniform random data in [0,1/cols] *)
  let adata = Matrix.(Dense.unwrap (unwrap mata)) in
  for j = 0 to cols - 1 do
    for k = 0 to rows - 1 do
      adata.{j, k} <- f_rand () /. f_rand_max /. float_of_int cols
    done
  done;

  (* Create anti-identity matrix *)
  let idata = Matrix.(Dense.unwrap (unwrap mati)) in
  for k = 0 to rows - 1 do
    idata.{cols - 1 - k, k} <- 1.
  done;

  (* Add anti-identity to ensure the solver needs to do row-swapping *)
  for k = 0 to rows - 1 do
    for j = 0 to cols - 1 do
      adata.{j, k} <- adata.{j, k} +. idata.{j, k}
    done
  done;

  (* Fill x vector with uniform random data in [0,1] *)
  let xdata = Nvector.unwrap x in
  for j = 0 to cols - 1 do
    xdata.{j} <- f_rand () /. f_rand_max
  done;

  (* copy A and x into B and y to print in case of solver failure *)
  Matrix.blit ~src:mata ~dst:matb;
  Nvector_serial.Ops.scale 1.0 x y;

  (* create right-hand side vector for linear solve *)
  (try
     Matrix.matvec mata x b
   with _ -> begin
     printf "FAIL: SUNLinSol SUNMatMatvec failure@\n";
     exit 1
   end);

  (* Create dense linear solver *)
  let ls = LS.Direct.lapack_dense x mata in

  (* Run Tests *)
  fails += Test.sunlinsolinitialize ls 0;
  fails += Test.sunlinsolsetup ls (Some mata) 0;
  fails += Test.sunlinsolsolve ls (Some mata) x b (100.0 *. unit_roundoff) true 0;

  fails += Test.sunlinsolgettype ls LS.Direct 0;
  fails += Test.sunlinsolgetid ls LS.LapackDense 0;
  fails += Test.sunlinsollastflag ls 0;
  fails += Test.sunlinsolspace ls 0;

  (* Print result *)
  if !fails > 0 then begin
    printf "FAIL: SUNLinSol module failed %d tests @\n @\n" !fails;
    printf "@\nanswer =@\n%a@\n" Nvector_serial.pp y;
    printf "computed =@\n%a@\n" Nvector_serial.pp x;
    Nvector_serial.Ops.linearsum 1.0 y (-1.0) x x;
    printf "diff (answer-computed) =@\n%a" Nvector_serial.pp x;
    if print_matrix_on_fail then begin
      printf "@\nA (original) =@\n%a@\n" Matrix.Dense.pp (Matrix.unwrap matb);
      printf "A (factored) =\n%a" Matrix.Dense.pp (Matrix.unwrap mata)
    end
  end else
    printf "SUCCESS: SUNLinSol module passed all tests @\n @\n";

  !fails

let _ = main ()

