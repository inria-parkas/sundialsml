(*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2018.
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNMatrix Dense module
 * implementation.
 * -----------------------------------------------------------------
 *)

module Matrix = Sundials.Matrix

let printf = Format.printf
let (+=) r x = r := !r + x

module Dense_tests =
struct
  type k = Matrix.standard
  type m = Matrix.Dense.t
  type nd = Nvector_serial.data
  type nk = Nvector_serial.kind

  type t = nk Matrix.dense
  type nvec = Nvector_serial.t

  let rewrap = Matrix.wrap_dense

  let check_matrix a b tol =
    let ca, cb = Matrix.(unwrap a, unwrap b) in
    let ((ma, na) as sa), ((mb, nb) as sb) = Matrix.Dense.(size ca, size cb) in
    let adata, bdata = Matrix.Dense.(unwrap ca, unwrap cb) in
    if sa <> sb then
      (printf ">>> ERROR: check_matrix: Different data array lengths @\n"; true)
    else begin
      let failure = ref 0 in
      for i = 0 to ma - 1 do
        for j = 0 to na - 1 do
          failure += Test_matrix.fneq adata.{j, i} bdata.{j, i} tol
        done
      done;
      !failure > 0
    end

  let check_matrix_entry a v tol =
    let ca = Matrix.unwrap a in
    let ma, na = Matrix.Dense.size ca in
    let adata = Matrix.Dense.unwrap ca in
    let failure = ref 0 in
    for i = 0 to ma - 1 do
      for j = 0 to na - 1 do
        failure += Test_matrix.fneq adata.{j, i} v tol
      done
    done;
    !failure > 0

  let is_square a =
    let ca = Matrix.unwrap a in
    let ma, na = Matrix.Dense.size ca in
    ma = na

  let check_vector x y tol =
    let xdata, ydata = Nvector.(unwrap x, unwrap y) in
    let nx, ny = Sundials.RealArray.(length xdata, length ydata) in
    if nx <> ny then
      (printf ">>> ERROR: check_vector: Different data array lengths @\n"; true)
    else begin
      let failure = ref 0 in
      for i = 0 to nx - 1 do
        failure += Test_matrix.fneq xdata.{i} ydata.{i} tol
      done;
      if !failure > 0 then begin
        printf "Check_vector failures:@\n";
        for i = 0 to nx - 1 do
          if Test_matrix.fneq xdata.{i} ydata.{i} tol <> 0 then
            printf "  xdata[%d] = %g != %g (err = %g)@\n" i
                   xdata.{i} ydata.{i} (abs_float (xdata.{i} -. ydata.{i}))
        done;
        true
      end
      else false
    end

end

module Test = Test_matrix.Test (Dense_tests) (Nvector_serial.Ops)

(* ----------------------------------------------------------------------
 * Main Matrix Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in

  (* check input and set vector length *)
  if Array.length Sys.argv < 4 then
    (printf "ERROR: THREE (3) Input required: matrix rows, matrix cols, print timing @\n";
     exit (-1));

  let matrows = int_of_string Sys.argv.(1) in
  if matrows <= 0 then
    (printf "ERROR: number of rows must be a positive integer @\n";
     exit (-1));

  let matcols = int_of_string Sys.argv.(2) in
  if matcols <= 0 then
    (printf "ERROR: number of cols must be a positive integer @\n";
     exit (-1));

  let print_timing = int_of_string Sys.argv.(3) in
  let _ = Test.set_timing (print_timing <> 0) in

  let square = (matrows = matcols) in
  printf "@\nDense matrix test: size %d by %d@\n@\n" matrows matcols;

  (* Create vectors and matrices *)
  let x = Nvector_serial.make matcols 0.0
  and y = Nvector_serial.make matrows 0.0
  and a = Matrix.dense ~m:matrows matcols
  and i = Matrix.dense matcols
  in

  let adata = Matrix.(Dense.unwrap (unwrap a)) in
  for j = 0 to matcols -1 do
    for i = 0 to matrows - 1 do
      adata.{j, i} <- float_of_int ((j + 1) * (i + j))
    done
  done;

  let idata = Matrix.(Dense.unwrap (unwrap i)) in
  if square then begin
    for i = 0 to matrows -1 do
      idata.{i, i} <- 1.0
    done
  end;

  let xdata = Nvector.unwrap x in
  for i = 0 to matcols - 1 do
    xdata.{i} <- 1.0 /. float_of_int (i + 1)
  done;

  let ydata = Nvector.unwrap y in
  for i = 0 to matrows - 1 do
    let m = float_of_int i in
    let n = m +. float_of_int (matcols - 1) in
    ydata.{i} <- 0.5 *. (n +. 1. -. m) *. (n +. m)
  done;

  (* SUNMatrix Tests *)
  fails += Test.test_sunmatgetid a Matrix.Dense 0;
  fails += Test.test_sunmatclone a 0;
  fails += Test.test_sunmatcopy a 0;
  fails += Test.test_sunmatzero a 0;
  fails += Test.test_sunmatscaleadd a i 0;
  if square then
    fails += Test.test_sunmatscaleaddi a i 0;
  fails += Test.test_sunmatmatvec a x y 0;
  fails += Test.test_sunmatspace a 0;

  (* Print result *)
  if !fails <> 0 then begin
    printf "FAIL: SUNMatrix module failed %d tests @\n @\n" !fails;
    printf "@\nA = %a@\n" Matrix.Dense.pp (Matrix.unwrap a);
    if square then
      printf "@\nI = %a@\n" Matrix.Dense.pp (Matrix.unwrap i);
    printf "@\nx = %a@\n" Nvector_serial.pp x;
    printf "@\ny = %a@\n" Nvector_serial.pp y
  end
  else
    printf "SUCCESS: SUNMatrix module passed all tests @\n @\n"

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

