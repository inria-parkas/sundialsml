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
 * This is the testing routine to check the SUNMatrix Band module
 * implementation.
 * -----------------------------------------------------------------
 *)

module Matrix = Sundials.Matrix

let printf = Format.printf
let (+=) r x = r := !r + x

module Band_tests =
struct
  type k = Matrix.standard
  type m = Matrix.Band.t
  type nd = Nvector_serial.data
  type nk = Nvector_serial.kind

  type t = nk Matrix.band
  type nvec = Nvector_serial.t

  let rewrap = Matrix.wrap_band

  let check_matrix a b tol =
    let ca, cb = Matrix.(unwrap a, unwrap b) in
    let open Matrix.Band in
    let { n = na; smu = smua; mu = mua; ml = mla } = dims ca in
    let { n = nb; smu = smub; mu = mub; ml = mlb } = dims cb in
    let adata, bdata = unwrap ca, unwrap cb in
    let d1a, d2a = Bigarray.Array2.(dim1 adata, dim2 adata) in
    let d1b, d2b = Bigarray.Array2.(dim1 adata, dim2 adata) in
    if na <> nb || mua <> mub || mla <> mlb || d1a <> d1b || d2a <> d2b
    then true
    else begin
      let failure = ref 0 in
      for j = 0 to na - 1 do
        (* compare entries in this column *)
        let istart = if j < mua then -j else -mua in
        let iend = if j > na - 1 - mla then na - 1 - j else mla in
        for i = istart to iend do
          failure += Test_matrix.fneq adata.{j, i + smua} bdata.{j, i + smub} tol
        done
      done;

      if !failure > 0 then begin
        printf "check_matrix failure, A = @\n%a@\n" pp ca;
        printf "B = @\n%a@\n" pp cb;
        true
      end else false
    end

  let check_matrix_entry a v tol =
    let ca = Matrix.unwrap a in
    let open Matrix.Band in
    let { n; smu; mu; ml } = dims ca in
    let adata = unwrap ca in
    let failure = ref 0 in
    for j = 0 to n - 1 do
      let istart = if j < mu then -j else -mu in
      let iend = if j > n - 1 - ml then n - 1 - j else ml in
      for i = istart to iend do
        if Test_matrix.fneq adata.{j, i + smu} v tol > 0 then begin
          incr failure;
          printf "j = %d, Acolj[%d] = %g, val = %g\n" j i adata.{j, i} v
        end
      done;
    done;
    !failure > 0

  let is_square a = true

  let check_vector x y tol =
    let xdata, ydata = Nvector.(unwrap x, unwrap y) in
    let nx = Sundials.RealArray.(length xdata) in
    let failure = ref 0 in
    for i = 0 to nx - 1 do
      failure += Test_matrix.fneq xdata.{i} ydata.{i} tol
    done;
    !failure > 0

  let nvec_pp = Nvector_serial.pp
end

module Test = Test_matrix.Test (Band_tests) (Nvector_serial.Ops)

(* ----------------------------------------------------------------------
 * Main Matrix Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in

  (* check input and set vector length *)
  if Array.length Sys.argv < 5 then
    (printf "ERROR: FOUR (4) Inputs required: matrix cols, matrix uband, matrix lband, print timing @\n";
     exit (-1));

  let cols = int_of_string Sys.argv.(1) in
  if cols <= 0 then
    (printf "ERROR: number of matrix columns must be a positive integer @\n";
     exit (-1));

  let uband = int_of_string Sys.argv.(2) in
  if uband <= 0 || uband >= cols then
    (printf "ERROR: matrix upper bandwidth must be a positive integer, less than number of columns @\n";
     exit (-1));

  let lband = int_of_string Sys.argv.(3) in
  if lband <= 0 || lband >= cols then
    (printf "ERROR: matrix lower bandwidth must be a positive integer, less than number of columns @\n";
     exit (-1));

  let print_timing = int_of_string Sys.argv.(4) in
  let _ = Test.set_timing (print_timing <> 0) in

  printf "@\nBand matrix test: size %d, bandwidths %d %d@\n@\n"
    cols uband lband;

  (* Create vectors and matrices *)
  let x = Nvector_serial.make cols 0.0
  and y = Nvector_serial.make cols 0.0
  and a = Matrix.band ~mu:uband ~ml:lband cols
  and i = Matrix.band ~mu:0 ~ml:0 cols
  in

  (* Fill matrices *)
  let adata = Matrix.(Band.unwrap (unwrap a))
  and idata = Matrix.(Band.unwrap (unwrap i))
  and xdata = Nvector.unwrap x
  and ydata = Nvector.unwrap y
  in
  for j = 0 to cols -1 do
    (* identity matrix *)
    idata.{j, 0} <- 1.0;

    (* A matrix *)
    let kstart = if j < uband then -j else -uband in
    let kend = if j > cols - 1 - lband then cols -1 - j else lband in
    for k = kstart to kend do
      adata.{j, k + (uband+lband)} <- float_of_int (j - k)
    done
  done;

  (* Fill vectors *)
  for i = 0 to cols - 1 do
    xdata.{i} <- float_of_int i;

    (* y vector *)
    ydata.{i} <- 0.0;
    let jstart = max 0 (i - lband) in
    let jend = min (cols - 1) (i + uband) in
    for j = jstart to jend do
      ydata.{i} <- ydata.{i} +. float_of_int ((j + j - i) * j)
    done
  done;

  (* SUNMatrix Tests *)
  fails += Test.test_sunmatgetid a Matrix.Band 0;
  fails += Test.test_sunmatclone a 0;
  fails += Test.test_sunmatcopy a 0;
  fails += Test.test_sunmatzero a 0;
  fails += Test.test_sunmatscaleadd a i 0;
  fails += Test.test_sunmatscaleaddi a i 0;
  fails += Test.test_sunmatmatvec a x y 0;
  fails += Test.test_sunmatspace a 0;

  (* Print result *)
  if !fails <> 0 then begin
    printf "FAIL: SUNMatrix module failed %d tests @\n @\n" !fails;
    printf "@\nA = %a@\n" Matrix.Band.pp (Matrix.unwrap a);
    printf "@\nI = %a@\n" Matrix.Band.pp (Matrix.unwrap i);
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

