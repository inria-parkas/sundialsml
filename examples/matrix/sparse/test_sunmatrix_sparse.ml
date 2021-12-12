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
 * This is the testing routine to check the SUNMatrix Sparse module
 * implementation.
 * -----------------------------------------------------------------
 *)

module Matrix = Sundials.Matrix

let printf = Format.printf
let (+=) r x = r := !r + x

let rand_float () =
  (float_of_int (Test_matrix.rand ()) /. float_of_int (Test_matrix.rand_max))

module Sparse_tests (S : sig type sformat end) =
struct
  type k = Matrix.standard
  type m = S.sformat Matrix.Sparse.t
  type nd = Nvector_serial.data
  type nk = Nvector_serial.kind

  type t = (S.sformat, nk) Matrix.sparse
  type nvec = Nvector_serial.t

  let rewrap = Matrix.wrap_sparse

  let check_matrix a b tol =
    let ac, bc = Matrix.(unwrap a, unwrap b) in
    let (am, an), (bm, bn) = Matrix.Sparse.(size ac, size bc) in
    let (annz, anp), (bnnz, _) = Matrix.Sparse.(dims ac, dims bc) in
    try
      if am <> bm then
        (printf ">>> ERROR: check_matrix: Different numbers of rows (%d vs %d)@\n"
           am bm; raise Exit);
      if an <> bn then
        (printf ">>> ERROR: check_matrix: Different numbers of columns (%d vs %d)@\n"
           an bn; raise Exit);
      if annz <> bnnz then
        (printf ">>> ERROR: check_matrix: Different numbers of nonzeros (%d vs %d)@\n"
           annz bnnz; raise Exit);
      let aindexvals, aindexptrs, adata = Matrix.Sparse.unwrap ac in
      let bindexvals, bindexptrs, bdata = Matrix.Sparse.unwrap bc in

      (* compare sparsity patterns *)
      for i = 0 to anp-1 do
        if aindexptrs.{i} <> bindexptrs.{i} then
          (printf ">>> ERROR: check_matrix: Different indexptrs @\n";
           raise Exit)
      done;
      for i = 0 to annz-1 do
        if aindexvals.{i} <> bindexvals.{i} then
          (printf ">>> ERROR: check_matrix: Different indexvals @\n";
           raise Exit)
      done;
      (* compare matrix values *)
      for i = 0 to annz - 1 do
        if Test_matrix.fneq adata.{i} bdata.{i} tol > 1 then
          (printf ">>> ERROR: check_matrix: Different entries @\n";
           raise Exit)
      done;
      false
    with Exit -> true

  let check_matrix_entry a v tol =
    let c = Matrix.unwrap a in
    let _, _, data = Matrix.Sparse.unwrap c in
    let nnz, _ = Matrix.Sparse.dims c in
    try
      for i = 0 to nnz - 1 do
        if Test_matrix.fneq data.{i} v tol > 0 then raise Exit
      done;
      false
    with Exit -> true

  let is_square a =
    let ca = Matrix.unwrap a in
    let ma, na = Matrix.Sparse.size ca in
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

  let nvec_pp = Nvector_serial.pp
end

(* ----------------------------------------------------------------------
 * Main Matrix Testing Routine
 * --------------------------------------------------------------------*)

let clone a = Matrix.(wrap_sparse ((get_ops a).m_clone (unwrap a)))

(* ----------------------------------------------------------------------
 * Extra ScaleAdd tests for sparse matrices:
 *    A and B should have different sparsity patterns, and neither should
 *      contain sufficient storage to for their sum
 *    y should already equal A*x
 *    z should already equal B*x
 * --------------------------------------------------------------------*)
let test_sunmatscaleadd2 check_vector a b x y z =
  let tol = match Sundials.Config.sundials_version with
            | 2,_,_ | 3,1,0 | 3,1,1 | 3,1,2 -> 1e-14
            | _ -> 100.0 *. Sundials.Config.unit_roundoff
  in
  try
    (* create clones for test *)
    let c = clone a in
    let u = Nvector_serial.Ops.clone y in
    let v = Nvector_serial.Ops.clone y in

    (* test 1: add A to B (output must be enlarged) *)
    (try Matrix.blit ~src:a ~dst:c
     with _ -> printf ">>> FAILED test -- SUNMatCopy returned 1 \n";
               raise Exit);
    (try Matrix.scale_add 1.0 c b
     with _ -> printf ">>> FAILED test -- SUNMatScaleAdd returned 1 \n";
               raise Exit);
    (try Matrix.matvec c x u
     with _ -> printf ">>> FAILED test -- SUNMatMatvec returned 1 \n";
               raise Exit);
    Nvector_serial.Ops.linearsum 1. y 1. z v;  (* v = y+z *)
    if not (check_vector u v tol)                 (* u ?= v *)
    then printf "    PASSED test -- SUNMatScaleAdd2 check 1 \n"
    else begin
      printf ">>> FAILED test -- SUNMatScaleAdd2 check 1 \n";
      printf "\nA =@.";
      Matrix.print_sparse a Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nB =@.";
      Matrix.print_sparse b Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nC =@.";
      Matrix.print_sparse c Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nx =\n%a@\n" Nvector_serial.pp x;
      printf "\ny =\n%a@\n" Nvector_serial.pp y;
      printf "\nz =\n%a@\n" Nvector_serial.pp z;
      printf "\nu =\n%a@\n" Nvector_serial.pp u;
      printf "\nv =\n%a@\n" Nvector_serial.pp v;
      raise Exit
    end;

    (* test 2: add A to a matrix with sufficient but misplaced storage *)
    let d = clone a in
    (try let a_nnz, _ = Matrix.(Sparse.dims (unwrap a)) in
         let b_nnz, _ = Matrix.(Sparse.dims (unwrap b)) in
         Matrix.(Sparse.resize ~nnz:(a_nnz+b_nnz) (unwrap d));
         Matrix.blit ~src:a ~dst:d (* D = A *)
     with _ -> printf ">>> FAILED test -- SUNMatCopy returned 1 @\n";
               raise Exit);
    (try Matrix.scale_add 1.0 d b (* D = A+B *)
     with _ -> printf ">>> FAILED test -- SUNMatScaleAdd returned 1 @\n";
               raise Exit);
    (try Matrix.matvec d x u (* u = Cx = Ax+Bx *)
     with _ -> printf ">>> FAILED test -- SUNMatMatvec returned 1 @\n";
               raise Exit);
    Nvector_serial.Ops.linearsum 1. y 1. z v; (* v = y+z *)
    if not (check_vector u v tol)                (* u ?= v *)
    then printf "    PASSED test -- SUNMatScaleAdd2 check 2 @\n"
    else begin
      printf ">>> FAILED test -- SUNMatScaleAdd2 check 2 @\n";
      printf "\nA =@.";
      Matrix.print_sparse a Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nB =@.";
      Matrix.print_sparse b Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nC =@.";
      Matrix.print_sparse c Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nx =\n%a@\n" Nvector_serial.pp x;
      printf "\ny =\n%a@\n" Nvector_serial.pp y;
      printf "\nz =\n%a@\n" Nvector_serial.pp z;
      printf "\nu =\n%a@\n" Nvector_serial.pp u;
      printf "\nv =\n%a@\n" Nvector_serial.pp v;
      raise Exit
    end;

    (* test 3: add A to a matrix with the appropriate structure already in place *)
    let e = clone c in
    (try Matrix.blit ~src:c ~dst:e  (* E = A + B *)
     with _ -> printf ">>> FAILED test -- SUNMatCopy returned 1 @\n";
               raise Exit);

    (try Matrix.scale_add (-1.0) e b (* E = -A *)
     with _ -> printf ">>> FAILED test -- SUNMatScaleAdd returned 1 @\n";
               raise Exit);
    (try Matrix.matvec e x u (* u = Ex = -Ax *)
     with _ -> printf ">>> FAILED test -- SUNMatMatvec returned 1 @\n";
               raise Exit);
    Nvector_serial.Ops.linearsum (-1.0) y 0.0 z v; (* v = -y *)
    if not (check_vector u v tol)                     (* v ?= u *)
    then printf "    PASSED test -- SUNMatScaleAdd2 check 3 @\n"
    else begin
      printf ">>> FAILED test -- SUNMatScaleAdd2 check 3 @\n";
      printf "\nA =@.";
      Matrix.print_sparse a Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nB =@.";
      Matrix.print_sparse b Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nC =@.";
      Matrix.print_sparse c Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nE =@.";
      Matrix.print_sparse e Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nx =\n%a@\n" Nvector_serial.pp x;
      printf "\ny =\n%a@\n" Nvector_serial.pp y;
      printf "\nu =\n%a@\n" Nvector_serial.pp u;
      printf "\nv =\n%a@\n" Nvector_serial.pp v;
      raise Exit
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * Extra ScaleAddI tests for sparse matrices:
 *    A should not contain values on the diagonal, nor should it contain
 *      sufficient storage to add those in
 *    y should already equal A*x
 * --------------------------------------------------------------------*)
let test_sunmatscaleaddi2 check_vector a x y =
  let tol = match Sundials.Config.sundials_version with
            | 2,_,_ | 3,1,0 | 3,1,1 | 3,1,2 -> 1e-14
            | _ -> 100.0 *. Sundials.Config.unit_roundoff
  in
  try
    (* create clones for test *)
    let b = clone a in
    let z = Nvector_serial.Ops.clone x in
    let w = Nvector_serial.Ops.clone x in

    (* test 1: add I to a matrix with insufficient storage *)
    (try Matrix.blit ~src:a ~dst:b
     with _ -> printf ">>> FAILED test -- SUNMatCopy returned 1 @\n";
               raise Exit);
    (try Matrix.scale_addi (-1.0) b (* B = I-A *)
     with _ -> printf ">>> FAILED test -- SUNMatScaleAddI returned 1 @\n";
               raise Exit);
    (try Matrix.matvec b x z
     with _ -> printf ">>> FAILED test -- SUNMatMatvec returned 1 @\n";
               raise Exit);
    Nvector_serial.Ops.linearsum 1. x (-1.0) y w;
    if not (check_vector z w tol)
    then printf "    PASSED test -- SUNMatScaleAddI2 check 1 @\n"
    else begin
      printf ">>> FAILED test -- SUNMatScaleAddI2 check 1 @\n";
      printf "\nA =@.";
      Matrix.print_sparse a Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nB =@.";
      Matrix.print_sparse b Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nz =\n%a@\n" Nvector_serial.pp z;
      printf "\nw =\n%a@\n" Nvector_serial.pp w;
      raise Exit
    end;

    (* test 2: add I to a matrix with sufficient but misplaced storage *)
    let c = clone a in
    (try let a_m       = Matrix.unwrap a in
         let a_nnz, _  = Matrix.Sparse.dims a_m in
         let a_rows, _ = Matrix.Sparse.size a_m in
         Matrix.(Sparse.resize ~nnz:(a_nnz+a_rows) (unwrap c));
         Matrix.blit ~src:a ~dst:c
     with _ -> printf ">>> FAILED test -- SUNMatCopy returned 1 @\n";
               raise Exit);
    (try Matrix.scale_addi (-1.0) c (* C = I-A *)
     with _ -> printf ">>> FAILED test -- SUNMatScaleAddI returned 1 @\n";
               raise Exit);
    (try Matrix.matvec c x z (* u = Cx = Ax+Bx *)
     with _ -> printf ">>> FAILED test -- SUNMatMatvec returned 1 @\n";
               raise Exit);
    Nvector_serial.Ops.linearsum 1. x (-1.0) y w;
    if not (check_vector z w tol)                 (* u ?= v *)
    then printf "    PASSED test -- SUNMatScaleAddI2 check 2 @\n"
    else begin
      printf ">>> FAILED test -- SUNMatScaleAddI2 check 2 @\n";
      printf "\nA =@.";
      Matrix.print_sparse a Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nC =@.";
      Matrix.print_sparse c Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nz =\n%a@\n" Nvector_serial.pp z;
      printf "\nw =\n%a@\n" Nvector_serial.pp w;
      raise Exit
    end;

    (* test 3: add I to a matrix with appropriate structure already in place *)
    let d = clone c in
    (try Matrix.blit ~src:c ~dst:d
     with _ -> printf ">>> FAILED test -- SUNMatCopy returned 1 @\n";
               raise Exit);
    (try Matrix.scale_addi (-1.0) d (* D = A *)
     with _ -> printf ">>> FAILED test -- SUNMatScaleAddI returned 1 @\n";
               raise Exit);
    (try Matrix.matvec d x z
     with _ -> printf ">>> FAILED test -- SUNMatMatvec returned 1 @\n";
               raise Exit);
    if not (check_vector z y tol)
    then printf "    PASSED test -- SUNMatScaleAddI2 check 3 @\n"
    else begin
      printf ">>> FAILED test -- SUNMatScaleAddI2 check 3 @\n";
      printf "\nA =@.";
      Matrix.print_sparse a Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nD =@.";
      Matrix.print_sparse d Sundials.Logfile.stdout;
      Sundials.Logfile.flush Sundials.Logfile.stdout;
      printf "\nz =\n%a@\n" Nvector_serial.pp z;
      printf "\ny =\n%a@\n" Nvector_serial.pp y;
      raise Exit
    end;
    0
  with Exit -> 1

let test_sunsparsematrix_convert (type s)
    (check_matrix : (s, 'a) Matrix.sparse -> (s, 'a) Matrix.sparse -> float -> bool)
    (am : (s, 'a) Matrix.sparse) =
  let tol = 200.0 *. Sundials.Config.unit_roundoff in
  let a = Matrix.unwrap am in
  match Matrix.Sparse.((sformat a : s sformat)) with
  | Matrix.Sparse.CSC ->
      let csr = Matrix.Sparse.copy_to_csr a in
      let csc = Matrix.(wrap_sparse (Sparse.copy_to_csc csr)) in
      if check_matrix am csc tol then
        (printf ">>> FAILED test --  Test_SUNSparseMatrixToCSR check_matrix failed\n"; 1)
      else (printf "    PASSED test -- SUNSparseMatrixToCSR\n"; 0)
  | Matrix.Sparse.CSR ->
      let csc = Matrix.Sparse.copy_to_csc a in
      let csr = Matrix.(wrap_sparse (Sparse.copy_to_csr csc)) in
      if check_matrix am csr tol then
        (printf ">>> FAILED test --  Test_SUNSparseMatrixToCSC check_martrix failed\n"; 1)
      else (printf "    PASSED test -- SUNSparseMatrixToCSC\n"; 0)

let int_of_mattype (type s) : s Matrix.Sparse.sformat -> int = function
  | Matrix.Sparse.CSC -> 0
  | Matrix.Sparse.CSR -> 1

let to_index = Sundials.Index.of_int

let rec main () =
  (* check input and set vector length *)
  if Array.length Sys.argv < 5 then
    (printf "ERROR: FOUR (4) Input required: matrix rows, matrix cols, matrix type (0/1), print timing @\n";
     exit (-1));

  let matrows = int_of_string Sys.argv.(1) in
  if matrows <= 0 then
    (printf "ERROR: number of rows must be a positive integer@\n";
     exit (-1));

  let matcols = int_of_string Sys.argv.(2) in
  if matcols <= 0 then
    (printf "ERROR: number of cols must be a positive integer@\n";
     exit (-1));

  let k = int_of_string Sys.argv.(3) in
  if k <> 0 && k <> 1 then
    (printf "ERROR: matrix type must be 0 or 1@\n";
     exit (-1));

  if k = 0
  then main_with_type matrows matcols Matrix.Sparse.CSC
  else main_with_type matrows matcols Matrix.Sparse.CSR

and main_with_type : type s. int -> int -> s Matrix.Sparse.sformat -> unit
  = fun matrows matcols mattype ->

  let module SparseTests =(Sparse_tests (struct type sformat = s end)) in
  let module Test = Test_matrix.Test (SparseTests) (Nvector_serial.Ops) in

  let fails = ref 0 in

  let print_timing = int_of_string Sys.argv.(4) in
  let _ = Test.set_timing (print_timing <> 0) in

  let square = (matrows = matcols) in
  printf "@\nSparse matrix test: size %d by %d, type = %d@\n@."
    matrows matcols (int_of_mattype mattype);

  (* check creating sparse matrix from dense matrix *)
  let b = Matrix.dense ~m:5 6 in
  let matdata = Matrix.(Dense.unwrap (unwrap b)) in

  matdata.{0, 2} <- 1.0;    (* [ 0 2 0 0 7 0 ] *)
  matdata.{1, 0} <- 2.0;    (* [ 0 0 4 0 8 0 ] *)
  matdata.{1, 4} <- 3.0;    (* [ 1 0 0 0 0 0 ] *)
  matdata.{2, 1} <- 4.0;    (* [ 0 0 5 6 0 0 ] *)
  matdata.{2, 3} <- 5.0;    (* [ 0 3 0 0 0 9 ] *)
  matdata.{3, 3} <- 6.0;
  matdata.{4, 0} <- 7.0;
  matdata.{4, 1} <- 8.0;
  matdata.{5, 4} <- 9.0;

  (match mattype with
  | Matrix.Sparse.CSR ->
      (* Check CSR *)
      let c = Matrix.Sparse.make mattype 5 6 9 in
      let colindices, rowptrs, matdata = Matrix.Sparse.unwrap c in
      rowptrs.{0} <- to_index 0;
      matdata.{0} <- 2.0;           colindices.{0} <- to_index 1;
      matdata.{1} <- 7.0;           colindices.{1} <- to_index 4;
      rowptrs.{1} <- to_index 2;
      matdata.{2} <- 4.0;           colindices.{2} <- to_index 2;
      matdata.{3} <- 8.0;           colindices.{3} <- to_index 4;
      rowptrs.{2} <- to_index 4;
      matdata.{4} <- 1.0;           colindices.{4} <- to_index 0;
      rowptrs.{3} <- to_index 5;
      matdata.{5} <- 5.0;           colindices.{5} <- to_index 2;
      matdata.{6} <- 6.0;           colindices.{6} <- to_index 3;
      rowptrs.{4} <- to_index 7;
      matdata.{7} <- 3.0;           colindices.{7} <- to_index 1;
      matdata.{8} <- 9.0;           colindices.{8} <- to_index 5;
      rowptrs.{5} <- to_index 9;

      let ca = Matrix.Sparse.from_dense mattype 0.0 (Matrix.unwrap b) in
      let a = Matrix.wrap_sparse ca in
      if SparseTests.check_matrix a (Matrix.wrap_sparse c) 1e-15 then
        (printf "FAIL: SUNMatrix SparseFromDense CSR conversion failed\n";
         exit 1)

  | Matrix.Sparse.CSC ->
      (* Check CSC *)
      let d = Matrix.Sparse.make mattype 5 6 9 in
      let rowindices, colptrs, matdata = Matrix.Sparse.unwrap d in
      colptrs.{0} <- to_index 0;
      matdata.{0} <- 1.0;           rowindices.{0} <- to_index 2;
      colptrs.{1} <- to_index 1;
      matdata.{1} <- 2.0;           rowindices.{1} <- to_index 0;
      matdata.{2} <- 3.0;           rowindices.{2} <- to_index 4;
      colptrs.{2} <- to_index 3;
      matdata.{3} <- 4.0;           rowindices.{3} <- to_index 1;
      matdata.{4} <- 5.0;           rowindices.{4} <- to_index 3;
      colptrs.{3} <- to_index 5;
      matdata.{5} <- 6.0;           rowindices.{5} <- to_index 3;
      colptrs.{4} <- to_index 6;
      matdata.{6} <- 7.0;           rowindices.{6} <- to_index 0;
      matdata.{7} <- 8.0;           rowindices.{7} <- to_index 1;
      colptrs.{5} <- to_index 8;
      matdata.{8} <- 9.0;           rowindices.{8} <- to_index 4;
      colptrs.{6} <- to_index 9;

      let ca = Matrix.Sparse.from_dense mattype 0.0 (Matrix.unwrap b) in
      let a = Matrix.wrap_sparse ca in
      if SparseTests.check_matrix a (Matrix.wrap_sparse d) 1e-15 then
        (printf "FAIL: SUNMatrix SparseFromDense CSC conversion failed\n";
         exit 1)
  );

  (* check creating sparse matrix from banded matrix *)
  let n = 7 in
  let uband = 1 in
  let lband = 2 in                                     (* B(i,j) = j + (j-i) *)
  let suband = 3 in
  let b = Matrix.band ~mu:uband ~smu:suband ~ml:lband n in
  let matdata = Matrix.(Band.unwrap (unwrap b))    (* B = [  0  2  0  0  0  0  0 ] *)
  in                                                   (* [ -1  1  3  0  0  0  0 ] *)
  for j = 0 to n - 1 do                                (* [ -2  0  2  4  0  0  0 ] *)
    let kstart = if j<uband then -j else -uband in     (* [  0 -1  1  3  5  0  0 ] *)
    let kend = if j>n-1-lband then n-1-j else lband in (* [  0  0  0  2  4  6  0 ] *)
    for k = kstart to kend do                          (* [  0  0  0  1  3  5  7 ] *)
      matdata.{j, k + suband} <- float_of_int (j - k); (* [  0  0  0  0  2  4  6 ] *)
    done
  done;

  (match mattype with
  | Matrix.Sparse.CSR ->
      (* CSR *)
      let c = Matrix.Sparse.make mattype 7 7 21 in
      let colindices, rowptrs, matdata = Matrix.Sparse.unwrap c in
      rowptrs.{ 0} <- to_index 0;
      matdata.{ 0} <- 2.0;          colindices.{ 0} <- to_index 1;
      rowptrs.{ 1} <- to_index 1;
      matdata.{ 1} <- -1.0;         colindices.{ 1} <- to_index 0;
      matdata.{ 2} <- 1.0;          colindices.{ 2} <- to_index 1;
      matdata.{ 3} <- 3.0;          colindices.{ 3} <- to_index 2;
      rowptrs.{ 2} <- to_index 4;
      matdata.{ 4} <- -2.0;         colindices.{ 4} <- to_index 0;
      matdata.{ 5} <- 2.0;          colindices.{ 5} <- to_index 2;
      matdata.{ 6} <- 4.0;          colindices.{ 6} <- to_index 3;
      rowptrs.{ 3} <- to_index 7;
      matdata.{ 7} <- -1.0;         colindices.{ 7} <- to_index 1;
      matdata.{ 8} <- 1.0;          colindices.{ 8} <- to_index 2;
      matdata.{ 9} <- 3.0;          colindices.{ 9} <- to_index 3;
      matdata.{10} <- 5.0;          colindices.{10} <- to_index 4;
      rowptrs.{ 4} <- to_index 11;
      matdata.{11} <- 2.0;          colindices.{11} <- to_index 3;
      matdata.{12} <- 4.0;          colindices.{12} <- to_index 4;
      matdata.{13} <- 6.0;          colindices.{13} <- to_index 5;
      rowptrs.{ 5} <- to_index 14;
      matdata.{14} <- 1.0;          colindices.{14} <- to_index 3;
      matdata.{15} <- 3.0;          colindices.{15} <- to_index 4;
      matdata.{16} <- 5.0;          colindices.{16} <- to_index 5;
      matdata.{17} <- 7.0;          colindices.{17} <- to_index 6;
      rowptrs.{ 6} <- to_index 18;
      matdata.{18} <- 2.0;          colindices.{18} <- to_index 4;
      matdata.{19} <- 4.0;          colindices.{19} <- to_index 5;
      matdata.{20} <- 6.0;          colindices.{20} <- to_index 6;
      rowptrs.{ 7} <- to_index 21;

      let ca = Matrix.Sparse.from_band mattype 0.0 (Matrix.unwrap b) in
      let a = Matrix.wrap_sparse ca in
      if SparseTests.check_matrix a (Matrix.wrap_sparse c) 1e-15 then
        (printf("FAIL: SUNMatrix SparseFromBand CSR conversion failed\n");
         exit 1)

  | Matrix.Sparse.CSC ->
      (* Check CSC *)
      let d = Matrix.Sparse.make mattype 7 7 21 in
      let rowindices, colptrs, matdata = Matrix.Sparse.unwrap d in
      colptrs.{ 0} <- to_index 0;
      matdata.{ 0} <- -1.0;         rowindices.{ 0} <- to_index 1;
      matdata.{ 1} <- -2.0;         rowindices.{ 1} <- to_index 2;
      colptrs.{ 1} <- to_index 2;
      matdata.{ 2} <- 2.0;          rowindices.{ 2} <- to_index 0;
      matdata.{ 3} <- 1.0;          rowindices.{ 3} <- to_index 1;
      matdata.{ 4} <- -1.0;         rowindices.{ 4} <- to_index 3;
      colptrs.{ 2} <- to_index 5;
      matdata.{ 5} <- 3.0;          rowindices.{ 5} <- to_index 1;
      matdata.{ 6} <- 2.0;          rowindices.{ 6} <- to_index 2;
      matdata.{ 7} <- 1.0;          rowindices.{ 7} <- to_index 3;
      colptrs.{ 3} <- to_index 8;
      matdata.{ 8} <- 4.0;          rowindices.{ 8} <- to_index 2;
      matdata.{ 9} <- 3.0;          rowindices.{ 9} <- to_index 3;
      matdata.{10} <- 2.0;          rowindices.{10} <- to_index 4;
      matdata.{11} <- 1.0;          rowindices.{11} <- to_index 5;
      colptrs.{ 4} <- to_index 12;
      matdata.{12} <- 5.0;          rowindices.{12} <- to_index 3;
      matdata.{13} <- 4.0;          rowindices.{13} <- to_index 4;
      matdata.{14} <- 3.0;          rowindices.{14} <- to_index 5;
      matdata.{15} <- 2.0;          rowindices.{15} <- to_index 6;
      colptrs.{ 5} <- to_index 16;
      matdata.{16} <- 6.0;          rowindices.{16} <- to_index 4;
      matdata.{17} <- 5.0;          rowindices.{17} <- to_index 5;
      matdata.{18} <- 4.0;          rowindices.{18} <- to_index 6;
      colptrs.{ 6} <- to_index 19;
      matdata.{19} <- 7.0;          rowindices.{19} <- to_index 5;
      matdata.{20} <- 6.0;          rowindices.{20} <- to_index 6;
      colptrs.{ 7} <- to_index 21;

      let ca = Matrix.Sparse.from_band mattype 1e-15 (Matrix.unwrap b) in
      let a = Matrix.wrap_sparse ca in
      if SparseTests.check_matrix a (Matrix.wrap_sparse d) 1e-15 then
        (printf("FAIL: SUNMatrix SparseFromBand CSC conversion failed\n");
         exit 1)
  );

  (* Create/fill I matrix *)
  let ci = Matrix.Sparse.make mattype matrows matcols matcols in
  let i = Matrix.wrap_sparse ci in
  let colindices, rowptrs, matdata = Matrix.Sparse.unwrap ci in
  for i = 0 to matrows - 1 do
    matdata.{i}    <- 1.0;
    colindices.{i} <- to_index i;
    rowptrs.{i}    <- to_index i
  done;
  rowptrs.{matrows} <- to_index matrows;

  (* Create/fill random dense matrices, create sparse from them *)
  let mc = Matrix.Dense.make matrows matcols 0.0 in
  let c = Matrix.wrap_dense mc in
  let md = Matrix.Dense.make matrows matcols 0.0 in
  let d = Matrix.wrap_dense md in
  for _ = 0 to 3*matrows - 1 do
    let i = Test_matrix.rand () mod matrows in
    let j = Test_matrix.rand () mod matcols in
    Matrix.Dense.set md i j (rand_float ())
  done;
  for _ = 0 to matrows - 1 do
    let i = Test_matrix.rand () mod matrows in
    let j = Test_matrix.rand () mod matcols in
    Matrix.Dense.set mc i j (rand_float ())
  done;
  let a = Matrix.wrap_sparse (Matrix.Sparse.from_dense mattype 0.0 mc) in
  let b = Matrix.wrap_sparse (Matrix.Sparse.from_dense mattype 0.0 md) in

  (* Create vectors and fill *)
  let x = Nvector_serial.make matcols 0.0
  and y = Nvector_serial.make matrows 0.0
  and z = Nvector_serial.make matrows 0.0 in
  let xdata = Nvector.unwrap x in
  for i = 0 to matcols - 1 do
    xdata.{i} <- rand_float ()
  done;

  (try Matrix.matvec c x y
    with _ -> printf "FAIL: SUNMatrix module Dense matvec failure \n \n";
              exit 1);
  (try Matrix.matvec d x z
    with _ -> printf "FAIL: SUNMatrix module Dense matvec failure \n \n";
              exit 1);

  (* SUNMatrix Tests *)
  fails += Test.test_sunmatgetid a Matrix.Sparse 0;
  fails += Test.test_sunmatclone a 0;
  fails += Test.test_sunmatcopy a 0;
  fails += Test.test_sunmatzero a 0;
  fails += Test.test_sunmatscaleadd a i 0;
  (match Sundials.Config.sundials_version with
   | 2,_,_ | 3,1,0 | 3,1,1 -> ()
   | _ -> fails += test_sunmatscaleadd2 SparseTests.check_vector a b x y z);
  if square then begin
    fails += Test.test_sunmatscaleaddi a i 0;
    (match Sundials.Config.sundials_version with
     | 2,_,_ | 3,1,0 | 3,1,1 -> ()
     | _ -> fails += test_sunmatscaleaddi2 SparseTests.check_vector a x y)
  end;
  fails += Test.test_sunmatmatvec a x y 0;
  fails += Test.test_sunmatspace a 0;
  (match Sundials.Config.sundials_version with
   | x, _, _ when x < 5 -> ()
   | _ -> fails += test_sunsparsematrix_convert SparseTests.check_matrix a);

  (* Print result *)
  if !fails <> 0 then begin
    printf "FAIL: SUNMatrix module failed %d tests @\n @\n" !fails;
    printf "@\nA = %a@\n" Matrix.Sparse.pp (Matrix.unwrap a);
    (match Sundials.Config.sundials_version with
     | 2,_,_ | 3,1,0 | 3,1,1 -> ()
     | _ -> printf "@\nB = %a@\n" Matrix.Sparse.pp (Matrix.unwrap b));
    if square then
      printf "@\nI = %a@\n" Matrix.Sparse.pp (Matrix.unwrap i);
    printf "@\nx = %a@\n" Nvector_serial.pp x;
    printf "@\ny = %a@\n" Nvector_serial.pp y;
    (match Sundials.Config.sundials_version with
     | 2,_,_ | 3,1,0 | 3,1,1 -> ()
     | _ -> printf "@\nz = %a@\n" Nvector_serial.pp z)
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
  for _ = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()

