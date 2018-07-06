(*
 * -----------------------------------------------------------------
 * $Revision: 4463 $
 * $Date: 2015-03-29 16:28:20 -0700 (Sun, 29 Mar 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: These testing routines are based on an
 *                   NVECTOR testing routine by Daniel R. Reynolds
 *                   @ SMU.
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, AIST, Jan 2016.
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * These test functions are designed to check an NVECTOR module
 * implementation.
 *
 * NOTE: Many of these tests rely on the N_VGetArrayPointer routine
 *       to get a pointer to the data component of an N_Vector. This
 *       assumes the internal data is stored in a contiguous
 *       realtype array.
 * -----------------------------------------------------------------
 *)

(* NOTE: a number of tests in the original C code validate memory
   management or functions that have no sensible OCaml counterpart.
   Those tests are omitted in this port.  These differences make the
   tests NOT suitable for performance comparison with C.  However, the
   output should still exactly agree with C if print_time is
   disabled.  *)

open Sundials

module type NVECTOR_OPS_EXT = sig
  include Nvector.NVECTOR_OPS
  val n_vgetarray : t -> RealArray.t
end

module Test (Nvector_ops : NVECTOR_OPS_EXT) =
struct (* Extends all the way to the end of file.  *)

let (+=) r x = r := !r + x
let printf = Printf.printf
let int_of_bool = function
  | true -> 1
  | false -> 0

exception TestFailed of int

(* define constatnts *)
let neg_two = -2.0
let neg_one = -1.0
let neg_half = -0.5
let zero = 0.0
let half = 0.5
let one = 1.0
let two = 2.0

(* NAN and floating point "equality" check, failure update macro *)
external stdc_version : unit -> int = "stdc_version"
let isnan (x : float) = x <> x

let fneq =
  if stdc_version () >= 199901 then
    (fun a b ->
      int_of_bool (isnan a || abs_float (a -. b) /. abs_float b > 1.0e-15))
  else
    (fun a b ->
      int_of_bool (abs_float (a -. b) /. abs_float b > 1.0e-15))


(* private functions *)

let print_time_flag = ref false

let print_time format time =
  if !print_time_flag then
    printf format time

external get_time : unit -> float = "get_time"
external set_timing : bool -> float = "SetTiming"
let set_timing b =
  print_time_flag := b;
  set_timing b

(* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*)
let check_ans ans x local_length =
  let failure = ref 0 in
  let xdata = Nvector_ops.n_vgetarray x in

  (* check vector data *)
  for i=0 to local_length-1 do
    failure += fneq xdata.{i} ans
  done;

  int_of_bool (!failure > 0)

(* ----------------------------------------------------------------------
 * N_VCloneVectorArray Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_n_vclonevectorarray count w local_length myid =
  (* clone array of vectors *)
  let start_time = get_time () in
  let _ = Array.init count (fun _ -> Nvector_ops.n_vclone w) in
  let stop_time = get_time () in

  (* check array of vectors *)
  (* check vectors in array *)

  if myid = 0 then (
    printf "    PASSED test -- N_VCloneVectorArray \n";
    print_time "    N_VCloneVectorArray Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  0

(* ----------------------------------------------------------------------
 * N_VCloneVectorArrayEmpty Test
 * --------------------------------------------------------------------*)
let test_n_vcloneemptyvectorarray count w myid =
  (* clone empty array *)
  let start_time = get_time () in
  let _ = Array.init count (fun _ -> Nvector_ops.n_vclone w) in
  let stop_time = get_time () in

  (* check array of vectors *)
  (* check vectors in array *)

  if myid = 0 then (
    printf "    PASSED test -- N_VCloneEmptyVectorArray \n";
    print_time "    N_VCloneEmptyVectorArray Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  0


(* ----------------------------------------------------------------------
 * N_VCloneEmpty Test
 * --------------------------------------------------------------------*)
let test_n_vcloneempty w myid =
  (* clone empty vector *)
  let start_time = get_time () in
  let _ = Nvector_ops.n_vclone w in
  let stop_time = get_time () in

  (* check vector *)
  (* check vector data *)
  if myid = 0 then (
    printf "    PASSED test -- N_VCloneEmpty \n";
    print_time "    N_VCloneEmpty Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  0


(* ----------------------------------------------------------------------
 * N_VClone Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_n_vclone w local_length myid =
  (* clone vector *)
  let start_time = get_time () in
  let x = Nvector_ops.n_vclone w in
  let stop_time = get_time () in

  (* check cloned vector *)

  (* check cloned vector data *)
  let _ = Nvector_ops.n_vgetarray x in

  Nvector_ops.n_vconst one x;
  let failure = check_ans one x local_length in
  if failure <> 0 then (
    printf ">>> FAILED test -- N_VClone Proc %d \n" myid;
    printf "    Failed N_VConst check \n \n";
    1
  ) else (
    if myid = 0 then (
      printf "    PASSED test -- N_VClone \n";
      print_time "    N_VClone Time: %22.15e \n \n" (stop_time -. start_time)
    );
    0
  )


(* ----------------------------------------------------------------------
 * N_VGetArrayPointer Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_n_vgetarraypointer w local_length myid =
  let failure = ref 0 in
  (* get vector data *)
  let start_time = get_time () in
  let wdata = Nvector_ops.n_vgetarray w in
  let stop_time = get_time () in

  (* check vector data *)
  Nvector_ops.n_vconst neg_half w;
  for i=0 to local_length-1 do
    failure += fneq wdata.{i} neg_half;
  done;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VGetArrayPointer Proc %d \n" myid;
    printf "    Failed N_VConst check \n \n";
    raise (TestFailed 1)
  );

  if myid = 0 then (
    printf "    PASSED test -- N_VGetArrayPointer \n";
    print_time "    N_VGetArrayPointer Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  0


(* ----------------------------------------------------------------------
 * N_VSetArrayPointer Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_n_vsetarraypointer w local_length myid =
  let failure = ref 0 in

  (* C version creates a buffer of the same length and sets that as
     the storage for w.  There's no counterpart in OCaml, so we instead
     clone w, throw that away, and get the storage for w.  *)

  (* create vector data *)
  let _ = Nvector_ops.n_vclone w in

  (* attach data to vector *)
  let start_time = get_time () in
  let wdata = Nvector_ops.n_vgetarray w in
  let stop_time = get_time () in

  (* check vector data *)
  Nvector_ops.n_vconst neg_half w;
  for i=0 to local_length-1 do
    failure += fneq wdata.{i} neg_half;
  done;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VSetArrayPointer Proc %d \n" myid;
    printf "    Failed N_VConst check \n \n";
    raise (TestFailed 1)
  );

  if myid = 0 then (
    printf "    PASSED test -- N_VSetArrayPointer \n";
    print_time "    N_VSetArrayPointer Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  0


(* ----------------------------------------------------------------------
 * N_VLinearSum Tests
 * --------------------------------------------------------------------*)
let test_n_vlinearsum x y z local_length myid =
  let fails = ref 0
  and failure = ref 0 in

  let xdata = Nvector_ops.n_vgetarray x in
  let ydata = Nvector_ops.n_vgetarray y in
  let zdata = Nvector_ops.n_vgetarray z in

  (*
   * Case 1a: y = x + y, (Vaxpy Case 1)
   *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- one;
    ydata.{i} <- neg_two;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum one x one y y;
  let stop_time = get_time () in

  (* y should be vector of -1 *)
  failure := check_ans neg_one y local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 1a Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 1a \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 1b: y = -x + y, (Vaxpy Case 2)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- one;
    ydata.{i} <- two;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum neg_one x one y y;
  let stop_time = get_time () in

  (* y should be vector of +1 *)
  failure := check_ans one y local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 1b Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 1b \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 1c: y = ax + y, (Vaxpy Case 3)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- neg_two;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum half x one y y;
  let stop_time = get_time () in

  (* y should be vector of -1 *)
  failure := check_ans neg_one y local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 1c Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 1c \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 2a: x = x + y, (Vaxpy Case 1)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- neg_one;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum one x one y x;
  let stop_time = get_time () in

  (* y should be vector of +1 *)
  failure := check_ans one x local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 2a Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 2a \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 2b: x = x - y, (Vaxpy Case 2)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- one;
    ydata.{i} <- two;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum one x neg_one y x;
  let stop_time = get_time () in

  (* y should be vector of -1 *)
  failure := check_ans neg_one x local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 2b Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 2b \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 2c: x = x + by, (Vaxpy Case 3)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- neg_half;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum one x two y x;
  let stop_time = get_time () in

  (* x should be vector of +1 *)
  failure := check_ans one x local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 2c Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 2c \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 3: z = x + y, (VSum)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_two;
    ydata.{i} <- one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum one x one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 3 Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 3 \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n" (stop_time -. start_time)
  );

  (*
   * Case 4a: z = x - y, (VDiff)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum one x neg_one y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  failure := check_ans one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 4a Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 4a \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 4b: z = -x + y, (VDiff)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum neg_one x one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 4b Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 4b \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 5a: z = x + by, (VLin1)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- neg_half;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum one x two y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  failure := check_ans one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 5a Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 5a \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 5b: z = ax + y, (VLin1)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- half;
    ydata.{i} <- neg_two;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum two x one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 5b Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 5b \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 6a: z = -x + by, (VLin2)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_two;
    ydata.{i} <- neg_half;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum neg_one x two y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  failure := check_ans one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 6a Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 6a \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 6b: z = ax - y, (VLin2)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- half;
    ydata.{i} <- two;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum two x neg_one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 6b Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 6b \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 7: z = a(x + y), (VScaleSum)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- one;
    ydata.{i} <- neg_half;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum two x two y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  failure := check_ans one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 7 Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 7 \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 8: z = a(x - y), (VScaleDiff)
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- half;
    ydata.{i} <- one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum two x neg_two y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 8 Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 8 \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 9: z = ax + by, All Other Cases
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- one;
    ydata.{i} <- neg_two;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vlinearsum two x half y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  failure := check_ans one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VLinearSum Case 9 Proc %d \n" myid;
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VLinearSum Case 9 \n";
    print_time "    N_VLinearSum Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VConst Test
 * --------------------------------------------------------------------*)
let test_n_vconst x local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.n_vgetarray x in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vconst one x;
  let stop_time = get_time () in

  (* x should be vector of +1 *)
  let failure = check_ans one x local_length in

  if failure <> 0 then (
    printf ">>> FAILED test -- N_VConst Proc %d \n" myid;
    print_time "    N_VConst Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VConst \n";
    print_time "    N_VConst Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VProd Test
 * --------------------------------------------------------------------*)
let test_n_vprod x y z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let ydata = Nvector_ops.n_vgetarray y in
  let zdata = Nvector_ops.n_vgetarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- neg_half;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vprod x y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VProd Proc %d \n" myid;
    print_time "    N_VProd Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VProd \n";
    print_time "    N_VProd Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VDiv Test
 * --------------------------------------------------------------------*)
let test_n_vdiv x y z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let ydata = Nvector_ops.n_vgetarray y in
  let zdata = Nvector_ops.n_vgetarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- one;
    ydata.{i} <- two;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vdiv x y z;
  let stop_time = get_time () in

  (* z should be vector of +1/2 *)
  failure := check_ans half z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VDiv Proc %d \n" myid;
    print_time "    N_VDiv Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VDiv \n";
    print_time "    N_VDiv Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VScale Tests
 * --------------------------------------------------------------------*)
let test_n_vscale x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let zdata = Nvector_ops.n_vgetarray z in

  (*
   * Case 1: x = cx, VScaleBy
   *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- half;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vscale two x x;
  let stop_time = get_time () in

  (* x should be vector of +1 *)
  failure := check_ans one x local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VScale Case 1 Proc %d \n" myid;
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VScale Case 1 \n";
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time)
  );

  (*
   * Case 2: z = x, VCopy
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vscale one x z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VScale Case 2 Proc %d \n" myid;
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VScale Case 2 \n";
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time)
  );

  (*
   * Case 3: z = -x, VNeg
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vscale neg_one x z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  failure := check_ans one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VScale Case 3 Proc %d \n" myid;
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VScale Case 3 \n";
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time)
  );

  (*
   * Case 4: z = cx, All other cases
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_half;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vscale two x z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VScale Case 4 Proc %d \n" myid;
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VScale Case 4 \n";
    print_time "    N_VScale Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VAbs Test
 * --------------------------------------------------------------------*)
let test_n_vabs x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let zdata = Nvector_ops.n_vgetarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vabs x z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  failure := check_ans one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VAbs Proc %d \n" myid;
    print_time "    N_VAbs Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VAbs \n";
    print_time "    N_VAbs Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VInv Test
 * --------------------------------------------------------------------*)
let test_n_vinv x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let zdata = Nvector_ops.n_vgetarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vinv x z;
  let stop_time = get_time () in

  (* z should be vector of +1/2 *)
  failure := check_ans half z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VInv Proc %d \n" myid;
    print_time "    N_VInv Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VInv \n";
    print_time "    N_VInv Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VAddConst Test
 * --------------------------------------------------------------------*)
let test_n_vaddconst x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let zdata = Nvector_ops.n_vgetarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- one;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  Nvector_ops.n_vaddconst x neg_two z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  failure := check_ans neg_one z local_length;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VAddConst Proc %d \n" myid;
    print_time "    N_VAddConst Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VAddConst \n";
    print_time "    N_VAddConst Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VDotProd Test
 * --------------------------------------------------------------------*)
let test_n_vdotprod x y local_length global_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let ydata = Nvector_ops.n_vgetarray y in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
    ydata.{i} <- half;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vdotprod x y in
  let stop_time = get_time () in

  (* ans should equal global vector length *)
  failure := fneq ans (float_of_int global_length);

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VDotProd Proc %d \n" myid;
    print_time "    N_VDotProd Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VDotProd \n";
    print_time "    N_VDotProd Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VMaxNorm Test
 * --------------------------------------------------------------------*)
let test_n_vmaxnorm x local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_one;
  done;
  xdata.{local_length-1} <- neg_two;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vmaxnorm x in
  let stop_time = get_time () in

  (* ans should equal 2 *)
  failure := if ans < zero then 1 else fneq ans two;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VMaxNorm Proc %d \n" myid;
    print_time "    N_VMaxNorm Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VMaxNorm \n";
    print_time "    N_VMaxNorm Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VWrmsNorm Test
 * --------------------------------------------------------------------*)
let test_n_vwrmsnorm x w local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let wdata = Nvector_ops.n_vgetarray w in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_half;
    wdata.{i} <- half;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vwrmsnorm x w in
  let stop_time = get_time () in

  (* ans should equal 1/4 *)
  failure := if ans < zero then 1 else fneq ans (half*.half);

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VWrmsNorm Proc %d \n" myid;
    print_time "    N_VWrmsNorm Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VWrmsNorm \n";
    print_time "    N_VWrmsNorm Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VWrmsNormMask Test
 * --------------------------------------------------------------------*)
let test_n_vwrmsnormmask x w id local_length global_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let wdata = Nvector_ops.n_vgetarray w in
  let id_data = Nvector_ops.n_vgetarray id in

  (*
   * Case 1: use all elements, id = 1
   *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_half;
    wdata.{i} <- half;
    id_data.{i} <- one;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vwrmsnormmask x w id in
  let stop_time = get_time () in

  (* ans equals 1/4 (same as wrms norm) *)
  failure := if ans < zero then 1 else fneq ans (half*.half);

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VWrmsNormMask Case 1 Proc %d \n" myid;
    print_time "    N_VWrmsNormMask Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VWrmsNormMask Case 1 \n";
    print_time "    N_VWrmsNormMask Time: %22.15e \n \n" (stop_time -. start_time)
  );

  (*
   * Case 2: use no elements, id = 0
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_half;
    wdata.{i} <- half;
    id_data.{i} <- zero;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vwrmsnormmask x w id in
  let stop_time = get_time () in

  (* ans equals 0 (skips all elements) *)
  failure := if ans < zero then 1 else fneq ans zero;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VWrmsNormMask Case 2 Proc %d \n" myid;
    print_time "    N_VWrmsNormMask Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VWrmsNormMask Case 2 \n";
    print_time "    N_VWrmsNormMask Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VMin Test
 * --------------------------------------------------------------------*)
let test_n_vmin x local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- two;
  done;
  xdata.{local_length-1} <- neg_one;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vmin x in
  let stop_time = get_time () in

  (* ans should equal -1 *)
  failure := fneq ans neg_one;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VMin Proc %d \n" myid;
    print_time "    N_VMin Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VMin \n";
    print_time "    N_VMin Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VWL2Norm Test
 * --------------------------------------------------------------------*)
let test_n_vwl2norm x w local_length global_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in
  let wdata = Nvector_ops.n_vgetarray w in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_half;
    wdata.{i} <- half;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vwl2norm x w in
  let stop_time = get_time () in

  (* ans should equal 1/4 * sqrt(global_length) *)
  failure :=
    if ans < zero then 1
    else fneq ans (half*.half*.sqrt(float_of_int global_length));

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VWL2Norm Proc %d \n" myid;
    print_time "    N_VWL2Norm Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VWL2Norm \n";
    print_time "    N_VWL2Norm Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VL1Norm Test
 * --------------------------------------------------------------------*)
let test_n_vl1norm x local_length global_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let xdata = Nvector_ops.n_vgetarray x in

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- neg_one;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vl1norm x in
  let stop_time = get_time () in

  (* ans should equal global_length *)
  failure := if ans < zero then 1 else fneq ans (float_of_int global_length);

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VL1Norm Proc %d \n" myid;
    print_time "    N_VL1Norm Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VL1Norm \n";
    print_time "    N_VL1Norm Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VCompare
 * --------------------------------------------------------------------*)
let test_n_vcompare x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  if local_length < 3 then (
    printf "Error Test_N_VCompare: Local vector length is %d length must be >= 3" local_length;
    raise (TestFailed (-1))
  );

  let xdata = Nvector_ops.n_vgetarray x in
  let zdata = Nvector_ops.n_vgetarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    zdata.{i} <- neg_one;

    match i mod 3 with
    | 0 ->
      (* abs(x[i]) < c *)
      xdata.{i} <- zero

    | 1 ->
      (* abs(x[i]) = c *)
      xdata.{i} <- neg_one

    | 2 ->
      (* abs(x[i]) > c *)
      xdata.{i} <- neg_two
    | _ -> assert false
  done;

  let start_time = get_time () in
  Nvector_ops.n_vcompare one x z;
  let stop_time = get_time () in

  (* check return vector *)
  for i=0 to local_length-1 do
    match i mod 3 with
    | 0 ->
      (* z[i] == 0 *)
      if zdata.{i} <> zero then
	failure := 1
    | 1 ->
      (* z[i] == 1 *)
      if zdata.{i} <> one then
	failure := 1
    | 2 ->
      (* z[i] == 1 *)
      if zdata.{i} <> one then
	failure := 1
    | _ -> assert false
  done;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VCompare Proc %d \n" myid;
    print_time "    N_VCompare Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VCompare \n";
    print_time "    N_VCompare Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VInvTest
 * --------------------------------------------------------------------*)
let test_n_vinvtest x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  if local_length < 2 then (
    printf "Error Test_N_VCompare: Local vector length is %d length must be >= 2" local_length;
    raise (TestFailed (-1))
  );

  let xdata = Nvector_ops.n_vgetarray x in
  let zdata = Nvector_ops.n_vgetarray z in

  (*
   * Case 1: All elements Nonzero, z[i] = 1/x[i], return True
   *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    xdata.{i} <- half;
    zdata.{i} <- zero;
  done;

  let start_time = get_time () in
  let test = Nvector_ops.n_vinvtest x z in
  let stop_time = get_time () in

  (* z should be vector of +2 *)
  failure := check_ans two z local_length;

  if !failure <> 0 || not test then (
    printf ">>> FAILED test -- N_VInvTest Case 1, Proc %d \n" myid;
    print_time "    N_VInvTest Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VInvTest Case 1 \n";
    print_time "    N_VInvTest Time: %22.15e \n \n" (stop_time -. start_time)
  );

  (*
   * Case 2: Some elements Zero, z[i] = 1/x[i] for x[i] != 0, return False
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    zdata.{i} <- zero;

    if i mod 2 <> 0 then
      xdata.{i} <- half
    else
      xdata.{i} <- zero
  done;

  let start_time = get_time () in
  let test = Nvector_ops.n_vinvtest x z in
  let stop_time = get_time () in

  (* check return vector *)
  for i=0 to local_length-1 do
    if i mod 2 <> 0 then (
      if zdata.{i} <> two then
	failure := 1
    ) else (
      if zdata.{i} <> zero then
	failure := 1
      );
  done;

  if !failure <> 0 || test then (
    printf ">>> FAILED test -- N_VInvTest Case 2 Proc %d \n" myid;
    print_time "    N_VInvTest Time: %22.15e \n \n" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VInvTest Case 2 \n";
    print_time "    N_VInvTest Time: %22.15e \n \n" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VConstrMask
 * --------------------------------------------------------------------*)
let test_n_vconstrmask c x m local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  if local_length < 7 then (
    printf "Error Test_N_VCompare: Local vector length is %d length must be >= 7" local_length;
    raise (TestFailed (-1))
  );

  let cdata = Nvector_ops.n_vgetarray c in
  let xdata = Nvector_ops.n_vgetarray x in
  let mdata = Nvector_ops.n_vgetarray m in

  (*
   * Case 1: Return True
   *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    mdata.{i} <- neg_one;

    match i mod 7 with
    | 0 ->
      (* c = -2, test for < 0*)
      cdata.{i} <- neg_two;
      xdata.{i} <- neg_two

    | 1 ->
      (* c = -1, test for <= 0 *)
      cdata.{i} <- neg_one;
      xdata.{i} <- neg_one

    | 2 ->
      (* c = -1, test for == 0 *)
      cdata.{i} <- neg_one;
      xdata.{i} <- zero

    | 3 ->
      (* c = 0, no test *)
      cdata.{i} <- zero;
      xdata.{i} <- half

    | 4 ->
      (* c = 1, test for == 0*)
      cdata.{i} <- one;
      xdata.{i} <- zero

    | 5 ->
      (* c = 1, test for >= 0*)
      cdata.{i} <- one;
      xdata.{i} <- one

    | 6 ->
      (* c = 2, test for > 0 *)
      cdata.{i} <- two;
      xdata.{i} <- two

    | _ -> assert false
  done;

  let start_time = get_time () in
  let test = Nvector_ops.n_vconstrmask c x m in
  let stop_time = get_time () in

  (* m should be vector of 0 *)
  failure := check_ans zero m local_length;

  if !failure <> 0 || not test then (
    printf ">>> FAILED test -- N_VConstrMask Case 1 Proc %d \n" myid;
    print_time "    N_VConstrMask Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VConstrMask Case 1 \n";
    print_time "    N_VConstrMask Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 2: Return False
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    mdata.{i} <- neg_one;

    match i mod 5 with
    | 0 ->
      (* c = -2, test for < 0*)
      cdata.{i} <- neg_two;
      xdata.{i} <- two

    | 1 ->
      (* c = -1, test for <= 0 *)
      cdata.{i} <- neg_one;
      xdata.{i} <- one

    | 2 ->
      (* c = 0, no test *)
      cdata.{i} <- zero;
      xdata.{i} <- half

    | 3 ->
      (* c = 1, test for >= 0*)
      cdata.{i} <- one;
      xdata.{i} <- neg_one

    | 4 ->
      (* c = 2, test for > 0 *)
      cdata.{i} <- two;
      xdata.{i} <- neg_two

    | _ -> assert false
  done;

  let start_time = get_time () in
  let test = Nvector_ops.n_vconstrmask c x m in
  let stop_time = get_time () in

  (* check mask vector *)
  for i=0 to local_length-1 do
    if i mod 5 = 2 then (
      if mdata.{i} <> zero then
	failure := 1
    ) else (
      if mdata.{i} <> one then
	failure := 1
    );
  done;

  if !failure <> 0 || test then (
    printf ">>> FAILED test -- N_VConstrMask Case 2, Proc %d \n" myid;
    print_time "    N_VConstrMask Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VConstrMask Case 2 \n";
    print_time "    N_VConstrMask Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VMinQuotient Test
 * --------------------------------------------------------------------*)
let test_n_vminquotient num denom local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  let num_data = Nvector_ops.n_vgetarray num in
  let denom_data = Nvector_ops.n_vgetarray denom in

  (*
   * Case 1: Pass
   *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    num_data.{i} <- two;
    denom_data.{i} <- two;
  done;
  num_data.{local_length-1} <- one;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vminquotient num denom in
  let stop_time = get_time () in

  (* ans should equal 1/2 *)
  failure := fneq ans half;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VMinQuotient Case 1 Proc %d \n" myid;
    print_time "    N_VMinQuotient Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VMinQuotient Case 1 \n";
    print_time "    N_VMinQuotient Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  (*
   * Case 2: Fail
   *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    num_data.{i} <- two;
    denom_data.{i} <- zero;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.n_vminquotient num denom in
  let stop_time = get_time () in

  (* ans should equal big_real *)
  failure := fneq ans Config.big_real;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VMinQuotient Case 2 Proc %d \n" myid;
    print_time "    N_VMinQuotient Time: %22.15e \n \n"
      (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    printf "    PASSED test -- N_VMinQuotient Case 2 \n";
    print_time "    N_VMinQuotient Time: %22.15e \n \n"
      (stop_time -. start_time)
  );

  !fails

end
