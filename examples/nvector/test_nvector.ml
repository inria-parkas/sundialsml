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
 * OCaml port: Timothy Bourke, Inria, Mar 2020.
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

let compat_ge400 =
  let n, _, _ = Config.sundials_version in
  n >= 4

let compat_ge500 =
  let n, _, _ = Config.sundials_version in
  n >= 5

let compat_ge540 =
  let n, m, _ = Config.sundials_version in
  (n = 5 && m >= 4) || n >= 6

let compat_ge570 =
  let n, m, _ = Config.sundials_version in
  (n = 5 && m >= 7) || n >= 6

let compat_ge600 =
  let n, _, _ = Config.sundials_version in
  n >= 6

let compat_neq600 =
  match Config.sundials_version with
  | 6, 0, 0 -> false
  | _ -> true

let print_passed s =
  if compat_ge400
  then Printf.printf "PASSED test -- %s \n" s
  else Printf.printf "    PASSED test -- %s \n" s

module type NVECTOR_OPS_EXT = sig
  include Nvector.NVECTOR_OPS
  type contents
  val get_id   : t -> Nvector.nvector_id
  val getarray : t -> contents

  val get : contents -> int -> float
  val set : contents -> int -> float -> unit

  val max_time : t -> float -> float
  val sync_device : unit -> unit
end

module type TEST = sig (* {{{ *)
  type t

  exception TestFailed of int

  val set_timing : bool -> bool -> unit
  val test_getvectorid : t -> Nvector.nvector_id -> int -> int
  val test_getlength                : t -> int -> int
  val test_getcommunicator          : t -> 'comm -> int -> int
  val test_clonevectorarray         : int -> t -> 'a -> int -> int
  val test_cloneemptyvectorarray    : int -> t -> int -> int
  val test_cloneempty               : t -> int -> int
  val test_clone                    : t -> int -> int -> int
  val test_getarraypointer          : t -> int -> int -> int
  val test_setarraypointer          : t -> int -> int -> int
  val test_linearsum                : t -> t -> t -> int -> int -> int
  val test_const                    : t -> int -> int -> int
  val test_prod                     : t -> t -> t -> int -> int -> int
  val test_div                      : t -> t -> t -> int -> int -> int
  val test_scale                    : t -> t -> int -> int -> int
  val test_abs                      : t -> t -> int -> int -> int
  val test_inv                      : t -> t -> int -> int -> int
  val test_addconst                 : t -> t -> int -> int -> int
  val test_dotprod                  : t -> t -> int -> int -> int -> int
  val test_maxnorm                  : t -> int -> int -> int
  val test_wrmsnorm                 : t -> t -> int -> int -> int
  val test_wrmsnormmask             : t -> t -> t -> int -> int -> int -> int
  val test_wrmsnormmask_lt400       : t -> t -> t -> int -> 'a -> int -> int
  val test_min                      : t -> int -> int -> int
  val test_wl2norm                  : t -> t -> int -> int -> int -> int
  val test_l1norm                   : t -> int -> int -> int -> int
  val test_compare                  : t -> t -> int -> int -> int
  val test_invtest                  : t -> t -> int -> int -> int
  val test_constrmask               : t -> t -> t -> int -> int -> int
  val test_minquotient              : t -> t -> int -> int -> int
  val test_linearcombination        : t -> int -> int -> int
  val test_scaleaddmulti            : t -> int -> int -> int
  val test_dotprodmulti             : t -> 'a -> int -> int -> int
  val test_linearsumvectorarray     : t -> int -> int -> int
  val test_scalevectorarray         : t -> int -> int -> int
  val test_constvectorarray         : t -> int -> int -> int
  val test_wrmsnormvectorarray      : t -> 'a -> int -> int
  val test_wrmsnormmaskvectorarray  : t -> int -> int -> int -> int
  val test_scaleaddmultivectorarray : t -> int -> int -> int
  val test_linearcombinationvectorarray : t -> int -> int -> int

  val test_dotprodlocal     : t -> t -> int -> int -> int
  val test_maxnormlocal     : t -> int -> int -> int
  val test_minlocal         : t -> int -> int -> int
  val test_l1normlocal      : t -> int -> int -> int
  val test_wsqrsumlocal     : t -> t -> int -> int -> int
  val test_wsqrsummasklocal : t -> t -> t -> int -> int -> int
  val test_invtestlocal     : t -> t -> int -> int -> int
  val test_constrmasklocal  : t -> t -> t -> int -> int -> int
  val test_minquotientlocal : t -> t -> int -> int -> int

  val test_dotprodmultilocal : t -> int -> int -> int
  val test_dotprodmultiallreduce : t -> int -> int -> int

  val test_bufsize   : t -> int -> int -> int
  val test_bufpack   : t -> int -> int -> int
  val test_bufunpack : t -> int -> int -> int
end (* }}} *)

module Test (Nvector_ops : NVECTOR_OPS_EXT) : TEST with type t = Nvector_ops.t =
struct (* Extends all the way to the end of file.  *)
type t = Nvector_ops.t

let (+=) r x = r := !r + x
let printf = Printf.printf

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

let fneqtol =
  if stdc_version () >= 199901 then
    (fun a b tol ->
      not (isnan a || abs_float (a -. b) /. abs_float b > tol))
  else
    (fun a b tol ->
      not (abs_float (a -. b) /. abs_float b > tol))

let fneq a b =
  fneqtol a b (if compat_ge500 then 10.0 *. Config.unit_roundoff else 1.0e-15)

(* private functions *)

let print_time_flag = ref false

let print_time name time =
  if !print_time_flag then
    if compat_ge400
    then printf "%s Time: %22.15e \n \n" name time
    else printf "    %s Time: %22.15e \n \n" name time

external get_time : unit -> float = "get_time"

(* May return: tick precision (in nanoseconds, as float) *)
external set_timing : bool -> (int * float) option = "SetTiming"

let set_timing b showres =
  (match set_timing b with
   | Some (n, f) ->
     if compat_ge400 && showres
     then printf "Timer resolution: %d ns = %g s\n" n f
   | None -> ());
  print_time_flag := b

(* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*)
let check_ans ans x local_length =
  let xdata = Nvector_ops.getarray x in

  (* check vector data *)
  try
    for i=0 to local_length-1 do
      if not (fneq (Nvector_ops.get xdata i) ans) then raise Exit
    done; true
  with Exit -> false

(* ----------------------------------------------------------------------
 * N_VGetVectorID Test
 * --------------------------------------------------------------------*)

let int_of_nvector_id = function
    Nvector.Serial        -> 0
  | Nvector.Parallel      -> 1
  | Nvector.OpenMP        -> 2
  | Nvector.Pthreads      -> 3
  | Nvector.ParHyp        -> 4
  | Nvector.PETSc         -> 5
  | Nvector.CUDA          -> 6
  | Nvector.HIP           -> 7
  | Nvector.SYCL          -> 8
  | Nvector.RAJA          -> 9
  | Nvector.Kokkos        -> 10
  | Nvector.OpenMPdev     -> 11
  | Nvector.Trilinos      -> 12
  | Nvector.ManyVector    -> 13
  | Nvector.MpiManyVector -> 14
  | Nvector.MpiPlusX      -> 15
  | Nvector.Custom        ->
      if Sundials_impl.Version.lt500 then 9 else 16

let test_getvectorid x id myid =
  if Nvector_ops.get_id x <> id then (
    printf ">>> FAILED test -- N_VGetVectorID, Proc %d \n" myid;
    printf "    Unrecognized vector type %d \n \n"
           (int_of_nvector_id (Nvector_ops.get_id x));
    1
  ) else (
    if myid = 0 then print_passed "N_VGetVectorID";
    0
  )

(* ----------------------------------------------------------------------
 * Test_N_VGetLength Test
 * --------------------------------------------------------------------*)
let test_getlength w myid =
  (* ask W for it's overall length *)
  let wlength = Nvector_ops.getlength w in
  (* use N_VConst and N_VDotProd to compute length *)
  Nvector_ops.const 1.0 w;
  let wlength2 = int_of_float (Nvector_ops.dotprod w w) in
  Nvector_ops.sync_device ();
  (* return error if lengths disagree *)
  if wlength <> wlength2 then (
    printf ">>> FAILED test -- N_VGetLength, Proc %d (%d != %d)\n"
      myid wlength wlength2;
    1
  ) else (
    if myid = 0 then printf "PASSED test -- N_VGetLength\n";
    0
  )

(* ----------------------------------------------------------------------
 * Test_N_VGetCommunicator Test (without MPI dependency)
 * --------------------------------------------------------------------*)
let test_getcommunicator _ _ myid =
  (* Cannot test without MPI *)
  if myid = 0 then printf "PASSED test -- N_VGetCommunicator\n";
  0

(* ----------------------------------------------------------------------
 * N_VCloneVectorArray Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_clonevectorarray count w _ myid =
  (* clone array of vectors *)
  let start_time = get_time () in
  let _ = Array.init count (fun _ -> Nvector_ops.clone w) in
  let stop_time = get_time () in

  (* check array of vectors *)
  (* check vectors in array *)

  if myid = 0 then (
    print_passed "N_VCloneVectorArray";
    print_time "N_VCloneVectorArray" (stop_time -. start_time)
  );

  0

(* ----------------------------------------------------------------------
 * N_VCloneVectorArrayEmpty Test
 * --------------------------------------------------------------------*)
let test_cloneemptyvectorarray count w myid =
  (* clone empty array *)
  let start_time = get_time () in
  let _ = Array.init count (fun _ -> Nvector_ops.clone w) in
  let stop_time = get_time () in

  (* check array of vectors *)
  (* check vectors in array *)

  if myid = 0 then (
    print_passed "N_VCloneEmptyVectorArray";
    print_time "N_VCloneEmptyVectorArray" (stop_time -. start_time)
  );

  0


(* ----------------------------------------------------------------------
 * N_VCloneEmpty Test
 * --------------------------------------------------------------------*)
let test_cloneempty w myid =
  (* clone empty vector *)
  let start_time = get_time () in
  let _ = Nvector_ops.clone w in
  let stop_time = get_time () in

  (* check vector *)
  (* check vector data *)
  if myid = 0 then (
    print_passed "N_VCloneEmpty";
    print_time "N_VCloneEmpty" (stop_time -. start_time)
  );

  0


(* ----------------------------------------------------------------------
 * N_VClone Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_clone w local_length myid =
  (* clone vector *)
  let start_time = get_time () in
  let x = Nvector_ops.clone w in
  let stop_time = get_time () in

  (* check cloned vector *)

  (* check cloned vector data *)
  let _ = Nvector_ops.getarray x in

  Nvector_ops.const one x;
  if not (check_ans one x local_length) then (
    printf ">>> FAILED test -- N_VClone Proc %d \n" myid;
    printf "    Failed N_VConst check \n \n";
    1
  ) else (
    if myid = 0 then (
      print_passed "N_VClone";
      print_time "N_VClone" (stop_time -. start_time)
    );
    0
  )


(* ----------------------------------------------------------------------
 * N_VGetArrayPointer Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_getarraypointer w local_length myid =
  (* get vector data *)
  let start_time = get_time () in
  let wdata = Nvector_ops.getarray w in
  let stop_time = get_time () in

  (* check vector data *)
  Nvector_ops.const neg_half w;
  try
    for i=0 to local_length-1 do
      if not (fneq (Nvector_ops.get wdata i) neg_half) then raise Exit
    done;
    if myid = 0 then (
      print_passed "N_VGetArrayPointer";
      print_time "N_VGetArrayPointer" (stop_time -. start_time)
    );
    0
  with Exit -> begin
    printf ">>> FAILED test -- N_VGetArrayPointer Proc %d \n" myid;
    printf "    Failed N_VConst check \n \n";
    raise (TestFailed 1)
  end

(* ----------------------------------------------------------------------
 * N_VSetArrayPointer Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*)
let test_setarraypointer w local_length myid =
  (* C version creates a buffer of the same length and sets that as
     the storage for w.  There's no counterpart in OCaml, so we instead
     clone w, throw that away, and get the storage for w.  *)

  (* create vector data *)
  let _ = Nvector_ops.clone w in

  (* attach data to vector *)
  let start_time = get_time () in
  let wdata = Nvector_ops.getarray w in
  let stop_time = get_time () in

  (* check vector data *)
  Nvector_ops.const neg_half w;
  try
    for i=0 to local_length-1 do
      if not (fneq (Nvector_ops.get wdata i) neg_half) then raise Exit
    done;
    if myid = 0 then (
      print_passed "N_VSetArrayPointer";
      print_time "N_VSetArrayPointer" (stop_time -. start_time)
    );
    0
  with Exit -> begin
    printf ">>> FAILED test -- N_VSetArrayPointer Proc %d \n" myid;
    printf "    Failed N_VConst check \n \n";
    raise (TestFailed 1)
  end

(* ----------------------------------------------------------------------
 * N_VLinearSum Tests
 * --------------------------------------------------------------------*)
let test_linearsum x y z local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x
  and ydata = Nvector_ops.getarray y
  and zdata = Nvector_ops.getarray z
  in

  (* Case 1a: y = x + y, (Vaxpy Case 1) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i one;
    Nvector_ops.set ydata i neg_two
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum one x one y y;
  let stop_time = get_time () in

  (* y should be vector of -1 *)
  if not (check_ans neg_one y local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 1a Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 1a";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 1b: y = -x + y, (Vaxpy Case 2) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i one;
    Nvector_ops.set ydata i two
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum neg_one x one y y;
  let stop_time = get_time () in

  (* y should be vector of +1 *)
  if not (check_ans one y local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 1b Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 1b";
    print_time "N_VLinearSum"
      (stop_time -. start_time)
  );

  (* Case 1c: y = ax + y, (Vaxpy Case 3) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set ydata i neg_two
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum half x one y y;
  let stop_time = get_time () in

  (* y should be vector of -1 *)
  if not (check_ans neg_one y local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 1c Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 1c";
    print_time "N_VLinearSum Time" (stop_time -. start_time)
  );

  (* Case 2a: x = x + y, (Vaxpy Case 1) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set ydata i neg_one
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum one x one y x;
  let stop_time = get_time () in

  (* y should be vector of +1 *)
  if not (check_ans one x local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 2a Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 2a";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 2b: x = x - y, (Vaxpy Case 2) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i one;
    Nvector_ops.set ydata i two
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum one x neg_one y x;
  let stop_time = get_time () in

  (* y should be vector of -1 *)
  if not (check_ans neg_one x local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 2b Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 2b";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 2c: x = x + by, (Vaxpy Case 3) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set ydata i neg_half
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum one x two y x;
  let stop_time = get_time () in

  (* x should be vector of +1 *)
  if not (check_ans one x local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 2c Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 2c";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 3: z = x + y, (VSum) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_two;
    Nvector_ops.set ydata i one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum one x one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 3 Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 3";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 4a: z = x - y, (VDiff) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set ydata i one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum one x neg_one y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  if not (check_ans one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 4a Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 4a";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 4b: z = -x + y, (VDiff) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set ydata i one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum neg_one x one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 4b Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 4b";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 5a: z = x + by, (VLin1) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set ydata i neg_half;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum one x two y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  if not (check_ans one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 5a Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 5a";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 5b: z = ax + y, (VLin1) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i half;
    Nvector_ops.set ydata i neg_two;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum two x one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 5b Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 5b";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 6a: z = -x + by, (VLin2) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_two;
    Nvector_ops.set ydata i neg_half;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum neg_one x two y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  if not (check_ans one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 6a Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 6a";
    print_time "N_VLinearSum Time" (stop_time -. start_time)
  );

  (* Case 6b: z = ax - y, (VLin2) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i half;
    Nvector_ops.set ydata i two;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum two x neg_one y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 6b Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 6b";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 7: z = a(x + y), (VScaleSum) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i one;
    Nvector_ops.set ydata i neg_half;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum two x two y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  if not (check_ans one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 7 Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 7";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 8: z = a(x - y), (VScaleDiff) *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i half;
    Nvector_ops.set ydata i one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum two x neg_two y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 8 Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 8";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  (* Case 9: z = ax + by, All Other Cases *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i one;
    Nvector_ops.set ydata i neg_two;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.linearsum two x half y z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  if not (check_ans one z local_length) then (
    printf ">>> FAILED test -- N_VLinearSum Case 9 Proc %d \n" myid;
    print_time "N_VLinearSum" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VLinearSum Case 9";
    print_time "N_VLinearSum" (stop_time -. start_time)
  );

  !fails

(* ----------------------------------------------------------------------
 * N_VConst Test
 * --------------------------------------------------------------------*)
let test_const x local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.const one x;
  let stop_time = get_time () in

  (* x should be vector of +1 *)
  if not (check_ans one x local_length) then (
    printf ">>> FAILED test -- N_VConst Proc %d \n" myid;
    print_time "N_VConst" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VConst";
    print_time "N_VConst" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VProd Test
 * --------------------------------------------------------------------*)
let test_prod x y z local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x
  and ydata = Nvector_ops.getarray y
  and zdata = Nvector_ops.getarray z
  in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set ydata i neg_half;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.prod x y z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VProd Proc %d \n" myid;
    print_time "N_VProd" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VProd";
    print_time "N_VProd" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VDiv Test
 * --------------------------------------------------------------------*)
let test_div x y z local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let ydata = Nvector_ops.getarray y in
  let zdata = Nvector_ops.getarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i one;
    Nvector_ops.set ydata i two;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.div x y z;
  let stop_time = get_time () in

  (* z should be vector of +1/2 *)
  if not (check_ans half z local_length) then (
    printf ">>> FAILED test -- N_VDiv Proc %d \n" myid;
    print_time "N_VDiv" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VDiv";
    print_time "N_VDiv" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VScale Tests
 * --------------------------------------------------------------------*)
let test_scale x z local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let zdata = Nvector_ops.getarray z in

  (* Case 1: x = cx, VScaleBy *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i half
  done;

  let start_time = get_time () in
  Nvector_ops.scale two x x;
  let stop_time = get_time () in

  (* x should be vector of +1 *)
  if not (check_ans one x local_length) then (
    printf ">>> FAILED test -- N_VScale Case 1 Proc %d \n" myid;
    print_time "N_VScale" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VScale Case 1";
    print_time "N_VScale" (stop_time -. start_time)
  );

  (* Case 2: z = x, VCopy *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.scale one x z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VScale Case 2 Proc %d \n" myid;
    print_time "N_VScale" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VScale Case 2";
    print_time "N_VScale" (stop_time -. start_time)
  );

  (* Case 3: z = -x, VNeg *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.scale neg_one x z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  if not (check_ans one z local_length) then (
    printf ">>> FAILED test -- N_VScale Case 3 Proc %d \n" myid;
    print_time "N_VScale" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VScale Case 3";
    print_time "N_VScale" (stop_time -. start_time)
  );

  (* Case 4: z = cx, All other cases *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_half;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.scale two x z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VScale Case 4 Proc %d \n" myid;
    print_time "N_VScale" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VScale Case 4";
    print_time "N_VScale" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VAbs Test
 * --------------------------------------------------------------------*)
let test_abs x z local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let zdata = Nvector_ops.getarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.abs x z;
  let stop_time = get_time () in

  (* z should be vector of +1 *)
  if not (check_ans one z local_length) then (
    printf ">>> FAILED test -- N_VAbs Proc %d \n" myid;
    print_time "N_VAbs" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VAbs";
    print_time "N_VAbs" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VInv Test
 * --------------------------------------------------------------------*)
let test_inv x z local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let zdata = Nvector_ops.getarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.inv x z;
  let stop_time = get_time () in

  (* z should be vector of +1/2 *)
  if not (check_ans half z local_length) then (
    printf ">>> FAILED test -- N_VInv Proc %d \n" myid;
    print_time "N_VInv" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VInv";
    print_time "N_VInv" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VAddConst Test
 * --------------------------------------------------------------------*)
let test_addconst x z local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let zdata = Nvector_ops.getarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i one;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  Nvector_ops.addconst x neg_two z;
  let stop_time = get_time () in

  (* z should be vector of -1 *)
  if not (check_ans neg_one z local_length) then (
    printf ">>> FAILED test -- N_VAddConst Proc %d \n" myid;
    print_time "N_VAddConst" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VAddConst";
    print_time "N_VAddConst" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VDotProd Test
 * --------------------------------------------------------------------*)
let test_dotprod x y _ global_length myid =
  let fails = ref 0 in

  (* fill vector data *)
  Nvector_ops.const two x;
  Nvector_ops.const half y;

  let start_time = get_time () in
  let ans = Nvector_ops.dotprod x y in
  let stop_time = get_time () in

  (* ans should equal global vector length *)
  if not (fneq ans (float_of_int global_length)) then (
    printf ">>> FAILED test -- N_VDotProd Proc %d \n" myid;
    print_time "N_VDotProd" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VDotProd";
    print_time "N_VDotProd" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VMaxNorm Test
 * --------------------------------------------------------------------*)
let test_maxnorm x local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in

  if Sundials_impl.Version.lt640 then begin

    (* fill vector data *)
    for i=0 to local_length-1 do
      Nvector_ops.set xdata i neg_one
    done;
    Nvector_ops.set xdata (local_length-1) neg_two;

    let start_time = get_time () in
    let ans = Nvector_ops.maxnorm x in
    let stop_time = get_time () in

    (* ans should equal 2 *)
    if (ans < zero || not (fneq ans two)) then (
      printf ">>> FAILED test -- N_VMaxNorm Proc %d \n" myid;
      print_time "N_VMaxNorm" (stop_time -. start_time);
      fails += 1
    ) else if myid = 0 then (
      print_passed "N_VMaxNorm";
      print_time "N_VMaxNorm" (stop_time -. start_time)
    )

  end else begin
    (* Case 1: zero vector *)

    (* fill vector data *)
    for i=0 to local_length-1 do
      Nvector_ops.set xdata i zero
    done;

    let start_time = get_time () in
    let ans = Nvector_ops.maxnorm x in
    let stop_time = get_time () in

    (* ans should equal 0 *)
    if (ans < zero || ans >= Sundials.Config.small_real) then (
      printf ">>> FAILED test -- N_VMaxNorm Case 1, Proc %d \n" myid;
      print_time "N_VMaxNorm" (stop_time -. start_time);
      fails += 1
    ) else if myid = 0 then (
      printf "PASSED test -- N_VMaxNorm Case 1\n";
      print_time "N_VMaxNorm" (stop_time -. start_time));

    (* Case 2: general vector *)

    (* reset failure *)

    (* fill vector data *)
    for i=0 to local_length-1 do
      Nvector_ops.set xdata i neg_half
    done;
    if myid = 0
    then Nvector_ops.set xdata (local_length - 1) neg_two
    else Nvector_ops.set xdata (local_length - 1) one;

    let start_time = get_time () in
    let ans = Nvector_ops.maxnorm x in
    let stop_time = get_time () in

    (* ans should equal 2 *)
    if ans < zero || not (fneq ans two) then (
      printf ">>> FAILED test -- N_VMaxNorm Case 2, Proc %d \n" myid;
      print_time "N_VMaxNorm" (stop_time -. start_time);
      fails += 1
    ) else if myid = 0 then (
      printf "PASSED test -- N_VMaxNorm Case 2\n";
      print_time "N_VMaxNorm" (stop_time -. start_time))
  end;

  !fails


(* ----------------------------------------------------------------------
 * N_VWrmsNorm Test
 * --------------------------------------------------------------------*)
let test_wrmsnorm x w local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let wdata = Nvector_ops.getarray w in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_half;
    Nvector_ops.set wdata i half;
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.wrmsnorm x w in
  let stop_time = get_time () in

  (* ans should equal 1/4 *)
  if ans < zero || not (fneq ans (half*.half)) then (
    printf ">>> FAILED test -- N_VWrmsNorm Proc %d \n" myid;
    print_time "N_VWrmsNorm" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VWrmsNorm";
    print_time "N_VWrmsNorm" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VWrmsNormMask Test
 * --------------------------------------------------------------------*)
let test_wrmsnormmask x w id local_length global_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let wdata = Nvector_ops.getarray w in
  let id_data = Nvector_ops.getarray id in

  (* factor used in checking solutions *)
  let fac = sqrt(float (global_length - 1) /. float global_length) in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_half;
    Nvector_ops.set wdata i half;
    Nvector_ops.set id_data i one
  done;
  if myid = 0 then Nvector_ops.set id_data (local_length-1) zero;

  let start_time = get_time () in
  let ans = Nvector_ops.wrmsnormmask x w id in
  let stop_time = get_time () in

  (* ans equals 1/4 (same as wrms norm) *)
  if ans < zero || not (fneq ans (fac*.half*.half)) then (
    printf ">>> FAILED test -- N_VWrmsNormMask, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VWrmsNormMask";
  );

  print_time "N_VWrmsNormMask" (stop_time -. start_time);
  !fails

let test_wrmsnormmask_lt400 x w id local_length _ myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let wdata = Nvector_ops.getarray w in
  let id_data = Nvector_ops.getarray id in

  (* Case 1: use all elements, id = 1 *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_half;
    Nvector_ops.set wdata i half;
    Nvector_ops.set id_data i one
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.wrmsnormmask x w id in
  let stop_time = get_time () in

  (* ans equals 1/4 (same as wrms norm) *)
  if ans < zero || not (fneq ans (half*.half)) then (
    printf ">>> FAILED test -- N_VWrmsNormMask Case 1 Proc %d \n" myid;
    print_time "N_VWrmsNormMask" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VWrmsNormMask Case 1";
    print_time "N_VWrmsNormMask" (stop_time -. start_time)
  );

  (* Case 2: use no elements, id = 0 *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_half;
    Nvector_ops.set wdata i half;
    Nvector_ops.set id_data i zero
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.wrmsnormmask x w id in
  let stop_time = get_time () in

  (* ans equals 0 (skips all elements) *)
  if ans < zero || not (fneq ans zero) then (
    printf ">>> FAILED test -- N_VWrmsNormMask Case 2 Proc %d \n" myid;
    print_time "N_VWrmsNormMask" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VWrmsNormMask Case 2";
    print_time "N_VWrmsNormMask" (stop_time -. start_time)
  );

  !fails

(* ----------------------------------------------------------------------
 * N_VMin Test
 * --------------------------------------------------------------------*)
let test_min x local_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i two;
  done;
  Nvector_ops.set xdata (local_length-1) neg_one;

  let start_time = get_time () in
  let ans = Nvector_ops.min x in
  let stop_time = get_time () in

  (* ans should equal -1 *)
  if not (fneq ans neg_one) then (
    printf ">>> FAILED test -- N_VMin Proc %d \n" myid;
    print_time "N_VMin" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VMin";
    print_time "N_VMin" (stop_time -. start_time)
  );

  if compat_ge570 then begin
    (* fill vector data *)
    Nvector_ops.const 2.0 x;
    let xdata = Nvector_ops.getarray x in
    if myid = 0
    then Nvector_ops.set xdata (local_length-1) neg_two
    else Nvector_ops.set xdata (local_length-1) neg_one;

    let start_time = get_time () in
    let ans = Nvector_ops.min x in
    let stop_time = get_time () in

    (* ans should equal -2 *)
    if not (fneq ans neg_two) then (
      printf ">>> FAILED test -- N_VMin, Proc %d \n" myid;
      printf "    min = %f, expected %f\n" ans (-2.0);
      fails += 1
    ) else if myid = 0 then (
      printf "PASSED test -- N_VMin \n";
    );

    (* find max time across all processes *)
    print_time "N_VMin"
      (Nvector_ops.max_time x (stop_time -. start_time));
  end;

  !fails


(* ----------------------------------------------------------------------
 * N_VWL2Norm Test
 * --------------------------------------------------------------------*)
let test_wl2norm x w local_length global_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in
  let wdata = Nvector_ops.getarray w in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_half;
    Nvector_ops.set wdata i half
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.wl2norm x w in
  let stop_time = get_time () in

  (* ans should equal 1/4 * sqrt(global_length) *)
  if ans < zero
     || not (fneq ans (half*.half*.sqrt(float_of_int global_length)))
  then (
    printf ">>> FAILED test -- N_VWL2Norm Proc %d \n" myid;
    print_time "N_VWL2Norm" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VWL2Norm";
    print_time "N_VWL2Norm" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VL1Norm Test
 * --------------------------------------------------------------------*)
let test_l1norm x local_length global_length myid =
  let fails = ref 0 in

  let xdata = Nvector_ops.getarray x in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i neg_one
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.l1norm x in
  let stop_time = get_time () in

  (* ans should equal global_length *)
  if ans < zero || not (fneq ans (float_of_int global_length)) then (
    printf ">>> FAILED test -- N_VL1Norm Proc %d \n" myid;
    print_time "N_VL1Norm" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VL1Norm";
    print_time "N_VL1Norm" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VCompare
 * --------------------------------------------------------------------*)
let test_compare x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  if local_length < 3 then (
    printf "Error Test_N_VCompare: Local vector length is %d length must be >= 3\n" local_length;
    raise (TestFailed (-1))
  );

  let xdata = Nvector_ops.getarray x in
  let zdata = Nvector_ops.getarray z in

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set zdata i neg_one;

    match i mod 3 with
    | 0 ->
      (* abs(x[i]) < c *)
      Nvector_ops.set xdata i zero

    | 1 ->
      (* abs(x[i]) = c *)
      Nvector_ops.set xdata i neg_one

    | 2 ->
      (* abs(x[i]) > c *)
      Nvector_ops.set xdata i neg_two
    | _ -> assert false
  done;

  let start_time = get_time () in
  Nvector_ops.compare one x z;
  let stop_time = get_time () in

  (* check return vector *)
  for i=0 to local_length-1 do
    match i mod 3 with
    | 0 ->
      (* z[i] == 0 *)
      if Nvector_ops.get zdata i <> zero then
        failure := 1
    | 1 ->
      (* z[i] == 1 *)
      if Nvector_ops.get zdata i <> one then
        failure := 1
    | 2 ->
      (* z[i] == 1 *)
      if Nvector_ops.get zdata i <> one then
        failure := 1
    | _ -> assert false
  done;

  if !failure <> 0 then (
    printf ">>> FAILED test -- N_VCompare Proc %d \n" myid;
    print_time "N_VCompare" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VCompare";
    print_time "N_VCompare" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VInvTest
 * --------------------------------------------------------------------*)
let test_invtest x z local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  if local_length < 2 then (
    printf "Error Test_N_VCompare: Local vector length is %d length must be >= 2\n" local_length;
    raise (TestFailed (-1))
  );

  let xdata = Nvector_ops.getarray x in
  let zdata = Nvector_ops.getarray z in

  (* Case 1: All elements Nonzero, z[i] = 1/x[i], return True *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set xdata i half;
    Nvector_ops.set zdata i zero
  done;

  let start_time = get_time () in
  let test = Nvector_ops.invtest x z in
  let stop_time = get_time () in

  (* z should be vector of +2 *)
  if not (check_ans two z local_length && test) then (
    printf ">>> FAILED test -- N_VInvTest Case 1, Proc %d \n" myid;
    print_time "N_VInvTest" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VInvTest Case 1";
    print_time "N_VInvTest" (stop_time -. start_time)
  );

  (* Case 2: Some elements Zero, z[i] = 1/x[i] for x[i] != 0, return False *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set zdata i zero;

    if i mod 2 <> 0 then
      Nvector_ops.set xdata i half
    else
      Nvector_ops.set xdata i zero
  done;

  let start_time = get_time () in
  let test = Nvector_ops.invtest x z in
  let stop_time = get_time () in

  (* check return vector *)
  for i=0 to local_length-1 do
    if i mod 2 <> 0 then (
      if Nvector_ops.get zdata i <> two then
        failure := 1
    ) else (
      if Nvector_ops.get zdata i <> zero then
        failure := 1
      );
  done;

  if !failure <> 0 || test then (
    printf ">>> FAILED test -- N_VInvTest Case 2 Proc %d \n" myid;
    print_time "N_VInvTest" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VInvTest Case 2";
    print_time "N_VInvTest" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VConstrMask
 * --------------------------------------------------------------------*)
let test_constrmask c x m local_length myid =
  let fails = ref 0
  and failure = ref 0
  in

  if local_length < 7 then (
    printf "Error Test_N_VCompare: Local vector length is %d length must be >= 7\n" local_length;
    raise (TestFailed (-1))
  );

  let cdata = Nvector_ops.getarray c in
  let xdata = Nvector_ops.getarray x in
  let mdata = Nvector_ops.getarray m in

  (* Case 1: Return True *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set mdata i neg_one;

    match i mod 7 with
    | 0 ->
      (* c = -2, test for < 0*)
      Nvector_ops.set cdata i neg_two;
      Nvector_ops.set xdata i neg_two

    | 1 ->
      (* c = -1, test for <= 0 *)
      Nvector_ops.set cdata i neg_one;
      Nvector_ops.set xdata i neg_one

    | 2 ->
      (* c = -1, test for == 0 *)
      Nvector_ops.set cdata i neg_one;
      Nvector_ops.set xdata i zero

    | 3 ->
      (* c = 0, no test *)
      Nvector_ops.set cdata i zero;
      Nvector_ops.set xdata i half

    | 4 ->
      (* c = 1, test for == 0*)
      Nvector_ops.set cdata i one;
      Nvector_ops.set xdata i zero

    | 5 ->
      (* c = 1, test for >= 0*)
      Nvector_ops.set cdata i one;
      Nvector_ops.set xdata i one

    | 6 ->
      (* c = 2, test for > 0 *)
      Nvector_ops.set cdata i two;
      Nvector_ops.set xdata i two

    | _ -> assert false
  done;

  let start_time = get_time () in
  let test = Nvector_ops.constrmask c x m in
  let stop_time = get_time () in

  (* m should be vector of 0 *)
  if not (check_ans zero m local_length && test) then (
    printf ">>> FAILED test -- N_VConstrMask Case 1 Proc %d \n" myid;
    print_time "N_VConstrMask" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VConstrMask Case 1";
    print_time "N_VConstrMask" (stop_time -. start_time)
  );

  (* Case 2: Return False *)

  (* reset failure *)
  failure := 0;

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set mdata i neg_one;

    match i mod 5 with
    | 0 ->
      (* c = -2, test for < 0*)
      Nvector_ops.set cdata i neg_two;
      Nvector_ops.set xdata i two

    | 1 ->
      (* c = -1, test for <= 0 *)
      Nvector_ops.set cdata i neg_one;
      Nvector_ops.set xdata i one

    | 2 ->
      (* c = 0, no test *)
      Nvector_ops.set cdata i zero;
      Nvector_ops.set xdata i half

    | 3 ->
      (* c = 1, test for >= 0*)
      Nvector_ops.set cdata i one;
      Nvector_ops.set xdata i neg_one

    | 4 ->
      (* c = 2, test for > 0 *)
      Nvector_ops.set cdata i two;
      Nvector_ops.set xdata i neg_two

    | _ -> assert false
  done;

  let start_time = get_time () in
  let test = Nvector_ops.constrmask c x m in
  let stop_time = get_time () in

  (* check mask vector *)
  for i=0 to local_length-1 do
    if i mod 5 = 2 then (
      if Nvector_ops.get mdata i <> zero then
        failure := 1
    ) else (
      if Nvector_ops.get mdata i <> one then
        failure := 1
    );
  done;

  if !failure <> 0 || test then (
    printf ">>> FAILED test -- N_VConstrMask Case 2, Proc %d \n" myid;
    print_time "N_VConstrMask" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VConstrMask Case 2";
    print_time "N_VConstrMask" (stop_time -. start_time)
  );

  !fails


(* ----------------------------------------------------------------------
 * N_VMinQuotient Test
 * --------------------------------------------------------------------*)
let test_minquotient num denom local_length myid =
  let fails = ref 0 in

  let num_data = Nvector_ops.getarray num in
  let denom_data = Nvector_ops.getarray denom in

  (* Case 1: Pass *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set num_data i two;
    Nvector_ops.set denom_data i two;
  done;
  Nvector_ops.set num_data (local_length-1) one;

  let start_time = get_time () in
  let ans = Nvector_ops.minquotient num denom in
  let stop_time = get_time () in

  (* ans should equal 1/2 *)
  if not (fneq ans half) then (
    printf ">>> FAILED test -- N_VMinQuotient Case 1 Proc %d \n" myid;
    print_time "N_VMinQuotient" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VMinQuotient Case 1";
    print_time "N_VMinQuotient" (stop_time -. start_time)
  );

  (* Case 2: Fail *)

  (* fill vector data *)
  for i=0 to local_length-1 do
    Nvector_ops.set num_data i two;
    Nvector_ops.set denom_data i zero
  done;

  let start_time = get_time () in
  let ans = Nvector_ops.minquotient num denom in
  let stop_time = get_time () in

  (* ans should equal big_real *)
  if not (fneq ans Config.big_real) then (
    printf ">>> FAILED test -- N_VMinQuotient Case 2 Proc %d \n" myid;
    print_time "N_VMinQuotient" (stop_time -. start_time);
    fails += 1
  ) else if myid = 0 then (
    print_passed "N_VMinQuotient Case 2";
    print_time "N_VMinQuotient" (stop_time -. start_time)
  );

  !fails

(* ----------------------------------------------------------------------
 * N_VLinearCombination Test
 * --------------------------------------------------------------------*)
let test_linearcombination x local_length myid =
  let fails = ref 0 in

  (* create vectors for testing *)
  (* set vectors in vector array *)
  let v = Array.init 3 (fun _ -> Nvector_ops.clone x) in
  let v_len1, v_len2, v_len3 = Array.sub v 0 1, Array.sub v 0 2, v in
  let y1, y2, y3 = v.(0), v.(1), v.(2) in

  (* initialize c values *)
  let c = RealArray.make 3 zero in

  (* Case 1a: v.(0) = a v.(0), N_VScale *)
  Nvector_ops.const two y1; (* fill vector data *)
  c.{0} <- half; (* set scaling factors *)

  let start_time = get_time () in
  Nvector_ops.linearcombination c v_len1 y1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* Y1 should be vector of +1 *)
  if not (check_ans one y1 local_length) then (
    printf ">>> FAILED test -- N_VLinearCombination Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombination Case 1a";

  (* find max time across all processes *)
  print_time "N_VLinearCombination"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 1b: X = a v.(0), N_VScale *)

  (* fill vector data and scaling factors *)
  Nvector_ops.const two y1;
  c.{0} <- half;
  Nvector_ops.const zero x;

  let start_time = get_time () in
  Nvector_ops.linearcombination c v_len1 x;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* X should be vector of +1 *)
  if not (check_ans one x local_length) then (
    printf ">>> FAILED test -- N_VLinearCombination Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombination Case 1b";

  (* find max time across all processes *)
  print_time "N_VLinearCombination"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2a: v.(0) = a v.(0) + b v.(1), N_VLinearSum *)
  Nvector_ops.const neg_two y1; (* fill vector data *)
  Nvector_ops.const one y2;
  c.{0} <- half;                    (* set scaling factors *)
  c.{1} <- two;

  let start_time = get_time () in
  Nvector_ops.linearcombination c v_len2 y1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* Y1 should be vector of +1 *)
  if not (check_ans one y1 local_length) then (
    printf ">>> FAILED test -- N_VLinearCombination Case 2a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombination Case 2a";

  (* find max time across all processes *)
  print_time "N_VLinearCombination"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2b: X = a v.(0) + b v.(1), N_VLinearSum *)
  Nvector_ops.const one y1;      (* fill vector data and scaling factors *)
  Nvector_ops.const neg_two y2;
  c.{0} <- two;
  c.{1} <- half;

  Nvector_ops.const zero x;

  let start_time = get_time () in
  Nvector_ops.linearcombination c v_len2 x;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* X should be vector of +1 *)
  if not (check_ans one y1 local_length) then (
    printf ">>> FAILED test -- N_VLinearCombination Case 2b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombination Case 2b";

  (* find max time across all processes *)
  print_time "N_VLinearCombination"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 3a: v.(0) = v.(0) + b v.(1) + c v.(2) *)
  Nvector_ops.const two y1;        (* fill vector data *)
  Nvector_ops.const neg_two y2;
  Nvector_ops.const neg_one y3;
  c.{0} <- one;                        (* set scaling factors *)
  c.{1} <- half;
  c.{2} <- neg_two;

  let start_time = get_time () in
  Nvector_ops.linearcombination c v_len3 y1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* Y1 should be vector of +3 *)
  if not (check_ans (two+.one) y1 local_length) then (
    printf ">>> FAILED test -- N_VLinearCombination Case 3a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombination Case 3a";

  (* find max time across all processes *)
  print_time "N_VLinearCombination"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 3b: v.(0) = a v.(0) + b v.(1) + c v.(2) *)
  Nvector_ops.const one y1;        (* fill vector data *)
  Nvector_ops.const neg_two y2;
  Nvector_ops.const neg_one y3;
  c.{0} <- two;                        (* set scaling factors *)
  c.{1} <- half;
  c.{2} <- neg_one;

  let start_time = get_time () in
  Nvector_ops.linearcombination c v_len3 y1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* Y1 should be vector of +2 *)
  if not (check_ans two y1 local_length) then (
    printf ">>> FAILED test -- N_VLinearCombination Case 3b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombination Case 3b";

  (* find max time across all processes *)
  print_time "N_VLinearCombination"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 3c: X = a v.(0) + b v.(1) + c v.(2) *)
  Nvector_ops.const one y1;    (* fill vector data and set scaling factors *)
  Nvector_ops.const neg_two y2;
  Nvector_ops.const neg_one y3;
  c.{0} <- two;
  c.{1} <- half;
  c.{2} <- neg_one;

  Nvector_ops.const zero x;

  let start_time = get_time () in
  Nvector_ops.linearcombination c v_len3 x;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* X should be vector of +2 *)
  if not (check_ans two x local_length) then (
    printf ">>> FAILED test -- N_VLinearCombination Case 3c, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombination Case 3c";

  (* find max time across all processes *)
  print_time "N_VLinearCombination"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VScaleaddmulti Test
 * --------------------------------------------------------------------*)
let test_scaleaddmulti x local_length myid =
  let fails = ref 0 in

  let avals = RealArray.make 3 zero in
  let avals1 = RealArray.sub avals 0 1
  and z = Array.init 3 (fun _ -> Nvector_ops.clone x)
  and v = Array.init 3 (fun _ -> Nvector_ops.clone x)
  in
  let v_len1, _, v_len3 = Array.sub v 0 1, Array.sub v 0 2, v
  and z_len1 = Array.sub z 0 1
  in

  (* Case 1a: v.(0) = a.(0) x + v.(0), N_VLinearSum *)
  Nvector_ops.const one x;           (* fill vector data *)
  Nvector_ops.const neg_one v.(0);
  avals.{0} <- two;                     (* set scaling factors *)

  let start_time = get_time () in
  Nvector_ops.scaleaddmulti avals1 x v_len1 v_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* v.(0) should be vector of +1 *)
  if not (check_ans one v.(0) local_length) then (
    printf ">>> FAILED test -- N_VScaleAddMulti Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMulti Case 1a";

  (* find max time across all processes *)
  print_time "N_VScaleAddMulti"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 1b: z.(0) = a.(0) x + v.(0), N_VLinearSum *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const one x;
  Nvector_ops.const neg_one v.(0);
  avals.{0} <- two;

  Nvector_ops.const zero z.(0);

  let start_time = get_time () in
  Nvector_ops.scaleaddmulti avals1 x v_len1 z_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(0) should be vector of +1 *)
  if not (check_ans one z.(0) local_length) then (
    printf ">>> FAILED test -- N_VScaleAddMulti Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMulti Case 1b";

  (* find max time across all processes *)
  print_time "N_VScaleAddMulti"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2a: v.(i) = a.(i) x + v.(i), N_VScaleAddMulti *)
  Nvector_ops.const one x;             (* fill vector data *)
  Nvector_ops.const neg_two v.(0);
  Nvector_ops.const two v.(1);
  Nvector_ops.const neg_one v.(2);
  avals.{0} <- one;                        (* set scaling factors *)
  avals.{1} <- neg_two;
  avals.{2} <- two;

  let start_time = get_time () in
  Nvector_ops.scaleaddmulti avals x v_len3 v_len3;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* v.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one v.(0) local_length
          && check_ans zero    v.(1) local_length
          && check_ans one     v.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMulti Case 2a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMulti Case 2a";

  (* find max time across all processes *)
  print_time "N_VScaleAddMulti"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2b: z.(i) = a.(i) x + v.(i), N_VScaleAddMulti *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const one x;
  Nvector_ops.const neg_two v.(0);
  Nvector_ops.const two v.(1);
  Nvector_ops.const neg_one v.(2);
  avals.{0} <- one;
  avals.{1} <- neg_two;
  avals.{2} <- two;

  Nvector_ops.const two z.(0);
  Nvector_ops.const two z.(1);
  Nvector_ops.const two z.(2);

  let start_time = get_time () in
  Nvector_ops.scaleaddmulti avals x v_len3 z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMulti Case 2b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMulti Case 2b";

  (* find max time across all processes *)
  print_time "N_VScaleAddMulti"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VDotProdMulti Test
 * --------------------------------------------------------------------*)
let test_dotprodmulti x _ global_length myid =
  let fails = ref 0
  and dotprods = RealArray.make 3 zero
  in
  let dotprods1 = RealArray.sub dotprods 0 1 in

  (* create vectors for testing *)
  let v = Array.init 3 (fun _ -> Nvector_ops.clone x) in
  let v_len1, _, v_len3 = Array.sub v 0 1, Array.sub v 0 2, v in

  (* Case 1: d.(0) = z . v.(0), N_VDotProd *)

  (* fill vector data *)
  Nvector_ops.const two x;
  Nvector_ops.const half v.(0);

  let start_time = get_time () in
  Nvector_ops.dotprodmulti x v_len1 dotprods1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod.(0) should equal the global vector length *)
  if not (fneq dotprods.{0} (float global_length)) then (
    printf ">>> FAILED test -- N_VDotProdMulti Case 1, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VDotProdMulti Case 1";

  (* find max time across all processes *)
  print_time "N_VDotProdMulti"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2: d.(i) = z . v.(i), N_VDotProd *)

  (* fill vector data *)
  Nvector_ops.const two x;
  Nvector_ops.const neg_half v.(0);
  Nvector_ops.const half v.(1);
  Nvector_ops.const one v.(2);

  let start_time = get_time () in
  Nvector_ops.dotprodmulti x v_len3 dotprods;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod.(i) should equal -1, +1, and 2 times the global vector length *)
  if not (   fneq dotprods.{0} (float (-1*global_length))
          && fneq dotprods.{1} (float (   global_length))
          && fneq dotprods.{2} (float ( 2*global_length)))
  then (
    printf ">>> FAILED test -- N_VDotProdMulti Case 2, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VDotProdMulti Case 2";

  (* find max time across all processes *)
  print_time "N_VDotProdMulti"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VLinearSumVectorArray Test
 * --------------------------------------------------------------------*)
let test_linearsumvectorarray v local_length myid =
  let fails = ref 0 in

  (* create vectors for testing *)
  let x = Array.init 3 (fun _ -> Nvector_ops.clone v) in
  let y = Array.init 3 (fun _ -> Nvector_ops.clone v) in
  let z = Array.init 3 (fun _ -> Nvector_ops.clone v) in
  let x_len1, y_len1, z_len1 = Array.sub x 0 1, Array.sub y 0 1, Array.sub z 0 1
  in

  (* Case 0: z.(0) = a x.(0) + b y.(0), N_VLinearSum *)

  (* fill vector data *)
  Nvector_ops.const neg_half x.(0);
  Nvector_ops.const two y.(0);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray two x_len1 half y_len1 z_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(0) should be a vector of 0 *)
  if not (check_ans zero z.(0) local_length) then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 0, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 0";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 1a: y.(i) = x.(i) + y.(i), (VaxpyVectorArray Case 1) *)

  (* fill vector data *)
  Nvector_ops.const neg_two x.(0);
  Nvector_ops.const one y.(0);

  Nvector_ops.const two x.(1);
  Nvector_ops.const neg_two y.(1);

  Nvector_ops.const two x.(2);
  Nvector_ops.const neg_one y.(2);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray one x one y y;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one y.(0) local_length
          && check_ans zero    y.(1) local_length
          && check_ans one     y.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 1a";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 1b: y = -x + y, (VaxpyVectorArray Case 2) *)

  (* fill vector data *)
  Nvector_ops.const two x.(0);
  Nvector_ops.const one y.(0);

  Nvector_ops.const neg_two x.(1);
  Nvector_ops.const neg_two y.(1);

  Nvector_ops.const neg_two x.(2);
  Nvector_ops.const neg_one y.(2);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray neg_one x one y y;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one y.(0) local_length
          && check_ans zero    y.(1) local_length
          && check_ans one     y.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 1b";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 1c: y = ax + y, (VaxpyVectorArray Case 3) *)

  (* fill vector data *)
  Nvector_ops.const two x.(0);
  Nvector_ops.const neg_two y.(0);

  Nvector_ops.const two x.(1);
  Nvector_ops.const neg_one y.(1);

  Nvector_ops.const neg_two x.(2);
  Nvector_ops.const two y.(2);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray half x one y y;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one y.(0) local_length
          && check_ans zero    y.(1) local_length
          && check_ans one     y.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 1c, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 1c";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 2a: x = x + y, (VaxpyVectorArray Case 1) *)

  (* fill vector data *)
  Nvector_ops.const neg_two x.(0);
  Nvector_ops.const one y.(0);

  Nvector_ops.const two x.(1);
  Nvector_ops.const neg_two y.(1);

  Nvector_ops.const two x.(2);
  Nvector_ops.const neg_one y.(2);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray one x one y x;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one x.(0) local_length
          && check_ans zero    x.(1) local_length
          && check_ans one     x.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 2a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 2a";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 2b: x = x - y, (VaxpyVectorArray Case 2) *)

  (* fill vector data *)
  Nvector_ops.const one x.(0);
  Nvector_ops.const two y.(0);

  Nvector_ops.const neg_two x.(1);
  Nvector_ops.const neg_two y.(1);

  Nvector_ops.const neg_one x.(2);
  Nvector_ops.const neg_two y.(2);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray one x neg_one y x;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one x.(0) local_length
          && check_ans zero    x.(1) local_length
          && check_ans one     x.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 2b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 2b";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 2c: x = x + by, (VaxpyVectorArray Case 3) *)

  (* fill vector data *)
  Nvector_ops.const neg_two x.(0);
  Nvector_ops.const two y.(0);

  Nvector_ops.const neg_one x.(1);
  Nvector_ops.const two y.(1);

  Nvector_ops.const two x.(2);
  Nvector_ops.const neg_two y.(2);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray one x half y x;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one x.(0) local_length
          && check_ans zero    x.(1) local_length
          && check_ans one     x.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 2c, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 2c";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 3: z = x + y, (VSumVectorArray) *)

  (* fill vector data *)
  Nvector_ops.const neg_two x.(0);
  Nvector_ops.const one y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const neg_one x.(1);
  Nvector_ops.const one y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const two x.(2);
  Nvector_ops.const neg_one y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray one x one y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 3, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 3";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 4a: z = x - y, (VDiffVectorArray) *)

  (* fill vector data *)
  Nvector_ops.const neg_two x.(0);
  Nvector_ops.const neg_one y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const neg_one x.(1);
  Nvector_ops.const neg_one y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const two x.(2);
  Nvector_ops.const one y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray one x neg_one y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 4a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 4a";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 4b: z = -x + y, (VDiffVectorArray) *)

  (* fill vector data *)
  Nvector_ops.const two x.(0);
  Nvector_ops.const one y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const neg_one x.(1);
  Nvector_ops.const neg_one y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const neg_two x.(2);
  Nvector_ops.const neg_one y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray neg_one x one y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 4b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 4b";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 5a: z = x + by, (VLin1VectorArray) *)

  (* fill vector data *)
  Nvector_ops.const neg_two x.(0);
  Nvector_ops.const two y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const one x.(1);
  Nvector_ops.const neg_two y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const half x.(2);
  Nvector_ops.const one y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray one x half y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 5a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 5a";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 5b: z = ax + y, (VLin1VectorArray) *)

  (* fill vector data *)
  Nvector_ops.const neg_two x.(0);
  Nvector_ops.const neg_two y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const one x.(1);
  Nvector_ops.const half y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const two x.(2);
  Nvector_ops.const two y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray neg_half x one y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 5b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 5b";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 6a: z = -x + by, (VLin2VectorArray) *)

  (* fill vector data *)
  Nvector_ops.const half x.(0);
  Nvector_ops.const neg_one y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const one x.(1);
  Nvector_ops.const two y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const neg_two x.(2);
  Nvector_ops.const neg_two y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray neg_one x half y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 6a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 6a";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 6b: z = ax - y, (VLin2VectorArray) *)

  (* fill vector data *)
  Nvector_ops.const half x.(0);
  Nvector_ops.const two y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const one x.(1);
  Nvector_ops.const two y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const neg_half x.(2);
  Nvector_ops.const neg_two y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray two x neg_one y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0) local_length
          && check_ans zero    z.(1) local_length
          && check_ans one     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 6b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 6b";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 7: z = a(x + y), (VScaleSumVectorArray) *)

  (* fill vector data *)
  Nvector_ops.const neg_one x.(0);
  Nvector_ops.const half y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const one x.(1);
  Nvector_ops.const half y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const one x.(2);
  Nvector_ops.const neg_half y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray two x two y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 3, +1 *)
  if not (   check_ans neg_one    z.(0) local_length
          && check_ans (two+.one) z.(1) local_length
          && check_ans one        z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 7, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 7";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
             (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 8: z = a(x - y), (VScaleDiffVectorArray) *)

  (* fill vector data *)
  Nvector_ops.const half x.(0);
  Nvector_ops.const one y.(0);
  Nvector_ops.const two z.(0);

  Nvector_ops.const two x.(1);
  Nvector_ops.const half y.(1);
  Nvector_ops.const two z.(0);

  Nvector_ops.const neg_half x.(2);
  Nvector_ops.const neg_one y.(2);
  Nvector_ops.const two z.(0);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray two x neg_two y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of -1, 3, +1 *)
  if not (   check_ans neg_one    z.(0) local_length
          && check_ans (two+.one) z.(1) local_length
          && check_ans one        z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 8, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 8";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 9: z = ax + by, All Other Cases *)

  (* fill vector data *)
  Nvector_ops.const neg_half x.(0);
  Nvector_ops.const two y.(0);

  Nvector_ops.const one x.(1);
  Nvector_ops.const neg_two y.(1);

  Nvector_ops.const half x.(2);
  Nvector_ops.const two y.(2);

  let start_time = get_time () in
  Nvector_ops.linearsumvectorarray two x half y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of 0, +1, +2 *)
  if not (   check_ans zero z.(0) local_length
          && check_ans one  z.(1) local_length
          && check_ans two  z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearSumVectorArray Case 9, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearSumVectorArray Case 9";

  (* find max time across all processes *)
  print_time "N_VLinearSumVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VScaleVectorArray Test
 * --------------------------------------------------------------------*)
let test_scalevectorarray x local_length myid =
  let fails = ref 0 in

  (* create vectors for testing *)
  let c = RealArray.make 3 0.0
  and y = Array.init 3 (fun _ -> Nvector_ops.clone x)
  and z = Array.init 3 (fun _ -> Nvector_ops.clone x) in
  let y_len1, z_len1 = Array.sub y 0 1, Array.sub z 0 1 in

  (* Case 1a: y.(0) = c.{0} y.(0), N_VScale *)

  (* fill vector data *)
  Nvector_ops.const half y.(0);
  c.{0} <- two;

  let start_time = get_time () in
  Nvector_ops.scalevectorarray c y_len1 y_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(0) should be a vector of +1 *)
  if not (check_ans one y.(0) local_length) then (
    printf ">>> FAILED test -- N_VScaleVectorArray Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleVectorArray Case 1a";

  (* find max time across all processes *)
  print_time "N_VScaleVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 1b: z.(0) = c.{0} y.(0), N_VScale *)

  (* fill vector data *)
  Nvector_ops.const half y.(0);
  c.{0} <- two;

  let start_time = get_time () in
  Nvector_ops.scalevectorarray c y_len1 z_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(0) should be a vector of +1 *)
  if not (check_ans one z.(0) local_length) then (
    printf ">>> FAILED test -- N_VScaleVectorArray Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleVectorArray Case 1b";

  (* find max time across all processes *)
  print_time "N_VScaleVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2a: y.(i) = c.{i} y.(i) *)

  (* fill vector data *)
  Nvector_ops.const half y.(0);
  Nvector_ops.const neg_two y.(1);
  Nvector_ops.const neg_one y.(2);

  c.{0} <- two;
  c.{1} <- half;
  c.{2} <- neg_two;

  let start_time = get_time () in
  Nvector_ops.scalevectorarray c y y;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(i) should be a vector of +1, -1, 2 *)
  if not (   check_ans one     y.(0) local_length
          && check_ans neg_one y.(1) local_length
          && check_ans two     y.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleVectorArray Case 2a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleVectorArray Case 2a";

  (* find max time across all processes *)
  print_time "N_VScaleVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2b: z.(i) = c.{i} y.(i) *)

  (* fill vector data *)
  Nvector_ops.const half y.(0);
  Nvector_ops.const neg_two y.(1);
  Nvector_ops.const neg_one y.(2);

  c.{0} <- two;
  c.{1} <- half;
  c.{2} <- neg_two;

  let start_time = get_time () in
  Nvector_ops.scalevectorarray c y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should be a vector of +1, -1, 2 *)
  if not (   check_ans one     z.(0) local_length
          && check_ans neg_one z.(1) local_length
          && check_ans two     z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleVectorArray Case 2b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleVectorArray Case 2b";

  (* find max time across all processes *)
  print_time "N_VScaleVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VConstVectorArray Test
 * --------------------------------------------------------------------*)
let test_constvectorarray x local_length myid =
  let fails = ref 0 in

  (* create vectors for testing *)
  let z = Array.init 3 (fun _ -> Nvector_ops.clone x) in
  let z_len1 = Array.sub z 0 1 in

  (* Case 1a: z.(0) = c, N_VConst *)

  (* fill vector data *)
  Nvector_ops.const zero z.(0);

  let start_time = get_time () in
  Nvector_ops.constvectorarray one z_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(0) should be a vector of 1 *)
  if not (check_ans one z.(0) local_length) then (
    printf ">>> FAILED test -- N_VConstVectorArray Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VConstVectorArray Case 1a";

  (* find max time across all processes *)
  print_time "N_VConstVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 1b: z.(i) = c *)

  (* fill vector data *)
  Nvector_ops.const zero z.(0);
  Nvector_ops.const zero z.(1);
  Nvector_ops.const zero z.(2);

  let start_time = get_time () in
  Nvector_ops.constvectorarray one z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(i) should be a vector of 1 *)
  if not (   check_ans one z.(0) local_length
          && check_ans one z.(1) local_length
          && check_ans one z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VConstVectorArray Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VConstVectorArray Case 1b";

  (* find max time across all processes *)
  print_time "N_VConstVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VWrmsNormVectorArray Test
 * --------------------------------------------------------------------*)
let test_wrmsnormvectorarray x _ myid =
  let fails = ref 0 in

  (* create vectors for testing *)
  let nrm = RealArray.make 3 neg_one
  and w = Array.init 3 (fun _ -> Nvector_ops.clone x)
  and z = Array.init 3 (fun _ -> Nvector_ops.clone x)
  in
  let w_len1, z_len1 = Array.sub w 0 1, Array.sub z 0 1 in

  (* Case 1a: nrm.(0) = ||z.(0)||, N_VWrmsNorm *)

  (* fill vector data *)
  Nvector_ops.const neg_half z.(0);
  Nvector_ops.const half w.(0);

  let start_time = get_time () in
  Nvector_ops.wrmsnormvectorarray z_len1 w_len1 nrm;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* nrm should equal 1/4 *)
  if not (nrm.{0} >= zero && fneq nrm.{0} (half*.half)) then (
    printf ">>> FAILED test -- N_VWrmsNormVectorArray Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VWrmsNormVectorArray Case 1a";

  (* find max time across all processes *)
  print_time "N_VWrmsNormVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 1b: nrm.(i) = ||z.(i)|| *)

  (* fill vector data *)
  Nvector_ops.const neg_half   z.(0);
  Nvector_ops.const (two*.two) z.(1);
  Nvector_ops.const half       z.(2);

  Nvector_ops.const half         w.(0);
  Nvector_ops.const (half*.half) w.(1);
  Nvector_ops.const one          w.(2);

  RealArray.fill nrm neg_one;

  let start_time = get_time () in
  Nvector_ops.wrmsnormvectorarray z w nrm;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal 1/4, 1, 1/2 *)
  if not (   (nrm.{0} >= zero && fneq nrm.{0} (half*.half))
          || (nrm.{1} >= zero && fneq nrm.{1} one)
          || (nrm.{2} >= zero && fneq nrm.{2} half))
  then (
    printf ">>> FAILED test -- N_VWrmsNormVectorArray Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VWrmsNormVectorArray Case 1b";

  (* find max time across all processes *)
  print_time "N_VWrmsNormVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VWrmsNormMaskVectorArray Test
 * --------------------------------------------------------------------*)
let test_wrmsnormmaskvectorarray x local_length global_length myid =
  let fails = ref 0 in
  let xdata = Nvector_ops.getarray x in

  (* factor used in checking solutions *)
  let fac = sqrt (float ((global_length - 1)/global_length)) in

  (* create vectors for testing *)
  let nrm = RealArray.make 3 neg_one
  and w = Array.init 3 (fun _ -> Nvector_ops.clone x)
  and z = Array.init 3 (fun _ -> Nvector_ops.clone x) in
  let w_len1, z_len1 = Array.sub w 0 1, Array.sub z 0 1 in

  (* Case 1: nrm.(0) = ||z.(0)|| *)

  (* fill vector data *)
  Nvector_ops.const neg_half z.(0);
  Nvector_ops.const half w.(0);

  (* use all elements except one *)
  Nvector_ops.const one x;
  if myid = 0 then Nvector_ops.set xdata (local_length - 1) zero;

  let start_time = get_time () in
  Nvector_ops.wrmsnormmaskvectorarray z_len1 w_len1 x nrm;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* nrm should equal fac/4 *)
  if nrm.{0} < zero || fneq nrm.{0} (fac*.half*.half) then (
    printf ">>> FAILED test -- N_VWrmsNormMaskVectorArray Case 1, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VWrmsNormMaskVectorArray Case 1";

  (* find max time across all processes *)
  print_time "N_VWrmsNormVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* Case 2: nrm.(i) = ||z.(i)|| *)

  (* fill vector data *)
  Nvector_ops.const neg_half   z.(0);
  Nvector_ops.const (two*.two) z.(1);
  Nvector_ops.const half       z.(2);

  Nvector_ops.const half         w.(0);
  Nvector_ops.const (half*.half) w.(1);
  Nvector_ops.const one          w.(2);

  (* use all elements except one *)
  Nvector_ops.const one x;
  if myid = 0 then Nvector_ops.set xdata (local_length - 1) zero;

  RealArray.fill nrm neg_one;

  let start_time = get_time () in
  Nvector_ops.wrmsnormmaskvectorarray z w x nrm;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal fac/4, fac, fac/2] *)
  if (   (nrm.{0} < zero || fneq nrm.{0} (fac*.half*.half))
      || (nrm.{1} < zero || fneq nrm.{1} fac)
      || (nrm.{2} < zero || fneq nrm.{2} (fac*.half)))
  then (
    printf ">>> FAILED test -- N_VWrmsNormMaskVectorArray Case 2, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VWrmsNormMaskVectorArray Case 2";

  (* find max time across all processes *)
  print_time "N_VWrmsNormVectorArray"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VScaleAddMultiVectorArray Test
 * --------------------------------------------------------------------*)
let test_scaleaddmultivectorarray v local_length myid =
  let fails = ref 0 in

  (* create vectors for testing *)
  let a = RealArray.make 3 zero
  and x = Array.init 3 (fun _ -> Nvector_ops.clone v)
  and y = Array.init 3 (fun _ -> Array.init 3 (fun _ -> Nvector_ops.clone v))
  and z = Array.init 3 (fun _ -> Array.init 3 (fun _ -> Nvector_ops.clone v))
  in
  let a_len1 = RealArray.sub a 0 1
  and x_len1 = Array.sub x 0 1
  and y_len1_1 = Array.init 1 (fun i -> Array.sub y.(i) 0 1)
  and z_len1_1 = Array.init 1 (fun i -> Array.sub z.(i) 0 1)
  and y_len1_3 = Array.sub y 0 1
  and z_len1_3 = Array.sub z 0 1
  and y_len3_1 = Array.map (fun yi -> Array.sub yi 0 1) y
  and z_len3_1 = Array.map (fun zi -> Array.sub zi 0 1) z
  in

  (* Case 1a (nvec = 1, nsum = 1):
   * z.(0).(0) = a.(0) x.(0) + y.(0).(0), N_VLinearSum *)

  (* fill scaling and vector data *)
  a.{0} <- two;

  Nvector_ops.const one x.(0);
  Nvector_ops.const neg_one y.(0).(0);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a_len1 x_len1 y_len1_1 y_len1_1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(0).(0) should be vector of +1 *)
  if not (check_ans one y.(0).(0) local_length) then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 1a";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 1b (nvec = 1, nsum = 1):
     z.(0).(0) = a.(0) x.(0) + y.(0).(0), N_VLinearSum *)

  (* fill scaling and vector data *)
  a.{0} <- two;

  Nvector_ops.const one     x.(0);
  Nvector_ops.const neg_one y.(0).(0);
  Nvector_ops.const zero    z.(0).(0);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a_len1 x_len1 y_len1_1 z_len1_1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(0).(0) should be vector of +1 *)
  if not (check_ans one z.(0).(0) local_length) then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 1b";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 2a (nvec = 1, nsum > 1):
     y.(j).(0) = a.(j) x.(0) + y.(j).(0), N_VScaleAddMulti *)

  (* fill scaling and vector data *)
  a.{0} <- one;
  a.{1} <- neg_two;
  a.{2} <- two;

  Nvector_ops.const one x.(0);

  Nvector_ops.const neg_two y.(0).(0);
  Nvector_ops.const two     y.(1).(0);
  Nvector_ops.const neg_one y.(2).(0);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a x_len1 y_len3_1 y_len3_1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(i).(0) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one y.(0).(0) local_length
          && check_ans zero    y.(1).(0) local_length
          && check_ans one     y.(2).(0) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 2a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 2a";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 2b (nvec = 1, nsum > 1):
     z.(j).(0) = a.(j) x.(0) + y.(j).(0), N_VScaleAddMulti *)

  (* fill scaling and vector data *)
  a.{0} <- one;
  a.{1} <- neg_two;
  a.{2} <- two;

  Nvector_ops.const one x.(0);

  Nvector_ops.const neg_two y.(0).(0);
  Nvector_ops.const two     y.(1).(0);
  Nvector_ops.const neg_one y.(2).(0);

  Nvector_ops.const zero z.(0).(0);
  Nvector_ops.const one  z.(1).(0);
  Nvector_ops.const two  z.(2).(0);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a x_len1 y_len3_1 z_len3_1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i).(0) should be a vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0).(0) local_length
          && check_ans zero    z.(1).(0) local_length
          && check_ans one     z.(2).(0) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 2b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 2b";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 3a (nvec > 1, nsum = 1):
     y.(0).(i) = a.(0) x.(i) + y.(0).(i), N_VLinearSumVectorArray *)

  (* fill scaling and vector data *)
  a.{0} <- two;

  Nvector_ops.const half    x.(0);
  Nvector_ops.const neg_one x.(1);
  Nvector_ops.const one     x.(2);

  Nvector_ops.const neg_two y.(0).(0);
  Nvector_ops.const two     y.(0).(1);
  Nvector_ops.const neg_one y.(0).(2);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a_len1 x y_len1_3 y_len1_3;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* y.(0).(i) should be vector of -1, 0, +1 *)
  if not (   check_ans neg_one y.(0).(0) local_length
          && check_ans zero    y.(0).(1) local_length
          && check_ans one     y.(0).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 3a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 3a";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 3b (nvec > 1, nsum = 1):
     z.(j).(0) = a.(j) x.(0) + y.(j).(0), N_VLinearSumVectorArray *)

  (* fill scaling and vector data *)
  a.{0} <- two;

  Nvector_ops.const half    x.(0);
  Nvector_ops.const neg_one x.(1);
  Nvector_ops.const one     x.(2);

  Nvector_ops.const neg_two y.(0).(0);
  Nvector_ops.const two     y.(0).(1);
  Nvector_ops.const neg_one y.(0).(2);

  Nvector_ops.const two z.(0).(0);
  Nvector_ops.const two z.(0).(1);
  Nvector_ops.const two z.(0).(2);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a_len1 x y_len1_3 z_len1_3;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(0).(i) should be vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0).(0) local_length
          && check_ans zero    z.(0).(1) local_length
          && check_ans one     z.(0).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 3b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 3b";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 4a (nvec > 1, nsum > 1):
     y.(j).(i) = a.(j) x.(i) + y.(j).(i), N_VScaleAddMultiVectorArray *)

  (* fill scaling and vector data *)
  a.{0} <- two;
  a.{1} <- one;
  a.{2} <- neg_two;

  Nvector_ops.const half     x.(0);
  Nvector_ops.const neg_two  y.(0).(0);
  Nvector_ops.const neg_half y.(1).(0);
  Nvector_ops.const two      y.(2).(0);

  Nvector_ops.const one     x.(1);
  Nvector_ops.const neg_one y.(0).(1);
  Nvector_ops.const neg_two y.(1).(1);
  Nvector_ops.const two     y.(2).(1);

  Nvector_ops.const neg_two        x.(2);
  Nvector_ops.const two            y.(0).(2);
  Nvector_ops.const (two*.two)     y.(1).(2);
  Nvector_ops.const (neg_two*.two) y.(2).(2);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a x y y;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

          (* y.(i).(0) should be vector of -1, 0, +1 *)
  if not (   check_ans neg_one y.(0).(0) local_length
          && check_ans zero    y.(1).(0) local_length
          && check_ans one     y.(2).(0) local_length

          (* y.(i).(1) should be vector of +1, -1, 0 *)
          && check_ans one     y.(0).(1) local_length
          && check_ans neg_one y.(1).(1) local_length
          && check_ans zero    y.(2).(1) local_length

          (* y.(i).(2) should be vector of -2, 2, 0 *)
          && check_ans neg_two y.(0).(2) local_length
          && check_ans two     y.(1).(2) local_length
          && check_ans zero    y.(2).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 4a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 4a";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 4b (nvec > 1, nsum > 1):
     z.(j).(i) = a.(j) x.(i) + y.(j).(i), N_VScaleAddMultiVectorArray *)

  (* fill scaling and vector data *)
  a.{0} <- two;
  a.{1} <- one;
  a.{2} <- neg_two;

  Nvector_ops.const half     x.(0);

  Nvector_ops.const neg_two  y.(0).(0);
  Nvector_ops.const neg_half y.(1).(0);
  Nvector_ops.const two      y.(2).(0);

  Nvector_ops.const half z.(0).(0);
  Nvector_ops.const half z.(1).(0);
  Nvector_ops.const half z.(2).(0);

  Nvector_ops.const one x.(1);

  Nvector_ops.const neg_one y.(0).(1);
  Nvector_ops.const neg_two y.(1).(1);
  Nvector_ops.const two     y.(2).(1);

  Nvector_ops.const half z.(0).(1);
  Nvector_ops.const half z.(1).(1);
  Nvector_ops.const half z.(2).(1);

  Nvector_ops.const neg_two x.(2);

  Nvector_ops.const two            y.(0).(2);
  Nvector_ops.const (two*.two)     y.(1).(2);
  Nvector_ops.const (neg_two*.two) y.(2).(2);

  Nvector_ops.const half z.(0).(2);
  Nvector_ops.const half z.(1).(2);
  Nvector_ops.const half z.(2).(2);

  let start_time = get_time () in
  Nvector_ops.scaleaddmultivectorarray a x y z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

          (* z.(i).(0) should be vector of -1, 0, +1 *)
  if not (   check_ans neg_one z.(0).(0) local_length
          && check_ans zero    z.(1).(0) local_length
          && check_ans one     z.(2).(0) local_length

          (* z.(i).(1) should be vector of +1, -1, 0 *)
          && check_ans one     z.(0).(1) local_length
          && check_ans neg_one z.(1).(1) local_length
          && check_ans zero    z.(2).(1) local_length

          (* z.(i).(2) should be vector of -2, 2, 0 *)
          && check_ans neg_two z.(0).(2) local_length
          && check_ans two     z.(1).(2) local_length
          && check_ans zero    z.(2).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VScaleAddMultiVectorArray Case 4b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VScaleAddMultiVectorArray Case 4b";

  (* find max time across all processes *)
  print_time "N_VScaleAddMultiVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VLinearCombinationVectorArray Test
 * --------------------------------------------------------------------*)
let test_linearcombinationvectorarray v local_length myid =
  let fails = ref 0 in

  (* create vectors for testing *)
  let c = RealArray.make 3 zero
  and z = Array.init 3 (fun _ -> Nvector_ops.clone v)
  and x = Array.init 3 (fun _ -> Array.init 3 (fun _ -> Nvector_ops.clone v))
  in
  let c_len1   = RealArray.sub c 0 1
  and c_len2   = RealArray.sub c 0 2
  and x_len1_1 = Array.init 1 (fun i -> Array.sub x.(i) 0 1)
  and x_len1_3 = Array.sub x 0 1
  and x_len2_1 = Array.init 2 (fun i -> Array.sub x.(i) 0 1)
  and x_len2_3 = Array.init 2 (fun i -> x.(i))
  and x_len3_1 = Array.map (fun xi -> Array.sub xi 0 1) x
  and z_len1   = Array.sub z 0 1
  in

  (* Case 1a: (nvec = 1, nsum = 1), N_VScale
     x.(0).(0) = c.{0} x.(0).(0) *)

  (* fill vector data and scaling factor *)
  Nvector_ops.const half x.(0).(0);
  c.{0} <- two;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len1 x_len1_1 x_len1_1.(0);
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(0) should equal +1 *)
  if not (check_ans one x.(0).(0) local_length) then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 1a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 1a";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 1b: (nvec = 1, nsum = 1), N_VScale
     z.(0) = c.{0} x.(0).(0) *)

  (* fill vector data and scaling factor *)
  Nvector_ops.const half x.(0).(0);
  Nvector_ops.const zero z.(0);
  c.{0} <- two;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len1 x_len1_3 z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(0) should equal +1 *)
  if not (check_ans one z.(0) local_length) then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 1b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 1b";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 2a: (nvec = 1, nsum = 2), N_VLinearSum
     x.(0).(0) = c.{0} x.(0).(0) + c.{1} x.(1).(0) *)

  (* fill vector data and scaling factor *)
  Nvector_ops.const half    x.(0).(0);
  Nvector_ops.const neg_one x.(1).(0);

  c.{0} <- two;
  c.{1} <- neg_one;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len2 x_len2_1 x_len2_1.(0);
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(0) should equal +2 *)
  if not (check_ans two x.(0).(0) local_length) then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 2a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 2a";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 2b: (nvec = 1, nsum = 2), N_VLinearSum
     z.(0) = c.{0} x.(0).(0) + c.{1} x.(1).(0) *)

  (* fill vector data and scaling factor *)
  Nvector_ops.const half    x.(0).(0);
  Nvector_ops.const neg_one x.(1).(0);

  c.{0} <- two;
  c.{1} <- neg_one;

  Nvector_ops.const zero z.(0);

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len2 x_len2_1 z_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(0) should equal +2 *)
  if not (check_ans two z.(0) local_length) then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 2b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 2b";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
             (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 3a: (nvec = 1, nsum > 2), N_VLinearCombination
     x.(0).(0) = c.{0} x.(0).(0) + c.{1} x.(1).(0) + c.{2} x.(2).(0) *)

  (* fill vector data *)
  Nvector_ops.const one x.(0).(0);
  Nvector_ops.const neg_two x.(1).(0);
  Nvector_ops.const neg_one x.(2).(0);

  (* set scaling factors *)
  c.{0} <- two;
  c.{1} <- half;
  c.{2} <- neg_one;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c x_len3_1 x_len3_1.(0);
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(0) should equal +2 *)
  if not (check_ans two x.(0).(0) local_length) then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 3a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 3a";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
             (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 3b: (nvec = 1, nsum > 2), N_VLinearCombination
     z.(0) = c.{0} x.(0).(0) + c.{1} x.(1).(0) + c.{2} x.(2).(0) *)

  (* fill vector data *)
  Nvector_ops.const one x.(0).(0);
  Nvector_ops.const neg_two x.(1).(0);
  Nvector_ops.const neg_one x.(2).(0);

  (* set scaling factors *)
  c.{0} <- two;
  c.{1} <- half;
  c.{2} <- neg_one;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c x_len3_1 z_len1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(0) should equal +2 *)
  if not (check_ans two z.(0) local_length) then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 3b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 3b";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 4a: (nvec > 1, nsum = 1), N_VScaleVectorArray
     x.(0).(i) = c.{0} x.(0).(i) *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const neg_two x.(0).(0);
  Nvector_ops.const neg_one x.(0).(1);
  Nvector_ops.const two     x.(0).(2);

  c.{0} <- half;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len1 x_len1_3 x.(0);
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(i) should equal to -1, -1/2, +1 *)
  if not (   check_ans neg_one  x.(0).(0) local_length
          && check_ans neg_half x.(0).(1) local_length
          && check_ans one      x.(0).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 4a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 4a";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 4b: (nvec > 1, nsum = 1), N_VScaleVectorArray
     z.(i) = c.{0} x.(0).(i) *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const neg_two x.(0).(0);
  Nvector_ops.const neg_one x.(0).(1);
  Nvector_ops.const two     x.(0).(2);

  c.{0} <- half;

  Nvector_ops.const zero z.(0);
  Nvector_ops.const zero z.(1);
  Nvector_ops.const zero z.(2);

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len1 x_len1_3 z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(i) should equal to -1, -1/2, +1 *)
  if not (   check_ans neg_one  z.(0) local_length
          && check_ans neg_half z.(1) local_length
          && check_ans one      z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 4b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 4b";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 5a: (nvec > 1, nsum = 2), N_VLinearSumVectorArray
     x.(0).(i) = c.{0} x.(0).(i) + c.{1} x.(1).(i) *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const neg_two x.(0).(0);
  Nvector_ops.const two x.(1).(0);

  Nvector_ops.const two x.(0).(1);
  Nvector_ops.const half x.(1).(1);

  Nvector_ops.const zero x.(0).(2);
  Nvector_ops.const half x.(1).(2);

  c.{0} <- half;
  c.{1} <- two;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len2 x_len2_3 x.(0);
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(i) should equal to +3, +2, +1 *)
  if not (   check_ans (one+.two) x.(0).(0) local_length
          && check_ans two        x.(0).(1) local_length
          && check_ans one        x.(0).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 5a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 5a";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 5b: (nvec > 1, nsum = 2), N_VLinearSumVectorArray
     z.(0) = c.{0} x.(0).(i) + c.{1} x.(1).(i) *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const neg_two x.(0).(0);
  Nvector_ops.const two     x.(1).(0);

  Nvector_ops.const two  x.(0).(1);
  Nvector_ops.const half x.(1).(1);

  Nvector_ops.const zero x.(0).(2);
  Nvector_ops.const half x.(1).(2);

  c.{0} <- half;
  c.{1} <- two;

  Nvector_ops.const zero z.(0);
  Nvector_ops.const zero z.(1);
  Nvector_ops.const zero z.(2);

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c_len2 x_len2_3 z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(i) should equal to +3, +2, +1 *)
  if not (   check_ans (one+.two) z.(0) local_length
          && check_ans two        z.(1) local_length
          && check_ans one        z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 5b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 5b";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 6a: (nvec > 1, nsum > 2)
     x.(0).(i) += c.{1} x.(1).(i) + c.{2} x.(2).(i) *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const two     x.(0).(0);
  Nvector_ops.const neg_two x.(1).(0);
  Nvector_ops.const neg_one x.(2).(0);

  Nvector_ops.const one x.(0).(1);
  Nvector_ops.const two x.(1).(1);
  Nvector_ops.const one x.(2).(1);

  Nvector_ops.const neg_one x.(0).(2);
  Nvector_ops.const two     x.(1).(2);
  Nvector_ops.const two     x.(2).(2);

  c.{0} <- one;
  c.{1} <- neg_half;
  c.{2} <- neg_one;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c x x.(0);
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(i) should equal to +4, -1, -4 *)
  if not (   check_ans (two+.two)   x.(0).(0) local_length
          && check_ans neg_one      x.(0).(1) local_length
          && check_ans (-.two-.two) x.(0).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 6a, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 6a";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 6b: (nvec > 1, nsum > 2)
     x.(0).(i) = c.{0} x.(0).(i) + c.{1} x.(1).(i) + c.{2} x.(2).(i) *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const one x.(0).(0);
  Nvector_ops.const neg_two x.(1).(0);
  Nvector_ops.const neg_one x.(2).(0);

  Nvector_ops.const neg_one x.(0).(1);
  Nvector_ops.const two x.(1).(1);
  Nvector_ops.const one x.(2).(1);

  Nvector_ops.const half x.(0).(2);
  Nvector_ops.const two x.(1).(2);
  Nvector_ops.const one x.(2).(2);

  c.{0} <- two;
  c.{1} <- half;
  c.{2} <- neg_one;

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c x x.(0);
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* x.(0).(i) should equal to +2, -2, +1 *)
  if not (   check_ans two      x.(0).(0) local_length
          && check_ans neg_two  x.(0).(1) local_length
          && check_ans one      x.(0).(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 6b, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 6b";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  (* Case 6c: (nvec > 1, nsum > 2)
     z.(i) = c.{0} x.(0).(i) + c.{1} x.(1).(i) + c.{2} x.(2).(i) *)

  (* fill vector data and set scaling factors *)
  Nvector_ops.const one     x.(0).(0);
  Nvector_ops.const neg_two x.(1).(0);
  Nvector_ops.const neg_one x.(2).(0);

  Nvector_ops.const neg_one x.(0).(1);
  Nvector_ops.const two     x.(1).(1);
  Nvector_ops.const one     x.(2).(1);

  Nvector_ops.const half x.(0).(2);
  Nvector_ops.const two  x.(1).(2);
  Nvector_ops.const one  x.(2).(2);

  c.{0} <- two;
  c.{1} <- half;
  c.{2} <- neg_one;

  Nvector_ops.const zero z.(0);
  Nvector_ops.const zero z.(1);
  Nvector_ops.const zero z.(2);

  let start_time = get_time () in
  Nvector_ops.linearcombinationvectorarray c x z;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* z.(i) should equal to +2, -2, +1 *)
  if not (   check_ans two      z.(0) local_length
          && check_ans neg_two  z.(1) local_length
          && check_ans one      z.(2) local_length)
  then (
    printf ">>> FAILED test -- N_VLinearCombinationVectorArray Case 6c, Proc %d \n" myid;
    fails += 1
  ) else if myid = 0 then print_passed "N_VLinearCombinationVectorArray Case 6c";

  (* find max time across all processes *)
  print_time "N_VLinearCombinationVectorArray"
    (Nvector_ops.max_time v (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VDotProdLocal test
 * --------------------------------------------------------------------*)
let test_dotprodlocal x y local_length myid =
  let fails = ref 0 in
  (* fill vector data *)
  let rmyid = float_of_int myid in
  let locleninv = 1.0 /. float_of_int local_length in
  Nvector_ops.const rmyid x;
  Nvector_ops.const locleninv y;

  let start_time = get_time () in
  let ans = Nvector_ops.Local.dotprod x y in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal rmyid *)
  if not (fneqtol ans rmyid (sqrt (Config.unit_roundoff)))
  then (printf ">>> FAILED test -- N_VDotProdLocal, Proc %d\n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VDotProdLocal\n";

  (* find max time across all processes *)
  print_time "N_VDotProdLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VMaxNormLocal test
 * --------------------------------------------------------------------*)
let test_maxnormlocal x local_length myid =
  let fails = ref 0 in
  (* fill vector data *)
  let xdata = Nvector_ops.getarray x in
  let myidp1 = float_of_int (myid+1) in
  Nvector_ops.const (-0.5) x;
  Nvector_ops.set xdata (local_length - 1) (-. myidp1);

  let start_time = get_time () in
  let ans = Nvector_ops.Local.maxnorm x in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal myidp1 *)
  if ans < 0.0 || not (fneqtol ans myidp1 (sqrt Config.unit_roundoff))
  then (printf ">>> FAILED test -- N_VMaxNormLocal, Proc %d\n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VMaxNormLocal\n";

  (* find max time across all processes *)
  print_time "N_VMaxNormLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails


(* ----------------------------------------------------------------------
 * N_VMinLocal test
 * --------------------------------------------------------------------*)
let test_minlocal x local_length myid =
  let fails = ref 0 in
  (* fill vector data *)
  let xdata = Nvector_ops.getarray x in
  let negmyid = float_of_int (-myid) in
  Nvector_ops.const 2.0 x;
  Nvector_ops.set xdata (local_length - 1) negmyid;

  let start_time = get_time () in
  let ans = Nvector_ops.Local.min x in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal negmyid *)
  if not (fneqtol ans negmyid (sqrt Config.unit_roundoff))
  then (printf ">>> FAILED test -- N_VMinLocal, Proc %d\n" myid; fails +=1)
  else if myid = 0 then printf "PASSED test -- N_VMinLocal\n";

  (* find max time across all processes *)
  print_time "N_VMinLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VL1NormLocal test
 * --------------------------------------------------------------------*)
let test_l1normlocal x local_length myid =
  let fails = ref 0 in
  (* fill vector data *)
  let v = -. (float_of_int myid) /. (float_of_int local_length) in
  Nvector_ops.const v x;

  let start_time = get_time () in
  let ans = Nvector_ops.Local.l1norm x in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal myid *)
  if ans < 0.0 || not (fneqtol ans (float_of_int myid) (sqrt Config.unit_roundoff))
  then (printf ">>> FAILED test -- N_VL1NormLocal, Proc %d\n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VL1NormLocal\n";

  (* find max time across all processes *)
  print_time "N_VL1NormLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VWSqrSumLocal test
 * --------------------------------------------------------------------*)
let test_wsqrsumlocal x w local_length myid =
  let fails = ref 0 in
  (* fill vector data *)
  let xval = sqrt (float_of_int myid) in
  let wval = 1.0 /. sqrt (float_of_int local_length) in
  Nvector_ops.const xval x;
  Nvector_ops.const wval w;

  let start_time = get_time () in
  let ans = Nvector_ops.Local.wsqrsum x w in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal myid *)
  if ans < 0.0 || not (fneqtol ans (float_of_int myid) (sqrt Config.unit_roundoff))
  then (printf ">>> FAILED test -- N_VWSqrSumLocal, Proc %d\n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VWSqrSumLocal\n";

  (* find max time across all processes *)
  print_time "N_VWSqrSumLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VWSqrSumMaskLocal test
 * --------------------------------------------------------------------*)
let test_wsqrsummasklocal x w id local_length myid =
  let fails = ref 0 in
  (* fill vector data *)
  let xval = sqrt (float_of_int myid) in
  let wval = 1.0 /. sqrt (float_of_int (local_length - 1)) in
  let id_data = Nvector_ops.getarray id in
  Nvector_ops.const xval x;
  Nvector_ops.const wval w;
  (* use all elements except one *)
  Nvector_ops.const 1.0 id;
  Nvector_ops.set id_data (local_length - 1) 0.0;

  let start_time = get_time () in
  let ans = Nvector_ops.Local.wsqrsummask x w id in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal myid *)
  if ans < 0.0 || not (fneqtol ans (float_of_int myid) (sqrt Config.unit_roundoff))
  then (printf ">>> FAILED test -- N_VWSqrSumMaskLocal, Proc %d\n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VWSqrSumMaskLocal\n";

  (* find max time across all processes *)
  print_time "N_VWSqrSumMaskLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VInvTestLocal test
 * --------------------------------------------------------------------*)
let test_invtestlocal x z local_length myid =
  let fails = ref 0 in

  if local_length < 2 then
    (printf "Error Test_N_VInvTestLocal: Local vector length is %d, length must be >= 2\n"
            local_length;
     raise (TestFailed (-1)));

  (*
   * Case 1: All elements Nonzero, z[i] = 1/x[i], return True
   *)

  (* fill vector data *)
  let xval = 1.0 /. float_of_int (myid + 2) in
  Nvector_ops.const xval x;
  Nvector_ops.const 0.0 z;

  let start_time = get_time () in
  let test = Nvector_ops.Local.invtest x z in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* Z should be vector of myid+2 *)
  if not (check_ans (float_of_int (myid + 2)) z local_length && test)
  then (printf ">>> FAILED test -- N_VInvTestLocal Case 1, Proc %d\n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VInvTestLocal Case 1\n";

  (* find max time across all processes *)
  print_time "N_VInvTestLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (*
   * Case 2: Some elements Zero, z[i] = 1/x[i] for x[i] != 0, return False
   *)

  (* reset failure *)

  (* fill vector data *)
  Nvector_ops.const 0.0 z;
  let xdata = Nvector_ops.getarray x in
  for i = 0 to local_length - 1 do
    Nvector_ops.set xdata i (if i mod 2 = 1 then 0.5 else 0.0)
  done;

  let start_time = get_time () in
  let test = Nvector_ops.Local.invtest x z in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* check return vector *)
  let zdata = Nvector_ops.getarray z in
  let failure =
    try
      for i = 0 to local_length - 1 do
        if Nvector_ops.get zdata i <> (if i mod 2 = 1 then 2.0 else 0.0)
        then raise Exit
      done;
      false
    with Exit -> true
  in

  if failure || test
  then (printf ">>> FAILED test -- N_VInvTestLocal Case 2, Proc %d\n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VInvTestLocal Case 2\n";

  (* find max time across all processes *)
  print_time "N_VInvTestLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VConstrMaskLocal test
 * --------------------------------------------------------------------*)
let test_constrmasklocal c x m local_length myid =
  let fails = ref 0 in
  if local_length < 7
  then (printf "Error Test_N_VConstrMaskLocal: Local vector length is %d, length must be >= 7\n" local_length;
        raise (TestFailed (-1)));

  (*
   * Case 1: Return True
   *)

  (* fill vector data *)
  let xdata = Nvector_ops.getarray x in
  let cdata = Nvector_ops.getarray c in
  Nvector_ops.const (-1.0) m;
  for i = 0 to local_length - 1 do
    let vc, vx = match i mod 7 with
      | 0 -> -2.0, -2.0 (* c = -2, test for < 0*)
      | 1 -> -1.0, -1.0 (* c = -1, test for <= 0 *)
      | 2 -> -1.0,  0.0 (* c = -1, test for == 0 *)
      | 3 ->  0.0,  0.5 (* c = 0, no test *)
      | 4 ->  1.0,  0.0 (* c = 1, test for == 0*)
      | 5 ->  1.0,  1.0 (* c = 1, test for >= 0*)
      | 6 ->  2.0,  2.0 (* c = 2, test for > 0 *)
      | _ -> assert false
    in
    Nvector_ops.set cdata i vc;
    Nvector_ops.set xdata i vx
  done;

  let start_time = get_time () in
  let test = Nvector_ops.Local.constrmask c x m in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* M should be vector of 0 *)
  if not (check_ans 0.0 m local_length && test)
  then (printf ">>> FAILED test -- N_VConstrMaskLocal Case 1, Proc %d \n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VConstrMaskLocal Case 1 \n";

  (* find max time across all processes *)
  print_time "N_VConstrMaskLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (*
   * Case 2: Return False
   *)

  (* fill vector data *)
  Nvector_ops.const (-1.0) m;
  for i = 0 to local_length - 1 do
    let vc, vx = match i mod 5 with
      | 0 -> -2.0,  2.0 (* c = -2, test for < 0*)
      | 1 -> -1.0,  1.0 (* c = -1, test for <= 0 *)
      | 2 ->  0.0,  0.5 (* c = 0, no test *)
      | 3 ->  1.0, -1.0 (* c = 1, test for >= 0*)
      | 4 ->  2.0, -2.0 (* c = 2, test for > 0 *)
      | _ -> assert false
    in
    Nvector_ops.set cdata i vc;
    Nvector_ops.set xdata i vx
  done;

  let start_time = get_time () in
  let test = Nvector_ops.Local.constrmask c x m in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* check mask vector *)
  let mdata = Nvector_ops.getarray m in
  let failure =
    try
      for i = 0 to local_length - 1 do
        if Nvector_ops.get mdata i <> (if i mod 5 = 2 then 0.0 else 1.0)
        then raise Exit
      done;
      false
    with Exit -> true
  in

  if failure || test
  then (printf ">>> FAILED test -- N_VConstrMaskLocal Case 2, Proc %d \n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VConstrMaskLocal Case 2 \n";

  (* find max time across all processes *)
  print_time "N_VConstrMaskLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VMinQuotientLocal test
 * --------------------------------------------------------------------*)
let test_minquotientlocal num denom _ myid =
  let fails = ref 0 in

  (*
   * Case 1: Pass
   *)

  (* fill vector data *)
  Nvector_ops.const 2.0 denom;
  Nvector_ops.const (2.0 *. float_of_int myid) num;

  let start_time = get_time () in
  let ans = Nvector_ops.Local.minquotient num denom in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal myid *)
  if not (fneqtol ans (float_of_int myid) (sqrt Config.unit_roundoff))
  then (printf ">>> FAILED test -- N_VMinQuotientLocal Case 1, Proc %d \n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VMinQuotientLocal Case 1 \n";

  (* find max time across all processes *)
  print_time "N_VMinQuotientLocal"
    (Nvector_ops.max_time num (stop_time -. start_time));

  (*
   * Case 2: Fail
   *)

  (* fill vector data *)
  Nvector_ops.const 2.0 num;
  Nvector_ops.const 0.0 denom;

  let start_time = get_time () in
  let ans = Nvector_ops.Local.minquotient num denom in
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* ans should equal BIG_REAL *)
  if not (fneqtol ans Config.big_real (sqrt Config.unit_roundoff))
  then (printf ">>> FAILED test -- N_VMinQuotientLocal Case 2, Proc %d \n" myid; fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VMinQuotientLocal Case 2 \n";

  (* find max time across all processes *)
  print_time "N_VMinQuotientLocal"
    (Nvector_ops.max_time num (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VDotProdMultiLocal Test
 * --------------------------------------------------------------------*)
let test_dotprodmultilocal x local_length myid =
  let fails = ref 0 in
  let dotprods = RealArray.make 3 0.0 in
  let dotprods1 = RealArray.make 1 0.0 in

  (* create vectors for testing *)
  let v = Array.init 3 (fun _ -> Nvector_ops.clone x) in
  let v1 = Array.sub v 0 1 in

  (*
   * Case 1: d[0] = z . V[0], N_VDotProd
   *)

  (* fill vector data *)
  Nvector_ops.const 2.0 x;
  Nvector_ops.const 0.5 v1.(0);

  let start_time = get_time () in
  Nvector_ops.Local.dotprodmulti x v1 dotprods1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod[0] should equal the local vector length *)
  if not (fneq dotprods1.{0} (float local_length))
  then (printf ">>> FAILED test -- N_VDotProdMultiLocal Case 1, Proc %d \n" myid;
        fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VDotProdMultiLocal Case 1 \n";

  (* find max time across all processes *)
  print_time "N_VDotProdMultiLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (*
   * Case 2: d[i] = z . V[i], N_VDotProd
   *)

  (* fill vector data *)
  Nvector_ops.const 2.0 x;
  Nvector_ops.const (-0.5) v.(0);
  Nvector_ops.const 0.5 v.(1);
  Nvector_ops.const 1.0 v.(2);

  let start_time = get_time () in
  Nvector_ops.Local.dotprodmulti x v dotprods;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod[i] should equal -1, +1, and 2 times the local vector length *)
  if ((if fneq dotprods.{0} (float (-1 * local_length)) then false
                                                        else (fails += 1; true))
     || (if fneq dotprods.{1} (float local_length) then false
                                                   else (fails += 1; true))

     || (if fneq dotprods.{2} (float (2 * local_length)) then false
                                                         else (fails += 1; true)))
  then printf ">>> FAILED test -- N_VDotProdMultiLocal Case 2, Proc %d \n" myid
  else if myid = 0 then printf "PASSED test -- N_VDotProdMultiLocal Case 2 \n";

  (* find max time across all processes *)
  print_time "N_VDotProdMultiLocal"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VDotProdMultiAllReduce Test
 * --------------------------------------------------------------------*)
let test_dotprodmultiallreduce x local_length myid =
  let fails = ref 0 in
  let dotprods = RealArray.make 3 0.0 in
  let dotprods1 = RealArray.make 1 0.0 in

  (* get global length *)
  let global_length = Nvector_ops.getlength x in

  (* create vectors for testing *)
  let v = Array.init 3 (fun _ -> Nvector_ops.clone x) in
  let v1 = Array.sub v 0 1 in

  (*
   * Case 1: d[0] = z . V[0], N_VDotProd
   *)

  (* fill vector data *)
  Nvector_ops.const 2.0 x;
  Nvector_ops.const 0.5 v1.(0);

  let start_time = get_time () in
  Nvector_ops.Local.dotprodmulti x v1 dotprods1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod[0] should equal the local vector length *)
  if not (fneq dotprods1.{0} (float local_length))
  then (printf ">>> FAILED test -- N_VDotProdMultiAllReduce Case 1, Proc %d \n" myid;
        fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VDotProdMultiAllReduce Case 1 \n";

  (* find max time across all processes *)
  print_time "N_VDotProdMultiAllReduce"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* perform the global reduction *)
  let start_time = get_time () in
  Nvector_ops.Local.dotprodmulti_allreduce x dotprods1;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod[0] should equal the global vector length *)
  if not (fneq dotprods1.{0} (float global_length))
  then (printf ">>> FAILED test -- N_VDotProdMultiAllReduce Case 1, Proc %d \n" myid;
        fails += 1)
  else if myid = 0 then printf "PASSED test -- N_VDotProdMultiAllReduce Case 1 \n";

  (* find max time across all processes *)
  print_time "N_VDotProdMultiAllReduce"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (*
   * Case 2: d[i] = z . V[i], N_VDotProd
   *)

  (* fill vector data *)
  Nvector_ops.const 2.0 x;
  Nvector_ops.const (-0.5) v.(0);
  Nvector_ops.const 0.5 v.(1);
  Nvector_ops.const 1.0 v.(2);

  let start_time = get_time () in
  Nvector_ops.Local.dotprodmulti x v dotprods;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod[i] should equal -1, +1, and 2 times the local vector length *)
  if ((if fneq dotprods.{0} (float (-1 * local_length)) then false
                                                        else (fails += 1; true))
     || (if fneq dotprods.{1} (float local_length) then false
                                                   else (fails += 1; true))

     || (if fneq dotprods.{2} (float (2 * local_length)) then false
                                                         else (fails += 1; true)))
  then printf ">>> FAILED test -- N_VDotProdMultiLocal Case 2, Proc %d \n" myid
  else if myid = 0 then printf "PASSED test -- N_VDotProdMultiLocal Case 2 \n";

  (* find max time across all processes *)
  print_time "N_VDotProdMultiAllReduce"
    (Nvector_ops.max_time x (stop_time -. start_time));

  (* perform the global reduction *)
  let start_time = get_time () in
  Nvector_ops.Local.dotprodmulti_allreduce x dotprods;
  Nvector_ops.sync_device ();
  let stop_time = get_time () in

  (* dotprod[i] should equal -1, +1, and 2 times the global vector length *)
  if ((if fneq dotprods.{0} (float (-1 * global_length)) then false
                                                        else (fails += 1; true))
     || (if fneq dotprods.{1} (float global_length) then false
                                                   else (fails += 1; true))

     || (if fneq dotprods.{2} (float (2 * global_length)) then false
                                                         else (fails += 1; true)))
  then printf ">>> FAILED test -- N_VDotProdMultiAllReduce Case 2, Proc %d \n" myid
  else if myid = 0 then printf "PASSED test -- N_VDotProdMultiAllReduce Case 2 \n";

  (* find max time across all processes *)
  print_time "N_VDotProdMultiAllReduce"
    (Nvector_ops.max_time x (stop_time -. start_time));

  !fails

(* ----------------------------------------------------------------------
 * N_VBufSize test
 * --------------------------------------------------------------------*)
let test_bufsize _ _ myid =
  (* Not implemented in Sundials/ML *)
  if myid = 0 then printf "PASSED test -- N_VBufSize\n"; 0

(* ----------------------------------------------------------------------
 * N_VBufPack test
 * --------------------------------------------------------------------*)
let test_bufpack _ _ myid =
  (* Not implemented in Sundials/ML *)
  if myid = 0 then printf "PASSED test -- N_VBufPack\n"; 0

(* ----------------------------------------------------------------------
 * N_VBufUnpack test
 * --------------------------------------------------------------------*)
let test_bufunpack _ _ myid =
  (* Not implemented in Sundials/ML *)
  if myid = 0 then printf "PASSED test -- N_VBufUnpack\n"; 0

end

