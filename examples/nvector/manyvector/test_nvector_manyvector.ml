(* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Sep 2021.
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the NVECTOR ManyVector
 * (serial) module implementation.
 * -----------------------------------------------------------------*)

module ROArray = Sundials.ROArray
module RealArray = Sundials.RealArray

module Nvector_manyvector_ops =
  struct
    include Nvector_many.Ops
    let get_id = Nvector.get_id
    type contents = Nvector.any ROArray.t * int
    let getarray = Nvector_many.unwrap

    let get (xxs, _) i =
      let xsub0 = Nvector_serial.Any.unwrap (ROArray.get xxs 0) in
      let xsub1 = Nvector_serial.Any.unwrap (ROArray.get xxs 1) in
      let x0len = RealArray.length xsub0 in
      if i < x0len then Sundials.RealArray.get xsub0 i
      else Sundials.RealArray.get xsub1 (i - x0len)

    let set (xxs, _) i v =
      let xsub0 = Nvector_serial.Any.unwrap (ROArray.get xxs 0) in
      let xsub1 = Nvector_serial.Any.unwrap (ROArray.get xxs 1) in
      let x0len = RealArray.length xsub0 in
      if i < x0len then RealArray.set xsub0 i v
      else RealArray.set xsub1 (i - x0len) v

    let max_time _ t = t
    let sync_device () = ()
  end

module Test =
  struct
  include Test_nvector.Test (Nvector_manyvector_ops)

  let make lens =
    let num = Array.length lens in
    Nvector_many.wrap
      (ROArray.init num (fun i -> Nvector_serial.Any.make lens.(i) 0.0))

  let id = Nvector.ManyVector
end

let printf = Printf.printf
let (+=) r x = r := !r + x

(* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in

  (* check input and set vector length *)
  if Array.length Sys.argv < 4 then
    (printf "ERROR: THREE (3) Inputs required: subvector 1 length, subvector 2 length, print timing \n";
     exit (-1))
  ;

  let len1 = int_of_string Sys.argv.(1) in
  if len1 <= 0 then
    (printf "ERROR: length of subvector 1 must be a positive integer \n";
     exit (-1))
  ;

  let len2 = int_of_string Sys.argv.(2) in
  if len2 <= 0 then
    (printf "ERROR: length of subvector 2 must be a positive integer \n";
     exit (-1))
  ;

  let length = len1 + len2 in

  let print_timing = int_of_string Sys.argv.(3) in
  let _ = Test.set_timing (print_timing <> 0) true in

  printf "Testing ManyVector (serial) N_Vector \n";
  printf "Vector lengths: %d %d \n" len1 len2;

  let x = Test.make [| len1; len2 |] in

  (* NVector Tests *)
  fails += Test.test_getvectorid x Test.id 0;
  fails += Test.test_getlength x 0;
  fails += Test.test_getcommunicator x 0 0;

  (* Test subvector accessors *)
  let xss, _ = Nvector_many.unwrap x in
  if ROArray.length xss <> 2 then (
    printf ">>> FAILED test -- N_VGetNumSubvectors_ManyVector\n";
    fails += 1);
  let x1 = Nvector_serial.Any.unwrap (ROArray.get xss 0) in
  if RealArray.length x1 <> len1 then (
    printf ">>> FAILED test -- N_VGetSubvector_ManyVector\n";
    fails += 1);
  let x2 = Nvector_serial.Any.unwrap (ROArray.get xss 1) in
  if RealArray.length x2 <> len2 then (
    printf ">>> FAILED test -- N_VGetSubvector_ManyVector\n";
    fails += 1);

  let y = Nvector_many.Ops.clone x in
  let z = Nvector_many.Ops.clone x in

  (* Standard vector operation tests *)
  printf "\nTesting standard vector operations:\n\n";

  fails += Test.test_const x length 0;
  fails += Test.test_linearsum x y z length 0;
  fails += Test.test_prod x y z length 0;
  fails += Test.test_div x y z length 0;
  fails += Test.test_scale x z length 0;
  fails += Test.test_abs x z length 0;
  fails += Test.test_inv x z length 0;
  fails += Test.test_addconst x z length 0;
  fails += Test.test_dotprod x y length length 0;
  fails += Test.test_maxnorm x length 0;
  fails += Test.test_wrmsnorm x y length 0;
  fails += Test.test_wrmsnormmask x y z length length 0;
  fails += Test.test_min x length 0;
  fails += Test.test_wl2norm x y length length 0;
  fails += Test.test_l1norm x length length 0;
  fails += Test.test_compare x z length 0;
  fails += Test.test_invtest x z length 0;
  fails += Test.test_constrmask x y z length 0;
  fails += Test.test_minquotient x y length 0;

  (* Fused and vector array operations tests (disabled) *)
  printf "\nTesting fused and vector array operations (disabled):\n\n";

  let u = Nvector_many.Ops.clone x in
  Nvector_many.enable ~with_fused_ops:false u;

  (* fused operations *)
  fails += Test.test_linearcombination u length 0;
  fails += Test.test_scaleaddmulti u length 0;
  fails += Test.test_dotprodmulti u length length 0;

  (* vector array operations *)
  fails += Test.test_linearsumvectorarray u length 0;
  fails += Test.test_scalevectorarray u length 0;
  fails += Test.test_constvectorarray u length 0;
  fails += Test.test_wrmsnormvectorarray u length 0;
  fails += Test.test_wrmsnormmaskvectorarray u length length 0;
  fails += Test.test_scaleaddmultivectorarray u length 0;
  fails += Test.test_linearcombinationvectorarray u length 0;

  (* Fused and vector array operations tests (enabled) *)
  printf "\nTesting fused and vector array operations (enabled):\n\n";

  let u = Nvector_many.Ops.clone x in
  Nvector_many.enable ~with_fused_ops:true u;

  (* fused operations *)
  fails += Test.test_linearcombination u length 0;
  fails += Test.test_scaleaddmulti u length 0;
  fails += Test.test_dotprodmulti u length length 0;

  (* vector array operations *)
  fails += Test.test_linearsumvectorarray u length 0;
  fails += Test.test_scalevectorarray u length 0;
  fails += Test.test_constvectorarray u length 0;
  fails += Test.test_wrmsnormvectorarray u length 0;
  fails += Test.test_wrmsnormmaskvectorarray u length length 0;
  fails += Test.test_scaleaddmultivectorarray u length 0;
  fails += Test.test_linearcombinationvectorarray u length 0;

  printf "\nTesting local reduction operations:\n\n";

  fails += Test.test_dotprodlocal x y length 0;
  fails += Test.test_maxnormlocal x length 0;
  fails += Test.test_minlocal x length 0;
  fails += Test.test_l1normlocal x length 0;
  fails += Test.test_wsqrsumlocal x y length 0;
  fails += Test.test_wsqrsummasklocal x y z length 0;
  fails += Test.test_invtestlocal x z length 0;
  fails += Test.test_constrmasklocal x y z length 0;
  fails += Test.test_minquotientlocal x y length 0;

  (* XBraid interface operations *)
  printf "\nTesting XBraid interface operations:\n\n";

  fails += Test.test_bufsize x length 0;
  fails += Test.test_bufpack x length 0;
  fails += Test.test_bufunpack x length 0;

  (* Print result *)
  if !fails <> 0 then printf "FAIL: NVector module failed %d tests \n\n" !fails
  else printf "SUCCESS: NVector module passed all tests \n\n"

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
