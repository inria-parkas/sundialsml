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
 * This is the testing routine to check the NVECTOR MPIPlusX
 * module implementation.
 * -----------------------------------------------------------------*)

module ROArray = Sundials.ROArray
module RealArray = Sundials.RealArray

module Nvector_mpiplusx_ops =
  struct
    include Nvector_mpiplusx.Ops
    let get_id = Nvector.get_id
    type contents = Nvector.any * Mpi.communicator
    let getarray = Nvector.unwrap

    let get (xs, _) i =
      let x = Nvector_serial.Any.unwrap xs in
      Sundials.RealArray.get x i

    let set (xs, _) i v =
      let x = Nvector_serial.Any.unwrap xs in
      Sundials.RealArray.set x i v

    let max_time x t =
      let _, comm = Nvector_mpiplusx.unwrap x in
      (* get max time across all MPI ranks *)
      Mpi.(reduce_float t Max 0 comm)

    let sync_device () = ()
  end

module Test =
  struct
  include Test_nvector.Test (Nvector_mpiplusx_ops)
  let id = Nvector.MpiPlusX
end

let printf = Printf.printf
let (+=) r x = r := !r + x

(* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in
  let comm = Mpi.comm_world in
  let nprocs = Mpi.comm_size comm in
  let myid = Mpi.comm_rank comm in

  (* check input and set vector length *)
  if Array.length Sys.argv < 3 then
    (printf "ERROR: TWO (2) Inputs required: vector local length, print timing \n";
     exit (-1))
  ;

  let local_length = int_of_string Sys.argv.(1) in
  if local_length <= 0 then
    (printf "ERROR: local vector length must be a positive integer \n";
     exit (-1))
  ;

  let global_length = nprocs * local_length in

  let print_timing = int_of_string Sys.argv.(2) in
  let _ = Test.set_timing (print_timing <> 0)
                          (Test_nvector.compat_neq600 && myid = 0)
  in

  if myid = 0 then begin
    printf "Testing the MPIPlusX N_Vector with X being a serial N_Vector\n";
    printf "Vector local length %d\n" local_length;
    printf "MPIPlusX vector global length %d\n" global_length;
    printf "MPI processes %d\n" nprocs
  end;

  let xlocal = Nvector_serial.Any.make local_length 0.0 in
  let x = Nvector_mpiplusx.wrap comm xlocal in

  (* NVector Tests *)
  fails += Test.test_getvectorid x Test.id myid;
  fails += Test.test_getlength x myid;
  fails += Test.test_getcommunicator x 0 myid;

  (* Test subvector accessors *)
  let u, _ = Nvector_mpiplusx.unwrap x in
  if Nvector.Ops.getlength u <> local_length then
    (printf ">>> FAILED test -- N_VGetLocalVector_MPIPlusX, Proc %d\n\n" myid;
     exit (-1));

  let y = Nvector_mpiplusx.Ops.clone x in
  let z = Nvector_mpiplusx.Ops.clone x in

  (* Standard vector operation tests *)
  if myid = 0 then printf "\nTesting standard vector operations:\n\n";

  fails += Test.test_const x local_length myid;
  fails += Test.test_linearsum x y z local_length myid;
  fails += Test.test_prod x y z local_length myid;
  fails += Test.test_div x y z local_length myid;
  fails += Test.test_scale x z local_length myid;
  fails += Test.test_abs x z local_length myid;
  fails += Test.test_inv x z local_length myid;
  fails += Test.test_addconst x z local_length myid;
  fails += Test.test_dotprod x y local_length global_length myid;
  fails += Test.test_maxnorm x local_length myid;
  fails += Test.test_wrmsnorm x y local_length myid;
  fails += Test.test_wrmsnormmask x y z local_length global_length myid;
  fails += Test.test_min x local_length myid;
  fails += Test.test_wl2norm x y local_length global_length myid;
  fails += Test.test_l1norm x local_length global_length myid;
  fails += Test.test_compare x z local_length myid;
  fails += Test.test_invtest x z local_length myid;
  fails += Test.test_constrmask x y z local_length myid;
  fails += Test.test_minquotient x y local_length myid;

  (* Fused and vector array operations tests (disabled) *)
  if myid = 0 then
    printf "\nTesting fused and vector array operations (disabled):\n\n";

  let u = Nvector_mpiplusx.Ops.clone x in
  let uany, _ = Nvector_mpiplusx.unwrap u in
  Nvector_serial.Any.enable ~with_fused_ops:false uany;

  (* fused operations *)
  fails += Test.test_linearcombination u local_length myid;
  fails += Test.test_scaleaddmulti u local_length myid;
  fails += Test.test_dotprodmulti u local_length global_length myid;

  (* vector array operations *)
  fails += Test.test_linearsumvectorarray u local_length myid;
  fails += Test.test_scalevectorarray u local_length myid;
  fails += Test.test_constvectorarray u local_length myid;
  fails += Test.test_wrmsnormvectorarray u local_length myid;
  fails += Test.test_wrmsnormmaskvectorarray u local_length global_length myid;
  fails += Test.test_scaleaddmultivectorarray u local_length myid;
  fails += Test.test_linearcombinationvectorarray u local_length myid;

  (* Fused and vector array operations tests (enabled) *)
  if myid = 0 then
    printf "\nTesting fused and vector array operations (enabled):\n\n";

  let u = Nvector_mpiplusx.Ops.clone x in
  let uany, _ = Nvector_mpiplusx.unwrap u in
  Nvector_serial.Any.enable ~with_fused_ops:true uany;

  (* fused operations *)
  fails += Test.test_linearcombination u local_length myid;
  fails += Test.test_scaleaddmulti u local_length myid;
  fails += Test.test_dotprodmulti u local_length global_length myid;

  (* vector array operations *)
  fails += Test.test_linearsumvectorarray u local_length myid;
  fails += Test.test_scalevectorarray u local_length myid;
  fails += Test.test_constvectorarray u local_length myid;
  fails += Test.test_wrmsnormvectorarray u local_length myid;
  fails += Test.test_wrmsnormmaskvectorarray u local_length global_length myid;
  fails += Test.test_scaleaddmultivectorarray u local_length myid;
  fails += Test.test_linearcombinationvectorarray u local_length myid;

  if myid = 0 then
    printf "\nTesting local reduction operations:\n\n";

  fails += Test.test_dotprodlocal x y local_length myid;
  fails += Test.test_maxnormlocal x local_length myid;
  fails += Test.test_minlocal x local_length myid;
  fails += Test.test_l1normlocal x local_length myid;
  fails += Test.test_wsqrsumlocal x y local_length myid;
  fails += Test.test_wsqrsummasklocal x y z local_length myid;
  fails += Test.test_invtestlocal x z local_length myid;
  fails += Test.test_constrmasklocal x y z local_length myid;
  fails += Test.test_minquotientlocal x y local_length myid;

  (* local fused reduction operations *)
  if Test_nvector.compat_ge600 then begin
    if myid = 0 then printf "\nTesting local fused reduction operations:\n\n";
    let v = Nvector_mpiplusx.Ops.clone x in
    let vany, _ = Nvector_mpiplusx.unwrap v in
    Nvector_serial.Any.enable ~with_fused_ops:true vany;
    fails += Test.test_dotprodmultilocal v local_length myid;
    fails += Test.test_dotprodmultiallreduce v local_length myid
  end;

  (* XBraid interface operations *)
  if myid = 0 then
    printf "\nTesting XBraid interface operations:\n\n";

  fails += Test.test_bufsize x local_length myid;
  fails += Test.test_bufpack x local_length myid;
  fails += Test.test_bufunpack x local_length myid;

  (* Print result *)
  if !fails <> 0 then printf "FAIL: NVector module failed %d tests \n\n" !fails
  else if myid = 0 then printf "SUCCESS: NVector module passed all tests \n\n"

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
