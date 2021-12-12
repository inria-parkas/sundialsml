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
 * This is the testing routine to check the NVECTOR MPIManyVector
 * (parallel, intercommunicator) module implementation.
 * -----------------------------------------------------------------*)

module ROArray = Sundials.ROArray
module RealArray = Sundials.RealArray

module Nvector_mpimanyvector_ops =
  struct
    include Nvector_mpimany.Ops
    let get_id = Nvector.get_id
    type contents = Nvector.any ROArray.t * Mpi.communicator * int
    let getarray = Nvector_mpimany.unwrap

    let get (xxs, _, _) i =
      let xsub0 = Nvector_serial.Any.unwrap (ROArray.get xxs 0) in
      let xsub1, _, _ = Nvector_parallel.Any.unwrap (ROArray.get xxs 1) in
      let x0len = RealArray.length xsub0 in
      if i < x0len then Sundials.RealArray.get xsub0 i
      else Sundials.RealArray.get xsub1 (i - x0len)

    let set (xxs, _, _) i v =
      let xsub0 = Nvector_serial.Any.unwrap (ROArray.get xxs 0) in
      let xsub1, _, _ = Nvector_parallel.Any.unwrap (ROArray.get xxs 1) in
      let x0len = RealArray.length xsub0 in
      if i < x0len then RealArray.set xsub0 i v
      else RealArray.set xsub1 (i - x0len) v

    let max_time x t =
      let _, comm, _ = Nvector_mpimany.unwrap x in
      (* get max time across all MPI ranks *)
      Mpi.(reduce_float t Max 0 comm)

    let sync_device () = ()
  end

module Test =
  struct
  include Test_nvector.Test (Nvector_mpimanyvector_ops)

  let id = Nvector.MpiManyVector
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

  (* check input and set vector globlen *)
  if Array.length Sys.argv < 4 then
    (printf "ERROR: THREE (3) Inputs required: subvector 1 local vector globlen, subvector 2 local vector globlen, print timing \n";
     exit (-1))
  ;

  let loclen1 = int_of_string Sys.argv.(1) in
  if loclen1 <= 0 then
    (printf "ERROR: local subvector 1 globlen must be a positive integer \n";
     exit (-1))
  ;

  let loclen2 = int_of_string Sys.argv.(2) in
  if loclen2 <= 0 then
    (printf "ERROR: local subvector 2 globlen must be a positive integer \n";
     exit (-1))
  ;

  (* Split main communicator into even/odd subcommunicators *)
  let subcomm = Mpi.comm_split comm (myid mod 2) 0 in
  let subprocs = Mpi.comm_size subcomm in
  let subid = Mpi.comm_rank subcomm in

  let globlen = subprocs * loclen2 in

  (* overall local length *)
  let local_length = loclen1 + loclen2 in

  (* overall global length *)
  let global_length = nprocs * (loclen1 + loclen2) in

  let print_timing = int_of_string Sys.argv.(3) in
  let _ = Test.set_timing (print_timing <> 0) (myid = 0) in

  if subid = 0 then begin
    if myid mod 2 = 0 then begin
      printf "Testing the ManyVector (parallel, custom comm) N_Vector\n";
      printf "  even subcomm: Vector 1 (serial) local length %d\n" loclen1;
      printf "  even subcomm: Vector 2 (parallel) global length %d\n" globlen;
      printf "  even subcomm processes %d\n" subprocs
    end else begin
      printf "  odd subcomm: Vector 1 (serial) local length %d\n" loclen1;
      printf "  odd subcomm: Vector 2 (parallel) global length %d\n" globlen;
      printf "  odd subcomm processes %d\n" subprocs;
    end
  end;

  if myid = 0 then begin
    printf "ManyVector global length %d\n" global_length;
  end;

  let x =
    Nvector_mpimany.wrap ~comm
      (ROArray.init 2 (fun i ->
         if i = 0 then (Nvector_serial.Any.make loclen1 0.0)
         else (Nvector_parallel.Any.make loclen2 globlen subcomm 0.0)))
  in

  (* NVector Tests *)
  fails += Test.test_getvectorid x Test.id myid;
  fails += Test.test_getlength x myid;
  fails += Test.test_getcommunicator x 0 myid;

  (* Test subvector accessors *)
  let xss, _, _ = Nvector_mpimany.unwrap x in
  if ROArray.length xss <> 2 then (
    printf ">>> FAILED test -- N_VGetNumSubvectors_MPIManyVector\n";
    fails += 1);
  let x1 = Nvector_serial.Any.unwrap (ROArray.get xss 0) in
  if RealArray.length x1 <> loclen1 then (
    printf ">>> FAILED test -- N_VGetSubvector_MPIManyVector\n";
    fails += 1);
  let x2, _, _ = Nvector_parallel.Any.unwrap (ROArray.get xss 1) in
  if RealArray.length x2 <> loclen2 then (
    printf ">>> FAILED test -- N_VGetSubvector_MPIManyVector\n";
    fails += 1);

  let y = Nvector_mpimany.Ops.clone x in
  let z = Nvector_mpimany.Ops.clone x in

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

  let u = Nvector_mpimany.Ops.clone x in
  Nvector_mpimany.enable ~with_fused_ops:false u;

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

  let u = Nvector_mpimany.Ops.clone x in
  Nvector_mpimany.enable ~with_fused_ops:true u;

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

  if myid = 0 then printf "\nTesting local reduction operations:\n\n";

  fails += Test.test_dotprodlocal x y local_length myid;
  fails += Test.test_maxnormlocal x local_length myid;
  fails += Test.test_minlocal x local_length myid;
  fails += Test.test_l1normlocal x local_length myid;
  fails += Test.test_wsqrsumlocal x y local_length myid;
  fails += Test.test_wsqrsummasklocal x y z local_length myid;
  fails += Test.test_invtestlocal x z local_length myid;
  fails += Test.test_constrmasklocal x y z local_length myid;
  fails += Test.test_minquotientlocal x y local_length myid;

  (* XBraid interface operations *)
  if myid = 0 then printf "\nTesting XBraid interface operations:\n\n";

  fails += Test.test_bufsize x local_length myid;
  fails += Test.test_bufpack x local_length myid;
  fails += Test.test_bufunpack x local_length myid;

  (* Print result *)
  if !fails <> 0 then printf "FAIL: NVector module failed %d tests \n\n" !fails
  else if myid = 0 then printf "SUCCESS: NVector module passed all tests \n\n";

  (* check if any other process failed *)
  if (Mpi.(allreduce_int (!fails) Max comm) <> 0) then failwith "Tests failed"

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
