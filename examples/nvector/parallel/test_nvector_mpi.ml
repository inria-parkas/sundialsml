(*
 * -----------------------------------------------------------------
 * $Revision: 4137 $
 * $Date: 2014-06-15 12:26:15 -0700 (Sun, 15 Jun 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, AIST, Mar 2016.
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
 * This is the testing routine to check the NVECTOR Parallel module
 * implementation.
 * -----------------------------------------------------------------
 *)
module Nvector_parallel_ops : Test_nvector.NVECTOR_OPS_EXT
  with type t = Nvector_parallel.t
  =
struct
  include Nvector_parallel.Ops
  let n_vgetarray = Nvector_parallel.local_array
end

module Test = Test_nvector.Test (Nvector_parallel_ops)

let printf = Printf.printf
let (+=) r x = r := !r + x


(* local vector length *)
let veclen = 10000

(* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in                  (* counter for test failures *)

  (* Get processor number and total number of processes *)
  let comm   = Mpi.comm_world in
  let nprocs = Mpi.comm_size comm in
  let myid   = Mpi.comm_rank comm in

  (* set local and global lengths *)
  let local_length = veclen in
  let global_length = nprocs*local_length in

  (* Create vectors *)
  let w = Nvector_parallel.make local_length global_length comm 0.0
  and x = Nvector_parallel.make local_length global_length comm 0.0
  and y = Nvector_parallel.make local_length global_length comm 0.0
  and z = Nvector_parallel.make local_length global_length comm 0.0
  in

  (* NVector Test *)
  fails += Test.test_n_vsetarraypointer w local_length myid;
  fails += Test.test_n_vgetarraypointer x local_length myid;
  fails += Test.test_n_vlinearsum x y z local_length myid;
  fails += Test.test_n_vconst x local_length myid;
  fails += Test.test_n_vprod x y z local_length myid;
  fails += Test.test_n_vdiv x y z local_length myid;
  fails += Test.test_n_vscale x z local_length myid;
  fails += Test.test_n_vabs x z local_length myid;
  fails += Test.test_n_vinv x z local_length myid;
  fails += Test.test_n_vaddconst x z local_length myid;
  fails += Test.test_n_vdotprod x y local_length global_length myid;
  fails += Test.test_n_vmaxnorm x local_length myid;
  fails += Test.test_n_vwrmsnorm x y local_length myid;
  fails += Test.test_n_vwrmsnormmask x y z local_length global_length myid;
  fails += Test.test_n_vmin x local_length myid;
  fails += Test.test_n_vwl2norm x y local_length global_length myid;
  fails += Test.test_n_vl1norm x local_length global_length myid;
  fails += Test.test_n_vcompare x z local_length myid;
  fails += Test.test_n_vinvtest x z local_length myid;
  fails += Test.test_n_vconstrmask x y z local_length myid;
  fails += Test.test_n_vminquotient x y local_length myid;
  fails += Test.test_n_vclonevectorarray 5 x local_length myid;
  fails += Test.test_n_vcloneemptyvectorarray 5 x myid;
  fails += Test.test_n_vcloneempty x myid;
  fails += Test.test_n_vclone x local_length myid;

  (* Print result *)
  if !fails <> 0 then
    printf "FAIL: NVector module failed %i tests, Proc %d \n \n" !fails myid
  else if myid = 0 then
    printf "SUCCESS: NVector module passed all tests, Proc %d \n \n" myid


(* Check environment variables for extra arguments.  *)
let reps =
  try int_of_string (Unix.getenv "NUM_REPS")
  with Not_found | Failure "int_of_string" -> 1
let gc_at_end =
  try int_of_string (Unix.getenv "GC_AT_END") <> 0
  with Not_found | Failure "int_of_string" -> false
let gc_each_rep =
  try int_of_string (Unix.getenv "GC_EACH_REP") <> 0
  with Not_found | Failure "int_of_string" -> false

(* Entry point *)
let _ =
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
