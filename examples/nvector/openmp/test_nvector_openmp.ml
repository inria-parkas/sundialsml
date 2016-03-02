(*
 * -----------------------------------------------------------------
 * $Revision: 4455 $
 * $Date: 2015-03-28 18:28:53 -0700 (Sat, 28 Mar 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine to check the open MP NVECTOR module
 * implementation.
 * -----------------------------------------------------------------
 *)
module Nvector_openmp_ops : Test_nvector.NVECTOR_OPS_EXT
  with type t = Nvector_openmp.t
  =
struct
  include Nvector_openmp.Ops
  let n_vgetarray = Nvector_openmp.unwrap
end

module Test = Test_nvector.Test (Nvector_openmp_ops)

let printf = Printf.printf
let (+=) r x = r := !r + x

(* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in

  (* check input and set vector length *)
  if Array.length Sys.argv < 4 then (
    printf "ERROR: THREE (3) arguments required: <vector length> <number of threads> <print timing>\n";
    exit (-1)
  );

  let veclen = int_of_string Sys.argv.(1) in
  if veclen <= 0 then (
    printf "ERROR: length of vector must be a positive integer \n";
    exit (-1);
  );

  let num_threads = int_of_string Sys.argv.(2) in
  if num_threads <= 0 then (
    printf "ERROR: numbber of threads must be a positive integer \n";
    exit (-1);
  );

  let print_timing = int_of_string Sys.argv.(3) in
  let _ = Test.set_timing (print_timing <> 0) in

  printf "\nRunning with vector length %d \n \n" veclen;
  printf "\nRunning with number of threads %d \n \n" num_threads;

  (* Create vectors *)
  let w = Nvector_openmp.make veclen num_threads 0.0
  and x = Nvector_openmp.make veclen num_threads 0.0
  and y = Nvector_openmp.make veclen num_threads 0.0
  and z = Nvector_openmp.make veclen num_threads 0.0
  in

  (* NVector Tests *)
  fails += Test.test_n_vsetarraypointer w veclen 0;
  fails += Test.test_n_vgetarraypointer x veclen 0;
  fails += Test.test_n_vlinearsum x y z veclen 0;
  fails += Test.test_n_vconst x veclen 0;
  fails += Test.test_n_vprod x y z veclen 0;
  fails += Test.test_n_vdiv x y z veclen 0;
  fails += Test.test_n_vscale x z veclen 0;
  fails += Test.test_n_vabs x z veclen 0;
  fails += Test.test_n_vinv x z veclen 0;
  fails += Test.test_n_vaddconst x z veclen 0;
  fails += Test.test_n_vdotprod x y veclen veclen 0;
  fails += Test.test_n_vmaxnorm x veclen 0;
  fails += Test.test_n_vwrmsnorm x y veclen 0;
  fails += Test.test_n_vwrmsnormmask x y z veclen veclen 0;
  fails += Test.test_n_vmin x veclen 0;
  fails += Test.test_n_vwl2norm x y veclen veclen 0;
  fails += Test.test_n_vl1norm x veclen veclen 0;
  fails += Test.test_n_vcompare x z veclen 0;
  fails += Test.test_n_vinvtest x z veclen 0;
  fails += Test.test_n_vconstrmask x y z veclen 0;
  fails += Test.test_n_vminquotient x y veclen 0;
  fails += Test.test_n_vclonevectorarray 5 x veclen 0;
  fails += Test.test_n_vcloneemptyvectorarray 5 x 0;
  fails += Test.test_n_vcloneempty x 0;
  fails += Test.test_n_vclone x veclen 0;

  (* Free vectors *)

  (* Print results *)
  if !fails <> 0 then
    printf "FAIL: NVector module failed %i tests \n \n" !fails
  else
    printf"SUCCESS: NVector module passed all tests \n \n"

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
