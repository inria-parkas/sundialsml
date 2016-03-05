(*
 * -----------------------------------------------------------------
 * $Revision: 4454 $
 * $Date: 2015-03-28 18:06:50 -0700 (Sat, 28 Mar 2015) $
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
 * This is the testing routine to check the POSIX Threads (Pthreads) 
 * NVECTOR module implementation which uses a LOCAL data struct to 
 * share data between threads. 
 * -----------------------------------------------------------------
 *)
module Nvector_pthreads_ops : Test_nvector.NVECTOR_OPS_EXT
  with type t = Nvector_pthreads.t
  =
struct
  include Nvector_pthreads.Ops
  let n_vgetarray = Nvector_pthreads.unwrap
end

module Test = Test_nvector.Test (Nvector_pthreads_ops)

let printf = Printf.printf
let (+=) r x = r := !r + x


(* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in

  (* check inputs, set vector length, and number of threads *)
  if Array.length Sys.argv < 4 then (
    printf "ERROR: THREE (3) Inputs required: vector length, number of threads, print timing \n";
    exit (-1)
  );

  let veclen = int_of_string Sys.argv.(1) in
  if veclen <= 0 then (
    printf("ERROR: length of vector must be a positive integer \n");
    exit (-1)
  );

  let nthreads = int_of_string Sys.argv.(2) in
  if nthreads < 1 then (
    printf "ERROR: number of threads must be at least 1 \n";
    exit (-1)
  );

  let print_timing = int_of_string Sys.argv.(3) in
  let _ = Test.set_timing (print_timing <> 0) in


  printf "\nRunning with %d threads and vector length %d \n \n" nthreads veclen;

  (* Create vectors *)
  let w = Nvector_pthreads.make nthreads veclen 0.0
  and x = Nvector_pthreads.make nthreads veclen 0.0
  and y = Nvector_pthreads.make nthreads veclen 0.0
  and z = Nvector_pthreads.make nthreads veclen 0.0
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

  (* Print result *)
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
