(*
 * -----------------------------------------------------------------
 * $Revision: 4137 $
 * $Date: 2014-06-15 12:26:15 -0700 (Sun, 15 Jun 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, AIST, Mar 2016.
 * OCaml port: Timothy Bourke, Inria, Jun 2020.
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
module Nvector_parallel_ops =
  functor (N : Nvector.NVECTOR with type data = Nvector_parallel.data) ->
struct
  include N.Ops
  let get_id = Nvector.get_id
  type data = Sundials.RealArray.t
  let n_vgetarray nv = let d, _, _ = Nvector.unwrap nv in d
  let get = Sundials.RealArray.get
  let set = Sundials.RealArray.set
  let max_time x t =
    let _, _, comm = Nvector.unwrap x in
    (* get max time across all MPI ranks *)
    Mpi.(reduce_float t Float_max 0 comm)
  let sync_device () = ()
end

(* Custom nvector with parallel operations reimplemented in OCaml *)
module Custom_parallel1 =
  Nvector_custom.MakeOps (struct
    type data = Nvector_parallel.data
    let ops = { (* {{{ *)
      Nvector_custom.n_vcheck
        = (fun (x, xg, xc) (y, yg, yc) ->
            Sundials.RealArray.(length x = length y)
            && (xg <= yg) && (xc == yc)
          );
      Nvector_custom.n_vclone        = Nvector_parallel.DataOps.n_vclone;
      Nvector_custom.n_vspace        = Some Nvector_parallel.DataOps.n_vspace;
      Nvector_custom.n_vlinearsum    = Nvector_parallel.DataOps.n_vlinearsum;
      Nvector_custom.n_vconst        = Nvector_parallel.DataOps.n_vconst;
      Nvector_custom.n_vprod         = Nvector_parallel.DataOps.n_vprod;
      Nvector_custom.n_vdiv          = Nvector_parallel.DataOps.n_vdiv;
      Nvector_custom.n_vscale        = Nvector_parallel.DataOps.n_vscale;
      Nvector_custom.n_vabs          = Nvector_parallel.DataOps.n_vabs;
      Nvector_custom.n_vinv          = Nvector_parallel.DataOps.n_vinv;
      Nvector_custom.n_vaddconst     = Nvector_parallel.DataOps.n_vaddconst;
      Nvector_custom.n_vmaxnorm      = Nvector_parallel.DataOps.n_vmaxnorm;
      Nvector_custom.n_vwrmsnorm     = Nvector_parallel.DataOps.n_vwrmsnorm;
      Nvector_custom.n_vmin          = Nvector_parallel.DataOps.n_vmin;
      Nvector_custom.n_vdotprod      = Nvector_parallel.DataOps.n_vdotprod;
      Nvector_custom.n_vcompare      = Nvector_parallel.DataOps.n_vcompare;
      Nvector_custom.n_vinvtest      = Nvector_parallel.DataOps.n_vinvtest;
      Nvector_custom.n_vwl2norm      = Some Nvector_parallel.DataOps.n_vwl2norm;
      Nvector_custom.n_vl1norm       = Some Nvector_parallel.DataOps.n_vl1norm;
      Nvector_custom.n_vwrmsnormmask = Some Nvector_parallel.DataOps.n_vwrmsnormmask;
      Nvector_custom.n_vconstrmask   = Some Nvector_parallel.DataOps.n_vconstrmask;
      Nvector_custom.n_vminquotient  = Some Nvector_parallel.DataOps.n_vminquotient;
      Nvector_custom.n_vlinearcombination
        = Some Nvector_parallel.DataOps.n_vlinearcombination;
      Nvector_custom.n_vscaleaddmulti
        = Some Nvector_parallel.DataOps.n_vscaleaddmulti;
      Nvector_custom.n_vdotprodmulti
        = Some Nvector_parallel.DataOps.n_vdotprodmulti;
      Nvector_custom.n_vlinearsumvectorarray
        = Some Nvector_parallel.DataOps.n_vlinearsumvectorarray;
      Nvector_custom.n_vscalevectorarray
        = Some Nvector_parallel.DataOps.n_vscalevectorarray;
      Nvector_custom.n_vconstvectorarray
        = Some Nvector_parallel.DataOps.n_vconstvectorarray;
      Nvector_custom.n_vwrmsnormvectorarray
        = Some Nvector_parallel.DataOps.n_vwrmsnormvectorarray;
      Nvector_custom.n_vwrmsnormmaskvectorarray
        = Some Nvector_parallel.DataOps.n_vwrmsnormmaskvectorarray;
      Nvector_custom.n_vscaleaddmultivectorarray
        = Some Nvector_parallel.DataOps.n_vscaleaddmultivectorarray;
      Nvector_custom.n_vlinearcombinationvectorarray
        = Some Nvector_parallel.DataOps.n_vlinearcombinationvectorarray;
    } (* }}} *)
  end)

(* Custom nvector with parallel operations reimplemented in OCaml using
   Nvector_parallel.MakeOps *)
module Custom_parallel2 =
  Nvector_custom.MakeOps (struct
    type data = Nvector_parallel.data
    module DataOps = Nvector_parallel.MakeOps (struct
      type local_data = Sundials.RealArray.t
      let get    = Sundials.RealArray.get
      let set    = Sundials.RealArray.set
      let fill   = Sundials.RealArray.fill
      let make   = Sundials.RealArray.make
      let clone  = Sundials.RealArray.copy
      let length = Sundials.RealArray.length
    end)
    let ops = { (* {{{ *)
      Nvector_custom.n_vcheck
        = (fun (x, xg, xc) (y, yg, yc) ->
            Sundials.RealArray.(length x = length y)
            && (xg <= yg) && (xc == yc)
          );
      Nvector_custom.n_vclone        = DataOps.n_vclone;
      Nvector_custom.n_vspace        = Some DataOps.n_vspace;
      Nvector_custom.n_vlinearsum    = DataOps.n_vlinearsum;
      Nvector_custom.n_vconst        = DataOps.n_vconst;
      Nvector_custom.n_vprod         = DataOps.n_vprod;
      Nvector_custom.n_vdiv          = DataOps.n_vdiv;
      Nvector_custom.n_vscale        = DataOps.n_vscale;
      Nvector_custom.n_vabs          = DataOps.n_vabs;
      Nvector_custom.n_vinv          = DataOps.n_vinv;
      Nvector_custom.n_vaddconst     = DataOps.n_vaddconst;
      Nvector_custom.n_vmaxnorm      = DataOps.n_vmaxnorm;
      Nvector_custom.n_vwrmsnorm     = DataOps.n_vwrmsnorm;
      Nvector_custom.n_vmin          = DataOps.n_vmin;
      Nvector_custom.n_vdotprod      = DataOps.n_vdotprod;
      Nvector_custom.n_vcompare      = DataOps.n_vcompare;
      Nvector_custom.n_vinvtest      = DataOps.n_vinvtest;
      Nvector_custom.n_vwl2norm      = Some DataOps.n_vwl2norm;
      Nvector_custom.n_vl1norm       = Some DataOps.n_vl1norm;
      Nvector_custom.n_vwrmsnormmask = Some DataOps.n_vwrmsnormmask;
      Nvector_custom.n_vconstrmask   = Some DataOps.n_vconstrmask;
      Nvector_custom.n_vminquotient  = Some DataOps.n_vminquotient;
      Nvector_custom.n_vlinearcombination = Some DataOps.n_vlinearcombination;
      Nvector_custom.n_vscaleaddmulti = Some DataOps.n_vscaleaddmulti;
      Nvector_custom.n_vdotprodmulti = Some DataOps.n_vdotprodmulti;
      Nvector_custom.n_vlinearsumvectorarray
        = Some DataOps.n_vlinearsumvectorarray;
      Nvector_custom.n_vscalevectorarray = Some DataOps.n_vscalevectorarray;
      Nvector_custom.n_vconstvectorarray = Some DataOps.n_vconstvectorarray;
      Nvector_custom.n_vwrmsnormvectorarray
        = Some DataOps.n_vwrmsnormvectorarray;
      Nvector_custom.n_vwrmsnormmaskvectorarray
        = Some DataOps.n_vwrmsnormmaskvectorarray;
      Nvector_custom.n_vscaleaddmultivectorarray
        = Some DataOps.n_vscaleaddmultivectorarray;
      Nvector_custom.n_vlinearcombinationvectorarray
        = Some DataOps.n_vlinearcombinationvectorarray;
    } (* }}} *)
  end)

module MakeTest =
  functor (N : Nvector.NVECTOR with type data = Nvector_parallel.data)
          (NI : sig
                  val id : Nvector.nvector_id
                  val comm : Mpi.communicator
                  val global_length : int
                end) ->
  struct
  include Test_nvector.Test (Nvector_parallel_ops (N))

  let make ?with_fused_ops n v =
    N.wrap ?with_fused_ops
           (Sundials.RealArray.make n v, NI.global_length, NI.comm)

  let id = NI.id
end

module type TEST = sig
  include Test_nvector.TEST
  val make : ?with_fused_ops:bool -> int -> float -> t
  val id   : Nvector.nvector_id
end

let printf = Printf.printf
let (+=) r x = r := !r + x

(* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*)
let main () =
  let fails = ref 0 in                  (* counter for test failures *)

  (* Get processor number and total number of processes *)
  let comm   = Mpi.comm_world in
  let nprocs = Mpi.comm_size comm in
  let myid   = Mpi.comm_rank comm in

  (* check input and set vector length *)
  if Array.length Sys.argv < 3 then begin
    if myid = 0 then
      printf "ERROR: TWO (2) Inputs required: vector length, print timing \n";
    exit (-1)
  end;

  let local_length = int_of_string Sys.argv.(1) in
  if local_length < 1 then begin
    if myid = 0 then
        printf "ERROR: local vector length must be a positive integer \n";
    exit (-1)
  end;

  (* global length *)
  let global_length = nprocs*local_length in

  let m =
    if (Array.length Sys.argv <= 3) then 0
    else int_of_string Sys.argv.(3)
  in
  let name, (module Test : TEST) =
    if m = 1
    then ("custom-1 parallel (MPI)", (module
      MakeTest (Custom_parallel1)
               (struct let id = Nvector.Custom
                       let comm = comm
                       let global_length = global_length
                end)))
    else if m =2
    then ("custom-2 parallel (MPI)", (module
      MakeTest (Custom_parallel2)
               (struct let id = Nvector.Custom
                       let comm = comm
                       let global_length = global_length
                end)))
    else ("parallel (MPI)", (module
      MakeTest (Nvector_parallel)
               (struct let id = Nvector.Parallel
                       let comm = comm
                       let global_length = global_length
                end)))
  in

  let print_timing = int_of_string Sys.argv.(2) in
  let _ = Test.set_timing (print_timing <> 0) (myid = 0) in

  (* Create vectors *)
  if Test_nvector.compat_ge400 && myid = 0 then begin
     printf "Testing the %S N_Vector \n" name;
     printf "Vector global length %d \n" global_length;
     printf "MPI processes %d \n" nprocs
  end;
  let w = Test.make local_length 0.0
  and x = Test.make local_length 0.0
  in

  (* NVector Test *)
  if Test_nvector.compat_ge400 then begin
    fails += Test.test_n_vgetvectorid x Test.id myid;
    fails += Test.test_n_vcloneempty x myid;
    fails += Test.test_n_vclone x local_length myid;
    fails += Test.test_n_vcloneemptyvectorarray 5 x myid;
    fails += Test.test_n_vclonevectorarray 5 x local_length myid
  end;

  (* Test setting/getting array data *)
  fails += Test.test_n_vsetarraypointer w local_length myid;
  fails += Test.test_n_vgetarraypointer x local_length myid;

  (* Clone additional vectors for testing *)
  let y = Test.make local_length 0.0
  and z = Test.make local_length 0.0
  in

  (* Standard vector operation tests *)
  if Test_nvector.compat_ge400 && myid = 0
  then printf "\nTesting standard vector operations:\n\n";

  if Test_nvector.compat_ge400 then begin
    fails += Test.test_n_vconst x local_length myid;
    fails += Test.test_n_vlinearsum x y z local_length myid;
  end else begin
    fails += Test.test_n_vlinearsum x y z local_length myid;
    fails += Test.test_n_vconst x local_length myid;
  end;
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
  if not Test_nvector.compat_ge400 then begin
    fails += Test.test_n_vclonevectorarray 5 x local_length myid;
    fails += Test.test_n_vcloneemptyvectorarray 5 x myid;
    fails += Test.test_n_vcloneempty x myid;
    fails += Test.test_n_vclone x local_length myid;
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (disabled) *)
    if myid = 0 then
      printf "\nTesting fused and vector array operations (disabled):\n\n";

    (* create vector and disable all fused and vector array operations *)
    let u = Test.make ~with_fused_ops:false local_length 0.0 in

    (* fused operations *)
    fails += Test.test_n_vlinearcombination u local_length myid;
    fails += Test.test_n_vscaleaddmulti u local_length myid;
    fails += Test.test_n_vdotprodmulti u local_length global_length myid;

    (* vector array operations *)
    fails += Test.test_n_vlinearsumvectorarray u local_length myid;
    fails += Test.test_n_vscalevectorarray u local_length myid;
    fails += Test.test_n_vconstvectorarray u local_length myid;
    fails += Test.test_n_vwrmsnormvectorarray u local_length myid;
    fails += Test.test_n_vwrmsnormmaskvectorarray u local_length global_length myid;
    fails += Test.test_n_vscaleaddmultivectorarray u local_length myid;
    fails += Test.test_n_vlinearcombinationvectorarray u local_length myid
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (enabled) *)
    if myid = 0 then
      printf "\nTesting fused and vector array operations (enabled):\n\n";

    let u = Test.make ~with_fused_ops:false local_length 0.0 in

    (* fused operations *)
    fails += Test.test_n_vlinearcombination u local_length myid;
    fails += Test.test_n_vscaleaddmulti u local_length myid;
    fails += Test.test_n_vdotprodmulti u local_length global_length myid;

    (* vector array operations *)
    fails += Test.test_n_vlinearsumvectorarray u local_length myid;
    fails += Test.test_n_vscalevectorarray u local_length myid;
    fails += Test.test_n_vconstvectorarray u local_length myid;
    fails += Test.test_n_vwrmsnormvectorarray u local_length myid;
    fails += Test.test_n_vwrmsnormmaskvectorarray u local_length global_length myid;
    fails += Test.test_n_vscaleaddmultivectorarray u local_length myid;
    fails += Test.test_n_vlinearcombinationvectorarray u local_length myid
  end;

  (* Print result *)
  if myid = 0 then begin
    if !fails <> 0 then
      printf "FAIL: NVector module failed %d tests \n%s\n" !fails
        (if Test_nvector.compat_ge400 then "" else " ")
    else
      printf "SUCCESS: NVector module passed all tests \n%s\n"
        (if Test_nvector.compat_ge400 then "" else " ")
  end;
  match Mpi.(allreduce_int (!fails) Int_max comm) with
  | 0 -> ();
  | n -> failwith "Tests failed"

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
