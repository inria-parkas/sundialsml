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
  type contents = Sundials.RealArray.t
  let getarray nv = let d, _, _ = Nvector.unwrap nv in d
  let get = Sundials.RealArray.get
  let set = Sundials.RealArray.set
  let max_time x t =
    let _, _, comm = Nvector.unwrap x in
    (* get max time across all MPI ranks *)
    Mpi.(reduce_float t Max 0 comm)
  let sync_device () = ()
end

module Nvector_generic_ops =
struct
  include Nvector.Ops
  let get_id = Nvector.get_id

  type contents = Sundials.RealArray.t
  let getarray x = match Nvector.unwrap x with
                   | Nvector_parallel.Par (xd, _, _) -> xd
                   | _ -> raise Nvector.BadGenericType

  let get = Sundials.RealArray.get
  let set = Sundials.RealArray.set

  let max_time xv t =
    match Nvector.unwrap xv with
    | Nvector_parallel.Par (_, _, comm) ->
      (* get max time across all MPI ranks *)
      Mpi.(reduce_float t Max 0 comm)
    | _ -> raise Nvector.BadGenericType
  let sync_device () = ()
end

(* Custom nvector with parallel operations reimplemented in OCaml *)
module Custom_parallel1 =
  Nvector_custom.MakeOps (struct
    type data = Nvector_parallel.data
    let ops = { (* {{{ *)
      Nvector_custom.check
        = (fun (x, xg, xc) (y, yg, yc) ->
            Sundials.RealArray.(length x = length y)
            && (xg <= yg) && (xc == yc)
          );
      Nvector_custom.clone        = Nvector_parallel.DataOps.clone;
      Nvector_custom.space        = Some Nvector_parallel.DataOps.space;
      Nvector_custom.getlength    = Nvector_parallel.DataOps.getlength;
      Nvector_custom.getcommunicator = None;
      Nvector_custom.linearsum    = Nvector_parallel.DataOps.linearsum;
      Nvector_custom.const        = Nvector_parallel.DataOps.const;
      Nvector_custom.prod         = Nvector_parallel.DataOps.prod;
      Nvector_custom.div          = Nvector_parallel.DataOps.div;
      Nvector_custom.scale        = Nvector_parallel.DataOps.scale;
      Nvector_custom.abs          = Nvector_parallel.DataOps.abs;
      Nvector_custom.inv          = Nvector_parallel.DataOps.inv;
      Nvector_custom.addconst     = Nvector_parallel.DataOps.addconst;
      Nvector_custom.maxnorm      = Nvector_parallel.DataOps.maxnorm;
      Nvector_custom.wrmsnorm     = Nvector_parallel.DataOps.wrmsnorm;
      Nvector_custom.min          = Nvector_parallel.DataOps.min;
      Nvector_custom.dotprod      = Nvector_parallel.DataOps.dotprod;
      Nvector_custom.compare      = Nvector_parallel.DataOps.compare;
      Nvector_custom.invtest      = Nvector_parallel.DataOps.invtest;
      Nvector_custom.wl2norm      = Some Nvector_parallel.DataOps.wl2norm;
      Nvector_custom.l1norm       = Some Nvector_parallel.DataOps.l1norm;
      Nvector_custom.wrmsnormmask = Some Nvector_parallel.DataOps.wrmsnormmask;
      Nvector_custom.constrmask   = Some Nvector_parallel.DataOps.constrmask;
      Nvector_custom.minquotient  = Some Nvector_parallel.DataOps.minquotient;
      Nvector_custom.linearcombination
        = Some Nvector_parallel.DataOps.linearcombination;
      Nvector_custom.scaleaddmulti
        = Some Nvector_parallel.DataOps.scaleaddmulti;
      Nvector_custom.dotprodmulti
        = Some Nvector_parallel.DataOps.dotprodmulti;
      Nvector_custom.linearsumvectorarray
        = Some Nvector_parallel.DataOps.linearsumvectorarray;
      Nvector_custom.scalevectorarray
        = Some Nvector_parallel.DataOps.scalevectorarray;
      Nvector_custom.constvectorarray
        = Some Nvector_parallel.DataOps.constvectorarray;
      Nvector_custom.wrmsnormvectorarray
        = Some Nvector_parallel.DataOps.wrmsnormvectorarray;
      Nvector_custom.wrmsnormmaskvectorarray
        = Some Nvector_parallel.DataOps.wrmsnormmaskvectorarray;
      Nvector_custom.scaleaddmultivectorarray
        = Some Nvector_parallel.DataOps.scaleaddmultivectorarray;
      Nvector_custom.linearcombinationvectorarray
        = Some Nvector_parallel.DataOps.linearcombinationvectorarray;

      Nvector_custom.dotprod_local     = Some Nvector_parallel.DataOps.Local.dotprod;
      Nvector_custom.maxnorm_local     = Some Nvector_parallel.DataOps.Local.maxnorm;
      Nvector_custom.min_local         = Some Nvector_parallel.DataOps.Local.min;
      Nvector_custom.l1norm_local      = Some Nvector_parallel.DataOps.Local.l1norm;
      Nvector_custom.invtest_local     = Some Nvector_parallel.DataOps.Local.invtest;
      Nvector_custom.constrmask_local  = Some Nvector_parallel.DataOps.Local.constrmask;
      Nvector_custom.minquotient_local = Some Nvector_parallel.DataOps.Local.minquotient;
      Nvector_custom.wsqrsum_local     = Some Nvector_parallel.DataOps.Local.wsqrsum;
      Nvector_custom.wsqrsummask_local = Some Nvector_parallel.DataOps.Local.wsqrsummask;
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
      let fill a v = Sundials.RealArray.fill a v
      let make   = Sundials.RealArray.make
      let clone  = Sundials.RealArray.copy
      let length = Sundials.RealArray.length
    end)
    let ops = { (* {{{ *)
      Nvector_custom.check
        = (fun (x, xg, xc) (y, yg, yc) ->
            Sundials.RealArray.(length x = length y)
            && (xg <= yg) && (xc == yc)
          );
      Nvector_custom.clone        = DataOps.clone;
      Nvector_custom.space        = Some DataOps.space;
      Nvector_custom.getlength    = Nvector_parallel.DataOps.getlength;
      Nvector_custom.getcommunicator = None;
      Nvector_custom.linearsum    = DataOps.linearsum;
      Nvector_custom.const        = DataOps.const;
      Nvector_custom.prod         = DataOps.prod;
      Nvector_custom.div          = DataOps.div;
      Nvector_custom.scale        = DataOps.scale;
      Nvector_custom.abs          = DataOps.abs;
      Nvector_custom.inv          = DataOps.inv;
      Nvector_custom.addconst     = DataOps.addconst;
      Nvector_custom.maxnorm      = DataOps.maxnorm;
      Nvector_custom.wrmsnorm     = DataOps.wrmsnorm;
      Nvector_custom.min          = DataOps.min;
      Nvector_custom.dotprod      = DataOps.dotprod;
      Nvector_custom.compare      = DataOps.compare;
      Nvector_custom.invtest      = DataOps.invtest;
      Nvector_custom.wl2norm      = Some DataOps.wl2norm;
      Nvector_custom.l1norm       = Some DataOps.l1norm;
      Nvector_custom.wrmsnormmask = Some DataOps.wrmsnormmask;
      Nvector_custom.constrmask   = Some DataOps.constrmask;
      Nvector_custom.minquotient  = Some DataOps.minquotient;
      Nvector_custom.linearcombination = Some DataOps.linearcombination;
      Nvector_custom.scaleaddmulti = Some DataOps.scaleaddmulti;
      Nvector_custom.dotprodmulti = Some DataOps.dotprodmulti;
      Nvector_custom.linearsumvectorarray
        = Some DataOps.linearsumvectorarray;
      Nvector_custom.scalevectorarray = Some DataOps.scalevectorarray;
      Nvector_custom.constvectorarray = Some DataOps.constvectorarray;
      Nvector_custom.wrmsnormvectorarray
        = Some DataOps.wrmsnormvectorarray;
      Nvector_custom.wrmsnormmaskvectorarray
        = Some DataOps.wrmsnormmaskvectorarray;
      Nvector_custom.scaleaddmultivectorarray
        = Some DataOps.scaleaddmultivectorarray;
      Nvector_custom.linearcombinationvectorarray
        = Some DataOps.linearcombinationvectorarray;

      Nvector_custom.dotprod_local     = Some Nvector_parallel.DataOps.Local.dotprod;
      Nvector_custom.maxnorm_local     = Some Nvector_parallel.DataOps.Local.maxnorm;
      Nvector_custom.min_local         = Some Nvector_parallel.DataOps.Local.min;
      Nvector_custom.l1norm_local      = Some Nvector_parallel.DataOps.Local.l1norm;
      Nvector_custom.invtest_local     = Some Nvector_parallel.DataOps.Local.invtest;
      Nvector_custom.constrmask_local  = Some Nvector_parallel.DataOps.Local.constrmask;
      Nvector_custom.minquotient_local = Some Nvector_parallel.DataOps.Local.minquotient;
      Nvector_custom.wsqrsum_local     = Some Nvector_parallel.DataOps.Local.wsqrsum;
      Nvector_custom.wsqrsummask_local = Some Nvector_parallel.DataOps.Local.wsqrsummask;
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

module MakeAnyTest =
  functor (NI : sig
                  val id : Nvector.nvector_id
                  val comm : Mpi.communicator
                  val global_length : int
                end) ->
  struct
  include Test_nvector.Test (Nvector_generic_ops)

  let make ?with_fused_ops n v =
    Nvector_parallel.Any.wrap ?with_fused_ops
      (Sundials.RealArray.make n v, NI.global_length, NI.comm)

  let id = Nvector.Parallel
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
    then ("generic parallel", (module MakeAnyTest
               (struct let id = Nvector.Custom
                       let comm = comm
                       let global_length = global_length
                end)))
    else if m =2
    then ("custom-1 parallel (MPI)", (module
      MakeTest (Custom_parallel1)
               (struct let id = Nvector.Custom
                       let comm = comm
                       let global_length = global_length
                end)))
    else if m =3
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
     printf "Testing the %s N_Vector \n" name;
     printf "Vector global length %d \n" global_length;
     printf "MPI processes %d \n" nprocs
  end;
  let w = Test.make local_length 0.0 in
  let x = Test.make local_length 0.0 in

  (* NVector Test *)
  if Test_nvector.compat_ge400 then begin
    fails += Test.test_getvectorid x Test.id myid;
    if Test_nvector.compat_ge500 then begin
      fails += Test.test_getlength x myid;
      fails += Test.test_getcommunicator x 0 myid;
    end;
    fails += Test.test_cloneempty x myid;
    fails += Test.test_clone x local_length myid;
    fails += Test.test_cloneemptyvectorarray 5 x myid;
    fails += Test.test_clonevectorarray 5 x local_length myid
  end;

  (* Test setting/getting array data *)
  fails += Test.test_setarraypointer w local_length myid;
  fails += Test.test_getarraypointer x local_length myid;

  (* Clone additional vectors for testing *)
  let y = Test.make local_length 0.0
  and z = Test.make local_length 0.0
  in

  (* Standard vector operation tests *)
  if Test_nvector.compat_ge400 && myid = 0
  then printf "\nTesting standard vector operations:\n\n";

  if Test_nvector.compat_ge400 then begin
    fails += Test.test_const x local_length myid;
    fails += Test.test_linearsum x y z local_length myid;
  end else begin
    fails += Test.test_linearsum x y z local_length myid;
    fails += Test.test_const x local_length myid;
  end;
  fails += Test.test_prod x y z local_length myid;
  fails += Test.test_div x y z local_length myid;
  fails += Test.test_scale x z local_length myid;
  fails += Test.test_abs x z local_length myid;
  fails += Test.test_inv x z local_length myid;
  fails += Test.test_addconst x z local_length myid;
  fails += Test.test_dotprod x y local_length global_length myid;
  fails += Test.test_maxnorm x local_length myid;
  fails += Test.test_wrmsnorm x y local_length myid;
  if Test_nvector.compat_ge400
  then fails += Test.test_wrmsnormmask x y z local_length global_length myid
  else fails += Test.test_wrmsnormmask_lt400 x y z local_length global_length myid;
  fails += Test.test_min x local_length myid;
  fails += Test.test_wl2norm x y local_length global_length myid;
  fails += Test.test_l1norm x local_length global_length myid;
  fails += Test.test_compare x z local_length myid;
  fails += Test.test_invtest x z local_length myid;
  fails += Test.test_constrmask x y z local_length myid;
  fails += Test.test_minquotient x y local_length myid;
  if not Test_nvector.compat_ge400 then begin
    fails += Test.test_clonevectorarray 5 x local_length myid;
    fails += Test.test_cloneemptyvectorarray 5 x myid;
    fails += Test.test_cloneempty x myid;
    fails += Test.test_clone x local_length myid;
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (disabled) *)
    if myid = 0 then
      printf "\nTesting fused and vector array operations (disabled):\n\n";

    (* create vector and disable all fused and vector array operations *)
    let u = Test.make ~with_fused_ops:false local_length 0.0 in

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
    fails += Test.test_linearcombinationvectorarray u local_length myid
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (enabled) *)
    if myid = 0 then
      printf "\nTesting fused and vector array operations (enabled):\n\n";

    let u = Test.make ~with_fused_ops:false local_length 0.0 in

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
    fails += Test.test_linearcombinationvectorarray u local_length myid
  end;

  (* local reduction operations *)
  if Test_nvector.compat_ge500 then begin
    if myid = 0 then printf "\nTesting local reduction operations:\n\n";

    fails += Test.test_dotprodlocal x y local_length myid;
    fails += Test.test_maxnormlocal x local_length myid;
    fails += Test.test_minlocal x local_length myid;
    fails += Test.test_l1normlocal x local_length myid;
    fails += Test.test_wsqrsumlocal x y local_length myid;
    fails += Test.test_wsqrsummasklocal x y z local_length myid;
    fails += Test.test_invtestlocal x z local_length myid;
    fails += Test.test_constrmasklocal x y z local_length myid;
    fails += Test.test_minquotientlocal x y local_length myid
  end;

  (* XBraid interface operations *)
  if Test_nvector.compat_ge540 then begin
    if myid = 0 then printf "\nTesting XBraid interface operations:\n\n";

    fails += Test.test_bufsize x local_length myid;
    fails += Test.test_bufpack x local_length myid;
    fails += Test.test_bufunpack x local_length myid
  end;

  (* Print result *)
  if myid = 0 then begin
    if !fails <> 0 then
      printf "FAIL: NVector module failed %d tests \n%s\n" !fails
        (if Test_nvector.compat_ge400 then "" else " ")
    else
      printf "SUCCESS: NVector module passed all tests%s\n%s\n"
        (if Test_nvector.compat_ge400 then " " else ", Proc 0 ")
        (if Test_nvector.compat_ge400 then "" else " ")
  end;
  match Mpi.(allreduce_int (!fails) Max comm) with
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
