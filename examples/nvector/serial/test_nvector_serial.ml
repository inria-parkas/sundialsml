(*
 * -----------------------------------------------------------------
 * $Revision: 4137 $
 * $Date: 2014-06-15 12:26:15 -0700 (Sun, 15 Jun 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, AIST, Mar 2016.
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
 * This is the testing routine to check the NVECTOR Serial module 
 * implementation. 
 * -----------------------------------------------------------------
 *)

module Nvector_serial_ops =
  functor (N : Nvector.NVECTOR with type data = Sundials.RealArray.t) ->
  (struct
    include N.Ops
    let get_id = Nvector.get_id
    type contents = Sundials.RealArray.t
    let getarray = Nvector.unwrap
    let get = Sundials.RealArray.get
    let set = Sundials.RealArray.set
    let max_time _ t = t
    let sync_device () = ()
  end)

module Nvector_generic_ops =
  struct
    type t = Nvector.any

    include Nvector.Ops
    let get_id = Nvector.get_id

    type contents = Sundials.RealArray.t
    let getarray x = match Nvector.unwrap x with
                     | Nvector.RA xa -> xa
                     | _ -> raise Nvector.BadGenericType

    let get = Sundials.RealArray.get
    let set = Sundials.RealArray.set

    let max_time _ t = t
    let sync_device () = ()
  end

module Nvector_array_ops =
  functor (N : Nvector.NVECTOR with type data = float array) ->
  (struct
    include N.Ops
    let get_id = Nvector.get_id
    type contents = float array
    let getarray = Nvector.unwrap
    let get = Array.get
    let set = Array.set
    let max_time _ t = t
    let sync_device () = ()
  end)

(* Custom nvector with serial operations reimplemented in OCaml *)
module Custom_serial =
  Nvector_custom.MakeOps (struct
    type data = Sundials.RealArray.t
    let ops = { (* {{{ *)
      Nvector_custom.check
        = (fun x y -> Sundials.RealArray.(length x = length y));
      Nvector_custom.clone        = Nvector_serial.DataOps.clone;
      Nvector_custom.space        = Some Nvector_serial.DataOps.space;
      Nvector_custom.getlength    = Nvector_serial.DataOps.getlength;
      Nvector_custom.getlocallength = Some Nvector_serial.DataOps.getlocallength;
      Nvector_custom.getcommunicator = None;
      Nvector_custom.print
        = Some (fun d logfile -> Nvector_serial.DataOps.print ?logfile d);
      Nvector_custom.linearsum    = Nvector_serial.DataOps.linearsum;
      Nvector_custom.const        = Nvector_serial.DataOps.const;
      Nvector_custom.prod         = Nvector_serial.DataOps.prod;
      Nvector_custom.div          = Nvector_serial.DataOps.div;
      Nvector_custom.scale        = Nvector_serial.DataOps.scale;
      Nvector_custom.abs          = Nvector_serial.DataOps.abs;
      Nvector_custom.inv          = Nvector_serial.DataOps.inv;
      Nvector_custom.addconst     = Nvector_serial.DataOps.addconst;
      Nvector_custom.maxnorm      = Nvector_serial.DataOps.maxnorm;
      Nvector_custom.wrmsnorm     = Nvector_serial.DataOps.wrmsnorm;
      Nvector_custom.min          = Nvector_serial.DataOps.min;
      Nvector_custom.dotprod      = Nvector_serial.DataOps.dotprod;
      Nvector_custom.compare      = Nvector_serial.DataOps.compare;
      Nvector_custom.invtest      = Nvector_serial.DataOps.invtest;
      Nvector_custom.wl2norm      = Some Nvector_serial.DataOps.wl2norm;
      Nvector_custom.l1norm       = Some Nvector_serial.DataOps.l1norm;
      Nvector_custom.wrmsnormmask = Some Nvector_serial.DataOps.wrmsnormmask;
      Nvector_custom.constrmask   = Some Nvector_serial.DataOps.constrmask;
      Nvector_custom.minquotient  = Some Nvector_serial.DataOps.minquotient;
      Nvector_custom.linearcombination
        = Some Nvector_serial.DataOps.linearcombination;
      Nvector_custom.scaleaddmulti
        = Some Nvector_serial.DataOps.scaleaddmulti;
      Nvector_custom.dotprodmulti
        = Some Nvector_serial.DataOps.dotprodmulti;
      Nvector_custom.linearsumvectorarray
        = Some Nvector_serial.DataOps.linearsumvectorarray;
      Nvector_custom.scalevectorarray
        = Some Nvector_serial.DataOps.scalevectorarray;
      Nvector_custom.constvectorarray
        = Some Nvector_serial.DataOps.constvectorarray;
      Nvector_custom.wrmsnormvectorarray
        = Some Nvector_serial.DataOps.wrmsnormvectorarray;
      Nvector_custom.wrmsnormmaskvectorarray
        = Some Nvector_serial.DataOps.wrmsnormmaskvectorarray;
      Nvector_custom.scaleaddmultivectorarray
        = Some Nvector_serial.DataOps.scaleaddmultivectorarray;
      Nvector_custom.linearcombinationvectorarray
        = Some Nvector_serial.DataOps.linearcombinationvectorarray;

      Nvector_custom.dotprod_local     = Some Nvector_serial.DataOps.Local.dotprod;
      Nvector_custom.maxnorm_local     = Some Nvector_serial.DataOps.Local.maxnorm;
      Nvector_custom.min_local         = Some Nvector_serial.DataOps.Local.min;
      Nvector_custom.l1norm_local      = Some Nvector_serial.DataOps.Local.l1norm;
      Nvector_custom.invtest_local     = Some Nvector_serial.DataOps.Local.invtest;
      Nvector_custom.constrmask_local  = Some Nvector_serial.DataOps.Local.constrmask;
      Nvector_custom.minquotient_local = Some Nvector_serial.DataOps.Local.minquotient;
      Nvector_custom.wsqrsum_local     = Some Nvector_serial.DataOps.Local.wsqrsum;
      Nvector_custom.wsqrsummask_local = Some Nvector_serial.DataOps.Local.wsqrsummask;

      Nvector_custom.dotprodmulti_local = Some Nvector_serial.DataOps.Local.dotprodmulti;
      Nvector_custom.dotprodmulti_allreduce = None;
    } (* }}} *)
  end)

(* Custom nvector with parallel operations reimplemented in OCaml *)
module Custom_array1 =
  Nvector_custom.MakeOps (struct
    type data = float array
    let ops = { (* {{{ *)
      Nvector_custom.check
        = (fun x y -> Array.(length x = length y));
      Nvector_custom.clone        = Nvector_array.DataOps.clone;
      Nvector_custom.space        = Some Nvector_array.DataOps.space;
      Nvector_custom.getlength    = Nvector_array.DataOps.getlength;
      Nvector_custom.getlocallength = Some Nvector_array.DataOps.getlocallength;
      Nvector_custom.getcommunicator = None;
      Nvector_custom.print
        = Some (fun d logfile -> Nvector_array.DataOps.print ?logfile d);
      Nvector_custom.linearsum    = Nvector_array.DataOps.linearsum;
      Nvector_custom.const        = Nvector_array.DataOps.const;
      Nvector_custom.prod         = Nvector_array.DataOps.prod;
      Nvector_custom.div          = Nvector_array.DataOps.div;
      Nvector_custom.scale        = Nvector_array.DataOps.scale;
      Nvector_custom.abs          = Nvector_array.DataOps.abs;
      Nvector_custom.inv          = Nvector_array.DataOps.inv;
      Nvector_custom.addconst     = Nvector_array.DataOps.addconst;
      Nvector_custom.maxnorm      = Nvector_array.DataOps.maxnorm;
      Nvector_custom.wrmsnorm     = Nvector_array.DataOps.wrmsnorm;
      Nvector_custom.min          = Nvector_array.DataOps.min;
      Nvector_custom.dotprod      = Nvector_array.DataOps.dotprod;
      Nvector_custom.compare      = Nvector_array.DataOps.compare;
      Nvector_custom.invtest      = Nvector_array.DataOps.invtest;
      Nvector_custom.wl2norm      = Some Nvector_array.DataOps.wl2norm;
      Nvector_custom.l1norm       = Some Nvector_array.DataOps.l1norm;
      Nvector_custom.wrmsnormmask = Some Nvector_array.DataOps.wrmsnormmask;
      Nvector_custom.constrmask   = Some Nvector_array.DataOps.constrmask;
      Nvector_custom.minquotient  = Some Nvector_array.DataOps.minquotient;
      Nvector_custom.linearcombination
        = Some Nvector_array.DataOps.linearcombination;
      Nvector_custom.scaleaddmulti
        = Some Nvector_array.DataOps.scaleaddmulti;
      Nvector_custom.dotprodmulti
        = Some Nvector_array.DataOps.dotprodmulti;
      Nvector_custom.linearsumvectorarray
        = Some Nvector_array.DataOps.linearsumvectorarray;
      Nvector_custom.scalevectorarray
        = Some Nvector_array.DataOps.scalevectorarray;
      Nvector_custom.constvectorarray
        = Some Nvector_array.DataOps.constvectorarray;
      Nvector_custom.wrmsnormvectorarray
        = Some Nvector_array.DataOps.wrmsnormvectorarray;
      Nvector_custom.wrmsnormmaskvectorarray
        = Some Nvector_array.DataOps.wrmsnormmaskvectorarray;
      Nvector_custom.scaleaddmultivectorarray
        = Some Nvector_array.DataOps.scaleaddmultivectorarray;
      Nvector_custom.linearcombinationvectorarray
        = Some Nvector_array.DataOps.linearcombinationvectorarray;

      Nvector_custom.dotprod_local     = Some Nvector_array.DataOps.Local.dotprod;
      Nvector_custom.maxnorm_local     = Some Nvector_array.DataOps.Local.maxnorm;
      Nvector_custom.min_local         = Some Nvector_array.DataOps.Local.min;
      Nvector_custom.l1norm_local      = Some Nvector_array.DataOps.Local.l1norm;
      Nvector_custom.invtest_local     = Some Nvector_array.DataOps.Local.invtest;
      Nvector_custom.constrmask_local  = Some Nvector_array.DataOps.Local.constrmask;
      Nvector_custom.minquotient_local = Some Nvector_array.DataOps.Local.minquotient;
      Nvector_custom.wsqrsum_local     = Some Nvector_array.DataOps.Local.wsqrsum;
      Nvector_custom.wsqrsummask_local = Some Nvector_array.DataOps.Local.wsqrsummask;
      Nvector_custom.dotprodmulti_local = Some Nvector_array.DataOps.Local.dotprodmulti;
      Nvector_custom.dotprodmulti_allreduce = None;
    } (* }}} *)
  end)

(* Custom nvector with parallel operations reimplemented in OCaml using
   Nvector_parallel.MakeOps *)
module Custom_array2 =
  Nvector_array.Make (struct
    type data  = float array
    let get    = Array.get
    let set    = Array.set
    let fill   = (fun a -> Array.fill a 0 (Array.length a))
    let make   = Array.make
    let clone  = Array.copy
    let length = Array.length
  end)

module MakeTest =
  functor (N : Nvector.NVECTOR with type data = Sundials.RealArray.t)
          (NI : sig val id : Nvector.nvector_id end) ->
  struct
  include Test_nvector.Test (Nvector_serial_ops (N))

  let make ?with_fused_ops n v =
    N.wrap ?with_fused_ops (Sundials.RealArray.make n v)

  let id = NI.id
end

module MakeAnyTest =
  struct
  include Test_nvector.Test (Nvector_generic_ops)

  let make ?with_fused_ops n v =
    Nvector_serial.Any.make ?with_fused_ops n v

  let id = Nvector.Serial
end

module MakeArrayTest =
  functor (N : Nvector.NVECTOR with type data = float array)
          (NI : sig val id : Nvector.nvector_id end) ->
  struct
  include Test_nvector.Test (Nvector_array_ops (N))

  let make ?with_fused_ops n v =
    N.wrap ?with_fused_ops (Array.make n v)

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
  let fails = ref 0 in

  (* check input and set vector length *)
  if Array.length Sys.argv < 3 then
    (printf "ERROR: ONE (1) Input required: vector length, print timing \n";
     exit (-1))
  ;

  let length = int_of_string Sys.argv.(1) in
  if length <= 0 then
    (printf "ERROR: length of vector must be a positive integer \n";
     exit (-1))
  ;

  let m =
    if (Array.length Sys.argv <= 3) then 0
    else int_of_string Sys.argv.(3)
  in
  let name, (module Test : TEST) =
    if m = 1
    then ("generic serial", (module MakeAnyTest))
    else if m = 2
    then ("custom serial", (module
      MakeTest (Custom_serial) (struct let id = Nvector.Custom end)))
    else if m = 3
    then ("custom array-1 (dataops)", (module
      MakeArrayTest (Custom_array1) (struct let id = Nvector.Custom end)))
    else if m = 4
    then ("custom array-2 (make)", (module
      MakeArrayTest (Custom_array2) (struct let id = Nvector.Custom end)))
    else ("serial", (module
      MakeTest (Nvector_serial) (struct let id = Nvector.Serial end)))
  in

  let print_timing = int_of_string Sys.argv.(2) in
  let _ = Test.set_timing (print_timing <> 0) Test_nvector.compat_neq600 in

  if Test_nvector.compat_ge400
  then printf "Testing %s N_Vector \nVector length %d \n" name length
  else printf "\nRunning with vector length %d \n \n" length;

  (* Create vectors *)
  let w = Test.make length 0.0
  and x = Test.make length 0.0
  and y = Test.make length 0.0
  and z = Test.make length 0.0
  in

  (* NVector Tests *)
  if Test_nvector.compat_ge400 then begin
    fails += Test.test_getvectorid x Test.id 0;
    if Test_nvector.compat_ge500 then begin
      fails += Test.test_getlength x 0;
      fails += Test.test_getcommunicator x 0 0;
    end;
    fails += Test.test_cloneempty x 0;
    fails += Test.test_clone x length 0;
    fails += Test.test_cloneemptyvectorarray 5 x 0;
    fails += Test.test_clonevectorarray 5 x length 0
  end;
  fails += Test.test_setarraypointer w length 0;
  fails += Test.test_getarraypointer x length 0;
  if Test_nvector.compat_ge400 then begin
    printf "\nTesting standard vector operations:\n\n";
    fails += Test.test_const x length 0
  end;
  fails += Test.test_linearsum x y z length 0;
  if not Test_nvector.compat_ge400 then begin
    fails += Test.test_const x length 0
  end;
  fails += Test.test_prod x y z length 0;
  fails += Test.test_div x y z length 0;
  fails += Test.test_scale x z length 0;
  fails += Test.test_abs x z length 0;
  fails += Test.test_inv x z length 0;
  fails += Test.test_addconst x z length 0;
  fails += Test.test_dotprod x y length length 0;
  fails += Test.test_maxnorm x length 0;
  fails += Test.test_wrmsnorm x y length 0;
  if Test_nvector.compat_ge400
  then fails += Test.test_wrmsnormmask x y z length length 0
  else fails += Test.test_wrmsnormmask_lt400 x y z length length 0;
  fails += Test.test_min x length 0;
  fails += Test.test_wl2norm x y length length 0;
  fails += Test.test_l1norm x length length 0;
  fails += Test.test_compare x z length 0;
  fails += Test.test_invtest x z length 0;
  fails += Test.test_constrmask x y z length 0;
  fails += Test.test_minquotient x y length 0;
  if not Test_nvector.compat_ge400 then begin
    fails += Test.test_clonevectorarray 5 x length 0;
    fails += Test.test_cloneemptyvectorarray 5 x 0;
    fails += Test.test_cloneempty x 0;
    fails += Test.test_clone x length 0
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (disabled) *)
    printf "\nTesting fused and vector array operations (disabled):\n\n";

    let u = Test.make ~with_fused_ops:false length 0.0 in

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
    fails += Test.test_linearcombinationvectorarray u length 0
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (enabled) *)
    printf "\nTesting fused and vector array operations (enabled):\n\n";

    let u = Test.make ~with_fused_ops:true length 0.0 in

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
  end;

  (* local reduction operations *)
  if Test_nvector.compat_ge500 then begin
    printf "\nTesting local reduction operations:\n\n";

    fails += Test.test_dotprodlocal x y length 0;
    fails += Test.test_maxnormlocal x length 0;
    fails += Test.test_minlocal x length 0;
    fails += Test.test_l1normlocal x length 0;
    fails += Test.test_wsqrsumlocal x y length 0;
    fails += Test.test_wsqrsummasklocal x y z length 0;
    fails += Test.test_invtestlocal x z length 0;
    fails += Test.test_constrmasklocal x y z length 0;
    fails += Test.test_minquotientlocal x y length 0
  end;

  (* local fused reduction operations *)
  if Test_nvector.compat_ge600 then begin
    printf "\nTesting local fused reduction operations:\n\n";
    let v = Test.make ~with_fused_ops:true length 0.0 in
    fails += Test.test_dotprodmultilocal v length 0
  end;

  (* XBraid interface operations *)
  if Test_nvector.compat_ge540 then begin
    printf "\nTesting XBraid interface operations:\n\n";

    fails += Test.test_bufsize x length 0;
    fails += Test.test_bufpack x length 0;
    fails += Test.test_bufunpack x length 0
  end;

  (* Free vectors *)

  (* Print result *)
  if !fails <> 0 then
    printf "FAIL: NVector module failed %d tests \n%s\n" !fails
      (if Test_nvector.compat_ge400 then "" else " ")
  else
    printf "SUCCESS: NVector module passed all tests \n%s\n"
      (if Test_nvector.compat_ge400 then "" else " ")

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
