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
    type data = Sundials.RealArray.t
    let n_vgetarray = Nvector.unwrap
    let get = Sundials.RealArray.get
    let set = Sundials.RealArray.set
    let max_time x t = t
    let sync_device () = ()
  end)

module Nvector_array_ops =
  functor (N : Nvector.NVECTOR with type data = float array) ->
  (struct
    include N.Ops
    let get_id = Nvector.get_id
    type data = float array
    let n_vgetarray = Nvector.unwrap
    let get = Array.get
    let set = Array.set
    let max_time x t = t
    let sync_device () = ()
  end)

(* Custom nvector with serial operations reimplemented in OCaml *)
module Custom_serial =
  Nvector_custom.MakeOps (struct
    type data = Sundials.RealArray.t
    let ops = { (* {{{ *)
      Nvector_custom.n_vcheck
        = (fun x y -> Sundials.RealArray.(length x = length y));
      Nvector_custom.n_vclone        = Nvector_serial.DataOps.n_vclone;
      Nvector_custom.n_vspace        = Some Nvector_serial.DataOps.n_vspace;
      Nvector_custom.n_vlinearsum    = Nvector_serial.DataOps.n_vlinearsum;
      Nvector_custom.n_vconst        = Nvector_serial.DataOps.n_vconst;
      Nvector_custom.n_vprod         = Nvector_serial.DataOps.n_vprod;
      Nvector_custom.n_vdiv          = Nvector_serial.DataOps.n_vdiv;
      Nvector_custom.n_vscale        = Nvector_serial.DataOps.n_vscale;
      Nvector_custom.n_vabs          = Nvector_serial.DataOps.n_vabs;
      Nvector_custom.n_vinv          = Nvector_serial.DataOps.n_vinv;
      Nvector_custom.n_vaddconst     = Nvector_serial.DataOps.n_vaddconst;
      Nvector_custom.n_vmaxnorm      = Nvector_serial.DataOps.n_vmaxnorm;
      Nvector_custom.n_vwrmsnorm     = Nvector_serial.DataOps.n_vwrmsnorm;
      Nvector_custom.n_vmin          = Nvector_serial.DataOps.n_vmin;
      Nvector_custom.n_vdotprod      = Nvector_serial.DataOps.n_vdotprod;
      Nvector_custom.n_vcompare      = Nvector_serial.DataOps.n_vcompare;
      Nvector_custom.n_vinvtest      = Nvector_serial.DataOps.n_vinvtest;
      Nvector_custom.n_vwl2norm      = Some Nvector_serial.DataOps.n_vwl2norm;
      Nvector_custom.n_vl1norm       = Some Nvector_serial.DataOps.n_vl1norm;
      Nvector_custom.n_vwrmsnormmask = Some Nvector_serial.DataOps.n_vwrmsnormmask;
      Nvector_custom.n_vconstrmask   = Some Nvector_serial.DataOps.n_vconstrmask;
      Nvector_custom.n_vminquotient  = Some Nvector_serial.DataOps.n_vminquotient;
      Nvector_custom.n_vlinearcombination
        = Some Nvector_serial.DataOps.n_vlinearcombination;
      Nvector_custom.n_vscaleaddmulti
        = Some Nvector_serial.DataOps.n_vscaleaddmulti;
      Nvector_custom.n_vdotprodmulti
        = Some Nvector_serial.DataOps.n_vdotprodmulti;
      Nvector_custom.n_vlinearsumvectorarray
        = Some Nvector_serial.DataOps.n_vlinearsumvectorarray;
      Nvector_custom.n_vscalevectorarray
        = Some Nvector_serial.DataOps.n_vscalevectorarray;
      Nvector_custom.n_vconstvectorarray
        = Some Nvector_serial.DataOps.n_vconstvectorarray;
      Nvector_custom.n_vwrmsnormvectorarray
        = Some Nvector_serial.DataOps.n_vwrmsnormvectorarray;
      Nvector_custom.n_vwrmsnormmaskvectorarray
        = Some Nvector_serial.DataOps.n_vwrmsnormmaskvectorarray;
      Nvector_custom.n_vscaleaddmultivectorarray
        = Some Nvector_serial.DataOps.n_vscaleaddmultivectorarray;
      Nvector_custom.n_vlinearcombinationvectorarray
        = Some Nvector_serial.DataOps.n_vlinearcombinationvectorarray;
    } (* }}} *)
  end)

(* Custom nvector with parallel operations reimplemented in OCaml *)
module Custom_array1 =
  Nvector_custom.MakeOps (struct
    type data = float array
    let ops = { (* {{{ *)
      Nvector_custom.n_vcheck
        = (fun x y -> Array.(length x = length y));
      Nvector_custom.n_vclone        = Nvector_array.DataOps.n_vclone;
      Nvector_custom.n_vspace        = Some Nvector_array.DataOps.n_vspace;
      Nvector_custom.n_vlinearsum    = Nvector_array.DataOps.n_vlinearsum;
      Nvector_custom.n_vconst        = Nvector_array.DataOps.n_vconst;
      Nvector_custom.n_vprod         = Nvector_array.DataOps.n_vprod;
      Nvector_custom.n_vdiv          = Nvector_array.DataOps.n_vdiv;
      Nvector_custom.n_vscale        = Nvector_array.DataOps.n_vscale;
      Nvector_custom.n_vabs          = Nvector_array.DataOps.n_vabs;
      Nvector_custom.n_vinv          = Nvector_array.DataOps.n_vinv;
      Nvector_custom.n_vaddconst     = Nvector_array.DataOps.n_vaddconst;
      Nvector_custom.n_vmaxnorm      = Nvector_array.DataOps.n_vmaxnorm;
      Nvector_custom.n_vwrmsnorm     = Nvector_array.DataOps.n_vwrmsnorm;
      Nvector_custom.n_vmin          = Nvector_array.DataOps.n_vmin;
      Nvector_custom.n_vdotprod      = Nvector_array.DataOps.n_vdotprod;
      Nvector_custom.n_vcompare      = Nvector_array.DataOps.n_vcompare;
      Nvector_custom.n_vinvtest      = Nvector_array.DataOps.n_vinvtest;
      Nvector_custom.n_vwl2norm      = Some Nvector_array.DataOps.n_vwl2norm;
      Nvector_custom.n_vl1norm       = Some Nvector_array.DataOps.n_vl1norm;
      Nvector_custom.n_vwrmsnormmask = Some Nvector_array.DataOps.n_vwrmsnormmask;
      Nvector_custom.n_vconstrmask   = Some Nvector_array.DataOps.n_vconstrmask;
      Nvector_custom.n_vminquotient  = Some Nvector_array.DataOps.n_vminquotient;
      Nvector_custom.n_vlinearcombination
        = Some Nvector_array.DataOps.n_vlinearcombination;
      Nvector_custom.n_vscaleaddmulti
        = Some Nvector_array.DataOps.n_vscaleaddmulti;
      Nvector_custom.n_vdotprodmulti
        = Some Nvector_array.DataOps.n_vdotprodmulti;
      Nvector_custom.n_vlinearsumvectorarray
        = Some Nvector_array.DataOps.n_vlinearsumvectorarray;
      Nvector_custom.n_vscalevectorarray
        = Some Nvector_array.DataOps.n_vscalevectorarray;
      Nvector_custom.n_vconstvectorarray
        = Some Nvector_array.DataOps.n_vconstvectorarray;
      Nvector_custom.n_vwrmsnormvectorarray
        = Some Nvector_array.DataOps.n_vwrmsnormvectorarray;
      Nvector_custom.n_vwrmsnormmaskvectorarray
        = Some Nvector_array.DataOps.n_vwrmsnormmaskvectorarray;
      Nvector_custom.n_vscaleaddmultivectorarray
        = Some Nvector_array.DataOps.n_vscaleaddmultivectorarray;
      Nvector_custom.n_vlinearcombinationvectorarray
        = Some Nvector_array.DataOps.n_vlinearcombinationvectorarray;
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
    then ("custom serial", (module
      MakeTest (Custom_serial) (struct let id = Nvector.Custom end)))
    else if m =2
    then ("custom array-1 (dataops)", (module
      MakeArrayTest (Custom_array1) (struct let id = Nvector.Custom end)))
    else if m =3
    then ("custom array-2 (make)", (module
      MakeArrayTest (Custom_array2) (struct let id = Nvector.Custom end)))
    else ("serial", (module
      MakeTest (Nvector_serial) (struct let id = Nvector.Serial end)))
  in

  let print_timing = int_of_string Sys.argv.(2) in
  let _ = Test.set_timing (print_timing <> 0) true in

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
    fails += Test.test_n_vgetvectorid x Test.id 0;
    fails += Test.test_n_vcloneempty x 0;
    fails += Test.test_n_vclone x length 0;
    fails += Test.test_n_vcloneemptyvectorarray 5 x 0;
    fails += Test.test_n_vclonevectorarray 5 x length 0
  end;
  fails += Test.test_n_vsetarraypointer w length 0;
  fails += Test.test_n_vgetarraypointer x length 0;
  if Test_nvector.compat_ge400 then begin
    printf "\nTesting standard vector operations:\n\n";
    fails += Test.test_n_vconst x length 0
  end;
  fails += Test.test_n_vlinearsum x y z length 0;
  if not Test_nvector.compat_ge400 then begin
    fails += Test.test_n_vconst x length 0
  end;
  fails += Test.test_n_vprod x y z length 0;
  fails += Test.test_n_vdiv x y z length 0;
  fails += Test.test_n_vscale x z length 0;
  fails += Test.test_n_vabs x z length 0;
  fails += Test.test_n_vinv x z length 0;
  fails += Test.test_n_vaddconst x z length 0;
  fails += Test.test_n_vdotprod x y length length 0;
  fails += Test.test_n_vmaxnorm x length 0;
  fails += Test.test_n_vwrmsnorm x y length 0;
  if Test_nvector.compat_ge400
  then fails += Test.test_n_vwrmsnormmask x y z length length 0
  else fails += Test.test_n_vwrmsnormmask_lt400 x y z length length 0;
  fails += Test.test_n_vmin x length 0;
  fails += Test.test_n_vwl2norm x y length length 0;
  fails += Test.test_n_vl1norm x length length 0;
  fails += Test.test_n_vcompare x z length 0;
  fails += Test.test_n_vinvtest x z length 0;
  fails += Test.test_n_vconstrmask x y z length 0;
  fails += Test.test_n_vminquotient x y length 0;
  if not Test_nvector.compat_ge400 then begin
    fails += Test.test_n_vclonevectorarray 5 x length 0;
    fails += Test.test_n_vcloneemptyvectorarray 5 x 0;
    fails += Test.test_n_vcloneempty x 0;
    fails += Test.test_n_vclone x length 0
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (disabled) *)
    printf "\nTesting fused and vector array operations (disabled):\n\n";

    let u = Test.make ~with_fused_ops:false length 0.0 in

    (* fused operations *)
    fails += Test.test_n_vlinearcombination u length 0;
    fails += Test.test_n_vscaleaddmulti u length 0;
    fails += Test.test_n_vdotprodmulti u length length 0;

    (* vector array operations *)
    fails += Test.test_n_vlinearsumvectorarray u length 0;
    fails += Test.test_n_vscalevectorarray u length 0;
    fails += Test.test_n_vconstvectorarray u length 0;
    fails += Test.test_n_vwrmsnormvectorarray u length 0;
    fails += Test.test_n_vwrmsnormmaskvectorarray u length length 0;
    fails += Test.test_n_vscaleaddmultivectorarray u length 0;
    fails += Test.test_n_vlinearcombinationvectorarray u length 0
  end;

  if Test_nvector.compat_ge400 then begin
    (* Fused and vector array operations tests (enabled) *)
    printf "\nTesting fused and vector array operations (enabled):\n\n";

    let u = Test.make ~with_fused_ops:true length 0.0 in

    (* fused operations *)
    fails += Test.test_n_vlinearcombination u length 0;
    fails += Test.test_n_vscaleaddmulti u length 0;
    fails += Test.test_n_vdotprodmulti u length length 0;

    (* vector array operations *)
    fails += Test.test_n_vlinearsumvectorarray u length 0;
    fails += Test.test_n_vscalevectorarray u length 0;
    fails += Test.test_n_vconstvectorarray u length 0;
    fails += Test.test_n_vwrmsnormvectorarray u length 0;
    fails += Test.test_n_vwrmsnormmaskvectorarray u length length 0;
    fails += Test.test_n_vscaleaddmultivectorarray u length 0;
    fails += Test.test_n_vlinearcombinationvectorarray u length 0
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
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
