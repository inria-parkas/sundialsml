(*
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2018.
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
 * These test functions are designed to check a SUNMatrix module
 * implementation.
 * -----------------------------------------------------------------
 *)

(* NOTE: a number of tests in the original C code validate memory
   management or functions that have no sensible OCaml counterpart.
   Those tests are omitted in this port. However, the output should
   still exactly agree with C if print_time is disabled. *)

module Matrix = Sundials.Matrix

module type MATRIX_TESTS = sig
  type k
  type m
  type nd
  type nk

  type t = (k, m, nd, nk) Matrix.t
  type nvec = (nd, nk) Nvector.t

  val rewrap : ?context:Sundials.Context.t -> m -> t
  val check_matrix : t -> t -> float -> bool
  val check_matrix_entry : t -> float -> float -> bool
  val is_square : t -> bool

  val check_vector : nvec -> nvec -> float -> bool

  val nvec_pp : Format.formatter -> nvec -> unit
end

let isnan (x : float) = x <> x

let fneq a b tol =
  if isnan a then 1
  else if (abs_float (a -. b) /. abs_float b) > tol then 1
  else 0

external rand : unit -> int = "ml_rand"
external rand_max : unit -> int = "ml_rand_max"

let rand_max = rand_max ()

module Test (M : MATRIX_TESTS) (NV : Nvector.NVECTOR_OPS with type t = M.nvec) =
struct (* Extends all the way to the end of file.  *)

(* private functions *)

let printf = Format.printf

let print_time_flag = ref false

let print_time format time =
  if !print_time_flag then
    printf format time

external get_time : unit -> float = "get_time"
external set_timing : bool -> float = "SetTiming"

let set_timing b =
  print_time_flag := b;
  set_timing b

let clone a = M.rewrap (Matrix.((get_ops a).m_clone (unwrap a)))

(* ----------------------------------------------------------------------
 * SUNMatGetID Test
 * --------------------------------------------------------------------*)
let test_sunmatgetid a sunid myid =
  let start_time = get_time () in
  let mysunid = Matrix.get_id a in
  let stop_time = get_time () in
  if sunid <> mysunid then begin
    printf ">>> FAILED test -- SUNMatGetID, Proc %d @\n" myid;
    print_time "    SUNMatGetID Time: %22.15e @\n @\n" (stop_time -. start_time)
  end else if myid = 0 then begin
    printf "    PASSED test -- SUNMatGetID @\n";
    print_time "    SUNMatGetID Time: %22.15e @\n @\n" (stop_time -. start_time)
  end;
  if sunid = mysunid then 0 else 1

(* ----------------------------------------------------------------------
 * N_VClone Test
 * --------------------------------------------------------------------*)
let test_sunmatclone a myid =
  let tol = 1e-15 in
  try
    (* clone vector *)
    let start_time = get_time () in
    let b = try clone a with _ -> begin
        printf ">>> FAILED test -- SUNMatClone, Proc %d @\n" myid;
        printf "    After SUNMatClone, B == NULL @\n @\n";
        raise Exit
      end
    in
    let stop_time = get_time () in

    (try Matrix.blit ~src:a ~dst:b
     with _ -> printf ">>> FAILED test -- SUNMatCopy, Proc %d @\n" myid;
               raise Exit);

    if M.check_matrix b a tol then begin
      printf ">>> FAILED test -- SUNMatClone, Proc %d @\n" myid;
      printf "    Failed SUNMatClone check @\n @\n";
      raise Exit
    end;

    if (myid == 0) then begin
      (match Sundials.Config.sundials_version with
       | 2,_,_ | 3,1,0 | 3,1,1 -> printf "    PASSED test -- N_VClone @\n"
       | _ -> printf "    PASSED test -- SUNMatClone @\n");
      print_time "    SUNMatClone Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNMatZero Test
 * --------------------------------------------------------------------*)
let test_sunmatzero a myid =
  let tol = 1e-15 in
  try
    (* protect A *)
    let b = clone a in
    (* set matrix data to zero *)
    let start_time = get_time () in
    (try Matrix.set_to_zero b
     with _ -> printf ">>> FAILED test -- SUNMatZero failed on Proc %d @\n" myid;
               raise Exit);
    let stop_time = get_time () in

    (* A data should be a vector of zeros *)
    if M.check_matrix_entry b 0.0 tol then begin
      printf ">>> FAILED test -- SUNMatZero check, Proc %d @\n" myid;
      print_time "    SUNMatZero Time: %22.15e @\n @\n" (stop_time -. start_time);
      raise Exit
    end
    else if myid = 0 then begin
      printf "    PASSED test -- SUNMatZero @\n";
      print_time "    SUNMatZero Time: %22.15e @\n @\n" (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNMatCopy Test
 * --------------------------------------------------------------------*)
let test_sunmatcopy a myid =
  let tol = 1e-15 in
  try
    let b = clone a in
    (* copy matrix data *)
    let start_time = get_time () in
    (try Matrix.blit ~src:a ~dst:b
     with _ -> printf ">>> FAILED test -- SUNMatZero failed on Proc %d @\n" myid;
               raise Exit);
    let stop_time = get_time () in

    (* check matrix entries *)
    if M.check_matrix b a tol then begin
      printf ">>> FAILED test -- SUNMatCopy check, Proc %d @\n" myid;
      print_time "    SUNMatCopy Time: %22.15e @\n @\n" (stop_time -. start_time);
      raise Exit
    end
    else if myid = 0 then begin
      (match Sundials.Config.sundials_version with
       | 2,_,_ | 3,1,0 | 3,1,1 -> printf "    PASSED test -- N_VConst @\n"
       | _ -> printf "    PASSED test -- SUNMatCopy @\n");
      print_time "    SUNMatCopy Time: %22.15e @\n @\n" (stop_time -. start_time);
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNMatScaleAdd Test: A = c * A + B
 * --------------------------------------------------------------------*)
let test_sunmatscaleadd a i myid =
  let tol = 1e-15 in
  try
    (*
     * Case 1: same sparsity/bandwith pattern
     *)

    (* protect A *)
    let b = clone a in
    (try Matrix.blit ~src:a ~dst:b
     with _ ->
       printf ">>> FAILED test -- SUNMatCopy failed on Proc %d @\n" myid;
       raise Exit);

    (* fill vector data *)
    let start_time = get_time () in
    (try Matrix.scale_add (-1.0) b b
     with _ ->
       printf ">>> FAILED test -- SUNMatScaleAdd failed on Proc %d @\n" myid;
       raise Exit);
    let stop_time = get_time () in

    (* check matrix entries *)
    if M.check_matrix_entry b 0.0 tol then begin
      printf ">>> FAILED test -- SUNMatScaleAdd case 1 check, Proc %d @\n" myid;
      print_time "    SUNMatScaleAdd Time: %22.15e @\n @\n"
        (stop_time -. start_time);
      raise Exit
    end
    else if (myid == 0) then begin
      printf "    PASSED test -- SUNMatScaleAdd case 1 @\n";
      print_time "    SUNMatScaleAdd Time: %22.15e @\n @\n"
        (stop_time -. start_time)
    end;

    (*
     * Case 2: different sparsity/bandwith patterns
     *)
    if M.is_square a then begin

      (* protect A and I *)
      let d = clone a in
      (try Matrix.blit ~src:a ~dst:d;
       with _ ->
         printf ">>> FAILED test -- SUNMatCopy failed on Proc %d @\n" myid;
         raise Exit);

      let c = clone i in
      (try Matrix.blit ~src:i ~dst:c
       with _ ->
         printf ">>> FAILED test -- SUNMatCopy failed on Proc %d @\n" myid;
         raise Exit);

      (* fill B and C *)
      let start_time = get_time () in
      (try Matrix.scale_add 1.0 d i
       with _ ->
         printf ">>> FAILED test -- SUNMatScaleAdd failed on Proc %d @\n" myid;
         raise Exit);
      (try Matrix.scale_add 1.0 c a with _ ->
         printf ">>> FAILED test -- SUNMatScaleAdd failed on Proc %d @\n" myid;
         raise Exit);
      let stop_time = get_time () in

      (* check matrix entries *)
      if M.check_matrix d c tol then begin
        printf ">>> FAILED test -- SUNMatScaleAdd case 2 check, Proc %d @\n" myid;
        print_time "    SUNMatScaleAdd Time: %22.15e @\n @\n"
          (stop_time -. start_time);
        raise Exit
      end
      else if myid = 0 then begin
        printf "    PASSED test -- SUNMatScaleAdd case 2 @\n";
        print_time "    SUNMatScaleAdd Time: %22.15e @\n @\n"
          (stop_time -. start_time);
      end
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNMatScaleAddI Tests
 * --------------------------------------------------------------------*)
let test_sunmatscaleaddi a i myid =
  let tol = 1e-15 in
  try
    (* protect A *)
    let b = clone a in
    (try Matrix.blit ~src:i ~dst:b
     with _ ->
       printf ">>> FAILED test -- SUNMatCopy failed on Proc %d @\n" myid;
       raise Exit);

    (* fill vector data *)
    let start_time = get_time () in
    (try Matrix.scale_addi (-1.0) b
     with _ ->
        printf ">>> FAILED test -- SUNMatScaleAddI failed on Proc %d @\n" myid;
        raise Exit);
    let stop_time = get_time () in

    (* check matrix *)
    if M.check_matrix_entry b 0.0 tol then begin
      printf ">>> FAILED test -- SUNMatScaleAddI check, Proc %d @\n" myid;
      print_time "    SUNMatScaleAddI Time: %22.15e @\n @\n"
        (stop_time -. start_time);
      raise Exit
    end
    else if (myid == 0) then begin
      printf "    PASSED test -- SUNMatScaleAddI @\n";
      print_time "    SUNMatScaleAddI Time: %22.15e @\n @\n"
        (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNMatMatvec Test (y should be correct A*x product)
 * --------------------------------------------------------------------*)
let test_sunmatmatvec a x y myid =
  let tol = 1e-14 in
  try
    (* harder tests for square matrices *)
    let failure, start_time, stop_time =
      if M.is_square a then begin
        (* protect A *)
        let b = clone a in
        (try Matrix.blit ~src:a ~dst:b
         with _ ->
           printf ">>> FAILED test -- SUNMatCopy failed on Proc %d @\n" myid;
           raise Exit);

        (* compute matrix vector product *)
        (try Matrix.scale_addi 3.0 b
         with _ ->
           printf ">>> FAILED test -- SUNMatScaleAddI failed on Proc %d @\n" myid;
           raise Exit);
        let z = NV.clone y in
        let w = NV.clone y in

        let start_time = get_time () in
        (try Matrix.matvec b x z
         with _ ->
            printf ">>> FAILED test -- SUNMatMatvec failed on Proc %d @\n" myid;
            raise Exit);
        let stop_time = get_time () in
        NV.linearsum 3.0 y 1.0 x w;
        M.check_vector w z tol, start_time, stop_time
      end
      else begin
        let z = NV.clone y in

        let start_time = get_time () in
        (try Matrix.matvec a x z
         with _ ->
           printf ">>> FAILED test -- SUNMatMatvec failed on Proc %d @\n" myid;
           raise Exit);
        let stop_time = get_time () in

        M.check_vector y z tol, start_time, stop_time
      end
    in
    if failure then begin
      printf ">>> FAILED test -- SUNMatMatvec check, Proc %d @\n" myid;
      print_time "    SUNMatMatvec Time: %22.15e @\n @\n"
        (stop_time -. start_time);
      raise Exit
    end
    else if myid = 0 then begin
      printf "    PASSED test -- SUNMatMatvec @\n";
      print_time "    SUNMatMatvec Time: %22.15e @\n @\n"
        (stop_time -. start_time)
    end;
    0
  with Exit -> 1

(* ----------------------------------------------------------------------
 * SUNMatSpace Test
 * --------------------------------------------------------------------*)
let test_sunmatspace ?cheat_leniw ?cheat_lenrw a myid =
  try
    let start_time = get_time () in
    let lenrw, leniw = try Matrix.space a with _ -> begin
        let stop_time = get_time () in
        printf ">>> FAILED test -- SUNMatSpace, Proc %d @\n" myid;
        print_time "    SUNMatSpace Time: %22.15e @\n @\n"
          (stop_time -. start_time);
        raise Exit
      end
    in
    let stop_time = get_time () in

    if myid = 0 then begin
      (* Hack around expected differences in the results of space *)
      let leniw =
        match cheat_leniw with
        | Some (actual, expected) when leniw = actual -> expected
        | _ -> leniw
      in
      let lenrw =
        match cheat_lenrw with
        | Some (actual, expected) when lenrw = actual -> expected
        | _ -> lenrw
      in
      if Sundials_impl.Version.lt500
      then printf "    PASSED test -- SUNMatSpace, lenrw = %d, leniw = %d@\n"
             lenrw leniw
      else printf "    PASSED test -- SUNMatSpace lenrw=%d leniw=%d@\n"
             lenrw leniw;
      print_time "    SUNMatSpace Time: %22.15e @\n @\n" (stop_time -. start_time);
    end;
    0
  with Exit -> 1

end

