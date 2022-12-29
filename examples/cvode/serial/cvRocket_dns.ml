(* -----------------------------------------------------------------
 * Programmer(s): Alan C. Hindmarsh @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Dec 2022.
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security and
 * Southern Methodist University. All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding needed
 * for its solution by CVODE. The problem is a simpliflied model of a
 * rocket, ascending vertically, with mass decreasing over time. The
 * system (of size 2) is given by
 *    y_1 = rocket height H, y_1(0) = 0,
 *    y_2 = rocket velocity v, y_2(0) = 0,
 *    dH/dt = v,
 *    dv/dt = a(t,v).
 * The upward acceleration a(t,v) is given by
 *    a(t,v) = F/(M_r + M_f) - Dv - g,
 * where F = engine thrust force (constant) M_r = rocket mass without
 * fuel, M_f = fuel mass = M_f0 - r*t, r = fuel burn rate,
 * D = drag coefficient, g = gravitational acceleration.
 * The engine force is reset to 0 when the fuel mass reaches 0, or
 * when H reaches a preset height H_c, whichever happens first.
 * Rootfinding is used to locate the time at which M_f = 0 or H =
 * H_c, and also the time at which the rocket reaches its maximum
 * height, given by the condition v = 0, t > 0.
 *
 * The problem is solved with the BDF method and Dense linear solver.
 *
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*)

open Sundials
let printf = Printf.printf

(* Problem Constants *)

let neq = 2 (* number of equations  *)

let force  = 2200.0 (* engine force *)
let massr  = 10.0   (* rocket mass (empty) *)
let massf0 = 1.0    (* initial fuel mass *)
let brate  = 0.1    (* fuel burn rate *)
let drag   = 0.3    (* Drag coefficient *)
let grav   = 32.0   (* acceleration due to gravity *)
let hcut   = 4000.0 (* height of engine cutoff *)

let y1     = 0.0    (* initial y components *)
let y2     = 0.0
let rtol   = 1.0e-5 (* scalar relative tolerance            *)
let atol1  = 1.0e-2 (* vector absolute tolerance components *)
let atol2  = 1.0e-1
let t0     = 0.0    (* initial time           *)
let t1     = 1.0    (* first output time      *)
let tinc   = 1.0    (* output time increment  *)
let nout   = 70     (* number of output times *)

let zero   = 0.0

type user_data = {
  mutable engine_on : bool;
}

(*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 *)

(* f routine. Compute function f(t,y). *)
let f ud t (y : RealArray.t) (ydot : RealArray.t) =
  let v = y.{1} in
  let acc =
    if ud.engine_on then force /. (massr +. massf0 -. brate *. t) else zero
  in
  ydot.{0} <- v;
  ydot.{1} <- acc -. drag *. v -. grav

(* Jacobian routine. Compute J(t,y) = df/dy. *)
let jac _ jmat =
  let jmatdata = Matrix.Dense.unwrap jmat in
  jmatdata.{0,1} <- 1.0;
  jmatdata.{1,1} <- -. drag

(* g routine. Compute functions g_i(t,y). *)
let g ud t (y : RealArray.t) (gout : RealArray.t) =
  if ud.engine_on then begin
    gout.{0} <- massf0 -. brate *. t;
    let h = y.{0} in
    gout.{1} <- h -. hcut
  end else begin
    let v = y.{1} in
    gout.{0} <- v
  end

(*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 *)

let print_output t y1 y2 =
  printf "At t = %0.4e      y =%14.6e  %14.6e\n" t y1 y2

let print_root_info root_f1 root_f2 numroot =
  if numroot = 2 then printf "    rootsfound[] = %3d %3d\n" root_f1 root_f2;
  if numroot = 1 then printf "    rootsfound[] = %3d\n" root_f1

(* Get and print some final statistics *)
let print_final_stats s =
  let open Cvode in
  let nst     = get_num_steps s
  and nfe     = get_num_rhs_evals s
  and nsetups = get_num_lin_solv_setups s
  and netf    = get_num_err_test_fails s
  and nni     = get_num_nonlin_solv_iters s
  and ncfn    = get_num_nonlin_solv_conv_fails s
  and nje     = Dls.get_num_jac_evals s
  and nfeLS   = Dls.get_num_lin_rhs_evals s
  and nge     = get_num_g_evals s
  in
  printf "\nFinal Statistics:\n";
  printf "nst = %-6d nfe  = %-6d nsetups = %-6d nfeLS = %-6d nje = % d\n"
         nst nfe nsetups nfeLS nje;
  printf "nni = %-6d ncfn = %-6d netf = %-6d nge = %d\n \n" nni ncfn netf nge

(*
 *-------------------------------
 * Main Program
 *-------------------------------
 *)

let main () =
  (* Create the SUNDIALS context *)
  let context = Sundials.Context.make () in

  (* Create serial vector of length NEQ for I.C. and abstol *)
  let y = Nvector_serial.make ~context neq 0.0 in
  (* Initialize y *)
  let ydata = Nvector.unwrap y in
  ydata.{0} <- y1;
  ydata.{1} <- y2;

  (* Set the vector absolute tolerance *)
  let abstol = Nvector_serial.wrap ~context (RealArray.of_list [atol1; atol2]) in
  (* Set the scalar relative tolerance *)
  let reltol = rtol in

  (* Provide sunbooleantype engine_on as user data for use in f and g routines *)
  let user_data = { engine_on = true } in

  (* Create dense SUNMatrix for use in linear solves *)
  let a = Matrix.dense ~context neq in

  (* Create dense SUNLinearSolver object for use by CVode *)
  (* Set the user-supplied Jacobian routine Jac *)
  let lsolver = Cvode.Dls.solver ~jac (LinearSolver.Direct.dense ~context y a) in

  (* Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula *)
  (* Call CVodeInit to initialize the integrator memory and specify the right-hand side function in
   * y'=f(t,y), the inital time T0, and the initial dependent variable vector y. *)
  (* Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tolerances *)
  (* Call CVodeRootInit to specify the root function g with 2 components *)
  (* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode *)
  let cvode_mem =
    Cvode.(init BDF ~lsolver (SVtolerances (reltol, abstol))
                (f user_data) ~roots:(2, g user_data) t0 y)
  in

  (* In loop, call CVode, print results, check for root stops, and test for error.  On the first
     root return, restart with engine turned off. Break out of loop when NOUT preset output times
     have been reached, or when the returned value of H is negative.  *)
  printf " \nAccelerating rocket problem\n\n";

  let rootsfound = Roots.create 2 in
  let r i = Roots.(int_of_root (get rootsfound i)) in

  let rec iterate iout tout numroot =
    if iout = nout || ydata.{0} < zero then ()
    else
      let t, flag = Cvode.solve_normal cvode_mem tout y in
      print_output t ydata.{0} ydata.{1};
      match flag with
      | Cvode.RootsFound when user_data.engine_on -> (* engine cutoff *)
          Cvode.get_root_info cvode_mem rootsfound;
          print_root_info (r 0) (r 1) numroot;
          user_data.engine_on <- false;

          (* Call CVodeRootInit to specify the root function g with 1 component *)
          (* Reinitialize the solver with current t and y values. *)
          Cvode.reinit cvode_mem ~roots:(1, g user_data) t y;
          iterate (iout + 1) (tout +. tinc) 1

      | Cvode.RootsFound -> (* max.  height *)
          Cvode.get_root_info cvode_mem rootsfound;
          print_root_info (r 0) (r 1) numroot;
          iterate iout tout numroot

      | Cvode.Success ->
          iterate (iout + 1) (tout +. tinc) numroot

      | Cvode.StopTimeReached -> ()
    in
  iterate 0 t1 2;

  (* Print some final statistics *)
  print_final_stats cvode_mem

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

