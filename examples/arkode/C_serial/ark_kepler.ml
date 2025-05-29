(* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *-----------------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2025.
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * We consider the Kepler problem. We choose one body to be the center of our
 * coordinate system and then we use the coordinates q = (q1, q2) to represent
 * the position of the second body relative to the first (center). This yields
 * the ODE:
 *    dq/dt = [ p1 ]
 *            [ p2 ]
 *    dp/dt = [ -q1 / (q1^2 + q2^2)^(3/2) ]
 *          = [ -q2 / (q1^2 + q2^2)^(3/2) ]
 * with the initial conditions
 *    q(0) = [ 1 - e ],  p(0) = [        0          ]
 *           [   0   ]          [ sqrt((1+e)/(1-e)) ]
 * where e = 0.6 is the eccentricity.
 *
 * The Hamiltonian for the system,
 *    H(p,q) = 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2)
 * is conserved as well as the angular momentum,
 *    L(p,q) = q1*p2 - q2*p1.
 *
 * By default we solve the problem by letting y = [ q, p ]^T then using a 4th
 * order symplectic integrator via the SPRKStep time-stepper of ARKODE with a
 * fixed time-step size.
 *
 * The rootfinding feature of SPRKStep is used to count the number of complete orbits.
 * This is done by defining the function,
 *    g(q) = q2
 * and providing it to SPRKStep as the function to find the roots for g(q).
 *
 * The program also accepts command line arguments to change the method
 * used and time-stepping strategy. The program has the following CLI arguments:
 *
 *   --step-mode <fixed, adapt>  should we use a fixed time-step or adaptive time-step (default fixed)
 *   --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep with an ERK method (default SPRK)
 *   --method <string>           which method to use (default ARKODE_SPRK_MCLACHLAN_4_4)
 *   --use-compensated-sums      turns on compensated summation in ARKODE where applicable
 *   --disable-tstop             turns off tstop mode
 *   --dt <Real>                 the fixed-time step size to use if fixed time stepping is turned on (default 0.01)
 *   --tf <Real>                 the final time for the simulation (default 100)
 *   --nout                      number of output times
 *   --count-orbits              use rootfinding to count the number of completed orbits
 *   --check-order               compute the order of the method used and check if it is within the expected range
 *
 * References:
 *    Ernst Hairer, Christain Lubich, Gerhard Wanner
 *    Geometric Numerical Integration: Structure-Preserving
 *    Algorithms for Ordinary Differential Equations
 *    Springer, 2006,
 *    ISSN 0179-3632
 * --------------------------------------------------------------------------*)

open Sundials
module ARKStep = Arkode.ARKStep
module SPRKStep = Arkode.SPRKStep
module SPRKTable = Arkode.SPRKStep.MethodTable

let printf = Printf.printf
let eprintf = Printf.eprintf
let fprintf = Printf.fprintf
let sprintf = Printf.sprintf

let num_dt = 8

type step_mode = Fixed | Adapt
type stepper = Sprk | Erk

let step_mode        = ref Fixed
let stepper          = ref Sprk
let method_name      = ref ""
let count_orbits     = ref false
let use_compsums     = ref false
let use_tstop        = ref true
let dt               = ref 1e-2
let tf               = ref 100.0
let check_order      = ref false
let num_output_times = ref 50

let set_step_mode = function
  | "fixed" -> step_mode := Fixed
  | "adapt" -> step_mode := Adapt
  | _ -> raise Arg.(Bad "ERROR: --step-mode must be 'fixed' or 'adapt'")

let set_stepper = function
  | "SPRK" -> stepper := Sprk
  | "ERK"  -> stepper := Erk
  | _ -> raise Arg.(Bad "ERROR: --stepper must be 'SPRK' or 'ERK'")

let parse_args () =
  Arg.parse [
      ("--step-mode <fixed, adapt>", Arg.String set_step_mode,
       "should we use a fixed time-step or adaptive time-step (default fixed)");

      ("--stepper <SPRK, ERK>", Arg.String set_stepper,
       "should we use SPRKStep or ARKStep with an ERK method (default SPRK)");

      ("--method <string>", Arg.Set_string method_name,
       "which method to use (default ARKODE_SPRK_MCLACHLAN_4_4)");

      ("--use-compensated-sums", Arg.Set use_compsums,
       "turns on compensated summation in ARKODE where applicable");

      ("--disable-tstop", Arg.Clear use_tstop,
       "turns off tstop mode");

      ("--dt <Real>", Arg.Set_float dt,
       "the fixed-time step size to use if fixed time stepping is turned on (default 0.01)");

      ("--tf <Real>", Arg.Set_float tf,
       "the final time for the simulation (default 100)");

      ("--nout <int>", Arg.Set_int num_output_times,
       "the number of output times (default 50)");

      ("--count-orbits", Arg.Set count_orbits,
       "use rootfinding to count the number of completed orbits");

      ("--check-order", Arg.Set check_order,
       "compute the order of the method used and check if it is within range of the expected");
    ] (fun _ -> raise Arg.(Bad "unrecognized option"))
      ("ark_kepler: an ARKODE example demonstrating the SPRKStep "
       ^ "time-stepping module solving the Kepler problem\n");
    if !method_name = "" then
      method_name := (match !stepper with
                      | Sprk -> "ARKODE_SPRK_MCLACHLAN_4_4"
                      | Erk  -> "ARKODE_ZONNEVELD_5_3_4")

let print_args () =
  printf "Problem Arguments:\n";
  printf "  stepper:              %d\n" (match !stepper with Sprk -> 0 | Erk -> 1);
  printf "  step mode:            %d\n" (match !step_mode with Fixed -> 0 | Adapt -> 0);
  printf "  use tstop:            %d\n" (if !use_tstop then 1 else 0);
  printf "  use compensated sums: %d\n" (if !use_compsums then 1 else 0);
  printf "  dt:                   %g\n" !dt;
  printf "  Tf:                   %g\n" !tf;
  printf "  nout:                 %d\n\n" !num_output_times

type convergence_args = {
    a11 : float;
    a12 : float;
    a21 : float;
    a22 : float;
    b1  : float;
    b2  : float;
    b1e : float;
    b2e : float;
  }
let convergence_zero =
  { a11 = 0.; a12 = 0.; a21 = 0.; a22 = 0.;
     b1 = 0.;  b2 = 0.; b1e = 0.; b2e = 0. }

let compute_convergence (orders : RealArray.t) _expected_order
                        { a11; a12; a21; a22; b1; b2; _ } =
  (* Compute/print overall estimated convergence rate *)
  let ord_sum, ord_max =
    RealArray.fold_left (fun (s, m) o -> (s +. o,  max m o)) (0., 0.) orders
  in
  let ord_avg = ord_sum /. (float (RealArray.length orders)) in
  let det = a11 *. a22 -. a12 *. a21 in
  let ord_est = (a11 *. b2 -. a21 *. b1) /. det in
  ord_avg, ord_max, ord_est

type user_data = { ecc : float }

(* RHS callback functions *)

let velocity _udata _t (y : RealArray.t) (ydot : RealArray.t) =
  let p1 = y.{2} in
  let p2 = y.{3} in
  ydot.{0} <- p1;
  ydot.{1} <- p2

let force _udata _t (y : RealArray.t) (ydot : RealArray.t) =
  let q1       = y.{0} in
  let q2       = y.{1} in
  let sqrt_qTq = sqrt(q1 *. q1 +. q2 *. q2) in
  ydot.{2} <- -.q1 /. Float.pow sqrt_qTq 3.0;
  ydot.{3} <- -.q2 /. Float.pow sqrt_qTq 3.0

let dydt udata t (y : RealArray.t) (ydot : RealArray.t) =
  force udata t y ydot;
  velocity udata t y ydot

(* g(q) callback function for rootfinding *)

let rootfn _udata _t y gout =
  let q2 = y.{1}  in
  gout.{0} <- q2

(* Helper functions *)

let initial_conditions (y0 : RealArray.t) ecc =
  y0.{0} <- 1.0 -. ecc;
  y0.{1} <- 0.0;
  y0.{2} <- 0.0;
  y0.{3} <- sqrt((1.0 +. ecc) /. (1.0 -. ecc))

let hamiltonian (y : RealArray.t) =
  let sqrt_qTq = sqrt(y.{0} *. y.{0} +. y.{1} *. y.{1}) in
  let pTp = y.{2} *. y.{2} +. y.{3} *. y.{3} in
  0.5 *. pTp -. 1.0 /. sqrt_qTq

let angular_momentum y =
  let q1 = y.{0} in
  let q2 = y.{1} in
  let p1 = y.{2} in
  let p2 = y.{3} in
  q1 *. p2 -. q2 *. p1

type solve_args = {
    count_orbits : bool;
    step_mode : step_mode;
    stepper : stepper;
    use_compsums : bool;
    num_output_times : int;
    method_name : string;
    dt : float;
    tf : float;
  }

let solve_problem { count_orbits; step_mode; stepper;
                    use_compsums; num_output_times;
                    method_name; dt; tf } context sol =

  (* Default problem parameters *)
  let t0    = 0.0 in
  let dTout = (tf -. t0) /. (float num_output_times) in
  let ecc   = 0.6 in

  printf "\n   Begin Kepler Problem\n\n";
  print_args ();

  (* Allocate and fill udata structure *)
  let udata = { ecc } in

  (* Allocate our state vector *)
  let y = Nvector_serial.make ~context 4 0.0 in
  let ydata = Nvector.unwrap y in

  (* Fill the initial conditions *)
  initial_conditions ydata ecc;

  (* Create SPRKStep integrator *)
  let session =
    match stepper with
    | Sprk -> begin
       if step_mode = Adapt then begin
          eprintf "ERROR: adaptive time-steps are not supported with SPRKStep\n";
          exit 1
       end;
       (* Optional: enable temporal root-finding *)
       let roots = if count_orbits then Some (1, rootfn udata) else None in
       let arkode_mem = SPRKStep.init ~context ~step:dt
                                      ~f1:(force udata) ~f2:(velocity udata)
                                      ?roots t0 y
       in
       SPRKStep.set_method_name arkode_mem method_name;
       SPRKStep.set_use_compensated_sums arkode_mem use_compsums;
       SPRKStep.set_max_num_steps arkode_mem (int_of_float (ceil(tf /. dt)) + 1);
       Arkode.SPRK arkode_mem
    end
    | Erk -> begin
       (* Optional: enable temporal root-finding *)
       let roots = if count_orbits then Some (1, rootfn udata) else None in
       let arkode_mem = ARKStep.(init ~context (explicit (dydt udata))
                                      (SStolerances (dt, dt)) ?roots t0 y)
       in
       ARKStep.set_table_name arkode_mem
         ~itable:"ARKODE_DIRK_NONE" ~etable:method_name ();
       ARKStep.set_max_num_steps arkode_mem (int_of_float (ceil(tf /. dt)) + 1);
       (if step_mode = Fixed then ARKStep.set_fixed_step arkode_mem (Some dt));
       Arkode.ARK arkode_mem
    end
  in

  (* Open output files *)
  let conserved_fp = Logfile.openfile ~trunc:true
                       (sprintf "ark_kepler_conserved_%s-dt-%.2e.txt" method_name dt)
  in
  let solution_fp = Logfile.openfile ~trunc:true
                      (sprintf "ark_kepler_solution_%s-dt-%.2e.txt" method_name dt)
  in
  let times_fp = Logfile.openfile ~trunc:true
                   (sprintf "ark_kepler_times_%s-dt-%.2e.txt" method_name dt)
  in

  (* Print out starting energy, momentum before integrating *)
  let tret = t0 in
  let tout = t0 +. dTout in
  let h0 = hamiltonian ydata in
  let l0 = angular_momentum ydata in
  printf "t = %.4f, H(p,q) = %.16f, L(p,q) = %.16f\n" tret h0 l0;
  Logfile.output_string times_fp (sprintf "%.16f\n" tret);
  Logfile.output_string conserved_fp (sprintf "%.16f, %.16f\n" h0 l0);
  Nvector.Ops.print ~logfile:solution_fp y;

  let rootsfound = Roots.create 1 in
  let loop set_stop_time evolve_normal get_root_info =
    let rec go iout tout num_orbits =
      if iout >= num_output_times then flush stdout
      else begin
        (* Optional: if the stop time is not set, then its possible that the
           exact requested output time will not be hit (even with a fixed
           time-step due to roundoff error accumulation) and interpolation will be
           used to get the solution at the output time. *)
        (if !use_tstop then set_stop_time tout);
        let tret, retval = evolve_normal tout y in
        (match retval with
         | Arkode.Common.RootsFound -> begin
             let num_orbits = num_orbits +. 0.5 in

             printf "ROOT RETURN:\t";
             get_root_info rootsfound;
             printf "  g.{0} = %3d, y.{0} = %3g, y.{1} = %3g, num. orbits is now %.2f\n"
                     (if Roots.detected rootsfound 0 then 1 else 0)
                     ydata.{0} ydata.{1} num_orbits;
             printf "t = %.4f, H(p,q)-H0 = %.16e, L(p,q)-L0 = %.16e\n" tret
                     (hamiltonian ydata -. h0) (angular_momentum ydata  -. l0);

             go iout tout num_orbits
           end
         | Arkode.Common.Success | Arkode.Common.StopTimeReached -> begin
             (* Output current integration status *)
             printf "t = %.4f, H(p,q)-H0 = %.16e, L(p,q)-L0 = %.16e\n"
                 tret (hamiltonian ydata -. h0)
                      (angular_momentum ydata -. l0);
             Logfile.output_string times_fp (sprintf "%.16f\n" tret);
             Logfile.output_string conserved_fp
               (sprintf "%.16f, %.16f\n" (hamiltonian ydata)
                                         (angular_momentum ydata));
             Nvector.Ops.print ~logfile:solution_fp y;

             go (iout + 1) (min (tout +. dTout) tf) num_orbits
         end)
      end
    in
    go 0 tout 0.
  in

  (* Do integration *)
  (match session with
   | Arkode.SPRK arkode_mem ->
       loop (SPRKStep.set_stop_time arkode_mem)
            (SPRKStep.evolve_normal arkode_mem)
            (SPRKStep.get_root_info arkode_mem);
       SPRKStep.print_all_stats arkode_mem OutputTable
   | Arkode.ARK arkode_mem ->
       loop (ARKStep.set_stop_time arkode_mem)
            (ARKStep.evolve_normal arkode_mem)
            (ARKStep.get_root_info arkode_mem);
       ARKStep.print_all_stats arkode_mem OutputTable
  | _ -> assert false);

  Nvector.Ops.scale 1.0 y sol;
  (* return the energy_error *)
  hamiltonian ydata -. h0

let main () =
  (* Create the SUNDIALS context object for this simulation *)
  let context = Sundials.Context.make () in

  (* Parse the command line arguments *)
  parse_args ();
  let solve_args = {
    count_orbits     = !count_orbits;
    step_mode        = !step_mode;
    stepper          = !stepper;
    use_compsums     = !use_compsums;
    num_output_times = !num_output_times;
    method_name      = !method_name;
    dt               = !dt;
    tf               = !tf;
  } in

  (* Allocate space for result variables *)
  let sol = Nvector_serial.make ~context 4 0.0 in

  if not !check_order then
    (* SolveProblem calls a stepper to evolve the problem to Tf *)
    ignore (solve_problem solve_args context sol)
  else begin
    (* Compute the order of accuracy of the method by testing
       it with different step sizes. *)
    let acc_orders = RealArray.make num_dt 0.0 in
    let con_orders = RealArray.make num_dt 0.0 in
    let acc_errors = RealArray.make num_dt 0.0 in
    let con_errors = RealArray.make num_dt 0.0 in
    let SPRKTable.{ method_order = expected_order; _ } =
      Option.get (SPRKTable.load_by_name !method_name) in

    let ref_sol = Nvector.clone sol in
    let error = Nvector.clone sol in
    let refine = 0.5 in
    let dt = if expected_order >= 3 then 1e-1 else 1e-3 in

    (* Create a reference solution using 8th order ERK with a small time step *)
    (* SolveProblem calls a stepper to evolve the problem to Tf *)
    ignore (solve_problem
              { solve_args with dt = 1e-3;
                                step_mode = Fixed;
                                stepper = Erk;
                                method_name = "ARKODE_ARK548L2SAb_ERK_8_4_5" }
              context ref_sol);
    let dts =
      RealArray.init num_dt (fun i -> dt *. Float.pow refine (float i))
    in
    let rec f i ({ a11; a12; a21; a22; b1; b2; b1e; b2e } as ab) =
      if i = num_dt then ab
      else begin
        (* Compute the error with various step sizes *)
        (* Set the dt to use for this solve *)
        (* SolveProblem calls a stepper to evolve the problem to Tf *)
        let energy_error = solve_problem { solve_args with dt = dts.{i} } context sol in
        printf "\n";

        (* Compute the error *)
        Nvector.Ops.linearsum 1.0 sol (-1.0) ref_sol error;
        acc_errors.{i} <-
          sqrt(Nvector.Ops.dotprod error error)
          /. (float (Nvector.Ops.getlength error));
        con_errors.{i} <- Float.abs(energy_error);


        if i >= 1 then begin
          acc_orders.{i - 1} <- log (acc_errors.{i} /. acc_errors.{i - 1})
                                /. log (dts.{i} /. dts.{i - 1});
          con_orders.{i - 1} <- log (con_errors.{i} /. con_errors.{i - 1})
                                /. log (dts.{i} /. dts.{i - 1})
        end;
        f (i + 1) { a11 = a11 +. 1.0;
                    a12 = a12 +. log dts.{i};
                    a21 = a21 +. log dts.{i};
                    a22 = a22 +. (log dts.{i} *. log dts.{i});
                    b1  =  b1 +. log acc_errors.{i};
                    b2  =  b2 +. (log acc_errors.{i} *. log dts.{i});
                    b1e = b1e +. log con_errors.{i};
                    b2e = b2e +. (log con_errors.{i} *. log dts.{i}) }
      end
    in
    let ab = f 0 convergence_zero in

    (* Compute the order of accuracy *)
    let ord_avg, ord_max_acc, ord_est =
      compute_convergence acc_orders expected_order ab
    in
    printf "Order of accuracy wrt solution:    expected = %d, max = %.4f,  avg = %.4f,  overall = %.4f\n"
           expected_order ord_max_acc ord_avg ord_est;

    (* Compute the order of accuracy with respect to conservation *)
    let ord_avg, ord_max_conv, ord_est =
      compute_convergence con_orders expected_order
                          { ab with b1 = ab.b1e; b2 = ab.b2e; }
    in
    printf "Order of accuracy wrt Hamiltonian: expected = %d, max = %.4f,  avg = %.4f,  overall = %.4f\n"
           expected_order ord_max_conv ord_avg ord_est;

    if ord_max_acc < (float expected_order -. 0.5) then begin
      printf ">>> FAILURE: computed order of accuracy wrt solution is below expected (%d)\n"
             expected_order;
      exit 1
    end;

    if ord_max_conv < (float expected_order -. 0.5) then begin
      printf ">>> FAILURE: computed order of accuracy wrt Hamiltonian is below expected (%d)\n"
             expected_order;
      exit 1
    end
  end

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

