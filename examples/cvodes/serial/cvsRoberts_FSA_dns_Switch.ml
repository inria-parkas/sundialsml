(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2010/12/01 22:58:00 $
 * -----------------------------------------------------------------
 * Programmers: Radu Serban and Alan Hindmarsh, and Cody Balos @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Feb 2019.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES. The problem is from chemical
 * kinetics, and consists of the following three rate equations:
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions y1 = 1.0, y2 = y3 = 0. The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODES dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 *)

open Sundials

module Sens = Cvodes.Sensitivity
module Densemat = Matrix.Dense
let unwrap = Nvector.unwrap

let clone = Nvector_serial.Ops.clone
let const = Nvector_serial.Ops.const

let printf = Printf.printf

(* Problem Constants *)

let mxsteps = 2000    (* output time factor *)
let neq     = 3       (* number of equations  *)
let t0      = 0.0     (* initial time *)
let t1      = 4.0e10  (* first output time *)
let zero  = 0.0

(* Type : UserData *)

type user_data = {
  mutable sensi   : bool; (* turn on (T) or off (F) sensitivity analysis    *)
  mutable errconS : bool; (* full (T) or partial error control (F)          *)
  mutable fsDQ    : bool; (* user provided r.h.s sensitivity analysis (T/F) *)
  mutable meth    : (RealArray.t, Nvector_serial.kind) Sens.sens_method;
                          (* sensitivity method                 *)
  p               : RealArray.t
}

(* Prototypes of functions by CVODES *)

(* f routine. Compute f(t,y). *)

let f data _ (y : RealArray.t) (ydot : RealArray.t) =
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2} in

  let yd1 = -.p1*.y.{0} +. p2*.y.{1}*.y.{2} in
  let yd3 = p3*.y.{1}*.y.{1} in
  ydot.{0} <- yd1;
  ydot.{1} <- (-.yd1 -. yd3);
  ydot.{2} <- yd3

(* Jacobian routine. Compute J(t,y). *)

let jac data { Cvode.jac_y = (y : RealArray.t) } jmat =
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}
  in
  let set = Matrix.Dense.set jmat in
  set 0 0 (-.p1);
  set 0 1 (p2*.y.{2});
  set 0 2 (p2*.y.{1});
  set 1 0 ( p1);
  set 1 1 (-.p2*.y.{2}-.2.0*.p3*.y.{1});
  set 1 2 (-.p2*.y.{1});
  set 2 1 (2.0*.p3*.y.{1})

(* fS routine. Compute sensitivity r.h.s. *)

let fS : user_data -> RealArray.t Sens.sensrhsfn1 =
  fun data iS { Sens.y = y } yS ySdot ->
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}
  and s1 = yS.{0}
  and s2 = yS.{1}
  and s3 = yS.{2}
  in
  let sd1 = -.p1*.s1 +. p2*.y.{2}*.s2 +. p2*.y.{1}*.s3 in
  let sd3 = 2.0*.p3*.y.{1}*.s2 in
  let sd2 = -.sd1-.sd3
  in
  let sd1, sd2, sd3 =
    (match iS with
     | 0 -> (sd1 -. y.{0}, sd2 +. y.{0}, sd3)
     | 1 -> (sd1 +. y.{1}*.y.{2}, sd2 -. y.{1}*.y.{2}, sd3)
     | 2 -> (sd1, sd2 -. y.{1}*.y.{1}, sd3 +. y.{1}*.y.{1})
     | _ -> assert false);
  in
  ySdot.{0} <- sd1;
  ySdot.{1} <- sd2;
  ySdot.{2} <- sd3

(* Print current t, step count, order, stepsize, and solution. *)

let print_output s t udata =
  let nst = Cvode.get_num_steps s
  and qu  = Cvode.get_last_order s
  and hu  = Cvode.get_last_step s
  in
  printf "%8.3e %2d  %8.3e %5d\n" t qu hu nst;
  print_string "                  Solution       ";
  printf "%12.4e %12.4e %12.4e \n" udata.{0} udata.{1} udata.{2}

(* Print sensitivities. *)

let print_output_s uS =
  let sdata = unwrap uS.(0) in
  print_string "                  Sensitivity 1  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2};
  let sdata = unwrap uS.(1) in
  print_string "                  Sensitivity 2  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2};
  let sdata = unwrap uS.(2) in
  print_string "                  Sensitivity 3  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2}

(* Print some final statistics from the CVODES memory. *)
(* For high NUM_REPS, the cost of OCaml printf becomes important! *)

let print_string_5d s i =
  print_string s;
  if i < 10 then print_string "    "
  else if i < 100 then print_string "   "
  else if i < 1000 then print_string "  "
  else if i < 10000 then print_string " ";
  print_int i

let print_final_stats s data =
  let open Cvode in
  let nst     = get_num_steps s
  and nfe     = get_num_rhs_evals s
  and nsetups = get_num_lin_solv_setups s
  and netf    = get_num_err_test_fails s
  and nni     = get_num_nonlin_solv_iters s
  and ncfn    = get_num_nonlin_solv_conv_fails s
  in
  print_string "Run statistics:\n";
  print_string_5d "   nst     = " nst;
  print_string_5d "\n   nfe     = " nfe;
  print_string_5d "\n   netf    = " netf;
  print_string_5d "    nsetups  = " nsetups;
  print_string_5d "\n   nni     = " nni;
  print_string_5d "    ncfn     = " ncfn;

  let njeD = Dls.get_num_jac_evals s
  and nfeD = Dls.get_num_lin_rhs_evals s
  in
  print_string_5d "\n   njeD    = " njeD;
  print_string_5d "    nfeD     = " nfeD;
  print_newline ();

  if data.sensi then begin
    let nfSe     = Sens.get_num_rhs_evals s
    and nfeS     = Sens.get_num_rhs_evals_sens s
    and nsetupsS = Sens.get_num_lin_solv_setups s
    and netfS    = if data.errconS then Sens.get_num_err_test_fails s else 0
    and nniS, ncfnS =
      match data.meth with
      | Sens.Staggered _ -> Sens.get_num_nonlin_solv_iters s,
                            Sens.get_num_nonlin_solv_conv_fails s
      | Sens.Staggered1 _ | Sens.Simultaneous _ -> 0,0
    in
    print_string "   -----------------------------------\n";
    print_string_5d "   nfSe    = " nfSe;
    print_string_5d "    nfeS     = " nfeS;
    print_string_5d "\n   netfs   = " netfS;
    print_string_5d "    nsetupsS = " nsetupsS;
    print_string_5d "\n   nniS    = " nniS;
    print_string_5d "    ncfnS    = " ncfnS;
    print_newline ()
  end

let print_header data =
  (* Print sensitivity control flags *)
  printf "Sensitivity: ";
  if data.sensi then begin
    printf "YES (%s + %s + %s)\n"
      (match data.meth with
       | Sens.Simultaneous _ -> "SIMULTANEOUS"
       | Sens.Staggered _    -> "STAGGERED"
       | Sens.Staggered1 _   -> "STAGGERED-1")
      (if data.errconS then "FULL ERROR CONTROL" else "PARTIAL ERROR CONTROL")
      (if data.fsDQ    then "DQ sensitivity RHS"
                       else "user-provided sensitivity RHS")
  end else printf "NO\n";

  (* Print current problem parameters *)
  let p = data.p in
  printf "Parameters: [%8.4e  %8.4e  %8.4e]\n" p.{0} p.{1} p.{2}

(* Runs integrator and prints final statistics when complete. *)
let run_cvode data cvode_mem y =
  (* Print header for current run *)
  print_header data;

  (* Call CVode in CV_NORMAL mode *)
  ignore (Cvode.solve_normal cvode_mem t1 y);

  (* Print final statistics *)
  print_final_stats cvode_mem data;
  printf "\n"

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* User data structure *)
  let data = {
    sensi   = true;                   (* sensitivity ON                *)
    meth    = Sens.Simultaneous None; (* simultaneous corrector method *)
    errconS = true;                   (* full error control            *)
    fsDQ    = false;                  (* user-provided sensitvity RHS  *)
    p       = RealArray.of_list [ 0.04; 1.0e4; 3.0e7 ]
  } in

  (* Initial conditions *)
  let y0data = RealArray.of_list [1.0; 0.0; 0.0] in
  let y0 = Nvector_serial.wrap y0data in
  let y = Nvector_serial.make neq 0.0 in

  (* Set integration tolerances *)
  let reltol = 1e-6 in
  let abstol = Nvector_serial.wrap (RealArray.of_list [1e-8; 1e-14; 1e-6]) in

  (* Create CVODES object *)
  let m = Matrix.dense neq in
  let cvode_mem =
    Cvode.(init BDF
                (SVtolerances (reltol, abstol))
                ~lsolver:Dls.(solver ~jac:(jac data) (dense y m))
                (f data) t0 y0)
  in
  Cvode.set_max_num_steps cvode_mem mxsteps;

  let ns = 3 in
  let pbar = Sundials.RealArray.copy data.p in
  let plist = Array.init ns (fun is -> is) in

  let yS0 = Array.init ns (fun _ -> let v = clone y in const 0.0 v; v) in

  let fSdata = Sens.OneByOne (Some (fS data)) in

  let sens_params = Sens.({
      pvals = Some data.p;
      pbar  = Some pbar;
      plist = Some plist;
    }) in

  (*
    Sensitivities are enabled
    Set full error control
    Set user-provided sensitivity RHS
    Run CVODES
  *)
  Sens.(init cvode_mem EEtolerances data.meth ~sens_params fSdata yS0);
  Sens.set_err_con cvode_mem data.errconS;
  run_cvode data cvode_mem y;

  (*
    Change parameters
    Toggle sensitivities OFF
    Reinitialize and run CVODES
  *)
  data.p.{0} <- 0.05;
  data.p.{1} <- 2.0e4;
  data.p.{2} <- 2.9e7;

  data.sensi <- false;
  Cvode.reinit cvode_mem t0 y0;
  Sens.toggle_off cvode_mem;
  run_cvode data cvode_mem y;

  (*
    Change parameters
    Switch to internal DQ sensitivity RHS function
    Toggle sensitivities ON (reinitialize sensitivities)
    Reinitialize and run CVODES
  *)
  data.p.{0} <- 0.06;
  data.p.{1} <- 3.0e4;
  data.p.{2} <- 2.8e7;

  data.sensi <- true;
  data.fsDQ  <- true;

  Cvode.reinit cvode_mem t0 y0;
  Sens.(init cvode_mem EEtolerances data.meth (OneByOne None) yS0);
  run_cvode data cvode_mem y;

  (*
    Switch to partial error control
    Switch back to user-provided sensitivity RHS
    Toggle sensitivities ON (reinitialize sensitivities)
    Change method to staggered
    Reinitialize and run CVODES
  *)

  data.sensi   <- true;
  data.errconS <- false;
  data.fsDQ    <- false;
  data.meth    <- Sens.Staggered None;

  Cvode.reinit cvode_mem t0 y0;
  Sens.set_err_con cvode_mem data.errconS;
  Sens.(init cvode_mem EEtolerances data.meth fSdata yS0);
  run_cvode data cvode_mem y;

  (*
    Free sensitivity-related memory
    (CVodeSensToggle is not needed, as CVodeSensFree toggles sensitivities OFF)
    Reinitialize and run CVODES
  *)
  data.sensi <- false;
  Sens.turn_off cvode_mem;
  Cvode.reinit cvode_mem t0 y0;
  run_cvode data cvode_mem y

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
