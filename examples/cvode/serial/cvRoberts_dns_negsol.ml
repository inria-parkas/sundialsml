(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/29 22:21:29 $
 * -----------------------------------------------------------------
 * Programmers: Radu Serban and Alan Hindmarsh @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Feb 2019.
 * -----------------------------------------------------------------
 * Modification of the CVODE example cvRoberts_dns to illustrate
 * the treatment of unphysical solution components through the RHS
 * function return flag.
 *
 * Note that, to make possible negative solution components, the
 * absolute tolerances had to be loosened a bit from their values
 * in cvRoberts_dns.
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * -----------------------------------------------------------------
 *)

open Sundials

let unwrap = Nvector.unwrap

let printf = Printf.printf

let ith (v : RealArray.t) i = v.{i - 1}
let set_ith (v : RealArray.t) i e = v.{i - 1} <- e

(* Problem Constants *)

let neq    = 3        (* number of equations  *)
let y1     = 1.0      (* initial y components *)
let y2     = 0.0
let y3     = 0.0
let rtol   = 1.0e-4   (* scalar relative tolerance            *)
let atol1  = 1.0e-7   (* vector absolute tolerance components *)
let atol2  = 1.0e-13
let atol3  = 1.0e-5
let t0     = 0.0      (* initial time           *)
let t1     = 0.4      (* first output time      *)
let tmult  = 10.0     (* output time factor     *)
let nout   = 14       (* number of output times *)
let nroots = 2        (* number of root functions *)

let f check_negative t (y : RealArray.t) (yd : RealArray.t) =
  if !check_negative && (y.{0} < 0.0 || y.{1} < 0.0 || y.{2} < 0.0)
  then raise RecoverableFailure;
  let yd1 = -0.04 *. y.{0} +. 1.0e4 *. y.{1} *. y.{2}
  and yd3 = 3.0e7 *. y.{1} *. y.{1}
  in
  yd.{0} <- yd1;
  yd.{1} <- (-. yd1 -. yd3);
  yd.{2} <- yd3

let print_output =
  printf "At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n"

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
  in
  printf "\nFinal Statistics:\n";
  printf "nst = %-6d nfe  = %-6d nsetups = %-6d nfeLS = %-6d nje = %d\n"
    nst nfe nsetups nfeLS nje;
  printf "nni = %-6d ncfn = %-6d netf = %-6d\n \n"
    nni ncfn netf

let main () =
  (* Create serial vector of length NEQ for I.C. and abstol *)
  let y = Nvector_serial.make neq 0.0
  and abstol = RealArray.create neq
  in
  let ydata = unwrap y in

  (* Initialize y *)
  set_ith ydata 1 y1;
  set_ith ydata 2 y2;
  set_ith ydata 3 y3;

  (* Set the vector absolute tolerance *)
  set_ith abstol 1 atol1;
  set_ith abstol 2 atol2;
  set_ith abstol 3 atol3;

  (* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration *)
  (* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. *)
  (* Call CVodeRootInit to specify the root function g with 2 components *)
  (* Call CVDense to specify the CVDENSE dense linear solver *)
  let m = Matrix.dense neq in
  let check_negative = ref false in
  let cvode_mem =
    Cvode.(init BDF ~lsolver:Dls.(solver (dense y m))
                (SVtolerances (rtol, (Nvector_serial.wrap abstol)))
                (f check_negative) t0 y)
  in

  (* Case 1: ignore negative solution components *)
  printf "Ignore negative solution components\n\n";
  (* In loop, call CVode in CV_NORMAL mode *)

  let tout = ref t1
  and iout = ref 0
  in
  while (!iout <> nout) do
    let (t, flag) = Cvode.solve_normal cvode_mem !tout y in
    print_output t (ith ydata 1) (ith ydata 2) (ith ydata 3);
    iout := !iout + 1;
    tout := !tout *. tmult
  done;

  (* Print some final statistics *)
  print_final_stats cvode_mem;

  (* Case 2: intercept negative solution components *)
  printf "Intercept negative solution components\n\n";
  check_negative := true;
  (* Reinitialize solver *)
  set_ith ydata 1 y1;
  set_ith ydata 2 y2;
  set_ith ydata 3 y3;
  Cvode.reinit cvode_mem t0 y;
  (* In loop, call CVode in CV_NORMAL mode *)

  let tout = ref t1
  and iout = ref 0
  in
  while (!iout <> nout) do
    let (t, flag) = Cvode.solve_normal cvode_mem !tout y in
    print_output t (ith ydata 1) (ith ydata 2) (ith ydata 3);
    iout := !iout + 1;
    tout := !tout *. tmult
  done;

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
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
