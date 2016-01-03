(*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2010/12/01 23:02:23 $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL  (based on idaRoberts_dns.c)
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Aug 2014.
 *             Timothy Bourke, Inria, Jan 2016.
 * -----------------------------------------------------------------
 * This simple example problem for IDA, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:

 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.
 *
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01.
 *
 * The problem is solved with IDA using IDAKLU for the linear
 * solver, with a user-supplied Jacobian. Output is printed at
 * t = .4, 4, 40, ..., 4e10.
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray
module Roots = Sundials.Roots

let printf = Printf.printf

(* Problem Constants *)

let neq    = 3        (* number of equations  *)
let y1     = 1.0      (* initial y components *)
let y2     = 0.0
let y3     = 0.0
let rtol   = 1.0e-4   (* scalar relative tolerance            *)
let atol1  = 1.0e-8   (* vector absolute tolerance components *)
let atol2  = 1.0e-14
let atol3  = 1.0e-6
let tmult  = 10.0     (* output time factor     *)
let nout   = 12       (* number of output times *)
let nroots = 2        (* number of root functions *)

let print_header rtol avtol yy =
  let open Printf in
  printf "\nidaRoberts_klu: Robertson kinetics DAE serial example problem for IDA.\n";
  printf "         Three equation chemical kinetics problem.\n\n";
  printf "Linear solver: IDAKLU, with user-supplied Jacobian.\n";

  printf "Tolerance parameters:  rtol = %g   atol = %g %g %g \n"
    rtol avtol.{0} avtol.{1} avtol.{2};
  printf "Initial conditions y0 = (%g %g %g)\n"
    yy.{0} yy.{1} yy.{2};

  printf "Constraints and id not used.\n\n";
  printf "-----------------------------------------------------------------------\n";
  printf "  t             y1           y2           y3";
  printf "      | nst  k      h\n";
  printf "-----------------------------------------------------------------------\n";

and print_output ida t y =
  let kused = Ida.get_last_order ida
  and nst = Ida.get_num_steps ida
  and hused = Ida.get_last_step ida
  in
  printf "%10.4e %12.4e %12.4e %12.4e | %3d  %1d %12.4e\n" 
    t y.{0} y.{1} y.{2} nst kused hused

and print_final_stats ida =
  let open Ida in
  let nst  = get_num_steps ida
  and nre  = get_num_res_evals ida
  and nje  = Ida_klu.get_num_jac_evals ida
  and nni  = get_num_nonlin_solv_iters ida
  and netf = get_num_err_test_fails ida
  and ncfn = get_num_nonlin_solv_conv_fails ida
  and nge  = get_num_g_evals ida
  in
  printf "\nFinal Run Statistics: \n\n";
  printf "Number of steps                    = %d\n" nst;
  printf "Number of residual evaluations     = %d\n" nre;
  printf "Number of Jacobian evaluations     = %d\n" nje;
  printf "Number of nonlinear iterations     = %d\n" nni;
  printf "Number of error test failures      = %d\n" netf;
  printf "Number of nonlinear conv. failures = %d\n" ncfn;
  printf "Number of root fn. evaluations     = %d\n" nge;

and print_root_info root_f1 root_f2 =
  (* For printing root_events.  Normally, string_of_root_event makes the output
   * easier to interpret, but we print them as int here in order to get the same
   * output as the C version of this example that comes with sundials.  *)
  let int_of_root_event = function
    | Roots.NoRoot -> 0
    | Roots.Rising -> 1
    | Roots.Falling -> -1
  in
  printf "    rootsfound[] = %3d %3d\n"
    (int_of_root_event root_f1)
    (int_of_root_event root_f2)

let resrob tres (y : RealArray.t) (yp : RealArray.t) (rr : RealArray.t) =
  rr.{0} <- -.0.04*.y.{0} +. 1.0e4*.y.{1}*.y.{2};
  rr.{1} <- -.rr.{0} -. 3.0e7*.y.{1}*.y.{1} -. yp.{1};
  rr.{0} <-  rr.{0} -. yp.{0};
  rr.{2} <-  y.{0} +. y.{1} +. y.{2} -. 1.0

let jacrob { Ida.jac_y = (yval : RealArray.t);
             Ida.jac_coef = cj }
           jacmat =
  let set_col = Sls.SparseMatrix.set_col jacmat in
  let set = Sls.SparseMatrix.set jacmat in

  set_col 0 0;
  set_col 1 3;
  set_col 2 6;
  set_col 3 9;

  set 0 0 (-0.04 -. cj);
  set 1 1 (0.04);
  set 2 2 (1.0);

  set 3 0 (1.0e4*.yval.{2});
  set 4 1 ((-1.0e4*.yval.{2}) -. (6.0e7*.yval.{1}) -. cj);
  set 5 2 (1.0);

  set 6 0 (1.0e4*.yval.{1});
  set 7 1 (-1.0e4*.yval.{1});
  set 8 2 (1.0)

and grob t (y : RealArray.t) y' (gout : RealArray.t) =
  let y1 = y.{0}
  and y3 = y.{2}
  in
  gout.{0} <- y1 -. 0.0001;
  gout.{1} <- y3 -. 0.01

let main () =
  (* Create and initialize y, y', and absolute tolerance vectors.  For
     larger vectors, you might want to use RealArray.create instead of
     RealArray.of_array to avoid making large temporary OCaml
     arrays.  *)
  let y = RealArray.of_array [|1.; 0.; 0.|]
  and y' = RealArray.of_array [|-0.04; 0.04; 0.|]
  and rtol = 1.0e-4
  and avtol = RealArray.of_array [|1.0e-8; 1.0e-6; 1.0e-6|] in
  (* Integration limits *)
  let t0 = 0.0
  and tout1 = 0.4
  in

  (* Wrap y and y' in nvectors.  Operations performed on the wrapped
     representation affect the originals y and y'.  *)
  let wy = Nvector_serial.wrap y
  and wy' = Nvector_serial.wrap y'
  in

  (* Print header information.  *)
  print_header rtol avtol y;

  (* Call IDACreate, IDAInit, and IDARootInit to initialize IDA memory with
   * a 2-component root function and the dense direct linear solver.  *)
  let nnz = neq * neq in
  let ida_mem =
    Ida.init (Ida_klu.klu jacrob nnz)
             (Ida.SVtolerances (rtol, Nvector_serial.wrap avtol))
             resrob ~roots:(nroots, grob) t0 wy wy'
  in
  (* In loop, call IDASolve, print results, and test for error.  Break out of
   * loop when NOUT preset output times have been reached. *)
  let iout = ref 0
  and tout = ref tout1
  in

  let roots = Roots.create nroots in
  let r = Roots.get roots in

  while (!iout <> nout) do
    let (t, flag) = Ida.solve_normal ida_mem !tout wy wy' in
    print_output ida_mem t y;
    match flag with
    | Ida.RootsFound ->
        Ida.get_root_info ida_mem roots;
        print_root_info (r 0) (r 1)

    | Ida.Success ->
        iout := !iout + 1;
        tout := !tout *. tmult

    | Ida.StopTimeReached ->
        iout := nout
  done;

  print_final_stats ida_mem

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
