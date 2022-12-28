(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2007/10/25 20:03:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Sep 2010.
 * -----------------------------------------------------------------
 * Example problem:
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
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a user-supplied function to compute the error weights
 * required for the WRMS norm calculations.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 *)

open Sundials

let unwrap = Nvector.unwrap

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
let t0     = 0.0      (* initial time           *)
let t1     = 0.4      (* first output time      *)
let tmult  = 10.0     (* output time factor     *)
let nout   = 12       (* number of output times *)
let nroots = 2        (* number of root functions *)

let f _ (y : RealArray.t) (yd : RealArray.t) =
  let yd1 = -0.04 *. y.{0} +. 1.0e4 *. y.{1} *. y.{2}
  and yd3 = 3.0e7 *. y.{1} *. y.{1}
  in
  yd.{0} <- yd1;
  yd.{1} <- (-. yd1 -. yd3);
  yd.{2} <- yd3

let g _ (y : RealArray.t) (gout : RealArray.t) =
  gout.{0} <- y.{0} -. 0.0001;
  gout.{1} <- y.{2} -. 0.01

let jac {Cvode.jac_y = (y : RealArray.t)} jmat =
  let set = Matrix.Dense.set jmat in
  set 0 0 (-0.04);
  set 0 1 (1.0e4 *. y.{2});
  set 0 2 (1.0e4 *. y.{1});
  set 1 0 (0.04);
  set 1 1 (-1.0e4 *. y.{2} -. 6.0e7 *. y.{1});
  set 1 2 (-1.0e4 *. y.{1});
  set 2 1 (6.0e7 *. y.{1})

let atol = [| atol1; atol2; atol3 |]
let ewt (y : RealArray.t) (w :RealArray.t) =
  for i = 0 to 2 do
    let yy = y.{i} in
    let ww = rtol *. abs_float(yy) +. atol.(i) in
    if (ww <= 0.0) then raise NonPositiveEwt;
    w.{i} <- (1.0 /. ww)
  done

let print_output =
  printf "At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n"

let print_root_info r1 r2 =
  printf "    rootsfound[] = %3d %3d\n"
    (Roots.int_of_root r1)
    (Roots.int_of_root r2)

let print_final_stats s =
  let open Cvode in
  let nst     = get_num_steps s
  and nfe     = get_num_rhs_evals s
  and nsetups = get_num_lin_solv_setups s
  and netf    = get_num_err_test_fails s
  and nni     = get_num_nonlin_solv_iters s
  and nnf     = get_num_nonlin_solv_conv_fails s
  and nje     = Dls.get_num_jac_evals s
  and nfeLS   = Dls.get_num_lin_rhs_evals s
  and nge     = get_num_g_evals s
  in
  printf "\nFinal Statistics:\n";
  printf "nst = %-6d nfe  = %-6d nsetups = %-6d nfeLS = %-6d nje = %d\n"
    nst nfe nsetups nfeLS nje;
  if Sundials_impl.Version.lt620
  then printf "nni = %-6d ncfn = %-6d netf = %-6d nge = %d\n \n"
              nni nnf netf nge
  else let ncfn = get_num_step_solve_fails s in
       printf "nni = %-6d nnf = %-6d netf = %-6d    ncfn = %-6d  nge = %d\n\n"
              nni nnf netf ncfn nge

let main () =
  (* Create serial vector of length NEQ for I.C. *)
  let y = Nvector_serial.make neq 0.0
  and roots = Roots.create nroots
  in
  let ydata = unwrap y in
  let r = Roots.get roots in

  (* Initialize y *)
  ydata.{0} <- y1;
  ydata.{1} <- y2;
  ydata.{2} <- y3;

  printf " \n3-species kinetics problem\n\n";

  (* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration *)
  (* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. *)
  (* Call CVodeRootInit to specify the root function g with 2 components *)
  (* Call CVDense to specify the CVDENSE dense linear solver *)
  (* Set the Jacobian routine to Jac (user-supplied) *)
  let m = Matrix.dense neq in
  let cvode_mem =
    Cvode.(init BDF ~lsolver:Dls.(solver ~jac (dense y m))
                (WFtolerances ewt) f ~roots:(nroots, g) t0 y)
  in

  (* In loop, call CVode, print results, and test for error.
  Break out of loop when NOUT preset output times have been reached.  *)

  let tout = ref t1
  and iout = ref 0
  in
  while (!iout <> nout) do

    let (t, flag) = Cvode.solve_normal cvode_mem !tout y
    in
    print_output t (ydata.{0}) (ydata.{1}) (ydata.{2});

    match flag with
    | Cvode.RootsFound ->
        Cvode.get_root_info cvode_mem roots;
        print_root_info (r 0) (r 1)

    | Cvode.Success ->
        iout := !iout + 1;
        tout := !tout *. tmult

    | Cvode.StopTimeReached ->
        iout := nout
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
  for _ = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
