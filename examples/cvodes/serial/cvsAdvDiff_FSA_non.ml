(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/31 00:04:42 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, George D. Byrne,
 *              and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2014.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODES. The problem is the semi-discrete form of
 * the advection-diffusion equation in 1-D:
 *   du/dt = q1 * d^2 u / dx^2 + q2 * du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is:
 *   u(x,y,t=0) = x(2-x)exp(2x).
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and functional iteration.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters q1 and q2.
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on FULL or PARTIAL,
 * respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % cvsAdvDiff_FSA_non -nosensi
 * If sensitivities are to be computed:
 *    % cvsAdvDiff_FSA_non -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 *)

open Sundials

module Sens = Cvodes.Sensitivity
let unwrap = Nvector.unwrap

let printf = Printf.printf
let vmax_norm = Nvector_serial.Ops.n_vmaxnorm

(* Problem Constants *)

let xmax  = 2.0   (* domain boundary           *)
let mx    = 10    (* mesh dimension            *)
let neq   = mx    (* number of equations       *)
let atol  = 1.e-5 (* scalar absolute tolerance *)
let t0    = 0.0   (* initial time              *)
let t1    = 0.5   (* first output time         *)
let dtout = 0.5   (* output time increment     *)
let nout  = 10    (* number of output times    *)

let np    = 2
let ns    = 2

let zero  = 0.0

(* Type : UserData
   contains problem parameters, grid constants, work array. *)

type user_data = {
  p  : RealArray.t;
  dx : float;
}

(* f routine. Compute f(t,u). *)

let f data t (u : RealArray.t) (udot : RealArray.t) =
  (* Extract needed problem constants from data *)
  let dx = data.dx in
  let hordc = data.p.{0}/.(dx*.dx) in
  let horac = data.p.{1}/.(2.0*.dx) in

  (* Loop over all grid points. *)
  for i=0 to (neq - 1) do

    (* Extract u at x_i and two neighboring points *)
    let ui = u.{i} in
    let ult = if i <> 0 then u.{i-1} else zero in
    let urt = if i <> neq-1 then u.{i+1} else zero in

    (* Set diffusion and advection terms and load into udot *)
    let hdiff = hordc*.(ult -. 2.0*.ui +. urt) in
    let hadv = horac*.(urt -. ult) in
    udot.{i} <- hdiff +. hadv
  done

(* Process and verify arguments. *)

let wrong_args name =
  printf "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n" name;
  printf "         sensi_meth = sim stg or stg1\n";
  printf "         err_con    = t or f\n";
  exit 0

let process_args u =
  let argv = Sys.argv in
  let argc = Array.length argv in
  if argc < 2 then wrong_args argv.(0);

  let sensi =
    if argv.(1) = "-nosensi" then false
    else if argv.(1) = "-sensi" then true
    else wrong_args argv.(0)
  in

  if not sensi then (None, false)
  else begin
    if argc <> 4 then wrong_args argv.(0);
    let sensi_meth =
      if argv.(2) = "sim" then Sens.Simultaneous
                          (Some (NonlinearSolver.FixedPoint.make_sens (ns+1) u))
      else if argv.(2) = "stg" then Sens.Staggered
                          (Some (NonlinearSolver.FixedPoint.make_sens ns u))
      else if argv.(2) = "stg1" then Sens.Staggered1
                          (Some (NonlinearSolver.FixedPoint.make u))
      else wrong_args argv.(0)
    in
    let err_con =
      if argv.(3) = "t" then true
      else if argv.(3) = "f" then false
      else wrong_args argv.(0)
    in
    (Some sensi_meth, err_con)
  end

(* Set initial conditions in u vector. *)

let set_ic u dx =
  (* Load initial profile into u vector *)
  for i=0 to neq - 1 do
    let x = float(i+1)*.dx in
    u.{i} <- x*.(xmax -. x) *. exp(2.0*.x);
  done

(* Print current t, step count, order, stepsize, and max norm of solution *)

let print_output cvode_mem t u =
  let nst = Cvode.get_num_steps cvode_mem in
  let qu  = Cvode.get_last_order cvode_mem in
  let hu  = Cvode.get_last_step cvode_mem in
  printf "%8.3e %2d  %8.3e %5d\n" t qu hu nst;
  print_string "                                Solution       ";
  printf "%12.4e \n" (vmax_norm u)

(* Print max norm of sensitivities *)

let print_output_s uS =
  print_string "                                Sensitivity 1  ";
  printf "%12.4e \n" (vmax_norm uS.(0));
  print_string "                                Sensitivity 2  ";
  printf "%12.4e \n" (vmax_norm uS.(1))

(* Print some final statistics located in the CVODES memory *)

let print_string_5d s i =
  print_string s;
  if i < 10 then print_string "    "
  else if i < 100 then print_string "   "
  else if i < 1000 then print_string "  "
  else if i < 10000 then print_string " ";
  print_int i

let print_string_3d s i =
  print_string s;
  if i < 10 then print_string "  "
  else if i < 100 then print_string " ";
  print_int i

let print_final_stats cvode_mem sensi =
  let open Cvode in
  let nst     = get_num_steps cvode_mem
  and nfe     = get_num_rhs_evals cvode_mem
  and nsetups = get_num_lin_solv_setups cvode_mem
  and netf    = get_num_err_test_fails cvode_mem
  and nni     = get_num_nonlin_solv_iters cvode_mem
  and ncfn    = get_num_nonlin_solv_conv_fails cvode_mem in
  print_string "\nFinal Statistics\n\n";
  print_string_5d "nst     = " nst;
  print_string_5d "\n\nnfe     = " nfe;
  print_string_5d "\nnetf    = " netf;
  print_string_5d "    nsetups  = " nsetups;
  print_string_5d "\nnni     = " nni;
  print_string_5d "    ncfn     = " ncfn;
  print_newline ();

  match sensi with
  | None -> ()
  | Some sensi_meth -> begin
      let open Sens in
      let nfSe     = get_num_rhs_evals cvode_mem
      and nfeS     = get_num_rhs_evals_sens cvode_mem
      and nsetupsS = get_num_lin_solv_setups cvode_mem
      and netfS    = get_num_err_test_fails cvode_mem
      in
      let nniS, ncfnS =
        match sensi_meth with
        | Staggered _ | Staggered1 _ ->
            get_num_nonlin_solv_iters cvode_mem,
            get_num_nonlin_solv_conv_fails cvode_mem
        | Simultaneous _ -> 0, 0
      in
      print_newline ();
      print_string_5d "nfSe    = " nfSe;
      print_string_5d "    nfeS     = " nfeS;
      print_string_5d "\nnetfs   = " netfS;
      print_string_5d "    nsetupsS = " nsetupsS;
      print_string_5d "\nnniS    = " nniS;
      print_string_5d "    ncfnS    = " ncfnS;
      print_newline ()
  end

let main () =
  let u = Nvector_serial.make neq 0.0 in

  (* Process arguments *)
  let sensi, err_con = process_args u in

  (* Set user data *)
  let dx = xmax/.(float (mx+1)) in
  let data = {
      p = RealArray.of_list [ 1.0; 0.5 ];
      dx = dx;
    } in

  (* Allocate and set initial states *)
  set_ic (unwrap u) dx;

  (* Set integration tolerances *)
  let reltol = zero in
  let abstol = atol in

  (* Create CVODES object *)
  let nlsolver = Sundials.NonlinearSolver.FixedPoint.make u in
  let cvode_mem = Cvode.(init
                           Adams
                           (SStolerances (reltol, abstol))
                           ~nlsolver
                           (f data)
                           t0
                           u)
  in
  print_string_3d "\n1-D advection-diffusion equation, mesh size =" mx;
  print_newline ();

  (* Sensitivity-related settings *)
  let print_sensi =
    match sensi with
    | None -> (print_string "Sensitivity: NO "; (fun _ -> ()))
    | Some sensi_meth -> begin
        let plist = Array.init ns (fun i -> i) in
        let pbar = RealArray.create ns in
        RealArray.mapi (fun is _ -> data.p.{plist.(is)}) pbar;

        let uS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0) in

        Sens.(init cvode_mem
                   EEtolerances
                   sensi_meth
                   ~sens_params:{ pvals = Some data.p;
                                  pbar = Some pbar;
                                  plist = Some plist; }
                   (OneByOne None)
                   uS);
        Sens.set_err_con cvode_mem err_con;
        Sens.set_dq_method cvode_mem Sens.DQCentered 0.0;

        print_string "Sensitivity: YES ";
        (match sensi_meth with
         | Sens.Simultaneous _ -> print_string "( SIMULTANEOUS +"
         | Sens.Staggered _    -> print_string "( STAGGERED +"
         | Sens.Staggered1 _   -> print_string "( STAGGERED1 +");
        print_string (if err_con then " FULL ERROR CONTROL )"
                                 else " PARTIAL ERROR CONTROL )");

        (fun s -> (ignore (Sens.get s uS); print_output_s uS))
      end
  in

  (* In loop over output points, call CVode, print results, test for error *)
  print_string "\n\n";
  print_string "============================================================\n";
  print_string "     T     Q       H      NST                    Max norm   \n";
  print_string "============================================================\n";

  let tout = ref t1 in
  for iout = 1 to nout do
    let t, _ = Cvode.solve_normal cvode_mem !tout u in
    print_output cvode_mem t u;
    print_sensi cvode_mem;
    print_string "------------------------------------------------------------\n";
    tout := !tout +. dtout
  done;

  (* Print final statistics *)
  print_final_stats cvode_mem sensi

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
