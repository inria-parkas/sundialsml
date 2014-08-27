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

module Sens = Cvodes.Sensitivity
module RealArray = Sundials.RealArray
let unvec = Sundials.unvec

(* As a slight deviation from the sundials/C code, we allow an extra
   argument to repeat the test, used to check that garbage collection
   works properly.  argv is updated to remove that extra argument so
   the rest of the test can exactly mirror the C code.  *)
let argv = ref Sys.argv

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

let f data t u udot =
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

let process_args () =
  let argv = !argv in
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
      else if argv.(2) = "stg" then Sens.Staggered
      else if argv.(2) = "stg1" then Sens.Staggered1
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
  printf "                                Solution       ";
  printf "%12.4e \n" (vmax_norm u)

(* Print max norm of sensitivities *)

let print_output_s uS =
  printf "                                Sensitivity 1  ";
  printf "%12.4e \n" (vmax_norm uS.(0));
  printf "                                Sensitivity 2  ";
  printf "%12.4e \n" (vmax_norm uS.(1))


(* Print some final statistics located in the CVODES memory *)

let print_final_stats cvode_mem sensi =
  let nst = Cvode.get_num_steps cvode_mem
  and nfe = Cvode.get_num_rhs_evals cvode_mem
  and nsetups = Cvode.get_num_lin_solv_setups cvode_mem
  and netf = Cvode.get_num_err_test_fails cvode_mem
  and nni = Cvode.get_num_nonlin_solv_iters cvode_mem
  and ncfn = Cvode.get_num_nonlin_solv_conv_fails cvode_mem in
  printf "\nFinal Statistics\n\n";
  printf "nst     = %5d\n\n" nst;
  printf "nfe     = %5d\n"   nfe;
  printf "netf    = %5d    nsetups  = %5d\n" netf nsetups;
  printf "nni     = %5d    ncfn     = %5d\n" nni ncfn;

  if sensi then begin
    let nfSe = Sens.get_num_rhs_evals cvode_mem
    and nfeS = Sens.get_num_rhs_evals_sens cvode_mem
    and nsetupsS = Sens.get_num_lin_solv_setups cvode_mem
    and netfS = Sens.get_num_err_test_fails cvode_mem
    and nniS = Sens.get_num_nonlin_solv_iters cvode_mem
    and ncfnS = Sens.get_num_nonlin_solv_conv_fails cvode_mem in
    printf "\n";
    printf "nfSe    = %5d    nfeS     = %5d\n" nfSe nfeS;
    printf "netfs   = %5d    nsetupsS = %5d\n" netfS nsetupsS;
    printf "nniS    = %5d    ncfnS    = %5d\n" nniS ncfnS
  end

let main () =
  (* Process arguments *)
  let sensi, err_con = process_args () in

  (* Set user data *)
  let dx = xmax/.(float (mx+1)) in
  let data = {
      p = RealArray.of_list [ 1.0; 0.5 ];
      dx = dx;
    } in

  (* Allocate and set initial states *)
  let u = Nvector_serial.make neq 0.0 in
  set_ic (unvec u) dx;

  (* Set integration tolerances *)
  let reltol = zero in
  let abstol = atol in

  (* Create CVODES object *)
  let cvode_mem = Cvode.init
                    Cvode.Adams
                    Cvode.Functional
                    (Cvode.SStolerances (reltol, abstol))
                    (f data)
                    ~t0:t0
                    u
  in
  printf "\n1-D advection-diffusion equation, mesh size =%3d\n" mx;

  (* Sensitivity-related settings *)
  let print_sensi =
    match sensi with
    | None -> (printf "Sensitivity: NO "; (fun _ -> ()))
    | Some sensi_meth -> begin
        let plist = Array.init ns (fun i -> i) in
        let pbar = RealArray.create ns in
        RealArray.mapi (fun is _ -> data.p.{plist.(is)}) pbar;

        let uS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0) in

        Sens.init cvode_mem
                         Sens.EEtolerances
                         sensi_meth
                         { Sens.pvals = Some data.p;
                           Sens.pbar = Some pbar;
                           Sens.plist = Some plist; }
                         (Sens.OneByOne None)
                         uS;
        Sens.set_err_con cvode_mem err_con;
        Sens.set_dq_method cvode_mem Sens.DQCentered 0.0;

        printf "Sensitivity: YES ";
        (match sensi_meth with
         | Sens.Simultaneous -> printf "( SIMULTANEOUS +"
         | Sens.Staggered    -> printf "( STAGGERED +"
         | Sens.Staggered1   -> printf "( STAGGERED1 +");
        printf (if err_con then " FULL ERROR CONTROL )"
                           else " PARTIAL ERROR CONTROL )");

        (fun s -> (ignore (Sens.get s uS); print_output_s uS))
      end
  in

  (* In loop over output points, call CVode, print results, test for error *)
  printf "\n\n";
  printf "============================================================\n";
  printf "     T     Q       H      NST                    Max norm   \n";
  printf "============================================================\n";

  let tout = ref t1 in
  for iout = 1 to nout do
    let t, _ = Cvode.solve_normal cvode_mem !tout u in
    print_output cvode_mem t u;
    print_sensi cvode_mem;
    printf "------------------------------------------------------------\n";
    tout := !tout +. dtout
  done;

  (* Print final statistics *)
  print_final_stats cvode_mem (sensi <> None)

(* Check if the last argument is a repetition count.  *)
let reps =
  let n = Array.length !argv in
  try
    if n >= 2 then
      let reps = int_of_string !argv.(n-1) in
      argv := Array.sub !argv 0 (n-1);
      reps
    else 1
  with _ -> 1
let _ = for i = 1 to reps do main () done
let _ = Gc.compact ()
