(*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:30 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, George D. Byrne,
 *                and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2014.
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
 * Note: This version uses MPI for user routines, and the CVODES
 *       solver. In what follows, N is the number of processors,
 *       N = NPEX*NPEY (see constants below) and it is assumed that
 *       the MPI script mpirun is used to run a parallel
 *       application.
 * If no sensitivities are desired:
 *    % mpirun -np N cvsAdvDiff_FSA_non_p -nosensi
 * If sensitivities are to be computed:
 *    % mpirun -np N cvsAdvDiff_FSA_non_p -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 *)

module Sens = Cvodes.Sensitivity
module RealArray = Sundials.RealArray
open Bigarray

let local_array = Nvector_parallel.local_array
let printf = Printf.printf
let eprintf = Printf.eprintf

let vmax_norm = Nvector_parallel.Ops.n_vmaxnorm

(* Problem Constants *)
let xmax =  2.0   (* domain boundary           *)
let mx =    10    (* mesh dimension            *)
let neq =   mx    (* number of equations       *)
let atol =  1.e-5 (* scalar absolute tolerance *)
let t0 =    0.0   (* initial time              *)
let t1 =    0.5   (* first output time         *)
let dtout = 0.5   (* output time increment     *)
let nout =  10    (* number of output times    *)

let np =    2
let ns =    2

let zero =  0.0

(* Type : UserData
   contains problem parameters, grid constants, work array. *)

type user_data = {
  p       : RealArray.t;      (* model parameters             *)

  dx      : float;            (* spatial discretization grid  *)

  npes    : int;              (* total number of processes    *)
  my_pe   : int;              (* current ID                   *)

  comm    : Mpi.communicator; (* MPI communicator             *)

  z       : float array;      (* work space                   *)
}

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(* Process and verify arguments to cvsfwdnonx_p. *)

let wrong_args my_pe name =
  if my_pe  = 0 then begin
    printf "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n" name;
    printf "         sensi_meth = sim, stg, or stg1\n";
    printf "         err_con    = t or f\n"
  end;
  exit 0

let process_args my_pe =
  let argv = Sys.argv in
  let argc = Array.length argv in
  if argc < 2 then wrong_args my_pe argv.(0);

  let sensi =
    if argv.(1) = "-nosensi" then false
    else if argv.(1) = "-sensi" then true
    else wrong_args my_pe argv.(0)
  in

  if not sensi then (None, false)
  else begin
    if argc <> 4 then wrong_args my_pe argv.(0);
    let sensi_meth =
      if argv.(2) = "sim" then Sens.Simultaneous
      else if argv.(2) = "stg" then Sens.Staggered
      else if argv.(2) = "stg1" then Sens.Staggered1
      else wrong_args my_pe argv.(0)
    in
    let err_con =
      if argv.(3) = "t" then true
      else if argv.(3) = "f" then false
      else wrong_args my_pe argv.(0)
    in
    (Some sensi_meth, err_con)
  end

(* Set initial conditions in u vector *)

let set_ic u dx my_length my_base =
  let udata = local_array u in

  (* Load initial profile into u vector *)
  for i=1 to my_length do
    let iglobal = my_base + i in
    let x = float iglobal *. dx in
    udata.{i - 1} <- x*.(xmax -. x) *. exp(2.0*.x)
  done

(* Print current t, step count, order, stepsize, and max norm of solution  *)

let print_output cvode_mem my_pe t u =
  let nst = Cvode.get_num_steps cvode_mem in
  let qu  = Cvode.get_last_order cvode_mem in
  let hu  = Cvode.get_last_step cvode_mem in
  let umax = vmax_norm u in

  if my_pe = 0 then begin
    printf "%8.3e %2d  %8.3e %5d\n" t qu hu nst;
    printf "                                Solution       ";
    printf "%12.4e \n" umax
  end

(* Print max norm of sensitivities *)

let print_output_s my_pe uS =
  let smax = vmax_norm uS.(0) in
  if my_pe = 0 then begin
    printf "                                Sensitivity 1  ";
    printf "%12.4e \n" smax
  end;
  let smax = vmax_norm uS.(1) in
  if my_pe = 0 then begin
    printf "                                Sensitivity 2  ";
    printf "%12.4e \n" smax
  end

(* Print some final statistics located in the iopt array *)

let print_final_stats cvode_mem sensi =
  let open Cvode in
  let nst     = get_num_steps cvode_mem
  and nfe     = get_num_rhs_evals cvode_mem
  and nsetups = get_num_lin_solv_setups cvode_mem
  and netf    = get_num_err_test_fails cvode_mem
  and nni     = get_num_nonlin_solv_iters cvode_mem
  and ncfn    = get_num_nonlin_solv_conv_fails cvode_mem in
  printf "\nFinal Statistics\n\n";
  printf "nst     = %5d\n\n" nst;
  printf "nfe     = %5d\n"   nfe;
  printf "netf    = %5d    nsetups  = %5d\n" netf nsetups;
  printf "nni     = %5d    ncfn     = %5d\n" nni ncfn;

  if sensi then begin
    let open Sens in
    let nfSe     = get_num_rhs_evals cvode_mem
    and nfeS     = get_num_rhs_evals_sens cvode_mem
    and nsetupsS = get_num_lin_solv_setups cvode_mem
    and netfS    = get_num_err_test_fails cvode_mem
    and nniS     = get_num_nonlin_solv_iters cvode_mem
    and ncfnS    = get_num_nonlin_solv_conv_fails cvode_mem in
    printf "\n";
    printf "nfSe    = %5d    nfeS     = %5d\n" nfSe nfeS;
    printf "netfs   = %5d    nsetupsS = %5d\n" netfS nsetupsS;
    printf "nniS    = %5d    ncfnS    = %5d\n" nniS ncfnS
  end

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 *)

(* f routine. Compute f(t,u). *)

let f data t (udata, _, _) (dudata, _, _) =

  (* Extract needed problem constants from data *)
  let dx = data.dx in
  let hordc = data.p.{0}/.(dx*.dx) in
  let horac = data.p.{1}/.(2.0*.dx) in

  (* Extract parameters for parallel computation. *)
  let comm = data.comm in
  let npes = data.npes in             (* Number of processes. *)
  let my_pe = data.my_pe in           (* Current process number. *)
  let my_length = Array1.dim udata in (* Number of local elements of u. *)
  let z = data.z in

  (* Compute related parameters. *)
  let my_pe_m1 = my_pe - 1 in
  let my_pe_p1 = my_pe + 1 in
  let last_pe  = npes - 1 in

  (* Store local segment of u in the working array z. *)
  for i = 1 to my_length do
    z.(i) <- udata.{i - 1};
  done;

  (* Pass needed data to processes before and after current process. *)
  if my_pe <> 0 then Mpi.send_float z.(1) my_pe_m1 0 comm;
  if my_pe <> last_pe then Mpi.send_float z.(my_length) my_pe_p1 0 comm;

  (* Receive needed data from processes before and after current process. *)
  z.(0) <- if my_pe <> 0 then Mpi.receive_float my_pe_m1 0 comm else zero;
  z.(my_length + 1) <-
    if my_pe <> last_pe then Mpi.receive_float my_pe_p1 0 comm else zero;

  (* Loop over all grid points in current process. *)
  for i=1 to my_length do
    (* Extract u at x_i and two neighboring points *)
    let ui = z.(i) in
    let ult = z.(i-1) in
    let urt = z.(i+1) in

    (* Set diffusion and advection terms and load into udot *)
    let hdiff = hordc*.(ult -. 2.0*.ui +. urt) in
    let hadv = horac*.(urt -. ult) in
    dudata.{i-1} <- hdiff +. hadv
  done

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* Get processor number, total number of pe's, and my_pe. *)
  let comm  = Mpi.comm_world in
  let npes  = Mpi.comm_size comm in
  let my_pe = Mpi.comm_rank comm in

  (* Process arguments *)
  let sensi, err_con = process_args my_pe in

  (* Set local vector length. *)
  let nperpe = neq/npes in
  let nrem = neq - npes*nperpe in
  let local_n = if my_pe < nrem then nperpe+1 else nperpe in
  let my_base = if my_pe < nrem then my_pe*local_n else my_pe*nperpe + nrem in

  (* USER DATA STRUCTURE *)
  let dx = xmax/.(float (mx+1)) in
  let p = RealArray.of_list [ 1.0; 0.5 ] in
  let data = {
      p       = p;
      dx      = dx;
      npes    = npes;
      my_pe   = my_pe;
      comm    = comm;
      z       = Array.make 100 0.0;
    } in

  (* INITIAL STATES *)
  let u = Nvector_parallel.make local_n neq comm 0.0 in
  set_ic u dx local_n my_base;

  (* TOLERANCES *)
  let reltol = zero in
  let abstol = atol in

  (* CVODE_CREATE & CVODE_MALLOC *)
  let cvode_mem = Cvode.(init Adams
                              Functional
                              (SStolerances (reltol, abstol))
                              (f data)
                              t0
                              u)
  in

  if my_pe = 0 then begin
    printf "\n1-D advection-diffusion equation, mesh size =%3d \n" mx;
    printf "\nNumber of PEs = %3d \n" npes;
  end;

  let print_sensi =
    match sensi with
    | None -> begin
          if my_pe = 0 then printf "Sensitivity: NO ";
          (fun _ -> ())
        end
    | Some sensi_meth -> begin
        (* sensitivity w.r.t. i-th parameter *)
        let plist = Array.init ns (fun i -> i) in
        let pbar = RealArray.create ns in
        RealArray.mapi (fun is _ -> data.p.{plist.(is)}) pbar;

        let uS = Array.init ns
                   (fun _ -> Nvector_parallel.make local_n neq comm 0.0)
        in

        Sens.(init cvode_mem
                         EEtolerances
                         sensi_meth
                         ~sens_params:{ pvals = Some data.p;
                                        pbar  = Some pbar;
                                        plist = Some plist; }
                         (OneByOne None)
                         uS);
        Sens.set_err_con cvode_mem err_con;
        Sens.(set_dq_method cvode_mem DQCentered 0.0);

        if my_pe = 0 then begin
          printf "Sensitivity: YES ";
          (match sensi_meth with
           | Sens.Simultaneous -> printf "( SIMULTANEOUS +"
           | Sens.Staggered    -> printf "( STAGGERED +"
           | Sens.Staggered1   -> printf "( STAGGERED1 +");
          if err_con then printf " FULL ERROR CONTROL )"
                     else printf " PARTIAL ERROR CONTROL )";
        end;
        (fun s -> (ignore (Sens.get s uS); print_output_s my_pe uS))
      end
  in

  (* In loop over output points, call CVode, print results, test for error *)

  if my_pe = 0 then begin
    printf "\n\n";
    printf "============================================================\n";
    printf "     T     Q       H      NST                    Max norm   \n";
    printf "============================================================\n";
  end;

  let tout = ref t1 in
  for iout = 1 to nout do
    let t, _ = Cvode.solve_normal cvode_mem !tout u in
    print_output cvode_mem my_pe t u;
    print_sensi cvode_mem;
    if my_pe = 0 then
      printf "------------------------------------------------------------\n";
    tout := !tout +. dtout
  done;

  (* Print final statistics *)
  if my_pe = 0 then print_final_stats cvode_mem (sensi <> None)

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
