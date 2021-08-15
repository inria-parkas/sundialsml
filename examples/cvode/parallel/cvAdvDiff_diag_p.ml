(*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:28 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, George Byrne,
 *                and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jul 2014.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODE. The problem is the semi-discrete
 * form of the advection-diffusion equation in 1-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is the following:
 *   u(x,t=0) = x(2-x)exp(2x) .
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the ADAMS integration method,
 * and with Newton iteration using diagonal approximate Jacobians.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * This version uses MPI for user routines.
 * Execute with Number of Processors = N,  with 1 <= N <= MX.
 * -----------------------------------------------------------------
 *)

open Sundials

open Bigarray
let local_array = Nvector_parallel.local_array

let printf = Printf.printf

let vmaxnorm = Nvector_parallel.Ops.maxnorm

(* Problem Constants *)

let zero =  0.0

let xmax =  2.0    (* domain boundary           *)
let mx =    10     (* mesh dimension            *)
let neq =   mx     (* number of equations       *)
let atol =  1.0e-5 (* scalar absolute tolerance *)
let t0 =    zero   (* initial time              *)
let t1 =    0.5    (* first output time         *)
let dtout = 0.5    (* output time increment     *)
let nout =  10     (* number of output times    *)

(* Type : UserData
   contains grid constants, parallel machine parameters, work array. *)

type user_data = {
  dx : float;
  hdcoef : float;
  hacoef : float;

  npes   : int;
  my_pe  : int;

  comm   : Mpi.communicator;

  z      : float array;
}

(************************ Private Helper Functions ***********************)

(* Set initial conditions in u vector *)

let set_ic u dx my_length my_base =
  (* Set pointer to data array and get local length of u. *)
  let udata = local_array u in
  let my_length = Array1.dim udata in

  (* Load initial profile into u vector *)
  for i=1 to my_length do
    let iglobal = float (my_base + i) in
    let x = iglobal*.dx in
    udata.{i-1} <- x*.(xmax -. x)*.exp(2.0*.x)
  done

(* Print problem introduction *)

let print_intro npes =
  printf "\n 1-D advection-diffusion equation, mesh size =%3d \n" mx;
  printf "\n Number of PEs = %3d \n"  npes;
  printf "\n Diagonal linear solver CVDiag \n\n"

(* Print data *)

let print_data t umax nst =
  printf "At t = %4.2f  max.norm(u) =%14.6e  nst =%4d \n" t umax nst

(* Print some final statistics located in the iopt array *)

let print_final_stats s =
  let open Cvode in
  let nst  = get_num_steps s
  and nfe  = get_num_rhs_evals s
  and netf = get_num_err_test_fails s
  and nni  = get_num_nonlin_solv_iters s
  and ncfn = get_num_nonlin_solv_conv_fails s
  in
  printf "\nFinal Statistics: \n\n";
  printf "nst = %-6d  nfe  = %-6d  " nst nfe;
  printf "nni = %-6d  ncfn = %-6d  netf = %d\n \n" nni ncfn netf

(***************** Function Called by the Solver ***********************)

(* f routine. Compute f(t,u). *)

let f data t (udata, _, _) (dudata, _, _) =
  (* Extract needed problem constants from data *)
  let hordc = data.hdcoef
  and horac = data.hacoef
  in

  (* Extract parameters for parallel computation. *)
  let comm      = data.comm in
  let npes      = data.npes in        (* Number of processes. *)
  let my_pe     = data.my_pe in       (* Current process number. *)
  let my_length = Array1.dim udata in (* Number of local elements of u. *)
  let z         = data.z in

  (* Compute related parameters. *)
  let my_pe_m1 = my_pe - 1 in
  let my_pe_p1 = my_pe + 1 in
  let last_pe = npes - 1 in

  (* Store local segment of u in the working array z. *)
  for i = 1 to my_length do
    z.(i) <- udata.{i - 1}
  done;

  (* Pass needed data to processes before and after current process. *)
  if my_pe <> 0 then Mpi.send_float z.(1) my_pe_m1 0 comm;
  if my_pe <> last_pe then Mpi.send_float z.(my_length) my_pe_p1 0 comm;

  (* Receive needed data from processes before and after current process. *)
  z.(0) <-
    if my_pe <> 0 then Mpi.receive_float my_pe_m1 0 comm else zero;
  z.(my_length + 1) <-
    if my_pe <> last_pe then Mpi.receive_float my_pe_p1 0 comm else zero;

  (* Loop over all grid points in current process. *)
  for i=1 to my_length do
    (* Extract u at x_i and two neighboring points *)
    let ui  = z.(i) in
    let ult = z.(i-1) in
    let urt = z.(i+1) in

    (* Set diffusion and advection terms and load into udot *)
    let hdiff = hordc*.(ult -. 2.0*.ui +. urt) in
    let hadv  = horac*.(urt -. ult) in
    dudata.{i-1} <- hdiff +. hadv
  done

(***************************** Main Program ******************************)

let main () =
  (* Get processor number, total number of pe's, and my_pe. *)
  let comm  = Mpi.comm_world in
  let npes  = Mpi.comm_size comm in
  let my_pe = Mpi.comm_rank comm in

  (* Set local vector length. *)
  let nperpe = neq/npes in
  let nrem = neq - npes*nperpe in
  let local_N = if my_pe < nrem then nperpe+1 else nperpe in
  let my_base = if my_pe < nrem then my_pe*local_N else my_pe*nperpe + nrem in

  let u = Nvector_parallel.make local_N neq comm 0.0 in (* Allocate u vector *)

  let reltol = zero  (* Set the tolerances *)
  and abstol = atol
  in

  let dx = xmax/.(float (mx+1)) in  (* Set grid coefficients in data *)

  let data = {
      dx     = dx;
      hdcoef = 1.0/.(dx*.dx);
      hacoef = 0.5/.(2.0*.dx);

      npes   = npes;
      my_pe  = my_pe;

      comm   = comm;

      z      = Array.make 100 0.0;
    }
  in
  set_ic u dx local_N my_base;  (* Initialize u vector *)

  let cvode_mem = Cvode.(init Adams
                              (SStolerances (reltol, abstol))
                              (f data) t0 u)
  in

  if my_pe = 0 then print_intro npes;

  let umax = vmaxnorm u in
  if my_pe = 0 then print_data t0 umax 0;

  (* In loop over output points, call CVode, print results, test for error *)

  let tout = ref t1 in
  for iout = 1 to nout do
    let (t, flag) = Cvode.solve_normal cvode_mem !tout u
    in
    let nst  = Cvode.get_num_steps cvode_mem in
    let umax = vmaxnorm u in
    if my_pe = 0 then print_data t umax nst;
    tout := !tout +. dtout
  done;

  (* Print some final statistics *)
  if my_pe = 0 then print_final_stats cvode_mem

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
