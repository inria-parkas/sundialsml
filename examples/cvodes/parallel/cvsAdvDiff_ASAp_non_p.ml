(*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:30 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2014.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODE. The problem is the semi-discrete form of
 * the advection-diffusion equation in 1-D:
 *   du/dt = p1 * d^2u / dx^2 + p2 * du / dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is:
 *   u(x,t=0) = x(2-x)exp(2x).
 * The nominal values of the two parameters are: p1=1.0, p2=0.5
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and functional iteration.
 * It uses scalar relative and absolute tolerances.
 *
 * In addition to the solution, sensitivities with respect to p1
 * and p2 as well as with respect to initial conditions are
 * computed for the quantity:
 *    g(t, u, p) = int_x u(x,t) at t = 5
 * These sensitivities are obtained by solving the adjoint system:
 *    dv/dt = -p1 * d^2 v / dx^2 + p2 * dv / dx
 * with homogeneous Ditrichlet boundary conditions and the final
 * condition:
 *    v(x,t=5) = 1.0
 * Then, v(x, t=0) represents the sensitivity of g(5) with respect
 * to u(x, t=0) and the gradient of g(5) with respect to p1, p2 is
 *    (dg/dp)^T = [  int_t int_x (v * d^2u / dx^2) dx dt ]
 *                [  int_t int_x (v * du / dx) dx dt     ]
 *
 * This version uses MPI for user routines.
 * Execute with Number of Processors = N,  with 1 <= N <= MX.
 * -----------------------------------------------------------------
 *)

module Adjoint = Cvodes.Adjoint
module RealArray = Sundials.RealArray
open Bigarray

let unvec = Sundials.unvec
let printf = Printf.printf
let eprintf = Printf.eprintf

let blit buf buf_offset dst dst_offset len =
  for i = 0 to len-1 do
    dst.(dst_offset + i) <- buf.{buf_offset + i}
  done

(* Problem Constants *)

let xmax  = 2.0   (* domain boundary            *)
let mx    = 20    (* mesh dimension             *)
let neq   = mx    (* number of equations        *)
let atol  = 1.e-5 (* scalar absolute tolerance  *)
let t0    = 0.0   (* initial time               *)
let tout  = 2.5   (* output time increment      *)

(* Adjoint Problem Constants *)

let np    = 2     (* number of parameters       *)
let steps = 200   (* steps between check points *)

let zero  = 0.0
let one   = 1.0
let two   = 2.0

(* Type : UserData *)

type user_data = {
  p       : float array;      (* model parameters             *)

  dx      : float;            (* spatial discretization grid  *)
  hdcoef  : float;            (* diffusion coefficient        *)
  hacoef  : float;            (* advection coefficient        *)

  local_n : int;

  npes    : int;              (* total number of processes    *)
  my_pe   : int;              (* current ID                   *)

  nperpe  : int;
  nrem    : int;

  comm    : Mpi.communicator; (* MPI communicator             *)

  z1      : float array;      (* work space                   *)
  z2      : float array;
}

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(* Set initial conditions in u vector *)

let set_ic u dx my_length my_base =
  (* Set pointer to data array and get local length of u. *)
  let udata, _, _ = unvec u in
  let my_length = Array1.dim udata in

  (* Load initial profile into u vector *)
  for i=1 to my_length do
    let iglobal = float (my_base + i) in
    let x = iglobal*.dx in
    udata.{i-1} <- x*.(xmax -. x)*.exp(two*.x)
  done

(* Set final conditions in uB vector *)

let set_ic_back uB my_base =
  (* Set pointer to data array and get local length of uB *)
  let uBdata, _, _ = unvec uB in

  (* Set adjoint states to 1.0 and quadrature variables to 0.0 *)
  Array1.fill uBdata (if my_base = -1 then zero else one)

(* Compute local value of the space integral int_x z(x) dx *)

let xintgr z l dx =
  let my_intgr = ref (0.5*.(z.{0} +. z.{l-1})) in
  for i = 1 to l-2 do
    my_intgr := !my_intgr +. z.{i}
  done;
  !my_intgr *. dx

let xintgr' z l dx =
  let my_intgr = ref (0.5*.(z.(0) +. z.(l-1))) in
  for i = 1 to l-2 do
    my_intgr := !my_intgr +. z.(i)
  done;
  !my_intgr *. dx

(* Compute value of g(u) *)

let compute_g data u =
  (* Extract MPI info. from data *)
  let comm  = data.comm in
  let npes  = data.npes in
  let my_pe = data.my_pe in
  let dx    = data.dx in

  if my_pe = npes then begin (* Loop over all other processes and sum *)
    let intgr = ref zero in
    for i=0 to npes-1 do
      intgr := !intgr +. Mpi.receive_float i 0 comm
    done;
    !intgr
  end else begin             (* Compute local portion of the integral *)
    let udata, _, _ = unvec u in
    let my_length = Array1.dim udata in
    let my_intgr = xintgr udata my_length dx in
    Mpi.send_float my_intgr npes 0 comm;
    my_intgr
  end

(* Print output after backward integration *)

let print_output data g_val uB =
  let comm    = data.comm in
  let npes    = data.npes in
  let my_pe   = data.my_pe in
  let nperpe  = data.nperpe in
  let nrem    = data.nrem in

  let uBdata, _, _ = unvec uB in

  if my_pe = npes then begin

    printf "\ng(tf) = %8e\n\n" g_val;
    printf "dgdp(tf)\n  [ 1]: %8e\n  [ 2]: %8e\n\n"
              (-.uBdata.{0}) (-.uBdata.{1});

    let mu = Array.make neq 0.0 in

    let indx = ref 0 in
    for i = 0 to npes - 1 do
      let ni = if i < nrem then nperpe+1 else nperpe in
      let buf = (Mpi.receive i 0 comm : RealArray.t) in
      blit buf 0 mu !indx ni;
      indx := !indx + ni
    done;

    printf "mu(t0)\n";
    Array.iteri (fun i v -> printf "  [%2d]: %8e\n" (i+1) v) mu

  end else Mpi.send uBdata npes 0 comm

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 *)

(* f routine. Compute f(t,u) for forward phase. *)

let f data t (udata, _, _) (dudata, _, _) =

  (* Extract MPI info. from data *)
  let comm      = data.comm in
  let npes      = data.npes in
  let my_pe     = data.my_pe in
  
  (* If this process is inactive, return now *)
  if my_pe <> npes then begin
    (* Extract problem constants from data *)
    let hordc = data.hdcoef in
    let horac = data.hacoef in

    (* Find related processes *)
    let my_pe_m1 = my_pe - 1 in
    let my_pe_p1 = my_pe + 1 in
    let last_pe  = npes - 1 in

    (* Obtain local arrays *)
    let my_length = Array1.dim udata in

    (* Pass needed data to processes before and after current process. *)
    if my_pe <> 0 then Mpi.send_float udata.{0} my_pe_m1 0 comm;
    if my_pe <> last_pe then Mpi.send_float udata.{my_length-1} my_pe_p1 0 comm;

    (* Receive needed data from processes before and after current process. *)
    let uLeft =
      if my_pe <> 0 then Mpi.receive_float my_pe_m1 0 comm else zero in
    let uRight =
      if my_pe <> last_pe then Mpi.receive_float my_pe_p1 0 comm else zero in

    (* Loop over all grid points in current process. *)
    for i=0 to my_length - 1 do
      (* Extract u at x_i and two neighboring points *)
      let ui = udata.{i} in
      let ult = if i = 0 then uLeft else udata.{i-1} in
      let urt = if i = my_length-1 then uRight else udata.{i+1} in

      (* Set diffusion and advection terms and load into udot *)
      let hdiff = hordc*.(ult -. two*.ui +. urt) in
      let hadv  = horac*.(urt -. ult) in
      dudata.{i} <- hdiff +. hadv
    done
  end

(* fB routine. Compute right hand side of backward problem *)

let fB data t (udata, _, _) (uBdata, _, _) (duBdata, _, _) =

  (* Extract MPI info. from data *)
  let comm      = data.comm in
  let npes      = data.npes in
  let my_pe     = data.my_pe in
  let my_length = Array1.dim uBdata in

  if my_pe = npes then begin (* This process performs the quadratures *)

    (* Loop over all other processes and load right hand side of quadrature eqs. *)
    duBdata.{0} <- zero;
    duBdata.{1} <- zero;
    for i=0 to npes - 1 do
      let intgr1 = Mpi.receive_float i 0 comm in
      duBdata.{0} <- duBdata.{0} +. intgr1;
      let intgr2 = Mpi.receive_float i 0 comm in
      duBdata.{1} <- duBdata.{1} +. intgr2
    done

  end else begin (* This process integrates part of the PDE *)

    (* Extract problem constants and work arrays from data *)
    let dx    = data.dx in
    let hordc = data.hdcoef in
    let horac = data.hacoef in
    let z1    = data.z1 in
    let z2    = data.z2 in

    (* Compute related parameters. *)
    let my_pe_m1 = my_pe - 1 in
    let my_pe_p1 = my_pe + 1 in
    let last_pe  = npes - 1 in

    (* Pass needed data to processes before and after current process. *)
    if my_pe <> 0 then begin
      let data_out = RealArray.of_list [ udata.{0}; uBdata.{0} ] in
      Mpi.send data_out my_pe_m1 0 comm
    end;
    if my_pe <> last_pe then begin
      let data_out = RealArray.of_list [ udata.{my_length-1}; uBdata.{my_length-1} ] in
      Mpi.send data_out my_pe_p1 0 comm
    end;
    
    (* Receive needed data from processes before and after current process. *)
    let uLeft, uBLeft =
      if my_pe <> 0 then
        let data_in = (Mpi.receive my_pe_m1 0 comm : RealArray.t) in
        data_in.{0}, data_in.{1}
      else zero, zero
    in

    let uRight, uBRight =
      if my_pe <> last_pe then
        let data_in = (Mpi.receive my_pe_p1 0 comm : RealArray.t) in
        data_in.{0}, data_in.{1}
      else zero, zero
    in

    (* Loop over all grid points in current process. *)
    for i=0 to my_length - 1 do
      
      (* Extract uB at x_i and two neighboring points *)
      let uBi = uBdata.{i} in
      let uBlt = if i = 0 then uBLeft else uBdata.{i-1} in
      let uBrt = if i = my_length-1 then uBRight else uBdata.{i+1} in
      
      (* Set diffusion and advection terms and load into udot *)
      let hdiff = hordc*.(uBlt -. two*.uBi +. uBrt) in
      let hadv = horac*.(uBrt -. uBlt) in
      duBdata.{i} <- -. hdiff +. hadv;

      (* Extract u at x_i and two neighboring points *)
      let ui = udata.{i} in
      let ult = if i = 0 then uLeft else udata.{i-1} in
      let urt = if i = my_length-1 then uRight else udata.{i+1} in

      (* Load integrands of the two space integrals *)
      z1.(i) <- uBdata.{i}*.(ult -. two*.ui +. urt)/.(dx*.dx);
      z2.(i) <- uBdata.{i}*.(urt -. ult)/.(two*.dx)
    done;

    (* Compute local integrals *)
    let intgr1 = xintgr' z1 my_length dx in
    let intgr2 = xintgr' z2 my_length dx in

    (* Send local integrals to 'quadrature' process *)
    Mpi.send_float intgr1 npes 0 comm;
    Mpi.send_float intgr2 npes 0 comm
  end

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (*------------------------------------------------------
    Initialize MPI and get total number of pe's, and my_pe
    ------------------------------------------------------*)
  let comm   = Mpi.comm_world in
  let nprocs = Mpi.comm_size comm in
  let my_pe  = Mpi.comm_rank comm in

  let npes = nprocs - 1 in (* pe's dedicated to PDE integration *)

  if npes <= 0 then begin
    if my_pe = npes then
      eprintf "\nMPI_ERROR(%d): number of processes must be >= 2\n\n" my_pe;
    exit 1
  end;

  (*-----------------------
    Set local vector length
    -----------------------*)
  let nperpe = neq/npes in
  let nrem = neq - npes*nperpe in

  let local_n, my_base =
    if my_pe < npes then
      (* PDE vars. distributed to this proccess *)
      let local_n = if my_pe < nrem then nperpe+1 else nperpe in
      local_n, (if my_pe < nrem then my_pe*local_n else my_pe*nperpe + nrem)
    else
      (* Make last process inactive for forward phase *)
      0, -1
  in

  (*-------------------------------------
    Allocate and load user data structure
    -------------------------------------*)
  let dx = xmax/.(float (mx+1)) in
  let p = [| one; 0.5 |] in
  let data = {
      p       = p;
      dx      = dx;
      hdcoef  = p.(0)/.(dx*.dx);
      hacoef  = p.(1)/.(two*.dx);
      local_n = local_n;
      npes    = npes;
      my_pe   = my_pe;
      nperpe  = nperpe;
      nrem    = nrem;
      comm    = comm;
      z1      = Array.make local_n 0.0;
      z2      = Array.make local_n 0.0;
    } in

  (*------------------------- 
    Forward integration phase
    -------------------------*)
  
  (* Set relative and absolute tolerances for forward phase *)
  let reltol = zero in
  let abstol = atol in

  (* Allocate and initialize forward variables *)
  let u = Nvector_parallel.make local_n neq comm 0.0 in
  set_ic u dx local_n my_base;

  (* Allocate CVODES memory for forward integration *)
  let cvode_mem = Cvode.init Cvode.Adams Cvode.Functional
                             (Cvode.SStolerances (reltol, abstol))
                             ~t0:t0 (f data) u
  in

  (* Allocate combined forward/backward memory *)
  Adjoint.init cvode_mem steps Adjoint.IHermite;

  (* Integrate to TOUT and collect check point information *)
  let _ = Adjoint.forward_normal cvode_mem tout u in

  (*---------------------------
    Compute and value of g(t_f)
    ---------------------------*)
  let g_val = compute_g data u in

  (*-------------------------- 
    Backward integration phase
    --------------------------*)
  (* Activate last process for integration of the quadrature equations *)
  let local_N = if my_pe = npes then np else local_n in

  (* Allocate and initialize backward variables *)
  let uB = Nvector_parallel.make local_N (neq+np) comm 0.0 in
  set_ic_back uB my_base;

  (* Allocate CVODES memory for the backward integration *)
  let bcvode_mem = Adjoint.init_backward cvode_mem
        Cvode.Adams
        Adjoint.Functional
        (Adjoint.SStolerances (reltol, abstol))
        (Adjoint.Basic (fB data)) tout uB in

  (* Integrate to T0 *)
  Adjoint.backward_normal cvode_mem t0;
  let _ = Adjoint.get bcvode_mem uB in

  (* Print results (adjoint states and quadrature variables) *)
  print_output data g_val uB

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
