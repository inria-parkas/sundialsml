(* ----------------------------------------------------------------------------- {{{
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Oct 2021.
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This demonstration problem simulates the advection and reaction of three
 * chemical species, u, v, and w, in a one dimensional domain. The reaction
 * mechanism is a variation of the Brusselator problem from chemical kinetics.
 * This is a PDE system with 3 components, Y = [u,v,w], satisfying the
 * equations,
 *
 *    u_t = -c * u_x + A - (w+1) * u + v * u^2
 *    v_t = -c * v_x + w * u - v * u^2
 *    w_t = -c * w_x + (B - w) / ep - w * u
 *
 * for t in [0, 10], x in [0, xmax] with periodic boundary conditions. The
 * initial condition is a Gaussian pertubation of the steady state
 * solution without advection
 *
 *    u(0,x) = k1 * A / k4 + p(x)
 *    v(0,x) = k2 * k4 * B / (k1 * k3 * A) + p(x)
 *    w(0,x) = 3.0 + p(x)
 *    p(x)   = alpha * e^( -(x - mu)^2 / (2*sigma^2) ).
 *
 * where alpha = 0.1, mu = xmax / 2.0, and sigma = xmax / 4.0.
 * The reaction rates are set so k_1 = k_2 = k_3 = k_4 = k, and k_5 = k_6 =
 * 1/5e-6. The spatial derivatives are discretized with first-order upwind
 * finite differences. An IMEX method is used to evolve the system in time with
 * the advection terms treated explicitly and the reaction terms implicitly. As
 * the reactions are purely local, the code uses a custom nonlinear solver to
 * exploit this locality so no parallel communication is needed in the implicit
 * solves. NOUT outputs are printed at equal intervals, and run statistics are
 * printed at the end.
 *
 * Command line options:
 *  --help           prints help message
 *  --printtime      print timing information
 *  --monitor        print solution information to screen (slower)
 *  --output-dir     the directory where all output files will be written
 *  --nout <int>     the number of output times
 *  --nx <int>       number of spatial mesh intervals
 *  --xmax <double>  maximum x value
 *  --explicit       use explicit method instead of IMEX
 *  --order <int>    method order
 *  --global-nls     use a global newton nonlinear solver instead of task-local
 *  --tf <double>    final time
 *  --A <double>     A parameter value
 *  --B <double>     B parameter value
 *  --k <double>     reaction rate
 *  --c <double>     advection speed
 *  --rtol <double>  relative tolerance
 *  --atol <double>  absolute tolerance
 * -------------------------------------------------------------------------- }}}*)

open Sundials

module ARKStep = Arkode.ARKStep
module ERKStep = Arkode.ERKStep
module NLS = NonlinearSolver

module MD = Matrix.Dense

let printf = Printf.printf
let eprintf = Printf.eprintf
let sprintf = Printf.sprintf
let fprintf = Printf.fprintf

(* Accessor macro:
   n = number of state variables
   i = mesh node index
   c = component *)
let idx n i c = n * i + c

let booli x = if x then 1 else 0

(*
 * User options structure
 *)

type user_options = {
  t0        : float;       (* initial time                 *)
  tf        : float;       (* final time                   *)
  rtol      : float;       (* relative tolerance           *)
  atol      : float;       (* absolute tolerance           *)
  order     : int;         (* method order                 *)
  explicit  : bool;        (* imex method or explicit      *)
  global    : bool;        (* use global nonlinear solve   *)
  fused     : bool;        (* use fused vector ops         *)
  nout      : int;         (* number of outputs            *)
  monitor   : bool;        (* print solution to screen     *)
  printtime : bool;        (* print timing information     *)
  tfid      : out_channel option; (* time output file pointer     *)
  ufid      : out_channel option; (* solution output file pointer *)
  vfid      : out_channel option;
  wfid      : out_channel option;
  outputdir : string;
}

(*
 * User data structure
 *)

type 't sunls =
  (Matrix.Dense.t, Nvector_serial.kind, 't) LinearSolver.serial_t

type 'k linear_system_data =
  | Explicit
  | MeshNode of {
      (* data structures for per mesh node linear system *)
      b_node   : Nvector_serial.t;
      jac_node : Nvector_serial.kind Matrix.dense;
      ls_node  : 'k sunls;
    }
  | TaskLocal of {
      (* data structures for task-local preconditioner *)
      pre   : Nvector_serial.kind Matrix.dense;
      prels : 'k sunls;
    }

type 'k user_data = {
  (* MPI data *)
  comm         : Mpi.communicator;
  myid         : int;
  nprocs       : int;
  mutable reqS : Mpi.request;
  mutable reqR : Mpi.request;
  wsend        : RealArray.t;
  esend        : RealArray.t;
  wrecv        : RealArray.t;
  erecv        : RealArray.t;

  linear_system: 'k linear_system_data;
  mutable nnlfi : int;

  (* solution masks *)
  umask : Nvector_mpiplusx.t;
  vmask : Nvector_mpiplusx.t;
  wmask : Nvector_mpiplusx.t;

  (* problem paramaters *)
  nvar  : int;     (* number of species            *)
  nx    : int;     (* number of intervals globally *)
  nxl   : int;     (* number of intervals locally  *)
  neq   : int;     (* number of equations locally  *)
  dx    : float;   (* mesh spacing                 *)
  xmax  : float;   (* maximum x value              *)
  a     : float;   (* concentration of species A   *)
  b     : float;   (* w source rate                *)
  k1    : float;   (* reaction rates               *)
  k2    : float;
  k3    : float;
  k4    : float;
  k5    : float;
  k6    : float;
  c     : float;    (* advection coefficient        *)

  (* integrator options *)
  uopt  : user_options;
}

(*
 * Definitions for a custom task local SUNNonlinearSolver
 *)

type local_nls_session =
  (Nvector_mpiplusx.data, Nvector_mpiplusx.kind) ARKStep.session
type local_nls =
  (Nvector.gdata, Nvector.gkind,
   local_nls_session, [`Nvec]) NonlinearSolver.t

type task_local_newton_content = {
  comm            : Mpi.communicator;
  myid            : Mpi.rank;
  nprocs          : int;
  mutable ncnf    : int;
  local_nls       : local_nls;
}

(* Compute the reaction term. *)
let reaction' { nvar; nxl = n; a; b; k1; k2; k3; k4; k5; k6;
               uopt = { explicit; _ }; _ } t (y  : RealArray.t) (dY : RealArray.t)
  =
  (* iterate over domain, computing reactions *)
  if explicit then
    (* when integrating explicitly, we add to ydot as we expect it
       to hold the advection term already *)
    for i = 0 to n - 1 do
      let u = y.{idx nvar i 0} in
      let v = y.{idx nvar i 1} in
      let w = y.{idx nvar i 2} in
      dY.{idx nvar i 0} <- dY.{idx nvar i 0}
               +. k1 *. a -. k2 *. w *. u +. k3 *. u *. u *. v -. k4 *. u;
      dY.{idx nvar i 1} <- dY.{idx nvar i 1}
               +. k2 *. w *. u -. k3 *. u *. u *. v;
      dY.{idx nvar i 2} <- dY.{idx nvar i 2}
               -. k2 *. w *. u +. k5 *. b -. k6 *. w
    done
  else begin
    (* set output to zero *)
    for i = 0 to n - 1 do
      let u = y.{idx nvar i 0} in
      let v = y.{idx nvar i 1} in
      let w = y.{idx nvar i 2} in
      dY.{idx nvar i 0} <-
              k1 *. a -. k2 *. w *. u +. k3 *. u *. u *. v -. k4 *. u;
      dY.{idx nvar i 1} <-
              k2 *. w *. u -. k3 *. u *. u *. v;
      dY.{idx nvar i 2} <-
              -. k2 *. w *. u +. k5 *. b -. k6 *. w;
    done
  end

let reaction udata t ((y, _)  : Nvector_mpiplusx.data)
                     ((dY, _) : Nvector_mpiplusx.data) =
  reaction' udata t (Nvector_serial.Any.unwrap y) (Nvector_serial.Any.unwrap dY)

(* --------------------------------------------------------------
 * (Non)linear system functions
 * --------------------------------------------------------------*)

let uw_any = function Nvector.RA d -> d | _ -> assert false

let task_local_nls_residual udata
                            (ycor : Nvector.gdata)
                            (f : Nvector.gdata)
                            (arkode_mem : local_nls_session)
  =
  let ARKStep.{ tcur; zpred; zi; fi; sdata; gamma }
    = ARKStep.get_nonlin_system_data arkode_mem in
  let ycor_d = uw_any ycor in
  let zpred_d = Nvector_serial.Any.unwrap (fst zpred) in
  let fi_d = Nvector_serial.Any.unwrap (fst fi) in
  let sdata_d = Nvector_serial.Any.unwrap (fst sdata) in
  let z_d = Nvector_serial.Any.unwrap (fst zi) in

  (* update 'z' value as stored predictor + current corrector *)
  Nvector_serial.(DataOps.linearsum 1.0 zpred_d 1.0 ycor_d z_d);

  (* compute implicit RHS and save for later *)
  reaction' udata tcur z_d fi_d;

  (* count calls to Fi as part of the nonlinear residual *)
  udata.nnlfi <- udata.nnlfi + 1;

  (* update with y, sdata, and gamma * fy *)
  Nvector_serial.DataOps.linearcombination
    (RealArray.of_array    [|    1.0;    -1.0; -. gamma |])
                           [| ycor_d; sdata_d;     fi_d |] (uw_any f)

let task_local_lsolve
      { nvar; nxl = n; k2; k3; k4; k6; linear_system; _ }
      delta (arkode_mem : local_nls_session)
  =
  let b_node, jac, ls = match linear_system with
    | MeshNode { b_node; jac_node; ls_node } -> b_node, jac_node, ls_node
    | _ -> assert false
  in
  let ARKStep.{ tcur; zpred; zi; fi; gamma; sdata } =
    ARKStep.get_nonlin_system_data arkode_mem
  in
  (* access solution array *)
  let bdata = uw_any delta in
  let zdata = Nvector_serial.Any.unwrap (fst zi) in

  (* solve the linear system at each mesh node *)
  for i = 0 to n - 1 do
    (* fill in Jacobian entries for this mesh node *)

    (* set nodal value shortcuts *)
    let u = zdata.{idx nvar i 0} in
    let v = zdata.{idx nvar i 1} in
    let w = zdata.{idx nvar i 2} in

    (* all vars wrt u *)
    let jdata = Matrix.unwrap jac in
    MD.set jdata 0 0 (-. k2 *. w +. 2.0 *. k3 *. u *. v -. k4);
    MD.set jdata 1 0 (   k2 *. w -. 2.0 *. k3 *. u *. v);
    MD.set jdata 2 0 (-. k2 *. w);

    (* all vars wrt v *)
    MD.set jdata 0 1 (   k3 *. u *. u);
    MD.set jdata 1 1 (-. k3 *. u *. u);
    MD.set jdata 2 1 (0.0);

    (* all vars wrt w *)
    MD.set jdata 0 2 (-. k2 *. u);
    MD.set jdata 1 2 (   k2 *. u);
    MD.set jdata 2 2 (-. k2 *. u -. k6);

    (* I - gamma*J *)
    MD.scale_addi (-. gamma) jdata;

    let b_node_data = Nvector_serial.unwrap b_node in

    (* grab just the portion of the vector 'b' for this mesh node *)
    for j = 0 to nvar - 1 do
      b_node_data.{j} <- bdata.{idx nvar i j}
    done;

    (* setup the linear system *)
    LinearSolver.setup ls jac;

    (* solve the linear system *)
    LinearSolver.solve ls jac b_node b_node 0.0;

    (* set just the portion of the vector 'b' for this mesh node *)
    for j = 0 to nvar - 1 do
      bdata.{idx nvar i j} <- b_node_data.{j}
    done
  done

(* SUNNonlinearSolver constructor *)

let task_local_newton_solve
    ({ comm; local_nls; _ } as content)
    ((y0, _)   : Nvector_mpiplusx.data)
    ((ycor, _) : Nvector_mpiplusx.data)
    ((w, _)    : Nvector_mpiplusx.data)
    tol
    callLSetup
    s
  =
  (* each tasks solves the local nonlinear system *)
  let solve_status =
    try
      NLS.solve local_nls ~y0 ~ycor ~w tol callLSetup s; 1
    with RecoverableFailure -> 2 | _ -> 0
  in
  (* if any process had a nonrecoverable failure, return it *)
  if Mpi.(allreduce_int solve_status Min comm) = 0
  then failwith "task_localnewton_solve: nonrecover";

  (* check if any process has a recoverable convergence failure *)
  if Mpi.(allreduce_int solve_status Max comm) = 2
  then (content.ncnf <- content.ncnf + 1; raise RecoverableFailure)

let rewrap_mpiplusx comm (y : Nvector.gdata) =
  match y with
  | Nvector.RA ydata -> Nvector_mpiplusx.wrap comm (Nvector_serial.Any.wrap ydata)
  | _ -> assert false

let task_local_newton_initialize udata { local_nls; comm; _ } () =
  NLS.set_sys_fn local_nls (task_local_nls_residual udata);
  NLS.set_lsolve_fn local_nls (task_local_lsolve udata)

let task_local_newton_setsysfn { local_nls; comm; _ }
      (sysfn : (Nvector_mpiplusx.t, local_nls_session) NLS.sysfn) =
  let rw = rewrap_mpiplusx comm in
  NLS.set_sys_fn local_nls (fun y fg -> sysfn (rw y) (rw fg))

let task_local_newton_setconvtestfn { local_nls; comm; _ }
      (ctestfn : (Nvector_mpiplusx.data, local_nls_session, [`Nvec]) NLS.convtestfn) =
  NLS.(set_convtest_fn local_nls
         ((assert_not_oconvtestfn ctestfn)
            : (Nvector.gdata, local_nls_session, [`Nvec]) NLS.convtestfn))

(* Pass via OCaml (for testing) *)
let task_local_newton_setconvtestfn' { local_nls; comm; _ }
      (ctestfn : (Nvector_mpiplusx.data, local_nls_session, [`Nvec]) NLS.convtestfn) =
  let rewrap = function
    | Nvector.RA a -> Nvector_serial.wrap a
    | _ -> assert false
  in
  let ocaml_convtestfn y del tol ewt mem =
    match ctestfn with
    | CConvTest cfn -> (Sundials.invoke cfn).f local_nls
                         (rewrap y) (rewrap del) tol (rewrap ewt) mem
    | _ -> assert false
  in
  NLS.(set_convtest_fn local_nls
         ((OConvTest ocaml_convtestfn)
            : (Nvector.gdata, local_nls_session, [`Nvec]) NLS.convtestfn))

let task_local_newton_getnumconvfails { ncnf; _ } () = ncnf

type nlsolver =
  (Nvector_mpiplusx.data, Nvector_mpiplusx.kind, local_nls_session, [`Nvec])
    NonlinearSolver.t

let task_local_newton udata (y : Nvector_mpiplusx.t) dfid =
  let yd, comm = Nvector_mpiplusx.unwrap y in
  let local_nls = NLS.Newton.make yd
  in
  let content = {
    comm;
    local_nls;
    myid      = Mpi.comm_rank comm;
    nprocs    = Mpi.comm_size comm;
    ncnf      = 0;
  } in
  (* Setup the local nonlinear solver monitoring *)
  (match dfid with None -> ()
   | Some dfid -> NLS.set_info_file local_nls dfid;
                  NLS.set_print_level local_nls true);
  (NLS.Custom.(make
                ~init:(task_local_newton_initialize udata content)
                ~set_convtest_fn:(task_local_newton_setconvtestfn content)
                ~get_num_conv_fails:(task_local_newton_getnumconvfails content)
                ~nls_type:RootFind
                ~solve:(task_local_newton_solve content)
                ~set_sys_fn:(task_local_newton_setsysfn content)
                ()) : nlsolver)

(* --------------------------------------------------------------
 * Right Hand Side (RHS) Functions
 * --------------------------------------------------------------*)

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 0) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 1) []) 0
  - header_and_empty_array_size

let bytes x = header_and_empty_array_size + x * float_cell_size

(* Starts the exchange of the neighbor information *)
let exchange_all_start
    ({ comm; c; nxl = n; nvar; myid; nprocs; esend; wsend; _ } as udata) ydata =

  (* shortcuts *)
  let first = 0 in
  let last = nprocs - 1 in
  let ipW = if myid = first then  last else myid - 1 in (* periodic BC *)
  let ipE = if myid = last  then first else myid + 1 in (* periodic BC *)

  if c > 0.0 then begin
    (* Right moving flow uses backward difference.
       Send from west to east (last processor sends to first) *)
    udata.reqR <- Mpi.ireceive (bytes nvar) ipW Mpi.any_tag comm; (* wrecv *)

    for var = 0 to nvar - 1 do
      esend.{idx nvar 0 var} <- ydata.{idx nvar (n - 1) var}
    done;

    udata.reqS <- Mpi.isend esend ipE 0 comm
  end else if c < 0.0 then begin
    (* Left moving flow uses forward difference.
       Send from east to west (first processor sends to last) *)
    udata.reqR <- Mpi.ireceive (bytes nvar) ipE Mpi.any_tag comm; (* erecv *)

    for var = 0 to nvar - 1 do
      wsend.{idx nvar 0 var} <- ydata.{idx nvar 0 var}
    done;

    udata.reqS <- Mpi.isend wsend ipW 0 comm
  end

(* Completes the exchange of the neighbor information *)
let exchange_all_end { c; reqR; reqS; wrecv; erecv; _ } =
  (* wait for exchange to finish *)
  if c > 0.0 then begin
    RealArray.blit ~src:(Mpi.wait_receive reqR) ~dst:wrecv;
    Mpi.wait reqS
  end else if c < 0.0 then begin
    RealArray.blit ~src:(Mpi.wait_receive reqR) ~dst:erecv;
    Mpi.wait reqS
  end

(* Compute the advection term. *)
let advection' ({ nvar; nxl = n; dx; c; erecv; wrecv; _ } as udata) t
               (y : RealArray.t) (dY : RealArray.t)
  =
  (* set output to zero *)
  RealArray.fill dY 0.0;

  (* begin exchanging boundary information *)
  exchange_all_start udata y;

  (* iterate over domain interior, computing advection *)
  let tmp = -. c /. dx in

  if c > 0.0 then
    (* right moving flow *)
    for i = 1 to n - 1 do
      for var = 0 to nvar - 1 do
        dY.{idx nvar i var} <-
          tmp *. (y.{idx nvar i var} -. y.{idx nvar (i-1) var})
      done
    done
  else if c < 0.0 then
    (* left moving flow *)
    for i = 0 to n - 2 do
      for var = 0 to nvar - 1 do
        dY.{idx nvar i var} <-
          tmp *. (y.{idx nvar (i + 1) var} -. y.{idx nvar i var})
      done
    done;

  (* finish exchanging boundary information *)
  exchange_all_end udata;

  (* compute advection at local boundaries *)
  if c > 0.0 then
    (* right moving flow (left boundary) *)
    for var = 0 to nvar - 1 do
      dY.{idx nvar 0 var} <-
        tmp *. (y.{idx nvar 0 var} -. wrecv.{idx nvar 0 var});
    done
  else if (c < 0.0) then
    (* left moving flow (right boundary) *)
    for var = 0 to nvar - 1 do
      dY.{idx nvar (n - 1) var} <-
        tmp *. (erecv.{idx nvar 0 var} -. y.{idx nvar (n - 1) var})
    done

let advection udata t ((y, _) : Nvector_mpiplusx.data) ((dY, _) : Nvector_mpiplusx.data)
  = advection' udata t (Nvector_serial.Any.unwrap y) (Nvector_serial.Any.unwrap dY)

(* Compute the RHS as Advection+Reaction. *)
let advection_reaction user_data t
    (y : Nvector_mpiplusx.data) (ydot : Nvector_mpiplusx.data) =
  (* NOTE: The order in which Advection and Reaction are
           called is critical here. Advection must be
           computed first. *)
  advection user_data t y ydot;
  reaction user_data t y ydot

(* --------------------------------------------------------------
 * Preconditioner functions
 * --------------------------------------------------------------*)

(* Sets P = I - gamma * J *)
let psetup { nvar; nxl = n; k2; k3; k4; k6; linear_system; _ }
           ARKStep.{ jac_y = ((y, _) : Nvector_mpiplusx.data); _ } jok gamma
  =
  let p, ls = match linear_system with
              | TaskLocal { pre; prels } -> pre, prels
              | _ -> assert false
  in
  if jok then false
  else begin
    let ydata = Nvector_serial.Any.unwrap y in
    let pdata = Matrix.unwrap p in
    (* setup the block diagonal preconditioner matrix *)
    for i = 0 to n - 1 do
      (* fill in Jacobian entries for this mesh node *)
      let blocki = nvar*i in

      (* set nodal value shortcuts *)
      let u = ydata.{idx nvar i 0} in
      let v = ydata.{idx nvar i 1} in
      let w = ydata.{idx nvar i 2} in

      (* all vars wrt u *)
      MD.set pdata (blocki  ) blocki (-. k2 *. w +. 2.0 *. k3 *. u *. v -. k4);
      MD.set pdata (blocki+1) blocki (   k2 *. w -. 2.0 *. k3 *. u *. v);
      MD.set pdata (blocki+2) blocki (-. k2 *. w);

      (* all vars wrt v *)
      MD.set pdata (blocki  ) (blocki+1) (   k3 *. u *. u);
      MD.set pdata (blocki+1) (blocki+1) (-. k3 *. u *. u);
      MD.set pdata (blocki+2) (blocki+1) 0.0;

      (* all vars wrt w *)
      MD.set pdata (blocki  ) (blocki+2) (-. k2 *. u);
      MD.set pdata (blocki+1) (blocki+2) (   k2 *. u);
      MD.set pdata (blocki+2) (blocki+2) (-. k2 *. u -. k6)
    done;
    MD.scale_addi (-. gamma) pdata;

    (* setup the linear system Pz = r *)
    LinearSolver.setup ls p;

    (* indicate that J is now current *)
    true
  end

(* Solves Pz = r *)
let psolve { linear_system; _ }
           ARKStep.{ jac_y = ((y, _) : Nvector_mpiplusx.data); _ }
           ARKStep.Spils.{ rhs = ((r, _) : Nvector_mpiplusx.data); _ }
           ((z, _) : Nvector_mpiplusx.data) =
  let p, ls = match linear_system with
              | TaskLocal { pre; prels } -> pre, prels
              | _ -> assert false
  in
  let z_local = Nvector_serial.(wrap (Any.unwrap z)) in
  let r_local = Nvector_serial.(wrap (Any.unwrap r)) in

  (* solve the task-local linear system Pz = r *)
  LinearSolver.solve ls p z_local r_local 0.0

(* --------------------------------------------------------------
 * Utility functions
 * --------------------------------------------------------------*)

(* Exchanges the periodic BCs only by sending the first
   mesh node to the last processor. *)
let exchange_bc_only ({ nvar; myid; nprocs; wsend; erecv; comm; _ } as udata)
    (y : Nvector_mpiplusx.t)  =
  (* shortcuts *)
  let first = 0 in
  let last = nprocs - 1 in

  (* extract the data *)
  let ydata = Nvector_serial.Any.unwrap (fst (Nvector_mpiplusx.unwrap y)) in

  (* open the East Irecv buffer *)
  if myid = last then
    udata.reqR <- Mpi.ireceive (bytes nvar) first Mpi.any_tag comm; (* erecv *)

  (* send first mesh node to the last processor *)
  if myid = first then begin
    for var = 0 to nvar - 1 do
      wsend.{idx nvar 0 var} <- ydata.{idx nvar 0 var}
    done;
    udata.reqS <- Mpi.isend wsend last 0 comm
  end;

  (* wait for exchange to finish *)
  if myid = last then
    RealArray.blit ~src:(Mpi.wait_receive udata.reqR) ~dst:erecv;

  if myid = first then Mpi.wait udata.reqS

let setup_problem () =
  let monitor    = ref false in
  let printtime  = ref false in
  let nout       = ref 40 in
  let nx         = ref 100 in
  let xmax       = ref 1.0 in
  let a          = ref 1.0 in
  let b          = ref 3.5 in
  let c          = ref 0.01 in
  let k1         = ref 1.0 in
  let k2         = ref 1.0 in
  let k3         = ref 1.0 in
  let k4         = ref 1.0 in
  let order      = ref 3 in
  let explicit   = ref false in
  let global     = ref false in
  let fused      = ref false in
  let tf         = ref 10.0 in
  let rtol       = ref 1.0e-6 in
  let atol       = ref 1.0e-9 in
  let outputdir  = ref "." in
  let args = Arg.[
    "--monitor", Set monitor,
    "print solution information to screen (slower)\n";

    "--output-dir", Set_string outputdir,
    "the directory where all output files will be written\n";

    "--nout", Set_int nout,
    "<int> number of output times\n";

    "--explicit", Set explicit,
    "use an explicit method instead of IMEX\n";

    "--global-nls", Set global,
    "use a global newton nonlinear solver instead of task-local (for IMEX only)\n";

    "--order", Set_int order,
    "<int> the method order to use\n";

    "--nx", Set_int nx,
    "<int> number of mesh points\n";

    "--xmax", Set_float xmax,
    "<double> maximum value of x (size of domain)\n";

    "--tf", Set_float tf,
    "<double> final time\n";

    "--A", Set_float a,
    "<double> A parameter value\n";

    "--B", Set_float b,
    "<double> B parameter value\n";

    "--k",
    Unit (fun () ->
            k1 := Float.of_string Sys.argv.(!Arg.current);
            k2 := Float.of_string Sys.argv.(!Arg.current + 1);
            k3 := Float.of_string Sys.argv.(!Arg.current + 2);
            k4 := Float.of_string Sys.argv.(!Arg.current + 3);
            Arg.current := !Arg.current + 4),
    "<double> reaction rate\n";

    "--c", Set_float c,
    "<double> advection speed\n";

    "--rtol", Set_float rtol,
    "<double> relative tolerance\n";

    "--atol", Set_float atol,
    "<double> absolute tolerance\n";

    "--printtime", Set printtime,
    "Print timing";

    "-fused", Set fused,
    "Enable the nvector fused operations";
  ] in
  Arg.parse args (fun _ -> ()) ("Command line options for " ^ Sys.argv.(0));

  let comm = Mpi.comm_world in
  let nvar = 3 in
  let nprocs = Mpi.comm_size comm in
  let myid = Mpi.comm_rank comm in
  let nx = !nx in
  let nxl = nx / nprocs in
  let neq = nvar * nxl in
  if nx mod nprocs <> 0 then failwith
    (sprintf "ERROR: The mesh size (nx = %d) must be divisible by the \
              number of processors (%d)\n" nx nprocs);
  let umask_nv = Nvector_serial.Any.make ~with_fused_ops:!fused neq 0.0 in
  let umask = Nvector_mpiplusx.wrap ~with_fused_ops:!fused comm umask_nv in
  let open_file name = open_out (sprintf "%s/%s.%06d.txt" !outputdir name myid) in
  let udata = {
    (* MPI data *)
    comm; myid; nprocs;
    reqS  = Mpi.null_request;
    reqR  = Mpi.null_request;
    wsend = RealArray.make nvar 0.0;
    esend = RealArray.make nvar 0.0;
    wrecv = RealArray.make nvar 0.0;
    erecv = RealArray.make nvar 0.0;

    linear_system =
      if !explicit then Explicit
      else if !global then begin
        (* Create MPI task-local data structures for preconditioning *)
        let pre = Matrix.dense neq in
        TaskLocal {
          pre;
          prels = LinearSolver.Direct.dense (Nvector_serial.make neq 0.0) pre
        }
      end else begin
        let b_node = Nvector_serial.make nvar 0.0 in
        let jac_node = Matrix.dense nvar in
        MeshNode {
          (* Create MPI task-local data structures for mesh node solves *)
          b_node;
          jac_node;
          ls_node  = LinearSolver.Direct.dense b_node jac_node
        }
      end;

    nnlfi = 0;

    (* solution masks *)
    umask;
    vmask = Nvector.clone umask;
    wmask = Nvector.clone umask;

    (* problem paramaters *)
    nvar; nx; nxl; neq;

    dx    = !xmax /. float nx;
    xmax  = !xmax;
    a     = !a;
    b     = !b;
    k1    = !k1;
    k2    = !k2;
    k3    = !k3;
    k4    = !k4;
    k5    = 1.0 /. 5.0e-6;
    k6    = 1.0 /. 5.0e-6;
    c     = !c;

    (* integrator options *)
    uopt = {
      t0        = 0.0;
      tf        = !tf;
      rtol      = !rtol;
      atol      = !atol;
      order     = !order;
      explicit  = !explicit;
      global    = !global;
      fused     = !fused;
      nout      = !nout;
      monitor   = !monitor;
      printtime = !printtime;
      tfid      = if !nout > 0 && myid = 0 then Some (open_file "t") else None;
      ufid      = if !nout > 0 then Some (open_file "u") else None;
      vfid      = if !nout > 0 then Some (open_file "v") else None;
      wfid      = if !nout > 0 then Some (open_file "w") else None;
      outputdir = !outputdir;
    }
  } in
  (* Create the solution masks *)
  let init_mask mask j =
    Nvector_mpiplusx.Ops.const 0.0 mask;
    let data = Nvector_serial.Any.unwrap (fst (Nvector_mpiplusx.unwrap mask)) in
    for i = 0 to udata.nxl - 1 do
      data.{idx nvar i j} <- 1.0
    done
  in
  init_mask udata.umask 0;
  init_mask udata.vmask 1;
  init_mask udata.wmask 2;

  (* Print problem setup *)
  let uopt = udata.uopt in
  if myid = 0 then begin
    printf "\n1D Advection-Reaction Test Problem\n\n";
    printf "Number of Processors = %d\n" udata.nprocs;
    printf "Mesh Info:\n";
    printf "  NX = %d, NXL = %d, dx = %.6f, xmax = %.6f\n"
           udata.nx udata.nxl udata.dx udata.xmax;
    printf "Problem Parameters:\n";
    printf "  A = %g\n" udata.a;
    printf "  B = %g\n" udata.b;
    printf "  k = %g\n" udata.k1;
    printf "  c = %g\n" udata.c;

    printf "Integrator Options:\n";
    printf "  order            = %d\n" uopt.order;
    printf "  method           = %s\n"
           (if uopt.explicit then "explicit" else "imex");
    printf "  fused vector ops = %d\n" (booli uopt.fused);
    printf "  t0               = %g\n" uopt.t0;
    printf "  tf               = %g\n"uopt.tf;
    printf "  reltol           = %.1e\n" uopt.rtol;
    printf "  abstol           = %.1e\n" uopt.atol;
    printf "  nout             = %d\n" uopt.nout;
    if uopt.explicit then printf "  nonlinear solver = none\n"
    else printf "  nonlinear solver = %s\n"
           (if uopt.global then "global" else "task local");
    printf "Output directory: %s\n" uopt.outputdir
  end;
  udata

let oout f = function None -> () | Some oc -> f oc

(* Write time and solution to disk *)
let write_output ({ myid; uopt; umask; vmask; wmask; nvar;
                    nx; nxl; nprocs; erecv; _ } as udata)
                 t (y : Nvector_mpiplusx.t) =
  let { monitor; nout; tfid; ufid; vfid; wfid; _ } = uopt in
  (* output current solution norm to screen *)
  if monitor then begin
    let u = Nvector_mpiplusx.Ops.wl2norm y umask in
    let u = sqrt (u *. u /. float nx) in
    let v = Nvector_mpiplusx.Ops.wl2norm y vmask in
    let v = sqrt (v *. v /. float nx) in
    let w = Nvector_mpiplusx.Ops.wl2norm y wmask in
    let w = sqrt (w *. w /. float nx) in
    if myid = 0 then printf "     %10.6f   %10.6f   %10.6f   %10.6f\n" t u v w
  end;

  if nout > 0 then begin
    (* get left end point for output *)
    exchange_bc_only udata y;

    (* get vector data array *)
    let data = Nvector_serial.Any.unwrap (fst (Nvector_mpiplusx.unwrap y)) in

    (* output the times to disk *)
    if myid = 0 then oout (fun o -> fprintf o " %.16e\n" t) tfid;

    (* output results to disk *)
    for i = 0 to nxl - 1 do
      oout (fun o -> fprintf o " %.16e" data.{idx nvar i 0}) ufid;
      oout (fun o -> fprintf o " %.16e" data.{idx nvar i 1}) vfid;
      oout (fun o -> fprintf o " %.16e" data.{idx nvar i 2}) wfid
    done;

    (* we have one extra output because of the periodic BCs *)
    if myid = nprocs - 1 then begin
      oout (fun o -> fprintf o " %.16e\n" erecv.{idx nvar 0 0}) ufid;
      oout (fun o -> fprintf o " %.16e\n" erecv.{idx nvar 0 1}) vfid;
      oout (fun o -> fprintf o " %.16e\n" erecv.{idx nvar 0 2}) wfid
    end else begin
      oout (fun o -> fprintf o "\n") ufid;
      oout (fun o -> fprintf o "\n") vfid;
      oout (fun o -> fprintf o "\n") wfid
    end
  end

(* Setup ARKODE and evolve problem in time with IMEX method*)
let evolve_problem_imex ({ myid; uopt; _ } as udata) (y : Nvector_mpiplusx.t) =
  let { rtol; atol; order; t0; tf; nout; monitor; global; outputdir; _ } = uopt
  in
  let dfid = if monitor
    then Some (Logfile.openfile
                 (sprintf "%s/diagnostics.%06d.txt" outputdir myid))
    else None
  in
  (* Create the (non)linear solver *)
  let (nlsolver : nlsolver), lsolver =
    if global then
      NonlinearSolver.Newton.make y,
      Some ARKStep.Spils.(solver
                            LinearSolver.Iterative.(spgmr y)
                            (prec_left ~setup:(psetup udata) (psolve udata)))
    else
      task_local_newton udata y dfid,
      (* The custom task-local nonlinear solver handles the linear solve
         as well, so we do not need a SUNLinearSolver *)
      None
  in
  (* Create the ARK timestepper module *)
  let arkode_mem =
    ARKStep.(init (imex ~nlsolver ?lsolver
                        ~fi:(reaction udata) (advection udata))
                  (SStolerances (rtol, atol))
                  ~order t0 y)
  in
  (* Increase the max number of steps allowed between outputs *)
  ARKStep.set_max_num_steps arkode_mem 100000;

  (* Open output file for integrator diagnostics *)
  (match dfid with None -> ()
   | Some dfid -> ARKStep.set_diagnostics ~logfile:dfid arkode_mem);

  (* Output initial condition *)
  if myid = 0 && monitor then begin
    printf "\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n";
    printf "   ----------------------------------------------------\n";
  end;
  write_output udata t0 y;

  (* Integrate to final time *)
  let dtout = if nout <> 0 then (tf -. t0) /. float nout else (tf -. t0) in
  let rec loop iout tout =
    (* Integrate to output time *)
    let t, _ = ARKStep.evolve_normal arkode_mem tout y in
    write_output udata t y;

    (* Update output time *)
    let tout = min tf (tout +. dtout) in

    if iout + 1 < nout then loop (iout + 1) tout
  in
  loop 0 (t0 +. dtout);

  (* close output stream *)
  (match dfid with None -> () | Some dfid -> Logfile.close dfid);

  (* Get final statistics *)
  let open ARKStep in
  let nst      = get_num_steps arkode_mem
  and nst_a    = get_num_step_attempts arkode_mem
  and nfe, nfi = get_num_rhs_evals arkode_mem
  and netf    = get_num_err_test_fails arkode_mem
  and nni     = get_num_nonlin_solv_iters arkode_mem
  and ncnf    = get_num_nonlin_solv_conv_fails arkode_mem
  in
  (* Print final statistics *)
  if myid = 0 then begin
    printf "\nFinal Solver Statistics (for processor 0):\n";
    printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
    printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe (nfi + udata.nnlfi);
    printf "   Total number of error test failures = %d\n" netf;
    printf "   Total number of nonlinear solver convergence failures = %d\n" ncnf;

    if global then begin
      let nli = Spils.get_num_lin_iters arkode_mem in
      let npre = Spils.get_num_prec_evals arkode_mem in
      let npsol = Spils.get_num_prec_solves arkode_mem in
      printf "   Total number of nonlinear iterations = %d\n" nni;
      printf "   Total number of linear iterations = %d\n" nli;
      printf "   Total number of preconditioner setups = %d\n" npre;
      printf "   Total number of preconditioner solves = %d\n" npsol
    end
  end

(* Setup ARKODE and evolve problem in time explicitly *)
let evolve_problem_explicit ({ myid; uopt; _ } as udata) (y : Nvector_mpiplusx.t) =
  let { rtol; atol; order; t0; tf; nout;
        monitor; global; outputdir; _ } = udata.uopt
  in
  let dfid = if monitor
    then Some (Logfile.openfile
                 (sprintf "%s/diagnostics.%06d.txt" outputdir myid))
    else None
  in

  (* Create the ERK timestepper module *)
  let arkode_mem = ERKStep.(init (SStolerances (rtol, atol))
                            ~order (advection_reaction udata) t0 y)
  in
  (* Increase the max number of steps allowed between outputs *)
  ERKStep.set_max_num_steps arkode_mem 1000000;

  (* Open output file for integrator diagnostics *)
  (match dfid with None -> ()
   | Some dfid -> ERKStep.set_diagnostics arkode_mem dfid);

  (* Output initial condition *)
  if myid = 0 && monitor then begin
    printf "\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n";
    printf "   ----------------------------------------------------\n";
  end;
  write_output udata t0 y;

  (* Integrate to final time *)
  let dtout = if nout <> 0 then (tf -. t0) /. float nout else (tf -. t0) in
  let rec loop iout tout =
    (* Integrate to output time *)
    let t, _ = ERKStep.evolve_normal arkode_mem tout y in
    write_output udata t y;

    (* Update output time *)
    let tout = min tf (tout +. dtout) in

    if iout + 1 < nout then loop (iout + 1) tout
  in
  loop 0 (t0 +. dtout);

  (* close output stream *)
  (match dfid with None -> () | Some dfid -> Logfile.close dfid);

  (* Get final statistics *)
  let open ERKStep in
  let nst     = get_num_steps arkode_mem
  and nst_a   = get_num_step_attempts arkode_mem
  and nfe     = get_num_rhs_evals arkode_mem
  and netf    = get_num_err_test_fails arkode_mem
  in
  (* Print final statistics *)
  if myid = 0 then begin
    printf "\nFinal Solver Statistics (for processor 0):\n";
    printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
    printf "   Total RHS evals:  Fe = %d" nfe;
    printf "   Total number of error test failures = %d\n" netf;
  end

(* Initial Condition Functions *)
let set_ic { nvar; nxl = n; dx; a; b; k1; k2; k3; k4; myid; xmax; _ }
           (y : Nvector_mpiplusx.t) =
  (* Gaussian distribution defaults *)
  let mu    = xmax /. 2.0 in
  let sigma = xmax /. 4.0 in
  let alpha = 0.1 in

  (* Access data array from NVector y *)
  let data = Nvector_serial.Any.unwrap (fst (Nvector.unwrap y)) in

  (* Steady state solution *)
  let us = k1 *. a /. k4 in
  let vs = k2 *. k4 *. b /. (k1 *. k3 *. a) in
  let ws = 3.0 in

  (* Gaussian perturbation of the steady state solution *)
  for i = 0 to n - 1 do
    let x = float (myid * n + i) *. dx in
    let p = alpha *. exp(-.((x -. mu) *. (x -. mu)) /. (2.0 *. sigma *. sigma)) in
    data.{idx nvar i 0} <- us +. p;
    data.{idx nvar i 1} <- vs +. p;
    data.{idx nvar i 2} <- ws +. p
  done

(* Main Program *)
let main () =
  (* Start timing *)
  let starttime = Mpi.wtime () in

  (* Process input args and setup the problem *)
  let { myid; neq; comm; dx; nx;
        uopt = { fused; explicit; nout; outputdir; printtime; _ } as uopt;
    _ } as udata = setup_problem ()
  in

  (* Create solution vector *)
  let ys = Nvector_serial.Any.make ~with_fused_ops:fused neq 0.0 in
  let y = Nvector_mpiplusx.wrap ~with_fused_ops:fused comm ys in

  (* Set the initial condition *)
  set_ic udata y;

  (* Output spatial mesh to disk (add extra point for periodic BC) *)
  if myid = 0 && nout > 0 then begin
    let mfid = open_out (sprintf "%s/mesh.txt" outputdir) in
    for i = 0 to nx do
      fprintf mfid "  %.16e\n" (dx *. float i)
    done;
    close_out mfid
  end;

  (* Integrate in time *)
  if explicit then evolve_problem_explicit udata y
  else evolve_problem_imex udata y;

  (* End timing *)
  let endtime = Mpi.wtime() in
  if myid = 0 && printtime then
    printf "\nTotal wall clock time: %.4f seconds\n" (endtime -. starttime)

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

