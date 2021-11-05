(*--------------------------------------------------------------- {{{
 * Programmer(s): Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0,
 *    v_t(t,0) = v_t(t,1) = 0,
 *    w_t(t,0) = w_t(t,1) = 0.
 * Note: these can also be implemented as Dirichlet boundary
 * conditions with values identical to the initial conditions.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * The data is stored using the ManyVector structure for a
 * structure-of-arrays format, i.e., each of u, v and w are
 * stored in separate serial vectors, and are attached together
 * using the ManyVector infrastructure.
 *
 * This program solves the problem with the ARK method, treating
 * only the reaction terms implicitly (diffusion is treated
 * explicitly), using a Newton iteration with the SUNSPGMR
 * iterative linear solver, and a user-supplied
 * Jacobian-vector-product routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *--------------------------------------------------------------- }}} *)

open Sundials
module ARKStep = Arkode.ARKStep

let printf = Printf.printf
let fprintf = Printf.fprintf
let unwrap = Nvector.unwrap

(* realtype constant macros *)
let zero = 0.0
let one  = 1.0
let two  = 2.0

(* user data structure *)
type user_data = {
    n  : int;     (* number of intervals     *)
    dx : float;   (* mesh spacing            *)
    a  : float;   (* constant forcing on u   *)
    b  : float;   (* steady-state value of w *)
    du : float;   (* diffusion coeff for u   *)
    dv : float;   (* diffusion coeff for v   *)
    dw : float;   (* diffusion coeff for w   *)
    ep : float;   (* stiffness parameter     *)
  }

(* fe routine to compute the diffusion portion of the ODE RHS. *)
let fe { n; du; dv; dw; dx; _ }
       t ((y : Nvector.any ROArray.t), _)
         (((dy : Nvector.any ROArray.t), _) as nvdy) =
  let y_u = Nvector_serial.Any.unwrap (ROArray.get y 0) in
  let y_v = Nvector_serial.Any.unwrap (ROArray.get y 1) in
  let y_w = Nvector_serial.Any.unwrap (ROArray.get y 2) in

  let f_u = Nvector_serial.Any.unwrap (ROArray.get dy 0) in
  let f_v = Nvector_serial.Any.unwrap (ROArray.get dy 1) in
  let f_w = Nvector_serial.Any.unwrap (ROArray.get dy 2) in

  Nvector_many.DataOps.const 0.0 nvdy; (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  let uconst = du /. dx /. dx in
  let vconst = dv /. dx /. dx in
  let wconst = dw /. dx /. dx in

  for i = 1 to n - 2 do
    (* Fill in ODE RHS for u *)
    f_u.{i} <- (y_u.{i-1} -. two*.y_u.{i} +. y_u.{i+1})*.uconst;

    (* Fill in ODE RHS for v *)
    f_v.{i} <- (y_v.{i-1} -. two*.y_v.{i} +. y_v.{i+1})*.vconst;

    (* Fill in ODE RHS for w *)
    f_w.{i} <- (y_w.{i-1} -. two*.y_w.{i} +. y_w.{i+1})*.wconst
  done;

  (* enforce stationary boundaries *)
  f_u.{0}   <- zero;
  f_v.{0}   <- zero;
  f_w.{0}   <- zero;
  f_u.{n-1} <- zero;
  f_v.{n-1} <- zero;
  f_w.{n-1} <- zero

(* fi routine to compute the reaction portion of the ODE RHS. *)
let fi { n; a; b; ep; _ }
       t ((y : Nvector.any ROArray.t), _)
         (((dy : Nvector.any ROArray.t), _) as nvdy) =
  let y_u = Nvector_serial.Any.unwrap (ROArray.get y 0) in
  let y_v = Nvector_serial.Any.unwrap (ROArray.get y 1) in
  let y_w = Nvector_serial.Any.unwrap (ROArray.get y 2) in

  let f_u = Nvector_serial.Any.unwrap (ROArray.get dy 0) in
  let f_v = Nvector_serial.Any.unwrap (ROArray.get dy 1) in
  let f_w = Nvector_serial.Any.unwrap (ROArray.get dy 2) in

  Nvector_many.DataOps.const 0.0 nvdy; (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  for i = 1 to n-2 do
    (* Fill in ODE RHS for u *)
    f_u.{i} <- a -. (y_w.{i}+.one)*.y_u.{i} +. y_v.{i}*.y_u.{i}*.y_u.{i};

    (* Fill in ODE RHS for v *)
    f_v.{i} <- y_w.{i}*.y_u.{i} -. y_v.{i}*.y_u.{i}*.y_u.{i};

    (* Fill in ODE RHS for w *)
    f_w.{i} <- (b-.y_w.{i})/.ep -. y_w.{i}*.y_u.{i};
  done;

  (* enforce stationary boundaries *)
  f_u.{0}   <- zero;
  f_v.{0}   <- zero;
  f_w.{0}   <- zero;
  f_u.{n-1} <- zero;
  f_v.{n-1} <- zero;
  f_w.{n-1} <- zero

(* Jacobian-vector product routine (implicit portion only) *)
let jacvi { n; ep; _ }
          { ARKStep.jac_y = ((y : Nvector.any ROArray.t), _) }
          ((v, _) : Nvector_many.data)
          (((jv, _) : Nvector_many.data) as nvjv) =
  let y_u = Nvector_serial.Any.unwrap (ROArray.get y 0) in
  let y_v = Nvector_serial.Any.unwrap (ROArray.get y 1) in
  let y_w = Nvector_serial.Any.unwrap (ROArray.get y 2) in

  let v_u = Nvector_serial.Any.unwrap (ROArray.get v 0) in
  let v_v = Nvector_serial.Any.unwrap (ROArray.get v 1) in
  let v_w = Nvector_serial.Any.unwrap (ROArray.get v 2) in

  let ju_v = Nvector_serial.Any.unwrap (ROArray.get jv 0) in
  let jv_v = Nvector_serial.Any.unwrap (ROArray.get jv 1) in
  let jw_v = Nvector_serial.Any.unwrap (ROArray.get jv 2) in

  Nvector_many.DataOps.const 0.0 nvjv; (* initialize ydot to zero *)

  (* iterate over domain, computing Jacobian-vector products *)
  for i = 1 to n - 2 do
    (* Fill in Jacobian-vector product for Ju*v *)
    ju_v.{i} <- -. v_w.{i} *. y_u.{i}
                -. y_w.{i} *. v_u.{i}
                -. v_u.{i}
                +. v_v.{i} *. y_u.{i} *. y_u.{i}
                +. two *. y_v.{i} *. y_u.{i} *. v_u.{i};

    (* Fill in Jacobian-vector product for Jv*v *)
    jv_v.{i} <- v_w.{i} *. y_u.{i}
                +. y_w.{i} *. v_u.{i}
                -. v_v.{i} *. y_u.{i} *. y_u.{i}
                -. two *. y_v.{i} *. y_u.{i} *. v_u.{i};

    (* Fill in Jacobian-vector product for Jw*v *)
    jw_v.{i} <- -. v_w.{i} /. ep
                -. v_w.{i} *. y_u.{i}
                -. y_w.{i} *. v_u.{i}
  done;

  (* enforce stationary boundaries *)
  ju_v.{0}   <- zero;
  jv_v.{0}   <- zero;
  jw_v.{0}   <- zero;
  ju_v.{n-1} <- zero;
  jv_v.{n-1} <- zero;
  jw_v.{n-1} <- zero

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in    (* initial time *)
  let tf = 10.0 in   (* final time *)
  let nt = 100 in    (* total number of output times *)
  let n = 201 in

  let reltol = 1.0e-6 in     (* tolerances *)
  let abstol = 1.0e-10 in

  (* store the inputs in the UserData structure *)
  let ud = {
    n;
    dx = 1.0 /. float(n-1);
    a = 0.6;
    b = 2.0;
    du = 0.001;
    dv = 0.001;
    dw = 0.001;
    ep = 1.0e-5         (* stiffness parameter *)
  } in

  (* Initial problem output *)
  printf "\n1D Brusselator PDE test problem:\n";
  printf "    N = %d\n" ud.n;
  printf "    problem parameters:  a = %g,  b = %g,  ep = %g\n" ud.a  ud.b  ud.ep;
  printf "    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n" ud.du ud.dv ud.dw;
  printf "    reltol = %.1e,  abstol = %.1e\n\n" reltol abstol;

  (* Initialize data structures *)
  let u = Nvector_serial.Any.make n 0.0 in
  let v = Nvector_serial.Any.make n 0.0 in
  let w = Nvector_serial.Any.make n 0.0 in

  (* Create manyvector for solution *)
  let y = Nvector_many.wrap (ROArray.of_list [u; v; w]) in

  let udata = Nvector_serial.Any.unwrap u in
  let vdata = Nvector_serial.Any.unwrap v in
  let wdata = Nvector_serial.Any.unwrap w in

  (* Set initial conditions into y *)
  let pi = 4.0*.atan(one) in
  for i=0 to n-1 do
    let fi = float i in
    udata.{i} <-    ud.a      +. 0.1 *. sin(pi *. fi *. ud.dx); (* u *)
    vdata.{i} <- ud.b /. ud.a +. 0.1 *. sin(pi *. fi *. ud.dx); (* v *)
    wdata.{i} <-    ud.b      +. 0.1 *. sin(pi *. fi *. ud.dx)  (* w *)
  done;

  (* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time t0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  let arkode_mem = ARKStep.(
    init
      (imex ~lsolver:Spils.(solver (spgmr ~maxl:10 y)
                                   ~jac_times_vec:(None, jacvi ud) prec_none)
            ~fi:(fi ud) (fe ud))
      (SStolerances (reltol, abstol))
      t0
      y
  ) in
  (* output spatial mesh to disk *)
  let fid = open_out "bruss_mesh.txt" in
  for i=0 to n-1 do
    fprintf fid "  %.16e\n" (ud.dx*.float i)
  done;
  close_out fid;

  (* Open output streams for results, access data array *)
  let ufid = open_out "bruss_u.txt" in
  let vfid = open_out "bruss_v.txt" in
  let wfid = open_out "bruss_w.txt" in

  (* output initial condition to disk *)
  for i=0 to n-1 do
    fprintf ufid " %.16e" udata.{i};
    fprintf vfid " %.16e" vdata.{i};
    fprintf wfid " %.16e" wdata.{i}
  done;
  fprintf ufid "\n";
  fprintf vfid "\n";
  fprintf wfid "\n";

  (* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let dTout = (tf-.t0)/. float nt in
  let tout = ref (t0+.dTout) in
  printf "        t      ||u||_rms   ||v||_rms   ||w||_rms\n";
  printf "   ----------------------------------------------\n";
  (try
     for iout=0 to nt-1 do
       (* call integrator *)
       let t, _ = ARKStep.evolve_normal arkode_mem !tout y in

       (* access/print solution statistics *)
       let unorm = sqrt (Nvector_serial.DataOps.dotprod udata udata /. float n) in
       let vnorm = sqrt (Nvector_serial.DataOps.dotprod vdata vdata /. float n) in
       let wnorm = sqrt (Nvector_serial.DataOps.dotprod wdata wdata /. float n) in
       printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" t unorm vnorm wnorm;
       (* successful solve: update output time *)
       tout := min (!tout +. dTout) tf;

       (* output results to disk *)
       for i=0 to n-1 do
         fprintf ufid " %.16e" udata.{i};
         fprintf vfid " %.16e" vdata.{i};
         fprintf wfid " %.16e" wdata.{i}
       done;
       fprintf ufid "\n";
       fprintf vfid "\n";
       fprintf wfid "\n"
     done
   with _ ->
     (* unsuccessful solve: break *)
     fprintf stderr "Solver failure, stopping integration\n");
  printf "   ----------------------------------------------\n";
  close_out ufid;
  close_out vfid;
  close_out wfid;

  (* Print some final statistics *)
  let open ARKStep in
  let nst      = get_num_steps arkode_mem in
  let nst_a    = get_num_step_attempts arkode_mem in
  let nfe, nfi = get_num_rhs_evals arkode_mem in
  let nsetups  = get_num_lin_solv_setups arkode_mem in
  let netf     = get_num_err_test_fails arkode_mem in
  let nni      = get_num_nonlin_solv_iters arkode_mem in
  let ncfn     = get_num_nonlin_solv_conv_fails arkode_mem in
  let nli      = Spils.get_num_lin_iters arkode_mem in
  let nlcf     = Spils.get_num_lin_conv_fails arkode_mem in
  let njv      = Spils.get_num_jtimes_evals arkode_mem in
  let nfels    = Spils.get_num_lin_rhs_evals arkode_mem in
  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe nfi;
  printf "   Total linear solver setups = %d\n" nsetups;
  printf "   Total linear iterations = %d\n" nli;
  printf "   Total linear convergence failures = %d\n" nlcf;
  printf "   Total J*v evaluations = %d\n" njv;
  printf "   Total RHS evals in linear solver = %d\n" nfels;
  printf "   Total number of Newton iterations = %d\n" nni;
  printf "   Total number of nonlinear solver convergence failures = %d\n" ncfn;
  printf "   Total number of error test failures = %d\n\n" netf

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

