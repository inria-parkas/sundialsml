(* -------------------------------------------------------------- {{{
 * Programmer(s): David J. Gardner @ LLNL
 * --------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Oct 2021.
 * --------------------------------------------------------------
 * Based on the ark_brusselator1D_omp.c ARKode example by
 * Daniel R. Reynolds @ SMU.
 * --------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * --------------------------------------------------------------
 * This program simulates 1D advection-reaction problem. The
 * brusselator problem from chemical kinetics is used for the
 * reaction terms. This is a PDE system with 3 components,
 * Y = [u,v,w], satisfying the equations,
 *
 *    u_t = -c*u_x + a - (w+1)*u + v*u^2
 *    v_t = -c*v_x + w*u - v*u^2
 *    w_t = -c*w_x + (b-w)/ep - w*u
 *
 * for t in [0, 10], x in [0, 1], with initial conditions
 *
 *    u(0,x) =  a  + 0.1*exp(-(x-0.5)^2 / 0.1)
 *    v(0,x) = b/a + 0.1*exp(-(x-0.5)^2 / 0.1)
 *    w(0,x) =  b  + 0.1*exp(-(x-0.5)^2 / 0.1),
 *
 * and with periodic boundary conditions i.e.,
 *
 *    u(t,0) = u(t,1),
 *    v(t,0) = v(t,1),
 *    w(t,0) = w(t,1).
 *
 * The spatial derivatives are computed using first-order
 * upwind differences for the advection terms. The data is
 * distributed over N points on a uniform spatial grid.
 *
 * This program use the MRIStep module with an explicit slow
 * method and an implicit fast method. The explicit method
 * uses a fixed step size and the implicit method uses adaptive
 * steps. Implicit systems are solved using a Newton iteration
 * with the band linear solver, and a user-supplied Jacobian
 * routine for the fast RHS.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 * -------------------------------------------------------------- }}} *)

open Sundials
module ARKStep = Arkode.ARKStep
module MRIStep = Arkode.MRIStep

let idx x v = 3 * x + v
let printf = Printf.printf
let fprintf = Printf.fprintf
let sqr x = x *. x

(* user data structure *)
type userdata = {
  n  : int;    (* number of intervals      *)
  dx : float;  (* mesh spacing             *)
  a  : float;  (* constant forcing on u    *)
  b  : float;  (* steady-state value of w  *)
  c  : float;  (* advection coefficient    *)
  ep : float;  (* stiffness parameter      *)
}

(* -----------------------------------
 * Functions called by the integrator
 * -----------------------------------*)

(* ff routine to compute the fast portion of the ODE RHS. *)
let ff { n; a; b; ep; _ } t (ydata : RealArray.t) (dydata : RealArray.t) =
  (* iterate over domain, computing reactions *)
  for i = 0 to n - 1 do

    (* set shortcuts *)
    let u = ydata.{idx i 0} in
    let v = ydata.{idx i 1} in
    let w = ydata.{idx i 2} in

    (* u_t = a - (w+1)*u + v*u^2 *)
    dydata.{idx i 0} <- a -. (w +. 1.0) *. u +. v *. u *. u;

    (* v_t = w*u - v*u^2 *)
    dydata.{idx i 1} <- w *. u -. v *. u *. u;

    (* w_t = (b-w)/ep - w*u *)
    dydata.{idx i 2} <- (b -. w) /. ep -. w *. u
  done

(* fs routine to compute the slow portion of the ODE RHS. *)
let fs { n; c; dx; _ } t (ydata : RealArray.t) (dydata : RealArray.t) =
  (* iterate over domain, computing advection *)
  let tmp = -. c /. dx in

  if c > 0.0 then begin
    (* right moving flow *)

    (* left boundary Jacobian entries *)
    dydata.{idx 0 0} <- tmp *. (ydata.{idx 0 0} -. ydata.{idx (n-1) 0});
    dydata.{idx 0 1} <- tmp *. (ydata.{idx 0 1} -. ydata.{idx (n-1) 1});
    dydata.{idx 0 2} <- tmp *. (ydata.{idx 0 2} -. ydata.{idx (n-1) 2});

    (* interior Jacobian entries *)
    for i = 1 to n-1 do
      dydata.{idx i 0} <- tmp *. (ydata.{idx i 0} -. ydata.{idx (i-1) 0});
      dydata.{idx i 1} <- tmp *. (ydata.{idx i 1} -. ydata.{idx (i-1) 1});
      dydata.{idx i 2} <- tmp *. (ydata.{idx i 2} -. ydata.{idx (i-1) 2})
    done

  end else if c < 0.0 then begin
    (* left moving flow *)

    (* interior Jacobian entries *)
    for i = 0 to n-2 do
      dydata.{idx i 0} <- tmp *. (ydata.{idx (i+1) 0} -. ydata.{idx i 0});
      dydata.{idx i 1} <- tmp *. (ydata.{idx (i+1) 1} -. ydata.{idx i 1});
      dydata.{idx i 2} <- tmp *. (ydata.{idx (i+1) 2} -. ydata.{idx i 2})
    done;

    (* right boundary Jacobian entries *)
    dydata.{idx (n-1) 0} <- tmp *. (ydata.{idx (n-1) 0} -. ydata.{idx 0 0});
    dydata.{idx (n-1) 1} <- tmp *. (ydata.{idx (n-1) 1} -. ydata.{idx 0 1});
    dydata.{idx (n-1) 2} <- tmp *. (ydata.{idx (n-1) 2} -. ydata.{idx 0 2})
  end


(* Js routine to compute the Jacobian of the fast portion of the ODE RHS. *)

let jf {n; ep; _} { MRIStep.jac_y = (ydata : RealArray.t) } jac =
  (* iterate over nodes, filling in Jacobian entries *)
  for i = 0 to n -1 do
    (* set nodal value shortcuts (shifted index due to start at first interior node) *)
    let u = ydata.{idx i 0} in
    let v = ydata.{idx i 1} in
    let w = ydata.{idx i 2} in

    (* all vars wrt u *)
    let open Matrix.Band in
    set jac (idx i 0) (idx i 0) (2.0 *. u *. v -. (w +. 1.0));
    set jac (idx i 1) (idx i 0) (w -. 2.0 *. u *. v);
    set jac (idx i 2) (idx i 0) (-. w);

    (* all vars wrt v *)
    set jac (idx i 0) (idx i 1) (u *. u);
    set jac (idx i 1) (idx i 1) (-. u *. u);

    (* all vars wrt w *)
    set jac (idx i 0) (idx i 2) (-. u);
    set jac (idx i 1) (idx i 2) (u);
    set jac (idx i 2) (idx i 2) (-1.0 /. ep -. u)
  done

(* ------------------------------
 * Private helper functions
 * ------------------------------*)

(* Set the initial condition *)
let setic { n; a; b; dx; _ } (y : Nvector_serial.t) =
  let data = Nvector.unwrap y in
  (* Set initial conditions into y *)
  for i = 0 to n - 1 do
    let x = float i *. dx in
    let p = 0.1 *. exp( -.(sqr (x -. 0.5)) /. 0.1) in
    data.{idx i 0} <-   a  +. p;
    data.{idx i 1} <- b/.a +. p;
    data.{idx i 2} <-   b  +. p
  done

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0     = 0.0 in      (* initial time                    *)
  let tf     = 10.0 in     (* final time                      *)
  let nt     = 100 in      (* total number of output times    *)
  let nvar   = 3 in        (* number of solution fields       *)
  let n      = 200 in      (* spatial mesh size (N intervals) *)
  let a      = 1.0 in      (* problem parameters              *)
  let b      = 3.5 in
  let c      = 0.25 in
  let ep     = 1.0e-6 in   (* stiffness parameter *)
  let reltol = 1.0e-6 in   (* tolerances          *)
  let abstol = 1.0e-10 in

  (* allocate udata structure *)
  let udata = {
    n; a; b; c; ep;
    dx = 1.0 /. float n   (* periodic BC, divide by N not N-1 *)
  } in

  (* set total allocated vector length *)
  let neq = nvar * udata.n in

  (* set the slow step size *)
  let hs = 0.5 *. (udata.dx /. abs_float c) in

  (* Initial problem output *)
  printf "\n1D Advection-Reaction example problem:\n";
  printf "    N = %d,  NEQ = %d\n" udata.n neq;
  printf "    problem parameters:  a = %g,  b = %g,  ep = %g\n"
         udata.a udata.b udata.ep;
  printf "    advection coefficient:  c = %g\n" udata.c;
  printf "    reltol = %.1e,  abstol = %.1e\n\n" reltol abstol;

  (* Create solution vector *)
  let y = Nvector_serial.make neq 0.0 in (* Create vector for solution *)

  (* Set initial condition *)
  setic udata y;

  (* Create vector masks *)
  let umask = Nvector_serial.make neq 0.0 in
  let vmask = Nvector_serial.make neq 0.0 in
  let wmask = Nvector_serial.make neq 0.0 in

  (* Set mask array values for each solution component *)
  let data = Nvector.unwrap umask in
  for i = 0 to n - 1 do
    data.{idx i 0} <- 1.0
  done;

  let data = Nvector.unwrap vmask in
  for i = 0 to n - 1 do
    data.{idx i 1} <- 1.0
  done;

  let data = Nvector.unwrap wmask in
  for i = 0 to n - 1 do
    data.{idx i 2} <- 1.0
  done;

  (*
   * Create the fast integrator and set options
   *)

  (* Initialize matrix and linear solver data structures *)
  let a = Matrix.band ~mu:4 ~ml:4 neq in
  let ls = LinearSolver.Direct.band y a in

  (* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, and the
     initial dependent variable vector y. *)
  (* Attach matrix and linear solver *)
  (* Set the Jacobian routine *)
  (* Specify fast tolerances *)
  (* Attach user data to fast integrator *)
  let inner_arkode_mem =
    ARKStep.(init (implicit ~lsolver:(Dls.solver ~jac:(jf udata) ls) (ff udata))
                  (SStolerances (reltol, abstol))
                  t0 y)
  in
  (* Set the fast method *)
  ARKStep.set_dirk_table_num inner_arkode_mem
             Arkode.ButcherTable.ARK324L2SA_DIRK_4_2_3;

  (*
   * Create the slow integrator and set options
   *)

  (* Initialize the slow integrator. Specify the slow right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, the
     initial dependent variable vector y, and the fast integrator. *)
  (* Pass udata to user functions *)
  (* Set the slow step size *)
  let arkode_mem = MRIStep.(init
                              InnerStepper.(from_arkstep inner_arkode_mem)
                               default_tolerances
                               (fs udata)
                               ~slowstep:hs
                               t0 y)
  in

  (* output spatial mesh to disk (add extra point for periodic BC) *)
  let fid = open_out "mesh.txt" in
  for i = 0 to n do
    fprintf fid "  %.16e\n" (udata.dx *. (float i));
  done;
  close_out fid;

  (* Open output stream for results, access data arrays *)
  let ufid = open_out "u.txt" in
  let vfid = open_out "v.txt" in
  let wfid = open_out "w.txt" in

  (* output initial condition to disk (extra output for periodic BC) *)
  let data = Nvector.unwrap y in

  for i = 0 to n - 1 do
    fprintf ufid " %.16e" data.{idx i 0}
  done;
  fprintf ufid " %.16e" data.{idx 0 0};
  fprintf ufid "\n";

  for i = 0 to n - 1 do
    fprintf vfid " %.16e" data.{idx i 1}
  done;
  fprintf vfid " %.16e" data.{idx 0 1};
  fprintf vfid "\n";

  for i = 0 to n - 1 do
    fprintf wfid " %.16e" data.{idx i 2}
  done;
  fprintf wfid " %.16e" data.{idx 0 2};
  fprintf wfid "\n";

  (* Main time-stepping loop: calls ARKStepEvolve to perform the integration,
     then prints results.  Stops when the final time has been reached *)
  let dTout = (tf -. t0) /. float nt in
  let tout = ref (t0 +. dTout) in
  printf "        t      ||u||_rms   ||v||_rms   ||w||_rms\n";
  printf "   ----------------------------------------------\n";
  for iout = 0 to nt - 1 do
    (* call integrator *)
    let t, r = MRIStep.evolve_normal arkode_mem !tout y in

    (* access/print solution statistics *)
    let u = Nvector_serial.Ops.wl2norm y umask in
    let u = sqrt (u *. u /. float n) in
    let v = Nvector_serial.Ops.wl2norm y vmask in
    let v = sqrt (v *. v /. float n) in
    let w = Nvector_serial.Ops.wl2norm y wmask in
    let w = sqrt (w *. w /. float n) in
    printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" t u v w;

    (* output results to disk (extr output for periodic BC) *)
    for i = 0 to n - 1 do
      fprintf ufid " %.16e" data.{idx i 0}
    done;
    fprintf ufid " %.16e" data.{idx 0 0};
    fprintf ufid "\n";

    for i = 0 to n - 1 do
      fprintf vfid " %.16e" data.{idx i 1}
    done;
    fprintf vfid " %.16e" data.{idx 0 1};
    fprintf vfid "\n";

    for i = 0 to n - 1 do
      fprintf wfid " %.16e" data.{idx i 2}
    done;
    fprintf wfid " %.16e" data.{idx 0 2};
    fprintf wfid "\n";

    (* successful solve: update output time *)
    tout := min tf (!tout +. dTout)
  done;
  printf "   ----------------------------------------------\n";
  close_out ufid;
  close_out vfid;
  close_out wfid;

  (* Get some slow integrator statistics *)
  let open MRIStep in
  let nsts  = get_num_steps arkode_mem in
  let nfs   = get_num_rhs_evals arkode_mem in

  (* Get some fast integrator statistics *)
  let open ARKStep in
  let nstf       = get_num_steps inner_arkode_mem in
  let nstf_a     = get_num_step_attempts inner_arkode_mem in
  let nffe, nffi = get_num_rhs_evals inner_arkode_mem in
  let nsetups    = get_num_lin_solv_setups inner_arkode_mem in
  let netf       = get_num_err_test_fails inner_arkode_mem in
  let nni        = get_num_nonlin_solv_iters inner_arkode_mem in
  let ncfn       = get_num_nonlin_solv_conv_fails inner_arkode_mem in
  let nje        = Dls.get_num_jac_evals inner_arkode_mem in
  let nfeLS      = Dls.get_num_lin_rhs_evals inner_arkode_mem in

  (* Print some final statistics *)
  printf "\nFinal Solver Statistics:\n";
  printf "   Slow Steps: nsts = %d\n" nsts;
  printf "   Fast Steps: nstf = %d (attempted = %d)\n" nstf nstf_a;
  printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfs nffi;
  printf "   Total number of fast error test failures = %d\n" netf;
  printf "   Total linear solver setups = %d\n" nsetups;
  printf "   Total RHS evals for setting up the linear system = %d\n" nfeLS;
  printf "   Total number of Jacobian evaluations = %d\n" nje;
  printf "   Total number of Newton iterations = %d\n" nni;
  printf "   Total number of nonlinear solver convergence failures = %d\n" ncfn

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
