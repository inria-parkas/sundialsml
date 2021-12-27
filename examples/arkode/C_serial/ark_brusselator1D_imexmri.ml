(* -----------------------------------------------------------------------------
 * Programmer(s): Rujeko Chinomona @SMU and @LLNL
 * -----------------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Dec 2021.
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    u_t = du*u_xx - au*u_x +  a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx - av*v_x +  w*u - v*u^2
 *    w_t = dw*w_xx - aw*w_x + (b-w)/ep - w*u
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
 * We use parameters:
 * du = dv = dw = 0.01 (diffusion coefficients)
 * au = av = aw = -0.001 (advection coefficients - velocity)
 * a  = 0.6
 * b  = 2
 * ep = 0.01
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 * Note: larger values of advection require advection schemes such as
 * upwinding not implemented here.
 *
 * This program solves the problem with multiple solvers listed below.
 * We select method to used based on solve_type input:
 * 0. MIS with third order dirk inner
 * 1. 5th order dirk method for reference solution
 * 2. MRI-GARK34a with erk inner
 * 3. MRI-GARK34a with dirk inner
 * 4. IMEX-MRI3b with erk inner
 * 5. IMEX-MRI3b with dirk inner
 * 6. IMEX-MRI4 with erk inner
 * 7. IMEX-MRI4 with dirk inner
 *
 *  We use Newton iteration with the SUNBAND linear solver and a user supplied
 * Jacobian routine for nonlinear solves.
 *
 * This program solves the problem with the MRI stepper. 10 outputs are printed
 * at equal intervals, and run statistics are printed at the end.
 * ---------------------------------------------------------------------------*)

open Sundials
module ARKStep = Arkode.ARKStep
module MRIStep = Arkode.MRIStep

let printf = Printf.printf
let fprintf = Printf.fprintf
let unwrap = Nvector_serial.unwrap

(* accessor macros between (x,v) location and 1D NVector array *)
let idx x v = 3*x+v

(* user data structure *)
type user_data = {
    n  : int;     (* number of intervals     *)
    dx : float;   (* mesh spacing            *)
    a  : float;   (* constant forcing on u   *)
    b  : float;   (* steady-state value of w *)
    pi : float;   (* value of pi             *)
    du : float;   (* diffusion coeff for u   *)
    dv : float;   (* diffusion coeff for v   *)
    dw : float;   (* diffusion coeff for w   *)
    au : float;   (* advection coeff for u   *)
    av : float;   (* advection coeff for v   *)
    aw : float;   (* advection coeff for w   *)
    ep : float;   (* stiffness parameter     *)
  }

(* f routine to compute the ODE RHS function f(t,y). *)
let f ud _ (y : RealArray.t) (dy : RealArray.t) =
  RealArray.fill dy 0.0; (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  let duconst = ud.du /. ud.dx /. ud.dx in
  let dvconst = ud.dv /. ud.dx /. ud.dx in
  let dwconst = ud.dw /. ud.dx /. ud.dx in
  let auconst = -. ud.au /. 2.0 /. ud.dx in
  let avconst = -. ud.av /. 2.0 /. ud.dx in
  let awconst = -. ud.aw /. 2.0 /. ud.dx in

  for i=1 to ud.n-1-1 do
    (* set shortcuts *)
    let u = y.{idx i 0} and ul = y.{idx (i-1) 0} and ur = y.{idx (i+1) 0} in
    let v = y.{idx i 1} and vl = y.{idx (i-1) 1} and vr = y.{idx (i+1) 1} in
    let w = y.{idx i 2} and wl = y.{idx (i-1) 2} and wr = y.{idx (i+1) 2} in

    (* Fill in ODE RHS for u *)
    dy.{idx i 0} <- (ul -. 2.0*.u +. ur)*.duconst
                        +. (ur -. ul)*.auconst
                        +. ud.a -. (w+.1.0)*.u +. v*.u*.u;

    (* Fill in ODE RHS for v *)
    dy.{idx i 1} <- (vl -. 2.0*.v +. vr)*.dvconst
                        +. (vr -. vl)*.avconst
                        +. w*.u -. v*.u*.u;

    (* Fill in ODE RHS for w *)
    dy.{idx i 2} <- (wl -. 2.0*.w +. wr)*.dwconst
                        +. (wr -. wl)*.awconst
                        +. (ud.b-.w)/.ud.ep -. w*.u
  done;

  (* enforce stationary boundaries *)
  dy.{idx 0 0} <- 0.0;
  dy.{idx 0 1} <- 0.0;
  dy.{idx 0 2} <- 0.0;
  dy.{idx (ud.n-1) 0} <- 0.0;
  dy.{idx (ud.n-1) 1} <- 0.0;
  dy.{idx (ud.n-1) 2} <- 0.0

(* Set the initial condition *)
let set_ic { n; a; b; dx; pi; _ } (y : Nvector_serial.t) =
  let data = Nvector.unwrap y in
  (* Set initial conditions into y *)
  for i = 0 to n - 1 do
    let p = 0.1 *. sin (pi *. float i *. dx) in
    data.{idx i 0} <-   a  +. p;  (* u *)
    data.{idx i 1} <- b/.a +. p;  (* v *)
    data.{idx i 2} <-   b  +. p;  (* w *)
  done

(* Routine to compute the Jacobian matrix from fse(t,y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there *)
let advection_jac { n; dx; au; av; aw; _ } c jac =
  let auconst = -. au /. 2.0 /. dx in
  let avconst = -. av /. 2.0 /. dx in
  let awconst = -. aw /. 2.0 /. dx in

  (* iterate over intervals, filling in Jacobian of (L*y) using SM_ELEMENT_B
     macro (see sunmatrix_band.h) *)
  let incr x y v = Matrix.Band.(set jac x y (get jac x y +. v)) in
  for i = 1 to n - 2 do
    incr (idx i 0) (idx (i - 1) 0) (-.c*.auconst);
    incr (idx i 1) (idx (i - 1) 1) (-.c*.avconst);
    incr (idx i 2) (idx (i - 1) 2) (-.c*.awconst);
    incr (idx i 0) (idx (i + 1) 0) (c*.auconst);
    incr (idx i 1) (idx (i + 1) 1) (c*.avconst);
    incr (idx i 2) (idx (i + 1) 2) (c*.awconst)
  done

(* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there *)
let laplace_matrix ud c jac =
  let inc_elem i j inc = Matrix.Band.update jac i j (fun v -> v +. inc) in

  (* iterate over intervals, filling in Jacobian of (L*y) *)
  for i=1 to ud.n-1-1 do
    inc_elem (idx i 0) (idx (i-1) 0) (c *. ud.du /. ud.dx /. ud.dx);
    inc_elem (idx i 1) (idx (i-1) 1) (c *. ud.dv /. ud.dx /. ud.dx);
    inc_elem (idx i 2) (idx (i-1) 2) (c *. ud.dw /. ud.dx /. ud.dx);
    inc_elem (idx i 0) (idx (i)   0) (-.c *.2.0 *.ud.du /. ud.dx /. ud.dx);
    inc_elem (idx i 1) (idx (i)   1) (-.c *.2.0 *.ud.dv /. ud.dx /. ud.dx);
    inc_elem (idx i 2) (idx (i)   2) (-.c *.2.0 *.ud.dw /. ud.dx /. ud.dx);
    inc_elem (idx i 0) (idx (i+1) 0) (c *. ud.du /. ud.dx /. ud.dx);
    inc_elem (idx i 1) (idx (i+1) 1) (c *. ud.dv /. ud.dx /. ud.dx);
    inc_elem (idx i 2) (idx (i+1) 2) (c *. ud.dw /. ud.dx /. ud.dx)
  done

(* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there *)
let reaction_jac ud c (y : RealArray.t) jac =
  let inc_elem i j inc = Matrix.Band.update jac i j (fun v -> v +. inc) in

  (* iterate over nodes, filling in Jacobian of reaction terms *)
  for i=1 to ud.n-1-1 do

    let u = y.{idx i 0} in   (* set nodal value shortcuts *)
    let v = y.{idx i 1} in
    let w = y.{idx i 2} in

    (* all vars wrt u *)
    inc_elem (idx i 0) (idx i 0) (c*.(2.0*.u*.v-.(w+.1.0)));
    inc_elem (idx i 1) (idx i 0) (c*.(w -. 2.0*.u*.v));
    inc_elem (idx i 2) (idx i 0) (c*.(-.w));

    (* all vars wrt v *)
    inc_elem (idx i 0) (idx i 1) (c*.(u*.u));
    inc_elem (idx i 1) (idx i 1) (c*.(-.u*.u));

    (* all vars wrt w *)
    inc_elem (idx i 0) (idx i 2) (c*.(-.u));
    inc_elem (idx i 1) (idx i 2) (c*.(u));
    inc_elem (idx i 2) (idx i 2) (c*.(-1.0/.ud.ep -. u))
  done

(* ff routine to compute the fast portion of the ODE RHS. *)
let ff { n; a; b; ep; _ } _ (y : RealArray.t) (ydot : RealArray.t) =
  RealArray.fill ydot 0.0;                    (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  for i = 1 to n - 2 do
    (* set shortcuts *)
    let u = y.{idx i 0} in
    let v = y.{idx i 1} in
    let w = y.{idx i 2} in

    (* Fill in ODE RHS for u *)
    ydot.{idx i 0} <- a -. (w+.1.0)*.u +. v*.u*.u;

    (* Fill in ODE RHS for v *)
    ydot.{idx i 1} <- w*.u -. v*.u*.u;

    (* Fill in ODE RHS for w *)
    ydot.{idx i 2} <- (b-.w)/.ep -. w*.u
  done;

  (* enforce stationary boundaries *)
  ydot.{idx 0 0} <- 0.0;
  ydot.{idx 0 1} <- 0.0;
  ydot.{idx 0 2} <- 0.0;
  ydot.{idx (n - 1) 0} <- 0.0;
  ydot.{idx (n - 1) 1} <- 0.0;
  ydot.{idx (n - 1) 2} <- 0.0

(* fse routine to compute the slow-explicit portion of the ODE RHS function. *)
let fse { n; au; av; aw; dx; _ } _ (y : RealArray.t) (ydot : RealArray.t) =
  RealArray.fill ydot 0.0;                    (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  let auconst = -. au /. 2.0 /. dx in
  let avconst = -. av /. 2.0 /. dx in
  let awconst = -. aw /. 2.0 /. dx in
  for i = 1 to n - 2 do
    (* set shortcuts *)
    let ul = y.{idx (i-1) 0} and ur = y.{idx (i+1) 0} in
    let vl = y.{idx (i-1) 1} and vr = y.{idx (i+1) 1} in
    let wl = y.{idx (i-1) 2} and wr = y.{idx (i+1) 2} in

    (* Fill in ODE RHS for u *)
    ydot.{idx i 0} <- (ur -. ul)*.auconst;

    (* Fill in ODE RHS for v *)
    ydot.{idx i 1} <- (vr -. vl)*.avconst;

    (* Fill in ODE RHS for w *)
    ydot.{idx i 2} <- (wr -. wl)*.awconst;
  done;

  (* enforce stationary boundaries *)
  ydot.{idx 0 0} <- 0.0;
  ydot.{idx 0 1} <- 0.0;
  ydot.{idx 0 2} <- 0.0;
  ydot.{idx (n - 1) 0} <- 0.0;
  ydot.{idx (n - 1) 1} <- 0.0;
  ydot.{idx (n - 1) 2} <- 0.0

(* fsi routine to compute the slow-implicit portion of the  ODE RHS. *)
let fsi { n; du; dv; dw; dx; _ } _ (y : RealArray.t) (ydot : RealArray.t) =
  RealArray.fill ydot 0.0;                    (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  let duconst = du/.dx/.dx in
  let dvconst = dv/.dx/.dx in
  let dwconst = dw/.dx/.dx in
  for i = 1 to n - 2 do
    (* set shortcuts *)
    let u = y.{idx i 0} and ul = y.{idx (i - 1) 0} and ur = y.{idx (i + 1) 0} in
    let v = y.{idx i 1} and vl = y.{idx (i - 1) 1} and vr = y.{idx (i + 1) 1} in
    let w = y.{idx i 2} and wl = y.{idx (i - 1) 2} and wr = y.{idx (i + 1) 2} in

    (* Fill in ODE RHS for u *)
    ydot.{idx i 0} <- (ul -. 2.0*.u +. ur)*.duconst;

    (* Fill in ODE RHS for v *)
    ydot.{idx i 1} <- (vl -. 2.0*.v +. vr)*.dvconst;

    (* Fill in ODE RHS for w *)
    ydot.{idx i 2} <- (wl -. 2.0*.w +. wr)*.dwconst
  done;

  (* enforce stationary boundaries *)
  ydot.{idx 0 0} <- 0.0;
  ydot.{idx 0 1} <- 0.0;
  ydot.{idx 0 2} <- 0.0;
  ydot.{idx (n - 1) 0} <- 0.0;
  ydot.{idx (n - 1) 1} <- 0.0;
  ydot.{idx (n - 1) 2} <- 0.0

(* fs routine to compute the slow portion of the ODE RHS. *)
let fs { n; du; dv; dw; au; av;aw; dx; _ } _ (y : RealArray.t) (ydot : RealArray.t) =
  RealArray.fill ydot 0.0;                    (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  let duconst = du/.dx/.dx in
  let dvconst = dv/.dx/.dx in
  let dwconst = dw/.dx/.dx in
  let auconst = -.au/.2.0/.dx in
  let avconst = -.av/.2.0/.dx in
  let awconst = -.aw/.2.0/.dx in
  for i = 1 to n - 2 do
    (* set shortcuts *)
    let u = y.{idx i 0} and ul = y.{idx (i - 1) 0} and ur = y.{idx (i + 1) 0} in
    let v = y.{idx i 1} and vl = y.{idx (i - 1) 1} and vr = y.{idx (i + 1) 1} in
    let w = y.{idx i 2} and wl = y.{idx (i - 1) 2} and wr = y.{idx (i + 1) 2} in

    (* Fill in ODE RHS for u *)
    ydot.{idx i 0} <- (ul -. 2.0*.u +. ur)*.duconst +. (ur -. ul)*.auconst;

    (* Fill in ODE RHS for v *)
    ydot.{idx i 1} <- (vl -. 2.0*.v +. vr)*.dvconst +. (vr -. vl)*.avconst;

    (* Fill in ODE RHS for w *)
    ydot.{idx i 2} <- (wl -. 2.0*.w +. wr)*.dwconst +. (wr -. wl)*.awconst;
  done;

  (* enforce stationary boundaries *)
  ydot.{idx 0 0} <- 0.0;
  ydot.{idx 0 1} <- 0.0;
  ydot.{idx 0 2} <- 0.0;
  ydot.{idx (n - 1) 0} <- 0.0;
  ydot.{idx (n - 1) 1} <- 0.0;
  ydot.{idx (n - 1) 2} <- 0.0

(* Placeholder function of zeroes *)
let f0 _ _ (_ : RealArray.t) (ydot : RealArray.t) =
  RealArray.fill ydot 0.0

(* Jf routine to compute Jacobian of the fast portion of the ODE RHS *)
let jf udata { MRIStep.jac_y = y; _ } j =
  Matrix.Band.set_to_zero j; (* Initialize Jacobian to zero *)
  (* Add in the Jacobian of the reaction terms matrix *)
  reaction_jac udata 1.0 y j

(* Jsi routine to compute the Jacobian of the slow-implicit portion of the ODE RHS. *)
let jsi udata _ j =
  Matrix.Band.set_to_zero j; (* Initialize Jacobian to zero *)
  (* Fill in the Laplace matrix *)
  laplace_matrix udata 1.0 j

(* Js routine to compute the Jacobian of the slow portion of ODE RHS. *)
let js udata _ j =
  Matrix.Band.set_to_zero j; (* Initialize Jacobian to zero *)
  (* Fill in the Laplace matrix *)
  laplace_matrix udata 1.0 j;
  (* Add Jacobian of the advection terms  *)
  advection_jac udata 1.0 j

(* Jac routine to compute the Jacobian of the full ODE RHS. *)
let jac udata { ARKStep.jac_y = y; _ } j =
  Matrix.Band.set_to_zero j; (* Initialize Jacobian to zero *)
  (* Fill in the Laplace matrix *)
  laplace_matrix udata 1.0 j;
  (* Add Jacobian of the advection terms  *)
  advection_jac udata 1.0 j;
  (* Add in the Jacobian of the reaction terms matrix *)
  reaction_jac udata 1.0 y j

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in                   (* initial time                 *)
  let tf = 10.0 in                  (* final time                   *)
  let nt = 10 in                    (* total number of output times *)
  let dtout = (tf-.t0)/.float nt in (* time between outputs         *)
  let nvar = 3 in                   (* number of solution fields    *)
  let m = 10.0 in                   (* time-scale separation factor *)
  let reltol = 1.0e-12 in           (* tolerances                   *)
  let abstol = 1.0e-14 in

  (* Create the SUNDIALS context object for this simulation. *)
  let context = Context.make () in

  (*
   * Initialization
   *)

  (* Retrieve the command-line options: solve_type h  *)
  if Array.length Sys.argv <= 2
  then (printf "ERROR: enter solve_type and hs \n"; exit (-1));
  let solve_type = int_of_string Sys.argv.(1) in
  let hs = float_of_string Sys.argv.(2) in

  (* Check arguments for validity *)
  (*   0 <= solve_type <= 7       *)
  (*   h > 0                      *)

  if solve_type < 0 || solve_type > 7
    then (printf "ERROR: solve_type be an integer in [0,7] \n"; exit (-1));
  let implicit_slow = solve_type > 1 in
  let imex_slow = solve_type > 3 in
  if hs <= 0.0 then (printf "ERROR: hs must be in positive\n"; exit (-1));

  (* store the inputs in the UserData structure *)
  let n = 101 in
  let hf = hs/.m in
  let neq = n * nvar in
  let udata = {
    n = n;                      (* spatial mesh size            *)
    a = 0.6;                    (* problem parameters           *)
    b = 2.0;
    du = 0.01;
    dv = 0.01;
    dw = 0.01;
    au = -0.001;
    av = -0.001;
    aw = -0.001;
    ep = 1.0e-2;                (* stiffness parameter          *)
    pi = 4.0*.atan 1.0;
    dx = 1.0/.float (n - 1);    (* set spatial mesh spacing     *)
  } in

  (* Initial problem output *)
  printf "\n1D Advection-Diffusion-Reaction (Brusselator) test problem:\n";
  printf "    time domain:  (%g,%g]\n" t0 tf;
  printf "    hs = %g\n" hs;
  printf "    hf = %g\n" hf;
  printf "    m  = %g\n" m;
  printf "    N  = %d,  NEQ = %d\n" n neq;
  printf "    dx = %g\n" udata.dx;
  printf "    problem parameters:  a = %g,  b = %g,  ep = %g\n"
         udata.a udata.b udata.ep;
  printf "    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n"
         udata.du udata.dv udata.dw;
  printf "    advection coefficients:  au = %g,  av = %g,  aw = %g\n"
         udata.au udata.av udata.aw;

  let reltol, abstol =
    match solve_type with
    | 0 -> begin
        (* reltol = SUNMAX(hs*hs*hs, 1e-10); *)
        (* abstol = 1e-11; *)
        printf "    solver: exp-3/dirk-3 (MIS / ESDIRK-3-3)\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | 1 -> begin
        let reltol = max (hs*.hs*.hs*.hs*.hs) 1e-14 in
        let abstol = 1e-14 in
        printf "    solver: none/dirk-5 (no slow, default 5th order dirk fast)\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | 2 -> begin
        (* reltol = SUNMAX(hs*hs*hs, 1e-10); *)
        (* abstol = 1e-11; *)
        printf "    solver: dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | 3 -> begin
        (* reltol = SUNMAX(hs*hs*hs, 1e-10); *)
        (* abstol = 1e-11; *)
        printf "    solver: dirk-3/dirk-3 (MRI-GARK-ESDIRK34a / ESDIRK-3-3) -- solve decoupled\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | 4 -> begin
        (* reltol = SUNMAX(hs*hs*hs, 1e-14); *)
        (* abstol = 1e-14; *)
        printf "    solver: ars343/exp-3 (IMEX-MRI3b / ERK-3-3) -- solve decoupled\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | 5 -> begin
        (* reltol = SUNMAX(hs*hs*hs, 1e-14); *)
        (* abstol = 1e-14; *)
        printf "    solver: ars343/dirk-3 (IMEX-MRI3b / ESDIRK-3-3) -- solve decoupled\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | 6 -> begin
        (* reltol = SUNMAX(hs*hs*hs*hs, 1e-14); *)
        (* abstol = 1e-14; *)
        printf "    solver: imexark4/exp-4 (IMEX-MRI4 / ERK-4-4) -- solve decoupled\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | 7 -> begin
        (* reltol = SUNMAX(hs*hs*hs*hs, 1e-14); *)
        (* abstol = 1e-14; *)
        printf "    solver: imexark4/dirk-4 (IMEX-MRI4 / CASH(5,3,4)-DIRK ) -- solve decoupled\n\n";
        printf "    reltol = %.2e,  abstol = %.2e\n\n" reltol abstol;
        reltol, abstol
      end
    | _ -> assert false
  in

  (* Create solution vector *)
  let y = Nvector_serial.make ~context neq 0.0 in (* Create vector for solution *)

  (* Set initial condition *)
  set_ic udata y;

  (* Create vector masks  *)
  let umask = Nvector.clone y in
  let vmask = Nvector.clone y in
  let wmask = Nvector.clone y in

  (* Set mask array values for each solution component *)
  Nvector_serial.Ops.const 0.0 umask;
  let data = Nvector.unwrap umask in
  for i = 0 to n - 1 do data.{idx i 0} <- 1.0 done;

  Nvector_serial.Ops.const 0.0 vmask;
  let data = Nvector.unwrap vmask in
  for i = 0 to n - 1 do data.{idx i 1} <- 1.0 done;

  Nvector_serial.Ops.const 0.0 wmask;
  let data = Nvector.unwrap wmask in
  for i = 0 to n - 1 do data.{idx i 2} <- 1.0 done;

  (*
   * Create the fast integrator and set options
   *)

  (* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the inital time T0,
     and the initial dependent variable vector y. *)
  let inner_arkode_mem =
    match solve_type with
    | 0 | 3 | 5 -> (* esdirk-3-3 fast solver *)
      let af = Matrix.band ~context ~mu:4 ~ml:4 neq in
      let lsf = LinearSolver.Direct.band ~context y af in
      let inner_arkode_mem =
        ARKStep.(init ~context
          (implicit ~lsolver:(Dls.solver ~jac:(jf udata) lsf) (ff udata))
          (SStolerances (reltol, abstol)) t0 y)
      in
      ARKStep.set_max_nonlin_iters inner_arkode_mem 10;
      let beta  = (sqrt 3.0) /. 6.0 +. 0.5 in
      let gamma = (-1.0 /. 8.0) *. ((sqrt 3.0) +. 1.0) in
      let implicit_table = Arkode.ButcherTable.{
        method_order = 3;
        stages = 3;
        stage_values = RealArray2.of_lists [
          [ 0.0; 4.0*.gamma +. 2.0*.beta; 0.5 -. beta -. gamma ];
          [ 0.0; 1.0 -. 4.0*.gamma -. 2.0*.beta; gamma ];
          [ 0.0; 0.0; beta ]
        ];
        stage_times = RealArray.of_list [ 0.0; 1.0; 0.5 ];
        coefficients = RealArray.of_list
                         [ 1.0 /. 6.0; 1.0 /. 6.0; 2.0 /. 3.0 ];
        embedding = None;
      } in
      ARKStep.set_tables inner_arkode_mem ~implicit_table ();
      inner_arkode_mem

    | 1 -> (* dirk 5th order fast solver (full problem) *)
      let af = Matrix.band ~context ~mu:4 ~ml:4 neq in
      let lsf = LinearSolver.Direct.band ~context y af in
      ARKStep.(init ~context
        (implicit ~lsolver:(Dls.solver ~jac:(jac udata) lsf) (f udata))
        (SStolerances (reltol, abstol)) ~order:5 t0 y)

    | 2 | 4 -> (* erk-3-3 fast solver *)
      let inner_arkode_mem =
        ARKStep.(init ~context (explicit (ff udata)) default_tolerances t0 y)
      in
      let explicit_table = Arkode.ButcherTable.{
        method_order = 3;
        stages = 3;
        stage_values = RealArray2.of_lists [
          [ 0.0; 0.5; -1.0 ];
          [ 0.0; 0.0;  2.0 ];
          [ 0.0; 0.0;  0.0 ]
        ];
        stage_times = RealArray.of_list [ 0.0; 0.5; 1.0 ];
        coefficients = RealArray.of_list
                         [ 1.0 /. 6.0; 2.0 /. 3.0; 1.0 /. 6.0 ];
        embedding = Some (2, RealArray.of_list [ 0.0; 1.0; 0.0 ] );
      } in
      ARKStep.set_tables inner_arkode_mem ~explicit_table ();
      inner_arkode_mem

    | 6 -> (* erk-4-4 fast solver *)
      let inner_arkode_mem =
        ARKStep.(init ~context (explicit (ff udata)) default_tolerances t0 y)
      in
      let explicit_table = Arkode.ButcherTable.{
        method_order = 4;
        stages = 4;
        stage_values = RealArray2.of_lists [
          [ 0.0; 0.5; 0.0; 0.0 ];
          [ 0.0; 0.0; 0.5; 0.0 ];
          [ 0.0; 0.0; 0.0; 1.0 ];
          [ 0.0; 0.0; 0.0; 0.0 ]
        ];
        stage_times = RealArray.of_list [ 0.0; 0.5; 0.5; 1.0 ];
        coefficients = RealArray.of_list
                         [ 1.0 /. 6.0; 1.0 /. 3.0; 1.0 /. 3.0; 1.0 /. 6.0 ];
        embedding = None;
      } in
      ARKStep.set_tables inner_arkode_mem ~explicit_table ();
      inner_arkode_mem

    | 7 -> (* Cash(5,3,4)-SDIRK fast solver *)
      let af = Matrix.band ~context ~mu:4 ~ml:4 neq in
      let lsf = LinearSolver.Direct.band ~context y af in
      let inner_arkode_mem =
        ARKStep.(init ~context
          (implicit ~lsolver:(Dls.solver ~jac:(jf udata) lsf) (ff udata))
          (SStolerances (reltol, abstol)) t0 y)
      in
      ARKStep.set_max_nonlin_iters inner_arkode_mem 10;
      ARKStep.set_dirk_table_num inner_arkode_mem Arkode.ButcherTable.Cash_5_3_4;
      inner_arkode_mem

    | _ -> assert false
  in
  (* Set the fast step size *)
  ARKStep.set_fixed_step inner_arkode_mem (Some hf);
  (* Create inner stepper *)
  let inner_stepper = MRIStep.InnerStepper.from_arkstep inner_arkode_mem in

  (*
   * Create the slow integrator and set options
   *)

  (* Initialize the slow integrator. Specify the slow right-hand side
     function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the inital time
     T0, the initial dependent variable vector y, and the fast integrator. *)
  let arkode_mem = match solve_type with
    | 0 -> (* use MIS outer integrator default for MRIStep *)
      MRIStep.(init ~context (explicit (fs udata))
                    default_tolerances inner_stepper ~slowstep:hs t0 y)

    | 1 -> (* no slow dynamics (use ERK-2-2) *)
      let b = Arkode.ButcherTable.{
        method_order = 2;
        stages = 2;
        stage_values = RealArray2.of_lists [
          [ 0.0; 2.0/.3.0 ];
          [ 0.0;   0.0    ];
        ];
        stage_times = RealArray.of_list [ 0.0; 2.0 /. 3.0 ];
        coefficients = RealArray.of_list [ 0.25; 0.75 ];
        embedding = None;
      } in
      let coupling =
        MRIStep.Coupling.mis_to_mri ~method_order:2 ~embedding_order:0 b
      in
      MRIStep.(init ~context (explicit (f0 udata))
                    default_tolerances inner_stepper ~coupling ~slowstep:hs t0 y)

    | 2 | 3 -> (* MRI-GARK-ESDIRK34a, solve-decoupled slow solver *)
      let a_s = Matrix.band ~context ~mu:4 ~ml:4 neq in
      let ls_s = LinearSolver.Direct.band ~context y a_s in
      let coupling = MRIStep.Coupling.(load_table GARK_ESDIRK34a) in
      MRIStep.(init ~context
        (implicit ~lsolver:(Dls.solver ~jac:(js udata) ls_s) (fs udata))
        (SStolerances (reltol, abstol))
        inner_stepper ~coupling ~slowstep:hs t0 y)

    | 4 | 5 -> (* IMEX-MRI-GARK3b, solve-decoupled slow solver *)
      let a_s = Matrix.band ~context ~mu:4 ~ml:4 neq in
      let ls_s = LinearSolver.Direct.band ~context y a_s in
      let coupling = MRIStep.Coupling.(load_table IMEX_GARK3b) in
      MRIStep.(init ~context
                    (imex ~lsolver:(Dls.solver ~jac:(jsi udata) ls_s)
                           ~fse:(fse udata) ~fsi:(fsi udata) ())
                    (SStolerances (reltol, abstol))
                    inner_stepper ~coupling ~slowstep:hs t0 y)

    | 6 | 7 -> (* IMEX-MRI-GARK4, solve-decoupled slow solver *)
      let a_s = Matrix.band ~context ~mu:4 ~ml:4 neq in
      let ls_s = LinearSolver.Direct.band ~context y a_s in
      let coupling = MRIStep.Coupling.(load_table IMEX_GARK4) in
      MRIStep.(init ~context
                    (imex ~lsolver:(Dls.solver ~jac:(jsi udata) ls_s)
                           ~fse:(fse udata) ~fsi:(fsi udata) ())
                    (SStolerances (reltol, abstol))
                    inner_stepper ~coupling ~slowstep:hs t0 y)

    | _ -> assert false
  in
  (* Set maximum number of steps taken by solver *)
  MRIStep.set_max_num_steps arkode_mem 1000000;

  (*
   * Integrate ODE
   *)

  (* output spatial mesh to disk *)
  let fid = open_out "bruss1D_mesh.txt" in
  for i = 0 to n - 1 do
    fprintf fid "  %.16e\n" (udata.dx *. float i)
  done;
  close_out fid;

  (* Open output stream for results, access data arrays *)
  let ufid = open_out (Printf.sprintf "bruss1D_u_%s_%s.txt"
                                      Sys.argv.(1) Sys.argv.(2))
  in
  let vfid = open_out (Printf.sprintf "bruss1D_v_%s_%s.txt"
                                      Sys.argv.(1) Sys.argv.(2))
  in
  let wfid = open_out (Printf.sprintf "bruss1D_w_%s_%s.txt"
                                      Sys.argv.(1) Sys.argv.(2))
  in

  (* output initial condition to disk *)
  let data = Nvector.unwrap y in

  for i =0 to n - 1 do
    fprintf ufid " %.16e" data.{idx i 0}
  done;
  for i =0 to n - 1 do
    fprintf vfid " %.16e" data.{idx i 1}
  done;
  for i =0 to n - 1 do
    fprintf wfid " %.16e" data.{idx i 2}
  done;
  fprintf ufid "\n";
  fprintf vfid "\n";
  fprintf wfid "\n";

  (* Main time-stepping loop: calls MRIStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  printf "        t      ||u||_rms   ||v||_rms   ||w||_rms\n";
  printf "   ----------------------------------------------\n";
  let n_f = float n in
  let rec loop iout tout =
    if iout = nt then ()
    else begin
      (* call integrator *)
      let t, _ = MRIStep.evolve_normal arkode_mem tout y in

      (* access/print solution statistics *)
      let u = Nvector_serial.Ops.wl2norm y umask in
      let u = sqrt (u *. u /. n_f) in
      let v = Nvector_serial.Ops.wl2norm y vmask in
      let v = sqrt (v *. v /. n_f) in
      let w = Nvector_serial.Ops.wl2norm y wmask in
      let w = sqrt (w *. w /. n_f) in
      printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" t u v w;

      (* output results to disk *)
      for i = 0 to n - 1 do
        fprintf ufid " %.16e" data.{idx i 0}
      done;
      for i = 0 to n - 1 do
        fprintf vfid " %.16e" data.{idx i 1}
      done;
      for i = 0 to n - 1 do
        fprintf wfid " %.16e" data.{idx i 2}
      done;
      fprintf ufid "\n";
      fprintf vfid "\n";
      fprintf wfid "\n";

      loop (iout + 1) (min (tout +. dtout) tf)
    end
  in
  loop 0 (t0 +. dtout);
  printf "   ----------------------------------------------\n";
  close_out ufid;
  close_out vfid;
  close_out wfid;

  (*
   * Finalize
   *)

  (* Get some slow integrator statistics *)
  let nsts = MRIStep.get_num_steps arkode_mem in
  let nfse, nfsi = MRIStep.get_num_rhs_evals arkode_mem in

  (* Get some fast integrator statistics *)
  let nstf = ARKStep.get_num_steps inner_arkode_mem in
  let nffe, nffi = ARKStep.get_num_rhs_evals inner_arkode_mem in

  (* Print some final statistics *)
  printf "\nFinal Solver Statistics:\n";
  printf "   Slow Steps: nsts = %d\n" nsts;
  printf "   Fast Steps: nstf = %d\n" nstf;
  if imex_slow then begin
    if solve_type=0 || solve_type=1 || solve_type=3 || solve_type=5 || solve_type=7
    then printf "   Total RHS evals:  Fse = %d, Fsi = %d,  Ff = %d\n" nfse nfsi nffi
    else printf "   Total RHS evals:  Fse = %d, Fsi = %d,  Ff = %d\n" nfse nfsi nffe
  end
  else if implicit_slow then begin
    if solve_type=0 || solve_type=1 || solve_type=3 || solve_type=5 || solve_type=7
    then printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfsi nffi
    else printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfsi nffe
  end
  else begin
    if solve_type=0 || solve_type=1 || solve_type=3 || solve_type=5 || solve_type=7
    then printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfse nffi
    else printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfse nffe
  end;

  (* Get/print slow integrator decoupled implicit solver statistics *)
  if solve_type > 1 then begin
    let nnis, nncs = MRIStep.get_nonlin_solv_stats arkode_mem in
    let njes = MRIStep.Dls.get_num_jac_evals arkode_mem in
    printf "   Slow Newton iters = %d\n" nnis;
    printf "   Slow Newton conv fails = %d\n" nncs;
    printf "   Slow Jacobian evals = %d\n" njes
  end;

  (* Get/print fast integrator implicit solver statistics *)
  if solve_type=0 || solve_type=1 || solve_type=3 || solve_type=5 || solve_type=7 then begin
    let nnif, nncf = ARKStep.get_nonlin_solv_stats inner_arkode_mem in
    let njef = ARKStep.Dls.get_num_jac_evals inner_arkode_mem in
    printf "   Fast Newton iters = %d\n" nnif;
    printf "   Fast Newton conv fails = %d\n" nncf;
    printf "   Fast Jacobian evals = %d\n" njef
  end

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

