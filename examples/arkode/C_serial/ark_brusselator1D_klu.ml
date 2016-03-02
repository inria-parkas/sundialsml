(*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jan 2016.
 *---------------------------------------------------------------
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 *---------------------------------------------------------------
 * Example problem:
 * 
 * The following test simulates a brusselator problem from chemical 
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w], 
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions
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
 * The number of spatial points N, the parameters a, b, du, dv, 
 * dw and ep, as well as the desired relative and absolute solver 
 * tolerances, are provided in the input file 
 * input_brusselator1D.txt.
 * 
 * This program solves the problem with the DIRK method, using a
 * Newton iteration.  The inner linear systems are solved using 
 * the ARKKLU linear solver.
 *
 * 100 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 *---------------------------------------------------------------*)

module RealArray = Sundials.RealArray
let printf = Printf.printf
let fprintf = Printf.fprintf
let unwrap = Nvector_serial.unwrap
let n_vwl2norm = Nvector_serial.Ops.n_vwl2norm

(* accessor macros between (x,v) location and 1D NVector array *)
let idx x v = 3*x+v

(* user data structure *)
type user_data = {
    n  : int;        (* number of intervals     *)
    dx : float;      (* mesh spacing            *)
    a  : float;      (* constant forcing on u   *)
    b  : float;      (* steady-state value of w *)
    du : float;      (* diffusion coeff for u   *)
    dv : float;      (* diffusion coeff for v   *)
    dw : float;      (* diffusion coeff for w   *)
    ep : float;      (* stiffness parameter     *)
    mutable r  : Sls.SparseMatrix.t option
                     (* temporary storage       *)
  }

(* f routine to compute the ODE RHS function f(t,y). *)
let f ud t y dy =
  RealArray.fill dy 0.0; (* initialize ydot to zero *)

  (* iterate over domain, computing all equations *)
  let uconst = ud.du /. ud.dx /. ud.dx in
  let vconst = ud.dv /. ud.dx /. ud.dx in
  let wconst = ud.dw /. ud.dx /. ud.dx in

  for i=1 to ud.n-1-1 do
    (* set shortcuts *)
    let u = y.{idx i 0} and ul = y.{idx (i-1) 0} and ur = y.{idx (i+1) 0} in
    let v = y.{idx i 1} and vl = y.{idx (i-1) 1} and vr = y.{idx (i+1) 1} in
    let w = y.{idx i 2} and wl = y.{idx (i-1) 2} and wr = y.{idx (i+1) 2} in

    (* u_t = du*u_xx + a - (w+1)*u + v*u^2 *)
    dy.{idx i 0} <- (ul -. 2.0*.u +. ur)*.uconst
                        +. ud.a -. (w+.1.0)*.u +. v*.u*.u;

    (* v_t = dv*v_xx + w*u - v*u^2 *)
    dy.{idx i 1} <- (vl -. 2.0*.v +. vr)*.vconst +. w*.u -. v*.u*.u;

    (* w_t = dw*w_xx + (b-w)/ep - w*u *)
    dy.{idx i 2} <- (wl -. 2.0*.w +. wr)*.wconst +. (ud.b-.w)/.ud.ep -. w*.u
  done;

  (* enforce stationary boundaries *)
  dy.{idx 0 0} <- 0.0;
  dy.{idx 0 1} <- 0.0;
  dy.{idx 0 2} <- 0.0;
  dy.{idx (ud.n-1) 0} <- 0.0;
  dy.{idx (ud.n-1) 1} <- 0.0;
  dy.{idx (ud.n-1) 2} <- 0.0

(* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there *)
let laplace_matrix ud jac =
  let nz = ref 0 in
  let set_col j = Sls.SparseMatrix.set_col jac j !nz in
  let set j v = (Sls.SparseMatrix.set jac !nz j v; incr nz) in

  (* set first column to zero *)
  set_col (idx 0 0);
  set_col (idx 0 1);
  set_col (idx 0 2);
  
  (* iterate over nodes, filling in Laplacian entries depending on these *)
  let uconst  = ud.du /. ud.dx /. ud.dx in
  let uconst2 = -2.0 *. uconst in
  let vconst  = ud.dv /. ud.dx /. ud.dx in
  let vconst2 = -2.0 *. vconst in
  let wconst  = ud.dw /. ud.dx /. ud.dx in
  let wconst2 = -2.0 *. wconst in

  for i=1 to ud.n-1-1 do
    (* dependence on u at this node *)
    set_col (idx i 0);
    if i>1 then set (idx (i-1) 0) uconst;      (* node to left *)
    set (idx (i) (0)) uconst2;                 (* self *)
    if i<ud.n-2 then set (idx (i+1) 0) uconst; (* node to right *)

    (* dependence on v at this node *)
    set_col (idx i 1);
    if i>1 then set (idx (i-1) 1) vconst;      (* node to left *)
    set (idx i 1) vconst2;                     (* self *)
    if i<ud.n-2 then set (idx (i+1) 1) vconst; (* node to right *)

    (* dependence on w at this node *)
    set_col (idx i 2);
    if i>1 then set (idx (i-1) 2) wconst;      (* node to left *)
    set (idx i 2) wconst2;                     (* self *)
    if i<ud.n-2 then set (idx (i+1) 2) wconst  (* node to right *)
  done;

  (* set last column to zero *)
  set_col (idx (ud.n-1) 0);
  set_col (idx (ud.n-1) 1);
  set_col (idx (ud.n-1) 2);
  
  (* end of data *)
  set_col ((idx (ud.n-1) 2)+1)

(* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there *)
let reaction_jac ud y jac =
  let nz = ref 0 in
  let set_col j = Sls.SparseMatrix.set_col jac j !nz in
  let set j v = (Sls.SparseMatrix.set jac !nz j v; incr nz) in

  (* set first matrix column to zero *)
  set_col (idx 0 0);
  set_col (idx 0 1);
  set_col (idx 0 2);
  
  (* iterate over interior nodes, filling in Jacobian entries *)
  for i=1 to ud.n-1-1 do
    (* set nodal value shortcuts *)
    let u = y.{idx i 0} in
    let v = y.{idx i 1} in
    let w = y.{idx i 2} in

    (* dependence on u at this node *)
    set_col (idx i 0);
    set (idx i 0) (2.0*.u*.v -. w -. 1.0);  (* fu wrt u *)
    set (idx i 1) (w -. 2.0*.u*.v);         (* fv wrt u *)
    set (idx i 2) (-.w);                    (* fw wrt u *)

    (* dependence on v at this node *)
    set_col (idx i 1);
    set (idx i 0) (u*.u);                   (* fu wrt v *)
    set (idx i 1) (-.u*.u);                 (* fv wrt v *)

    (* dependence on w at this node *)
    set_col (idx i 2);
    set (idx i 0) (-.u);                    (* fu wrt w *)
    set (idx i 1) (u);                      (* fv wrt w *)
    set (idx i 2) (-1.0/.ud.ep -. u)           (* fw wrt w *)
  done;

  (* set last matrix column to zero *)
  set_col (idx (ud.n-1) 0);
  set_col (idx (ud.n-1) 1);
  set_col (idx (ud.n-1) 2);

  (* end of data *)
  set_col ((idx (ud.n-1) 2)+1)

(* Jacobian routine to compute J(t,y) = df/dy. *)
let jac ud { Arkode.jac_y = y } j =
  let m, n, nnz = Sls.SparseMatrix.size j in

  (* ensure that Jac is the correct size *)
  if (m <> ud.n*3) || (n <> ud.n*3) then
    (printf "Jacobian calculation error: matrix is the wrong size!\n";
     raise Sundials.RecoverableFailure);
  
  (* Fill in the Laplace matrix *)
  laplace_matrix ud j;

  (* Create empty reaction Jacobian matrix (if not done already) *)
  (match ud.r with
   | Some _ -> ()
   | None ->
      try ud.r <- Some (Sls.SparseMatrix.make m n nnz)
      with _ ->
        (printf "Jacobian calculation error in allocating R matrix!\n";
         raise Sundials.RecoverableFailure));
      
  (* Add in the Jacobian of the reaction terms matrix *)
  (match ud.r with
   | None -> raise Sundials.RecoverableFailure
   | Some r -> begin
       reaction_jac ud y r;
       (* Add R to J *)
       try Sls.SparseMatrix.add j r
       with _ ->
         (printf "Jacobian calculation error in adding sparse matrices!\n";
          raise Sundials.RecoverableFailure)
     end)

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0 = 0.0 in    (* initial time *)
  let tf = 10.0 in   (* final time *)
  let nt = 10 in     (* total number of output times *)
  let nvar = 3 in

  let reltol = 1.0e-6 in     (* tolerances *)
  let abstol = 1.0e-10 in

  (* store the inputs in the UserData structure *)
  let n_mesh = 201 in
  let udata = {
    n = n_mesh;     (* spatial mesh size *)
    dx = 1.0/.float(n_mesh-1);    (* set spatial mesh spacing *)
    a = 0.6;
    b = 2.0;
    du = 0.025;
    dv = 0.025;
    dw = 0.025;
    ep = 1.0e-5;        (* stiffness parameter *)
    r  = None;
  } in

  (* set total allocated vector length *)
  let neq = nvar*n_mesh in

  (* Initial problem output *)
  printf("\n1D Brusselator PDE test problem (KLU solver):\n");
  printf "    N = %d,  NEQ = %d\n" udata.n neq;
  printf "    problem parameters:  a = %g,  b = %g,  ep = %g\n"
                                                    udata.a  udata.b  udata.ep;
  printf "    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n"
                                                    udata.du udata.dv udata.dw;
  printf "    reltol = %.1e,  abstol = %.1e\n\n" reltol abstol;

  (* Initialize data structures *)
  let data = RealArray.create neq in  (* Access data array for new NVector y *)
  let y = Nvector_serial.wrap data in (* Create serial vector for solution *)

  (* Set initial conditions into y *)
  let pi = 4.0*.atan(1.0) in
  for i=0 to n_mesh-1 do
    let fi = float i in
    data.{idx i 0} <-      udata.a       +. 0.1*.sin(pi*.fi*.udata.dx);  (* u *)
    data.{idx i 1} <- udata.b /. udata.a +. 0.1*.sin(pi*.fi*.udata.dx);  (* v *)
    data.{idx i 2} <-      udata.b       +. 0.1*.sin(pi*.fi*.udata.dx)   (* w *)
  done;

  (* Set mask array values for each solution component *)
  let data = RealArray.make neq 0.0 in
  let umask = Nvector_serial.wrap data in
  for i=0 to n_mesh-1 do
    data.{idx i 0} <- 1.0
  done;

  let data = RealArray.make neq 0.0 in
  let vmask = Nvector_serial.wrap data in
  for i=0 to n_mesh-1 do
    data.{idx i 1} <- 1.0
  done;

  let data = RealArray.make neq 0.0 in
  let wmask = Nvector_serial.wrap data in
  for i=0 to n_mesh-1 do
    data.{idx i 2} <- 1.0
  done;

  (* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time t0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. *)
  let nnz = 5*neq in
  let arkode_mem = Arkode.(
    init
      (Implicit (f udata, Newton (Arkode_klu.klu (jac udata) nnz), Nonlinear))
      (SStolerances (reltol, abstol))
      t0
      y
  ) in
  (* output spatial mesh to disk *)
  let fid = open_out "bruss_mesh.txt" in
  for i=0 to n_mesh-1 do
    fprintf fid "  %.16e\n" (udata.dx*.float i)
  done;
  close_out fid;

  (* Open output stream for results, access data arrays *)
  let ufid = open_out "bruss_u.txt" in
  let vfid = open_out "bruss_v.txt" in
  let wfid = open_out "bruss_w.txt" in

  (* output initial condition to disk *)
  let data = unwrap y in
  for i=0 to n_mesh-1 do
    fprintf ufid " %.16e" data.{idx i 0};
    fprintf vfid " %.16e" data.{idx i 1};
    fprintf wfid " %.16e" data.{idx i 2}
  done;
  fprintf ufid "\n";
  fprintf vfid "\n";
  fprintf wfid "\n";

  (* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let dTout = tf /. float nt in
  let tout = ref (t0+.dTout) in
  printf "        t      ||u||_rms   ||v||_rms   ||w||_rms\n";
  printf "   ----------------------------------------------\n";
  (try
     for iout=0 to nt-1 do
       (* call integrator *)
       let t, _ = Arkode.solve_normal arkode_mem !tout y in
 
       (* access/print solution statistics *)
       let u = n_vwl2norm y umask in
       let u = sqrt(u*.u/. float n_mesh) in
       let v = n_vwl2norm y vmask in
       let v = sqrt(v*.v/. float n_mesh) in
       let w = n_vwl2norm y wmask in
       let w = sqrt(w*.w/. float n_mesh) in
       printf "  %10.6f  %10.6f  %10.6f  %10.6f\n" t u v w;
       (* successful solve: update output time *)
       tout := min (!tout +. dTout) tf;
 
       (* output results to disk *)
       for i=0 to n_mesh-1 do
         fprintf ufid " %.16e" data.{idx i 0};
         fprintf vfid " %.16e" data.{idx i 1};
         fprintf wfid " %.16e" data.{idx i 2}
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
  let open Arkode in
  let nst      = get_num_steps arkode_mem in
  let nst_a    = get_num_step_attempts arkode_mem in
  let nfe, nfi = get_num_rhs_evals arkode_mem in
  let nsetups  = get_num_lin_solv_setups arkode_mem in
  let netf     = get_num_err_test_fails arkode_mem in
  let nni      = get_num_nonlin_solv_iters arkode_mem in
  let ncfn     = get_num_nonlin_solv_conv_fails arkode_mem in
  let nje      = Arkode_klu.get_num_jac_evals arkode_mem in

  printf "\nFinal Solver Statistics:\n";
  printf "   Internal solver steps = %d (attempted = %d)\n" nst nst_a;
  printf "   Total RHS evals:  Fe = %d,  Fi = %d\n" nfe nfi;
  printf "   Total linear solver setups = %d\n" nsetups;
  printf "   Total number of Jacobian evaluations = %d\n" nje;
  printf "   Total number of nonlinear iterations = %d\n" nni;
  printf "   Total number of nonlinear solver convergence failures = %d\n" ncfn;
  printf "   Total number of error test failures = %d\n" netf

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
