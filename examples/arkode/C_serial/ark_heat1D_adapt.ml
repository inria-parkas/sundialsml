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
 * The following test simulates a simple 1D heat equation,
 *    u_t = k*u_xx + f
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  0
 * Dirichlet boundary conditions, i.e. 
 *    u_t(t,0) = u_t(t,1) = 0,
 * and a heating term of the form
 *    f = 2*exp(-200*(x-0.25)*(x-0.25))
 *        - exp(-400*(x-0.7)*(x-0.7))
 *        + exp(-500*(x-0.4)*(x-0.4))
 *        - 2*exp(-600*(x-0.55)*(x-0.55));
 * 
 * The spatial derivatives are computed using a three-point 
 * centered stencil (second order for a uniform mesh).  The data
 * is initially uniformly distributed over N points in the interval
 * [0, 1], but as the simulation proceeds the mesh is adapted.
 *
 * This program solves the problem with a DIRK method, solved with 
 * a Newton iteration and PCG linear solver, with a user-supplied 
 * Jacobian-vector product routine.
 *---------------------------------------------------------------*)

module RealArray = Sundials.RealArray
module LintArray = Sundials.LintArray
let printf = Printf.printf
let fprintf = Printf.fprintf
let n_vdotprod = Nvector_serial.Ops.n_vdotprod

exception IllegalMeshCreated

(* user data structure *)
type user_data = {
    mutable n  : int;         (* number of intervals   *)
    mutable x  : RealArray.t; (* current mesh *)
    k          : float;       (* diffusion coefficient *)
    refine_tol : float;       (* adaptivity tolerance *)
  }

(* Adapts the current mesh, using a simple adaptivity strategy of 
   refining when an approximation of the scaled second-derivative is 
   too large.  We only do this in one sweep, so no attempt is made to 
   ensure the resulting mesh meets these same criteria after adaptivity:
      y [input] -- the current solution vector
      nnew [output] -- the size of the new mesh
      udata [input] -- the current system information 
   The return for this function is a pointer to the new mesh. *)
let adapt_mesh { x = xold; n; refine_tol } (y : RealArray.t) =
  (* create marking array *)
  let marks = Array.make (n-1) 0 in

  (* /\* perform marking:  *)
  (*     0 -> leave alone *)
  (*     1 -> refine *)
  (* realtype ymax, ymin; *)
  (* for (i=0; i<(n-1); i++) { *)

  (*   /\* check for refinement *\/ *)
  (*   if (fabs(Y[i+1] - Y[i]) > refine_tol) { *)
  (*     marks[i] = 1; *)
  (*     continue; *)
  (*   } *)
  (* } *)

  (* perform marking: 
      0 -> leave alone
      1 -> refine *)
  for i=1 to n-1-1 do
    (* approximate scaled second-derivative *)
    let ydd = y.{i-1} -. 2.0*.y.{i} +. y.{i+1} in

    (* check for refinement *)
    if abs_float ydd > refine_tol then (marks.(i-1) <- 1;
                                        marks.(i) <- 1)
  done;

  (* allocate new mesh *)
  let nnew = n + Array.fold_left (fun s v->if v = 1 then s + 1 else s) 0 marks
  in
  let xnew = RealArray.create nnew in

  (* fill new mesh *)
  xnew.{0} <- xold.{0};    (* store endpoints *)
  xnew.{nnew-1} <- xold.{n-1};
  let j= ref 1 in
  (* iterate over old intervals *)
  for i=0 to n-1-1 do
    (* if mark is 0, reuse old interval *) 
    if marks.(i) = 0 then (xnew.{!j} <- xold.{i+1}; incr j)
    else if marks.(i) = 1 then begin
      (* if mark is 1, refine old interval *)
      xnew.{!j} <- 0.5 *. (xold.{i} +. xold.{i+1});
      incr j;
      xnew.{!j} <- xold.{i+1};
      incr j
    end
  done;

  (* verify that new mesh is legal *)
  for i=0 to nnew-1-1 do
    if xnew.{i+1} <= xnew.{i} then
      (fprintf stderr "adapt_mesh error: illegal mesh created\n";
       raise IllegalMeshCreated)
  done;
  (nnew, xnew)

(* Projects one vector onto another:
      Nold [input] -- the size of the old mesh
      xold [input] -- the old mesh
      yold [input] -- the vector defined over the old mesh
      nnew [input] -- the size of the new mesh
      xnew [input] -- the new mesh
      ynew [output] -- the vector defined over the new mesh
                       (allocated prior to calling project) *)
let project nold (xold : RealArray.t) (yold : RealArray.t)
            nnew (xnew : RealArray.t) (ynew : RealArray.t) =
  (* loop over new mesh, finding corresponding interval within old mesh, 
     and perform piecewise linear interpolation from yold to ynew *)
  let iv= ref 0 in
  for i=0 to nnew - 1 do
    (* find old interval, start with previous value since sorted *)
    (try
       for j=(!iv) to nold-1-1 do
         if xnew.{i} >= xold.{j} && xnew.{i} <= xold.{j+1} then
           (iv := j; raise Exit);
       done;
       iv := nold-1     (* just in case it wasn't found above *)
     with Exit -> ());

    (* perform interpolation *) 
    ynew.{i} <-
         yold.{!iv}  *.(xnew.{i}-.xold.{!iv+1})/.(xold.{!iv}  -.xold.{!iv+1})
      +. yold.{!iv+1}*.(xnew.{i}-.xold.{!iv})  /.(xold.{!iv+1}-.xold.{!iv});
  done

(* f routine to compute the ODE RHS function f(t,y). *)
let f { n; k; x } t (y : RealArray.t) (ydot : RealArray.t) =
  (* iterate over domain, computing all equations *)
  ydot.{0} <- 0.0;           (* left boundary condition *)
  for i=1 to n-1-1 do        (* interior *)
    let dxL = x.{i} -. x.{i-1} in
    let dxR = x.{i+1} -. x.{i} in
    ydot.{i} <- y.{i-1} *. k *. 2.0 /. (dxL *. (dxL +. dxR))
             -. y.{i}   *. k *. 2.0 /. (dxL *. dxR)
             +. y.{i+1} *. k *. 2.0 /. (dxR *. (dxL +. dxR))
  done;
  ydot.{n-1} <- 0.0;         (* right boundary condition *)

  (* source term *)
  (* RealArray.mapi (fun i ydot ->
    ydot +. (2.0 *. exp(-200.0 *. (x.{i} -. 0.25) *. (x.{i} -. 0.25))
                 -. exp(-400.0 *. (x.{i} -. 0.70) *. (x.{i} -. 0.70))
                 +. exp(-500.0 *. (x.{i} -. 0.40) *. (x.{i} -. 0.40))
          -. 2.0 *. exp(-600.0 *. (x.{i} -. 0.55) *. (x.{i} -. 0.55)))) ydot *)
  for i=0 to n-1-1 do
    ydot.{i} <- ydot.{i} +.
            (2.0 *. exp(-200.0 *. (x.{i} -. 0.25) *. (x.{i} -. 0.25))
                 -. exp(-400.0 *. (x.{i} -. 0.7 ) *. (x.{i} -. 0.7 ))
                 +. exp(-500.0 *. (x.{i} -. 0.4 ) *. (x.{i} -. 0.4 ))
          -. 2.0 *. exp(-600.0 *. (x.{i} -. 0.55) *. (x.{i} -. 0.55)))
  done

(* Jacobian routine to compute J(t,y) = df/dy. *)
let jac { n; k; x } _ (v : RealArray.t) (jv : RealArray.t) =
  (* iterate over domain, computing all Jacobian-vector products *)
  jv.{0} <- 0.0;
  for i=1 to n-1-1 do
    let dxL = x.{i} -. x.{i-1} in
    let dxR = x.{i+1} -. x.{i} in
    jv.{i} <- v.{i-1} *. k *. 2.0 /. (dxL *. (dxL +. dxR)) 
           -. v.{i}   *. k *. 2.0 /. (dxL *. dxR)
           +. v.{i+1} *. k *. 2.0 /. (dxR *. (dxL +. dxR))
  done;
  jv.{n-1} <- 0.0

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0     = 0.0 in       (* initial time *)
  let tf     = 1.0 in       (* final time *)
  let rtol   = 1.e-3 in     (* relative tolerance *)
  let atol   = 1.e-10 in    (* absolute tolerance *)
  let hscale = 1.0 in       (* time step change factor on resizes *)
  let n_mesh = 21 in        (* initial spatial mesh size *)

  (* allocate and fill initial udata structure *)
  let udata = {
    n = n_mesh;
    x = RealArray.init n_mesh (fun i-> float i /. float (n_mesh-1));
    k = 0.5;             (* heat conductivity *)
    refine_tol = 3.e-3;  (* adaptivity refinement tolerance *)
  } in

  (* Initial problem output *)
  printf "\n1D adaptive Heat PDE test problem:\n";
  printf "  diffusion coefficient:  k = %g\n" udata.k;
  printf "  initial N = %d\n" n_mesh;

  (* Initialize data structures *)

  (* Create initial serial vector for solution *)
  (* Set initial conditions *)
  let y = Nvector_serial.make n_mesh 0.0 in

  (* output mesh to disk *)
  let xfid = open_out "heat_mesh.txt" in

  (* output initial mesh to disk *)
  RealArray.iter (fun d -> fprintf xfid " %.16e" d) udata.x;
  fprintf xfid "\n";

  (* Open output stream for results, access data array *)
  let ufid = open_out "heat1D.txt" in

  (* output initial condition to disk *)
  let data = Nvector.unwrap y in
  RealArray.iter (fun d -> fprintf ufid " %.16e" d) data;
  fprintf ufid "\n";

  (* Initialize the integrator memory *)
  let jac = jac udata in
  let arkode_mem = Arkode.(
    init
      (Implicit (f udata,
                 Newton Spils.(pcg ~jac_times_vec:jac prec_none),
                 Nonlinear))
      (SStolerances (rtol, atol))
      t0
      y
  ) in
  Arkode.set_max_num_steps arkode_mem 10000;      (* Increase max num steps  *)
  Arkode.(set_adaptivity_method arkode_mem
            (Icontroller { ks = None; method_order = false}));
  Arkode.(set_predictor_method arkode_mem TrivialPredictor);

  (* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached *)
  let border = String.make 88 '-' in
  printf "  iout          dt_old                 dt_new               ||u||_rms       N   NNI  NLI\n";
  printf " %s\n" border;
  printf " %4d  %19.15e  %19.15e  %19.15e  %d   %2d  %3d\n"
         0 0.0 0.0 (sqrt(n_vdotprod y y /. float udata.n)) udata.n 0 0;
  let rec loop t newdt y iout nni_cur nni_tot nli_tot =
    if t >= tf then iout, nni_tot, nli_tot
    else begin
      (* "set" routines *)
      Arkode.set_stop_time arkode_mem tf;
      Arkode.set_init_step arkode_mem newdt;

      (* call integrator *)
      let t, _ = Arkode.solve_one_step arkode_mem tf y in

      (* "get" routines *)
      let olddt = Arkode.get_last_step arkode_mem in
      let newdt = Arkode.get_current_step arkode_mem in
      let nni   = Arkode.get_num_nonlin_solv_iters arkode_mem in
      let nli   = Arkode.Spils.get_num_lin_iters arkode_mem in

      (* print current solution stats *)
      printf " %4d  %19.15e  %19.15e  %19.15e  %d   %2d  %3d\n"
             (iout + 1) olddt newdt (sqrt(n_vdotprod y y /. float udata.n))
             udata.n (nni-nni_cur) nli;

      (* output results and current mesh to disk *)
      let data = Nvector.unwrap y in
      RealArray.iter (fun d -> fprintf ufid " %.16e" d) data;
      fprintf ufid "\n";
      RealArray.iter (fun d -> fprintf xfid " %.16e" d) udata.x;
      fprintf xfid "\n";

      (* adapt the spatial mesh *)
      let nnew, xnew = adapt_mesh udata data in

      (* create N_Vector of new length *)
      let y2 = Nvector_serial.make nnew 0.0 in
      
      (* project solution onto new mesh *)
      project udata.n udata.x data nnew xnew (Nvector.unwrap y2);

      (* swap x and xnew so that new mesh is stored in udata structure *)
      udata.x <- xnew;
      udata.n <- nnew;   (* store size of new mesh *)

      (* call ARKodeResize to notify integrator of change in mesh *)
      Arkode.(resize arkode_mem
                     ~linsolv:Spils.(pcg ~jac_times_vec:jac prec_none)
                     (SStolerances (rtol, atol))
                     hscale y2 t);

      loop t newdt y2 (iout + 1) nni nni (nli_tot + nli)
    end
  in
  let iout, nni_tot, nli_tot = loop t0 0.0 y 0 0 0 0 in
  printf " %s\n" border;

  (* print some final statistics *)
  printf " Final solver statistics:\n";
  printf "   Total number of time steps = %d\n" iout;
  printf "   Total nonlinear iterations = %d\n" nni_tot;
  printf "   Total linear iterations    = %d\n\n" nli_tot

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
