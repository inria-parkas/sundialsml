(*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2011/11/23 23:53:02 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Ting Yan @ SMU
 *      Based on cvsRoberts_ASAi_dns.c and modified to use KLU
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Dec 2016.
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Adjoint sensitivity example problem.
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES. The problem is from chemical
 * kinetics, and consists of the following three rate equations.
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The reaction rates are:
 * p1=0.04, p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODE KLU linear solver, and a user-supplied
 * Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * 
 * Optionally, CVODES can compute sensitivities with respect to
 * the problem parameters p1, p2, and p3 of the following quantity:
 *   G = int_t0^t1 g(t,p,y) dt
 * where
 *   g(t,p,y) = y3
 *        
 * The gradient dG/dp is obtained as:
 *   dG/dp = int_t0^t1 (g_p - lambda^T f_p ) dt - lambda^T(t0)*y0_p
 *         = - xi^T(t0) - lambda^T(t0)*y0_p
 * where lambda and xi are solutions of:
 *   d(lambda)/dt = - (f_y)^T * lambda - (g_y)^T
 *   lambda(t1) = 0
 * and
 *   d(xi)/dt = - (f_p)^T * lambda + (g_p)^T
 *   xi(t1) = 0
 * 
 * During the backward integration, CVODES also evaluates G as
 *   G = - phi(t0)
 * where
 *   d(phi)/dt = g(t,y,p)
 *   phi(t1) = 0
 * -----------------------------------------------------------------
 *)

module Quad = Cvodes.Quadrature
module Adj = Cvodes.Adjoint
module QuadAdj = Cvodes.Adjoint.Quadrature
module RealArray = Sundials.RealArray
module Densemat = Dls.DenseMatrix
let unwrap = Nvector.unwrap

let printf = Printf.printf

(* Accessor macros *)

let ith (v : RealArray.t) i = v.{i - 1}
let set_ith (v : RealArray.t) i e = v.{i - 1} <- e

(* Problem Constants *)

let neq   = 3       (* number of equations                  *)

let rtol  = 1.0e-6  (* scalar relative tolerance            *)

let atol1 = 1.0e-8  (* vector absolute tolerance components *)
let atol2 = 1.0e-14
let atol3 = 1.0e-6

let atoll = 1.0e-8  (* absolute tolerance for adjoint vars. *)
let atolq = 1.0e-6  (* absolute tolerance for quadratures   *)

let t0    = 0.0     (* initial time                         *)
let tout  = 4.0e7   (* final time                           *)

let tb1   = 4.0e7   (* starting point for adjoint problem   *)
let tb2   = 50.0    (* starting point for adjoint problem   *)
let tBout1= 40.0    (* intermediate t for adjoint problem   *)

let steps = 150     (* number of steps between check points *)

let np    = 3       (* number of problem parameters         *)

let zero  = 0.0

(* Type : UserData *)

type user_data = { p : float array }

(* f routine. Compute f(t,y). *)

let f data t (y : RealArray.t) (ydot : RealArray.t) =
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2) in

  let yd1 = -.p1*.y.{0} +. p2*.y.{1}*.y.{2} in
  let yd3 = p3*.y.{1}*.y.{1} in
  ydot.{0} <- yd1;
  ydot.{1} <- (-.yd1 -. yd3);
  ydot.{2} <- yd3

(* Jacobian routine. Compute J(t,y). *)

let jac data {Cvode.jac_y = (y : RealArray.t)} smat =
  let set_col = Sls.SparseMatrix.set_col smat in
  let set = Sls.SparseMatrix.set smat in
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  in
  Sls.SparseMatrix.set_to_zero smat;

  set_col 0 0;
  set_col 1 3;
  set_col 2 6;
  set_col 3 9;

  set 0 0 (-.p1);
  set 1 1  p1;
  set 2 2  0.00;

  set 3 0 (p2 *. y.{2});
  set 4 1 (-.p2 *. y.{2} -. 2.0 *. p3 *. y.{1});
  set 5 2 (2.0 *. y.{1});

  set 6 0 (p2 *. y.{1});
  set 7 1 (-.p2 *. y.{1});
  set 8 2 0.00

(* fQ routine. Compute fQ(t,y). *)

let fQ data t (y : RealArray.t) (qdot : RealArray.t) = qdot.{0} <- y.{2}
 
(* EwtSet function. Computes the error weights at the current solution. *)

let atol = Array.of_list [ atol1; atol2; atol3 ]
let ewt data (y : RealArray.t) (w : RealArray.t) =
  let rtol = rtol;
  in
  for i = 0 to 2 do
    let ww = rtol *. (abs_float y.{i}) +. atol.(i) in
    if ww <= 0.0 then raise Sundials.NonPositiveEwt;
    w.{i} <- (1.0/.ww)
  done

(* fB routine. Compute fB(t,y,yB). *)

let fB : user_data -> RealArray.t Adj.brhsfn_no_sens =
  fun data args yBdot ->
  let y = args.Adj.y
  and yB = args.Adj.yb
  in
  let p1 = data.p.(0) (* The p vector *)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  and l1 = yB.{0}  (* The lambda vector *)
  and l2 = yB.{1}
  and l3 = yB.{2}
  in
  (* Temporary variables *)
  let l21 = l2-.l1
  and l32 = l3-.l2
  in
  (* Load yBdot *)
  yBdot.{0} <- -. p1*.l21;
  yBdot.{1} <- p2*.y.{2}*.l21 -. 2.0*.p3*.y.{1}*.l32;
  yBdot.{2} <- p2*.y.{1}*.l21 -. 1.0

(* JacB routine. Compute JB(t,y,yB). *)

let jacb data { Adj.jac_y = (y : RealArray.t) } smat =
  let set_col = Sls.SparseMatrix.set_col smat in
  let set = Sls.SparseMatrix.set smat in
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  in
  Sls.SparseMatrix.set_to_zero smat;

  set_col 0 0;
  set_col 1 3;
  set_col 2 6;
  set_col 3 9;

  set 0 0  p1;
  set 1 1  (-.p2 *. y.{2});
  set 2 2  (-.p2 *. y.{1});

  set 3 0 (-.p1);
  set 4 1 (p2 *. y.{2} +. 2.0 *. p3 *. y.{1});
  set 5 2 (p2 *. y.{1});

  set 6 0 0.0;
  set 7 1 (-2.0 *. p3 *. y.{1});
  set 8 2 0.00

(* fQB routine. Compute integrand for quadratures *)

let fQB : user_data -> RealArray.t QuadAdj.bquadrhsfn_no_sens =
  fun data args qBdot ->
  let y = args.QuadAdj.y
  and yB = args.QuadAdj.yb
  in
  let l1 = yB.{0}  (* The lambda vector *)
  and l2 = yB.{1}
  and l3 = yB.{2}
  in
  (* Temporary variables *)
  let l21 = l2-.l1
  and l32 = l3-.l2
  and y23 = y.{1}*.y.{2}
  in
  qBdot.{0} <- y.{0}*.l21;
  qBdot.{1} <- -. y23*.l21;
  qBdot.{2} <- y.{1}*.y.{1}*.l32

(* Print intermediate results during backward integration *)

let print_output1 time t y yB =
  printf "--------------------------------------------------------\n";
  printf "returned t: %12.4e\n" time;
  printf "tout:       %12.4e\n" t;
  printf "lambda(t):  %12.4e %12.4e %12.4e\n" (ith yB 1) (ith yB 2) (ith yB 3);
  printf "y(t):       %12.4e %12.4e %12.4e\n" (ith y 1) (ith y 2) (ith y 3);
  printf "--------------------------------------------------------\n\n"

(* Print results after backward integration *)

let print_output tfinal y yB qB =
  printf "--------------------------------------------------------\n";
  (match Sundials.sundials_version with
   | 2,5,_ ->
       printf "tB0:        %12.4e\n" tfinal;
       printf "dG/dp:      %12.4e %12.4e %12.4e\n"
                                   (-.(ith qB 1)) (-.(ith qB 2)) (-.(ith qB 3));
       printf "lambda(t0): %12.4e %12.4e %12.4e\n"
                                   (ith yB 1) (ith yB 2) (ith yB 3)
   | _ ->
       printf "returned t: %12.4e\n" tfinal;
       printf "lambda(t0): %12.4e %12.4e %12.4e\n"
                                 (ith yB 1) (ith yB 2) (ith yB 3);
       printf "y(t0):      %12.4e %12.4e %12.4e\n"
                                 (ith y 1) (ith y 2) (ith y 3);
       printf "dG/dp:      %12.4e %12.4e %12.4e\n"
                                 (-.(ith qB 1)) (-.(ith qB 2)) (-.(ith qB 3)));
  printf "--------------------------------------------------------\n\n"

let print_head tB0 =
  printf "Backward integration from tB0 = %12.4e\n\n" tB0

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* Print problem description *)
  printf "\nAdjoint Sensitivity Example for Chemical Kinetics\n";
  printf "-------------------------------------------------\n\n";
  printf "ODE: dy1/dt = -p1*y1 + p2*y2*y3\n";
  printf "     dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2\n";
  printf "     dy3/dt =  p3*(y2)^2\n\n";
  printf "Find dG/dp for\n";
  printf "     G = int_t0^tB0 g(t,p,y) dt\n";
  printf "     g(t,p,y) = y3\n\n\n";

  (* User data structure *)
  let data = { p = Array.of_list [ 0.04; 1.0e4; 3.0e7 ] } in

  (* Initialize y *)
  let ydata = RealArray.of_list [ 1.0; zero; zero ] in
  let y = Nvector_serial.wrap ydata in

  (* Initialize q *)
  let qdata = RealArray.of_list [ zero ] in
  let q = Nvector_serial.wrap qdata in

  (* Set the scalar realtive and absolute tolerances reltolQ and abstolQ *)
  let reltolQ = rtol
  and abstolQ = atolq in

  (* Create and allocate CVODES memory for forward run *)
  printf "Create and allocate CVODES memory for forward runs\n";

  let nnz = neq * neq in
  let cvode_mem =
    Cvode.(init BDF (Newton Sls.Klu.(solver_csc (jac data) nnz))
                (WFtolerances (ewt data)) (f data) t0 y)
  in

  Quad.init cvode_mem (fQ data) q;
  Quad.(set_tolerances cvode_mem (SStolerances (reltolQ, abstolQ)));

  (* Allocate global memory *)
  Adj.(init cvode_mem steps IHermite) (* IPolynomial *);

  (* Perform forward run *)
  printf "Forward integration ... ";
  let _, ncheck, _ = Adj.forward_normal cvode_mem tout y in
  let nst = Cvode.get_num_steps cvode_mem in
  
  printf "done ( nst = %d )\n" nst;
  printf "\nncheck = %d\n\n"  ncheck;
  ignore (Quad.get cvode_mem q);

  printf "--------------------------------------------------------\n";
  printf "G:          %12.4e \n" (ith qdata 1);
  printf "--------------------------------------------------------\n\n";

  (* Initialize yB *)
  let yBdata = RealArray.of_list [ zero; zero; zero ] in
  let yB = Nvector_serial.wrap yBdata in

  (* Initialize qB *)
  let qBdata = RealArray.of_list [ zero; zero; zero ] in
  let qB = Nvector_serial.wrap qBdata in

  (* Set the scalar relative tolerance reltolB *)
  let reltolB = rtol in

  (* Set the scalar absolute tolerance abstolB *)
  let abstolB = atoll in

  (* Set the scalar absolute tolerance abstolQB *)
  let abstolQB = atolq in

  (* Create and allocate CVODES memory for backward run *)
  printf "Create and allocate CVODES memory for backward run\n";

  let cvode_memB =
    Adj.(init_backward
          cvode_mem Cvode.BDF
                    (Newton Sls.(Klu.solver_csc (NoSens (jacb data)) nnz))
                    (SStolerances (reltolB, abstolB))
                    (NoSens (fB data))
                    tb1 yB)
  in
  QuadAdj.(init cvode_memB (NoSens (fQB data)) qB);
  QuadAdj.(set_tolerances cvode_memB (SStolerances (reltolB, abstolQB)));

  (* Backward Integration *)

  (* First get results at t = tBout1 *)
  (match Sundials.sundials_version with
   | 2,5,_ -> printf "Backward integration ... "
   | _ ->
      print_head tb1;
      Adj.backward_normal cvode_mem tBout1;
      let time = Adj.get cvode_memB yB in
      Adj.get_y cvode_mem y tBout1;
      print_output1 time tBout1 ydata yBdata);

  (* Then at t = T0 *)

  Adj.backward_normal cvode_mem t0;
  let nstB = Adj.get_num_steps cvode_memB in

  (match Sundials.sundials_version with
   | 2,5,_ -> printf "done ( nst = %d )\n"  nstB
   | _     -> printf "Done ( nst = %d )\n"  nstB);

  ignore (Adj.get cvode_memB yB);
  let time = QuadAdj.get cvode_memB qB in
  Adj.get_y cvode_mem y t0;

  (match Sundials.sundials_version with
   | 2,5,_ -> print_output tb1 ydata yBdata qBdata
   | _     -> print_output time ydata yBdata qBdata);

  (* Reinitialize backward phase (new tB0) *)

  set_ith yBdata 1 zero;
  set_ith yBdata 2 zero;
  set_ith yBdata 3 zero;

  set_ith qBdata 1 zero;
  set_ith qBdata 2 zero;
  set_ith qBdata 3 zero;

  printf "Re-initialize CVODES memory for backward run\n";

  Adj.reinit cvode_memB tb2 yB;
  QuadAdj.reinit cvode_memB qB;

  (* First get results at t = tBout1 *)
  (match Sundials.sundials_version with
   | 2,5,_ -> printf "Backward integration ... "
   | _     ->
      print_head tb2;
      Adj.backward_normal cvode_mem tBout1;
      let time = Adj.get cvode_memB yB in
      Adj.get_y cvode_mem y tBout1;
      print_output1 time tBout1 ydata yBdata);

  (* Then at t = T0 *)

  Adj.backward_normal cvode_mem t0;
  let nstB = Adj.get_num_steps cvode_memB in

  (match Sundials.sundials_version with
   | 2,5,_ -> printf "done ( nst = %d )\n"  nstB
   | _     -> printf "Done ( nst = %d )\n"  nstB);

  ignore (Adj.get cvode_memB yB);
  let time = QuadAdj.get cvode_memB qB in
  Adj.get_y cvode_mem y t0;

  (match Sundials.sundials_version with
   | 2,5,_ -> print_output tb2 ydata yBdata qBdata
   | _     -> print_output time ydata yBdata qBdata);

  printf "Free memory\n\n"

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
