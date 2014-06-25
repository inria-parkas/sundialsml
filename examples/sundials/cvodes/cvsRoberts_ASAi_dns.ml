(*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2011/11/23 23:53:02 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2014.
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
 * iteration with the CVODE dense linear solver, and a user-supplied
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

module Cvode = Cvode_serial
module Quad = Cvodes_serial.Quadrature
module Adj = Cvodes_serial.Adjoint
module QuadAdj = Cvodes_serial.Adjoint.Quadrature
module RealArray = Cvode.RealArray
module Densemat = Dls.DenseMatrix

let printf = Printf.printf

(* Accessor macros *)

let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

let ijth v i j       = Densemat.get v (i - 1) (j - 1)
let set_ijth v i j e = Densemat.set v (i - 1) (j - 1) e

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

let steps = 150     (* number of steps between check points *)

let np    = 3       (* number of problem parameters         *)

let zero  = 0.0

(* Type : UserData *)

type user_data = { p : float array }

(* f routine. Compute f(t,y). *)

let f data t y ydot =
  let y1 = ith y 1
  and y2 = ith y 2
  and y3 = ith y 3
  and p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2) in

  let yd1 = -.p1*.y1 +. p2*.y2*.y3 in
  let yd3 = p3*.y2*.y2 in
  set_ith ydot 1 yd1;
  set_ith ydot 3 yd3; 
  set_ith ydot 2 (-.yd1 -. yd3)

(* Jacobian routine. Compute J(t,y). *)

let jac data { Cvode.jac_y = y } jmat =
  let y2 = ith y 2
  and y3 = ith y 3
  and p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  in
  set_ijth jmat 1 1 (-.p1);
  set_ijth jmat 1 2 (p2*.y3);
  set_ijth jmat 1 3 (p2*.y2);
  set_ijth jmat 2 1 ( p1);
  set_ijth jmat 2 2 (-.p2*.y3-.2.0*.p3*.y2);
  set_ijth jmat 2 3 (-.p2*.y2);
  set_ijth jmat 3 2 (2.0*.p3*.y2)

(* fQ routine. Compute fQ(t,y). *)

let fQ data t y qdot = set_ith qdot 1 (ith y 3)
 
(* EwtSet function. Computes the error weights at the current solution. *)

let ewt data y w =
  let rtol = rtol;
  and atol = Array.of_list [ atol1; atol2; atol3 ]
  in
  for i = 1 to 3 do
    let yy = ith y i in
    let ww = rtol *. (abs_float yy) +. atol.(i-1) in
    if ww <= 0.0 then failwith "ewt: ww negative";
    set_ith w i (1.0/.ww)
  done

(* fB routine. Compute fB(t,y,yB). *)

let fB data t y yB yBdot =
  let p1 = data.p.(0) (* The p vector *)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  and y2 = ith y 2
  and y3 = ith y 3
  and l1 = ith yB 1  (* The lambda vector *)
  and l2 = ith yB 2
  and l3 = ith yB 3
  in
  (* Temporary variables *)
  let l21 = l2-.l1
  and l32 = l3-.l2
  in
  (* Load yBdot *)
  set_ith yBdot 1 (-. p1*.l21);
  set_ith yBdot 2 (p2*.y3*.l21 -. 2.0*.p3*.y2*.l32);
  set_ith yBdot 3 (p2*.y2*.l21 -. 1.0)

(* JacB routine. Compute JB(t,y,yB). *)

let jacb data { Adj.jac_u = y } jbmat =
  let p1 = data.p.(0) (* The p vector *)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  and y2 = ith y 2
  and y3 = ith y 3
  in
  (* Load JB *)
  set_ijth jbmat 1 1 (p1);
  set_ijth jbmat 1 2 (-.p1); 
  set_ijth jbmat 2 1 (-.p2*.y3);
  set_ijth jbmat 2 2 (p2*.y3+.2.0*.p3*.y2);
  set_ijth jbmat 2 3 (-.2.0*.p3*.y2);
  set_ijth jbmat 3 1 (-.p2*.y2);
  set_ijth jbmat 3 2 (p2*.y2)

(* fQB routine. Compute integrand for quadratures *)

let fQB data t y yB qBdot =
  let y1 = ith y 1   (* The y vector *)
  and y2 = ith y 2
  and y3 = ith y 3
  and l1 = ith yB 1  (* The lambda vector *)
  and l2 = ith yB 2
  and l3 = ith yB 3
  in
  (* Temporary variables *)
  let l21 = l2-.l1
  and l32 = l3-.l2
  and y23 = y2*.y3
  in
  set_ith qBdot 1 (y1*.l21);
  set_ith qBdot 2 (-. y23*.l21);
  set_ith qBdot 3 (y2*.y2*.l32)

(* Print results after backward integration *)

let print_output tfinal yB qB =
  printf "--------------------------------------------------------\n";
  printf "tB0:        %12.4e\n" tfinal;
  printf "dG/dp:      %12.4e %12.4e %12.4e\n"
                                  (-.(ith qB 1)) (-.(ith qB 2)) (-.(ith qB 3));
  printf "lambda(t0): %12.4e %12.4e %12.4e\n" (ith yB 1) (ith yB 2) (ith yB 3);
  printf "--------------------------------------------------------\n\n"

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
  let y = RealArray.of_list [ 1.0; zero; zero ] in

  (* Initialize q *)
  let q = RealArray.of_list [ zero ] in

  (* Set the scalar realtive and absolute tolerances reltolQ and abstolQ *)
  let reltolQ = rtol
  and abstolQ = atolq in

  (* Create and allocate CVODES memory for forward run *)
  printf "Create and allocate CVODES memory for forward runs\n";

  let cvode_mem =
    Cvode.init Cvode.BDF (Cvode.Newton (Cvode.Dense (Some (jac data))))
      (Cvode.WFtolerances (ewt data)) (f data) ~t0:t0 y
  in

  Quad.init cvode_mem (fQ data) q;
  Quad.set_tolerances cvode_mem (Quad.SStolerances (reltolQ, abstolQ));

  (* Allocate global memory *)
  Adj.init cvode_mem steps Adj.IHermite (* Adj.IPolynomial *);

  (* Perform forward run *)
  printf "Forward integration ... ";
  let _, ncheck, _ = Adj.forward_normal cvode_mem tout y in
  let nst = Cvode.get_num_steps cvode_mem in
  
  printf "done ( nst = %d )\n" nst;
  printf "\nncheck = %d\n\n"  ncheck;
  ignore (Quad.get cvode_mem q);

  printf "--------------------------------------------------------\n";
  printf "G:          %12.4e \n" (ith q 1);
  printf "--------------------------------------------------------\n\n";

  (* Initialize yB *)
  let yB = RealArray.of_list [ zero; zero; zero ] in

  (* Initialize qB *)
  let qB = RealArray.of_list [ zero; zero; zero ] in

  (* Set the scalar relative tolerance reltolB *)
  let reltolB = rtol in

  (* Set the scalar absolute tolerance abstolB *)
  let abstolB = atoll in

  (* Set the scalar absolute tolerance abstolQB *)
  let abstolQB = atolq in

  (* Create and allocate CVODES memory for backward run *)
  printf "Create and allocate CVODES memory for backward run\n";

  let cvode_memB =
    Adj.init_backward cvode_mem Cvode.BDF
                                (Adj.Newton (Adj.Dense (Some (jacb data))))
                                (Adj.SStolerances (reltolB, abstolB))
                                (Adj.Basic (fB data))
                                tb1 yB
  in
  QuadAdj.init cvode_memB (QuadAdj.Basic (fQB data)) qB;
  QuadAdj.set_tolerances cvode_memB (QuadAdj.SStolerances (reltolB, abstolQB));

  (* Backward Integration *)
  printf "Backward integration ... ";
  Adj.backward_normal cvode_mem t0;
  let nstB = Adj.get_num_steps cvode_memB in

  printf "done ( nst = %d )\n"  nstB;

  ignore (Adj.get cvode_memB yB);
  ignore (QuadAdj.get cvode_memB qB);

  print_output tb1 yB qB;

  (* Reinitialize backward phase (new tB0) *)

  set_ith yB 1 zero;
  set_ith yB 2 zero;
  set_ith yB 3 zero;

  set_ith qB 1 zero;
  set_ith qB 2 zero;
  set_ith qB 3 zero;

  printf "Re-initialize CVODES memory for backward run\n";

  Adj.reinit cvode_memB tb2 yB;
  QuadAdj.reinit cvode_memB qB;

  printf "Backward integration ... ";

  Adj.backward_normal cvode_mem t0;
  let nstB = Adj.get_num_steps cvode_memB in

  printf "done ( nst = %d )\n"  nstB;

  ignore (Adj.get cvode_memB yB);
  ignore (QuadAdj.get cvode_memB qB);

  print_output tb2 yB qB

let _ = main ()
let _ = printf "Free memory\n\n"; Gc.compact ()

