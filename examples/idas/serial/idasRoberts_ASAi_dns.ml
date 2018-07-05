(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2010/12/01 23:05:10 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Jul 2014.
 * -----------------------------------------------------------------
 * Copyright (c) 2007, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Adjoint sensitivity example problem.
 *
 * This simple example problem for IDAS, due to Robertson,
 * is from chemical kinetics, and consists of the following three
 * equations:
 *
 *      dy1/dt + p1*y1 - p2*y2*y3            = 0
 *      dy2/dt - p1*y1 + p2*y2*y3 + p3*y2**2 = 0
 *                 y1  +  y2  +  y3  -  1    = 0
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7
 *
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 *
 * IDAS can also compute sensitivities with respect to
 * the problem parameters p1, p2, and p3 of the following quantity:
 *   G = int_t0^t1 g(t,p,y) dt
 * where
 *   g(t,p,y) = y3
 *
 * The gradient dG/dp is obtained as:
 *   dG/dp = int_t0^t1 (g_p - lambda^T F_p ) dt -
 *           lambda^T*F_y'*y_p | _t0^t1
 *         = int_t0^t1 (lambda^T*F_p) dt
 * where lambda and are solutions of the adjoint system:
 *   d(lambda^T * F_y' )/dt -lambda^T F_y = -g_y
 *
 * During the backward integration, IDAS also evaluates G as
 *   G = - phi(t0)
 * where
 *   d(phi)/dt = g(t,y,p)
 *   phi(t1) = 0
 * -----------------------------------------------------------------
 *)

open Sundials

module Quad = Idas.Quadrature
module Sens = Idas.Sensitivity
module Adjoint = Idas.Adjoint
module AdjQuad = Adjoint.Quadrature
module VarId = Ida.VarId

let printf = Printf.printf

let nvconst = Nvector_serial.DataOps.n_vconst
let nvscale = Nvector_serial.DataOps.n_vscale

let unwrap = Nvector.unwrap
let wrap = Nvector_serial.wrap

(* Problem Constants *)

let neq =      3             (* number of equations                  *)

let rtol =     1e-06 (* scalar relative tolerance            *)

let atol1 =    1e-08 (* vector absolute tolerance components *)
let atol2 =    1e-12
let atol3 =    1e-08

let atola =    1e-08 (* absolute tolerance for adjoint vars. *)
let atolq =    1e-06 (* absolute tolerance for quadratures   *)

let t0 =       0.0   (* initial time                         *)
let tout =     4e10  (* final time                           *)

let tb1 =      50.0  (* starting point for adjoint problem   *)
let tb2 =      tout          (* starting point for adjoint problem   *)

let t1b =      49.0  (* for IDACalcICB                       *)

let steps =    100           (* number of steps between check points *)

let np =       3             (* number of problem parameters         *)


(* Type : UserData *)

type user_data = { p : RealArray.t }

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDAS
 *--------------------------------------------------------------------
 *)

(*
 * f routine. Compute f(t,y).
*)

let res data t (yy : RealArray.t) (yp : RealArray.t) (rval : RealArray.t) =
  let y1  = yy.{0}
  and y2  = yy.{1}
  and y3  = yy.{2}
  and yp1 = yp.{0}
  and yp2 = yp.{1}

  and p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}
  in

  rval.{0} <- p1*.y1-.p2*.y2*.y3;
  rval.{1} <- -.rval.{0} +. p3*.y2*.y2 +. yp2;
  rval.{0} <- rval.{0} +. yp1;
  rval.{2} <- y1+.y2+.y3-.1.0

(*
 * Jacobian routine. Compute J(t,y).
*)

let jac data jac_arg j =
  let cj = jac_arg.Ida.jac_coef
  and (yy : RealArray.t) = jac_arg.Ida.jac_y
  in

  let y2 = yy.{1}
  and y3 = yy.{2}

  and p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}
  in
  let set = Matrix.Dense.set j in
  set 0 0 (p1+.cj);
  set 1 0 (-.p1);
  set 2 0 (1.0);

  set 0 1 (-.p2*.y3);
  set 1 1 (p2*.y3+.2.0*.p3*.y2+.cj);
  set 2 1 (1.0);

  set 0 2 (-.p2*.y2);
  set 1 2 (p2*.y2);
  set 2 2 (1.0)

(*
 * rhsQ routine. Compute fQ(t,y).
*)

let rhsQ data t (yy : RealArray.t) yp (qdot : RealArray.t) =
  qdot.{0} <- yy.{2}

(* EwtSet function. Computes the error weights at the current solution. *)

let ewt data (y : RealArray.t) (w : RealArray.t) =
  let atol = [|atol1; atol2; atol3|] in

  for i = 1 to 3 do
    let yy = y.{i-1} in
    let ww = rtol *. abs_float yy +. atol.(i-1) in
    if ww <= 0.0 then raise NonPositiveEwt;
    w.{i-1} <- 1.0/.ww
  done


(* resB routine. *)

let resB : user_data -> RealArray.t Adjoint.bresfn_no_sens =
  fun data { Adjoint.y = yy; Adjoint.yb = yyB; Adjoint.yb' = ypB } rrB ->

  (* The p vector *)
  let p1 = data.p.{0} and p2 = data.p.{1} and p3 = data.p.{2} in

  (* The y  vector *)
  let y2 = yy.{1} and y3 = yy.{2} in

  (* The lambda vector *)
  let l1 = yyB.{0} and l2 = yyB.{1} and l3 = yyB.{2} in

  (* The lambda dot vector *)
  let lp1 = ypB.{0} and lp2 = ypB.{1} in

  (* Temporary variables *)
  let l21 = l2-.l1 in

  (* Load residual. *)
  rrB.{0} <- lp1 +. p1*.l21 -. l3;
  rrB.{1} <- lp2 -. p2*.y3*.l21 -. 2.0*.p3*.y2*.l2-.l3;
  rrB.{2} <- -. p2*.y2*.l21 -.l3 +. 1.0

(*Jacobian for backward problem. *)
let jacB data { Adjoint.jac_coef = cj; Adjoint.jac_y = (yy : RealArray.t) } jB =
  let y2 = yy.{1} and y3 = yy.{2} in

  let p1 = data.p.{0} and p2 = data.p.{1} and p3 = data.p.{2} in

  let set = Matrix.Dense.set jB in
  set 0 0 (-.p1+.cj);
  set 0 1 (p1);
  set 0 2 (-.1.0);

  set 1 0 (p2*.y3);
  set 1 1 (-.(p2*.y3+.2.0*.p3*.y2)+.cj);
  set 1 2 (-.1.0);

  set 2 0 (p2*.y2);
  set 2 1 (-.p2*.y2);
  set 2 2 (-.1.0)

let rhsQB : user_data -> RealArray.t AdjQuad.bquadrhsfn_no_sens =
  fun data { AdjQuad.y = yy; AdjQuad.yb = yyB } rrQB ->

  (* The y vector *)
  let y1 = yy.{0} and y2 = yy.{1} and y3 = yy.{2} in

  (* The lambda vector *)
  let l1 = yyB.{0} and l2 = yyB.{1}
  in

  (* Temporary variables *)
  let l21 = l2-.l1 in

  rrQB.{0} <- y1*.l21;
  rrQB.{1} <- -.y3*.y2*.l21;
  rrQB.{2} <- -.y2*.y2*.l2

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(* Print results after backward integration *)

let print_output tfinal yB ypB qB =
  printf "--------------------------------------------------------\n";
  printf "tB0:        %12.4e\n" tfinal;
  printf "dG/dp:      %12.4e %12.4e %12.4e\n"
         (-.qB.{0}) (-.qB.{1}) (-.qB.{2});
  printf "lambda(t0): %12.4e %12.4e %12.4e\n"
         yB.{0} yB.{1} yB.{2};
  printf "--------------------------------------------------------\n\n"


(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* Print problem description *)
  print_string "\nAdjoint Sensitivity Example for Chemical Kinetics\n";
  print_string "-------------------------------------------------\n\n";
  print_string "DAE: dy1/dt + p1*y1 - p2*y2*y3 = 0\n";
  print_string "     dy2/dt - p1*y1 + p2*y2*y3 + p3*(y2)^2 = 0\n";
  print_string "               y1  +  y2  +  y3 = 0\n\n";
  print_string "Find dG/dp for\n";
  print_string "     G = int_t0^tB0 g(t,p,y) dt\n";
  print_string "     g(t,p,y) = y3\n\n\n";

  (* User data structure *)
  let data = { p = RealArray.of_array [|0.04; 1.0e4; 3.0e7|] } in

  (* Initialize y *)
  let yy = RealArray.of_array [|1.0; 0.0; 0.0|] in

  (* Initialize yprime *)
  let yp = RealArray.of_array [|-0.04; 0.04; 0.0|] in

  (* Initialize q *)
  let q = RealArray.of_array [|0.0|] in

  (* Wrap arrays into nvectors. *)
  let wyy = wrap yy
  and wyp = wrap yp
  and wq  = wrap q
  in

  (* Set the scalar realtive and absolute tolerances reltolQ and abstolQ *)
  let reltolQ = rtol
  and abstolQ = atolq
  in

  (* Create and allocate IDAS memory for forward run *)
  print_string "Create and allocate IDAS memory for forward runs\n";

  let m = Matrix.dense neq in
  let ida_mem =
    Ida.(init Dls.(solver ~jac:(jac data) (dense wyy m))
              (WFtolerances (ewt data))
              (res data)
              t0
              wyy wyp)
  in

  Quad.init ida_mem (rhsQ data) wq;
  Quad.(set_tolerances ida_mem (SStolerances (reltolQ,abstolQ)));

  (* Allocate global memory *)

  Adjoint.(init ida_mem steps IHermite);

  (* Perform forward run *)
  print_string "Forward integration ... ";

  (* Integrate till TB1 and get the solution (y, y') at that time. *)
  let _ = Adjoint.forward_normal ida_mem tb1 wyy wyp in

  let yyTB1 = RealArray.copy yy
  and ypTB1 = RealArray.copy yp
  in
  (* Save the states at t=TB1. *)
  nvscale 1.0 yy yyTB1;
  nvscale 1.0 yp ypTB1;

  (* Continue integrating till TOUT is reached. *)
  let _ = Adjoint.forward_normal ida_mem tout wyy wyp in

  let nst = Ida.get_num_steps ida_mem in

  printf "done ( nst = %d )\n" nst;

  let _ = Quad.get ida_mem wq in

  print_string "--------------------------------------------------------\n";
  printf "G:          %12.4e \n" q.{0};
  print_string "--------------------------------------------------------\n\n";

  (* Create BACKWARD problem. *)

  (* Allocate yB (i.e. lambda_0). *)
  (* Consistently initialize yB. *)
  let yB = RealArray.of_array [|0.0;0.0;1.0|] in

  (* Allocate ypB (i.e. lambda'_0). *)
  (* Consistently initialize ypB. *)
  let ypB = RealArray.of_array [|1.0;1.0;0.0|] in

  (* Set the scalar relative tolerance reltolB *)
  let reltolB = rtol in

  (* Set the scalar absolute tolerance abstolB *)
  let abstolB = atola in

  (* Set the scalar absolute tolerance abstolQB *)
  let abstolQB = atolq in

  (* Create and allocate IDAS memory for backward run *)
  print_string "Create and allocate IDAS memory for backward run\n";

  let wyB  = wrap yB
  and wypB = wrap ypB
  in

  let m = Matrix.dense neq in
  let indexB =
    Adjoint.(init_backward ida_mem
               Dls.(solver ~jac:(NoSens (jacB data)) (dense wyB m))
               (SStolerances (reltolB, abstolB))
               (NoSens (resB data))
               tb2 wyB wypB)
  in
  Adjoint.set_max_num_steps indexB 1000;

  (* Quadrature for backward problem. *)

  (* Initialize qB *)
  let qB = RealArray.of_array [|0.0; 0.0; 0.0|] in
  let wqB = wrap qB in

  AdjQuad.(init indexB (NoSens (rhsQB data)) wqB);
  (* Include quadratures in error control. *)
  AdjQuad.(set_tolerances indexB (SStolerances (reltolB, abstolQB)));

  (* Backward Integration *)
  print_string "Backward integration ... ";

  Adjoint.backward_normal ida_mem t0;

  let nstB = Adjoint.get_num_steps indexB in
  printf "done ( nst = %d )\n" nstB;

  let _ = Adjoint.get indexB wyB wypB in

  let _ = AdjQuad.get indexB wqB in

  print_output tb2 yB ypB qB;


  (* Reinitialize backward phase and start from a different time (TB1). *)
  print_string "Re-initialize IDAS memory for backward run\n";

  (* Both algebraic part from y and the entire y' are computed by IDACalcIC. *)
  yB.{0} <- 0.0;
  yB.{1} <- 0.0;
  yB.{2} <- 0.50; (* not consistent *)

  (* Rough guess for ypB. *)
  ypB.{0} <- 0.80;
  ypB.{1} <- 0.75;
  ypB.{2} <- 0.0;

  (* Initialize qB *)
  qB.{0} <- 0.0;
  qB.{1} <- 0.0;
  qB.{2} <- 0.0;

  Adjoint.reinit indexB tb1 wyB wypB;

  (* Also reinitialize quadratures. *)
  AdjQuad.(init indexB (NoSens (rhsQB data)) wqB);

  (* Use IDACalcICB to compute consistent initial conditions
     for this backward problem. *)

  let wyyTB1 = wrap yyTB1
  and wypTB1 = wrap ypTB1
  in

  let id = RealArray.of_array [|VarId.differential;
                                VarId.differential;
                                VarId.algebraic|]
  in

  (* Specify which variables are differential (1) and which algebraic (0).*)
  (* Get the consistent IC found by IDAS. *)
  Adjoint.set_id indexB (wrap id);
  Adjoint.calc_ic indexB t1b wyyTB1 wypTB1 ~yb:wyyTB1 ~yb':wypTB1;

  print_string "Backward integration ... ";

  let _ = Adjoint.backward_normal ida_mem t0 in

  let nstB = Adjoint.get_num_steps indexB in

  printf "done ( nst = %d )\n" nstB;

  let _ = Adjoint.get indexB wyB wypB in

  let _ = AdjQuad.get indexB wqB in

  print_output tb1 yB ypB qB;

  print_string "Free memory\n\n"


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
