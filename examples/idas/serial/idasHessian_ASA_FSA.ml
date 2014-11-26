(*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:39 $
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
 *
 * Hessian using adjoint sensitivity example problem. 
 * 
 * This simple example problem for IDAS, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:
 *
 *   [ y1' + p1 * y1 - p2 * y2 * y3              = 0 
 *   [ y2' - p1 * y1 + p2 * y2 * y3 + p3 * y2^2  = 0 
 *   [ y1 + y2 + y3 -1                               = 0 
 * 
 *        [1]        [-p1]
 *   y(0)=[0]  y'(0)=[ p1]   p1 = 0.04   p2 = 1e4   p3 = 1e07   
 *        [0]        [ 0 ]
 *
 *       80
 *      / 
 *  G = | 0.5 * (y1^2 + y2^2 + y3^2) dt
 *      /
 *      0
 * Compute the gradient (using FSA and ASA) and Hessian (FSA over ASA)
 * of G with respect to parameters p1 and p2.
 *
 * Reference: D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * Error handling was suppressed for code readibility reasons.
*)
module RealArray = Sundials.RealArray
module Quad = Idas.Quadrature
module Sens = Idas.Sensitivity
module QuadSens = Sens.Quadrature
module Adjoint = Idas.Adjoint
module AdjQuad = Adjoint.Quadrature

let printf = Printf.printf

let nvconst = Nvector_serial.DataOps.n_vconst
let nvscale = Nvector_serial.DataOps.n_vscale

(* Problem Constants *)
let neq =      3             (* number of equations                  *)
let np =       2             (* number of sensitivities              *)

let t0 =       0.0   (* Initial time. *)
let tf =       80.0  (* Final time. *)

(* Tolerances *)
let rtol =     1e-08 (* scalar relative tolerance            *)
let atol =     1e-10 (* vector absolute tolerance components *)
let rtola =    1e-08 (* for adjoint integration              *)
let atola =    1e-08 (* for adjoint integration              *)

(* Parameters *)
let p1 = 0.04
let p2 = 1.0e4
let p3 = 3.0e7

(* User defined struct *)
type user_data = { p : RealArray.t;     (* size 3 *) }

(* residual for forward problem *)
let res data tres (yy : RealArray.t) (yp : RealArray.t) (rr : RealArray.t) =
  let y1  = yy.{0}
  and y2  = yy.{1}
  and y3  = yy.{2}

  and yp1 = yp.{0}
  and yp2 = yp.{1}

  and p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}
  in

  rr.{0} <- p1*.y1-.p2*.y2*.y3;
  rr.{1} <- -.rr.{0} +. p3*.y2*.y2 +. yp2;
  rr.{0} <- rr.{0} +. yp1;
  rr.{2} <- y1+.y2+.y3-.1.0

let resS : user_data -> RealArray.t Sens.sensresfn =
  fun data { Sens.y = yy; Sens.s = yyS; Sens.s' = ypS } resvalS ->
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}

  and y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}

  in

  for is = 0 to np-1 do
    let s1 = yyS.(is).{0}
    and s2 = yyS.(is).{1}
    and s3 = yyS.(is).{2}

    and sd1 = ypS.(is).{0}
    and sd2 = ypS.(is).{1}
    in

    let rs1 = sd1 +. p1*.s1 -. p2*.y3*.s2 -. p2*.y2*.s3
    and rs2 = sd2 -. p1*.s1 +. p2*.y3*.s2 +. p2*.y2*.s3 +. 2.0*.p3*.y2*.s2
    and rs3 = s1 +. s2 +. s3
    in

    let rs1,rs2 =
      match is with
       | 0 -> (rs1 +. y1,
               rs2 -. y1)
       | 1 -> (rs1 -. y2*.y3,
               rs2 +. y2*.y3)
       | _ -> (rs1, rs2)
    in

    resvalS.(is).{0} <- rs1;
    resvalS.(is).{1} <- rs2;
    resvalS.(is).{2} <- rs3
  done

let rhsQ data t (yy : RealArray.t) (yp : RealArray.t) (qdot : RealArray.t) =
  let y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}
  in
  qdot.{0} <- 0.5*.(y1*.y1+.y2*.y2+.y3*.y3)

let rhsQS : user_data -> RealArray.t QuadSens.quadsensrhsfn =
  fun data { QuadSens.y = yy; QuadSens.s = yyS } rhsQS ->
  let y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}
  in

  (* 1st sensitivity RHS *)
  let s1 = yyS.(0).{0}
  and s2 = yyS.(0).{1}
  and s3 = yyS.(0).{2}
  in
  rhsQS.(0).{0} <- y1*.s1 +. y2*.s2 +. y3*.s3;

  (* 2nd sensitivity RHS *)
  let s1 = yyS.(1).{0}
  and s2 = yyS.(1).{1}
  and s3 = yyS.(1).{2}
  in
  rhsQS.(1).{0} <- y1*.s1 +. y2*.s2 +. y3*.s3

(* Residuals for adjoint model. *)
let resBS1 : user_data -> RealArray.t Adjoint.bresfn_with_sens =
  fun data args yyS ypS rrBS ->
  let yy = args.Adjoint.y
  and yyB = args.Adjoint.yb
  and ypB = args.Adjoint.yb'
  in
  (* The parameters. *)
  (* Note: constants P1,P2,P3 from the original C source have names
     that clash with these local variables, but the constants are not
     used in this function.  *)
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}

  (* The y vector. *)
  and y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}

  (* The lambda vector. *)
  and l1 = yyB.{0}
  and l2 = yyB.{1}
  and l3 = yyB.{2}
  (* The mu vector. *)
  and m1 = yyB.{3}
  and m2 = yyB.{4}
  and m3 = yyB.{5}

  (* The lambda dot vector. *)
  and lp1 = ypB.{0}
  and lp2 = ypB.{1}
  (* The mu dot vector. *)
  and mp1 = ypB.{3}
  and mp2 = ypB.{4}

  (* The sensitivity with respect to p1 *)
  and s1 = yyS.(0).{0}
  and s2 = yyS.(0).{1}
  and s3 = yyS.(0).{2}
  in

  (* Temporary variables *)
  let l21 = l2-.l1 in

  rrBS.{0} <- lp1 +. p1*.l21 -. l3 +. y1;
  rrBS.{1} <- lp2 -. p2*.y3*.l21 -. 2.0*.p3*.y2*.l2 -. l3 +. y2;
  rrBS.{2} <- -.p2*.y2*.l21 -. l3 +. y3;

  rrBS.{3} <- mp1 +. p1*.(-.m1+.m2) -. m3 +. l21 +. s1;
  rrBS.{4} <- mp2 +. p2*.y3*.m1 -. (p2*.y3+.2.0*.p3*.y2)*.m2 -. m3
                +. p2*.s3*.l1 -. (2.0*.p3*.s2+.p2*.s3)*.l2 +. s2;
  rrBS.{5} <- p2*.y2*.(m1-.m2) -. m3 -. p2*.s2*.l21 +. s3

let rhsQBS1 : user_data -> RealArray.t AdjQuad.bquadrhsfn_with_sens =
  fun data { AdjQuad.y = yy; AdjQuad.yb = yyB } yyS ypS rhsBQS ->

  (* The y vector *)
  let y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}

  (* The lambda vector. *)
  and l1 = yyB.{0}
  and l2 = yyB.{1}
  (* The mu vector. *)
  and m1 = yyB.{3}
  and m2 = yyB.{4}

  (* The sensitivity with respect to p1 *)
  and s1 = yyS.(0).{0}
  and s2 = yyS.(0).{1}
  and s3 = yyS.(0).{2}
  in

  (* Temporary variables *)
  let l21 = l2-.l1 in

  rhsBQS.{0} <- -.y1*.l21;
  rhsBQS.{1} <- y2*.y3*.l21;

  rhsBQS.{2} <- y1*.(m1-.m2) -. s1*.l21;
  rhsBQS.{3} <- y2*.y3*.(m2-.m1) +. (y3*.s2+.y2*.s3)*.l21

let resBS2 : user_data -> RealArray.t Adjoint.bresfn_with_sens =
  fun data { Adjoint.y = yy; Adjoint.yb = yyB; Adjoint.yb' = ypB } yyS ypS rrBS
  ->
  (* The parameters. *)
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2};

  (* The y vector. *)
  and y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}

  (* The lambda vector. *)
  and l1 = yyB.{0}
  and l2 = yyB.{1}
  and l3 = yyB.{2}
  (* The mu vector. *)
  and m1 = yyB.{3}
  and m2 = yyB.{4}
  and m3 = yyB.{5}

  (* The lambda dot vector. *)
  and lp1 = ypB.{0}
  and lp2 = ypB.{1}
  (* The mu dot vector. *)
  and mp1 = ypB.{3}
  and mp2 = ypB.{4}

  (* The sensitivity with respect to p2 *)
  and s1 = yyS.(1).{0}
  and s2 = yyS.(1).{1}
  and s3 = yyS.(1).{2}
  in

  (* Temporary variables *)
  let l21 = l2-.l1 in

  rrBS.{0} <- lp1 +. p1*.l21 -. l3 +. y1;
  rrBS.{1} <- lp2 -. p2*.y3*.l21 -. 2.0*.p3*.y2*.l2 -. l3 +. y2;
  rrBS.{2} <- -.p2*.y2*.l21 -. l3 +. y3;

  rrBS.{3} <- mp1 +. p1*.(-.m1+.m2) -. m3 +. s1;
  rrBS.{4} <- mp2 +. p2*.y3*.m1 -. (p2*.y3+.2.0*.p3*.y2)*.m2 -. m3 +. (y3+.p2*.s3)*.l1 -. (y3+.2.0*.p3*.s2+.p2*.s3)*.l2 +. s2;
  rrBS.{5} <- p2*.y2*.(m1-.m2) -. m3 -. (y2+.p2*.s2)*.l21 +. s3

let rhsQBS2 : user_data -> RealArray.t AdjQuad.bquadrhsfn_with_sens =
  fun data {AdjQuad.y = yy; AdjQuad.yb = yyB } yyS ypS rhsBQS ->
  (* The y vector *)
  let y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}

  (* The lambda vector. *)
  and l1 = yyB.{0}
  and l2 = yyB.{1}
  (* The mu vector. *)
  and m1 = yyB.{3}
  and m2 = yyB.{4}

  (* The sensitivity with respect to p2 *)
  and s1 = yyS.(1).{0}
  and s2 = yyS.(1).{1}
  and s3 = yyS.(1).{2}
  in

  (* Temporary variables *)
  let l21 = l2-.l1 in

  rhsBQS.{0} <- -.y1*.l21;
  rhsBQS.{1} <-  y2*.y3*.l21;

  rhsBQS.{2} <- y1*.(m1-.m2) -. s1*.l21;
  rhsBQS.{3} <- y2*.y3*.(m2-.m1) +. (y3*.s2+.y2*.s3)*.l21

let main () =
  (* Print problem description *)
  printf("\nAdjoint Sensitivity Example for Chemical Kinetics\n");
  printf("---------------------------------------------------------\n");
  printf("DAE: dy1/dt + p1*y1 - p2*y2*y3 = 0\n");
  printf("     dy2/dt - p1*y1 + p2*y2*y3 + p3*(y2)^2 = 0\n");
  printf("               y1  +  y2  +  y3 = 0\n\n");
  printf("Find dG/dp and d^2G/dp^2, where p=[p1,p2] for\n");
  printf("     G = int_t0^tB0 g(t,p,y) dt\n");
  printf("     g(t,p,y) = y3\n\n\n");

  (* Alocate and initialize user data. *)
  let data = { p = RealArray.of_array [|p1;p2;p3|] } in

  (* Consistent IC *)
  let yy = RealArray.create neq
  and yp = RealArray.create neq
  in
  yy.{0} <- 1.0;
  yy.{1} <- 0.0;
  yy.{2} <- 0.0;
  yp.{0} <- -. p1;
  yp.{1} <- p1;
  yp.{2} <- 0.0;

  let q = RealArray.create 1 in
  nvconst 0.0 q;

  (* Wrap arrays in nvectors.  *)
  let wyy = Nvector_serial.wrap yy
  and wyp = Nvector_serial.wrap yp
  and wq  = Nvector_serial.wrap q
  in

  let yyS = Array.init np (fun _ -> Nvector_serial.wrap (RealArray.copy yy))
  and ypS = Array.init np (fun _ -> Nvector_serial.wrap (RealArray.copy yp))
  in
  nvconst 0.0 (Nvector.unwrap yyS.(0));
  nvconst 0.0 (Nvector.unwrap yyS.(1));
  nvconst 0.0 (Nvector.unwrap ypS.(0));
  nvconst 0.0 (Nvector.unwrap ypS.(1));

  let qS = Array.init np (fun _ -> Nvector_serial.wrap (RealArray.copy q)) in
  nvconst 0.0 (Nvector.unwrap qS.(0));

  (* Forward problem's setup. *)
  let ti = t0 in
  let ida_mem =
    Ida.init (Ida.Dls.dense ())
      (Ida.SStolerances (rtol,atol))
      (res data)
      ti
      wyy wyp
  in
  Ida.set_max_num_steps ida_mem 1500;

  (* Quadrature's setup. *)
  Quad.init ida_mem (rhsQ data) wq;

  Quad.set_tolerances ida_mem (Quad.SStolerances (rtol,atol));

  (* Sensitivity's setup. *)
  Sens.init ida_mem Sens.EEtolerances Sens.Simultaneous
    ~fs:(resS data) yyS ypS;
  Sens.set_err_con ida_mem true;

  (* Setup of quadrature's sensitivities *)
  QuadSens.init ida_mem ~fqs:(rhsQS data) qS;
  QuadSens.set_tolerances ida_mem QuadSens.EEtolerances;

  (* Initialize ASA. *)
  Adjoint.init ida_mem 100 Adjoint.IHermite;

  printf "---------------------------------------------------------\n";
  printf "Forward integration\n";
  printf "---------------------------------------------------------\n\n";

  let _ = Adjoint.forward_normal ida_mem tf wyy wyp in
  let _ = Quad.get ida_mem wq in
  let g = q.{0} in
  printf "     G:    %12.4e\n" g;

  (* Sensitivities are needed for IC of backward problems. *)
  Sens.get_dky ida_mem yyS tf 0;
  Sens.get_dky ida_mem ypS tf 1;

  let _ = QuadSens.get ida_mem qS in
  printf "   dG/dp:  %12.4e %12.4e\n"
    (Nvector.unwrap qS.(0)).{0}
    (Nvector.unwrap qS.(1)).{0};
  printf "\n";
  (******************************
  * BACKWARD PROBLEM #1
  *******************************)

  (* Consistent IC. *)
  let yyB1 = RealArray.create (2*neq)
  and ypB1 = RealArray.create (2*neq)
  in

  nvconst 0.0 yyB1;
  yyB1.{2} <- yy.{2};
  yyB1.{5} <- (Nvector.unwrap yyS.(0)).{2};

  nvconst 0.0 ypB1;
  ypB1.{0} <- yy.{2} -. yy.{0};
  ypB1.{1} <- yy.{2} -. yy.{1};
  ypB1.{3} <- (Nvector.unwrap yyS.(0)).{2} -. (Nvector.unwrap yyS.(0)).{0};
  ypB1.{4} <- (Nvector.unwrap yyS.(0)).{2} -. (Nvector.unwrap yyS.(0)).{1};

  let qB1 = RealArray.create (2*np) in
  nvconst 0.0 qB1;

  let wyyB1 = Nvector_serial.wrap yyB1
  and wypB1 = Nvector_serial.wrap ypB1
  and wqB1  = Nvector_serial.wrap qB1
  in

  let indexB1 =
    Adjoint.init_backward ida_mem (Adjoint.Dls.dense ())
      (Adjoint.SStolerances (rtola, atola))
      (Adjoint.WithSens (resBS1 data))
      tf wyyB1 wypB1
  in
  Adjoint.set_max_num_steps indexB1 5000;
  AdjQuad.init indexB1 (AdjQuad.WithSens (rhsQBS1 data)) wqB1;

  (******************************
  * BACKWARD PROBLEM #2  
  *******************************)

  (* Consistent IC. *)
  let yyB2 = RealArray.create (2*neq)
  and ypB2 = RealArray.create (2*neq)
  in

  nvconst 0.0 yyB2;
  yyB2.{2} <- yy.{2};
  yyB2.{5} <- (Nvector.unwrap yyS.(1)).{2};

  nvconst 0.0 ypB2;
  ypB2.{0} <- yy.{2}-.yy.{0};
  ypB2.{1} <- yy.{2}-.yy.{1};
  ypB2.{3} <- (Nvector.unwrap yyS.(1)).{2} -. (Nvector.unwrap yyS.(1)).{0};
  ypB2.{4} <- (Nvector.unwrap yyS.(1)).{2} -. (Nvector.unwrap yyS.(1)).{1};

  let qB2 = RealArray.create (2*np) in
  nvconst 0.0 qB2;

  let wyyB2 = Nvector_serial.wrap yyB2
  and wypB2 = Nvector_serial.wrap ypB2
  and wqB2  = Nvector_serial.wrap qB2
  in

  let indexB2 =
    Adjoint.init_backward ida_mem (Adjoint.Dls.dense ())
      (Adjoint.SStolerances (rtola, atola))
      (Adjoint.WithSens (resBS2 data))
      tf
      wyyB2
      wypB2
  in
  Adjoint.set_max_num_steps indexB2 2500;
  AdjQuad.init indexB2 (AdjQuad.WithSens (rhsQBS2 data)) wqB2;

  (* Integrate backward problems. *)
  printf "---------------------------------------------------------\n";
  printf "Backward integration \n";
  printf "---------------------------------------------------------\n\n";

  Adjoint.backward_normal ida_mem ti;
  let _ = Adjoint.get indexB1 wyyB1 wypB1 in

  let _ = AdjQuad.get indexB1 wqB1 in
  let _ = AdjQuad.get indexB2 wqB2 in
  printf "   dG/dp:  %12.4e %12.4e   (from backward pb. 1)\n" qB1.{0} qB1.{1};
  printf "   dG/dp:  %12.4e %12.4e   (from backward pb. 2)\n" qB2.{0} qB2.{1};

  printf "\n";
  printf "   H = d2G/dp2:\n";
  printf "        (1)            (2)\n";
  printf "  %12.4e  %12.4e\n" qB1.{2} qB2.{2};
  printf "  %12.4e  %12.4e\n" qB1.{3} qB2.{3};

  (*********************************
  * Use Finite Differences to verify
  **********************************)

  (* Perturbations are of different magnitudes as p1 and p2 are. *)
  let dp1 = 1.0e-3
  and dp2 = 2.5e+2
  in

  printf "\n";
  printf "---------------------------------------------------------\n";
  printf "Finite Differences ( dp1=%6.1e and dp2 = %6.1e )\n" dp1 dp2;
  printf "---------------------------------------------------------\n\n";



  (********************
  * Forward FD for p1
  ********************)
  let rtolFD = 1.0e-12
  and atolFD = 1.0e-14
  in

  let ti = t0 in

  data.p.{0} <- data.p.{0} +. dp1;
  yy.{0} <- 1.0;
  yy.{1} <- 0.0;
  yy.{2} <- 0.0;
  yp.{0} <- -.data.p.{0};
  yp.{1} <- -.yp.{0};
  yp.{2} <- 0.0;
  nvconst 0.0 q;

  let ida_mem =
    Ida.init (Ida.Dls.dense ())
      (Ida.SStolerances (rtolFD, atolFD))
      (res data)
      ti wyy wyp
  in
  Ida.set_max_num_steps ida_mem 10000;

  Quad.init ida_mem (rhsQ data) wq;
  Quad.set_tolerances ida_mem (Quad.SStolerances (rtolFD,atolFD));

  let _ = Ida.solve_normal ida_mem tf wyy wyp in
  let _ = Quad.get ida_mem wq in
  let gp = q.{0} in

  (********************
  * Backward FD for p1
  ********************)
  data.p.{0} <- data.p.{0} -. 2.0*.dp1;

  yy.{0} <- 1.0; yy.{1} <- 0.0; yy.{2} <- 0.0;
  yp.{0} <- -.data.p.{0}; yp.{1} <- -.yp.{0}; yp.{2} <- 0.0;
  nvconst 0.0 q;

  Ida.reinit ida_mem ti wyy wyp;
  Quad.reinit ida_mem wq;

  let _ = Ida.solve_normal ida_mem tf wyy wyp in
  let _ = Quad.get ida_mem wq in
  let gm = q.{0} in

  (* Compute FD for p1. *)
  let grdG_fwd = RealArray.create 2
  and grdG_bck = RealArray.create 2
  and grdG_cntr = RealArray.create 2
  in
  grdG_fwd.{0} <- (gp-.g)/.dp1;
  grdG_bck.{0} <- (g-.gm)/.dp1;
  grdG_cntr.{0} <- (gp-.gm)/.(2.0*.dp1);
  let h11 = (gp -. 2.0*.g +. gm) /. (dp1*.dp1) in

  (********************
  * Forward FD for p2
  ********************)
  (*restore p1*)
  data.p.{0} <- data.p.{0} +. dp1; 
  data.p.{1} <- data.p.{1} +. dp2;

  yy.{0} <- 1.0; yy.{1} <- 0.0; yy.{2} <- 0.0;
  yp.{0} <- -.data.p.{0}; yp.{1} <- -.yp.{0}; yp.{2} <- 0.0;
  nvconst 0.0 q;

  Ida.reinit ida_mem ti wyy wyp;
  Quad.reinit ida_mem wq;

  let _ = Ida.solve_normal ida_mem tf wyy wyp in
  let _ = Quad.get ida_mem wq in
  let gp = q.{0} in

  (********************
  * Backward FD for p2
  ********************)
  data.p.{1} <- data.p.{1} -. 2.0*.dp2;

  yy.{0} <- 1.0; yy.{1} <- 0.0; yy.{2} <- 0.0;
  yp.{0} <- -.data.p.{0}; yp.{1} <- -.yp.{0}; yp.{2} <- 0.0;
  nvconst 0.0 q;

  Ida.reinit ida_mem ti wyy wyp;
  Quad.reinit ida_mem wq;

  let _ = Ida.solve_normal ida_mem tf wyy wyp in
  let _ = Quad.get ida_mem wq in
  let gm = q.{0} in

  (* Compute FD for p2. *)
  grdG_fwd.{1} <- (gp-.g)/.dp2;
  grdG_bck.{1} <- (g-.gm)/.dp2;
  grdG_cntr.{1} <- (gp-.gm)/.(2.0*.dp2);
  let h22 = (gp -. 2.0*.g +. gm) /. (dp2*.dp2) in

  printf "\n";
  printf "   dG/dp:  %12.4e  %12.4e   (fwd FD)\n"  grdG_fwd.{0}  grdG_fwd.{1};
  printf "           %12.4e  %12.4e   (bck FD)\n"  grdG_bck.{0}  grdG_bck.{1};
  printf "           %12.4e  %12.4e   (cntr FD)\n" grdG_cntr.{0} grdG_cntr.{1};
  printf "\n";
  printf "  H(1,1):  %12.4e\n" h11;
  printf "  H(2,2):  %12.4e\n" h22


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
