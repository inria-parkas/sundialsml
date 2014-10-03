(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2010/12/01 22:57:59 $
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
 *
 * Hessian through adjoint sensitivity example problem.
 *
 *        [ - p1 * y1^2 - y3 ]           [ 1 ]
 *   y' = [    - y2          ]    y(0) = [ 1 ]
 *        [ -p2^2 * y2 * y3  ]           [ 1 ]
 *
 *   p1 = 1.0
 *   p2 = 2.0
 *
 *           2
 *          /
 *   G(p) = |  0.5 * ( y1^2 + y2^2 + y3^2 ) dt
 *          /
 *          0
 *
 * Compute the gradient (ASA) and Hessian (FSA over ASA) of G(p).
 *
 * See D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module Quad = Cvodes.Quadrature
module Sens = Cvodes.Sensitivity
module QuadSens = Cvodes.Sensitivity.Quadrature
module Adj = Cvodes.Adjoint
module QuadAdj = Cvodes.Adjoint.Quadrature
open Bigarray
let unwrap = Nvector.unwrap

let printf = Printf.printf

let ith (v : RealArray.t) i = v.{i - 1}
let set_ith (v : RealArray.t) i e = v.{i - 1} <- e

let zero = 0.0
let one  = 1.0

type user_data = {
  mutable p1 : float;
  mutable p2 : float;
}

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 *)

let f data t (y : RealArray.t) (ydot : RealArray.t) =
  let p1 = data.p1
  and p2 = data.p2 in 
  ydot.{0} <- -.p1*.y.{0}*.y.{0} -. y.{2};
  ydot.{1} <- -.y.{1};
  ydot.{2} <- -.p2*.p2*.y.{1}*.y.{2}

let fQ data t (y : RealArray.t) (qdot : RealArray.t) =
  qdot.{0} <- 0.5 *. ( y.{0}*.y.{0} +. y.{1}*.y.{1} +. y.{2}*.y.{2})

let fS data t (y : RealArray.t) (ydot : RealArray.t)
              (yS : RealArray.t array) (ySdot : RealArray.t array) tmp1 tmp2 =
  let p1 = data.p1
  and p2 = data.p2 in
  (* 1st sensitivity RHS *)
  let s1 = yS.(0).{0} in
  let s2 = yS.(0).{1} in
  let s3 = yS.(0).{2} in
  let fys1 = -. 2.0*.p1*.y.{0} *. s1 -. s3
  and fys2 = -. s2
  and fys3 = -. p2*.p2*.y.{2} *. s2 -. p2*.p2*.y.{1} *. s3
  in
  ySdot.(0).{0} <- (fys1 -. y.{0}*.y.{0});
  ySdot.(0).{1} <- fys2;
  ySdot.(0).{2} <- fys3;

  (* 2nd sensitivity RHS *)
  let s1 = yS.(1).{0}
  and s2 = yS.(1).{1}
  and s3 = yS.(1).{2}
  in
  let fys1 = -. 2.0*.p1*.y.{0} *. s1 -. s3
  and fys2 = -. s2
  and fys3 = -. p2*.p2*.y.{2} *. s2 -. p2*.p2*.y.{1} *. s3
  in
  ySdot.(1).{0} <- fys1;
  ySdot.(1).{1} <- fys2;
  ySdot.(1).{2} <- (fys3 -. 2.0*.p2*.y.{1}*.y.{2})

let fQS data t (y : RealArray.t) (yS : RealArray.t array)
               yQdot (yQSdot : RealArray.t array) tmp tmpQ =
  (* 1st sensitivity RHS *)
  let s1 = yS.(0).{0}
  and s2 = yS.(0).{1}
  and s3 = yS.(0).{2}
  in
  yQSdot.(0).{0} <- y.{0}*.s1 +. y.{1}*.s2 +. y.{2}*.s3;

  (* 1st sensitivity RHS *)
  let s1 = yS.(1).{0}
  and s2 = yS.(1).{1}
  and s3 = yS.(1).{2}
  in
  yQSdot.(1).{0} <- y.{0}*.s1 +. y.{1}*.s2 +. y.{2}*.s3

let fB1 data t (y : RealArray.t) (yS : RealArray.t array)
               (yB : RealArray.t) (yBdot : RealArray.t) =
  let p1 = data.p1 
  and p2 = data.p2
  in
  let s1 = yS.(0).{0} (* sensitivity 1 *)
  and s2 = yS.(0).{1}
  and s3 = yS.(0).{2}
  in
  let l1 = yB.{0}     (* lambda *)
  and l2 = yB.{1}
  and l3 = yB.{2}
  in
  let m1 = yB.{3}     (* mu *)
  and m2 = yB.{4}
  and m3 = yB.{5}
  in
  yBdot.{0} <- 2.0*.p1*.y.{0} *. l1      -. y.{0};
  yBdot.{1} <- l2 +. p2*.p2*.y.{2} *. l3 -. y.{1};
  yBdot.{2} <- l1 +. p2*.p2*.y.{1} *. l3 -. y.{2};

  yBdot.{3} <- 2.0*.p1*.y.{0} *. m1      +. l1 *. 2.0*.(y.{0} +. p1*.s1) -. s1;
  yBdot.{4} <- m2 +. p2*.p2*.y.{2} *. m3 +. l3 *. p2*.p2*.s3             -. s2;
  yBdot.{5} <- m1 +. p2*.p2*.y.{1} *. m3 +. l3 *. p2*.p2*.s2             -. s3

let fQB1 data t (y : RealArray.t) (yS : RealArray.t array)
                                  (yB : RealArray.t) (qBdot : RealArray.t) =
  let p2 = data.p2
  in
  let s1 = yS.(0).{0} (* sensitivity 1 *)
  and s2 = yS.(0).{1}
  and s3 = yS.(0).{2}
  in
  let l1 = yB.{0}     (* lambda *)
  and l3 = yB.{2}
  in
  let m1 = yB.{3}     (* mu *)
  and m3 = yB.{5}
  in
  qBdot.{0} <- -.y.{0}*.y.{0} *. l1;
  qBdot.{1} <- -.2.0*.p2*.y.{1}*.y.{2} *. l3;

  qBdot.{2} <- -.y.{0}*.y.{0} *. m1          -. l1 *. 2.0*.y.{0}*.s1;
  qBdot.{3} <- -.2.0*.p2*.y.{1}*.y.{2} *. m3 -.
                                    l3 *. 2.0*.(p2*.y.{2}*.s2 +. p2*.y.{1}*.s3)

let fB2 data t (y : RealArray.t) (yS : RealArray.t array)
               (yB : RealArray.t) (yBdot : RealArray.t) =
  let p1 = data.p1 
  and p2 = data.p2
  in
  let s1 = yS.(1).{0} (* sensitivity 2 *)
  and s2 = yS.(1).{1}
  and s3 = yS.(1).{2}
  in
  let l1 = yB.{0}     (* lambda *)
  and l2 = yB.{1}
  and l3 = yB.{2}
  in
  let m1 = yB.{3}     (* mu *)
  and m2 = yB.{4}
  and m3 = yB.{5}
  in
  yBdot.{0} <- 2.0*.p1*.y.{0} *. l1      -. y.{0};
  yBdot.{1} <- l2 +. p2*.p2*.y.{2} *. l3 -. y.{1};
  yBdot.{2} <- l1 +. p2*.p2*.y.{1} *. l3 -. y.{2};

  yBdot.{3} <- 2.0*.p1*.y.{0} *. m1      +.
                                     l1 *. 2.0*.p1*.s1                    -. s1;
  yBdot.{4} <- m2 +. p2*.p2*.y.{2} *. m3 +.
                                     l3 *. (2.0*.p2*.y.{2} +. p2*.p2*.s3) -. s2;
  yBdot.{5} <- m1 +. p2*.p2*.y.{1} *. m3 +.
                                     l3 *. (2.0*.p2*.y.{2} +. p2*.p2*.s2) -. s3

let fQB2 data t (y : RealArray.t) (yS : RealArray.t array)
                (yB : RealArray.t) (qBdot : RealArray.t) =
  let p2 = data.p2
  in
  let s1 = yS.(1).{0} (* sensitivity 2 *)
  and s2 = yS.(1).{1}
  and s3 = yS.(1).{2}
  in
  let l1 = yB.{0}     (* lambda *)
  and l3 = yB.{2}
  in
  let m1 = yB.{3}     (* mu *)
  and m3 = yB.{5}
  in
  qBdot.{0} <- -.y.{0}*.y.{0} *. l1;
  qBdot.{1} <- -.2.0*.p2*.y.{1}*.y.{2} *. l3;

  qBdot.{2} <- -.y.{0}*.y.{0} *. m1          -. l1 *. 2.0*.y.{0}*.s1;
  qBdot.{3} <- -.2.0*.p2*.y.{1}*.y.{2} *. m3 -.
                    l3 *. 2.0*.(p2*.y.{2}*.s2 +. p2*.y.{1}*.s3 +. y.{1}*.y.{2})

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

let print_fwd_stats cvode_mem =
  let { Cvode.num_steps = nst;
        Cvode.num_rhs_evals = nfe;
        Cvode.num_lin_solv_setups = nsetups;
        Cvode.num_err_test_fails = netf;
        Cvode.last_order = qlast;
        Cvode.current_order = qcur;
        Cvode.actual_init_step = h0u;
        Cvode.last_step = hlast;
        Cvode.current_step = hcur;
        Cvode.current_time = tcur;
    } = Cvode.get_integrator_stats cvode_mem
  in
  let nni, ncfn = Cvode.get_nonlin_solv_stats cvode_mem
  and nfQe, netfQ = Quad.get_stats cvode_mem
  and { Sens.num_rhs_evals = nfSe;
        Sens.num_sens_evals = nfeS;
        Sens.num_err_test_fails = netfS;
        Sens.num_lin_solv_setups = nsetupsS;
      } = Sens.get_stats cvode_mem
  and nniS, ncfnS = Sens.get_nonlin_solv_stats cvode_mem
  and nfQSe, netfQS = QuadSens.get_stats cvode_mem
  in
  printf " Number steps: %5d\n\n" nst;
  printf " Function evaluations:\n";
  printf "  f:        %5d\n  fQ:       %5d\n  fS:       %5d\n  fQS:      %5d\n" 
         nfe nfQe nfSe nfQSe;
  printf " Error test failures:\n";
  printf "  netf:     %5d\n  netfQ:    %5d\n  netfS:    %5d\n  netfQS:   %5d\n" 
         netf netfQ netfS netfQS;
  printf " Linear solver setups:\n";
  printf "  nsetups:  %5d\n  nsetupsS: %5d\n" nsetups nsetupsS;
  printf " Nonlinear iterations:\n";
  printf "  nni:      %5d\n  nniS:     %5d\n" nni nniS;
  printf " Convergence failures:\n";
  printf "  ncfn:     %5d\n  ncfnS:    %5d\n" ncfn ncfnS;
  printf "\n"

let print_bck_stats cvode_mem =
  let { Cvode.num_steps = nst;
        Cvode.num_rhs_evals = nfe;
        Cvode.num_lin_solv_setups = nsetups;
        Cvode.num_err_test_fails = netf;
        Cvode.last_order = qlast;
        Cvode.current_order = qcur;
        Cvode.actual_init_step = h0u;
        Cvode.last_step = hlast;
        Cvode.current_step = hcur;
        Cvode.current_time = tcur;
    } = Adj.get_integrator_stats cvode_mem
  in
  let nni, ncfn = Adj.get_nonlin_solv_stats cvode_mem
  and nfQe, netfQ = QuadAdj.get_stats cvode_mem
  in
  printf " Number steps: %5d\n\n" nst;
  printf " Function evaluations:\n";
  printf "  f:        %5d\n  fQ:       %5d\n" nfe nfQe;
  printf " Error test failures:\n";
  printf "  netf:     %5d\n  netfQ:    %5d\n" netf netfQ;
  printf " Linear solver setups:\n";
  printf "  nsetups:  %5d\n" nsetups;
  printf " Nonlinear iterations:\n";
  printf "  nni:      %5d\n" nni;
  printf " Convergence failures:\n";
  printf "  ncfn:     %5d\n" ncfn;
  printf "\n"

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  let grdG_fwd  = Array.make 2 0.0 in
  let grdG_bck  = Array.make 2 0.0 in
  let grdG_cntr = Array.make 2 0.0 in

  (* User data structure *)
  let data = { p1 = 1.0; p2 = 2.0 } in

  (* Problem size, integration interval, and tolerances *)
  let neq = 3 in
  let np  = 2 in
  let np2 = 2*np in

  let t0  = 0.0 in
  let tf  = 2.0 in

  let reltol   = 1.0e-8 in

  let abstol   = 1.0e-8 in
  let abstolQ  = 1.0e-8 in

  let abstolB  = 1.0e-8 in
  let abstolQB = 1.0e-8 in

  (* Initializations for forward problem *)
  let y = Nvector_serial.make neq one in
  let ydata = unwrap y in
  let yQ = Nvector_serial.make 1 zero in
  let yS = Array.init np (fun _ -> Nvector_serial.make neq zero) in
  let yQS = Array.init np (fun _ -> Nvector_serial.make 1 zero) in

  (* Create and initialize forward problem *)
  let cvode_mem =
    Cvode.init
        Cvode.BDF
        (Cvode.Newton (Cvode.Dls.dense None))
        (Cvode.SStolerances (reltol, abstol))
        (f data)
        t0
        y
  in
  Quad.init cvode_mem (fQ data) yQ;
  Quad.set_tolerances cvode_mem (Quad.SStolerances (reltol, abstolQ));

  Sens.init cvode_mem Sens.EEtolerances Sens.Simultaneous
                      Sens.no_sens_params (Sens.AllAtOnce (Some (fS data))) yS;
  Sens.set_err_con cvode_mem true;

  QuadSens.init cvode_mem (Some (fQS data)) yQS;
  QuadSens.set_tolerances cvode_mem (QuadSens.EEtolerances);

  (* Initialize ASA *)
  let steps = 100 in
  Adj.init cvode_mem steps Adj.IPolynomial;

  (* Forward integration *)
  printf "-------------------\n";
  printf "Forward integration\n";
  printf "-------------------\n\n";
  let time, ncheck, _ = Adj.forward_normal cvode_mem tf y in
  ignore (Quad.get cvode_mem yQ);
  let g = ith (unwrap yQ) 1 in

  ignore (Sens.get cvode_mem yS);
  ignore (QuadSens.get cvode_mem yQS);

  printf "ncheck = %d\n"  ncheck;
  printf "\n";
  printf "     y:    %12.4e %12.4e %12.4e"
                    (ith ydata 1) (ith ydata 2) (ith ydata 3);
  printf "     G:    %12.4e\n"  (ith (unwrap yQ) 1);
  printf "\n";
  printf "     yS1:  %12.4e %12.4e %12.4e\n"
                                (ith (unwrap yS.(0)) 1)
                                (ith (unwrap yS.(0)) 2)
                                (ith (unwrap yS.(0)) 3);
  printf "     yS2:  %12.4e %12.4e %12.4e\n"
                                (ith (unwrap yS.(1)) 1)
                                (ith (unwrap yS.(1)) 2)
                                (ith (unwrap yS.(1)) 3);
  printf "\n";
  printf "   dG/dp:  %12.4e %12.4e\n" (ith (unwrap yQS.(0)) 1)
                                      (ith (unwrap yQS.(1)) 1);
  printf "\n";

  printf "Final Statistics for forward pb.\n";
  printf "--------------------------------\n";
  print_fwd_stats cvode_mem;

  (* Initializations for backward problems *)
  let yB1  = Nvector_serial.make (2 * neq) zero in
  let yQB1 = Nvector_serial.make np2 zero in
  let yB2  = Nvector_serial.make (2 * neq) zero in
  let yQB2 = Nvector_serial.make np2 zero in

  (* Create and initialize backward problems (one for each column of the Hessian) *)
  let cvode_memB1 =
    Adj.init_backward cvode_mem Cvode.BDF
                                (Adj.Newton (Adj.Dls.dense None))
                                (Adj.SStolerances (reltol, abstolB))
                                (Adj.WithSens (fB1 data))
                                tf yB1
  in
  QuadAdj.init cvode_memB1 (QuadAdj.WithSens (fQB1 data)) yQB1;
  QuadAdj.set_tolerances cvode_memB1 (QuadAdj.SStolerances (reltol, abstolQB));

  let cvode_memB2 =
    Adj.init_backward cvode_mem Cvode.BDF
                                (Adj.Newton (Adj.Dls.dense None))
                                (Adj.SStolerances (reltol, abstolB))
                                (Adj.WithSens (fB2 data))
                                tf yB2
  in
  QuadAdj.init cvode_memB2 (QuadAdj.WithSens (fQB2 data)) yQB2;
  QuadAdj.set_tolerances cvode_memB2 (QuadAdj.SStolerances (reltol, abstolB));

  (* Backward integration *)
  printf "---------------------------------------------\n";
  printf "Backward integration ... (2 adjoint problems)\n";
  printf "---------------------------------------------\n\n";

  Adj.backward_normal cvode_mem t0;

  ignore (Adj.get cvode_memB1 yB1);
  ignore (QuadAdj.get cvode_memB1 yQB1);

  ignore (Adj.get cvode_memB2 yB2);
  ignore (QuadAdj.get cvode_memB2 yQB2);

  printf "   dG/dp:  %12.4e %12.4e   (from backward pb. 1)\n"
                                                  (-.ith (unwrap yQB1) 1)
                                                  (-.ith (unwrap yQB1) 2);
  printf "           %12.4e %12.4e   (from backward pb. 2)\n"
                                                 (-.ith (unwrap yQB2) 1)
                                                 (-.ith (unwrap yQB2) 2);
  printf "\n";
  printf "   H = d2G/dp2:\n";
  printf "        (1)            (2)\n";
  printf "  %12.4e   %12.4e\n" (-.ith (unwrap yQB1) 3) (-.ith (unwrap yQB2) 3);
  printf "  %12.4e   %12.4e\n" (-.ith (unwrap yQB1) 4) (-.ith (unwrap yQB2) 4);
  printf "\n";

  printf "Final Statistics for backward pb. 1\n";
  printf "-----------------------------------\n";
  print_bck_stats cvode_memB1;

  printf "Final Statistics for backward pb. 2\n";
  printf "-----------------------------------\n";
  print_bck_stats cvode_memB2;

  (* Free CVODES memory *)

  (* Finite difference tests *)
  let dp = 1.0e-2 in

  printf "-----------------------\n";
  printf "Finite Difference tests\n";
  printf "-----------------------\n\n";

  printf "del_p = %g\n\n" dp;

  RealArray.fill (unwrap y) one;
  let cvode_mem =
    Cvode.init
        Cvode.BDF
        (Cvode.Newton (Cvode.Dls.dense None))
        (Cvode.SStolerances (reltol, abstol))
        (f data)
        t0
        y
  in
  RealArray.fill (unwrap yQ) zero;
  Quad.init cvode_mem (fQ data) yQ;
  Quad.set_tolerances cvode_mem (Quad.SStolerances (reltol, abstolQ));

  data.p1 <- data.p1 +. dp;

  ignore (Cvode.solve_normal cvode_mem tf y);
  ignore (Quad.get cvode_mem yQ);

  let gp = ith (unwrap yQ) 1 in
  printf "p1+  y:   %12.4e %12.4e %12.4e" (ith (unwrap y) 1)
                                          (ith (unwrap y) 2)
                                          (ith (unwrap y) 3);
  printf "     G:   %12.4e\n" (ith (unwrap yQ) 1);

  data.p1 <- data.p1 -. 2.0*.dp;

  RealArray.fill (unwrap y) one;
  RealArray.fill (unwrap yQ) zero;

  Cvode.reinit cvode_mem t0 y;
  Quad.reinit cvode_mem yQ;

  ignore (Cvode.solve_normal cvode_mem tf y);
  ignore (Quad.get cvode_mem yQ);

  let gm = ith (unwrap yQ) 1 in
  printf "p1-  y:   %12.4e %12.4e %12.4e" (ith (unwrap y) 1)
                                          (ith (unwrap y) 2)
                                          (ith (unwrap y) 3);
  printf "     G:   %12.4e\n" (ith (unwrap yQ) 1);
 
  data.p1 <- data.p1 +. dp;

  grdG_fwd.(0)  <- (gp-.g)/.dp;
  grdG_bck.(0)  <- (g-.gm)/.dp;
  grdG_cntr.(0) <- (gp-.gm)/.(2.0*.dp);
  let h11 = (gp -. 2.0*.g +. gm) /. (dp*.dp) in

  data.p2 <- data.p2 +. dp;

  RealArray.fill (unwrap y) one;
  RealArray.fill (unwrap yQ) zero;

  Cvode.reinit cvode_mem t0 y;
  Quad.reinit cvode_mem yQ;

  ignore (Cvode.solve_normal cvode_mem tf y);
  ignore (Quad.get cvode_mem yQ);

  let gp = ith (unwrap yQ) 1 in
  printf "p2+  y:   %12.4e %12.4e %12.4e" (ith (unwrap y) 1)
                                          (ith (unwrap y) 2)
                                          (ith (unwrap y) 3);
  printf "     G:   %12.4e\n" (ith (unwrap yQ) 1);
 
  data.p2 <- data.p2 -. 2.0*.dp;

  RealArray.fill (unwrap y) one;
  RealArray.fill (unwrap yQ) zero;

  Cvode.reinit cvode_mem t0 y;
  Quad.reinit cvode_mem yQ;

  ignore (Cvode.solve_normal cvode_mem tf y);
  ignore (Quad.get cvode_mem yQ);

  let gm = ith (unwrap yQ) 1 in
  printf "p2-  y:   %12.4e %12.4e %12.4e" (ith (unwrap y) 1)
                                          (ith (unwrap y) 2)
                                          (ith (unwrap y) 3);
  printf "     G:   %12.4e\n" (ith (unwrap yQ) 1);

  data.p2 <- data.p2 +. dp;

  grdG_fwd.(1)  <- (gp-.g)/.dp;
  grdG_bck.(1)  <- (g-.gm)/.dp;
  grdG_cntr.(1) <- (gp-.gm)/.(2.0*.dp);
  let h22 = (gp -. 2.0*.g +. gm) /. (dp*.dp) in

  printf "\n";

  printf "   dG/dp:  %12.4e %12.4e   (fwd FD)\n"  grdG_fwd.(0)  grdG_fwd.(1);
  printf "           %12.4e %12.4e   (bck FD)\n"  grdG_bck.(0)  grdG_bck.(1);
  printf "           %12.4e %12.4e   (cntr FD)\n" grdG_cntr.(0) grdG_cntr.(1);
  printf "\n";
  printf "  H(1,1):  %12.4e\n"  h11;
  printf "  H(2,2):  %12.4e\n"  h22


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
