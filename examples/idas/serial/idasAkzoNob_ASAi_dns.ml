(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/04/17 20:12:55 $
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
 * Adjoint sensitivity example problem
 *
 * This IVP is a stiff system of 6 non-linear DAEs of index 1. The
 * problem originates from Akzo Nobel Central research in Arnhern,
 * The Netherlands, and describes a chemical process in which 2
 * species are mixed, while carbon dioxide is continuously added.
 * See http://pitagora.dm.uniba.it/~testset/report/chemakzo.pdf
 *
 * IDAS also computes the sensitivities with respect to initial
 * conditions of the following quantity:
 *   G = int_t0^t1 y1 dt
 * The sensitivity of G is the solution of the adjoint system at t0.
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray
module Quad = Idas.Quadrature
module Sens = Idas.Sensitivity
module Adjoint = Idas.Adjoint

let printf = Printf.printf

let nvconst = Nvector_serial.DataOps.n_vconst
let nvscale = Nvector_serial.DataOps.n_vscale

(* Problem Constants *)
let neq = 6
let t0 = 0.0
let tf = 180.0

let rtol =  1.0e-08
let atol =  1.0e-10
let rtolb = 1.0e-06
let atolb = 1.0e-08
let rtolq = 1.0e-10
let atolq = 1.0e-12

let steps = 150

type user_data = { k1 : float;
                   k2 : float;
                   k3 : float;
                   k4 : float;
                   k : float;
                   klA : float;
                   ks : float;
                   pCO2 : float;
                   h : float; }

let r_power_i base exponent =
  let go prod expt =
    let r = ref 1.0 in
    for i = 0 to expt - 1 do
      r := !r *. base
    done;
    !r
  in
  if exponent < 0 then 1. /.go  1.0 (- exponent)
  else go 1.0 exponent

let res data t (yy : RealArray.t) (yd : RealArray.t) (res : RealArray.t) =
  let k1 = data.k1
  and k2 = data.k2
  and k3 = data.k3
  and k4 = data.k4
  and k = data.k
  and klA = data.klA
  and ks = data.ks
  and pCO2 = data.pCO2
  and h = data.h

  and y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}
  and y4 = yy.{3}
  and y5 = yy.{4}
  and y6 = yy.{5}

  and yd1 = yd.{0}
  and yd2 = yd.{1}
  and yd3 = yd.{2}
  and yd4 = yd.{3}
  and yd5 = yd.{4}
  in

  let r1 = k1 *. r_power_i y1 4 *. sqrt y2
  and r2 = k2 *. y3 *. y4
  and r3 = k2/.k *. y1 *. y5
  and r4 = k3 *. y1 *. y4 *. y4
  and r5 = k4 *. y6 *. y6 *. sqrt y2
  and fin = klA *. ( pCO2/.h -. y2 )
  in

  res.{0} <- yd1 +. 2.0*.r1 -. r2 +. r3 +. r4;
  res.{1} <- yd2 +. 0.5*.r1 +. r4 +. 0.5*.r5 -. fin;
  res.{2} <- yd3 -. r1 +. r2 -. r3;
  res.{3} <- yd4 +. r2 -. r3 +. 2.0*.r4;
  res.{4} <- yd5 -. r2 +. r3 -. r5;
  res.{5} <- ks*.y1*.y4 -. y6

(*
 * rhsQ routine. Computes quadrature(t,y).
 *)
let rhsQ data t yy yp qdot =
  qdot.{0} <- yy.{0}

(*
 * resB routine. Residual for adjoint system.
 *)
let resB : user_data -> RealArray.t Adjoint.bresfn_no_sens =
  fun data { Adjoint.y = yy; Adjoint.yb = yyB; Adjoint.yb' = ypB } rrB ->
  let k1 = data.k1
  and k2 = data.k2
  and k3 = data.k3
  and k4 = data.k4
  and k = data.k
  and klA = data.klA
  and ks = data.ks

  and y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}
  and y4 = yy.{3}
  and y5 = yy.{4}
  and y6 = yy.{5}

  and yB1 = yyB.{0}
  and yB2 = yyB.{1}
  and yB3 = yyB.{2}
  and yB4 = yyB.{3}
  and yB5 = yyB.{4}
  and yB6 = yyB.{5}

  and ypB1 = ypB.{0}
  and ypB2 = ypB.{1}
  and ypB3 = ypB.{2}
  and ypB4 = ypB.{3}
  and ypB5 = ypB.{4}
  in

  let y2tohalf = sqrt y2
  and y1to3 = y1*.y1*.y1
  and k2overk = k2/.k
  in

  let tmp1 = k1*. y1to3 *. y2tohalf
  and tmp2 = k3*.y4*.y4
  in
  rrB.{0} <- 1.0 +.  ypB1 -. (8.0*.tmp1 +. k2overk*.y5 +. tmp2)*.yB1
    -. (2.0*.tmp1+.tmp2)*.yB2 +. (4.0*.tmp1+.k2overk*.y5)*.yB3
    +. k2overk*.y5*.(yB4-.yB5) -. 2.0*.tmp2*.yB4 +. ks*.y4*.yB6;

  let tmp1 = k1 *. y1*.y1to3 *. (y2tohalf/.y2)
  and tmp2 = k4 *. y6*.y6 *. (y2tohalf/.y2)
  in
  rrB.{1} <- ypB2 -. tmp1*.yB1 -. (0.25*.tmp1 +. 0.25*.tmp2 +. klA)*.yB2
    +. 0.5*.tmp1*.yB3 +. 0.5*.tmp2*.yB5;

  rrB.{2} <- ypB3 +. k2*.y4*.(yB1-.yB3-.yB4+.yB5);

  let tmp1 = k3*.y1*.y4
  and tmp2 = k2*.y3
  in
  rrB.{3} <- ypB4 +. (tmp2-.2.0*.tmp1)*.yB1 -. 2.0*.tmp1*.yB2 -. tmp2*.yB3
    -. (tmp2+.4.0*.tmp1)*.yB4 +. tmp2*.yB5 +. ks*.y1*.yB6;

  rrB.{4} <- ypB5 -. k2overk*.y1*.(yB1-.yB3-.yB4+.yB5);

  rrB.{5} <- k4*.y6*.y2tohalf*.(2.0*.yB5-.yB2) -. yB6

(*
 * Print results after backward integration
 *)
let print_output tfinal yB ypB =
  printf "dG/dy0: \t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n"
         yB.{0} yB.{1} yB.{2} yB.{3} yB.{4} yB.{5};
  printf "--------------------------------------------------------\n\n"

(* Main program *)
let main () =
  printf "\nAdjoint Sensitivity Example for Akzo-Nobel Chemical Kinetics\n";
  printf "-------------------------------------------------------------\n";
  printf "Sensitivity of G = int_t0^tf (y1) dt with respect to IC.\n";
  printf "-------------------------------------------------------------\n\n";
  (* Fill user's data with the appropriate values for coefficients. *)
  let data = { k1 = 18.7;
               k2 = 0.58;
               k3 = 0.09;
               k4 = 0.42;
               k = 34.4;
               klA = 3.3;
               ks = 115.83;
               pCO2 = 0.9;
               h = 737.0;
             }
  in

  (* Allocate N-vectors. *)
  let yy = RealArray.create neq
  and yp = RealArray.create neq
  in

  (* Consistent IC for  y, y'. *)
  let y01 = 0.444
  and y02 = 0.00123
  and y03 = 0.00
  and y04 = 0.007
  and y05 = 0.0
  in
  yy.{0} <- y01;
  yy.{1} <- y02;
  yy.{2} <- y03;
  yy.{3} <- y04;
  yy.{4} <- y05;
  yy.{5} <- data.ks *. y01 *. y04;

  (* Get y' = - res(t0, y, 0) *)
  nvconst 0. yp;

  let rr = RealArray.create neq in
  res data t0 yy yp rr;
  nvscale (-1.0) rr yp;

  (* Create and initialize q0 for quadratures. *)
  let q = RealArray.create 1 in
  q.{0} <- 0.0;

  (* Wrap arrays in nvectors.  *)
  let wyy = Nvector_serial.wrap yy
  and wyp = Nvector_serial.wrap yp
  and wq  = Nvector_serial.wrap q
  in

  (* Call IDACreate and IDAInit to initialize IDA memory *)
  let m = Matrix.dense neq in
  let mem = Ida.(init Dls.(solver Lsolver.Direct.(dense wyy m))
                      (SStolerances (rtol,atol))
                      (res data) t0 wyy wyp)
  in

  (* Initialize QUADRATURE(S). *)
  Quad.init mem (rhsQ data) wq;

  (* Set tolerances and error control for quadratures. *)
  Quad.(set_tolerances mem (SStolerances (rtolq,atolq)));

  (* Prepare ADJOINT. *)
  Adjoint.(init mem steps IHermite);

  (* FORWARD run. *)
  printf "Forward integration ... ";
  let _ = Adjoint.forward_normal mem tf wyy wyp in
  let nst = Ida.get_num_steps mem in

  printf "done ( nst = %d )\n" nst;

  let _ = Quad.get mem wq in

  printf "G:          %24.16f \n" q.{0};
  printf "--------------------------------------------------------\n\n";


  (* BACKWARD run *)

  (* Initialize yB *)
  let yB = RealArray.create neq in
  nvconst 0. yB;

  let ypB = RealArray.create neq in
  nvconst 0. ypB;
  ypB.{0} <- - 1.0;

  let wyB = Nvector_serial.wrap yB
  and wypB = Nvector_serial.wrap ypB
  in

  let m = Matrix.dense neq in
  let indexB = Adjoint.(init_backward mem Dls.(solver
                                               Lsolver.Direct.(dense wyB m))
                                      (SStolerances (rtolb, atolb))
                                      (NoSens (resB data))
                                      tf wyB wypB)
  in
  Adjoint.set_max_num_steps indexB 1000;

  printf "Backward integration ... ";

  Adjoint.backward_normal mem t0;

  let nstB = Adjoint.get_num_steps indexB in

  printf "done ( nst = %d )\n" nstB;

  let time = Adjoint.get indexB wyB wypB in

  print_output time yB ypB

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
