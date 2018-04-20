(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2009/09/30 23:33:29 $
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
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray
module Quad = Idas.Quadrature
module Sens = Idas.Sensitivity
module QuadSens = Idas.Sensitivity.Quadrature

let printf = Printf.printf

let nvconst = Nvector_serial.DataOps.n_vconst
let nvscale = Nvector_serial.DataOps.n_vscale

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

(* Problem Constants *)
let neq = 6

let t0 = 0.0

let t1 = 1e-8  (* first time for output *)

let tf = 180.0 (* Final time. *)
let nf = 25    (* Total number of outputs. *)

let rtol =  1.0e-08
let atol =  1.0e-10
let rtolq = 1.0e-10
let atolq = 1.0e-12


type user_data = { k1 : float;
                   k2 : float;
                   k3 : float;
                   k4 : float;
                   k : float;
                   klA : float;
                   ks : float;
                   pCO2 : float;
                   h : float;
                 }

let res data t (y : RealArray.t) (yd : RealArray.t) (res : RealArray.t) =
  let k1 = data.k1
  and k2 = data.k2
  and k3 = data.k3
  and k4 = data.k4
  and k = data.k
  and klA = data.klA
  and ks = data.ks
  and pCO2 = data.pCO2
  and h = data.h
  in

  let r1 = k1 *. (r_power_i y.{0} 4) *. sqrt y.{1}
  and r2 = k2 *. y.{2} *. y.{3}
  and r3 = k2/.k *. y.{0} *. y.{4}
  and r4 = k3 *. y.{0} *. y.{3} *. y.{3}
  and r5 = k4 *. y.{5} *. y.{5} *. sqrt y.{1}
  and fin = klA *. ( pCO2/.h -. y.{1} )
  in

  res.{0} <- yd.{0} +. 2.0*.r1 -. r2 +. r3 +. r4;
  res.{1} <- yd.{1} +. 0.5*.r1 +. r4 +. 0.5*.r5 -. fin;
  res.{2} <- yd.{2} -. r1 +. r2 -. r3;
  res.{3} <- yd.{3} +. r2 -. r3 +. 2.0*.r4;
  res.{4} <- yd.{4} -. r2 +. r3 -. r5;
  res.{5} <- ks*.y.{0}*.y.{3} -. y.{5}

(*
 * rhsQ routine. Computes quadrature(t,y).
 *)
let rhsQ data t (yy : RealArray.t) yp (qdot : RealArray.t) =
  qdot.{0} <- yy.{0}

let print_header rtol avtol y =
  print_string "\nidasAkzoNob_dns: Akzo Nobel chemical kinetics DAE serial example problem for IDAS\n";
  print_string "Linear solver: IDADENSE, Jacobian is computed by IDAS.\n";
  printf "Tolerance parameters:  rtol = %g   atol = %g\n"
         rtol avtol;
  print_string "---------------------------------------------------------------------------------\n";
  print_string "   t        y1        y2       y3       y4       y5";
  print_string "      y6    | nst  k      h\n";
  print_string "---------------------------------------------------------------------------------\n"


let print_output mem t y =
  let kused = Ida.get_last_order mem
  and nst   = Ida.get_num_steps mem
  and hused = Ida.get_last_step mem
  in
  printf "%8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e | %3d  %1d %8.2e\n"
         t y.{0} y.{1} y.{2} y.{3} y.{4} y.{5} nst kused hused


let print_final_stats mem =
  let open Ida in
  let nst   = get_num_steps mem
  and nre   = get_num_res_evals mem
  and nje   = Dls.get_num_jac_evals mem
  and nni   = get_num_nonlin_solv_iters mem
  and netf  = get_num_err_test_fails mem
  and ncfn  = get_num_nonlin_solv_conv_fails mem
  and nreLS = Dls.get_num_res_evals mem
  in

  print_string "\nFinal Run Statistics: \n\n";
  print_string "Number of steps                    = ";   print_int nst;
  print_string "\nNumber of residual evaluations     = "; print_int (nre+nreLS);
  print_string "\nNumber of Jacobian evaluations     = "; print_int nje;
  print_string "\nNumber of nonlinear iterations     = "; print_int nni;
  print_string "\nNumber of error test failures      = "; print_int netf;
  print_string "\nNumber of nonlinear conv. failures = "; print_int ncfn;
  print_newline ()

(* Main program *)
let main () =

  (* Fill user's data with the appropriate values for coefficients. *)
  let data = { k1 = 18.7;
               k2 = 0.58;
               k3 = 0.09;
               k4 = 0.42;
               k = 34.4;
               klA = 3.3;
               ks = 115.83;
               pCO2 = 0.9;
               h = 737.0; }
  in

  (* Allocate N-vectors. *)
  let yy = RealArray.create neq
  and yp = RealArray.create neq
  in

  (* Consistent IC for  y, y'. *)
  let y01 =0.444
  and y02 =0.00123
  and y03 =0.00
  and y04 =0.007
  and y05 =0.0
  in
  yy.{0} <- y01;
  yy.{1} <- y02;
  yy.{2} <- y03;
  yy.{3} <- y04;
  yy.{4} <- y05;
  yy.{5} <- data.ks *. y01 *. y04;

  (* Get y' = - res(t0, y, 0) *)
  nvconst 0.0 yp;

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
  let mem = Ida.(init Dls.(solver Direct.(dense wyy m) m)
                      (SStolerances (rtol,atol))
                      (res data) t0 wyy wyp)
  in

  (* Initialize QUADRATURE(S). *)
  Quad.init mem (rhsQ data) wq;

  (* Set tolerances and error control for quadratures. *)
  Quad.(set_tolerances mem (SStolerances (rtolq,atolq)));

  print_header rtol atol yy;
  (* Print initial states *)
  print_output mem 0.0 yy;

  let tout = ref t1
  and incr = (tf/.t1) ** (1.0 /. float_of_int nf)
  in

  (* FORWARD run. *)
  for nout = 0 to nf do

    let (time, _) = Ida.solve_normal mem !tout wyy wyp in
    print_output mem time yy;

    tout := !tout *. incr;
  done;

  let _ = Quad.get mem wq in

  print_string "\n--------------------------------------------------------\n";
  printf "G:          %24.16f \n" q.{0};
  print_string "--------------------------------------------------------\n\n";

  print_final_stats mem

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
