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

let nvconst = Nvector_array.Bigarray.array_nvec_ops.Nvector_custom.nvconst
let nvscale = Nvector_array.Bigarray.array_nvec_ops.Nvector_custom.nvscale

let r_power_i base exponent =
  let rec go prod expt =
    if expt = 0 then prod
    else go (prod *. base) (expt - 1)
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

let res data t yy yd res =
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

  let r1 = k1 *. (r_power_i y1 4) *. sqrt y2
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

let print_header rtol avtol y =
  printf "\nidasAkzoNob_dns: Akzo Nobel chemical kinetics DAE serial example problem for IDAS\n";
  printf "Linear solver: IDADENSE, Jacobian is computed by IDAS.\n";
  printf "Tolerance parameters:  rtol = %g   atol = %g\n"
         rtol avtol;
  printf "---------------------------------------------------------------------------------\n";
  printf "   t        y1        y2       y3       y4       y5";
  printf "      y6    | nst  k      h\n";
  printf "---------------------------------------------------------------------------------\n"


let print_output mem t y =
  let kused = Ida.get_last_order mem
  and nst = Ida.get_num_steps mem
  and hused = Ida.get_last_step mem
  in
  printf "%8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e | %3d  %1d %8.2e\n"
         t y.{0} y.{1} y.{2} y.{3} y.{4} y.{5} nst kused hused


let print_final_stats mem =
  let nst = Ida.get_num_steps mem
  and nre = Ida.get_num_res_evals mem
  and nje = Ida.Dls.get_num_jac_evals mem
  and nni = Ida.get_num_nonlin_solv_iters mem
  and netf = Ida.get_num_err_test_fails mem
  and ncfn = Ida.get_num_nonlin_solv_conv_fails mem
  and nreLS = Ida.Dls.get_num_res_evals mem
  in

  printf "\nFinal Run Statistics: \n\n";
  printf "Number of steps                    = %d\n" nst;
  printf "Number of residual evaluations     = %d\n" (nre+nreLS);
  printf "Number of Jacobian evaluations     = %d\n" nje;
  printf "Number of nonlinear iterations     = %d\n" nni;
  printf "Number of error test failures      = %d\n" netf;
  printf "Number of nonlinear conv. failures = %d\n" ncfn

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
  let mem = Ida.init (Ida.Dls.dense None) (Ida.SStolerances (rtol,atol))
                     (res data) ~t0:t0 wyy wyp
  in

  (* Initialize QUADRATURE(S). *)
  Quad.init mem (rhsQ data) wq;

  (* Set tolerances and error control for quadratures. *)
  Quad.ss_tolerances mem rtolq atolq;
  Quad.set_err_con mem true;

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

  printf "\n--------------------------------------------------------\n";
  printf "G:          %24.16f \n" q.{0};
  printf "--------------------------------------------------------\n\n";

  print_final_stats mem

(* Check if the last argument is a repetition count.  *)
let reps =
  match Sys.argv with
  | [|_; n|] -> int_of_string n
  | _ -> 1
let _ = for i = 1 to reps do main () done
let _ = Gc.full_major ()

