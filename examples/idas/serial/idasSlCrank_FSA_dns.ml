(*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2009/04/29 20:40:07 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban and Cosmin Petra @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Jul 2014.
 * -----------------------------------------------------------------
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The mechanism moves under the action of a constant horizontal 
 * force applied to the connecting rod and a spring-damper connecting 
 * the crank and connecting rod.
 *
 * The equations of motion are formulated as a system of stabilized
 * index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
 *
 * IDAS also computes sensitivities with respect to the problem
 * parameters k (spring constant) and c (damper constant) of the
 * kinetic energy:
 *   G = int_t0^tend g(t,y,p) dt,
 * where
 *   g(t,y,p) = 0.5*J1*v1^2 + 0.5*J2*v3^2 + 0.5*m2*v2^2
 *
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray
module Quad = Idas.Quadrature
module Sens = Idas.Sensitivity
module QuadSens = Sens.Quadrature
module Adjoint = Idas.Adjoint
module AdjQuad = Adjoint.Quadrature
module VarId = Ida.VarId

let printf = Printf.printf

let nvconst = Nvector_serial.DataOps.n_vconst

(* Problem Constants *)

let neq =   10
let np =     2

let tbegin =  0.0
let tend =    10.000

let rtolf =   1.0e-06
let atolf =   1.0e-07

let rtolq =   1.0e-06
let atolq =   1.0e-08

let rtolfd =  1.0e-06
let atolfd =  1.0e-08

type user_data = { a : float;
                   j1 : float;
                   j2 : float;
                   m1 : float;
                   m2 : float;
                   l0 : float;
                   params : RealArray.t; (* size 2 *)
                   f : float;
                 }

let force data (yy : RealArray.t) (_Q : RealArray.t) =

  let a = data.a in
  let k = data.params.{0} in
  let c = data.params.{1} in
  let l0 = data.l0 in
  let _F = data.f in

  let q = yy.{0} in
  let x = yy.{1} in
  let p = yy.{2} in

  let qd = yy.{3} in
  let xd = yy.{4} in
  let pd = yy.{5} in

  let s1 = sin(q) in
  let c1 = cos(q) in
  let s2 = sin(p) in
  let c2 = cos(p) in
  let s21 = s2*.c1 -. c2*.s1 in
  let c21 = c2*.c1 +. s2*.s1 in

  let l2 = x*.x -. x*.(c2+.a*.c1) +. (1.0 +. a*.a)/.4.0 +. a*.c21/.2.0 in
  let l = sqrt l2 in
  let ld = (2.0*.x*.xd -. xd*.(c2+.a*.c1) +. x*.(s2*.pd+.a*.s1*.qd) -. a*.s21*.(pd-.qd)/.2.0) in
  let ld = ld /. (2.0*.l)
  in

  let f = k*.(l-.l0) +. c*.ld in
  let fl = f/.l in

  _Q.{0} <- -. fl *. a *. (s21/.2.0 +. x*.s1) /. 2.0;
  _Q.{1} <- fl *. (c2/.2.0 -. x +. a*.c1/.2.0) +. _F;
  _Q.{2} <- -. fl *. (x*.s2 -. a*.s21/.2.0) /. 2.0 -. _F*.s2

let set_ic data yy yp =
  nvconst 0.0 yy;
  nvconst 0.0 yp;

  let pi = 4.0*.atan(1.0) in

  let a = data.a in
  let j1 = data.j1 in
  let m2 = data.m2 in
  let j2 = data.j2 in

  let q = pi/.2.0 in
  let p = asin(-.a) in
  let x = cos(p) in

  yy.{0} <- q;
  yy.{1} <- x;
  yy.{2} <- p;

  let _Q = RealArray.create 3 in

  force data yy _Q;

  yp.{3} <- _Q.{0}/.j1;
  yp.{4} <- _Q.{1}/.m2;
  yp.{5} <- _Q.{2}/.j2

let ressc data tres (yval : RealArray.t) (ypval : RealArray.t)
                    (rval : RealArray.t) =
  let a  = data.a in
  let j1 = data.j1 in
  let m2 = data.m2 in
  let j2 = data.j2 in

  let q = yval.{0} in
  let x = yval.{1} in
  let p = yval.{2} in

  let qd = yval.{3} in
  let xd = yval.{4} in
  let pd = yval.{5} in

  let lam1 = yval.{6} in
  let lam2 = yval.{7} in

  let mu1 = yval.{8} in
  let mu2 = yval.{9} in

  let s1 = sin(q) in
  let c1 = cos(q) in
  let s2 = sin(p) in
  let c2 = cos(p) in

  let _Q = RealArray.create 3 in

  force data yval _Q;

  rval.{0} <- ypval.{0} -. qd +. a*.s1*.mu1 -. a*.c1*.mu2;
  rval.{1} <- ypval.{1} -. xd +. mu1;
  rval.{2} <- ypval.{2} -. pd +. s2*.mu1 -. c2*.mu2; 

  rval.{3} <- j1*.ypval.{3} -. _Q.{0} +. a*.s1*.lam1 -. a*.c1*.lam2;
  rval.{4} <- m2*.ypval.{4} -. _Q.{1} +. lam1;
  rval.{5} <- j2*.ypval.{5} -. _Q.{2} +. s2*.lam1 -. c2*.lam2; 

  rval.{6} <- x -. c2 -. a*.c1;
  rval.{7} <- -.s2 -. a*.s1;

  rval.{8} <- a*.s1*.qd +. xd +. s2*.pd;
  rval.{9} <- -.a*.c1*.qd -. c2*.pd

let rhsQ data t (yy : RealArray.t) (yp : RealArray.t) (qdot : RealArray.t) =
  let j1 = data.j1 in
  let m2 = data.m2 in
  let j2 = data.j2 in

  let v1 = yy.{3} in
  let v2 = yy.{4} in
  let v3 = yy.{5} in

  qdot.{0} <- 0.5*.(j1*.v1*.v1 +. m2*.v2*.v2 +. j2*.v3*.v3)

let rhsQS : user_data -> RealArray.t QuadSens.quadsensrhsfn =
  fun data args rhsQS ->
  let yy = args.QuadSens.y
  and yyS = args.QuadSens.yS
  in
  let j1 = data.j1 in
  let m2 = data.m2 in
  let j2 = data.j2 in

  let v1 = yy.{3} in
  let v2 = yy.{4} in
  let v3 = yy.{5} in

  (* Sensitivities of v. *)
  let s1 = yyS.(0).{3} in
  let s2 = yyS.(0).{4} in
  let s3 = yyS.(0).{5} in

  rhsQS.(0).{0} <- j1*.v1*.s1 +. m2*.v2*.s2 +. j2*.v3*.s3;

  let s1 = yyS.(1).{3} in
  let s2 = yyS.(1).{4} in
  let s3 = yyS.(1).{5} in

  rhsQS.(1).{0} <- j1*.v1*.s1 +. m2*.v2*.s2 +. j2*.v3*.s3

let print_final_stats mem =
  let nst = Ida.get_num_steps mem in
  let nre = Ida.get_num_res_evals mem in
  let nje = Ida.Dls.get_num_jac_evals mem in
  let nni = Ida.get_num_nonlin_solv_iters mem in
  let netf = Ida.get_num_err_test_fails mem in
  let ncfn = Ida.get_num_nonlin_solv_conv_fails mem in
  let nreLS = Ida.Dls.get_num_res_evals mem in

  printf "\nFinal Run Statistics: \n\n";
  printf "Number of steps                    = %d\n" nst;
  printf "Number of residual evaluations     = %d\n" (nre+nreLS);
  printf "Number of Jacobian evaluations     = %d\n" nje;
  printf "Number of nonlinear iterations     = %d\n" nni;
  printf "Number of error test failures      = %d\n" netf;
  printf "Number of nonlinear conv. failures = %d\n" ncfn

(*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 *)

let main () =
  let pbar = RealArray.create 2
  and gm = RealArray.create 2
  and gp = RealArray.create 2
  and atolS = RealArray.create np
  in

  let id = RealArray.create neq in
  let yy = RealArray.create neq in
  let yp = RealArray.create neq in
  let q = RealArray.create 1 in

  let yyS= Array.init np (fun _ -> Nvector_serial.wrap (RealArray.copy yy)) in
  let ypS= Array.init np (fun _ -> Nvector_serial.wrap (RealArray.copy yp)) in
  let qS = Array.init np (fun _ -> Nvector_serial.wrap (RealArray.copy q)) in

  let data = { a = 0.5;   (* half-length of crank *)
               j1 = 1.0;  (* crank moment of inertia *)
               m2 = 1.0;  (* mass of connecting rod *)
               m1 = 1.0;
               j2 = 2.0;  (* moment of inertia of connecting rod *)
               params = RealArray.of_array [|1.0;   (* spring constant *)
                                             1.0|]; (* damper constant *)
               l0 = 1.0;  (* spring free length *)
               f = 1.0;   (* external constant force *)
             }
  in

  nvconst VarId.differential id;
  id.{9} <- VarId.algebraic;
  id.{8} <- VarId.algebraic;
  id.{7} <- VarId.algebraic;
  id.{6} <- VarId.algebraic;

  printf "\nSlider-Crank example for IDAS:\n";

  (* Consistent IC*)
  set_ic data yy yp;

  for is = 0 to np - 1 do
    nvconst 0.0 (Nvector.unwrap yyS.(is));
    nvconst 0.0 (Nvector.unwrap ypS.(is));
  done;

  (* Wrap arrays in nvectors.  Operations performed on the wrapped
     representation affect the original bigarrays.  *)
  let wyy = Nvector_serial.wrap yy
  and wyp = Nvector_serial.wrap yp
  and wid = Nvector_serial.wrap id
  and wq  = Nvector_serial.wrap q
  in

  (* IDA initialization *)
  (* Call IDADense and set up the linear solver. *)
  let mem =
    Ida.init (Ida.Dls.dense ())
      (Ida.SStolerances (rtolf, atolf))
      (ressc data)
      ~varid:wid
      tbegin
      wyy wyp
  in
  Ida.set_suppress_alg mem true;
  Ida.set_max_num_steps mem 20000;

  pbar.{0} <- data.params.{0}; pbar.{1} <- data.params.{1};
  Sens.init mem Sens.EEtolerances Sens.Simultaneous
    { Sens.pvals = Some data.params;
      Sens.pbar = Some pbar;
      Sens.plist = None }
    yyS
    ypS;
  Sens.set_err_con mem true;

  nvconst 0.0 q;
  Quad.init mem (rhsQ data) wq;
  Quad.set_tolerances mem (Quad.SStolerances (rtolq, atolq)) ;

  nvconst 0.0 (Nvector.unwrap qS.(0));
  QuadSens.init mem ~fQS:(rhsQS data) qS;
  atolS.{0} <- atolq; atolS.{1} <- atolq;
  QuadSens.set_tolerances mem (QuadSens.SStolerances (rtolq, atolS));

  (* Perform forward run *)
  printf "\nForward integration ... ";

  let _ = Ida.solve_normal mem tend wyy wyp in

  printf "done!\n";

  print_final_stats mem;

  let _ = Quad.get mem wq in
  printf "--------------------------------------------\n";
  printf "  G = %24.16f\n" q.{0};
  printf "--------------------------------------------\n\n";

  let _ = QuadSens.get mem qS in
  printf "-------------F O R W A R D------------------\n";
  printf "   dG/dp:  %12.4e %12.4e\n" (Nvector.unwrap qS.(0)).{0}
                                      (Nvector.unwrap qS.(1)).{0};
  printf "--------------------------------------------\n\n";



  (* Finite differences for dG/dp *)
  let dp = 0.00001 in
  data.params.{0} <- 1.0;
  data.params.{1} <- 1.0;

  set_ic data yy yp;

  (* Call IDADense and set up the linear solver. *)
  let mem =
    Ida.init (Ida.Dls.dense ())
      (Ida.SStolerances (rtolfd, atolfd))
      (ressc data)
      ~varid:wid
      tbegin wyy wyp
  in
  Ida.set_suppress_alg mem true;

  nvconst 0.0 q;
  Quad.init mem (rhsQ data) wq;
  Quad.set_tolerances mem (Quad.SStolerances (rtolq, atolq));

  let _ = Ida.solve_normal mem tend wyy wyp in

  let _ = Quad.get mem wq in
  let g = q.{0} in
  (*printf "  G  =%12.6e\n" q.{0};*)

  (******************************
  * BACKWARD for k
  ******************************)
  data.params.{0} <- data.params.{0} -. dp;
  set_ic data yy yp;

  Ida.reinit mem tbegin wyy wyp;

  nvconst 0.0 q;
  Quad.reinit mem wq;

  let _ = Ida.solve_normal mem tend wyy wyp in
  let _ = Quad.get mem wq in
  gm.{0} <- q.{0};
  (*printf "Gm.{0}=%12.6e\n" q.{0};*)

  (****************************
  * FORWARD for k *
  ****************************)
  data.params.{0} <- data.params.{0} +. (2.0*.dp);
  set_ic data yy yp;
  Ida.reinit mem tbegin wyy wyp;

  nvconst 0.0 q;
  Quad.reinit mem wq;

  let _ = Ida.solve_normal mem tend wyy wyp in
  let _ = Quad.get mem wq in
  gp.{0} <- q.{0};
  (*printf "Gp.{0}=%12.6e\n" q.{0};*)


  (* Backward for c *)
  data.params.{0} <- 1.0;
  data.params.{1} <- data.params.{1} -. dp;
  set_ic data yy yp;
  Ida.reinit mem tbegin wyy wyp;

  nvconst 0.0 q;
  Quad.reinit mem wq;

  let _ = Ida.solve_normal mem tend wyy wyp in
  let _ = Quad.get mem wq in
  gm.{1} <- q.{0};

  (* Forward for c *)
  data.params.{1} <- data.params.{1} +. (2.0*.dp);
  set_ic data yy yp;
  Ida.reinit mem tbegin wyy wyp;

  nvconst 0.0 q;
  Quad.reinit mem wq;

  let _ = Ida.solve_normal mem tend wyy wyp in
  let _ = Quad.get mem wq in
  gp.{1} <- q.{0};

  printf "\n\n   Checking using Finite Differences \n\n";

  printf "---------------BACKWARD------------------\n";
  printf "   dG/dp:  %12.4e %12.4e\n" ((g-.gm.{0})/.dp) ((g-.gm.{1})/.dp);
  printf "-----------------------------------------\n\n";

  printf "---------------FORWARD-------------------\n";
  printf "   dG/dp:  %12.4e %12.4e\n" ((gp.{0}-.g)/.dp) ((gp.{1}-.g)/.dp);
  printf "-----------------------------------------\n\n";

  printf "--------------CENTERED-------------------\n";
  printf "   dG/dp:  %12.4e %12.4e\n" ((gp.{0}-.gm.{0})/.(2.0*.dp)) ((gp.{1}-.gm.{1})/.(2.0*.dp));
  printf "-----------------------------------------\n\n"


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
