open Sundials

(* TODO: write a proper and working example... *)
(* TODO: include the quadrature module in this skeleton. *)

(* Initialize a session [s] per the skeleton at {!Ida.init} *)
module Adj = Idas.Adjoint
...

(* Adding quadrature variables using {!Quadrature.init} if desired.
   Initial value correction, if needed, should be held off
   until sensitivity calculations are activated below. *)

(* Initialize the adjoint computation *)
init s nsteps IHermite

(* Integrate forward problem *)
let t, ncheck, r = forward_normal s tout y0

(* Setup the backward problem and attach a linear solver *)
let yB0  = RealArray.of_list [0.0; 0.0; ...]
let yB'0 = RealArray.of_list [0.0; 0.0; ...]
let bs = Adj.init_backward s (Spils.spgmr ...) (SStolerances ...) (NoSens fB) tB0 yB0 yB'0

(* Set optional inputs *)
Adj.set_max_ord bs ...

(* Initialize quadrature calculation *)
Adj.Quadrature.init bs fQb yQB0

(* Integrate backward problem *)
Adj.backward_normal s tB

(* Extract quadrature variables *)
let t = Adj.Quadrature.get s yQS

