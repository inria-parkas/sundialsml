(* TODO: write a proper and working example... *)
(* TODO: include the quadrature module in this skeleton. *)

(* Initialize a session [s] per the skeleton at {!Cvode.init}
   Adding quadrature variables using {!Quadrature.init} if desired. *)
module Adj = Cvodes.Adjoint
...

(* Initialize the adjoint computation *)
Adj.init s nsteps Adj.IHermite

(* Integrate forward problem *)
let t, ncheck, r = Adj.forward_normal s tout y0

(* Setup the backward problem and attach a linear solver *)
let yB0 = RealArray.of_list [0.0; 0.0; ...]
let bs = Adj.init_backward s lmm (Newton ...) (SStolerances ...) fB tB0 yB0

(* Set optional inputs *)
Adj.set_max_ord bs ...

(* Initialize quadrature calculation *)
Quadrature.init bs fQb yQB0

(* Integrate backward problem *)
Adj.backward_normal s tB

(* Extract quadrature variables *)
let t = Quadrature.get s yQS

