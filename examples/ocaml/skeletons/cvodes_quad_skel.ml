(* TODO: write a proper and working example... *)

(* Initialize a session [s] per the skeleton at {!Cvode.session}}
   The vector of initial values should not include values for the
   quadrature variables. *)
...
module Quad = Cvode.Quadrature

(* Set vector of quadrature variables
   The length of this vector determines the number of quadrature variables. *)
let yQ = Cvode.RealArray.of_array [| 0.0; 0.0 |]

(* Initialize quadrature integration *)
Quad.init s fQ yQ

(* Specify integration tolerances (optional)}, e.g. *)
Quad.set_tolerances s SStolerances (reltol, abstol)

(* Advance the solution in time as per normal *)

(* Extract quadrature variables *)
Quad.get s yQ

(* Get quadrature optional outputs
   Call any of the [get_*] functions to examine solver statistics. *)
let nre = Quad.get_num_rhs_evals s in ...

