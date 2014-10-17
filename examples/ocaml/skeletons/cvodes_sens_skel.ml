(* TODO: write a proper and working example... *)
(* TODO: include the quadrature module in this skeleton. *)

(* Initialize a session [s] per the skeleton at
   {!Cvode.session} or {!Quadrature.init} *)
module Sens = Cvodes.Sensitivity
...

(* Define the sensitivity problem *)
let p = Cvode.RealArray.make np in
let sp = { pvals = Some p; pbar = ...; plist = ... }

(* Set sensitivity initial conditions *)
let yS0 = Array.init ns (fun _ -> RealArray.init neq 0.0)

(* Activate sensitivity calculations *)
Sens.init s (Sens.SStolerances ...) Sens.Simultaneous sp fS yS0;

(* Set optional inputs *)
Sens.set_dq_method s ...

(* Advance the solution in time as per normal *)
(* Extract sensitivity solution *)
let t = Sens.get s yS

