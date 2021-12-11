open Sundials

(* TODO: write a proper and working example... *)
(* TODO: include the quadrature module in this skeleton. *)

(* Initialize a session [s] per the skeleton at
   {!Ida.session} or {!Idas.Quadrature.init} *)
module Sens = Idas.Sensitivity
...

(* Initial value correction, if needed, should be held off
  until sensitivity calculations are activated below. *)

(* Define the sensitivity problem *)
let p = Ida.RealArray.make np in
let sp = { pvals = Some p; pbar = ...; plist = ... }

(* Set sensitivity initial conditions *)
let yS0 = Array.init ns (fun _ -> RealArray.init neq 0.0)
let yS'0 = Array.init ns (fun _ -> RealArray.init neq 0.0)

(* Activate sensitivity calculations} *)
Sens.init s (SStolerances ...) Simultaneous sp fS yS0 y'S0

(* Correct initial values (optional) *)
Sens.calc_ic_y s tout1

(* Set optional inputs *)
Sens.set_dq_method s ...

(* Advance the solution in time as per normal *)

(* Extract sensitivity solution *)
let t = Sens.get s yS

