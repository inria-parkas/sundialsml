(* Set vector of initial values
   The length of this vector determines the problem size. *)
let y = Nvector_array.wrap [| 0.0; 0.0; 0.0 |]

(* Create and initialize a solver session
   This will initialize a specific linear solver and the root-finding
   mechanism, if necessary. *)
let s = Ida.init (Ida.Spgmr spils_no_precond) tols f ~roots:(2, g) (3, y)

(* Specify integration tolerances (optional)}, e.g. *)
set_tolerances s SStolerances (reltol, abstol)

(* Set optional inputs}, e.g.
  Call any of the [set_*] functions to change solver parameters from their
  defaults. *)
set_stop_time s 10.0; ...

(* Advance solution in time}, e.g.,
   Repeatedly call either [solve_normal] or [solve_one_step] to
   advance the simulation. *)
let (t', result) = Ida.solve_normal s !t y in
...
t := t' + 0.1

(* Get optional outputs
   Call any of the [get_*] functions to examine solver statistics. *)
let stats = get_integrator_stats s in ...
