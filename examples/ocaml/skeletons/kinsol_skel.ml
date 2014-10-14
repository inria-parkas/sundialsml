(* Set vector with initial guess
   The length of this vector determines the problem size. *)
let u = Nvector_array.wrap [| 0.0; 0.0; 0.0 |]
(* Create and initialize a solver session
   This will also initialize the specified linear solver. *)
let s = Kinsol.init (Kinsol.Spgmr callbacks) f u
(* Set optional inputs,
   Call any of the [set_*] functions to change solver and linear solver
   parameters from their defaults. *)
set_num_max_iters s 500;

(* Solve problem *)
let result = Kinsol.solve s u u_scale f_scale
(* Get optional outputs
   Call any of the [get_*] functions to examine solver and linear solver
   statistics. *)
let fnorm = Kinsol.get_func_norm s in ...
