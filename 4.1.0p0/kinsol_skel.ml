open Sundials

(* 1. Define a system function. *)
let pi, e = 4. *. atan 1., exp 1.
let sysf u r =
  r.{0} <- 0.5 *. sin(u.{0}*.u.{1}) -. 0.25 *. u.{1} /. pi -. 0.5 *. u.{0};
  r.{1} <- (1. -. 0.25/.pi)*.(exp(2.*.u.{0})-.e) +. e*.u.{1}/.pi -. 2.*.e*.u.{0};
  r.{2} <- u.{2} -. u.{0} +. 0.25;
  r.{3} <- u.{3} -. u.{0} +. 1.0;
  r.{4} <- u.{4} -. u.{1} +. 1.5;
  r.{5} <- u.{5} -. u.{1} +. 2.0 *. pi

(* 2. Set vector with initial guess.
      The length of this vector determines the problem size. *)
let ud = RealArray.of_list [ 0.25; 1.5; 0.; -0.75; 0.; 1.5 -. 2. *.pi ]
let u = Nvector_serial.wrap ud

(* 3. Create and initialize a solver session.
      This will initialize a specific linear solver. *)
let m = Matrix.dense 6
let s = Kinsol.(init ~lsolver:(Dls.(solver (dense u m))) sysf u);;

(* 4. Set optional inputs, e.g.,
      call [set_*] functions to change solver parameters. *)
let c = RealArray.of_list [
   Constraint.unconstrained;
   Constraint.unconstrained;
   Constraint.geq_zero;
   Constraint.leq_zero;
   Constraint.geq_zero;
   Constraint.leq_zero ] in
Kinsol.set_constraints s (Nvector_serial.wrap c);
Kinsol.set_func_norm_tol s 1.0e-5;
Kinsol.set_scaled_step_tol s 1.0e-5;
Kinsol.set_max_setup_calls s 1;;

(* 5. Solve the problem. *)
let snv = Nvector_serial.make (RealArray.length ud) 1.0 in
ignore (Kinsol.(solve s u Newton snv snv));

Printf.printf "%8.5g %8.6g\n" ud.{0} ud.{1};;

(* 6. Get optional outputs,
      call the [get_*] functions to examine solver statistics. *)
let fnorm = Kinsol.get_func_norm s;;
