(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/17 19:38:48 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2014.
 * -----------------------------------------------------------------
 * Example (serial):
 *
 * This example solves a nonlinear system from.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 4 from Section 14.1, Chapter 14: Ferraris and Tronconi
 * 
 * This problem involves a blend of trigonometric and exponential terms.
 *    0.5 sin(x1 x2) - 0.25 x2/pi - 0.5 x1 = 0
 *    (1-0.25/pi) ( exp(2 x1)-e ) + e x2 / pi - 2 e x1 = 0
 * such that
 *    0.25 <= x1 <=1.0
 *    1.5 <= x2 <= 2 pi
 * 
 * The treatment of the bound constraints on x1 and x2 is done using
 * the additional variables
 *    l1 = x1 - x1_min >= 0
 *    L1 = x1 - x1_max <= 0
 *    l2 = x2 - x2_min >= 0
 *    L2 = x2 - x2_max >= 0
 * 
 * and using the constraint feature in KINSOL to impose
 *    l1 >= 0    l2 >= 0
 *    L1 <= 0    L2 <= 0
 * 
 * The Ferraris-Tronconi test problem has two known solutions.
 * The nonlinear system is solved by KINSOL using different 
 * combinations of globalization and Jacobian update strategies 
 * and with different initial guesses (leading to one or the other
 * of the known solutions).
 *
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray

let printf = Printf.printf
let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

(* Problem Constants *)

let nvar   = 2
let neq    = 3 * nvar

let ftol   = 1.0e-5 (* function tolerance *)
let stol   = 1.0e-5 (* step tolerance *)

let zero   = 0.0
let pt25   = 0.25
let pt5    = 0.5
let one    = 1.0
let onept5 = 1.5
let two    = 2.0

let pi     = 3.1415926
let e      = 2.7182818

(* Accessor macros *)
let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

(* System function for predator-prey system *)
let lb = RealArray.make nvar
let ub = RealArray.make nvar

let func udata fdata =
  let x1  = udata.{0} in
  let x2  = udata.{1} in
  let l1  = udata.{2} in
  let l1' = udata.{3} in
  let l2  = udata.{4} in
  let l2' = udata.{5} in

  fdata.{0} <- pt5 *. sin(x1*.x2) -. pt25 *. x2 /. pi -. pt5 *. x1;
  fdata.{1} <- (one -. pt25/.pi)*.(exp(two*.x1)-.e) +. e*.x2/.pi -. two*.e*.x1;
  fdata.{2} <- l1  -. x1 +. lb.{0};
  fdata.{3} <- l1' -. x1 +. ub.{0};
  fdata.{4} <- l2  -. x2 +. lb.{1};
  fdata.{5} <- l2' -. x2 +. ub.{1}

(* Initial guesses *)
let set_initial_guess1 udata =
  (* There are two known solutions for this problem *)
  (* this init. guess should take us to (0.29945; 2.83693) *)
  let x1 = lb.{0} in
  let x2 = lb.{1} in
  udata.{0} <- x1;
  udata.{1} <- x2;
  udata.{2} <- x1 -. lb.{0};
  udata.{3} <- x1 -. ub.{0};
  udata.{4} <- x2 -. lb.{1};
  udata.{5} <- x2 -. ub.{1}

let set_initial_guess2 udata =
  (* There are two known solutions for this problem *)
  (* this init. guess should take us to (0.5; 3.1415926) *)
  let x1 = pt5 *. (lb.{0} +. ub.{0}) in
  let x2 = pt5 *. (lb.{1} +. ub.{1}) in

  udata.{0} <- x1;
  udata.{1} <- x2;
  udata.{2} <- x1 -. lb.{0};
  udata.{3} <- x1 -. ub.{0};
  udata.{4} <- x2 -. lb.{1};
  udata.{5} <- x2 -. ub.{1}

(* Print first lines of output (problem description) *)
let print_header fnormtol scsteptol =
  printf "\nFerraris and Tronconi test problem\n";
  printf "Tolerance parameters:\n";
  printf "  fnormtol  = %10.6g\n  scsteptol = %10.6g\n" fnormtol scsteptol

(* Print solution *)
let print_output u = printf " %8.6g  %8.6g\n" (ith u 1) (ith u 2)

(* Print final statistics contained in iopt *)
let print_final_stats kmem =
  let nni  = Kinsol.get_num_nonlin_solv_iters kmem in
  let nfe  = Kinsol.get_num_func_evals kmem in
  let nje  = Kinsol.Dls.get_num_jac_evals kmem in
  let nfeD = Kinsol.Dls.get_num_func_evals kmem in
  printf "Final Statistics:\n";
  printf "  nni = %5d    nfe  = %5d \n" nni nfe;
  printf "  nje = %5d    nfeD = %5d \n" nje nfeD

(* MAIN PROGRAM *)
let solve_it kmem u s glstr mset =
  printf("\n");
  if mset==1 then printf "Exact Newton" else printf "Modified Newton";
  if not glstr then printf "\n" else printf " with line search\n";
  Kinsol.set_max_setup_calls kmem (Some mset);
  ignore (Kinsol.solve kmem u glstr s s);
  printf "Solution:\n  [x1,x2] = ";
  print_output (Sundials.unvec u);
  print_final_stats kmem

let main () =
  lb.{0} <- pt25;
  lb.{1} <- onept5;
  ub.{0} <- one;
  ub.{1} <- two*.pi;

  (* Create serial vectors of length NEQ *)
  let u1 = RealArray.make neq in
  let u2 = RealArray.make neq in
  let u  = RealArray.make neq in
  let u_nvec = Nvector_serial.wrap u in
  let s_nvec = Nvector_serial.make neq one in (* no scaling *)
  let c = RealArray.of_list [
     zero; (* no constraint on x1 *)
     zero; (* no constraint on x2 *)
      one; (* l1 = x1 - x1_min >= 0 *)
    -.one; (* L1 = x1 - x1_max <= 0 *)
      one; (* l2 = x2 - x2_min >= 0 *)
    -.one; (* L2 = x2 - x22_min <= 0 *)
  ] in
  let c_nvec = Nvector_serial.wrap c in
  set_initial_guess1 u1;
  set_initial_guess2 u2;

  let fnormtol = ftol in
  let scsteptol = stol in

  (* Call KINDense to specify the linear solver *)
  let kmem = Kinsol.init (Kinsol.Dls.dense None) func u_nvec in
  Kinsol.set_constraints kmem c_nvec;
  Kinsol.set_func_norm_tol kmem (Some fnormtol);
  Kinsol.set_scaled_step_tol kmem (Some scsteptol);

  (* Print out the problem size, solution parameters, initial guess. *)
  print_header fnormtol scsteptol;

  (* --------------------------- *)

  printf "\n------------------------------------------\n";
  printf "\nInitial guess on lower bounds\n";
  printf "  [x1,x2] = ";
  print_output u1;

  RealArray.blit u1 u;
  solve_it kmem u_nvec s_nvec false 1;

  (* --------------------------- *)

  RealArray.blit u1 u;
  solve_it kmem u_nvec s_nvec true 1;

  (* --------------------------- *)

  RealArray.blit u1 u;
  solve_it kmem u_nvec s_nvec false 0;

  (* --------------------------- *)

  RealArray.blit u1 u;
  solve_it kmem u_nvec s_nvec true 0;

  (* --------------------------- *)

  printf "\n------------------------------------------\n";
  printf "\nInitial guess in middle of feasible region\n";
  printf "  [x1,x2] = ";
  print_output u2;

  RealArray.blit u2 u;
  solve_it kmem u_nvec s_nvec false 1;

  (* --------------------------- *)

  RealArray.blit u2 u;
  solve_it kmem u_nvec s_nvec true 1;

  (* --------------------------- *)

  RealArray.blit u2 u;
  solve_it kmem u_nvec s_nvec false 0;

  (* --------------------------- *)

  RealArray.blit u2 u;
  solve_it kmem u_nvec s_nvec true 0

let _ = main ()
let _ = Gc.compact ()

