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

(* System function for predator-prey system *)
let lb = RealArray.create nvar
let ub = RealArray.create nvar

let jac { Kinsol.jac_u = (y : RealArray.t) } jacmat =
  let set_row = Sls.SparseMatrix.set_row jacmat in
  let set = Sls.SparseMatrix.set jacmat in
  Sls.SparseMatrix.set_to_zero jacmat;

  set_row 0  0;
  set_row 1  2;
  set_row 2  4;
  set_row 3  6;
  set_row 4  8;
  set_row 5 10;
  set_row 6 12;

  set  0 0 (pt5 *. cos(y.{0}*.y.{1}) *. y.{1} -. pt5);
  set  1 1 (pt5 *. cos(y.{0}*.y.{1}) *. y.{0} -. pt25/.pi);

  set  2 0 (two *. (one -. pt25/.pi) *. (exp(two*.y.{0}) -. e));
  set  3 1 (e /. pi);

  set  4 0 (-.one);
  set  5 2 (one);
 
  set  6 0 (-.one);
  set  7 3 (one);

  set  8 1 (-.one);
  set  9 4 (one);

  set 10 1 (-.one);
  set 11 5 (one)

let func (udata : RealArray.t) (fdata : RealArray.t) =
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
let set_initial_guess1 (udata : RealArray.t) =
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

let set_initial_guess2 (udata : RealArray.t) =
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
  print_string "\nFerraris and Tronconi test problem\n";
  print_string "Tolerance parameters:\n";
  printf "  fnormtol  = %10.6g\n  scsteptol = %10.6g\n" fnormtol scsteptol

(* Print solution *)
let print_output (u : RealArray.t) =
  printf " %8.6g  %8.6g\n" u.{0} u.{1}

(* Print final statistics contained in iopt *)
(* For high NUM_REPS, the cost of OCaml printf becomes important! *)
let print_final_stats kmem =
  let open Kinsol in
  let nni  = get_num_nonlin_solv_iters kmem in
  let nfe  = get_num_func_evals kmem in
  let nje  = Sls.Klu.get_num_jac_evals kmem in
  print_string "Final Statistics:\n";
  printf "  nni = %5d    nfe  = %5d \n  nje = %5d    \n" nni nfe nje

(* MAIN PROGRAM *)
let solve_it kmem u s glstr mset =
  print_newline ();
  print_string (if mset==1 then "Exact Newton" else "Modified Newton");
  if glstr then print_string " with line search\n" else print_newline ();
  Kinsol.set_max_setup_calls kmem mset;
  ignore Kinsol.(solve kmem u (if glstr then LineSearch else Newton) s s);
  print_string "Solution:\n  [x1,x2] = ";
  print_output (Nvector.unwrap u);
  print_final_stats kmem

let main () =
  lb.{0} <- pt25;
  lb.{1} <- onept5;
  ub.{0} <- one;
  ub.{1} <- two*.pi;

  (* Create serial vectors of length NEQ *)
  let u1 = RealArray.create neq in
  let u2 = RealArray.create neq in
  let u  = RealArray.create neq in
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
  let nnz = 12 in

  (* Call KINKlu to specify the linear solver *)
  let kmem = Kinsol.(init ~lsolver:(Sls.Klu.solver_csr jac nnz) func u_nvec) in
  Kinsol.set_constraints kmem c_nvec;
  Kinsol.set_func_norm_tol kmem fnormtol;
  Kinsol.set_scaled_step_tol kmem scsteptol;

  (* Print out the problem size, solution parameters, initial guess. *)
  print_header fnormtol scsteptol;

  (* --------------------------- *)

  print_string "\n------------------------------------------\n";
  print_string "\nInitial guess on lower bounds\n";
  print_string "  [x1,x2] = ";
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

  print_string "\n------------------------------------------\n";
  print_string "\nInitial guess in middle of feasible region\n";
  print_string "  [x1,x2] = ";
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

(* Check environment variables for extra arguments.  *)
let reps =
  try int_of_string (Unix.getenv "NUM_REPS")
  with Not_found | Failure _ -> 1
let gc_at_end =
  try int_of_string (Unix.getenv "GC_AT_END") <> 0
  with Not_found | Failure _ -> false
let gc_each_rep =
  try int_of_string (Unix.getenv "GC_EACH_REP") <> 0
  with Not_found | Failure _ -> false

(* Entry point *)
let _ =
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
