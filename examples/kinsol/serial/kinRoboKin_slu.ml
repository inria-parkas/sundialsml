(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2010/12/01 23:08:49 $
 * -----------------------------------------------------------------
 * Programmer(s): Carol S. Woodward @ LLNL.  Adapted from the file
 *    kinRoboKin_dns.c by Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, December 2015.
 * -----------------------------------------------------------------
 * This example solves a nonlinear system from robot kinematics.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 6 from Section 14.1, Chapter 14
 *
 * The nonlinear system is solved by KINSOL using the SUPERLU_MT linear
 * solver.
 *
 * Constraints are imposed to make all components of the solution
 * be within .{-1,1}.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray

let printf = Printf.printf

(* Problem Constants *)

let nvar  = 8              (* variables *)
let neq   = 3*nvar         (* equations + bounds *)

let ftol  = 1.e-5 (* function tolerance *)
let stol  = 1.e-5 (* step tolerance *)

let zero  = 0.0
let one   = 1.0
let two   = 2.0

let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

(* System function *)
let func (yd : RealArray.t) (fd : RealArray.t) =
  let x1 = yd.{0} and l1 = yd.{ 8} and u1 = yd.{16}
  and x2 = yd.{1} and l2 = yd.{ 9} and u2 = yd.{17}
  and x3 = yd.{2} and l3 = yd.{10} and u3 = yd.{18}
  and x4 = yd.{3} and l4 = yd.{11} and u4 = yd.{19}
  and x5 = yd.{4} and l5 = yd.{12} and u5 = yd.{20}
  and x6 = yd.{5} and l6 = yd.{13} and u6 = yd.{21}
  and x7 = yd.{6} and l7 = yd.{14} and u7 = yd.{22}
  and x8 = yd.{7} and l8 = yd.{15} and u8 = yd.{23} in

  (* Nonlinear equations *)
  let eq1 = - 0.1238*.x1 +. x7 -. 0.001637*.x2
    -. 0.9338*.x4 +. 0.004731*.x1*.x3 -. 0.3578*.x2*.x3 -. 0.3571 in
  let eq2 = 0.2638*.x1 -. x7 -. 0.07745*.x2
    -. 0.6734*.x4 +. 0.2238*.x1*.x3 +. 0.7623*.x2*.x3 -. 0.6022 in
  let eq3 = 0.3578*.x1 +. 0.004731*.x2 +. x6*.x8 in
  let eq4 = -. 0.7623*.x1 +. 0.2238*.x2 +. 0.3461 in
  let eq5 = x1*.x1 +. x2*.x2 -. 1.0 in
  let eq6 = x3*.x3 +. x4*.x4 -. 1.0 in
  let eq7 = x5*.x5 +. x6*.x6 -. 1.0 in
  let eq8 = x7*.x7 +. x8*.x8 -. 1.0 in

  (* Lower bounds ( l_i = 1 + x_i >= 0)*)
  let lb1 = l1 -. 1.0 -. x1 in
  let lb2 = l2 -. 1.0 -. x2 in
  let lb3 = l3 -. 1.0 -. x3 in
  let lb4 = l4 -. 1.0 -. x4 in
  let lb5 = l5 -. 1.0 -. x5 in
  let lb6 = l6 -. 1.0 -. x6 in
  let lb7 = l7 -. 1.0 -. x7 in
  let lb8 = l8 -. 1.0 -. x8 in

  (* Upper bounds ( u_i = 1 - x_i >= 0)*)
  let ub1 = u1 -. 1.0 +. x1 in
  let ub2 = u2 -. 1.0 +. x2 in
  let ub3 = u3 -. 1.0 +. x3 in
  let ub4 = u4 -. 1.0 +. x4 in
  let ub5 = u5 -. 1.0 +. x5 in
  let ub6 = u6 -. 1.0 +. x6 in
  let ub7 = u7 -. 1.0 +. x7 in
  let ub8 = u8 -. 1.0 +. x8 in

  fd.{0} <- eq1; fd.{ 8} <- lb1; fd.{16} <- ub1;
  fd.{1} <- eq2; fd.{ 9} <- lb2; fd.{17} <- ub2;
  fd.{2} <- eq3; fd.{10} <- lb3; fd.{18} <- ub3;
  fd.{3} <- eq4; fd.{11} <- lb4; fd.{19} <- ub4;
  fd.{4} <- eq5; fd.{12} <- lb5; fd.{20} <- ub5;
  fd.{5} <- eq6; fd.{13} <- lb6; fd.{21} <- ub6;
  fd.{6} <- eq7; fd.{14} <- lb7; fd.{22} <- ub7;
  fd.{7} <- eq8; fd.{15} <- lb8; fd.{23} <- ub8

(* System Jacobian *)
let jac { Kinsol.jac_u = (yd : RealArray.t); Kinsol.jac_fu  = f } smat =
  let set_col = Matrix.Sparse.set_col smat in
  let set = Matrix.Sparse.set smat in

  let x1 = yd.{0}
  and x2 = yd.{1}
  and x3 = yd.{2}
  and x4 = yd.{3}
  and x5 = yd.{4}
  and x6 = yd.{5}
  and x7 = yd.{6}
  and x8 = yd.{7} in

  set_col 0 0;
  set_col 1 7;
  set_col 2 14;
  set_col 3 19;
  set_col 4 24;
  set_col 5 27;
  set_col 6 31;
  set_col 7 36;
  set_col 8 40;
  set_col 9 41;
  set_col 10 42;
  set_col 11 43;
  set_col 12 44;
  set_col 13 45;
  set_col 14 46;
  set_col 15 47;
  set_col 16 48;
  set_col 17 49;
  set_col 18 50;
  set_col 19 51;
  set_col 20 52;
  set_col 21 53;
  set_col 22 54;
  set_col 23 55;
  set_col 24 56;

  (* Nonlinear equations *)

  (*
     - 0.1238*x1 + x7 - 0.001637*x2
     - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571
  *)
  set  0 0 (-. 0.1238 +. 0.004731*.x3);
  set  7 0 (-. 0.001637 -. 0.3578*.x3);
  set 14 0 (0.004731*.x1 -. 0.3578*.x2);
  set 19 0 (-. 0.9338);
  set 31 0 (1.0);

  (*
    0.2638*x1 - x7 - 0.07745*x2
    - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022
  *)
  set  1 1 (0.2638 +. 0.2238*.x3);
  set  8 1 (-. 0.07745 +. 0.7623*.x3);
  set 15 1 (0.2238*.x1 +. 0.7623*.x2);
  set 20 1 (-. 0.6734);
  set 32 1 (-.1.0);

  (*
    0.3578*x1 + 0.004731*x2 + x6*x8
  *)
  set  2 2 (0.3578);
  set  9 2 (0.004731);
  set 27 2 (x8);
  set 36 2 (x6);

  (*
    - 0.7623*x1 + 0.2238*x2 + 0.3461
  *)
  set  3 3 (-. 0.7623);
  set 10 3 (0.2238);

  (*
    x1*x1 + x2*x2 - 1
  *)
  set  4 4 (2.0*.x1);
  set 11 4 (2.0*.x2);

  (*
    x3*x3 + x4*x4 - 1
  *)
  set 16 5 (2.0*.x3);
  set 21 5 (2.0*.x4);

  (*
    x5*x5 + x6*x6 - 1
  *)
  set 24 6 (2.0*.x5);
  set 28 6 (2.0*.x6);

  (*
    x7*x7 + x8*x8 - 1
  *)
  set 33 7 (2.0*.x7);
  set 37 7 (2.0*.x8);

  (*
    Lower bounds ( l_i = 1 + x_i >= 0)
    l_i - 1.0 - x_i
   *)
  set  5  8 (-1.0);
  set 12  9 (-1.0);
  set 17 10 (-1.0);
  set 22 11 (-1.0);
  set 25 12 (-1.0);
  set 29 13 (-1.0);
  set 34 14 (-1.0);
  set 38 15 (-1.0);

  set 40  8 (1.0);
  set 41  9 (1.0);
  set 42 10 (1.0);
  set 43 11 (1.0);
  set 44 12 (1.0);
  set 45 13 (1.0);
  set 46 14 (1.0);
  set 47 15 (1.0);

  (*
    Upper bounds ( u_i = 1 - x_i >= 0)
    u_i - 1.0 + x_i
   *)
  set  6 16 (1.0);
  set 13 17 (1.0);
  set 18 18 (1.0);
  set 23 19 (1.0);
  set 26 20 (1.0);
  set 30 21 (1.0);
  set 35 22 (1.0);
  set 39 23 (1.0);

  set 48 16 (1.0);
  set 49 17 (1.0);
  set 50 18 (1.0);
  set 51 19 (1.0);
  set 52 20 (1.0);
  set 53 21 (1.0);
  set 54 22 (1.0);
  set 55 23 (1.0)

(* Print solution *)
let print_output y =
  printf "     l=x+1          x         u=1-x\n";
  printf "   ----------------------------------\n";

  for i=1 to nvar do
    printf " %10.6g   %10.6g   %10.6g\n"
      (ith y (i+nvar))
      (ith y i)
      (ith y (i+2*nvar))
  done

(* Print final statistics *)
let print_final_stats kmem =
  let nni = Kinsol.get_num_nonlin_solv_iters kmem in
  let nfe = Kinsol.get_num_func_evals kmem in
  let nje = Kinsol.Dls.get_num_jac_evals kmem in
  printf "\nFinal Statistics.. \n";
  printf "nni    = %5d    nfe   = %5d \n" nni nfe;
  printf "nje    = %5d \n" nje

(* MAIN PROGRAM *)
let main () =
  printf "\nRobot Kinematics Example\n";
  printf "8 variables; -1 <= x_i <= 1\n";
  printf "KINSOL problem size: 8 + 2*8 = 24 \n\n";

  (* Create vectors for solution, scales, and constraints *)
  let y = Nvector_serial.make neq one in
  let ydata = Nvector.unwrap y in
  for i = 1 to nvar do
    set_ith ydata i (sqrt(two) /. two)
  done;
  let scale = Nvector_serial.make neq one in

  (* Initialize and allocate memory for KINSOL *)
  (* Attach dense linear solver *)
  let m = Matrix.sparse_csc ~nnz:56 neq in
  let kmem = Kinsol.(init
              ~linsolv:Dls.(solver Direct.(superlumt ~nthreads:2 y m) ~jac m)
              func y)
  in

  (* Set optional inputs *)
  let constraints = RealArray.make neq zero in
  for i=nvar+1 to neq do
    set_ith constraints i one
  done;
  Kinsol.set_constraints kmem (Nvector_serial.wrap constraints);
  Kinsol.set_func_norm_tol kmem ftol;
  Kinsol.set_scaled_step_tol kmem stol;

  (* Indicate exact Newton *)
  Kinsol.set_max_setup_calls kmem 1;

  (* Initial guess *)
  printf "Initial guess:\n";
  print_output ydata;

  (* Call KINSol to solve problem *)
  ignore Kinsol.(solve
                    kmem        (* KINSol memory block *)
                    y           (* initial guess on input; solution vector *)
                    LineSearch  (* global strategy choice *)
                    scale       (* scaling vector, for the variable cc *)
                    scale);     (* scaling vector for function values fval *)

  printf "\nComputed solution:\n";
  print_output ydata;

  (* Print final statistics and free memory *)
  print_final_stats kmem

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
