(*
 * -----------------------------------------------------------------
 * $Revision: 4074 $
 * $Date: 2014-04-23 14:13:52 -0700 (Wed, 23 Apr 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2015.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by the accelerated fixed point solver in
 * KINSOL.
 * The problem is from chemical kinetics, and consists of solving
 * the first time step in a Backward Euler solution for the
 * following three rate equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2
 *    dy3/dt = 3.e2*(y2)^2
 * on the interval from t = 0.0 to t = 0.1, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
let unwrap = Nvector.unwrap
let printf = Printf.printf

(* Problem Constants *)

let neq   = 3      (* number of equations  *)
let y10   = 1.0    (* initial y components *)
let y20   = 0.0
let y30   = 0.0
let tol   = 1.e-10 (* function tolerance *)
let dstep = 0.1    (* Size of the single time step used *)

let priors = 2

(* User-defined vector accessor macro: Ith *)

let ith (v : RealArray.t) i = v.{i - 1}
let set_ith (v : RealArray.t) i e = v.{i - 1} <- e

(* System function *)

let func_roberts (y : RealArray.t) (g : RealArray.t) =
  let y1 = ith y 1 in
  let y2 = ith y 2 in
  let y3 = ith y 3 in

  let yd1 = dstep *. ( -0.04*.y1 +. 1.0e4*.y2*.y3 ) in
  let yd3 = dstep *. 3.0e2*.y2*.y2 in

  set_ith g 1 (yd1 +. y10);
  set_ith g 2 (-.yd1 -. yd3 +. y20);
  set_ith g 3 (yd3 +. y30)

(* Print solution at selected points *)

let print_output (y : RealArray.t) =
  let y1 = ith y 1 in
  let y2 = ith y 2 in
  let y3 = ith y 3 in
  printf "y =%14.6e  %14.6e  %14.6e\n" y1 y2 y3

(* Print final statistics *)

let print_final_stats kmem =
  let nni = Kinsol.get_num_nonlin_solv_iters kmem in
  let nfe = Kinsol.get_num_func_evals kmem in
  printf "\nFinal Statistics.. \n\n";
  printf "nni      = %6d    nfe     = %6d \n" nni nfe

let main () =
  (* Print problem description *)

  printf "Example problem from chemical kinetics solving\n";
  printf "the first time step in a Backward Euler solution for the\n";
  printf "following three rate equations:\n";
  printf "    dy1/dt = -.04*y1 + 1.e4*y2*y3\n";
  printf "    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2\n";
  printf "    dy3/dt = 3.e2*(y2)^2\n";
  printf "on the interval from t = 0.0 to t = 0.1, with initial\n";
  printf "conditions: y1 = 1.0, y2 = y3 = 0.\n";
  printf "Solution method: Anderson accelerated fixed point iteration.\n";

  (* Create vectors for solution and scales *)
  let y = Nvector_serial.make neq 0.0 in

  (* No scaling used *)
  let scale = Nvector_serial.make neq 1.0 in

  (* Initialize and allocate memory for KINSOL *)

  (* y is used as a template *)
  (* Set number of prior residuals used in Anderson acceleration *)
  let kmem = Kinsol.init ~maa:priors func_roberts y in

  (* Set optional inputs *)

  (* Specify stopping tolerance based on residual *)
  let fnormtol  = tol in
  Kinsol.set_func_norm_tol kmem fnormtol;

  (* Initial guess *)
  set_ith (unwrap y) 1 1.0;

  (* Call KINSol to solve problem *)

  (* Call main solver *)
  ignore Kinsol.(solve
                    kmem       (* KINSol memory block *)
                    y          (* initial guess on input; solution vector *)
                    FixedPoint (* global strategy choice *)
                    scale      (* scaling vector, for the variable cc *)
                    scale);    (* scaling vector for function values fval *)

  (* Print solution and solver statistics *)

  (* Get scaled norm of the system function *)
  let fnorm = Kinsol.get_func_norm kmem in

  printf "\nComputed solution (||F|| = %g):\n\n" fnorm;
  print_output (unwrap y);

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

