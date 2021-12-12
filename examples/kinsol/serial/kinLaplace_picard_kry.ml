(* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * Based on kinLaplace_picard_bnd.c by Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, November 2021.
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This example solves a 2D elliptic PDE
 *
 *    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
 *
 * subject to homogeneous Dirichlet boundary conditions.
 * The PDE is discretized on a uniform NX+2 by NY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving a system of size NEQ = NX*NY.
 * The nonlinear system is solved by KINSOL using the Picard
 * iteration and the SPGMR linear solver.
 * -----------------------------------------------------------------
 *)

open Sundials

let unwrap = Nvector.unwrap
let printf = Printf.printf

(* Problem Constants *)

let nx   = 31             (* no. of points in x direction *)
let ny   = 31             (* no. of points in y direction *)
let neq  = nx*ny          (* problem dimension *)

let skip = 3              (* no. of points skipped for printing *)

let ftol = 1.e-12 (* function tolerance *)

let zero = 0.0
let one  = 1.0
let two  = 2.0

(* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage.
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= NX, 1 <= j <= NY.
   The vdata array is obtained via the call vdata = N_VGetArrayPointer(v),
   where v is an N_Vector.
   The variables are ordered by the y index j, then by the x index i. *)
let ijth (v : RealArray.t) i j = v.{(j - 1) + (i - 1)*ny}
let set_ijth (v : RealArray.t) i j e = v.{(j - 1) + (i - 1)*ny} <- e

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(*
 * System function
 *)

let func (u : RealArray.t) (f : RealArray.t) =
  let dx = one/.float (nx+1) in
  let dy = one/.float (ny+1) in
  let hdc = one/.(dx*.dx) in
  let vdc = one/.(dy*.dy) in

  for j=1 to ny do
    for i=1 to nx do
      (* Extract u at x_i, y_j and four neighboring points *)
      let uij = ijth u i j in
      let udn = if j = 1  then zero else ijth u i (j-1) in
      let uup = if j = ny then zero else ijth u i (j+1) in
      let ult = if i = 1  then zero else ijth u (i-1) j in
      let urt = if i = nx then zero else ijth u (i+1) j in

      (* Evaluate diffusion components *)
      let hdiff = hdc*.(ult -. two*.uij +. urt) in
      let vdiff = vdc*.(uup -. two*.uij +. udn) in

      (* Set residual at x_i, y_j *)
      set_ijth f i j (hdiff +. vdiff +. uij -. uij*.uij*.uij +. 2.0)
    done
  done

(*
 * Jacobian vector product function
 *)

let jactimes (v : RealArray.t) (jv : RealArray.t) _ _ =
  let dx  = one /. float (nx+1) in
  let dy  = one /. float (ny+1) in
  let hdc = one /. (dx *. dx) in
  let vdc = one /. (dy *. dy) in

  for j=1 to ny do
    for i=1 to nx do

      (* Extract v at x_i, y_j and four neighboring points *)
      let vij = ijth v i j in
      let vdn = if j = 1  then zero else ijth v i (j-1) in
      let vup = if j = ny then zero else ijth v i (j+1) in
      let vlt = if i = 1  then zero else ijth v (i-1) j in
      let vrt = if i = nx then zero else ijth v (i+1) j in

      (* Evaluate diffusion components *)

      let hdiff = hdc *. (vlt -. two *. vij +. vrt) in
      let vdiff = vdc *. (vup -. two *. vij +. vdn) in

      (* Set Jv at x_i, y_j *)
      set_ijth jv i j (hdiff +. vdiff)
    done
  done;
  false

(*
 * Print solution at selected points
 *)

let prloop max f =
  let rec go i =
    if i <= max then (f i; go (i + skip)) else ()
  in go 1

let print_output u =
  let dx = one/.float(nx+1) in
  let dy = one/.float(ny+1) in
  printf "            ";

  prloop nx (fun i -> printf "%-8.5f " (float i *. dx));
  printf("\n\n");

  prloop ny (fun j ->
    printf "%-8.5f    " (float j *. dy);
    prloop nx (fun i -> printf "%-8.5f " (ijth u i j));
    printf("\n")
  )

(*
 * Print final statistics
 *)

let print_final_stats kmem =
  let open Kinsol in
  (* Main solver statistics *)
  let nni = get_num_nonlin_solv_iters kmem in
  let nfe = get_num_func_evals kmem in

  (* Linear solver statistics *)
  let nli = Spils.get_num_lin_iters kmem in
  let nfeLS = Spils.get_num_lin_func_evals kmem in
  let ncfl = Spils.get_num_lin_conv_fails kmem in
  let njvevals = Spils.get_num_jtimes_evals kmem in
  let npe = Spils.get_num_prec_evals kmem in
  let nps = Spils.get_num_prec_solves kmem in

  (* Main solver workspace size *)
  let lenrw, leniw = get_work_space kmem in

  (* Linear solver workspace size *)
  let lenrwLS, leniwLS = Spils.get_work_space kmem in

  printf "\nFinal Statistics.. \n\n";
  printf "nni = %6d  nli   = %6d  ncfl = %6d\n" nni nli ncfl;
  printf "nfe = %6d  nfeLS = %6d  njt  = %6d\n" nfe nfeLS njvevals;
  printf "npe = %6d  nps   = %6d\n" npe nps;
  printf "\n";
  printf "lenrw   = %6d  leniw   = %6d\n" lenrw leniw;
  printf "lenrwLS = %6d  leniwLS = %6d\n" lenrwLS leniwLS

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* -------------------------
   * Print problem description
   * ------------------------- *)

  printf "\n2D elliptic PDE on unit square\n";
  printf "   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0\n";
  printf " + homogeneous Dirichlet boundary conditions\n\n";
  printf "Solution method: Anderson accelerated Picard iteration with SPGMR linear solver.\n";
  printf "Problem size: %2d x %2d = %4d\n" nx ny neq;

  (* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- *)
  let y = Nvector_serial.make neq zero in

  (* No scaling used *)
  let scale = Nvector_serial.make neq one in

  (* ----------------------
   * Create SUNLinearSolver
   * ---------------------- *)
  let ls = LinearSolver.Iterative.spgmr ~maxl:10 y in

  (* ----------------------------------------------------------------------------------
   * Attach linear solver
   * Set Jacobian vector product function
   *
   * Initialize and allocate memory for KINSOL, set parametrs for Anderson acceleration
   * ---------------------------------------------------------------------------------- *)

  (* y is used as a template *)
  (* Use acceleration with up to 3 prior residuals *)
  let kmem = Kinsol.(init ~maa:3
                       ~lsolver:Spils.(solver ~jac_times_vec:jactimes
                                              ls prec_none)
                       func y)
  in

  (* -------------------
   * Set optional inputs
   * ------------------- *)

  (* Specify stopping tolerance based on residual *)

  let fnormtol  = ftol in
  Kinsol.set_func_norm_tol kmem fnormtol;

  (* Set information file *)
  let infofp = Logfile.openfile "KINSOL.log" in
  Kinsol.(set_info_file kmem ~print_level:ShowGlobalValues infofp);

  (* -------------
   * Initial guess
   * ------------- *)
  Nvector_serial.Ops.const zero y;
  set_ijth (Nvector.unwrap y) 2 2 one;

  (* ----------------------------
   * Call KINSol to solve problem
   * ---------------------------- *)

  (* Call main solver *)
  ignore (Kinsol.solve
            kmem              (* KINSol memory block *)
            y                 (* initial guess on input; solution vector *)
            Kinsol.Picard     (* global strategy choice *)
            scale             (* scaling vector, for the variable cc *)
            scale);           (* scaling vector for function values fval *)

  (* ------------------------------------
   * Print solution and solver statistics
   * ------------------------------------ *)

  (* Get scaled norm of the system function *)

  let fnorm = Kinsol.get_func_norm kmem in
  printf "\nComputed solution (||F|| = %g):\n\n" fnorm;

  print_output (Nvector.unwrap y);

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
  for _ = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
