(*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:41 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2014.
 * -----------------------------------------------------------------
 * This example solves a 2D elliptic PDE
 *
 *    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0
 *
 * subject to homogeneous Dirichelt boundary conditions.
 * The PDE is discretized on a uniform NX+2 by NY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an system of size NEQ = NX*NY.
 * The nonlinear system is solved by KINSOL using the BAND linear
 * solver.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
let unvec = Sundials.unvec

let printf = Printf.printf
let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

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
   The vdata array is obtained via the macro call vdata = NV_DATA_S(v),
   where v is an N_Vector. 
   The variables are ordered by the y index j, then by the x index i. *)
let ijth v i j = v.{(j - 1) + (i - 1)*ny}
let set_ijth v i j e = v.{(j - 1) + (i - 1)*ny} <- e

(* System function *)
let func u f =
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

(* Print solution at selected points *)
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

(* Print final statistics *)
let print_final_stats kmem =
  (* Main solver statistics *)
  let nni = Kinsol.get_num_nonlin_solv_iters kmem in
  let nfe = Kinsol.get_num_func_evals kmem in

  (* Linesearch statistics *)
  let nbcfails = Kinsol.get_num_beta_cond_fails kmem in
  let nbacktr = Kinsol.get_num_backtrack_ops kmem in

  (* Main solver workspace size *)
  let lenrw, leniw = Kinsol.get_work_space kmem in

  (* Band linear solver statistics *)
  let nje = Kinsol.Dls.get_num_jac_evals kmem in
  let nfeD = Kinsol.Dls.get_num_func_evals kmem in

  (* Band linear solver workspace size *)
  let lenrwB, leniwB = Kinsol.Dls.get_work_space kmem in

  printf "\nFinal Statistics.. \n\n";
  printf "nni      = %6d    nfe     = %6d \n" nni nfe;
  printf "nbcfails = %6d    nbacktr = %6d \n" nbcfails nbacktr;
  printf "nje      = %6d    nfeB    = %6d \n" nje nfeD;
  printf "\n";
  printf "lenrw    = %6d    leniw   = %6d \n" lenrw leniw;
  printf "lenrwB   = %6d    leniwB  = %6d \n" lenrwB leniwB

(* MAIN PROGRAM *)
let main () =
  (* -------------------------
   * Print problem description
   * ------------------------- *)
  printf "\n2D elliptic PDE on unit square\n";
  printf "   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0\n";
  printf " + homogeneous Dirichlet boundary conditions\n\n";
  printf "Solution method: Modified Newton with band linear solver\n";
  printf "Problem size: %2d x %2d = %4d\n" nx ny neq;

  (* -------------
   * Initial guess 
   * ------------- *)
  let y = Nvector_serial.make neq zero in

  (* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * y is used as a template
   * Attach band linear solver 
   * ----------------------------------------- *)
  let kmem = Kinsol.init
                (Kinsol.Dls.band {Kinsol.mupper=ny; Kinsol.mlower=ny} None)
                func y
  in
  (* -------------------
   * Set optional inputs 
   * ------------------- *)

  (* Specify stopping tolerance based on residual *)
  Kinsol.set_func_norm_tol kmem (Some ftol);

  (* ------------------------------
   * Parameters for Modified Newton
   * ------------------------------ *)

  (* Force a Jacobian re-evaluation every mset iterations *)
  Kinsol.set_max_setup_calls kmem (Some 100);

  (* Every msubset iterations, test if a Jacobian evaluation
     is necessary *)
  Kinsol.set_max_sub_setup_calls kmem (Some 1);

  (* ----------------------------
   * Call KINSol to solve problem 
   * ---------------------------- *)

  (* No scaling used *)
  let scale = Nvector_serial.make neq one in

  (* Call main solver *)
  ignore (Kinsol.solve
            kmem           (* KINSol memory block *)
            y              (* initial guess on input; solution vector *)
            true           (* global strategy choice *)
            scale          (* scaling vector, for the variable cc *)
            scale);        (* scaling vector for function values fval *)


  (* ------------------------------------
   * Print solution and solver statistics 
   * ------------------------------------ *)

  (* Get scaled norm of the system function *)
  let fnorm = Kinsol.get_func_norm kmem in
  printf "\nComputed solution (||F|| = %g):\n\n" fnorm;

  print_output (unvec y);

  print_final_stats kmem

let _ = main ()
let _ = Gc.compact ()

