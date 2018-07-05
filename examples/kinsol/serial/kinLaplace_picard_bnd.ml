(*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:41 $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2015.
 * -----------------------------------------------------------------
 * This example solves a 2D elliptic PDE
 *
 *    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0
 *
 * subject to homogeneous Dirichelt boundary conditions.
 * The PDE is discretized on a uniform NX+2 by NY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving a system of size NEQ = NX*NY.
 * The nonlinear system is solved by KINSOL using the Picard
 * iteration and the BAND linear solver.
 *
 * This file is strongly based on the kinLaplace_bnd.c file
 * developed by Radu Serban.
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
   The vdata array is obtained via the macro call vdata = NV_DATA_S(v),
   where v is an N_Vector.
   The variables are ordered by the y index j, then by the x index i. *)
let ijth (v : RealArray.t) i j = v.{(j - 1) + (i - 1)*ny}
let set_ijth (v : RealArray.t) i j e = v.{(j - 1) + (i - 1)*ny} <- e

let bandset bm i j v = Matrix.Band.set bm i j v

(* Jacobian function *)
let jac { Kinsol.jac_u = (yd : RealArray.t); Kinsol.jac_fu  = f } jac =

  let dx  = 1.0 /. float(nx+1) in
  let dy  = 1.0 /. float(ny+1) in
  let hdc = 1.0 /. (dx*.dx) in
  let vdc = 1.0 /. (dy*.dy) in

  (* The components of f(t,u) which depend on u_{i,j} are
     f_{i,j}, f_{i-1,j}, f_{i+1,j}, f_{i,j+1}, and f_{i,j-1}.
     Thus, a column of the Jacobian will contain an entry from
     each of these equations exception the ones on the boundary.

     f_{i,j}   = hdc*(u_{i-1,j}  -2u_{i,j}  +u_{i+1,j})
                                + vdc*(u_{i,j-1}  -2u_{i,j}  +u_{i,j+1})
     f_{i-1,j} = hdc*(u_{i-2,j}  -2u_{i-1,j}+u_{i,j})
                                + vdc*(u_{i-1,j-1}-2u_{i-1,j}+u_{i-1,j+1})
     f_{i+1,j} = hdc*(u_{i,j}    -2u_{i+1,j}+u_{i+2,j})
                                + vdc*(u_{i+1,j-1}-2u_{i+1,j}+u_{i+1,j+1})
     f_{i,j-1} = hdc*(u_{i-1,j-1}-2u_{i,j-1}+u_{i+1,j-1})
                                + vdc*(u_{i,j-2}  -2u_{i,j-1}+u_{i,j})
     f_{i,j+1} = hdc*(u_{i-1,j+1}-2u_{i,j+1}+u_{i+1,j+1})
                                + vdc*(u_{i,j}    -2u_{i,j+1}+u_{i,j+2}) *)
  for j = 0 to ny - 1 do
    for i =0 to nx - 1 do
      (* Evaluate diffusion coefficients *)
      let k = i + j*nx in
      bandset jac k k (-2.0*.hdc -. 2.0*.vdc);
      if ( i <> (nx-1) ) then bandset jac k (k+1)  hdc;
      if ( i <> 0 )      then bandset jac k (k-1)  hdc;
      if ( j <> (ny-1) ) then bandset jac k (k+nx) vdc;
      if ( j <> 0 )      then bandset jac k (k-nx) vdc
    done
  done

(* System function *)
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
  let open Kinsol in
  (* Main solver statistics *)
  let nni = get_num_nonlin_solv_iters kmem in
  let nfe = get_num_func_evals kmem in

  (* Band linear solver statistics *)
  let nje  = Dls.get_num_jac_evals kmem in
  let nfeD = Dls.get_num_func_evals kmem in

  (* Band linear solver workspace size *)
  let lenrwB, leniwB = Dls.get_work_space kmem in

  printf "\nFinal Statistics.. \n\n";
  printf "nni      = %6d    nfe     = %6d \n" nni nfe;
  printf "nje      = %6d    nfeB    = %6d \n" nje nfeD;
  printf "\n";
  printf "lenrwB   = %6d    leniwB  = %6d \n" lenrwB leniwB

(* MAIN PROGRAM *)
let main () =
  (* -------------------------
   * Print problem description
   * ------------------------- *)
  printf "\n2D elliptic PDE on unit square\n";
  printf "   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0\n";
  printf " + homogeneous Dirichlet boundary conditions\n\n";
  printf "Solution method: %s.\n"
         "Anderson accelerated Picard iteration with band linear solver";
  printf "Problem size: %2d x %2d = %4d\n" nx ny neq;

  (* -------------
   * Initial guess
   * ------------- *)
  let y = Nvector_serial.make neq zero in

  (* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * set parameters for Anderson acceleration
   * y is used as a template
   * Attach band linear solver
   * Use acceleration with up to 3 prior residuals
   * ----------------------------------------- *)
  let m = Matrix.band ~smu:(2*nx) ~mu:nx ~ml:nx neq in
  let kmem = Kinsol.(init ~maa:3
                          ~linsolv:Dls.(solver ~jac:jac (band y m))
                          func y)
  in

  (* -------------------
   * Set optional inputs
   * ------------------- *)

  (* Specify stopping tolerance based on residual *)
  Kinsol.set_func_norm_tol kmem ftol;

  (* ----------------------------
   * Call KINSol to solve problem
   * ---------------------------- *)

  (* No scaling used *)
  let scale = Nvector_serial.make neq one in
  set_ijth (unwrap y) 2 2 1.0;

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
