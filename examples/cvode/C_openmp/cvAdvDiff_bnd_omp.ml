(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/29 22:21:28 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * Programmer(s): Ting Yan @ SMU
 *     Based on cvAdvDiff_bnd.c and parallelized with OpenMP
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Sep 2010.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODE.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVBAND band linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Optionally, we can set the number of threads from environment
 * variable or command line. To check the current value for number
 * of threads from environment:
 *      % echo $OMP_NUM_THREADS
 *
 * Execution:
 *
 * If the user want to use the default value or the number of threads
 * from environment value:
 *      % ./cvAdvDiff_bnd_omp
 * If the user want to specify the number of threads to use
 *      % ./cvAdvDiff_bnd_omp num_threads
 * where num_threads is the number of threads the user want to use
 * -----------------------------------------------------------------
 *)

open Sundials

let unwrap = Nvector.unwrap

let printf = Printf.printf
let vmax_norm = Nvector_openmp.Ops.n_vmaxnorm

(* Header files with a description of contents used in cvbanx.c *)

let set bm i j v = Matrix.Band.set bm i j v

(* Problem Constants *)

let xmax   = 2.0        (* domain boundaries         *)
let ymax   = 1.0
let mx     = 10         (* mesh dimensions           *)
let my     = 5
let neq    = mx * my    (* number of equations       *)
let atol   = 1.0e-5     (* scalar absolute tolerance *)
let t0     = 0.0        (* initial time              *)
let t1     = 0.1        (* first output time         *)
let dtout  = 0.1        (* output time increment     *)
let nout   = 10         (* number of output times    *)

let zero  = 0.0
let half  = 0.5
let one   = 1.0
let two   = 2.0
let five  = 5.0

(* User-defined vector access macro IJth *)

(* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage.
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The vdata array is obtained via the macro call vdata = NV_DATA_S(v),
   where v is an N_Vector.
   The variables are ordered by the y index j, then by the x index i. *)
let set_ijth (vdata : RealArray.t) i j e = vdata.{(j-1) + (i-1)*my} <- e

(* Type : UserData (contains grid constants) *)

type user_data = {
  dx : float;
  dy : float;
  hdcoef : float;
  hacoef : float;
  vdcoef : float;
}

(* f routine. Compute f(t,u). *)

let f data t (udata : RealArray.t) (dudata : RealArray.t) =
  (* Extract needed constants from data *)
  let hordc = data.hdcoef
  and horac = data.hacoef
  and verdc = data.vdcoef
  in

  (* Loop over all grid points. *)
  for j = 1 to my do
    for i = 1 to mx do

      (* Extract u at x_i, y_j and four neighboring points *)
      let uij = udata.{j-1 + (i-1)*my};
      and udn = (if j = 1  then zero else udata.{(j-2) + (i-1)*my})
      and uup = (if j = my then zero else udata.{j     + (i-1)*my})
      and ult = (if i = 1  then zero else udata.{(j-1) + (i-2)*my})
      and urt = (if i = mx then zero else udata.{(j-1) + i*my})
      in

      (* Set diffusion and advection terms and load into udot *)
      let hdiff = hordc *. (ult -. two *. uij +. urt)
      and hadv  = horac *. (urt -. ult)
      and vdiff = verdc *. (uup -. two *. uij +. udn)
      in
      dudata.{(j-1)+(i-1)*my} <- hdiff +. hadv +. vdiff

    done
  done

(* Jacobian routine. Compute J(t,u). *)

let jac data arg jmat =
  (*
    The components of f = udot that depend on u(i,j) are
    f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
      df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
      df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
      df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
      df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
      df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
  *)

  let hordc = data.hdcoef
  and horac = data.hacoef
  and verdc = data.vdcoef
  in

  (* set non-zero Jacobian etnries *)
  for j = 1 to my do
    for i = 1 to mx do
      (* set the kth column of jmat *)
      let k = j - 1 + (i - 1) * my in
      set jmat k k (-. two *. (verdc +. hordc));
      if (i <> 1)  then set jmat (k - my) k (hordc +. horac);
      if (i <> mx) then set jmat (k + my) k (hordc -. horac);
      if (j <> 1)  then set jmat (k - 1)  k verdc;
      if (j <> my) then set jmat (k + 1)  k verdc
    done
  done

(* Set initial conditions in u vector *)

let set_ic u data =
  (* Extract needed constants from data *)
  let dx = data.dx
  and dy = data.dy
  in

  (* Load initial profile into u vector *)
  for j = 1 to my do
    let y = float(j) *. dy in
    for i = 1 to mx do
      let x = float(i) *. dx in
      set_ijth u i j
        (x *. (xmax -. x) *. y *. (ymax -. y) *. exp(five *. x *. y))
    done
  done

(* Print first lines of output (problem description) *)

let print_header reltol abstol umax =
  printf "\n2-D Advection-Diffusion Equation\n";
  printf "Mesh dimensions = %d X %d\n" mx my;
  printf "Total system size = %d\n" neq;

  printf "Tolerance parameters: reltol = %g   abstol = %g\n\n" reltol abstol;
  printf "At t = %g      max.norm(u) =%14.6e \n" t0 umax

(* Print current value *)

let print_output = printf "At t = %4.2f   max.norm(u) =%14.6e   nst = %4d\n"

let print_final_stats s =
  let open Cvode in
  let nst     = get_num_steps s
  and nfe     = get_num_rhs_evals s
  and nsetups = get_num_lin_solv_setups s
  and netf    = get_num_err_test_fails s
  and nni     = get_num_nonlin_solv_iters s
  and ncfn    = get_num_nonlin_solv_conv_fails s
  and nje     = Dls.get_num_jac_evals s
  and nfeLS   = Dls.get_num_lin_rhs_evals s
  in
  printf "\nFinal Statistics:\n";
  printf "nst = %-6d nfe  = %-6d nsetups = %-6d nfeLS = %-6d nje = %d\n"
  nst nfe nsetups nfeLS nje;
  printf "nni = %-6d ncfn = %-6d netf = %d\n"
  nni ncfn netf

let main () =
  (* Set the number of threads to use *)
  let num_threads =
    if Array.length Sys.argv > 1
    then int_of_string Sys.argv.(1)
    else 1
  in

  (* Create a serial vector *)
  let u = Nvector_openmp.make num_threads neq 0.0 in (* Allocate u vector *)

  let reltol = zero  (* Set the tolerances *)
  and abstol = atol
  in

  let dx = xmax /. float(mx + 1); (* Set grid coefficients in data *)
  and dy = ymax /. float(my + 1);
  in
  let data = {
    dx = dx; dy = dy;
    hdcoef = one /. (dx *. dx);
    hacoef = half /. (two *. dx);
    vdcoef = one /. (dy *. dy);
  } in

  set_ic (unwrap u) data;  (* Initialize u vector *)

  (* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration *)
  (* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. *)
  (* Call CVBand to specify the CVBAND band linear solver *)
  (* Set the user-supplied Jacobian routine Jac *)
  let m = Matrix.band ~smu:(2*my) ~mu:my ~ml:my neq in
  let cvode_mem = Cvode.(init BDF
                    (SStolerances (reltol, abstol))
                    ~lsolver:(Dls.(solver ~jac:(jac data) (band u m)))
                    (f data) t0 u)
  in

  (* In loop over output points: call CVode, print results, test for errors *)

  print_header reltol abstol (vmax_norm u);

  let tout = ref t1 in
  for iout = 1 to nout do
    let (t, flag) = Cvode.solve_normal cvode_mem !tout u
    in
    let nst = Cvode.get_num_steps cvode_mem in

    print_output t (vmax_norm u) nst;
    tout := !tout +. dtout
  done;

  print_final_stats cvode_mem;  (* Print some final statistics   *)
  printf "num_threads = %d\n\n" num_threads

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
