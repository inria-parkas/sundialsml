(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2009/02/17 02:48:46 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
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
 * iteration with the LAPACK band linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module Roots = Sundials.Roots
module Col = Dls.BandMatrix.Col
module Dls = Cvode.Dls
let unvec = Sundials.unvec

let printf = Printf.printf
let vmax_norm = Nvector_serial.Ops.n_vmaxnorm

let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

(* Header files with a description of contents used in cvbanx.c *)


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

let ijth vdata i j = vdata.{(j-1) + (i-1)*my}
let set_ijth vdata i j e = vdata.{(j-1) + (i-1)*my} <- e

(* Type : UserData (contains grid constants) *)

type user_data = {
  dx : float;
  dy : float;
  hdcoef : float;
  hacoef : float;
  vdcoef : float;
}

(* f routine. Compute f(t,u). *)

let f data t udata dudata =
  (* Extract needed constants from data *)
  let hordc = data.hdcoef
  and horac = data.hacoef
  and verdc = data.vdcoef
  in

  (* Loop over all grid points. *)
  let ijth = ijth udata in

  for j = 1 to my do
    for i = 1 to mx do

      (* Extract u at x_i, y_j and four neighboring points *)
      let uij = ijth i j;
      and udn = (if j = 1  then zero else ijth i (j-1))
      and uup = (if j = my then zero else ijth i (j+1))
      and ult = (if i = 1  then zero else ijth (i-1) j)
      and urt = (if i = mx then zero else ijth (i+1) j)
      in

      (* Set diffusion and advection terms and load into udot *)
      let hdiff = hordc *. (ult -. two *. uij +. urt)
      and hadv  = horac *. (urt -. ult)
      and vdiff = verdc *. (uup -. two *. uij +. udn)
      in
      set_ijth dudata i j (hdiff +. hadv +. vdiff)

    done
  done

(* Jacobian routine. Compute J(t,u). *)

let jac data {Cvode.mupper=mupper; Cvode.mlower=mlower} arg jmat =
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
  
      let k = j - 1 + (i - 1) * my in
      let kthCol = Col.get_col jmat k in

      (* set the kth column of jmat *)
      Col.set kthCol k k (-. two *. (verdc +. hordc));
      if (i <> 1)  then Col.set kthCol (k - my) k (hordc +. horac);
      if (i <> mx) then Col.set kthCol (k + my) k (hordc -. horac);
      if (j <> 1)  then Col.set kthCol (k - 1)  k verdc;
      if (j <> my) then Col.set kthCol (k + 1)  k verdc
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
  let nst = Cvode.get_num_steps s
  and nfe = Cvode.get_num_rhs_evals s
  and nsetups = Cvode.get_num_lin_solv_setups s
  and netf = Cvode.get_num_err_test_fails s
  and nni = Cvode.get_num_nonlin_solv_iters s
  and ncfn = Cvode.get_num_nonlin_solv_conv_fails s
  and nje = Cvode.Dls.get_num_jac_evals s
  and nfeLS = Cvode.Dls.get_num_rhs_evals s
  in
  printf "\nFinal Statistics:\n";
  printf "nst = %-6d nfe  = %-6d nsetups = %-6d nfeLS = %-6d nje = %d\n"
  nst nfe nsetups nfeLS nje;
  printf "nni = %-6d ncfn = %-6d netf = %d\n \n"
  nni ncfn netf

let main () =
  (* Create a serial vector *)
  let u = Nvector_serial.make neq 0.0 in (* Allocate u vector *)

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

  set_ic (unvec u) data;  (* Initialize u vector *)

  (* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration *)
  (* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. *)
  (* Call CVLapackBand to specify the CVBAND band linear solver *)
  (* Set the user-supplied Jacobian routine Jac *)
  let solver = Cvode.Dls.lapack_band {Cvode.mupper = my; Cvode.mlower = my}
                              (Some (jac data))
  in
  let cvode_mem = Cvode.init Cvode.BDF (Cvode.Newton solver)
                             (Cvode.SStolerances (reltol, abstol))
                             ~t0:t0 (f data) u
  in
  Gc.compact ();

  (* In loop over output points: call CVode, print results, test for errors *)

  print_header reltol abstol (vmax_norm (unvec u));

  let tout = ref t1 in
  for iout = 1 to nout do
    let (t, flag) = Cvode.solve_normal cvode_mem !tout u
    in
    let nst = Cvode.get_num_steps cvode_mem in

    print_output t (vmax_norm (unvec u)) nst;
    tout := !tout +. dtout
  done;

  print_final_stats cvode_mem  (* Print some final statistics   *)

let n =
  match Sys.argv with
  | [|_; n|] -> int_of_string n
  | _ -> 1
let _ = for i = 1 to n do main () done


let _ = Gc.compact ()

