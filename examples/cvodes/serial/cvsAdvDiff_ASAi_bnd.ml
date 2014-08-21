(*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2011/11/23 23:53:02 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2014.
 * -----------------------------------------------------------------
 * Adjoint sensitivity example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODES.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is the following:
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODE band linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Additionally, CVODES integrates backwards in time the
 * the semi-discrete form of the adjoint PDE:
 *   d(lambda)/dt = - d^2(lambda) / dx^2 + 0.5 d(lambda) / dx
 *                  - d^2(lambda) / dy^2 - 1.0
 * with homogeneous Dirichlet boundary conditions and final
 * conditions:
 *   lambda(x,y,t=t_final) = 0.0
 * whose solution at t = 0 represents the sensitivity of
 *   G = int_0^t_final int_x int _y u(t,x,y) dx dy dt
 * with respect to the initial conditions of the original problem.
 * -----------------------------------------------------------------
 *)

module Adjoint = Cvodes.Adjoint
module RealArray = Sundials.RealArray
module Col = Dls.BandMatrix.Col
module Dls = Cvode.Dls
let unvec = Sundials.unvec

let printf = Printf.printf

let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

(* Header files with a description of contents used in cvbanx.c *)


(* Problem Constants *)

let xmax   = 2.0        (* domain boundaries         *)
let ymax   = 1.0
let mx     = 40         (* mesh dimensions           *)
let my     = 20
let neq    = mx * my    (* number of equations       *)
let atol   = 1.0e-5     (* scalar absolute tolerance *)
let rtolb  = 1.0e-6
let t0     = 0.0        (* initial time              *)
let t1     = 0.1        (* first output time         *)
let dtout  = 0.1        (* output time increment     *)
let nout   = 10         (* number of output times    *)

let tout   = 1.0        (* final time                    *)
let nstep  = 50         (* check point saved every NSTEP *)

let zero  = 0.0
let one   = 1.0
let two   = 2.0
let five  = 5.0

(* User-defined vector access macro IJth *)

(* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. 
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The vdata array is obtained via the macro call vdata = N_VDATA_S(v),
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

(* f routine. right-hand side of forward ODE. *)

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

(* Jac function. Jacobian of forward ODE. *)

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

(* fB function. Right-hand side of backward ODE. *)

let fB data tB u uB uBdot =
  (* Extract needed constants from data *)
  let hordc = data.hdcoef in
  let horac = data.hacoef in
  let verdc = data.vdcoef in

  (* Loop over all grid points. *)
  for j=1 to my do
    for i=1 to mx do
      (* Extract u at x_i, y_j and four neighboring points *)
      let uBij = ijth uB i j in
      let uBdn = if j = 1  then zero else ijth uB i (j-1) in
      let uBup = if j = my then zero else ijth uB i (j+1) in
      let uBlt = if i = 1  then zero else ijth uB (i-1) j in
      let uBrt = if i = mx then zero else ijth uB (i+1) j in

      (* Set diffusion and advection terms and load into udot *)

      let hdiffB = hordc*.(-. uBlt +. two*.uBij -. uBrt) in
      let hadvB  = horac*.(uBrt -. uBlt) in
      let vdiffB = verdc*.(-. uBup +. two*.uBij -. uBdn) in
      set_ijth uBdot i j (hdiffB +. hadvB +. vdiffB -. one)
    done
  done

(* JacB function. Jacobian of backward ODE. *)

let jacb data { Adjoint.mupper = muB; Adjoint.mlower = mlB }
              { Adjoint.jac_t = tB;
                Adjoint.jac_y = u;
                Adjoint.jac_yb = uB;
                Adjoint.jac_fyb = fuB;
                Adjoint.jac_tmp = (tmp1B, tmp2B, tmp3B) } jb =

  (* The Jacobian of the adjoint system is: JB = -J^T *)
  let hordc = data.hdcoef in
  let horac = data.hacoef in
  let verdc = data.vdcoef in

  for j=1 to my do
    for i=1 to mx do
      let k = j-1 + (i-1)*my in
      let kthCol = Col.get_col jb k in

      (* set the kth column of J *)
      Col.set kthCol k k (two*.(verdc+.hordc));
      if i != 1  then Col.set kthCol (k-my) k (-. hordc +. horac);
      if i != mx then Col.set kthCol (k+my) k (-. hordc -. horac);
      if j != 1  then Col.set kthCol (k-1) k  (-. verdc);
      if j != my then Col.set kthCol (k+1) k  (-. verdc)
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

(* Print results after backward integration *)

let print_output uB data =
  let x = ref 0.0 in
  let y = ref 0.0 in

  let dx = data.dx in
  let dy = data.dy in

  let uBmax = ref 0.0 in

  for j=1 to my do
    for i=1 to mx do
      let uBij = ijth uB i j in
      if (abs_float uBij > !uBmax) then begin
        uBmax := uBij;
        x := float i *.dx;
        y := float j *.dy
      end
    done
  done;

  printf "\nMaximum sensitivity\n";
  printf "  lambda max = %e\n" !uBmax;
  printf "at\n";
  printf "  x = %e\n  y = %e\n" !x !y

let main () =
  (* Create a serial vector *)
  let u = RealArray.create neq in  (* Allocate u vector *)
  let u_nvec = Nvector_serial.wrap u in

  let reltol = zero  (* Set the tolerances *)
  and abstol = atol
  in

  let dx = xmax /. float(mx + 1); (* Set grid coefficients in data *)
  and dy = ymax /. float(my + 1);
  in
  let data = {
    dx = dx; dy = dy;
    hdcoef = one /. (dx *. dx);
    hacoef = 1.5 /. (two *. dx);
    vdcoef = one /. (dy *. dy);
  } in

  set_ic u data;  (* Initialize u vector *)

  (* Create and allocate CVODES memory for forward run *)
  printf "\nCreate and allocate CVODES memory for forward runs\n";

  (* Call CVBand with  bandwidths ml = mu = MY, *)
  let solver = Cvode.Dls.band {Cvode.mupper = my; Cvode.mlower = my}
                              (Some (jac data))
  in
  let cvode_mem = Cvode.init Cvode.BDF (Cvode.Newton solver)
                             (Cvode.SStolerances (reltol, abstol))
                             ~t0:t0 (f data) u_nvec
  in
  Gc.compact ();

  (* Allocate global memory *)
  printf "\nAllocate global memory\n";

  Adjoint.init cvode_mem nstep Adjoint.IHermite;

  (* Perform forward run *)
  printf "\nForward integration\n";
  let t, ncheck, _ = Adjoint.forward_normal cvode_mem tout u_nvec in
  printf "\nncheck = %d\n" ncheck;

  (* Allocate uB *)
  let uB = Nvector_serial.make neq 0.0 in

  (* Create and allocate CVODES memory for backward run *)
  printf "\nCreate and allocate CVODES memory for backward run\n";

  let bsolver = Adjoint.Dls.band {Adjoint.mupper = my; Adjoint.mlower = my}
                                 (Some (jacb data)) in
  let bcvode_mem = Adjoint.init_backward cvode_mem
        Cvode.BDF
        (Adjoint.Newton bsolver)
        (Adjoint.SStolerances (rtolb, atol))
        (Adjoint.Basic (fB data)) tout uB in

  (* Perform backward integration *)
  printf "\nBackward integration\n";
  Adjoint.backward_normal cvode_mem t0;
  let _ = Adjoint.get bcvode_mem uB in

  print_output (unvec uB) data


let n =
  match Sys.argv with
  | [|_; n|] -> int_of_string n
  | _ -> 1
let _ = for i = 1 to n do main (); Gc.compact () done
let _ = Gc.full_major ()

