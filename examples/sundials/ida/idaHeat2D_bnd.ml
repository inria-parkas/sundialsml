(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2009/09/30 23:25:59 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, serial, banded.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the band solver IDABand, and IDACalcIC.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform M x M
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = M^2. Here M = 10.
 *
 * The system is solved with IDA using the banded linear system
 * solver, half-bandwidths equal to M, and default
 * difference-quotient Jacobian. For purposes of illustration,
 * IDACalcIC is called to compute correct values at the boundary,
 * given incorrect values as input initial guesses. The constraints
 * u >= 0 are posed for all components. Output is taken at
 * t = 0, .01, .02, .04, ..., 10.24. (Output at t = 0 is for
 * IDACalcIC cost statistics only.)
 * -----------------------------------------------------------------
 *)
module Ida = Ida_serial;;
module RealArray = Ida.RealArray
module Roots = Ida.Roots
module Dls = Ida.Dls
module Constraints = Ida.Constraints

let printf = Printf.printf
let vmax_norm = Nvector_array.Bigarray.array_nvec_ops.Nvector.Mutable.nvmaxnorm
let vscale = Nvector_array.Bigarray.array_nvec_ops.Nvector.Mutable.nvscale


(* Problem Constants *)
let nout = 11
let mgrid = 10                          (* mesh grid size *)
let neq = mgrid * mgrid
let bval = 0.1

type user_data = { mm : int; dx : float; coeff : float }

(*
 * heatres: heat equation system residual function                       
 * This uses 5-point central differencing on the interior points, and    
 * includes algebraic equations for the boundary values.                 
 * So for each interior point, the residual component has the form       
 *    res_i = u'_i - (central difference)_i                              
 * while for each boundary point, it is res_i = u_i.                     
 *)
let heatres t u u' resval data =
  let mm = data.mm
  and coeff = data.coeff
  in
  (* Initialize resval to u, to take care of boundary equations. *)
  RealArray.blit u resval;

  (* Loop over interior points; set res = u' - (central difference). *)
  for j = 1 to mm-2 do
    let offset = mm*j in
    for i = 1 to mm-2 do
      let loc = offset + i in
      resval.{loc} <-
        u'.{loc} -.
        coeff *.
        (u.{loc-1} +. u.{loc+1} +. u.{loc-mm} +. u.{loc+mm} -. 4. *. u.{loc})
    done
  done

let set_initial_profile data u u' id res =
  (* Initialize id to Differential. *)
  Ida.Id.fill id Ida.Id.Differential;

  let mm = data.mm in
  let mm1 = mm - 1 in

  (* Initialize u on all grid points. *)
  for j = 0 to mm-1 do
    let yfact = data.dx *. float_of_int j
    and offset = mm*j in
    for i = 0 to mm-1 do
      let xfact = data.dx *. float_of_int i
      and loc = offset + i in
      u.{loc} <- 16.0 *. xfact *. (1. -. xfact) *. yfact *. (1. -. yfact)
    done
  done;

  (* Initialize u' vector to 0. *)
  RealArray.fill u' 0.;

  (* heatres sets res to negative of ODE RHS values at interior points. *)
  heatres 0. u u' res data;

  (* Copy -res into u' to get correct interior initial u' values. *)
  vscale (-1.) res u';

  (* Finally, set values of u, u', and id at boundary points. *)
  for j = 0 to mm-1 do
    let offset = mm*j in
    for i = 0 to mm-1 do
      let loc = offset + i in
      if j = 0 || j = mm1 || i = 0 || i = mm1
      then (u.{loc} <- bval;
            u'.{loc} <- 0.;
            Ida.Id.set id loc Ida.Id.Algebraic)
    done
  done

let print_header rtol atol =
  printf "\nidaHeat2D_bnd: Heat equation, serial example problem for IDA\n";
  printf "          Discretized heat equation on 2D unit square.\n";
  printf "          Zero boundary conditions,";
  printf " polynomial initial conditions.\n";
  printf "          Mesh dimensions: %d x %d" mgrid mgrid;
  printf "        Total system size: %d\n\n" neq;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Constraints set to force all solution components >= 0. \n";
  printf "Linear solver: IDABAND, banded direct solver \n";
  printf "       difference quotient Jacobian, half-bandwidths = %d \n" mgrid;
  printf "IDACalcIC called with input boundary values = %g \n" bval;
  (* Print output table heading and initial line of table.  *)
  printf "\n   Output Summary (umax = max-norm of solution) \n\n";
  printf "  time       umax     k  nst  nni  nje   nre   nreLS    h      \n" ;
  printf " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n"

let print_output mem t u =
  let umax = vmax_norm u in
  let kused = Ida.get_last_order mem
  and nst = Ida.get_num_steps mem
  and nni = Ida.get_num_nonlin_solv_iters mem
  and nre = Ida.get_num_res_evals mem
  and hused = Ida.get_last_step mem
  and nje = Ida.Dls.get_num_jac_evals mem
  and nreLS = Ida.Dls.get_num_res_evals mem in
  
  printf " %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e \n"
         t umax kused nst nni nje nre nreLS hused

let main () =
  (* Create vectors uu, up, res, constraints, id.  *)
  let u = RealArray.make neq
  and u' = RealArray.make neq
  and res = RealArray.make neq
  and constraints = Constraints.create neq
  and id = Ida.Id.create neq in

  (* Create and load problem data block.  *)
  let data =
    let mm = mgrid in
    let dx = 1. /. (float_of_int mgrid -. 1.) in
    let coeff = 1. /. ( dx *. dx ) in
    { mm=mm; dx=dx; coeff=coeff }
  in

  (* Initialize u, u', id.  *)
  set_initial_profile data u u' id res;

  (* Set constraints to all NonNegative for nonnegative solution values.  *)
  Constraints.fill constraints Constraints.NonNegative;

  (* Set remaining input parameters.  *)
  let t0 = 0.
  and t1 = 0.01
  and rtol = 0.
  and atol = 1.0e-3
  in

  (* Call IDACreate and IDAMalloc to initialize solution, and call IDABand to
     specify the linear solver.  *)
  let mu = mgrid and ml = mgrid in
  let mem =
    Ida.init (Ida.Band ({ Ida.mupper=mu; Ida.mlower=ml }, None))
             (Ida.SStolerances (rtol, atol))
             (fun t u u' r -> heatres t u u' r data)
             ~t0:t0 u u'
  in
  Ida.set_constraints mem constraints;

  (* Call IDASetId and IDACalcIC to correct the initial values.  *)
  Ida.calc_ic_ya_yd' mem id t1;

  (* Print output heading. *)
  print_header rtol atol;

  print_output mem t0 u;

  (* Loop over output times, call IDASolve, and print results. *)
  let tout = ref t1 in
  for iout = 1 to nout do
    let (tret, flag) = Ida.solve_normal mem !tout u u' in
    print_output mem tret u;
    tout := 2. *. !tout
  done;

  (* Print remaining counters. *)
  let netf = Ida.get_num_err_test_fails mem
  and ncfn = Ida.get_num_nonlin_solv_conv_fails mem
  in
  printf "\n netf = %d,   ncfn = %d \n" netf ncfn
;;

let _ = main ()
let _ = Gc.compact ()
