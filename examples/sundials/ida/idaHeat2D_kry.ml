(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2009/09/30 23:25:59 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, serial, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver IDASpgmr.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
 * PDE is treated with central differences on a uniform M x M grid.
 * The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = M^2. Here M = 10.
 *
 * The system is solved with IDA using the Krylov linear solver
 * IDASPGMR. The preconditioner uses the diagonal elements of the
 * Jacobian only. Routines for preconditioning, required by
 * IDASPGMR, are supplied here. The constraints u >= 0 are posed
 * for all components. Output is taken at t = 0, .01, .02, .04,
 * ..., 10.24. Two cases are run -- with the Gram-Schmidt type
 * being Modified in the first case, and Classical in the second.
 * The second run uses IDAReInit and IDAReInitSpgmr.
 * -----------------------------------------------------------------
 *)
module Ida = Ida_serial
module Carray = Ida.Carray

(* Problem Constants *)
let nout  = 11
let mgrid = 10
let neq   = mgrid*mgrid

(* Shorthands *)
let nvscale = Nvector_array.Bigarray.array_nvec_ops.Nvector.Mutable.nvscale
let nvprod = Nvector_array.Bigarray.array_nvec_ops.Nvector.Mutable.nvprod
let nvmaxnorm = Nvector_array.Bigarray.array_nvec_ops.Nvector.Mutable.nvmaxnorm
let printf = Printf.printf

(* User data *)
type user_data =
  {  
    mm    : int;                        (* number of grid points *)
    dx    : float;
    coeff : float;
    pp    : Ida.nvec;                   (* vector of prec. diag. elements *)
  }

(*
 * resHeat: heat equation system residual function (user-supplied)      
 * This uses 5-point central differencing on the interior points, and   
 * includes algebraic equations for the boundary values.                
 * So for each interior point, the residual component has the form      
 *    res_i = u'_i - (central difference)_i                             
 * while for each boundary point, it is res_i = u_i.                     
 *)
let res_heat data t u u' r =
  let coeff = data.coeff
  and mm    = data.mm in
  
  (* Initialize r to u, to take care of boundary equations. *)
  Carray.blit u r;
  
  (* Loop over interior points; set res = up - (central difference).  *)
  for j = 1 to mgrid-2 do
    let offset = mm*j in
    for i = 1 to mm-2 do
      let loc = offset + i in
      let dif1 = u.{loc-1} +. u.{loc+1} -. 2. *. u.{loc}
      and dif2 = u.{loc-mm} +. u.{loc+mm} -. 2. *. u.{loc} in
      r.{loc} <- u'.{loc} -. coeff *. (dif1 +. dif2)
    done
  done

(*
 * p_setup_heat: setup for diagonal preconditioner for idaHeat2D_kry.   
 *                                                                 
 * The optional user-supplied functions p_setup_heat and          
 * p_solve_heat together must define the left preconditoner        
 * matrix P approximating the system Jacobian matrix               
 *                   J = dF/du + cj*dF/du'                         
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   
 * systems P z = r.   This is done in this case by keeping only    
 * the diagonal elements of the J matrix above, storing them as    
 * inverses in a vector pp, when computed in p_setup_heat, for    
 * subsequent use in p_solve_heat.                                 
 *                                                                 
 * In this instance, only cj and data (user data structure, with    
 * pp etc.) are used from the p_setup_heat argument list.         
 *)
let p_setup_heat data jac =
  let pp = data.pp
  and mm = data.mm
  and c_j = jac.Ida.jac_coef
  in

  (* Initialize the entire vector to 1., then set the interior points to the
     correct value for preconditioning. *)
  Carray.fill pp 1.;
  
  (* Compute the inverse of the preconditioner diagonal elements. *)
  let pelinv = 1. /. (c_j +. 4.*.data.coeff) in

  for j = 1 to mm-2 do
    let offset = mm*j in
    for i = 1 to mm-2 do
      let loc = offset + i in
      pp.{loc} <- pelinv
    done
  done

(*
 * p_solve_heat: solve preconditioner linear system.              
 * This routine multiplies the input vector rvec by the vector pp 
 * containing the inverse diagonal Jacobian elements (previously  
 * computed in PrecondHeateq), returning the result in zvec.      
 *)
let p_solve_heat data jac rvec zvec delta =
  nvprod data.pp rvec zvec

(*
 * set_initial_profile: routine to initialize u and u' vectors.
 *)

let set_initial_profile data u u' res =
  let mm = data.mm in

  (* Initialize uu on all grid points. *)
  let mm1 = mm - 1 in
  for j = 0 to mm-1 do
    let yfact = data.dx *. float_of_int j
    and offset = mm*j in
    for i = 0 to mm-1 do
      let xfact = data.dx *. float_of_int i
      and loc = offset + i in
      u.{loc} <- 16.0 *. xfact *. (1. -. xfact) *. yfact *. (1. -. yfact)
    done
  done;

  (* Initialize up vector to 0. *)
  Carray.fill u' 0.;

  (* res_heat sets res to negative of ODE RHS values at interior points. *)
  res_heat data 0. u u' res;

  (* Copy -res into up to get correct interior initial up values. *)
  nvscale (-1.) res u';

  (* Set up at boundary points to zero. *)
  for j = 0 to mm-1 do
    let offset = mm*j in
    for i = 0 to mm-1 do
      let loc = offset + i in
      if j = 0 || j = mm1 || i = 0 || i = mm1 then
        u'.{loc} <- 0.
    done
  done

(*
 * Print first lines of output (problem description)
 *)
let print_header rtol atol =
  printf "\nidaHeat2D_kry: Heat equation, serial example problem for IDA \n";
  printf "         Discretized heat equation on 2D unit square. \n";
  printf "         Zero boundary conditions,";
  printf " polynomial initial conditions.\n";
  printf "         Mesh dimensions: %d x %d" mgrid mgrid;
  printf "        Total system size: %d\n\n" neq;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Constraints set to force all solution components >= 0. \n";
  printf "Linear solver: IDASPGMR, preconditioner using diagonal elements. \n"

(*
 * print_output: print max norm of solution and current solver statistics
 *)
let print_output mem t u =
  let umax = nvmaxnorm u;
  and kused = Ida.get_last_order mem
  and nst = Ida.get_num_steps mem
  and nni = Ida.get_num_nonlin_solv_iters mem
  and nre = Ida.get_num_res_evals mem
  and hused = Ida.get_last_step mem
  and nje = Ida.Spils.get_num_jtimes_evals mem
  and nreLS = Ida.Spils.get_num_res_evals mem
  and npe = Ida.Spils.get_num_prec_evals mem
  and nps = Ida.Spils.get_num_prec_solves mem in
  printf " %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e  %3d %3d\n"
         t umax kused nst nni nje nre nreLS hused npe nps


let main () =
  (* Allocate N-vectors and the user data structure. *)
  let u = Carray.create neq
  and u' = Carray.create neq
  and res = Carray.create neq
  and constraints = Ida.Constraints.create neq in

  let dx = 1. /. float_of_int (mgrid - 1) in
  let data = { pp = Carray.create neq;
               mm = mgrid;
               dx = dx;
               coeff = 1. /. (dx *. dx);
             }
  in

  (* Initialize u, u'. *)
  set_initial_profile data u u' res;

  (* Set constraints to all NonNegative. *)
  Ida.Constraints.fill constraints Ida.Constraints.NonNegative;

  (* Assign various parameters. *)

  let t0   = 0.0
  and t1   = 0.01
  and rtol = 0.0
  and atol = 1.0e-3 in

  (* Call IDACreate to initialize solution with SPGMR linear solver.  *)

  let solver = Ida.Spgmr { Ida.maxl = Some 5;
                           Ida.prec_setup_fn = Some (p_setup_heat data);
                           Ida.prec_solve_fn = Some (p_solve_heat data);
                           Ida.jac_times_vec_fn = None;
                         }
  in
  let mem = Ida.init solver (Ida.SStolerances (rtol, atol))
                     (res_heat data) ~t0:t0 u u' in
  Ida.set_constraints mem constraints;

  (* Print output heading. *)
  print_header rtol atol;
  
  (* 
   * -------------------------------------------------------------------------
   * CASE I 
   * -------------------------------------------------------------------------
   *)
  
  (* Print case number, output table heading, and initial line of table. *)

  printf "\n\nCase 1: gsytpe = MODIFIED_GS\n";
  printf "\n   Output Summary (umax = max-norm of solution) \n\n";
  printf "  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n" ;
  printf "----------------------------------------------------------------------\n";

  (* Loop over output times, call IDASolve, and print results. *)

  let tout = ref t1 in
  for iout = 1 to nout do
    let (tret, flag) = Ida.solve_normal mem !tout u u' in
    print_output mem tret u;
    tout := !tout *. 2.
  done;

  (* Print remaining counters. *)

  let netf = Ida.get_num_err_test_fails mem
  and ncfn = Ida.get_num_nonlin_solv_conv_fails mem
  and ncfl = Ida.Spils.get_num_conv_fails mem in
  printf "\nError test failures            = %d\n" netf;
  printf "Nonlinear convergence failures = %d\n" ncfn;
  printf "Linear convergence failures    = %d\n" ncfl;

  (*
   * -------------------------------------------------------------------------
   * CASE II
   * -------------------------------------------------------------------------
   *)
  
  (* Re-initialize u, u'. *)

  set_initial_profile data u u' res;
  
  (* Re-initialize IDA and IDASPGMR *)

  Ida.reinit mem t0 u u';
  Ida.Spils.set_gs_type mem Spils.ClassicalGS;
  
  (* Print case number, output table heading, and initial line of table. *)
  printf "\n\nCase 2: gstype = CLASSICAL_GS\n";
  printf "\n   Output Summary (umax = max-norm of solution) \n\n";
  printf "  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n" ;
  printf "----------------------------------------------------------------------\n";

  (* Loop over output times, call IDASolve, and print results. *)
  let tout = ref t1 in
  for iout = 1 to nout do
    let (tret, flag) = Ida.solve_normal mem !tout u u' in
    print_output mem tret u;
    tout := !tout *. 2.
  done;

  (* Print remaining counters. *)

  let netf = Ida.get_num_err_test_fails mem
  and ncfn = Ida.get_num_nonlin_solv_conv_fails mem
  and ncfl = Ida.Spils.get_num_conv_fails mem in
  printf "\nError test failures            = %d\n" netf;
  printf "Nonlinear convergence failures = %d\n" ncfn;
  printf "Linear convergence failures    = %d\n" ncfl

let _ = main ()
let _ = Gc.compact ()
