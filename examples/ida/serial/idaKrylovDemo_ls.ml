(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2009/09/30 23:25:59 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Aug 2014.
 * -----------------------------------------------------------------
 *
 * This example loops through the available iterative linear solvers:
 * SPGMR, SPBCG and SPTFQMR.
 *
 * Example problem for IDA: 2D heat equation, serial, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version loops through the Krylov solvers IDASpgmr, IDASpbcg
 * and IDASptfqmr.
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
 * The system is solved with IDA using the following Krylov
 * linear solvers: IDASPGMR, IDASPBCG and IDASPTFQMR. The
 * preconditioner uses the diagonal elements of the Jacobian only.
 * Routines for preconditioning, required by IDASP*, are supplied
 * here. The constraints u >= 0 are posed for all components. Output
 * is taken at t = 0, .01, .02, .04,..., 10.24.
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray

let nvscale = Nvector_serial.DataOps.n_vscale
and nvprod = Nvector_serial.DataOps.n_vprod
and nvmaxnorm = Nvector_serial.DataOps.n_vmaxnorm

(* Problem Constants *)
let nout  = 11
let mgrid = 10
let neq   = mgrid*mgrid

(* Linear Solver Loop Constants *)

type solver_to_use = USE_SPGMR | USE_SPBCG | USE_SPTFQMR

(* User data type *)

type user_data = {
  mm : int;                             (* number of grid points *)
  dx : float;
  coeff : float;
  pp : RealArray.t;                        (* vector of prec. diag. elements *)
}

(* Output functions *)
let printf = Printf.printf

(*
 * Print first lines of output (problem description)
 *)

let print_header rtol atol linsolver =
  printf "\nidaKrylovDemo_ls: Heat equation, serial example problem for IDA\n";
  printf "               Discretized heat equation on 2D unit square.\n";
  printf "               Zero boundary conditions,";
  printf " polynomial initial conditions.\n";
  printf "         Mesh dimensions: %d x %d" mgrid mgrid;
  printf "       Total system size: %d\n\n" neq;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Constraints set to force all solution components >= 0. \n";

  match linsolver with
  | USE_SPGMR ->
    printf "Linear solver: IDASPGMR, preconditioner using diagonal elements. \n"
  | USE_SPBCG ->
    printf "Linear solver: IDASPBCG, preconditioner using diagonal elements. \n"
  | USE_SPTFQMR ->
    printf "Linear solver: IDASPTFQMR, preconditioner using diagonal elements. \n"

(*
 * print_output: print max norm of solution and current solver statistics
 *)

let print_output mem t u linsolver =
  let umax = nvmaxnorm u

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


(*
 * resHeat: heat equation system residual function (user-supplied)      
 * This uses 5-point central differencing on the interior points, and   
 * includes algebraic equations for the boundary values.                
 * So for each interior point, the residual component has the form      
 *    res_i = u'_i - (central difference)_i                             
 * while for each boundary point, it is res_i = u_i.                     
 *)
let res_heat data t u (u' : RealArray.t) res =
  let coeff = data.coeff
  and mm    = data.mm in

  (* Initialize res to u, to take care of boundary equations. *)
  RealArray.blit_all u res;

  (* Loop over interior points; set res = up - (central difference). *)
  for j = 1 to mgrid-2 do
    let offset = mm*j in
    for i = 1 to mm-2 do
      let loc = offset + i in
      let dif1 = u.{loc-1}  +. u.{loc+1}  -. 2. *. u.{loc}
      and dif2 = u.{loc-mm} +. u.{loc+mm} -. 2. *. u.{loc} in
      res.{loc} <- u'.{loc} -. coeff *. ( dif1 +. dif2 )
    done
  done

(*
 * PsetupHeat: setup for diagonal preconditioner.   
 *                                                                 
 * The optional user-supplied functions PsetupHeat and          
 * PsolveHeat together must define the left preconditoner        
 * matrix P approximating the system Jacobian matrix               
 *                   J = dF/du + cj*dF/du'                         
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   
 * systems P z = r.   This is done in this case by keeping only    
 * the diagonal elements of the J matrix above, storing them as    
 * inverses in a vector pp, when computed in PsetupHeat, for    
 * subsequent use in PsolveHeat.                                 
 *                                                                 
 * In this instance, only cj and data (user data structure, with    
 * pp etc.) are used from the PsetupdHeat argument list.         
 *)
let p_setup_heat data jac =
  let c_j = jac.Ida.jac_coef
  and mm = data.mm
  and pp = data.pp in

  (* Initialize the entire vector to 1., then set the interior points to the
     correct value for preconditioning. *)
  RealArray.fill data.pp 1.;
    
  (* Compute the inverse of the preconditioner diagonal elements. *)
  let pelinv = 1./.(c_j +. 4.*.data.coeff) in

  for j = 1 to mm-2 do
    let offset = mm * j in
    for i = 1 to mm-2 do
      let loc = offset + i in
      pp.{loc} <- pelinv
    done
  done

(*
 * PsolveHeat: solve preconditioner linear system.              
 * This routine multiplies the input vector rvec by the vector pp 
 * containing the inverse diagonal Jacobian elements (previously  
 * computed in PrecondHeateq), returning the result in zvec.      
 *)
let p_solve_heat data jac r z delta =
  nvprod data.pp r z

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
  RealArray.fill u' 0.;

  (* res_heat sets res to negative of ODE RHS values at interior points. *)
  res_heat data 0. u u' res;

  (* Copy -res into up to get correct interior initial up values. *)
  nvscale (-1.) res u';

  (* Set up at boundary points to zero. *)
  for j = 0 to mm-1 do
    let offset = mm*j in
    for i = 0 to mm-1 do
      let loc = offset + i in
      if j == 0 || j == mm1 || i == 0 || i == mm1 then u'.{loc} <- 0.
    done
  done

let main() =
  (* Allocate N-vectors and the user data structure. *)
  let u = RealArray.create neq
  and u' = RealArray.create neq
  and res = RealArray.create neq
  and constraints = RealArray.create neq in
  let dx = 1. /. float_of_int (mgrid - 1) in
  let data =
    {
      mm = mgrid;
      dx = dx;
      coeff = 1. /. (dx *. dx);
      pp = RealArray.create neq;
    }
  in

  (* Initialize u, u'. *)

  set_initial_profile data u u' res;

  (* Set constraints to all nonnegative solution values. *)
  RealArray.fill constraints Ida.Constraint.non_negative;

  (* Assign various parameters. *)

  let t0   = 0.
  and t1   = 0.01
  and rtol = 0.
  and atol = 1.0e-3 in

  (* Wrap u and u' in nvectors.  Operations performed on the wrapped
     representation affect the originals u and u'.  *)
  let wu = Nvector_serial.wrap u
  and wu' = Nvector_serial.wrap u'
  in

  (* Call IDACreate with dummy linear solver *)

  let mem = Ida.init (Ida.Dls.dense None) (Ida.SStolerances (rtol, atol))
                     (res_heat data) ~t0:t0 wu wu' in
  Ida.set_constraints mem (Nvector_serial.wrap constraints);

  (* START: Loop through SPGMR, SPBCG and SPTFQMR linear solver modules *)
  let solvers = [|USE_SPGMR; USE_SPBCG; USE_SPTFQMR|] in
  for i = 0 to Array.length solvers - 1 do
    let linsolver = solvers.(i) in

    (* Note: the original C version of this example reinitializes the linear
       solver only when i > 0, but the OCaml interface prohibits setting the
       linear solver without a reinit (which isn't a sensible thing to do
       anyway).  So we reinit each time here.  *)

    if i <> 0 then
      (* Re-initialize uu, up. *)
      set_initial_profile data u u' res;

    (* Print header and reinit with a new solver module *)
    let spils_init = { Ida.Spils.prec_setup_fn = Some (p_setup_heat data);
                       Ida.Spils.prec_solve_fn = Some (p_solve_heat data);
                       Ida.Spils.jac_times_vec_fn = None;
                     }
    in
    begin
      match linsolver with
      | USE_SPGMR -> (printf " -------";
                      printf " \n| SPGMR |\n";
                      printf " -------\n";
                      flush stdout;
                      Ida.reinit mem ~linsolv:(Ida.Spils.spgmr None spils_init)
                        t0 wu wu')
      | USE_SPBCG -> (printf " -------";
                      printf " \n| SPBCG |\n";
                      printf " -------\n";
                      flush stdout;
                      Ida.reinit mem ~linsolv:(Ida.Spils.spbcg None spils_init)
                        t0 wu wu')
      | USE_SPTFQMR -> (printf " ---------";
                        printf " \n| SPTFQMR |\n";
                        printf " ---------\n";
                      flush stdout;
                        Ida.reinit mem ~linsolv:(Ida.Spils.sptfqmr None
                                                   spils_init)
                          t0 wu wu')
    end;

    (* Print output heading. *)
    print_header rtol atol linsolver;

    (* Print output table heading, and initial line of table. *)

    printf "\n   Output Summary (umax = max-norm of solution) \n\n";
    printf "  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n";
    printf "----------------------------------------------------------------------\n";

    (* Loop over output times, call IDASolve, and print results. *)

    let tout = ref t1 in
    for iout = 1 to nout do
      let (tret, flag) = Ida.solve_normal mem !tout wu wu' in
      print_output mem tret u linsolver;
      tout := !tout *. 2.
    done;

    (* Print remaining counters. *)
    let netf = Ida.get_num_err_test_fails mem
    and ncfn = Ida.get_num_nonlin_solv_conv_fails mem
    and ncfl = Ida.Spils.get_num_conv_fails mem in

    printf "\nError test failures            = %d\n" netf;
    printf "Nonlinear convergence failures = %d\n" ncfn;
    printf "Linear convergence failures    = %d\n" ncfl;

    if i < Array.length solvers - 1 then
      printf "\n======================================================================\n\n";

  done


(* Check environment variables for extra arguments.  *)
let reps =
  try int_of_string (Unix.getenv "NUM_REPS")
  with Not_found | Failure "int_of_string" -> 1
let gc_at_end =
  try int_of_string (Unix.getenv "GC_AT_END") <> 0
  with Not_found | Failure "int_of_string" -> false
let gc_each_rep =
  try int_of_string (Unix.getenv "GC_EACH_REP") <> 0
  with Not_found | Failure "int_of_string" -> false

(* Entry point *)
let _ =
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
