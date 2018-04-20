(*
 * -----------------------------------------------------------------
 * $Revision:  $
 * $Date:  $
 * -----------------------------------------------------------------
 * Programmer(s): Ting Yan @ SMU
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Dec 2016.
 * -----------------------------------------------------------------
 * Example program for IDAS: Food web problem, OpenMP, GMRES,
 * user-supplied preconditioner
 *
 * This example program uses the IDASPGMR as the linear
 * solver, and IDACalcIC for initial condition calculation.
 *
 * The mathematical problem solved in this example is a DAE system
 * that arises from a system of partial differential equations after
 * spatial discretization. The PDE system is a food web population
 * model, with predator-prey interaction and diffusion on the unit
 * square in two dimensions. The dependent variable vector is:
 *
 *         1   2         ns
 *   c = (c , c ,  ..., c  ) , ns = 2 * np
 *
 * and the PDE's are as follows:
 *
 *     i             i      i
 *   dc /dt = d(i)*(c    + c  )  +  R (x,y,c)   (i = 1,...,np)
 *                   xx     yy       i
 *
 *              i      i
 *   0 = d(i)*(c    + c  )  +  R (x,y,c)   (i = np+1,...,ns)
 *              xx     yy       i
 *
 *   where the reaction terms R are:
 *
 *                   i             ns         j
 *   R  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *  a(i,i) = -AA   (all i)
 *  a(i,j) = -GG   (i <= np , j >  np)
 *  a(i,j) =  EE   (i >  np, j <= np)
 *  all other a(i,j) = 0
 *  b(i) = BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y)) (i <= np)
 *  b(i) =-BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y)) (i  > np)
 *  d(i) = DPREY   (i <= np)
 *  d(i) = DPRED   (i > np)
 *
 * The various scalar parameters required are set using '#define'
 * statements or directly in routine InitUserData. In this program,
 * np = 1, ns = 2. The boundary conditions are homogeneous Neumann:
 * normal derivative = 0.
 *
 * A polynomial in x and y is used to set the initial values of the
 * first np variables (the prey variables) at each x,y location,
 * while initial values for the remaining (predator) variables are
 * set to a flat value, which is corrected by IDACalcIC.
 *
 * The PDEs are discretized by central differencing on a MX by MY
 * mesh.
 *
 * The DAE system is solved by IDAS using the IDABAND linear solver.
 * Output is printed at t = 0, .001, .01, .1, .4, .7, 1.
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
 *      % ./idasFoodWeb_kry_omp
 * If the user want to specify the number of threads to use
 *      % ./idasFoodWeb_kry_omp num_threads
 * where num_threads is the number of threads the user want to use
 *
 * -----------------------------------------------------------------
 * References:
 * [1] Peter N. Brown and Alan C. Hindmarsh,
 *     Reduced Storage Matrix Methods in Stiff ODE systems, Journal
 *     of Applied Mathematics and Computation, Vol. 31 (May 1989),
 *     pp. 40-91.
 *
 * [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Using Krylov Methods in the Solution of Large-Scale
 *     Differential-Algebraic Systems, SIAM J. Sci. Comput., 15
 *     (1994), pp. 1467-1488.
 *
 * [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Consistent Initial Condition Calculation for Differential-
 *     Algebraic Systems, SIAM J. Sci. Comput., 19 (1998),
 *     pp. 1495-1512.
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray
module RealArray2 = Sundials.RealArray2
module Roots = Sundials.Roots
module LintArray = Sundials.LintArray

let printf = Printf.printf

(* Problem Constants. *)
let nprey       = 1                     (* No. of prey (= no. of predators). *)
let num_species = 2*nprey
let pi          = 3.1415926535898
let fourpi      = 4.0*.pi
let mx          = 20                    (* MX = number of x mesh points *)
let my          = 20                    (* MY = number of y mesh points *)
let nsmx        = num_species * mx
let neq         = num_species * mx * my
let aa          = 1.0                   (* Coefficient in above eqns. for a *)
let ee          = 10000.                (* Coefficient in above eqns. for a *)
let gg          = 0.5e-6                (* Coefficient in above eqns. for a *)
let bb          = 1.0                   (* Coefficient in above eqns. for b *)
let dprey       = 1.0                   (* Coefficient in above eqns. for d *)
let dpred       = 0.05                  (* Coefficient in above eqns. for d *)
let alpha       = 50.                   (* Coefficient alpha in above eqns. *)
let beta        = 1000.                 (* Coefficient beta in above eqns. *)
let ax          = 1.0                   (* Total range of x variable *)
let ay          = 1.0                   (* Total range of y variable *)
let rtol        = 1.e-5                 (* Relative tolerance *)
let atol        = 1.e-5                 (* Absolute tolerance *)
let nout        = 6                     (* Number of output times *)
let tmult       = 10.0                  (* Multiplier for tout values *)
let tadd        = 0.3                   (* Increment for tout values *)

(* User-defined vector and accessor functions: index and ij_v.  These are
 * defined in order to express the underlying 3-D structure of the dependent
 * variable vector from its underlying 1-D storage (an N_Vector).  index i j k
 * returns the index corresponding to the location in vv corresponding to
 * species index is = k, x-index ix = i, and y-index jy = j.  ij_v returns a
 * slice of vv corresponding to x-index ix = i, and y-index jy = j, and species
 * index anywhere between 0 and num_species.  *)
let index i j k = i*num_species + j*nsmx + k
let ij_v (v : RealArray.t) i j =
  (* v is type annotated so that the layout of the bigarray is fixed.  This
     ensures that the compiler can inline accesses to the bigarray elements. *)
  let i0 = index i j 0
  and iend = index i j num_species in
  Bigarray.Array1.sub v i0 (iend - i0)

type user_data =
  { neq   : int;
    ns    : int;
    np    : int;
    mx    : int;
    my    : int;
    dx    : float;
    dy    : float;
    acoef : Matrix.ArrayDense.t;
    cox   : float array;
    coy   : float array;
    bcoef : float array;
    rates : RealArray.t;
    pp    : Matrix.ArrayDense.t array array;
    pivot : LintArray.t array array;
    ewt   : Nvector_openmp.t;
    mutable ida_mem :
      (Nvector_openmp.data, Nvector_openmp.kind) Ida.session option;
  }

(* init_user_data: Load problem constants in webdata (of type user_data).  *)
let init_user_data num_threads =
  let webdata =
    { neq   = neq;
      mx    = mx;
      my    = my;
      ns    = num_species;
      np    = nprey;
      dx    = ax /. float_of_int (mx-1);
      dy    = ay /. float_of_int (my-1);
      acoef = Matrix.ArrayDense.create num_species num_species;
      cox   = Array.make num_species 0.;
      coy   = Array.make num_species 0.;
      bcoef = Array.make num_species 0.;
      rates = RealArray.create neq;
      pp    = Array.init mx (fun jx ->
                Array.init my (fun jy ->
                  Matrix.ArrayDense.create num_species num_species));
      pivot = Array.init mx (fun jx ->
                Array.init my (fun jy ->
                  LintArray.make num_species 0));
      ewt   = Nvector_openmp.make num_threads neq 0.0;
      ida_mem = None;
    } in
  (* Set up shorthands. *)
  let acoef = webdata.acoef
  and bcoef = webdata.bcoef
  and cox = webdata.cox
  and coy = webdata.coy
  and dx = webdata.dx
  and dy = webdata.dy
  and np = webdata.np
  in
  (* Set up the coefficients a and b, and others found in the equations. *)
  let np = np
  and dx2 = dx *. dx
  and dy2 = dy *. dy in

  let open Matrix.ArrayDense in
  for i = 0 to np-1 do
    (* Fill in the portion of acoef in the four quadrants, row by row. *)
    for j = 0 to np-1 do
      set acoef (np+j) i      (-. gg);
      set acoef j      (i+np)      ee;
      set acoef j      i           0.;
      set acoef (np+j) (i+np)      0.;
    done;

    (* Reset the diagonal elements of acoef to -AA.  *)
    set acoef i i (-. aa);
    set acoef (i+np) (i+np) (-. aa);

    (* Set coefficients for b and diffusion terms.  *)
    bcoef.(i) <- bb; bcoef.(i+np) <- -. bb;
    cox.(i) <- dprey /. dx2; cox.(i+np) <- dpred /. dx2;
    coy.(i) <- dprey /. dy2; coy.(i+np) <- dpred /. dy2
  done;
  webdata

(* web_rates: Evaluate reaction rates at a given spatial point.
 * At a given (x,y), evaluate the array of ns reaction terms R. *)
let web_rates webdata x y ((cxy : RealArray.t), cxy_off)
                          ((ratesxy : RealArray.t), ratesxy_off) =
  let acoef = RealArray2.unwrap webdata.acoef
  and bcoef = webdata.bcoef in

  for is = 0 to num_species-1 do
    (* ratesxy.{is} <- dotprod cxy (Directdensematrix.column acoef is) *)
    ratesxy.{ratesxy_off + is} <- 0.;
    for j = 0 to num_species-1 do
      ratesxy.{ratesxy_off + is} <- ratesxy.{ratesxy_off + is}
                                  +. cxy.{cxy_off + j} *. acoef.{is, j}
    done
  done;

  let fac = 1. +. alpha*.x*.y +. beta*.sin(fourpi*.x)*.sin(fourpi*.y) in

  for is = 0 to num_species-1 do
    ratesxy.{ratesxy_off + is} <- cxy.{cxy_off + is}*.( bcoef.(is)*.fac
                                    +. ratesxy.{ratesxy_off + is} )
  done

(* fweb: Rate function for the food-web problem.
 * This routine computes the right-hand sides of the system equations,
 * consisting of the diffusion term and interaction term.
 * The interaction term is computed by the function WebRates.  *)
let fweb webdata t c (crate : RealArray.t) =
  let cox = webdata.cox
  and coy = webdata.coy in
  (* Loop over grid points, evaluate interaction vector (length ns), form
     diffusion difference terms, and load crate. *)
  for jy = 0 to my-1 do
    let yy = webdata.dy *. float_of_int jy
    and idyu = if jy <> my-1 then nsmx else -nsmx
    and idyl = if jy <> 0    then nsmx else -nsmx
    in
    for jx = 0 to mx-1 do
      let xx = webdata.dx *. float_of_int jx
      and idxu = if jx <> mx-1 then num_species else -num_species
      and idxl = if jx <> 0    then num_species else -num_species in

      (* Get interaction vector at this grid point. *)
      web_rates webdata xx yy (c, index jx jy 0)
                              (webdata.rates, index jx jy 0);

      (* Loop over species, do differencing, load crate segment. *)
      let off = index jx jy 0 in
      for is = 0 to num_species-1 do
        (* Differencing in y. *)
        let dcyli = c.{off + is} -. c.{off + (is - idyl)}
        and dcyui = c.{off + idyu + is} -. c.{off + is} in
        (* Differencing in x. *)
        let dcxli = c.{off + is} -. c.{off + (is - idxl)}
        and dcxui = c.{off + (is + idxu)} -. c.{off + is} in
        (* Compute the crate values at (xx,yy). *)
        crate.{off + is} <- (coy.(is) *. (dcyui -. dcyli)
                            +. cox.(is) *. (dcxui -. dcxli)
                            +. webdata.rates.{index jx jy is})
      done
    done
  done

(* System residual function for predator-prey system.  This routine calls fweb
 * to get all the right-hand sides of the equations, then loads the residual
 * vector accordingly, using c' in the case of prey species.  *)
let resweb webdata t c (c' : RealArray.t) res =
  let np = webdata.np in

  (* Call Fweb to set res to vector of right-hand sides. *)
  fweb webdata t c res;

  (* Loop over all grid points, setting residual values appropriately
     for differential or algebraic components. *)
  for jy = 0 to my-1 do
    let yloc = nsmx * jy in
    for jx = 0 to mx-1 do
      let loc = yloc + num_species * jx in
      for is = 0 to num_species-1 do
        if is < np
        then res.{loc+is} <- c'.{loc+is} -. res.{loc+is}
        else res.{loc+is} <- -. res.{loc+is}
      done
    done
  done

let precond webdata jac =
  let cc = jac.Ida.jac_y
  and cp = jac.Ida.jac_y'
  and cj = jac.Ida.jac_coef
  in

  let mx =    webdata.mx in
  let my =    webdata.my in
  let np =    webdata.np in
  let dx =    webdata.dx in
  let dy =    webdata.dy in
  let rates = webdata.rates in
  let ns =    webdata.ns in

  let perturb_rates = RealArray.create num_species in

  let uround = Sundials.unit_roundoff in
  let sqru = sqrt uround in
  let mem =
    match webdata.ida_mem with
    | Some m -> m
    | None -> invalid_arg "Internal error: webdata.ida_mem not set"
  in
  let ewt = Nvector_openmp.unwrap webdata.ewt in
  Ida.get_err_weights mem webdata.ewt;
  let hh = Ida.get_current_step mem in

  for jy = 0 to my-1 do
    let yy = float_of_int jy *. dy in

    for jx = 0 to mx-1 do
      let xx = float_of_int jx *. dx in
      let pxy = webdata.pp.(jx).(jy) in
      let off = index jx jy 0 in

      for js = 0 to ns-1 do
        let inc = sqru*.(max (abs_float (cc.{off+js}))
                             (max (hh*.abs_float (cp.{off+js}))
                                  (1.0/.ewt.{off+js}))) in
        let cctemp = cc.{off+js} in(* Save the (js,jx,jy) element of cc. *)
        cc.{off+js} <- cc.{off+js} +. inc;    (* Perturb the (js,jx,jy) element of cc. *)
        let fac = -.1.0/.inc in

        web_rates webdata xx yy (cc, off) (perturb_rates, 0);

        let pxycol = RealArray2.col pxy js in

        for is = 0 to ns-1 do
          pxycol.{is} <- (perturb_rates.{is} -. rates.{off+is})*.fac
        done;

        (* Add partial with respect to cp. *)
        if js < np then pxycol.{js} <- pxycol.{js} +. cj;

        cc.{off+js} <- cctemp; (* Restore (js,jx,jy) element of cc. *)

      done; (* End of js loop. *)

      (* Do LU decomposition of matrix block for grid point (jx,jy). *)
      (try Matrix.ArrayDense.getrf pxy webdata.pivot.(jx).(jy)
       with Matrix.ZeroDiagonalElement _ -> raise Sundials.RecoverableFailure)

    done (* End of jx loop. *)
  done (* End of jy loop. *)

let psolve webdata jac rvec zvec delta =
  RealArray.blit rvec zvec;

  (* Loop through subgrid and apply preconditioner factors at each point. *)
  for jx = 0 to webdata.mx-1 do
    for jy = 0 to webdata.my-1 do

      (* For grid point (jx,jy), do backsolve on local vector.
         zxy is the address of the local portion of zvec, and
         Pxy is the address of the corresponding block of PP.  *)
      let zxy = ij_v zvec jx jy in
      let pxy = webdata.pp.(jx).(jy) in
      let pivot = webdata.pivot.(jx).(jy) in
      Matrix.ArrayDense.getrs pxy pivot zxy
    done (* End of jy loop. *)
  done (* End of jx loop. *)

let set_initial_profiles webdata c c' id =
  (* Loop over grid, load c values and id values.  *)
  for jy = 0 to my-1 do
    let y = webdata.dy *. float_of_int jy
    and yloc = nsmx * jy in
    for jx = 0 to mx-1 do
      let x = webdata.dx *. float_of_int jx in
      let xyfactor = 16.0*.x*.(1.-.x)*.y*.(1.-.y) in
      let xyfactor = xyfactor *. xyfactor
      and loc = yloc + num_species*jx in
      for is = 0 to num_species-1 do
        if is < webdata.np
        then (c.{loc+is} <- 10.0 +. float_of_int (is+1) *. xyfactor;
              id.{loc+is} <- Ida.VarId.differential)
        else (c.{loc+is} <- 1.0e5;
              id.{loc+is} <- Ida.VarId.algebraic)
      done
    done
  done;

  (* Set c' for the prey by calling the function fweb. *)
  fweb webdata 0. c c';

  (* Set c' for predators to 0. *)
  for jy = 0 to my-1 do
    let yloc = nsmx * jy in
    for jx = 0 to mx-1 do
      let loc = yloc + num_species * jx in
      for is = webdata.np to num_species-1 do
        c'.{loc+is} <- 0.0
      done
    done
  done

(* Print first lines of output (problem description) *)
let print_header maxl rtol atol =
  printf "\nidasFoodWeb_kry_omp: Predator-prey DAE OpenMP example problem for IDAS \n\n";
  printf "Number of species ns: %d" num_species;
  printf "     Mesh dimensions: %d x %d" mx my;
  printf "     System size: %d\n" neq;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Linear solver: IDASpgmr,  Spgmr parameters maxl = %d\n" maxl;
  printf "CalcIC called to correct initial predator concentrations.\n\n";
  printf "-----------------------------------------------------------\n";
  printf "  t        bottom-left  top-right";
  printf "    | nst  k      h\n";
  printf "-----------------------------------------------------------\n\n";
;;
(*
 * PrintOutput: Print output values at output time t = tt.
 * Selected run statistics are printed.  Then values of the concentrations
 * are printed for the bottom left and top right grid points only.
 *)

let print_output mem c t =
  let kused = Ida.get_last_order mem
  and nst   = Ida.get_num_steps mem
  and hused = Ida.get_last_step mem
  and c_bl = ij_v c 0 0
  and c_tr = ij_v c (mx-1) (my-1) in

  printf "%8.2e %12.4e %12.4e   | %3d  %1d %12.4e\n"
         t c_bl.{0} c_tr.{0} nst kused hused;
  for i = 1 to num_species-1 do
    printf "         %12.4e %12.4e   |\n" c_bl.{i} c_tr.{i}
  done;
  printf "\n"

(*
 * PrintFinalStats: Print final run data contained in iopt.
 *)

let print_final_stats mem =
  let open Ida in
  let nst     = get_num_steps mem
  and sli     = Spils.get_num_lin_iters mem
  and nre     = get_num_res_evals mem
  and netf    = get_num_err_test_fails mem
  and nps     = Spils.get_num_prec_solves mem
  and npevals = Spils.get_num_prec_evals mem
  in
  printf "-----------------------------------------------------------\n";
  printf "Final run statistics: \n\n";
  printf "Number of steps                       = %d\n" nst;
  printf "Number of residual evaluations        = %d\n" nre;
  printf "Number of Preconditioner evaluations  = %d\n" npevals;
  printf "Number of linear iterations           = %d\n" sli;
  printf "Number of error test failures         = %d\n" netf;
  printf "Number of precond solve fun called    = %d\n" nps

let main () =
  (* Set the number of threads to use *)
  let num_threads =
    if Array.length Sys.argv > 1
    then int_of_string Sys.argv.(1)
    else 1
  in

  let webdata = init_user_data num_threads
  and cc = RealArray.create neq
  and cp = RealArray.create neq
  and id = RealArray.create neq in
  set_initial_profiles webdata cc cp id;

  (* Set remaining inputs to IDAInit. *)
  let t0 = 0. in

  (* Wrap c and c' in nvectors.  Operations performed on the wrapped
     representation affect the originals c and c'.  *)
  let wcc = Nvector_openmp.wrap num_threads cc
  and wcp = Nvector_openmp.wrap num_threads cp
  in

  (* Call IDACreate and IDABand to initialize IDA including the linear
     solver. *)
  let maxl = 16 in
  let mem =
    Ida.(init Spils.(solver Iterative.(spgmr ~maxl wcc)
                            (prec_left ~setup:(precond webdata)
                                       (psolve webdata)))
      (SStolerances (rtol, atol))
      (resweb webdata) t0 wcc wcp)
  in
  webdata.ida_mem <- Some mem;
  let tout1 = 0.001 in
  Ida.calc_ic_ya_yd' mem ~varid:(Nvector_openmp.wrap num_threads id) tout1;

  (* Print heading, basic parameters, and initial values. *)
  print_header maxl rtol atol;
  print_output mem cc 0.;

  (* Loop over iout, call IDASolve (normal mode), print selected output. *)
  let tout = ref tout1 in
  for iout = 1 to nout do
    let (tret, retval) = Ida.solve_normal mem !tout wcc wcp in
    print_output mem cc tret;
    if iout < 3 then tout := !tout *. tmult
    else tout := !tout +. tadd
  done;

  print_final_stats mem;
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
