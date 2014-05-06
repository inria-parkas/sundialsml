(*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/01 23:08:49 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, May 2014.
 * -----------------------------------------------------------------
 * Example (serial):
 *
 * This example solves a nonlinear system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion
 * on the unit square in two dimensions. The dependent variable
 * vector is the following:
 * 
 *       1   2         ns
 * c = (c , c ,  ..., c  )     (denoted by the variable cc)
 * 
 * and the PDE's are as follows:
 *
 *                    i       i
 *         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)
 *                    xx      yy         i
 *
 *   where
 *
 *                   i             ns         j
 *   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The number np is both the
 * number of prey and predator species. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *   a(i,i) = -AA   (all i)
 *   a(i,j) = -GG   (i <= np , j >  np)
 *   a(i,j) =  EE   (i >  np,  j <= np)
 *   b(i) = BB * (1 + alpha * x * y)   (i <= np)
 *   b(i) =-BB * (1 + alpha * x * y)   (i >  np)
 *   d(i) = DPREY   (i <= np)
 *   d(i) = DPRED   ( i > np)
 *
 * The various scalar parameters are set earlier.
 *
 * The boundary conditions are: normal derivative = 0, and the
 * initial guess is constant in x and y, but the final solution
 * is not.
 *
 * The PDEs are discretized by central differencing on an MX by
 * MY mesh.
 * 
 * The nonlinear system is solved by KINSOL using the method
 * specified in local variable globalstrat.
 *
 * The preconditioner matrix is a block-diagonal matrix based on
 * the partial derivatives of the interaction terms f only.
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 * References:
 *
 * 1. Peter N. Brown and Youcef Saad,
 *    Hybrid Krylov Methods for Nonlinear Systems of Equations
 *    LLNL report UCRL-97645, November 1987.
 *
 * 2. Peter N. Brown and Alan C. Hindmarsh,
 *    Reduced Storage Matrix Methods in Stiff ODE systems,
 *    Lawrence Livermore National Laboratory Report  UCRL-95088,
 *    Rev. 1, June 1987, and  Journal of Applied Mathematics and
 *    Computation, Vol. 31 (May 1989), pp. 40-91. (Presents a
 *    description of the time-dependent version of this test
 *    problem.)
 * -----------------------------------------------------------------
 *)

module Kinsol = Kinsol_serial
module Carray = Sundials.Carray
module Dense = Dls.ArrayDenseMatrix

type real_array = Sundials.real_array

let printf = Printf.printf
let subarray = Bigarray.Array1.sub
let slice_left = Bigarray.Array2.slice_left
let unwrap = Sundials.Realarray2.unwrap
let nvwl2norm =
  match Nvector_array.Bigarray.array_nvec_ops.Nvector.Mutable.nvwl2norm with
  | Some fn -> fn
  | None -> failwith "nvwl2norm not found!"

let ith v i = v.{i - 1}
let set_ith v i e = v.{i - 1} <- e

(* Problem Constants *)

let num_species =   6  (* must equal 2*(number of prey or predators)
                              number of prey = number of predators       *) 

let pi          = 3.1415926535898   (* pi *) 

let mx          = 8                 (* MX = number of x mesh points *)
let my          = 8                 (* MY = number of y mesh points *)
let nsmx        = (num_species * mx)
let neq         = (nsmx * my)       (* number of equations in the system *)
let aa          = 1.0    (* value of coefficient AA in above eqns *)
let ee          = 10000. (* value of coefficient EE in above eqns *)
let gg          = 0.5e-6 (* value of coefficient GG in above eqns *)
let bb          = 1.0    (* value of coefficient BB in above eqns *)
let dprey       = 1.0    (* value of coefficient dprey above *)
let dpred       = 0.5    (* value of coefficient dpred above *)
let alpha       = 1.0    (* value of coefficient alpha above *)
let ax          = 1.0    (* total range of x variable *)
let ay          = 1.0    (* total range of y variable *)
let ftol        = 1.e-7  (* ftol tolerance *)
let stol        = 1.e-13 (* stol tolerance *)
let thousand    = 1000.0 (* one thousand *)
let zero        = 0.0    (* 0. *)
let one         = 1.0    (* 1. *)
let two         = 2.0    (* 2. *)
let preyin      = 1.0    (* initial guess for prey concentrations. *)
let predin      = 30000.0(* initial guess for predator concs.      *)

(* User-defined vector access macro: IJ_Vptr *)

(* ij_vptr is defined in order to translate from the underlying 3D structure
   of the dependent variable vector to the 1D storage scheme for an N-vector.
   ij_vptr vv i j  returns a pointer to the location in vv corresponding to 
   indices is = 0, jx = i, jy = j.    *)
let ij_vptr_idx i j = i*num_species + j*nsmx
let ij_vptr vv i j = subarray vv (ij_vptr_idx i j) num_species

(* Type : UserData 
   contains preconditioner blocks, pivot arrays, and problem constants *)

let p =
  Array.init mx (fun jx ->
    Array.init my (fun jy ->
      Dense.make num_species num_species
    ))

let pivot =
  Array.init mx (fun jx ->
    Array.init my (fun jy ->
      let v = Sundials.make_lint_array num_species in
      Bigarray.Array1.fill v 0;
      v
    ))

let acoef = Sundials.make_real_array2 num_species num_species
let bcoef = Carray.create num_species
let cox = Carray.create num_species
let coy = Carray.create num_species

let rates = Carray.create neq

(* Load problem constants in data *)

let dx = ax /. float(mx-1)
let dy = ay /. float(my-1)
let uround = Sundials.unit_roundoff
let sqruround = sqrt(uround)

let init_user_data =
  let np = num_species/2 in
  let dx2 = dx*.dx in
  let dy2 = dy*.dy in

  (* Fill in the portion of acoef in the four quadrants, row by row *)
  for i = 0 to np - 1 do
    let a1 = subarray (slice_left acoef i) np np in
    let a2 = slice_left acoef (i + np) in
    let a3 = slice_left acoef i in
    let a4 = subarray (slice_left acoef (i + np)) np np in

    for j = 0 to np - 1 do
      a1.{j} <- -.gg;
      a2.{j} <-   ee;
      a3.{j} <-   zero;
      a4.{j} <-   zero
    done;

    (* and then change the diagonal elements of acoef to -AA *)
    acoef.{i, i} <- -.aa;
    acoef.{i+np, i+np} <- -.aa;

    bcoef.{i} <- bb;
    bcoef.{i+np} <- -.bb;

    cox.{i} <- dprey/.dx2;
    cox.{i+np} <- dpred/.dx2;

    coy.{i} <- dprey/.dy2;
    coy.{i+np} <- dpred/.dy2
  done

(* Dot product routine for realtype arrays *)
let dot_prod size x1 x2 =
  let temp =ref zero in
  for i = 0 to size - 1 do
    temp := !temp +. x1.{i} *. x2.{i}
  done;
  !temp

(* Interaction rate function routine *)

let web_rate xx yy cxy ratesxy =
  for i = 0 to num_species - 1 do
    ratesxy.{i} <- dot_prod num_species cxy (slice_left acoef i)
  done;
  
  let fac = one +. alpha *. xx *. yy in
  for i = 0 to num_species - 1 do
    ratesxy.{i} <- cxy.{i} *. (bcoef.{i} *. fac +. ratesxy.{i})
  done

(* System function for predator-prey system *)

let func cc fval =
  let delx = dx in
  let dely = dy in
  
  (* Loop over all mesh points, evaluating rate array at each point*)
  for jy = 0 to my - 1 do
    let yy = dely *. float(jy) in

    (* Set lower/upper index shifts, special at boundaries. *)
    let idyl = if jy <> 0    then nsmx else -nsmx in
    let idyu = if jy <> my-1 then nsmx else -nsmx in
    
    for jx = 0 to mx - 1 do
      let xx = delx *. float(jx) in

      (* Set left/right index shifts, special at boundaries. *)
      let idxl = if jx <>  0   then num_species else -num_species in
      let idxr = if jx <> mx-1 then num_species else -num_species in

      let cxy = ij_vptr cc jx jy in
      let cxy_ij = ij_vptr_idx jx jy in
      let rxy = ij_vptr rates jx jy in
      let fxy = ij_vptr fval jx jy in

      (* Get species interaction rate array at (xx,yy) *)
      web_rate xx yy cxy rxy;

      for is = 0 to num_species - 1 do
        (* Differencing in x direction *)
        let dcyli = cxy.{is} -. cc.{cxy_ij - idyl + is} in
        let dcyui = cc.{cxy_ij + idyu + is} -. cxy.{is} in
        
        (* Differencing in y direction *)
        let dcxli = cxy.{is} -. cc.{cxy_ij - idxl + is} in
        let dcxri = cc.{cxy_ij + idxr+is} -. cxy.{is} in
        
        (* Compute the total rate value at (xx,yy) *)
        fxy.{is} <- coy.{is} *. (dcyui -. dcyli)
                    +. cox.{is} *. (dcxri -. dcxli) +. rxy.{is}

      done (* end of is loop *)
    done (* end of jx loop *)
  done (* end of jy loop *)

(* Preconditioner setup routine. Generate and preprocess P. *)
let prec_setup_bd { Kinsol.jac_u=cc;
                    Kinsol.jac_fu=fval;
                    Kinsol.jac_tmp=(vtemp1, vtemp2)}
                  { Kinsol.uscale=cscale;
                    Kinsol.fscale=fscale } =
  let perturb_rates = Sundials.make_real_array num_species in
  
  let delx = dx in
  let dely = dy in

  let fac = nvwl2norm fval fscale in
  let r0 = thousand *. uround *. fac *. float(neq) in
  let r0 = if r0 = zero then one else r0 in
  
  (* Loop over spatial points; get size NUM_SPECIES Jacobian block at each *)
  for jy = 0 to my - 1 do
    let yy = float(jy) *. dely in
    
    for jx = 0 to mx - 1 do
      let xx = float(jx) *. delx in
      let pxy = p.(jx).(jy) in
      let cxy = ij_vptr cc jx jy in
      let scxy = ij_vptr cscale jx jy in
      let ratesxy = ij_vptr rates jx jy in
      
      (* Compute difference quotients of interaction rate fn. *)
      for j = 0 to num_species - 1 do
        let csave = cxy.{j} in  (* Save the j,jx,jy element of cc *)
        let r = max (sqruround *. abs_float csave) (r0/.scxy.{j}) in
        cxy.{j} <- cxy.{j} +. r; (* Perturb the j,jx,jy element of cc *)
        let fac = one/.r in
        web_rate xx yy cxy perturb_rates;
        
        (* Restore j,jx,jy element of cc *)
        cxy.{j} <- csave;
        
        (* Load the j-th column of difference quotients *)
        let pxycol = slice_left (unwrap pxy) j in
        for i = 0 to num_species - 1 do
          pxycol.{i} <- (perturb_rates.{i} -. ratesxy.{i}) *. fac
        done
      done; (* end of j loop *)
      
      (* Do LU decomposition of size NUM_SPECIES preconditioner block *)
      Dense.getrf pxy pivot.(jx).(jy)
    done (* end of jx loop *)
  done (* end of jy loop *)
  
(* Preconditioner solve routine *)
let prec_solve_bd { Kinsol.jac_u=cc;
                    Kinsol.jac_fu=fval;
                    Kinsol.jac_tmp=ftem}
                  { Kinsol.uscale=cscale;
                    Kinsol.fscale=fscale }
                  vv =
  for jx = 0 to mx - 1 do
    for jy = 0 to my - 1 do
      (* For each (jx,jy), solve a linear system of size NUM_SPECIES.
         vxy is the address of the corresponding portion of the vector vv;
         Pxy is the address of the corresponding block of the matrix P;
         piv is the address of the corresponding block of the array pivot. *)
      let vxy = ij_vptr vv jx jy in
      let pxy = p.(jx).(jy) in
      let piv = pivot.(jx).(jy) in
      Dense.getrs pxy piv vxy
    done (* end of jy loop *)
  done (* end of jx loop *)

(* Set initial conditions in cc *)
let set_initial_profiles cc sc =
  (* Load initial profiles into cc and sc vector. *)
  for jy = 0 to my - 1 do
    for jx = 0 to mx - 1 do
      let cloc = ij_vptr cc jx jy in
      let sloc = ij_vptr sc jx jy in
      for i = 0 to num_species/2 - 1 do
        cloc.{i} <- preyin;
        sloc.{i} <- one
      done;
      for i = num_species/2 to num_species - 1 do
        cloc.{i} <- predin;
        sloc.{i} <- 0.00001
      done
    done
  done

(* Print first lines of output (problem description) *)
let print_header globalstrategy maxl maxlrst fnormtol scsteptol =
  printf "\nPredator-prey test problem --  KINSol (serial version)\n\n";
  printf "Mesh dimensions = %d X %d\n" mx my;
  printf "Number of species = %d\n" num_species;
  printf "Total system size = %d\n\n" neq;
  printf "Flag globalstrategy = %d (0 = None, 1 = Linesearch)\n"
         (if globalstrategy then 1 else 0);
  printf "Linear solver is SPGMR with maxl = %d, maxlrst = %d\n" maxl maxlrst;
  printf "Preconditioning uses interaction-only block-diagonal matrix\n";
  printf "Positivity constraints imposed on all components \n";
  printf "Tolerance parameters:  fnormtol = %g   scsteptol = %g\n"
         fnormtol scsteptol;
  printf "\nInitial profile of concentration\n";
  printf "At all mesh points:  %g %g %g   %g %g %g\n"
         preyin preyin preyin predin predin predin

(* Print sampled values of current cc *)
let print_output cc =
  let jy = 0 in
  let jx = 0 in
  let ct = ij_vptr cc jx jy in
  printf "\nAt bottom left:";

  (* Print out lines with up to 6 values per line *)
  for is = 0 to num_species - 1 do
    if ((is mod 6)*6 = is) then printf "\n";
    printf " %g" ct.{is}
  done;
  
  let jy = my-1 in
  let jx = mx-1 in
  let ct = ij_vptr cc jx jy in
  printf("\n\nAt top right:");

  (* Print out lines with up to 6 values per line *)
  for is = 0 to num_species - 1 do
    if ((is mod 6)*6 = is) then printf("\n");
    printf " %g" ct.{is}
  done;
  printf("\n\n")

(* Print final statistics contained in iopt *)
let print_final_stats kmem =
  let nni = Kinsol.get_num_nonlin_solv_iters kmem in
  let nfe = Kinsol.get_num_func_evals kmem in
  let nli = Kinsol.Spils.get_num_lin_iters kmem in
  let npe = Kinsol.Spils.get_num_prec_evals kmem in
  let nps = Kinsol.Spils.get_num_prec_solves kmem in
  let ncfl = Kinsol.Spils.get_num_conv_fails kmem in
  let nfeSG = Kinsol.Spils.get_num_func_evals kmem in
  printf "Final Statistics.. \n";
  printf "nni    = %5d    nli   = %5d\n" nni nli;
  printf "nfe    = %5d    nfeSG = %5d\n" nfe nfeSG;
  printf "nps    = %5d    npe   = %5d     ncfl  = %5d\n" nps npe ncfl

(* MAIN PROGRAM *)
let main () =
  let globalstrategy = false in

  (* Create serial vectors of length NEQ *)
  let cc = Carray.create neq in
  let sc = Carray.create neq in
  set_initial_profiles cc sc;

  let fnormtol  = ftol in
  let scsteptol = stol in

  let maxl = 15 in
  let maxlrst = 2 in

  (* Call KINCreate/KINInit to initialize KINSOL using the linear solver
     KINSPGMR with preconditioner routines prec_setup_bd
     and prec_solve_bd. *)
  let kmem = Kinsol.init
              (Kinsol.Spgmr (Some maxl, Some maxlrst, {
                             Kinsol.prec_setup_fn=Some prec_setup_bd;
                             Kinsol.prec_solve_fn=Some prec_solve_bd;
                             Kinsol.jac_times_vec_fn=None; }))
              func cc in
  Kinsol.set_constraints kmem (Carray.init neq two);
  Kinsol.set_func_norm_tol kmem (Some fnormtol);
  Kinsol.set_scaled_step_tol kmem (Some scsteptol);

  (* Print out the problem size, solution parameters, initial guess. *)
  print_header globalstrategy maxl maxlrst fnormtol scsteptol;

  (* Call KINSol and print output concentration profile *)
  ignore (Kinsol.solve kmem           (* KINSol memory block *)
                 cc             (* initial guess on input; solution vector *)
                 globalstrategy (* global stragegy choice *)
                 sc             (* scaling vector, for the variable cc *)
                 sc);           (* scaling vector for function values fval *)

  printf("\n\nComputed equilibrium species concentrations:\n");
  print_output cc;

  (* Print final statistics and free memory *)  
  print_final_stats kmem

let _ = main ()
let _ = Gc.compact ()

