(*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2010/12/01 23:09:24 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem for KINSOL (parallel machine case) using the BBD
 * preconditioner.
 *
 * This example solves a nonlinear system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions. The dependent variable vector
 * is the following:
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
 * prey and the last np being predators. The number np is both
 * the number of prey and predator species. The coefficients a(i,j),
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
 * The various scalar parameters are set using define's or in
 * routine InitUserData.
 *
 * The boundary conditions are: normal derivative = 0, and the
 * initial guess is constant in x and y, although the final
 * solution is not.
 *
 * The PDEs are discretized by central differencing on a MX by
 * MY mesh.
 *
 * The nonlinear system is solved by KINSOL using the method
 * specified in the local variable globalstrat.
 *
 * The preconditioner matrix is a band-block-diagonal matrix
 * using the KINBBDPRE module. The half-bandwidths are as follows:
 *
 *   Difference quotient half-bandwidths mldq = mudq = 2*ns - 1
 *   Retained banded blocks have half-bandwidths mlkeep = mukeep = ns.
 *
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
 * ----------------------------------------------------------------------
 *  Run command line: mpirun -np N -machinefile machines kinFoodWeb_kry_bbd_p
 *  where N = NPEX * NPEY is the number of processors.
 * ----------------------------------------------------------------------
 *)

module Nvector = Nvector_parallel
module RealArray = Sundials.RealArray
module RealArray2 = Sundials.RealArray2
module LintArray = Sundials.LintArray
module Dense = Dls.ArrayDenseMatrix
module Bbd = Kinsol_bbd
open Bigarray
let local_array = Nvector_parallel.local_array

let printf = Printf.printf
let eprintf = Printf.eprintf
let subarray = Array1.sub
let slice_left = Array2.slice_left
let slice = Array1.sub
let unwrap = RealArray2.unwrap

let sqr x = x *. x
let nvwl2norm = Nvector.DataOps.n_vwl2norm

(* Problem Constants *)

let num_species =   6  (* must equal 2*(number of prey or predators)
                              number of prey = number of predators       *) 

let pi          = 3.1415926535898   (* pi *) 

let npex        = 2            (* number of processors in the x-direction  *)
let npey        = 2            (* number of processors in the y-direction  *)
let mxsub       = 10           (* number of x mesh points per subgrid      *)
let mysub       = 10           (* number of y mesh points per subgrid      *)
let mx          = npex*mxsub   (* number of mesh points in the x-direction *)
let my          = npey*mysub   (* number of mesh points in the y-direction *)
let nsmxsub     = num_species * mxsub
let nsmxsub2    = num_species * (mxsub+2)
let nsmysub     = num_species * mysub
let neq         = num_species*mx*my  (* number of equations in the system *)
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
let preyin      = 1.0    (* initial guess for prey concentrations. *)
let predin      = 30000.0(* initial guess for predator concs.      *)

(* User-defined vector access macro: IJ_Vptr *)

(* ij_vptr is defined in order to translate from the underlying 3D structure
   of the dependent variable vector to the 1D storage scheme for an N-vector.
   ij_vptr vv i j  returns a pointer to the location in vv corresponding to 
   indices is = 0, jx = i, jy = j.    *)
let ij_vptr_idx i j = i*num_species + j*nsmxsub
let ij_vptr vv i j = subarray vv (ij_vptr_idx i j) num_species

(* Type : UserData
   contains problem constants and extended array *)

type user_data = {
  acoef : RealArray2.data;
  bcoef : RealArray.t;
  cox   : RealArray.t;
  coy   : RealArray.t;
  rates : RealArray.t;
  cext  : RealArray.t;

  my_pe : int;
  comm  : Mpi.communicator;
  isubx : int;
  isuby : int;
}

(* Load problem constants in data *)

let dx = ax /. float(mx-1)
let dy = ay /. float(my-1)

let init_user_data my_pe comm =
  let np = num_species/2 in
  let dx2 = dx*.dx in
  let dy2 = dy*.dy in

  let isuby = my_pe/npex in
  let isubx = my_pe - isuby*npex in
  let data = {
      acoef = RealArray2.make_data num_species num_species;
      bcoef = RealArray.create num_species;
      cox = RealArray.create num_species;
      coy = RealArray.create num_species;

      rates = RealArray.create neq;
      cext  = RealArray.create (num_species * (mxsub+2)*(mysub+2));

      my_pe = my_pe;
      comm  = comm;
      isubx = isubx;
      isuby = isuby;
    } in

  (* Fill in the portion of acoef in the four quadrants, row by row *)
  for i = 0 to np - 1 do
    let a1 = subarray (slice_left data.acoef i) np np in
    let a2 = slice_left data.acoef (i + np) in
    let a3 = slice_left data.acoef i in
    let a4 = subarray (slice_left data.acoef (i + np)) np np in

    for j = 0 to np - 1 do
      a1.{j} <- -.gg;
      a2.{j} <-   ee;
      a3.{j} <-   zero;
      a4.{j} <-   zero
    done;

    (* and then change the diagonal elements of acoef to -AA *)
    data.acoef.{i, i} <- -.aa;
    data.acoef.{i+np, i+np} <- -.aa;

    data.bcoef.{i} <- bb;
    data.bcoef.{i+np} <- -.bb;

    data.cox.{i} <- dprey/.dx2;
    data.cox.{i+np} <- dpred/.dx2;

    data.coy.{i} <- dprey/.dy2;
    data.coy.{i+np} <- dpred/.dy2
  done;
  data

let blit (buf : RealArray.t) buf_offset (dst : RealArray.t) dst_offset len =
  for i = 0 to len-1 do
    dst.{dst_offset + i} <- buf.{buf_offset + i}
  done

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_string (RealArray.create 0) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_string (RealArray.create 1) []) 0
  - header_and_empty_array_size

let bytes x = header_and_empty_array_size + x * float_cell_size

(* Routine to send boundary data to neighboring PEs *)

let bsend comm my_pe isubx isuby dsizex dsizey udata =
  let buf = RealArray.create (num_species*mysub) in

  (* If isuby > 0, send data from bottom x-line of u *)
  if isuby <> 0 then Mpi.send (slice udata 0 dsizex) (my_pe-npex) 0 comm;

  (* If isuby < NPEY-1, send data from top x-line of u *)
  if isuby <> npey-1 then begin
    let offsetu = (mysub-1)*dsizex in
    Mpi.send (slice udata offsetu dsizex) (my_pe+npex) 0 comm
  end;

  (* If isubx > 0, send data from left y-line of u (via bufleft) *)
  if isubx <> 0 then begin
    for ly = 0 to mysub-1 do
      blit udata (ly*dsizex) buf (ly*num_species) num_species
    done;
    Mpi.send buf (my_pe-1) 0 comm
  end;

  (* If isubx < NPEX-1, send data from right y-line of u (via bufright) *)
  if isubx <> npex-1 then begin
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*num_species in
      let offsetu = offsetbuf*mxsub + (mxsub-1)*num_species in
      blit udata offsetu buf offsetbuf num_species
    done;
    Mpi.send buf (my_pe+1) 0 comm
  end
 
(* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*num_species*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. *)

let brecvpost comm my_pe isubx isuby dsizex dsizey =
  (* If isuby > 0, receive data for bottom x-line of cext *)
  let r0 = if isuby <> 0
           then Mpi.ireceive (bytes dsizex) (my_pe-npex) 0 comm
           else Mpi.null_request
  in
  (* If isuby < NPEY-1, receive data for top x-line of cext *)
  let r1 = if isuby <> npey-1
           then Mpi.ireceive (bytes dsizex) (my_pe+npex) 0 comm
           else Mpi.null_request
  in
  (* If isubx > 0, receive data for left y-line of cext (via bufleft) *)
  let r2 = if isubx <> 0
           then Mpi.ireceive (bytes dsizey) (my_pe-1) 0 comm
           else Mpi.null_request
  in
  (* If isubx < NPEX-1, receive data for right y-line of cext (via bufright) *)
  let r3 = if isubx <> npex-1
           then Mpi.ireceive (bytes dsizey) (my_pe+1) 0 comm
           else Mpi.null_request
  in
  Array.of_list [r0; r1; r2; r3]

(* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*num_species*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. *)

let brecvwait request isubx isuby dsizex cext =
  let dsizex2 = dsizex + 2*num_species in

  (* If isuby > 0, receive data for bottom x-line of cext *)
  if isuby <> 0 then begin
    let buf = (Mpi.wait_receive request.(0) : RealArray.t) in
    blit buf 0 cext num_species dsizex
  end;

  (* If isuby < NPEY-1, receive data for top x-line of cext *)
  if isuby <> npey-1 then begin
    let buf = (Mpi.wait_receive request.(1) : RealArray.t) in
    blit buf 0 cext (num_species*(1 + (mysub+1)*(mxsub+2))) dsizex
  end;

  (* If isubx > 0, receive data for left y-line of cext (via bufleft) *)
  if isubx <> 0 then begin
    let bufleft = (Mpi.wait_receive request.(2) : RealArray.t) in
    (* Copy the buffer to cext *)
    for ly = 0 to mysub - 1 do
      let offsetbuf = ly*num_species in
      let offsetue = (ly+1)*dsizex2 in
      blit bufleft offsetbuf cext offsetue num_species
    done
  end;

  (* If isubx < NPEX-1, receive data for right y-line of cext (via bufright) *)
  if isubx <> npex-1 then begin
    let bufright = (Mpi.wait_receive request.(3) : RealArray.t) in
    (* Copy the buffer to cext *)
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*num_species in
      let offsetue = (ly+2)*dsizex2 - num_species in
      blit bufright offsetbuf cext offsetue num_species
    done
  end

(* ccomm routine.  This routine performs all communication 
   between processors of data needed to calculate f. *)

let ccomm data udata =
  let comm    = data.comm
  and my_pe   = data.my_pe
  and isubx   = data.isubx
  and isuby   = data.isuby
  and cext    = data.cext
  in
  (* Start receiving boundary data from neighboring PEs *)
  let request = brecvpost comm my_pe isubx isuby nsmxsub nsmysub in

  (* Send data from boundary of local grid to neighboring PEs *)
  bsend comm my_pe isubx isuby nsmxsub nsmysub udata;

  (* Finish receiving boundary data from neighboring PEs *)
  brecvwait request isubx isuby nsmxsub cext
(* Interaction rate function routine *)

let web_rate data xx yy ((cxy : RealArray.t), cxy_off)
                        ((ratesxy : RealArray.t), ratesxy_off) =
  let acoef = data.acoef in
  let bcoef = data.bcoef in

  for i = 0 to num_species - 1 do
    ratesxy.{ratesxy_off + i} <- 0.0;
    for j = 0 to num_species - 1 do
      ratesxy.{ratesxy_off + i} <- ratesxy.{ratesxy_off + i}
                                      +. cxy.{cxy_off + j} *. acoef.{i, j}
    done
  done;
  
  let fac = one +. alpha *. xx *. yy in
  for i = 0 to num_species - 1 do
    ratesxy.{ratesxy_off + i} <- cxy.{cxy_off + i} *. (bcoef.{i}
                                    *. fac +. ratesxy.{ratesxy_off + i})
  done

(* System function for predator-prey system - calculation part *)

let func_local data (cdata, _, _) (fval, _, _) =
  
  (* Get subgrid indices, data sizes, extended work array cext *)
  let isubx = data.isubx in
  let isuby = data.isuby in
  let cext = data.cext in
  let rates = data.rates in

  (* Copy local segment of cc vector into the working extended array cext *)
  let offsetc = ref 0 in
  let offsetce = ref (nsmxsub2 + num_species) in
  for ly = 0 to mysub - 1 do
    for i = 0 to nsmxsub - 1 do
      cext.{!offsetce+i} <- cdata.{!offsetc+i}
    done;
    offsetc := !offsetc + nsmxsub;
    offsetce := !offsetce + nsmxsub2
  done;

  (* To facilitate homogeneous Neumann boundary conditions, when this is a
  boundary PE, copy data from the first interior mesh line of cc to cext   *)

  (* If isuby = 0, copy x-line 2 of cc to cext *)
  if isuby = 0 then
    for i = 0 to nsmxsub - 1 do
      cext.{num_species+i} <- cdata.{nsmxsub+i}
    done;

  (* If isuby = NPEY-1, copy x-line MYSUB-1 of cc to cext *)
  if isuby = npey-1 then begin
    let offsetc = ref ((mysub-2)*nsmxsub) in
    let offsetce = ref ((mysub+1)*nsmxsub2 + num_species) in
    for i = 0 to nsmxsub - 1 do
      cext.{!offsetce+i} <- cdata.{!offsetc+i}
    done
  end;

  (* If isubx = 0, copy y-line 2 of cc to cext *)
  if isubx = 0 then
    for ly = 0 to mysub - 1 do
      let offsetc = ref (ly*nsmxsub + num_species) in
      let offsetce = ref ((ly+1)*nsmxsub2) in
      for i = 0 to num_species - 1 do
        cext.{!offsetce+i} <- cdata.{!offsetc+i}
      done
    done;

  (* If isubx = NPEX-1, copy y-line MXSUB-1 of cc to cext *)
  if isubx = npex-1 then
    for ly = 0 to mysub - 1 do
      let offsetc = ref ((ly+1)*nsmxsub - 2*num_species) in
      let offsetce = ref ((ly+2)*nsmxsub2 - num_species) in
      for i = 0 to num_species - 1 do
        cext.{!offsetce+i} <- cdata.{!offsetc+i}
      done
    done;
  
  (* Loop over all mesh points, evaluating rate arra at each point *)
  let delx = dx in
  let dely = dy in
  let shifty = (mxsub+2)*num_species in
  
  for jy = 0 to mysub - 1 do
    let yy = dely*.float (jy + isuby * mysub) in

    for jx = 0 to mxsub - 1 do
      let xx = delx *. float (jx + isubx * mxsub) in

      let off = ij_vptr_idx jx jy in
      web_rate data xx yy (cdata, off) (rates, off);

      let offsetc = (jx+1)*num_species + (jy+1)*nsmxsub2 in
      let offsetcd = offsetc - shifty in
      let offsetcu = offsetc + shifty in
      let offsetcl = offsetc - num_species in
      let offsetcr = offsetc + num_species in
      
      for is = 0 to num_species - 1 do
        
        (* differencing in x *)
        let dcydi = cext.{offsetc+is}  -. cext.{offsetcd+is} in
        let dcyui = cext.{offsetcu+is} -. cext.{offsetc+is} in
        
        (* differencing in y *)
        let dcxli = cext.{offsetc+is}  -. cext.{offsetcl+is} in
        let dcxri = cext.{offsetcr+is} -. cext.{offsetc+is} in
        
        (* compute the value at xx , yy *)
        fval.{off + is} <- data.coy.{is} *. (dcyui -. dcydi)
                           +. data.cox.{is} *. (dcxri -. dcxli)
                           +. rates.{off + is}
      done (* end of is loop *)
    done (* end of jx loop *)
  done (* end of jy loop *)

(*
 * System function routine.  Evaluate f(cc).  First call ccomm to do
 * communication of subgrid boundary data into cext.  Then calculate f
 * by a call to func_local.
 *)

let func data ((cdata, _, _) as cc) ((fvdata, _, _) as fv) =
  (* Call ccomm to do inter-processor communicaiton *)
  ccomm data cdata;

  (* Call func_local to calculate all right-hand sides *)
  func_local data cc fv

(* Set initial conditions in cc *)
let set_initial_profiles cc sc =
  (* Load initial profiles into cc and sc vector. *)
  for jy = 0 to mysub - 1 do
    for jx = 0 to mxsub - 1 do
      let idx = ij_vptr_idx jx jy in
      for i = 0 to num_species/2 - 1 do
        cc.{idx + i} <- preyin;
        sc.{idx + i} <- one
      done;
      for i = num_species/2 to num_species - 1 do
        cc.{idx + i} <- predin;
        sc.{idx + i} <- 0.00001
      done
    done
  done

(* Print first lines of output (problem description) *)
let print_header globalstrategy maxl maxlrst
                 mudq mldq mukeep mlkeep fnormtol scsteptol =
  printf "\nPredator-prey test problem--  KINSol (parallel-BBD version)\n\n";
  printf "Mesh dimensions = %d X %d\n" mx my;
  printf "Number of species = %d\n" num_species;
  printf "Total system size = %d\n\n" neq;
  printf "Subgrid dimensions = %d X %d\n" mxsub mysub;
  printf "Processor array is %d X %d\n\n" npex npey;
  printf "Flag globalstrategy = %d (0 = None, 1 = Linesearch)\n"
         (if globalstrategy = Kinsol.LineSearch then 1 else 0);
  printf "Linear solver is SPGMR with maxl = %d, maxlrst = %d\n" maxl maxlrst;
  printf "Preconditioning uses band-block-diagonal matrix from KINBBDPRE\n";
  printf "  Difference quotient half-bandwidths: mudq = %d, mldq = %d\n"
                                                                     mudq mldq;
  printf "  Retained band block half-bandwidths: mukeep = %d, mlkeep = %d\n"
                                                                mukeep mlkeep;
  printf "Tolerance parameters:  fnormtol = %g   scsteptol = %g\n"
         fnormtol scsteptol;
  printf "\nInitial profile of concentration\n";
  printf "At all mesh points:  %g %g %g   %g %g %g\n"
         preyin preyin preyin predin predin predin

(* Print sample of current cc values *)
let print_output my_pe comm cc =
  let npelast = npex*npey - 1 in
  let ct = local_array cc in
  let i0 = num_species*(mxsub*mysub-1) in
  
  (* Send the cc values (for all species) at the top right mesh point to PE 0 *)
  if my_pe = npelast then begin
    if npelast <> 0 then Mpi.send (slice ct i0 num_species) 0 0 comm
  end;
  
  (* On PE 0, receive the cc values at top right, then print performance data 
     and sampled solution values *)
  if my_pe = 0 then begin
    let tempc =
      if npelast <> 0 then (Mpi.receive npelast 0 comm : RealArray.t)
      else RealArray.init num_species (fun is -> ct.{i0 + is})
    in
    
    printf "\nAt bottom left:";
    for is = 0 to num_species - 1 do
      if (is mod 6)*6 = is then printf "\n";
      printf " %g" ct.{is}
    done;
    
    printf "\n\nAt top right:";
    for is = 0 to num_species - 1 do
      if (is mod 6)*6 = is then printf "\n";
      printf " %g" tempc.{is}
    done;
    printf "\n\n"
  end

(* Print final statistics contained in iopt *)
let print_final_stats kmem =
  let open Kinsol in
  let nni   = get_num_nonlin_solv_iters kmem in
  let nfe   = get_num_func_evals kmem in
  let nli   = Spils.get_num_lin_iters kmem in
  let npe   = Spils.get_num_prec_evals kmem in
  let nps   = Spils.get_num_prec_solves kmem in
  let ncfl  = Spils.get_num_conv_fails kmem in
  let nfeSG = Spils.get_num_func_evals kmem in
  printf "Final Statistics.. \n";
  printf "nni    = %5d    nli   = %5d\n" nni nli;
  printf "nfe    = %5d    nfeSG = %5d\n" nfe nfeSG;
  printf "nps    = %5d    npe   = %5d     ncfl  = %5d\n" nps npe ncfl

(* MAIN PROGRAM *)
let main () =
  let globalstrategy = Kinsol.Newton in

  let comm    = Mpi.comm_world in
  let my_pe   = Mpi.comm_rank comm in
  let npes    = Mpi.comm_size comm in
  let npelast = npex*npey-1 in

  if npes <> npex*npey then begin
    if my_pe = 0 then
      eprintf "\nMPI_ERROR(0): npes=%d is not equal to NPEX*NPEY=%d\n"
              npes (npex*npey);
    exit 1
  end;

  (* Allocate memory, and set problem data, initial values, tolerances *)

  (* Set local length *)
  let local_N = num_species*mxsub*mysub in

  let data = init_user_data my_pe comm in

  (* Create serial vectors of length NEQ *)
  let cc = Nvector.make local_N neq comm 0.0 in
  let sc = Nvector.make local_N neq comm 0.0 in
  set_initial_profiles (local_array cc) (local_array sc);

  let fnormtol  = ftol in
  let scsteptol = stol in

  let maxl = 20 in
  let maxlrst = 2 in

  (* Call KINBBDPrecInit to initialize and allocate memory for the
     band-block-diagonal preconditioner, and specify the local and
     communication functions func_local and gcomm=NULL (all communication
     needed for the func_local is already done in func). *)
  let mudq = 2*num_species - 1 in
  let mldq = 2*num_species - 1 in
  let mukeep = num_species in
  let mlkeep = num_species in
  let kmem =
    Kinsol.(init
        ~linsolv:(Spils.spgmr ~maxl:maxl ~max_restarts:maxlrst
                              (Kinsol_bbd.prec_right
                                  Bbd.({ mudq; mldq; mukeep; mlkeep; })
                                  (func_local data)))
        (func data) cc) in
  Kinsol.set_constraints kmem (Nvector.make local_N neq comm 0.0);
  Kinsol.set_func_norm_tol kmem fnormtol;
  Kinsol.set_scaled_step_tol kmem scsteptol;

  (* Print out the problem size, solution parameters, initial guess. *)
  if my_pe = 0 then print_header globalstrategy maxl maxlrst
                                 mudq mldq mukeep mlkeep fnormtol scsteptol;

  (* Call KINSol and print output concentration profile *)
  ignore (Kinsol.solve kmem           (* KINSol memory block *)
                 cc             (* initial guess on input; solution vector *)
                 globalstrategy (* global stragegy choice *)
                 sc             (* scaling vector, for the variable cc *)
                 sc);           (* scaling vector for function values fval *)

  if my_pe = 0 then printf("\n\nComputed equilibrium species concentrations:\n");
  if my_pe = 0 || my_pe = npelast then print_output my_pe comm cc;

  (* Print final statistics and free memory *)  
  if my_pe = 0 then print_final_stats kmem

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
