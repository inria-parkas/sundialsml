(*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/01 23:03:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Aug 2014.
 * -----------------------------------------------------------------
 * Example program for IDA: Food web, parallel, GMRES, user
 * preconditioner.
 *
 * This example program for IDA uses IDASPGMR as the linear solver.
 * It is written for a parallel computer system and uses a
 * block-diagonal preconditioner (setup and solve routines) for the
 * IDASPGMR package. It was originally run on a Sun SPARC cluster
 * and used MPICH.
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
 *   0 = d(i)*(c   +  c  )  +  R  (x,y,c)   (i = np+1,...,ns)
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
 *   a(i,i) = -AA  (all i)
 *   a(i,j) = -GG  (i <= np , j >  np)
 *   a(i,j) =  EE  (i >  np,  j <= np)
 *   all other a(i,j) = 0
 *   b(i) = BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))  (i <= np)
 *   b(i) =-BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))  (i  > np)
 *   d(i) = DPREY  (i <= np)
 *   d(i) = DPRED  (i > np)
 *
 * Note: The above equations are written in 1-based indices,
 * whereas the code has 0-based indices, being written in C.
 *
 * The various scalar parameters required are set using '#define'
 * statements or directly in routine InitUserData. In this program,
 * np = 1, ns = 2. The boundary conditions are homogeneous Neumann:
 * normal derivative  =  0.
 *
 * A polynomial in x and y is used to set the initial values of the
 * first np variables (the prey variables) at each x,y location,
 * while initial values for the remaining (predator) variables are
 * set to a flat value, which is corrected by IDACalcIC.
 *
 * The PDEs are discretized by central differencing on a MX by MY
 * mesh, and so the system size Neq is the product
 * MX * MY * NUM_SPECIES. The system is actually implemented on
 * submeshes, processor by processor, with an MXSUB by MYSUB mesh
 * on each of NPEX * NPEY processors.
 *
 * The DAE system is solved by IDA using the IDASPGMR linear
 * solver, which uses the preconditioned GMRES iterative method to
 * solve linear systems. The precondtioner supplied to IDASPGMR is
 * the block-diagonal part of the Jacobian with ns by ns blocks
 * arising from the reaction terms only. Output is printed at
 * t = 0, .001, .01, .1, .4, .7, 1.
 * -----------------------------------------------------------------
 * References:
 * [1] Peter N. Brown and Alan C. Hindmarsh,
 *     Reduced Storage Matrix Methods in Stiff ODE systems,
 *     Journal of Applied Mathematics and Computation, Vol. 31
 *     (May 1989), pp. 40-91.
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
module LintArray = Sundials.LintArray

let fprintf = Printf.fprintf
let printf = Printf.printf

let vconst = Nvector_parallel.DataOps.n_vconst
let vscale = Nvector_parallel.DataOps.n_vscale

let slice = Bigarray.Array1.sub

let blit (buf : RealArray.t) buf_offset (dst : RealArray.t) dst_offset len =
  for i = 0 to len-1 do
    dst.{dst_offset + i} <- buf.{buf_offset + i}
  done

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 0) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 1) []) 0
  - header_and_empty_array_size

let bytes x = header_and_empty_array_size + x * float_cell_size

(* Problem Constants. *)

let nprey =       1        (* Number of prey (= number of predators). *)
let num_species = 2*nprey

let pi =          3.1415926535898   (* pi *)
let fourpi =      (4.0*.pi)          (* 4 pi *)

let mxsub =       10    (* Number of x mesh points per processor subgrid *)
let mysub =       10    (* Number of y mesh points per processor subgrid *)
let npex =        2     (* Number of subgrids in the x direction *)
let npey =        2     (* Number of subgrids in the y direction *)
let mx =          (mxsub*npex)      (* mx = number of x mesh points *)
let my =          (mysub*npey)      (* my = number of y mesh points *)
let nsmxsub =     (num_species * mxsub)
let neq =         (num_species*mx*my) (* Number of equations in system *)
let aa =          1.0    (* Coefficient in above eqns. for a *)
let ee =          10000. (* Coefficient in above eqns. for a *)
let gg =          0.5e-6 (* Coefficient in above eqns. for a *)
let bb =          1.0    (* Coefficient in above eqns. for b *)
let dprey =       1.0    (* Coefficient in above eqns. for d *)
let dpred =       0.05   (* Coefficient in above eqns. for d *)
let alpha =       50.    (* Coefficient alpha in above eqns. *)
let beta =        1000.  (* Coefficient beta in above eqns. *)
let ax =          1.0    (* Total range of x variable *)
let ay =          1.0    (* Total range of y variable *)
let rtol =        1.e-5  (*  rtol tolerance *)
let atol =        1.e-5  (*  atol tolerance *)
let nout =        6
let tmult =       10.0   (* Multiplier for tout values *)
let tadd =        0.3    (* Increment for tout values *)

(* User-defined vector accessor macro IJ_Vptr. *)

(* IJ_Vptr is defined in order to express the underlying 3-d structure of the
   dependent variable vector from its underlying 1-d storage (an N_Vector).
   IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to
   species index is = 0, x-index ix = i, and y-index jy = j.                *)

let ij_vptr_idx i j = i*num_species + j*nsmxsub

let ij_vptr (local, _, _) i j =
  let offset = ij_vptr_idx i j in
  slice local offset (RealArray.length local - offset)

(* Type: UserData.  Contains problem constants, preconditioner data, etc. *)

type user_data =
  {
    ns : int;
    np : int;
    thispe : int;
    npes : int;
    ixsub : int; jysub : int;
    npex : int; npey : int;
    mxsub : int; mysub : int;
    nsmxsub : int; nsmxsub2 : int;
    dx : float; dy : float;
    acoef : RealArray2.t;
    cox : RealArray.t;                  (* size = num_species *)
    coy : RealArray.t;                  (* size = num_species *)
    bcoef : RealArray.t;                (* size = num_species *)
    rhs : RealArray.t;                  (* size = num_species *)
    cext : RealArray.t;                 (* size = (mxsub+2)*(mysub+2)
                                                  *num_species *)
    comm : Mpi.communicator;
    rates : Nvector_parallel.data;
    pp : RealArray2.t array array; (* size = mxsub, mysub, variable *)
    pivot : LintArray.t array array; (* size = mxsub, mysub, variable *)
    ewt : Nvector_parallel.t;
    mutable ida_mem :
      (Nvector_parallel.data, Nvector_parallel.kind) Ida.session option;
  }

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(*
 * AllocUserData: Allocate memory for data structure of type UserData.
 *)
(*
 * InitUserData: Load problem constants in webdata (of type UserData).
 *)

let alloc_init_user_data comm local_N system_size thispe npes =
  let rates = Nvector_parallel.make local_N system_size comm 0.
  and pp =
    Array.init mxsub (fun ix ->
        Array.init mysub (fun jy ->
            RealArray2.create num_species num_species))
  and pivot =
    Array.init mxsub (fun ix ->
        Array.init mysub (fun jy ->
            LintArray.create num_species))
  in

  let acoef = RealArray2.create num_species num_species
  and ewt = Nvector_parallel.make local_N system_size comm 0.
  in

  let jysub = thispe / npex in
  let ixsub = thispe - jysub*npex in
  let ns = num_species in
  let np = nprey in
  let dx = ax/.float_of_int(mx-1) in
  let dy = ay/.float_of_int(my-1) in
  let nsmxsub = mxsub * num_species in
  let nsmxsub2 = (mxsub+2)*num_species in

  (* Set up the coefficients a and b plus others found in the equations. *)
  let dx2 = dx*.dx and dy2 = dy*.dy in

  let bcoef = RealArray.create num_species in
  let cox   = RealArray.create num_species in
  let coy   = RealArray.create num_species in
  let rhs   = RealArray.create num_species in
  let cext  = RealArray.create ((mxsub+2)*(mysub+2)*num_species) in

  for i = 0 to np-1 do
    (*  Fill in the portion of acoef in the four quadrants, row by row. *)
    for j = 0 to np-1 do
      RealArray2.set acoef (np+j) i      (-.gg);
      RealArray2.set acoef j      (i+np) ee;
      RealArray2.set acoef j      i      0.0;
      RealArray2.set acoef (np+j) (i+np) 0.0
    done;

    (* Reset the diagonal elements of acoef to -aa. *)
    RealArray2.set acoef i      i      (-.aa);
    RealArray2.set acoef (i+np) (i+np) (-.aa);

    (* Set coefficients for b and diffusion terms. *)
    bcoef.{i} <- bb; bcoef.{i+np} <- (-.bb);
    cox.{i} <- dprey/.dx2; cox.{i+np} <- dpred/.dx2;
    coy.{i} <- dprey/.dy2; coy.{i+np} <- dpred/.dy2
  done;

  {
    ns = ns; np = np; thispe = thispe; npes = npes;
    ixsub = ixsub; jysub = jysub;
    npex = npex; npey = npey; mxsub = mxsub; mysub = mysub;
    nsmxsub = nsmxsub; nsmxsub2 = nsmxsub2;
    dx = dx; dy = dy; acoef = acoef; cox = cox; coy = coy; bcoef = bcoef;
    rhs = rhs; cext = cext;
    comm = comm; rates = Nvector.unwrap rates;
    pp = pp; pivot = pivot; ewt = ewt;
    ida_mem = None;
  }

(*
 * Print first lines of output (problem description)
 *)

let spgmr = match Sundials.sundials_version with
  | 2,_,_ -> "IDASPGMR"
  | _ -> "SUNSPGMR"


let print_header system_size maxl rtol atol =
  printf "\nidaFoodWeb_kry_p: Predator-prey DAE parallel example problem for IDA \n\n";
  printf "Number of species ns: %d" num_species;
  printf "     Mesh dimensions: %d x %d" mx my;
  printf "     Total system size: %d\n" system_size;
  printf "Subgrid dimensions: %d x %d" mxsub mysub;
  printf "     Processor array: %d x %d\n" npex npey;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Linear solver: %s     Max. Krylov dimension maxl: %d\n" spgmr maxl;
  printf "Preconditioner: block diagonal, block size ns,";
  printf " via difference quotients\n";
  printf "CalcIC called to correct initial predator concentrations \n\n";

  printf "-----------------------------------------------------------\n";
  printf "  t        bottom-left  top-right";
  printf "    | nst  k      h\n";
  printf "-----------------------------------------------------------\n\n"

(*
 * PrintOutput: Print output values at output time t = tt.
 * Selected run statistics are printed.  Then values of c1 and c2
 * are printed for the bottom left and top right grid points only.
 * (NOTE: This routine is specific to the case NUM_SPECIES = 2.)
 *)

let print_output webdata comm mem ((cdata : RealArray.t),_,_) tt =
  let thispe = webdata.thispe in
  let npelast = webdata.npes - 1 in

  let clast = RealArray.create 2 in

  (* Send conc. at top right mesh point from PE npes-1 to PE 0. *)
  if thispe = npelast then
    let ilast = num_species*mxsub*mysub - 2 in
    if npelast <> 0
    then Mpi.send (slice cdata ilast 2) 0 0 comm
    else (clast.{0} <- cdata.{ilast};
          clast.{1} <- cdata.{ilast+1})

  (* On PE 0, receive conc. at top right from PE npes - 1.
     Then print performance data and sampled solution values. *)

  else if (thispe = 0) then begin

    if npelast <> 0 then begin
      let buf = (Mpi.receive npelast 0 comm : RealArray.t) in
      blit buf 0 clast 0 2
    end;

    let kused = Ida.get_last_order mem in
    let nst   = Ida.get_num_steps mem in
    let hused = Ida.get_last_step mem in

    printf "%8.2e %12.4e %12.4e   | %3d  %1d %12.4e\n"
           tt cdata.{0} clast.{0} nst kused hused;
    for i = 1 to num_species - 1 do
      printf "         %12.4e %12.4e   |\n" cdata.{i} clast.{i}
    done;
    printf "\n"
  end

(*
 * PrintFinalStats: Print final run data contained in iopt.
 *)

let print_final_stats mem =
  let open Ida in
  let nst  = get_num_steps mem in
  let nre  = get_num_res_evals mem in
  let netf = get_num_err_test_fails mem in
  let ncfn = get_num_nonlin_solv_conv_fails mem in
  let nni  = get_num_nonlin_solv_iters mem in

  let ncfl  = Spils.get_num_conv_fails mem in
  let nli   = Spils.get_num_lin_iters mem in
  let npe   = Spils.get_num_prec_evals mem in
  let nps   = Spils.get_num_prec_solves mem in
  let nreLS = Spils.get_num_res_evals mem in

  printf "-----------------------------------------------------------\n";
  printf "\nFinal statistics: \n\n";

  printf "Number of steps                    = %d\n" nst;
  printf "Number of residual evaluations     = %d\n" (nre+nreLS);
  printf "Number of nonlinear iterations     = %d\n" nni;
  printf "Number of error test failures      = %d\n" netf;
  printf "Number of nonlinear conv. failures = %d\n\n" ncfn;

  printf "Number of linear iterations        = %d\n" nli;
  printf "Number of linear conv. failures    = %d\n\n" ncfl;

  printf "Number of preconditioner setups    = %d\n" npe;
  printf "Number of preconditioner solves    = %d\n" nps

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA & SUPPORTING FUNCTIONS
 *--------------------------------------------------------------------
 *)
(*
 * BSend: Send boundary data to neighboring PEs.
 * This routine sends components of cc from internal subgrid boundaries
 * to the appropriate neighbor PEs.
 *)

let bsend comm my_pe isubx isuby dsizex dsizey udata =
  let bufleft = RealArray.create (num_species*mysub)
  and bufright = RealArray.create (num_species*mysub)
  in

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
      blit udata (ly*dsizex) bufleft (ly*num_species) num_species
    done;
    Mpi.send (slice bufleft 0 dsizey) (my_pe-1) 0 comm
  end;

  (* If isubx < NPEX-1, send data from right y-line of u (via bufright) *)
  if isubx <> npex-1 then begin
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*num_species in
      let offsetu = offsetbuf*mxsub + (mxsub-1)*num_species in
      blit udata offsetu bufright offsetbuf num_species
    done;
    Mpi.send (slice bufright 0 dsizey) (my_pe+1) 0 comm
  end


(*
 * BRecvPost: Start receiving boundary data from neighboring PEs.
 * (1) buffer should be able to hold 2*NUM_SPECIES*MYSUB realtype entries,
 *     should be passed to both the BRecvPost and BRecvWait functions, and
 *     should not be manipulated between the two calls.
 * (2) request should have 4 entries, and is also passed in both calls.
 *)

let brecvpost comm my_pe isubx isuby dsizex dsizey =
  (* If isuby > 0, receive data for bottom x-line of uext *)
  let r0 = if isuby <> 0
           then Mpi.ireceive (bytes dsizex) (my_pe-npex) 0 comm
           else Mpi.null_request
  in
  (* If isuby < NPEY-1, receive data for top x-line of uext *)
  let r1 = if isuby <> npey-1
           then Mpi.ireceive (bytes dsizex) (my_pe+npex) 0 comm
           else Mpi.null_request
  in
  (* If isubx > 0, receive data for left y-line of uext (via bufleft) *)
  let r2 = if isubx <> 0
           then Mpi.ireceive (bytes dsizey) (my_pe-1) 0 comm
           else Mpi.null_request
  in
  (* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) *)
  let r3 = if isubx <> npex-1
           then Mpi.ireceive (bytes dsizey) (my_pe+1) 0 comm
           else Mpi.null_request
  in
  [|r0; r1; r2; r3|]

(*
 * BRecvWait: Finish receiving boundary data from neighboring PEs.
 * (1) buffer should be able to hold 2*NUM_SPECIES*MYSUB realtype entries,
 *     should be passed to both the BRecvPost and BRecvWait functions, and
 *     should not be manipulated between the two calls.
 * (2) request should have 4 entries, and is also passed in both calls.
 *)

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

(*
 * WebRates: Evaluate reaction rates at a given spatial point.
 * At a given (x,y), evaluate the array of ns reaction terms R.
 *)

let web_rates webdata xx yy ((cxy : RealArray.t), cxy_off)
                            ((ratesxy : RealArray.t), ratesxy_off) =
  let acoef = RealArray2.unwrap webdata.acoef in

  for is = 0 to num_species-1 do
    ratesxy.{ratesxy_off + is} <- 0.;
    for j = 0 to num_species - 1 do
      ratesxy.{ratesxy_off + is} <- ratesxy.{ratesxy_off + is} +.
                      cxy.{cxy_off+j} *. acoef.{is, j}
    done
  done;

  let fac = 1. +. alpha*.xx*.yy +. beta*.sin(fourpi*.xx)*.sin(fourpi*.yy) in

  for is = 0 to num_species-1 do
    ratesxy.{ratesxy_off+is} <- cxy.{cxy_off+is}
                       *.( webdata.bcoef.{is}*.fac +. ratesxy.{ratesxy_off+is} )
  done

(*
 * reslocal: Compute res = F(t,cc,cp).
 * This routine assumes that all inter-processor communication of data
 * needed to calculate F has already been done.  Components at interior
 * subgrid boundaries are assumed to be in the work array cext.
 * The local portion of the cc vector is first copied into cext.
 * The exterior Neumann boundary conditions are explicitly handled here
 * by copying data from the first interior mesh line to the ghost cell
 * locations in cext.  Then the reaction and diffusion terms are
 * evaluated in terms of the cext array, and the residuals are formed.
 * The reaction terms are saved separately in the vector webdata.rates
 * for use by the preconditioner setup routine.
 *)

let reslocal webdata tt cc cp res =
  let mxsub =      (webdata.mxsub) in
  let mysub =      (webdata.mysub) in
  let npex =       (webdata.npex) in
  let npey =       (webdata.npey) in
  let ixsub =      (webdata.ixsub) in
  let jysub =      (webdata.jysub) in
  let nsmxsub =    (webdata.nsmxsub) in
  let nsmxsub2 =   (webdata.nsmxsub2) in
  let np =         (webdata.np) in
  let dx =         (webdata.dx) in
  let dy =         (webdata.dy) in
  let cox =        (webdata.cox) in
  let coy =        (webdata.coy) in
  let rhs =        (webdata.rhs) in
  let cext =       (webdata.cext) in
  let rates,_,_ =  (webdata.rates) in

  let (ccdata : RealArray.t),_,_ = cc in
  let (cpdata : RealArray.t),_,_ = cp in
  let (resdata : RealArray.t), _, _ = res in

  (* Get data pointers, subgrid data, array sizes, work array cext. *)

  (* Copy local segment of cc vector into the working extended array cext. *)
  let locc = ref 0
  and locce = ref (nsmxsub2 + num_species) in
  for jy = 0 to mysub-1 do
    for i = 0 to nsmxsub-1 do
      cext.{!locce+i} <- ccdata.{!locc+i}
    done;
    locc := !locc + nsmxsub;
    locce := !locce + nsmxsub2
  done;

  (* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of cc to cext. *)

  (* If jysub = 0, copy x-line 2 of cc to cext. *)
  if jysub = 0 then
    for i = 0 to nsmxsub-1 do
      cext.{num_species+i} <- ccdata.{nsmxsub+i}
    done
  ;

  (* If jysub = npey-1, copy x-line mysub-1 of cc to cext. *)
  if jysub = npey-1 then begin
    let locc = (mysub-2)*nsmxsub
    and locce = (mysub+1)*nsmxsub2 + num_species in
    for i = 0 to nsmxsub-1 do
      cext.{locce+i} <- ccdata.{locc+i}
    done
  end;

  (* If ixsub = 0, copy y-line 2 of cc to cext. *)
  if ixsub = 0 then begin
    for jy = 0 to mysub-1 do
      let locc = jy*nsmxsub + num_species
      and locce = (jy+1)*nsmxsub2 in
      for i = 0 to num_species - 1 do
        cext.{locce+i} <- ccdata.{locc+i}
      done
    done
  end;

  (* If ixsub = npex-1, copy y-line mxsub-1 of cc to cext. *)
  if ixsub = npex-1 then begin
    for jy = 0 to mysub - 1 do
      let locc  = (jy+1)*nsmxsub - 2*num_species
      and locce = (jy+2)*nsmxsub2 - num_species in
      for i = 0 to num_species-1 do
        cext.{locce+i} <- ccdata.{locc+i}
      done
    done
  end;

  (* Loop over all grid points, setting local array rates to right-hand sides.
     Then set res values appropriately for prey/predator components of F. *)
  for jy = 0 to mysub-1 do
    let ylocce = (jy+1)*nsmxsub2
    and yy     = float_of_int(jy+jysub*mysub)*.dy in

    for ix = 0 to mxsub-1 do
      let locce = ylocce + (ix+1)*num_species
      and xx = float_of_int(ix + ixsub*mxsub)*.dx in

      let off = ij_vptr_idx ix jy in
      web_rates webdata xx yy (cext, locce) (rates, off);

      for is = 0 to num_species-1 do
        let dcyli = cext.{locce+is}          -. cext.{locce+is-nsmxsub2} in
        let dcyui = cext.{locce+is+nsmxsub2} -. cext.{locce+is} in

        let dcxli = cext.{locce+is} -. cext.{locce+is-num_species} in
        let dcxui = cext.{locce+is+num_species} -. cext.{locce+is} in

        rhs.{is} <- cox.{is}*.(dcxui-.dcxli)
                    +. coy.{is}*.(dcyui-.dcyli)
                    +. rates.{off+is};

        if is < np
        then resdata.{off+is} <- cpdata.{off+is} -. rhs.{is}
        else resdata.{off+is} <-                 -. rhs.{is}
      done (* End of is (species) loop. *)
    done (* End of ix loop. *)
  done (* End of jy loop. *)

(*
 * rescomm: Communication routine in support of resweb.
 * This routine performs all inter-processor communication of components
 * of the cc vector needed to calculate F, namely the components at all
 * interior subgrid boundaries (ghost cell data).  It loads this data
 * into a work array cext (the local portion of c, extended).
 * The message-passing uses blocking sends, non-blocking receives,
 * and receive-waiting, in routines BRecvPost, BSend, BRecvWait.
 *)

let rescomm webdata cc cp =
  (* Get comm, thispe, subgrid indices, data sizes, extended array cext. *)
  let comm = webdata.comm   and thispe = webdata.thispe
  and ixsub = webdata.ixsub and jysub = webdata.jysub
  and cext = webdata.cext
  and nsmxsub = webdata.nsmxsub and nsmysub = (webdata.ns)*(webdata.mysub) in
  let (cdata : RealArray.t), _, _ = cc in

  (* Start receiving boundary data from neighboring PEs. *)
  let rs = brecvpost comm thispe ixsub jysub nsmxsub nsmysub in

  (* Send data from boundary of local grid to neighboring PEs. *)
  bsend comm thispe ixsub jysub nsmxsub nsmysub cdata;

  (* Finish receiving boundary data from neighboring PEs. *)
  brecvwait rs ixsub jysub nsmxsub cext

(*
 * resweb: System residual function for predator-prey system.
 * To compute the residual function F, this routine calls:
 *    rescomm, for needed communication, and then
 *    reslocal, for computation of the residuals on this processor.
 *)

let resweb webdata tt cc cp res =
  (* Call rescomm to do inter-processor communication. *)
  rescomm webdata cc cp;
  (* Call reslocal to calculate the local portion of residual vector. *)
  reslocal webdata tt cc cp res

(*
 * SetInitialProfiles: Set initial conditions in cc, cp, and id.
 * A polynomial profile is used for the prey cc values, and a constant
 * (1.0e5) is loaded as the initial guess for the predator cc values.
 * The id values are set to 1 for the prey and 0 for the predators.
 * The prey cp values are set according to the given system, and
 * the predator cp values are set to zero.
 *)

let set_initial_profiles webdata cc cp id res =
  let ixsub = webdata.ixsub in
  let jysub = webdata.jysub in
  let mxsub = webdata.mxsub in
  let mysub = webdata.mxsub in
  let dx = webdata.dx in
  let dy = webdata.dy in
  let np = webdata.np in

  let (ccdata : RealArray.t), _, _ = cc in
  let (cpdata : RealArray.t), _, _ = cp in
  let (iddata : RealArray.t), _, _ = id in

  (* Loop over grid, load cc values and id values. *)
  for jy = 0 to mysub-1 do
    let yy = float_of_int (jy + jysub*mysub) *. dy in
    for ix = 0 to mxsub-1 do
      let xx = float_of_int (ix + ixsub*mxsub) *. dx in
      let xyfactor = 16.0*.xx*.(1.0 -. xx)*.yy*.(1.0 -. yy) in
      let xyfactor = xyfactor *. xyfactor in

      let off = ij_vptr_idx ix jy in
      for is = 0 to num_species-1 do
        if is < np
        then (ccdata.{off+is} <- 10.0 +. float_of_int(is+1)*.xyfactor;
              iddata.{off+is} <- 1.0)
        else (ccdata.{off+is} <- 1.0e5;
              iddata.{off+is} <- 0.0)
      done
    done
  done;

  (* Set c' for the prey by calling the residual function with cp = 0. *)
  vconst 0.0 cp;
  resweb webdata 0.0 cc cp res;
  vscale (-.1.0) res cp;

  (* Set c' for predators to 0. *)
  for jy = 0 to mysub-1 do
    for ix = 0 to mxsub-1 do
      let off = ij_vptr_idx ix jy in
      for is = np to num_species-1 do
        cpdata.{off+is} <- 0.0
      done
    done
  done

(*
 * Preconbd: Preconditioner setup routine.
 * This routine generates and preprocesses the block-diagonal
 * preconditoner pp.  At each spatial point, a block of pp is computed
 * by way of difference quotients on the reaction rates R.
 * The base value of R are taken from webdata.rates, as set by webres.
 * Each block is LU-factored, for later solution of the linear systems.
 *)

let precondbd webdata jac =
  let cc,_,_ = jac.Ida.jac_y
  and cp,_,_ = jac.Ida.jac_y'
  and cj = jac.Ida.jac_coef
  in

  let mxsub =      webdata.mxsub in
  let mysub =      webdata.mysub in
  let ixsub =      webdata.ixsub in
  let jysub =      webdata.jysub in
  let np =         webdata.np in
  let dx =         webdata.dx in
  let dy =         webdata.dy in
  let rates,_,_ =  webdata.rates in
  let ns =         webdata.ns in

  let perturb_rates = RealArray.create num_species in

  let uround = Sundials.unit_roundoff in
  let sqru = sqrt uround in
  let mem =
    match webdata.ida_mem with
    | Some m -> m
    | None -> invalid_arg "Internal error: webdata.ida_mem not set"
  in
  let ewt = Nvector_parallel.local_array webdata.ewt in
  Ida.get_err_weights mem webdata.ewt;
  let hh = Ida.get_current_step mem in

  for jy = 0 to mysub-1 do
    let yy = float_of_int(jy + jysub*mysub)*.dy in

    for ix = 0 to mxsub-1 do
      let xx = float_of_int(ix+ ixsub*mxsub)*.dx in
      let pxy = webdata.pp.(ix).(jy) in
      let off = ij_vptr_idx ix jy in

      for js = 0 to ns-1 do
        let inc = sqru*.(max (abs_float (cc.{off+js}))
                             (max (hh*.abs_float (cp.{off+js}))
                                  (1.0/.ewt.{off+js}))) in
        let cctemp = cc.{off+js} in(* Save the (js,ix,jy) element of cc. *)
        cc.{off+js} <- cc.{off+js} +. inc;    (* Perturb the (js,ix,jy) element of cc. *)
        let fac = -.1.0/.inc in

        web_rates webdata xx yy (cc, off) (perturb_rates, 0);

        let pxycol = RealArray2.col pxy js in

        for is = 0 to ns-1 do
          pxycol.{is} <- (perturb_rates.{is} -. rates.{off+is})*.fac
        done;

        (* Add partial with respect to cp. *)
        if js < np then pxycol.{js} <- pxycol.{js} +. cj;

        cc.{off+js} <- cctemp; (* Restore (js,ix,jy) element of cc. *)

      done; (* End of js loop. *)

      (* Do LU decomposition of matrix block for grid point (ix,jy). *)
      (try Matrix.ArrayDense.getrf pxy webdata.pivot.(ix).(jy)
       with Matrix.ZeroDiagonalElement _ -> raise Sundials.RecoverableFailure)

    done (* End of ix loop. *)
  done (* End of jy loop. *)

(*
 * PSolvebd: Preconditioner solve routine.
 * This routine applies the LU factorization of the blocks of the
 * preconditioner PP, to compute the solution of PP * zvec = rvec.
 *)

let psolvebd webdata jac rvec zvec delta =
  vscale 1.0 rvec zvec;

  (* Loop through subgrid and apply preconditioner factors at each point. *)
  for ix = 0 to mxsub-1 do
    for jy = 0 to mysub-1 do

      (* For grid point (ix,jy), do backsolve on local vector.
         zxy is the address of the local portion of zvec, and
         Pxy is the address of the corresponding block of PP.  *)
      let zxy = ij_vptr zvec ix jy in
      let pxy = webdata.pp.(ix).(jy) in
      let pivot = webdata.pivot.(ix).(jy) in
      Matrix.ArrayDense.getrs pxy pivot zxy;
    done (* End of jy loop. *)
  done (* End of ix loop. *)

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* Set communicator, and get processor number and total number of PE's. *)

  let comm = Mpi.comm_world in
  let thispe = Mpi.comm_rank comm in
  let npes = Mpi.comm_size comm in

  if npes <> npex*npey then begin
    if thispe = 0 then
      fprintf stderr
        "\nMPI_ERROR(0): npes = %d not equal to NPEX*NPEY = %d\n"
        npes (npex*npey)
    ;
    exit 1
  end;

  (* Set local length (local_N) and global length (SystemSize). *)

  let local_N = mxsub*mysub*num_species in
  let system_size = neq in

  (* Set up user data block webdata. *)

  let webdata = alloc_init_user_data comm local_N system_size thispe npes in

  (* Create needed vectors, and load initial values.
     The vector res is used temporarily only.        *)

  let cc = Nvector_parallel.make local_N system_size comm 0. in
  let cp = Nvector_parallel.make local_N system_size comm 0. in
  let res = Nvector_parallel.make local_N system_size comm 0. in
  let id = Nvector_parallel.make local_N system_size comm 0. in

  set_initial_profiles webdata
    (Nvector.unwrap cc) (Nvector.unwrap cp)
    (Nvector.unwrap id) (Nvector.unwrap res);

  (* Set remaining inputs to IDAMalloc. *)

  let t0 = 0.0 in

  (* Call IDACreate and IDAMalloc to initialize IDA.
     A pointer to IDA problem memory is returned and stored in idamem. *)

  (* Call IDASpgmr to specify the IDA linear solver IDASPGMR and specify
     the preconditioner routines supplied (Precondbd and PSolvebd).
     maxl (max. Krylov subspace dim.) is set to 16. *)
  let maxl = 16 in
  let lsolver = Ida.Spils.(spgmr ~maxl cc) in
  let mem =
    Ida.(init
      Spils.(solver lsolver
                    (prec_left ~setup:(precondbd webdata) (psolvebd webdata)))
      (SStolerances (rtol, atol))
      (resweb webdata)
      t0 cc cp)
  in
  (match Sundials.sundials_version with
   | 2,_,_ -> () | _ -> Ida.Spils.set_max_restarts lsolver 5);
  webdata.ida_mem <- Some mem;

  (* Call IDACalcIC (with default options) to correct the initial values. *)

  let tout = ref 0.001 in
  Ida.calc_ic_ya_yd' mem ~varid:id !tout;

  (* On PE 0, print heading, basic parameters, initial values. *)

  if thispe = 0 then
    print_header system_size maxl rtol atol
  ;
  print_output webdata comm mem (Nvector.unwrap cc) t0;

  (* Loop over iout, call IDASolve (normal mode), print selected output. *)

  for iout = 1 to nout do

    let (tret, _) = Ida.solve_normal mem !tout cc cp in

    print_output webdata comm mem (Nvector.unwrap cc) tret;

    if iout < 3 then tout := !tout *. tmult
    else             tout := !tout +. tadd;

  done;

  (* On PE 0, print final set of statistics. *)
  if thispe = 0 then print_final_stats mem

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
