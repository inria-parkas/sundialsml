(* * -----------------------------------------------------------------
 * $Revision:
 * $Date:
 * -----------------------------------------------------------------
 * Programmer(s): Cosmin Petra and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Aug 2014.
 * -----------------------------------------------------------------
 * Example program for IDA: Brusselator, parallel, GMRES, IDABBD
 * preconditioner, FSA.
 *
 * This example program for IDAS uses IDASPGMR as the linear solver.
 * It is written for a parallel computer system and uses the
 * IDABBDPRE band-block-diagonal preconditioner module for the
 * IDASPGMR package.
 *
 * The mathematical problem solved in this example is a DAE system
 * that arises from a system of partial differential equations after
 * spatial discretization.
 *
 * The PDE system is a two-species time-dependent PDE known as
 * Brusselator PDE and models a chemically reacting system.
 *
 *
 *  du/dt = eps1(u  + u  ) + u^2 v -(B+1)u + A
 *                xx   yy
 *                                              domain [0,L]X[0,L]
 *  dv/dt = eps2(v  + v  ) - u^2 v + Bu
 *                xx   yy
 *
 *  B.C. Neumann
 *  I.C  u(x,y,t0) = u0(x,y) =  1  - 0.5*cos(pi*y/L)
 *       v(x,y,t0) = v0(x,y) = 3.5 - 2.5*cos(pi*x/L)
 *
 * The PDEs are discretized by central differencing on a MX by MY
 * mesh, and so the system size Neq is the product MX*MY*NUM_SPECIES.
 * The system is actually implemented on submeshes, processor by
 * processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
 * processors.
 *
 * The average of the solution u at final time is also computed.
 *            / /
 *        g = | | u(x,y,tf) dx dy
 *            / /
 * Also the sensitivities of g with respect to parameters eps1 and
 * eps2 are computed.
 *                  / /
 *       dg/d eps = | | u  (x,y,tf)  dx dy
 *                  / /  eps
 *)

module RealArray = Sundials.RealArray
module RealArray2 = Sundials.RealArray2
module LintArray = Sundials.LintArray
let unvec = Sundials.unvec
module Sens = Idas.Sensitivity
open Nvector_parallel.DataOps

let fprintf = Printf.fprintf
let printf = Printf.printf

let slice = Bigarray.Array1.sub

let blit (buf : RealArray.t) buf_offset (dst : RealArray.t) dst_offset len =
  for i = 0 to len-1 do
    dst.{dst_offset + i} <- buf.{buf_offset + i}
  done

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_string (RealArray.empty) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_string (RealArray.create 1) []) 0
  - header_and_empty_array_size

let bytes x = header_and_empty_array_size + x * float_cell_size

(* Problem Constants *)
let num_species = 2
let ctL =         1.0    (* Domain =[0,L]^2 *)
let ctA =         1.0
let ctB =         3.4
let ctEps =       2.0e-3
let num_sens =    2

let pi =          3.1415926535898 (* pi *)

let mxsub =       41    (* Number of x mesh points per processor subgrid *)
let mysub =       41    (* Number of y mesh points per processor subgrid *)
let npex =        2     (* Number of subgrids in the x direction *)
let npey =        2     (* Number of subgrids in the y direction *)
let mx =          (mxsub*npex)      (* MX = number of x mesh points *)
let my =          (mysub*npey)      (* MY = number of y mesh points *)
let nsmxsub =     (num_species * mxsub)
let neq =         (num_species*mx*my) (* Number of equations in system *)


let rtol =        1.e-5  (*  rtol tolerance *)
let atol =        1.e-5  (*  atol tolerance *)
let nout =        6
let tmult =       10.0   (* Multiplier for tout values *)
let tadd =        0.3    (* Increment for tout values *)

let zero =        0.0
let half =        0.5
let one =         1.0


(* User-defined vector accessor macro IJ_Vptr. *)

(*
 * IJ_Vptr is defined in order to express the underlying 3-d structure of the
 * dependent variable vector from its underlying 1-d storage (an N_Vector).
 * IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to
 * species index is = 0, x-index ix = i, and y-index jy = j.
 *)
let index i j = i*num_species + j*nsmxsub
let ij_vptr (local,_,_) i j =
  let offset = index i j in
  slice local offset (RealArray.length local - offset)

type user_data =
  {
    ns : int;
    thispe : int;
    npes : int;
    ixsub : int;
    jysub : int;
    npex : int;
    npey : int;
    mxsub : int;
    mysub : int;
    nsmxsub : int;
    nsmxsub2 : int;
    a : float;
    b : float;
    l : float;
    eps : RealArray.t;                  (* size = num_species *)
    dx : float;
    dy : float;
    cox : RealArray.t;                  (* size = num_species *)
    coy : RealArray.t;                  (* size = num_species *)
    gridext : RealArray.t;              (* size = (mxsub+2)*(mysub+2)*num_species *)
    rhs : RealArray.t;                  (* size = num_species *)
    comm : Mpi.communicator;
    rates : RealArray.t;                (* size = 2 *)
    n_local : int;
  }

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA & SUPPORTING FUNCTIONS
 *--------------------------------------------------------------------
 *)

(*
 * BRecvPost: Start receiving boundary data from neighboring PEs.
 *)

let brecvpost comm my_pe ixsub jysub dsizex dsizey cext =
  (* If jysub > 0, receive data for bottom x-line of cext. *)
  let r0 = if jysub <> 0
           then Mpi.ireceive (bytes dsizex) (my_pe-npex) 0 comm
           else Mpi.null_request
  in

  (* If jysub < npey-1, receive data for top x-line of cext. *)
  let r1 = if jysub <> npey-1
           then Mpi.ireceive (bytes dsizex) (my_pe+npex) 0 comm
           else Mpi.null_request
  in

  (* If ixsub > 0, receive data for left y-line of cext (via bufleft). *)
  let r2 = if ixsub <> 0
           then Mpi.ireceive (bytes dsizey) (my_pe-1) 0 comm
           else Mpi.null_request
  in

  (* If ixsub < npex-1, receive data for right y-line of cext (via bufright). *)
  let r3 = if ixsub <> npex-1
           then Mpi.ireceive (bytes dsizey) (my_pe+1) 0 comm
           else Mpi.null_request
  in
  [|r0;r1;r2;r3|]

(*
 * BRecvWait: Finish receiving boundary data from neighboring PEs.
 * (1) buffer should be able to hold 2*num_species*mysub realtype entries,
 *     should be passed to both the brecvpost and BRecvWait functions, and
 *     should not be manipulated between the two calls.
 * (2) request should have 4 entries, and is also passed in both calls.
 *)

let brecvwait request ixsub jysub dsizex cext =
  let dsizex2 = dsizex + 2*num_species in

  (* If jysub > 0, receive data for bottom x-line of cext. *)
  if jysub <> 0 then begin
    let buf = (Mpi.wait_receive request.(0) : RealArray.t) in
    blit buf 0 cext num_species dsizex
  end;

  (* If jysub < npey-1, receive data for top x-line of cext. *)
  if jysub <> npey-1 then begin
    let buf = (Mpi.wait_receive request.(1) : RealArray.t) in
    let offsetce = num_species*(1 + (mysub+1)*(mxsub+2)) in
    blit buf 0 cext offsetce dsizex
  end;

  (* If ixsub > 0, receive data for left y-line of cext (via bufleft). *)
  if ixsub <> 0 then begin
    let bufleft = (Mpi.wait_receive request.(2) : RealArray.t) in

    (* Copy the buffer to cext *)
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*num_species in
      let offsetce = (ly+1)*dsizex2 in
      for i = 0 to num_species-1 do
        cext.{offsetce+i} <- bufleft.{offsetbuf+i}
      done
    done
  end;

  (* If ixsub < npex-1, receive data for right y-line of cext (via bufright). *)
  if ixsub <> npex-1 then begin
    let bufright = (Mpi.wait_receive request.(3) : RealArray.t) in

    (* Copy the buffer to cext *)
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*num_species in
      let offsetce = (ly+2)*dsizex2 - num_species in
      for i = 0 to num_species-1 do
        cext.{offsetce+i} <- bufright.{offsetbuf+i}
      done
    done
  end

(*
 * BSend: Send boundary data to neighboring PEs.
 * This routine sends components of uv from internal subgrid boundaries
 * to the appropriate neighbor PEs.
 *)
let bsend comm my_pe ixsub jysub dsizex dsizey (cdata : RealArray.t) =
  let bufleft = RealArray.create (num_species * mysub)
  and bufright = RealArray.create (num_species * mysub)
  in
  (* If jysub > 0, send data from bottom x-line of uv. *)

  if jysub <> 0 then
    Mpi.send (slice cdata 0 dsizex) (my_pe-npex) 0 comm
  ;

  (* If jysub < npey-1, send data from top x-line of uv. *)

  if jysub <> npey-1 then begin
    let offsetc = (mysub-1)*dsizex in
    Mpi.send (slice cdata offsetc dsizex) (my_pe+npex) 0 comm
  end;

  (* If ixsub > 0, send data from left y-line of uv (via bufleft). *)

  if ixsub <> 0 then begin
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*num_species in
      let offsetc = ly*dsizex in
      for i = 0 to num_species-1 do
        bufleft.{offsetbuf+i} <- cdata.{offsetc+i}
      done
    done;
    Mpi.send (slice bufleft 0 dsizey) (my_pe-1) 0 comm
  end;

  (* If ixsub < npex-1, send data from right y-line of uv (via bufright). *)

  if ixsub <> npex-1 then begin
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*num_species in
      let offsetc = offsetbuf*mxsub + (mxsub-1)*num_species in
      for i = 0 to num_species-1 do
        bufright.{offsetbuf+i} <- cdata.{offsetc+i}
      done
    done;
    Mpi.send (slice bufright 0 dsizey) (my_pe+1) 0 comm
  end

(*
 * ReactRates: Evaluate reaction rates at a given spatial point.
 * At a given (x,y), evaluate the array of ns reaction terms R.
 *)

let react_rates data xx yy ((uvval : RealArray.t), uvval_off)
                           (rates : RealArray.t) =
  let a = data.a and b = data.b in

  rates.{0} <- uvval.{uvval_off}*.uvval.{uvval_off}*.uvval.{uvval_off + 1};
  rates.{1} <- -. rates.{0};

  rates.{0} <- rates.{0} +. (a-.(b+.1.0)*.uvval.{uvval_off});
  rates.{1} <- rates.{1} +. b*.uvval.{uvval_off}

(*
 * reslocal: Compute res = F(t,uv,uvp).
 * This routine assumes that all inter-processor communication of data
 * needed to calculate F has already been done.  Components at interior
 * subgrid boundaries are assumed to be in the work array cext.
 * The local portion of the uv vector is first copied into cext.
 * The exterior Neumann boundary conditions are explicitly handled here
 * by copying data from the first interior mesh line to the ghost cell
 * locations in cext.  Then the reaction and diffusion terms are
 * evaluated in terms of the cext array, and the residuals are formed.
 * The reaction terms are saved separately in the vector data.rates
 * for use by the preconditioner setup routine.
 *)

let reslocal data tt ((uv : RealArray.t), _, _)
                     ((uvp : RealArray.t), _, _)
                     ((rr : RealArray.t), _, _) =
  let mxsub =      data.mxsub in
  let mysub =      data.mysub in
  let npex =       data.npex in
  let npey =       data.npey in
  let ixsub =      data.ixsub in
  let jysub =      data.jysub in
  let nsmxsub =    data.nsmxsub in
  let nsmxsub2 =   data.nsmxsub2 in
  let dx =         data.dx in
  let dy =         data.dy in
  let gridext =    data.gridext in
  let eps =        data.eps in
  (* Get data pointers, subgrid data, array sizes, work array cext. *)
  let rates = RealArray.create 2 in

  let dx2 = dx *. dx in
  let dy2 = dy *. dy in

  (* Copy local segment of uv vector into the working extended array gridext. *)
  let locc = ref 0 in
  let locce = ref (nsmxsub2 + num_species) in
  for jy = 0 to mysub-1 do
    for i = 0 to nsmxsub-1 do
      gridext.{!locce+i} <- uv.{!locc+i}
    done;
    locc := !locc + nsmxsub;
    locce := !locce + nsmxsub2;
  done;

  (* To facilitate homogeneous Neumann boundary conditions, when this is
     a boundary PE, copy data from the first interior mesh line of uv to gridext. *)

  (* If jysub = 0, copy x-line 2 of uv to gridext. *)
  if jysub = 0 then
    for i = 0 to nsmxsub-1 do
      gridext.{num_species+i} <- uv.{nsmxsub+i}
    done
  ;

  (* If jysub = npey-1, copy x-line mysub-1 of uv to gridext. *)
  if jysub = npey-1 then begin
    let locc = (mysub-2)*nsmxsub in
    let locce = (mysub+1)*nsmxsub2 + num_species in
    for i = 0 to nsmxsub-1 do
      gridext.{locce+i} <- uv.{locc+i}
    done
  end;


  (* If ixsub = 0, copy y-line 2 of uv to gridext. *)
  if ixsub = 0 then begin
    for jy = 0 to mysub-1 do
      let locc = jy*nsmxsub + num_species in
      let locce = (jy+1)*nsmxsub2 in
      for i = 0 to num_species-1 do
        gridext.{locce+i} <- uv.{locc+i}
      done
    done
  end;


  (* If ixsub = npex-1, copy y-line mxsub-1 of uv to gridext. *)
  if ixsub = npex-1 then begin
    for jy = 0 to mysub-1 do
      let locc = (jy+1)*nsmxsub - 2*num_species in
      let locce = (jy+2)*nsmxsub2 - num_species in
      for i = 0 to num_species-1 do
        gridext.{locce+i} <- uv.{locc+i}
      done
    done
  end;

  (* Loop over all grid points, setting local array rates to right-hand sides.
     Then set rr values appropriately (ODE in the interior and DAE on the boundary)*)
  let ixend = if ixsub=npex-1 then 1 else 0 in
  let ixstart = if ixsub=0 then 1 else 0 in
  let jystart = if jysub=0 then 1 else 0 in
  let jyend = if jysub=npey-1 then 1 else 0 in

  for jy = jystart to mysub-jyend-1 do
    let ylocce = (jy+1)*nsmxsub2 in
    let yy = float_of_int(jy+jysub*mysub)*.dy in

    for ix = ixstart to mxsub-ixend-1 do
      let locce = ylocce + (ix+1)*num_species in
      let xx = float_of_int(ix + ixsub*mxsub)*.dx in

      react_rates data xx yy (gridext, locce) rates;

      let off = index ix jy in

      for is = 0 to num_species-1 do
        let dcyli = gridext.{locce+is}          -. gridext.{locce+is-nsmxsub2} in
        let dcyui = gridext.{locce+is+nsmxsub2} -. gridext.{locce+is} in

        let dcxli = gridext.{locce+is}             -. gridext.{locce+is-num_species} in
        let dcxui = gridext.{locce+is+num_species} -. gridext.{locce+is} in

        rr.{off + is} <- uvp.{off + is}
                      -. eps.{is}*.( (dcxui-.dcxli)/.dx2 +. (dcyui-.dcyli)/.dy2 )
                      -. rates.{is};
      done
    done
  done;

  (* Algebraic equation correspoding to boundary mesh point. *)
  if jysub=0 then begin
    for ix = 0 to mxsub-1 do
      let locce = nsmxsub2 + num_species * (ix+1) in
      let rr_off = index ix 0 in

      for is = 0 to num_species-1 do
        rr.{rr_off + is} <- gridext.{locce+is+nsmxsub2} -. gridext.{locce+is}
      done
    done
  end;

  if ixsub=npex-1 then begin
    for jy = 0 to mysub-1 do
      let locce = (jy+1)*nsmxsub2 + nsmxsub2-num_species in
      let rr_off = index (mxsub-1) jy in

      for is = 0 to num_species-1 do
        rr.{rr_off + is} <- gridext.{locce+is-num_species} -. gridext.{locce+is}
      done
    done
  end;

  if ixsub=0 then begin
    for jy = 0 to mysub-1 do
      let locce = (jy+1)*nsmxsub2 + num_species in
      let rr_off = index 0 jy in

      for is = 0 to num_species-1 do
        rr.{rr_off + is} <- gridext.{locce+is-num_species} -. gridext.{locce+is}
      done
    done
  end;

  if jysub=npey-1 then begin
    for ix = 0 to mxsub-1 do
      let locce = nsmxsub2*mysub + (ix+1)*num_species in
      let rr_off = index ix (mysub - 1) in

      for is = 0 to num_species-1 do
        rr.{rr_off + is} <- gridext.{locce+is-nsmxsub2} -. gridext.{locce+is}
      done
    done
  end

(*
 * rescomm: Communication routine in support of resweb.
 * This routine performs all inter-processor communication of components
 * of the uv vector needed to calculate F, namely the components at all
 * interior subgrid boundaries (ghost cell data).  It loads this data
 * into a work array cext (the local portion of c, extended).
 * The message-passing uses blocking sends, non-blocking receives,
 * and receive-waiting, in routines BRecvPost, BSend, BRecvWait.
 *)
let rescomm data tt uv uvp =
  let cdata,_,_ = uv in

  (* Get comm, thispe, subgrid indices, data sizes, extended array cext. *)
  let comm = data.comm in
  let thispe = data.thispe in

  let ixsub = data.ixsub in
  let jysub = data.jysub in
  let gridext = data.gridext in
  let nsmxsub = data.nsmxsub in
  let nsmysub = (data.ns)*(data.mysub) in

  (* Start receiving boundary data from neighboring PEs. *)
  let request = brecvpost comm thispe ixsub jysub nsmxsub nsmysub gridext
  in

  (* Send data from boundary of local grid to neighboring PEs. *)
  bsend comm thispe ixsub jysub nsmxsub nsmysub cdata;

  (* Finish receiving boundary data from neighboring PEs. *)
  brecvwait request ixsub jysub nsmxsub gridext

(*
 * res: System residual function
 *
 * To compute the residual function F, this routine calls:
 * rescomm, for needed communication, and then
 * reslocal, for computation of the residuals on this processor.
 *)

let res data tt uv uvp rr =
  (* Call rescomm to do inter-processor communication. *)
  rescomm data tt uv uvp;

  (* Call reslocal to calculate the local portion of residual vector. *)
  reslocal data tt uv uvp rr


(* Integrate over the spatial domain. Each process computes the integral on its
   grid. Then processes call MPI_REDUCE to compute sum of the local values. *)
let integr comm uv data =
  (* compute the integral on the (local) grid *)
  let uvdata = Nvector_parallel.unwrap uv in
  let buf = RealArray.create 2 in

  let intval = ref 0.0 in

  for jy = 1 to mysub-1 do
    for ix = 1 to mxsub-1 do
      (* consider only u *)
      intval := !intval +. uvdata.{ix*num_species + jy*nsmxsub};
    done
  done;
  intval := !intval *. (data.dx *. data.dy);

  buf.{0} <- !intval;

  (* Sum local values and get the result on all processors. *)
  buf.{1} <- Mpi.allreduce_float buf.{0} Mpi.Float_sum comm;

  intval := buf.{1};
  !intval

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(*
 * InitUserData: Load problem constants in data (of type UserData).
 *)
let init_user_data thispe npes comm =
  let jysub = thispe / npex in
  let ixsub = thispe - (jysub)*npex in
  let ns = num_species in
  let dx = ctL/.float_of_int(mx-1) in
  let dy = ctL/.float_of_int(my-1) in
  let thispe = thispe in
  let npes = npes in
  let nsmxsub = mxsub * num_species in
  let nsmxsub2 = (mxsub+2)*num_species in
  let n_local = mxsub*mysub*num_species in
  let a = ctA in
  let b = ctB in
  let l = ctL in
  let eps = RealArray.of_array [|ctEps; ctEps|] in
  {
    ns = ns;
    thispe = thispe;
    npes = npes;
    ixsub = ixsub;
    jysub = jysub;
    npex = npex;
    npey = npey;
    mxsub = mxsub;
    mysub = mysub;
    nsmxsub = nsmxsub;
    nsmxsub2 = nsmxsub2;
    a = a;
    b = b;
    l = l;
    eps = eps;
    dx = dx;
    dy = dy;
    cox = RealArray.create num_species;
    coy = RealArray.create num_species;
    gridext = RealArray.create ((mxsub+2)*(mysub+2)*num_species);
    rhs = RealArray.create num_species;
    comm = comm;
    rates = RealArray.create 2;
    n_local = n_local;
  }


(*
 * SetInitialProfiles: Set initial conditions in uv, uvp, and id.
 *)

let set_initial_profiles data uv uvp id resid =
  let ixsub = data.ixsub in
  let jysub = data.jysub in
  let mxsub = data.mxsub in
  let mysub = data.mysub in
  let npex = data.npex in
  let npey = data.npey in
  let dx = data.dx in
  let dy = data.dy in
  let l = data.l in

  n_vconst 0. uv;

  (* Loop over grid, load uv values and id values. *)
  for jy = 0 to mysub-1 do
    let y = float_of_int(jy + jysub*mysub) *. dy in
    for ix = 0 to mxsub-1 do

      let x = float_of_int(ix + ixsub*mxsub) *. dx in
      let uvxy = ij_vptr uv ix jy in

      uvxy.{0} <- 1.0 -. half*.cos(pi*.y/.l);
      uvxy.{1} <- 3.5 -. 2.5*.cos(pi*.x/.l);
    done
  done;

  n_vconst one id;

  if jysub = 0 then begin
    for ix = 0 to mxsub-1 do
      let idxy = ij_vptr id ix 0 in
      idxy.{0} <- zero;
      idxy.{1} <- zero;

      let uvxy = ij_vptr uv ix 0 in
      let uvxy1 = ij_vptr uv ix 1 in
      uvxy.{0} <- uvxy1.{0};
      uvxy.{1} <- uvxy1.{1};
    done
  end;

  if ixsub = npex-1 then begin
    for jy = 0 to mysub-1 do
      let idxy = ij_vptr id (mxsub-1) jy in
      idxy.{0} <- zero;
      idxy.{1} <- zero;

      let uvxy = ij_vptr uv (mxsub-1) jy in
      let uvxy1 = ij_vptr uv (mxsub-2) jy in
      uvxy.{0} <- uvxy1.{0};
      uvxy.{1} <- uvxy1.{1};

    done
  end;

  if ixsub = 0 then begin
    for jy = 0 to mysub-1 do
      let idxy = ij_vptr id 0 jy in
      idxy.{0} <- zero;
      idxy.{1} <- zero;

      let uvxy = ij_vptr uv 0 jy in
      let uvxy1 = ij_vptr uv 1 jy in
      uvxy.{0} <- uvxy1.{0};
      uvxy.{1} <- uvxy1.{1}
    done
  end;

  if jysub = npey-1 then begin
    for ix = 0 to mxsub-1 do
      let idxy = ij_vptr id ix jysub in
      idxy.{0} <- zero;
      idxy.{1} <- zero;

      let uvxy = ij_vptr uv ix (mysub-1) in
      let uvxy1 = ij_vptr uv ix (mysub-2) in
      uvxy.{0} <- uvxy1.{0};
      uvxy.{1} <- uvxy1.{1}
    done
  end;

  (* Derivative found by calling the residual function with uvp = 0. *)
  n_vconst zero uvp;
  res data zero uv uvp resid;
  n_vscale (-.one) resid uvp

(*
 * Print first lines of output (problem description)
 * and table headerr
 *)

let print_header system_size maxl mudq mldq mukeep mlkeep rtol atol =
  printf "\n Brusselator PDE -  DAE parallel example problem for IDA \n\n";
  printf "Number of species ns: %d" num_species;
  printf "     Mesh dimensions: %d x %d\n" mx my;
  printf "Total system size: %d\n" system_size;
  printf "Subgrid dimensions: %d x %d" mxsub mysub;
  printf "     Processor array: %d x %d\n" npex npey;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Linear solver: IDASPGMR     Max. Krylov dimension maxl: %d\n" maxl;
  printf "Preconditioner: band-block-diagonal (IDABBDPRE), with parameters\n";
  printf "     mudq = %d,  mldq = %d,  mukeep = %d,  mlkeep = %d\n"
         mudq mldq mukeep mlkeep;
  printf "CalcIC called to correct initial concentrations \n\n";
  printf "-----------------------------------------------------------\n";
  printf "  t        bottom-left  top-right";
  printf "    | nst  k      h\n";
  printf "-----------------------------------------------------------\n\n"

(*
 * PrintOutput: Print output values at output time t = tt.
 * Selected run statistics are printed.  Then values of c1 and c2
 * are printed for the bottom left and top right grid points only.
 *)

let print_output mem uv tt data comm =
  let thispe = data.thispe in
  let npelast = data.npes - 1 in
  let cdata = Nvector_parallel.unwrap uv in
  let clast = RealArray.create 2 in

  (* Send conc. at top right mesh point from PE npes-1 to PE 0. *)
  if thispe = npelast then begin
    let ilast = num_species*mxsub*mysub - 2 in
    if npelast <> 0
    then Mpi.send (slice cdata ilast 2) 0 0 comm
    else (clast.{0} <- cdata.{ilast}; clast.{1} <- cdata.{ilast+1})
  end;

  (* On PE 0, receive conc. at top right from PE npes - 1.
     Then print performance data and sampled solution values. *)
  if thispe = 0 then begin

    if npelast <> 0 then begin
      let buf = (Mpi.receive npelast 0 comm : RealArray.t) in
      blit buf 0 clast 0 2
    end;

    let kused = Ida.get_last_order mem in
    let nst = Ida.get_num_steps mem in
    let hused = Ida.get_last_step mem in

    printf "%8.2e %12.4e %12.4e   | %3d  %1d %12.4e\n"
         tt cdata.{0} clast.{0} nst kused hused;
    for i = 1 to num_species-1 do
      printf "         %12.4e %12.4e   |\n" cdata.{i} clast.{i}
    done;
    printf "\n"
  end

(*
 * PrintSol the PE's portion of the solution to a file.
 *)
let print_sol data mem uv uvp comm =
  let thispe = data.thispe in
  let szFilename = Printf.sprintf "ysol%d.txt" thispe in

  (* NB: the original C code opens with "w+" instead of "w", but seems
     to only write to the file. *)
  let fout =
    try open_out szFilename
    with e ->
      printf "PE[% 2d] is unable to write solution to disk!\n" thispe;
      raise e
  in

  let mxsub = data.mxsub in
  let mysub = data.mysub in

  for jy = 0 to mysub-1 do
    for ix = 0 to mxsub-1 do
      let uvxy = ij_vptr uv ix jy in
      fprintf fout "%g\n%g\n" uvxy.{0} uvxy.{1};
    done
  done;
  close_out fout




(*
 * PrintFinalStats: Print final run data contained in iopt.
 *)

let print_final_stats mem =
  let nst = Ida.get_num_steps mem in
  let nre = Ida.get_num_res_evals mem in
  let netf = Ida.get_num_err_test_fails mem in
  let ncfn = Ida.get_num_nonlin_solv_conv_fails mem in
  let nni = Ida.get_num_nonlin_solv_iters mem in

  let ncfl = Ida.Spils.get_num_conv_fails mem in
  let nli = Ida.Spils.get_num_lin_iters mem in
  let npe = Ida.Spils.get_num_prec_evals mem in
  let nps = Ida.Spils.get_num_prec_solves mem in
  let nreLS = Ida.Spils.get_num_res_evals mem in

  let nge = Ida_bbd.get_num_gfn_evals mem in

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
  printf "Number of preconditioner solves    = %d\n" nps;
  printf "Number of local residual evals.    = %d\n" nge

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

  (* Set local length (local_N) and global length (system_size). *)
  let local_N = mxsub*mysub*num_species in
  let system_size = neq in

  (* Set up user data block data. *)
  let data = init_user_data thispe npes comm in

  (* Create needed vectors, and load initial values.
     The vector resid is used temporarily only.        *)

  let uv = Nvector_parallel.make local_N system_size comm 0. in

  let uvp = Nvector_parallel.make local_N system_size comm 0. in

  let resid = Nvector_parallel.make local_N system_size comm 0. in

  let id = Nvector_parallel.make local_N system_size comm 0. in

  let uvS = Array.init num_sens (fun _ -> Nvector_parallel.clone uv) in

  let uvpS = Array.init num_sens (fun _ -> Nvector_parallel.clone uv) in

  set_initial_profiles data (unvec uv) (unvec uvp)
    (unvec id) (unvec resid);

  (* Set remaining inputs to IDAS. *)
  let t0 = zero in

  (* Call IDACreate and IDAInit to initialize solution *)
  (* Call IDASpgmr to specify the IDAS LINEAR SOLVER IDASPGMR *)
  (* Call IDABBDPrecInit to initialize the band-block-diagonal preconditioner.
     The half-bandwidths for the difference quotient evaluation are exact
     for the system Jacobian, but only a 5-diagonal band matrix is retained. *)
  let mudq = nsmxsub in
  let mldq = nsmxsub in
  let mukeep = 2 in
  let mlkeep = 2 in
  let maxl = 16 in
  let linsolv =
    Ida_bbd.spgmr (Some maxl)
      { Ida_bbd.mudq = mudq;
        Ida_bbd.mldq = mldq;
        Ida_bbd.mukeep = mukeep;
        Ida_bbd.mlkeep = mlkeep;
      }
      (Some zero)
      {
        Ida_bbd.local_fn = reslocal data;
        Ida_bbd.comm_fn = None;
      }
  in
  let mem =
    Ida.init linsolv (Ida.SStolerances (rtol,atol))
      (res data)
      ~t0:t0 uv uvp
  in

  (* Enable forward sensitivity analysis. *)
  let sparams =
    {
      Sens.pvals = Some data.eps;
      Sens.pbar  = Some (RealArray.clone data.eps);
      Sens.plist = None;
    }
  in
  Sens.init mem Sens.EEtolerances Sens.Simultaneous sparams None uvS uvpS;
  Sens.set_err_con mem true;

  (* Call IDACalcIC (with default options) to correct the initial values. *)
  let tout = ref 0.001 in
  Ida.calc_ic_ya_yd' mem id !tout ~y:uv ~y':uvp;

  (* On PE 0, print heading, basic parameters, initial values. *)
  if thispe = 0 then
    print_header system_size maxl mudq mldq mukeep mlkeep rtol atol
  ;

  print_output mem uv t0 data comm;


  (* Call IDAS in tout loop, normal mode, and print selected output. *)
  for iout = 1 to nout do

    let (tret, _) = Ida.solve_normal mem !tout uv uvp in

    print_output mem uv tret data comm;

    if iout < 3 then tout := !tout *. tmult
    else             tout := !tout +. tadd
  done;
  (* Print each PE's portion of the solution in a separate file. *)
  (* print_sol(mem, uv, uvp, data, comm); *)


  (* On PE 0, print final set of statistics. *)
  if thispe = 0 then
    print_final_stats mem
  ;

  (* calculate integral of u over domain. *)
  let intval = integr comm uv data in
  if thispe = 0 then
    printf "\n\nThe average of u on the domain:\ng = %g\n" intval
  ;

  (* integrate the sensitivities of u over domain. *)
  let _ = Sens.get mem uvS in
  if thispe = 0 then
    printf "\nSensitivities of g:\n"
  ;

  for is = 0 to num_sens-1 do
    let intval = integr comm uvS.(is) data in
    if thispe = 0 then
      printf "w.r.t. eps%d = %14.10f\n" is intval
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
