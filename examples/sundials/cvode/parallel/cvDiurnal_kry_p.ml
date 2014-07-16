(*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2010/12/01 22:52:20 $
 * -----------------------------------------------------------------
 * Programmer(s): S. D. Cohen, A. C. Hindmarsh, M. R. Wittman, and
 *                Radu Serban  @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2014.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(y) = Kv0*exp(y/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * mesh, with simple polynomial initial profiles.
 *
 * The problem is solved by CVODE on NPE processors, treated
 * as a rectangular process grid of size NPEX by NPEY, with
 * NPE = NPEX*NPEY. Each processor contains a subgrid of size MXSUB
 * by MYSUB of the (x,y) mesh.  Thus the actual mesh sizes are
 * MX = MXSUB*NPEX and MY = MYSUB*NPEY, and the ODE system size is
 * neq = 2*MX*MY.
 *
 * The solution is done with the BDF/GMRES method (i.e. using the
 * CVSPGMR linear solver) and the block-diagonal part of the
 * Newton matrix as a left preconditioner. A copy of the
 * block-diagonal part of the Jacobian is saved and conditionally
 * reused within the preconditioner routine.
 *
 * Performance data and sampled solution values are printed at
 * selected output times, and all performance counters are printed
 * on completion.
 *
 * This version uses MPI for user routines.
 * 
 * Execution: mpirun -np N cvDiurnal_kry_p   with N = NPEX*NPEY
 * (see constants below).
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module Roots  = Sundials.Roots
module Direct = Dls.ArrayDenseMatrix
module BandedSpils = Cvode.Spils.Banded
open Bigarray

let unvec = Sundials.unvec
let slice = Array1.sub
let printf = Printf.printf
let eprintf = Printf.eprintf

let blit buf buf_offset dst dst_offset len =
  for i = 0 to len-1 do
    dst.{dst_offset + i} <- buf.{buf_offset + i}
  done

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_string (RealArray.empty) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_string (RealArray.make 1) []) 0
  - header_and_empty_array_size

let bytes x = header_and_empty_array_size + x * float_cell_size

(* Problem Constants *)

let nvars =    2            (* number of species         *)
let kh =       4.0e-6       (* horizontal diffusivity Kh *)
let vel =      0.001        (* advection velocity V      *)
let kv0 =      1.0e-8       (* coefficient in Kv(y)      *)
let q1 =       1.63e-16     (* coefficients q1, q2, c3   *) 
let q2 =       4.66e-16
let c3 =       3.7e16
let a3 =       22.62        (* coefficient in expression for q3(t) *)
let a4 =       7.601        (* coefficient in expression for q4(t) *)
let c1_scale = 1.0e6        (* coefficients in initial profiles    *)
let c2_scale = 1.0e12

let t0 =       0.0          (* initial time *)
let nout =     12           (* number of output times *)
let twohr =    7200.0       (* number of seconds in two hours  *)
let halfday =  4.32e4       (* number of seconds in a half day *)
let pi =       3.1415926535898  (* pi *) 

let xmin =     0.0          (* grid boundaries in x  *)
let xmax =     20.0           
let ymin =     30.0         (* grid boundaries in y  *)
let ymax =     50.0

let npex =     2            (* no. PEs in x direction of PE array *)
let npey =     2            (* no. PEs in y direction of PE array *)
                            (* Total no. PEs = NPEX*NPEY *)
let mxsub =    5            (* no. x points per subgrid *)
let mysub =    5            (* no. y points per subgrid *)

let mx =       npex*mxsub   (* MX = number of x mesh points *)
let my =       npey*mysub   (* MY = number of y mesh points *)
                            (* Spatial mesh is MX by MY *)

(* CVodeMalloc Constants *)

let rtol =     1.0e-5       (* scalar relative tolerance *)
let floor =    100.0        (* value of C1 or C2 at which tolerances *)
                            (* change from relative to absolute      *)
let atol =     rtol*.floor  (* scalar absolute tolerance *)

(* User-defined matrix accessor macro: IJth *)

(* IJth is defined in order to write code which indexes into dense
   matrices with a (row,column) pair, where 1 <= row,column <= NVARS.   

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NVARS. The small matrix routines in sundials_dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. *)

let ijth v i j       = Direct.get v (i - 1) (j - 1)
let set_ijth v i j e = Direct.set v (i - 1) (j - 1) e

(* Type : UserData 
   contains problem constants, preconditioner blocks, pivot arrays, 
   grid constants, and processor indices, as well as data needed
   for the preconditiner *)

type user_data = {

        mutable q4 : float;
        om         : float;
        dx         : float;
        dy         : float;
        hdco       : float;
        haco       : float;
        vdco       : float;

        uext       : RealArray.t;
        
        my_pe      : int;
        isubx      : int;
        isuby      : int;

        nvmxsub    : int;
        nvmxsub2   : int;

        comm       : Mpi.communicator;

        (* For preconditioner *)
        p          : Direct.t array array;
        jbd        : Direct.t array array;
        pivot      : Sundials.LintArray.t array array;

    }

(*********************** Private Helper Functions ************************)

(* Load constants in data *)

let sqr x = x ** 2.0

let init_user_data my_pe comm =
  let new_dmat _ = Direct.make nvars nvars in
  let new_int1 _  = Sundials.LintArray.make nvars in
  let new_y_arr elinit _ = Array.init mysub elinit in
  let new_xy_arr elinit  = Array.init mxsub (new_y_arr elinit) in

  let dx    = (xmax-.xmin)/.(float (mx-1)) in
  let dy    = (ymax-.ymin)/.(float (my-1)) in
  let isuby = my_pe/npex in

  {
    q4       = 0.0; (* set later *)

    (* Set problem constants *)
    om       = pi/.halfday;
    dx       = dx;
    dy       = dy;
    hdco     = kh/.sqr(dx);
    haco     = vel/.(2.0*.dx);
    vdco     = (1.0/.sqr(dy))*.kv0;

    (* Set machine-related constants *)
    comm     = comm;
    my_pe    = my_pe;

    (* isubx and isuby are the PE grid indices corresponding to my_pe *)
    isuby    = isuby;
    isubx    = my_pe - isuby*npex;

    uext     = RealArray.init (nvars*(mxsub+2)*(mysub+2)) 0.0;

    (* Set the sizes of a boundary x-line in u and uext *)
    nvmxsub  = nvars*mxsub;
    nvmxsub2 = nvars*(mxsub+2);

    (* Preconditioner-related fields *)
    p     = new_xy_arr new_dmat;
    jbd   = new_xy_arr new_dmat;
    pivot = new_xy_arr new_int1;
  }

(* Set initial conditions in u *)

let set_initial_profiles data u =
  (* Set pointer to data array in vector u *)
  let udata, _, _ = unvec u in

  (* Get mesh spacings, and subgrid indices for this PE *)
  let dx = data.dx
  and dy = data.dy
  and isubx = data.isubx
  and isuby = data.isuby
  in
  (* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. *)
  let offset = ref 0 in
  let xmid = 0.5*.(xmin +. xmax) in
  let ymid = 0.5*.(ymin +. ymax) in
  for ly = 0 to mysub-1 do
    let jy = ly + isuby*mysub in
    let y = ymin +. (float jy)*.dy in
    let cy = sqr(0.1*.(y -. ymid)) in
    let cy = 1.0 -. cy +. 0.5*.(sqr cy) in
    for lx = 0 to mxsub-1 do
      let jx = lx + isubx*mxsub in
      let x  = xmin +. (float jx)*.dx in
      let cx = sqr(0.1*.(x -. xmid)) in
      let cx = 1.0 -. cx +. 0.5*.(sqr cx) in
      udata.{!offset  } <- c1_scale *. cx *. cy; 
      udata.{!offset+1} <- c2_scale *. cx *. cy;
      offset := !offset + 2
    done
  done

(* Print current t, step count, order, stepsize, and sampled c1,c2 values *)

let print_output s my_pe comm u t =
  let npelast = npex*npey - 1 in
  let tempu = RealArray.make 2 in
  let udata, _, _ = unvec u in

  (* Send c1,c2 at top right mesh point to PE 0 *)
  if my_pe = npelast then begin
    let i0 = nvars*mxsub*mysub - 2 in
    let i1 = i0 + 1 in
    if npelast <> 0 then Mpi.send (slice udata i0 2) 0 0 comm
    else (tempu.{0} <- udata.{i0}; tempu.{1} <- udata.{i1})
  end;

  (* On PE 0, receive c1,c2 at top right, then print performance data
     and sampled solution values *) 
  if my_pe = 0 then begin
    if npelast <> 0 then begin
      let buf = (Mpi.receive npelast 0 comm : RealArray.t) in
      blit buf 0 tempu 0 2
    end;

    let nst = Cvode.get_num_steps s
    and qu  = Cvode.get_last_order s
    and hu  = Cvode.get_last_step s
    in
    printf "t = %.2e   no. steps = %d   order = %d   stepsize = %.2e\n" 
                                                                  t nst qu hu;
    printf "At bottom left:  c1, c2 = %12.3e %12.3e \n" udata.{0} udata.{1};
    printf "At top right:    c1, c2 = %12.3e %12.3e \n\n" tempu.{0} tempu.{1}
  end

(* Print final statistics contained in iopt *)

let print_final_stats s =
  let lenrw, leniw = Cvode.get_work_space s
  and nst          = Cvode.get_num_steps s
  and nfe          = Cvode.get_num_rhs_evals s
  and nsetups      = Cvode.get_num_lin_solv_setups s
  and netf         = Cvode.get_num_err_test_fails s
  and nni          = Cvode.get_num_nonlin_solv_iters s
  and ncfn         = Cvode.get_num_nonlin_solv_conv_fails s
  in
  let lenrwLS, leniwLS = Cvode.Spils.get_work_space s
  and nli   = Cvode.Spils.get_num_lin_iters s
  and npe   = Cvode.Spils.get_num_prec_evals s
  and nps   = Cvode.Spils.get_num_prec_solves s
  and ncfl  = Cvode.Spils.get_num_conv_fails s
  and nfeLS = Cvode.Spils.get_num_rhs_evals s
  in
  printf "\nFinal Statistics: \n\n";
  printf "lenrw   = %5d     leniw   = %5d\n"   lenrw leniw;
  printf "lenrwls = %5d     leniwls = %5d\n"   lenrwLS leniwLS;
  printf "nst     = %5d\n"                      nst;
  printf "nfe     = %5d     nfels   = %5d\n"   nfe nfeLS;
  printf "nni     = %5d     nli     = %5d\n"   nni nli;
  printf "nsetups = %5d     netf    = %5d\n"   nsetups netf;
  printf "npe     = %5d     nps     = %5d\n"   npe nps;
  printf "ncfn    = %5d     ncfl    = %5d\n\n" ncfn ncfl
 
(* Routine to send boundary data to neighboring PEs *)

let bsend comm my_pe isubx isuby dsizex dsizey udata =
  let buf = RealArray.make (nvars*mysub) in

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
      blit udata (ly*dsizex) buf (ly*nvars) nvars
    done;
    Mpi.send buf (my_pe-1) 0 comm
  end;

  (* If isubx < NPEX-1, send data from right y-line of u (via bufright) *)
  if isubx <> npex-1 then begin
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*nvars in
      let offsetu = offsetbuf*mxsub + (mxsub-1)*nvars in
      blit udata offsetu buf offsetbuf nvars
    done;
    Mpi.send buf (my_pe+1) 0 comm
  end
 
(* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. *)

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
  Array.of_list [r0; r1; r2; r3]

(* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. *)

let brecvwait request isubx isuby dsizex uext =
  let dsizex2 = dsizex + 2*nvars in

  (* If isuby > 0, receive data for bottom x-line of uext *)
  if isuby <> 0 then begin
    let buf = (Mpi.wait_receive request.(0) : RealArray.t) in
    blit buf 0 uext nvars dsizex
  end;

  (* If isuby < NPEY-1, receive data for top x-line of uext *)
  if isuby <> npey-1 then begin
    let buf = (Mpi.wait_receive request.(1) : RealArray.t) in
    blit buf 0 uext (nvars*(1 + (mysub+1)*(mxsub+2))) dsizex
  end;

  (* If isubx > 0, receive data for left y-line of uext (via bufleft) *)
  if isubx <> 0 then begin
    let bufleft = (Mpi.wait_receive request.(2) : RealArray.t) in
    (* Copy the buffer to uext *)
    for ly = 0 to mysub - 1 do
      let offsetbuf = ly*nvars in
      let offsetue = (ly+1)*dsizex2 in
      blit bufleft offsetbuf uext offsetue nvars
    done
  end;

  (* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) *)
  if isubx <> npex-1 then begin
    let bufright = (Mpi.wait_receive request.(3) : RealArray.t) in
    (* Copy the buffer to uext *)
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*nvars in
      let offsetue = (ly+2)*dsizex2 - nvars in
      blit bufright offsetbuf uext offsetue nvars
    done
  end

(* ucomm routine.  This routine performs all communication 
   between processors of data needed to calculate f. *)

let ucomm data t udata =
  let comm    = data.comm
  and my_pe   = data.my_pe
  and isubx   = data.isubx
  and isuby   = data.isuby
  and nvmxsub = data.nvmxsub
  and nvmysub = nvars*mysub
  and uext    = data.uext
  in
  (* Start receiving boundary data from neighboring PEs *)
  let request = brecvpost comm my_pe isubx isuby nvmxsub nvmysub in

  (* Send data from boundary of local grid to neighboring PEs *)
  bsend comm my_pe isubx isuby nvmxsub nvmysub udata;

  (* Finish receiving boundary data from neighboring PEs *)
  brecvwait request isubx isuby nvmxsub uext

(* fcalc routine. Compute f(t,y).  This routine assumes that communication 
   between processors of data needed to calculate f has already been done,
   and this data is in the work array uext. *)

let fcalc data t udata dudata =
  (* Get subgrid indices, data sizes, extended work array uext *)
  let isubx    = data.isubx
  and isuby    = data.isuby
  and nvmxsub  = data.nvmxsub
  and nvmxsub2 = data.nvmxsub2
  and uext     = data.uext
  in
  (* Copy local segment of u vector into the working extended array uext *)
  for ly = 0 to mysub-1 do
    blit udata (ly*nvmxsub) uext ((ly + 1)*nvmxsub2 + nvars) nvmxsub;
  done;

  (* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of u to uext *)

  (* If isuby = 0, copy x-line 2 of u to uext *)
  if isuby = 0 then blit udata nvmxsub uext nvars nvmxsub;

  (* If isuby = NPEY-1, copy x-line MYSUB-1 of u to uext *)
  if isuby = npey-1 then blit udata ((mysub-2)*nvmxsub)
                              uext  ((mysub+1)*nvmxsub2 + nvars) nvmxsub;

  (* If isubx = 0, copy y-line 2 of u to uext *)
  if isubx = 0 then
    for ly = 0 to mysub-1 do
      blit udata (ly*nvmxsub + nvars) uext ((ly+1)*nvmxsub2) nvars
    done;

  (* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext *)
  if isubx = npex-1 then
    for ly = 0 to mysub-1 do
      blit udata ((ly+1)*nvmxsub - 2*nvars) uext ((ly+2)*nvmxsub2 - nvars) nvars
    done;

  (* Make local copies of problem variables, for efficiency *)
  let dely   = data.dy
  and verdco = data.vdco
  and hordco = data.hdco
  and horaco = data.haco
  in
  (* Set diurnal rate coefficients as functions of t, and save q4 in 
  data block for use by preconditioner evaluation routine *)
  let s = sin(data.om *. t) in
  let q3, q4coef =
    if s > 0.0  then (exp(-.a3/.s), exp(-.a4/.s)) else (0.0, 0.0)
  in
  data.q4 <- q4coef;

  (* Loop over all grid points in local subgrid *)
  for ly = 0 to mysub-1 do
    let jy = ly + isuby*mysub in

    (* Set vertical diffusion coefficients at jy +- 1/2 *)
    let ydn = ymin +. (float jy -. 0.5)*.dely in
    let yup = ydn +. dely in
    let cydn = verdco*.exp(0.2*.ydn) in
    let cyup = verdco*.exp(0.2*.yup) in

    for lx = 0 to mxsub-1 do
      (* Extract c1 and c2, and set kinetic rate terms *)
      let offsetue = (lx+1)*nvars + (ly+1)*nvmxsub2 in
      let c1 = uext.{offsetue} in
      let c2 = uext.{offsetue+1} in
      let qq1 = q1*.c1*.c3 in
      let qq2 = q2*.c1*.c2 in
      let qq3 = q3*.c3 in
      let qq4 = q4coef*.c2 in
      let rkin1 = -.qq1 -. qq2 +. 2.0*.qq3 +. qq4 in
      let rkin2 = qq1 -. qq2 -. qq4 in

      (* Set vertical diffusion terms *)
      let c1dn = uext.{offsetue-nvmxsub2} in
      let c2dn = uext.{offsetue-nvmxsub2+1} in
      let c1up = uext.{offsetue+nvmxsub2} in
      let c2up = uext.{offsetue+nvmxsub2+1} in
      let vertd1 = cyup*.(c1up -. c1) -. cydn*.(c1 -. c1dn) in
      let vertd2 = cyup*.(c2up -. c2) -. cydn*.(c2 -. c2dn) in

      (* Set horizontal diffusion and advection terms *)
      let c1lt = uext.{offsetue-2} in
      let c2lt = uext.{offsetue-1} in
      let c1rt = uext.{offsetue+2} in
      let c2rt = uext.{offsetue+3} in
      let hord1 = hordco*.(c1rt -. 2.0*.c1 +. c1lt) in
      let hord2 = hordco*.(c2rt -. 2.0*.c2 +. c2lt) in
      let horad1 = horaco*.(c1rt -. c1lt) in
      let horad2 = horaco*.(c2rt -. c2lt) in
      (* Load all terms into dudata *)
      let offsetu = lx*nvars + ly*nvmxsub in
      dudata.{offsetu}   <- vertd1 +. hord1 +. horad1 +. rkin1; 
      dudata.{offsetu+1} <- vertd2 +. hord2 +. horad2 +. rkin2
    done
  done

(***************** Functions Called by the Solver *************************)

(* f routine.  Evaluate f(t,y).  First call ucomm to do communication of 
   subgrid boundary data into uext.  Then calculate f by a call to fcalc. *)

let f data t ((udata : RealArray.t),_,_) ((dudata : RealArray.t),_,_) =
  (* Call ucomm to do inter-processor communication *)
  ucomm data t udata;
  (* Call fcalc to calculate all right-hand sides *)
  fcalc data t udata dudata

(* Preconditioner setup routine. Generate and preprocess P. *)
let precond data jacarg jok gamma =
  let { Cvode.jac_t   = tn;
        Cvode.jac_y   = (udata, _, _);
      } = jacarg
  in
  (* Make local copies of pointers in user_data, and of pointer to u's data *)
  let p       = data.p
  and jbd     = data.jbd
  and pivot   = data.pivot
  and isuby   = data.isuby
  and nvmxsub = data.nvmxsub
  in

  let r =
    if jok then begin
      (* jok = TRUE: Copy Jbd to P *)
      for ly = 0 to mysub - 1 do
        for lx = 0 to mxsub - 1 do
          Direct.copy jbd.(lx).(ly) p.(lx).(ly)
        done
      done;
      false
    end
    else begin
      (* jok = FALSE: Generate Jbd from scratch and copy to P *)
      (* Make local copies of problem variables, for efficiency. *)
      let q4coef = data.q4
      and dely   = data.dy
      and verdco = data.vdco
      and hordco = data.hdco
      in
      (* Compute 2x2 diagonal Jacobian blocks (using q4 values 
         computed on the last f call).  Load into P. *)
      for ly = 0 to mysub - 1 do
        let jy  = ly + isuby*mysub in
        let ydn = ymin +. (float jy -. 0.5) *. dely in
        let yup = ydn +. dely in
        let cydn = verdco *. exp(0.2 *. ydn)
        and cyup = verdco *. exp(0.2 *. yup)
        in
        let diag = -. (cydn +. cyup +. 2.0 *. hordco) in

        for lx = 0 to mxsub - 1 do
          let offset = lx*nvars + ly*nvmxsub in
          let c1     = udata.{offset} in
          let c2     = udata.{offset+1} in
          let j      = jbd.(lx).(ly) in
          let a      = p.(lx).(ly) in
          set_ijth j 1 1 ((-. q1 *. c3 -. q2 *. c2) +. diag);
          set_ijth j 1 2 (-. q2 *. c1 +. q4coef);
          set_ijth j 2 1 (q1 *. c3 -. q2 *. c2);
          set_ijth j 2 2 ((-. q2 *. c1 -. q4coef) +. diag);
          Direct.copy j a
        done
      done;
      true
    end
  in

  (* Scale by -gamma *)
  for ly = 0 to mysub - 1 do
    for lx = 0 to mxsub - 1 do
      Direct.scale (-. gamma) p.(lx).(ly)
    done
  done;
  
  (* Add identity matrix and do LU decompositions on blocks in place. *)
  for lx = 0 to mxsub - 1 do
    for ly = 0 to mysub - 1 do
      Direct.add_identity p.(lx).(ly);
      Direct.getrf p.(lx).(ly) pivot.(lx).(ly)
    done
  done;
  r

(* Preconditioner solve routine *)

let psolve data jac_arg solve_arg (zdata, _, _) =
  let { Cvode.Spils.rhs = (r, _, _);
        Cvode.Spils.gamma = gamma;
        Cvode.Spils.delta = delta;
        Cvode.Spils.left = lr } = solve_arg
  in
  (* Extract the P and pivot arrays from user_data. *)
  let p       = data.p
  and pivot   = data.pivot
  and nvmxsub = data.nvmxsub
  in

  Bigarray.Array1.blit r zdata;

  (* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. *)
  for lx = 0 to mxsub - 1 do
    for ly = 0 to mysub - 1 do
      Direct.getrs p.(lx).(ly) pivot.(lx).(ly)
        (Array1.sub zdata (lx*nvars + ly*nvmxsub) nvars)
    done
  done

(***************************** Main Program ******************************)

let main () =
  (* Set problem size neq *)
  let neq = nvars*mx*my in

  (* Get processor number and total number of pe's *)
  let comm  = Mpi.comm_world in
  let npes  = Mpi.comm_size comm in
  let my_pe = Mpi.comm_rank comm in

  if npes <> npex*npey then begin
    if my_pe = 0 then
      eprintf "\nMPI_ERROR(0: npes = %d is not equal to NPEX*NPEY = %d\n\n"
                                                              npes (npex*npey);
    exit 1
  end;

  (* Set local length *)
  let local_N = nvars*mxsub*mysub in

  (* Allocate and load user data block; allocate preconditioner block *)
  let data = init_user_data my_pe comm in

  (* Allocate u, and set initial values and tolerances *) 
  let u = Nvector_parallel.make local_N neq comm 0.0 in
  set_initial_profiles data u;
  let abstol = atol
  and reltol = rtol
  in
  (* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration *)
  let cvode_mem =
    Cvode.init Cvode.BDF
      (Cvode.Newton
        (Cvode.Spils.spgmr
                None
                Spils.PrecLeft
                { Cvode.Spils.prec_setup_fn = Some (precond data);
                  Cvode.Spils.prec_solve_fn = Some (psolve data);
                  Cvode.Spils.jac_times_vec_fn = None }))
      (Cvode.SStolerances (reltol, abstol))
      (f data) ~t0:t0 u
  in
    
  if my_pe = 0 then
    printf "\n2-species diurnal advection-diffusion problem\n\n";

  (* In loop over output points, call CVode, print results, test for error *)
  let tout = ref twohr in
  for iout=1 to nout do
    let (t, flag) = Cvode.solve_normal cvode_mem !tout u in
    print_output cvode_mem my_pe comm u t;
    tout := !tout +. twohr
  done;

  (* Print final statistics *)  
  if my_pe = 0 then print_final_stats cvode_mem

let n =
  match Sys.argv with
  | [|_; n|] -> int_of_string n
  | _ -> 1
let _ = for i = 1 to n do main () done

let _ = Gc.full_major ()

