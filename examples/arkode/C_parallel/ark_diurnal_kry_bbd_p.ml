(*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jan 2016.
 *---------------------------------------------------------------
 * Copyright (c) 2015, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 *---------------------------------------------------------------
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
 * The problem is solved by ARKODE on NPE processors, treated
 * as a rectangular process grid of size NPEX by NPEY, with
 * NPE = NPEX*NPEY. Each processor contains a subgrid of size MXSUB
 * by MYSUB of the (x,y) mesh. Thus the actual mesh sizes are
 * MX = MXSUB*NPEX and MY = MYSUB*NPEY, and the ODE system size is
 * neq = 2*MX*MY.
 *
 * The solution is done with the DIRK/GMRES method (i.e. using the
 * ARKSPGMR linear solver) and a block-diagonal matrix with banded
 * blocks as a preconditioner, using the ARKBBDPRE module.
 * Each block is generated using difference quotients, with
 * half-bandwidths mudq = mldq = 2*MXSUB, but the retained banded
 * blocks have half-bandwidths mukeep = mlkeep = 2.
 * A copy of the approximate Jacobian is saved and conditionally
 * reused within the preconditioner routine.
 *
 * The problem is solved twice -- with left and right preconditioning.
 *
 * Performance data and sampled solution values are printed at
 * selected output times, and all performance counters are printed
 * on completion.
 *
 * This version uses MPI for user routines.
 * Execute with number of processors = NPEX*NPEY (see constants below).
 *-----------------------------------------------------------------
 *)

open Sundials
module ARKStep = Arkode.ARKStep

module BBD = Arkode_bbd
open Bigarray

let lt600 =
  let n, _, _ = Sundials.Config.sundials_version in
  n < 6

let ge670 = match Sundials.Config.version with
            | 6, m, _, _ -> m >= 7
            | m, _, _, _ -> m > 6

let local_array = Nvector_parallel.local_array
let slice = Array1.sub
let printf = Printf.printf
let eprintf = Printf.eprintf

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 0) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 1) []) 0
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

(* ARKodeInit Constants *)

let rtol =     1.0e-5       (* scalar relative tolerance *)
let floor =    100.0        (* value of C1 or C2 at which tolerances *)
                            (* change from relative to absolute      *)
let atol =     rtol*.floor  (* scalar absolute tolerance *)

(* Type : UserData
   contains problem constants, extended dependent variable array,
   grid constants, processor indices, MPI communicator *)

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

    }

(*********************** Private Helper Functions ************************)

(* Load constants in data *)

let sqr x = x *. x

let init_user_data my_pe comm =
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

    uext     = RealArray.make (nvars*(mxsub+2)*(mysub+2)) 0.0;

    (* Set the sizes of a boundary x-line in u and uext *)
    nvmxsub  = nvars*mxsub;
    nvmxsub2 = nvars*(mxsub+2);
  }

(* Set initial conditions in u *)

let set_initial_profiles data u =
  (* Set pointer to data array in vector u *)
  let udata = local_array u in

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

(* Print problem introduction *)

let print_intro npes mudq mldq mukeep mlkeep =
  printf "\n2-species diurnal advection-diffusion problem\n";
  printf "  %d by %d mesh on %d processors\n" mx my npes;
  printf "  Using ARKBBDPRE preconditioner module\n";
  printf "    Difference-quotient half-bandwidths are";
  printf " mudq = %d,  mldq = %d\n" mudq mldq;
  printf "    Retained band block half-bandwidths are";
  printf " mukeep = %d,  mlkeep = %d" mukeep mlkeep

(* Print current t, step count, order, stepsize, and sampled c1,c2 values *)

let print_output s my_pe comm u t =
  let npelast = npex*npey - 1 in
  let tempu = RealArray.create 2 in
  let udata = local_array u in

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
      RealArray.blitn ~src:buf ~dst:tempu 2
    end;

    let nst = ARKStep.get_num_steps s
    and hu  = ARKStep.get_last_step s
    in
    printf "t = %.2e   no. steps = %d   stepsize = %.2e\n" t nst hu;
    printf "At bottom left:  c1, c2 = %12.3e %12.3e \n" udata.{0} udata.{1};
    printf "At top right:    c1, c2 = %12.3e %12.3e \n\n" tempu.{0} tempu.{1}
  end

(* Print final statistics contained in iopt *)

let print_final_stats s =
  let open ARKStep in
  let lenrw, leniw = get_work_space s
  and nst          = get_num_steps s
  and nfe,nfi      = get_num_rhs_evals s
  and nsetups      = get_num_lin_solv_setups s
  and netf         = get_num_err_test_fails s
  and nni          = get_num_nonlin_solv_iters s
  and ncfn         = get_num_nonlin_solv_conv_fails s
  in
  let lenrwLS, leniwLS = Spils.get_work_space s
  and nli   = Spils.get_num_lin_iters s
  and npe   = Spils.get_num_prec_evals s
  and nps   = Spils.get_num_prec_solves s
  and ncfl  = Spils.get_num_lin_conv_fails s
  and nfeLS = Spils.get_num_lin_rhs_evals s
  in
  printf "\nFinal Statistics: \n\n";
  printf "lenrw   = %5d     leniw   = %5d\n"   lenrw leniw;
  printf "lenrwls = %5d     leniwls = %5d\n"   lenrwLS leniwLS;
  printf "nst     = %5d     nfe     = %5d\n"   nst nfe;
  printf "nfe     = %5d     nfels   = %5d\n"   nfi nfeLS;
  printf "nni     = %5d     nli     = %5d\n"   nni nli;
  printf "nsetups = %5d     netf    = %5d\n"   nsetups netf;
  printf "npe     = %5d     nps     = %5d\n"   npe nps;
  printf "ncfn    = %5d     ncfl    = %5d\n\n" ncfn ncfl;

  let lenrwBBDP, leniwBBDP = BBD.get_work_space s in
  let ngevalsBBDP = BBD.get_num_gfn_evals s in
  printf "In ARKBBDPRE: real/integer local work space sizes = %d, %d\n"
                                                          lenrwBBDP leniwBBDP;
  printf "             no. flocal evals. = %d\n" ngevalsBBDP

(* Routine to send boundary data to neighboring PEs *)

let bsend comm my_pe isubx isuby dsizex _ (udata : RealArray.t) =
  let buf = RealArray.create (nvars*mysub) in

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
      RealArray.blitn ~src:udata ~spos:(ly*dsizex)
                      ~dst:buf ~dpos:(ly*nvars)
                      nvars
    done;
    Mpi.send buf (my_pe-1) 0 comm
  end;

  (* If isubx < NPEX-1, send data from right y-line of u (via bufright) *)
  if isubx <> npex-1 then begin
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*nvars in
      let offsetu = offsetbuf*mxsub + (mxsub-1)*nvars in
      RealArray.blitn ~src:udata ~spos:offsetu
                      ~dst:buf ~dpos:offsetbuf
                      nvars
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

let brecvwait request isubx isuby dsizex (uext : RealArray.t) =
  let dsizex2 = dsizex + 2*nvars in

  (* If isuby > 0, receive data for bottom x-line of uext *)
  if isuby <> 0 then begin
    let buf = (Mpi.wait_receive request.(0) : RealArray.t) in
    RealArray.blitn ~src:buf ~dst:uext ~dpos:nvars dsizex
  end;

  (* If isuby < NPEY-1, receive data for top x-line of uext *)
  if isuby <> npey-1 then begin
    let buf = (Mpi.wait_receive request.(1) : RealArray.t) in
    RealArray.blitn ~src:buf
                    ~dst:uext ~dpos:(nvars*(1 + (mysub+1)*(mxsub+2)))
                    dsizex
  end;

  (* If isubx > 0, receive data for left y-line of uext (via bufleft) *)
  if isubx <> 0 then begin
    let bufleft = (Mpi.wait_receive request.(2) : RealArray.t) in
    (* Copy the buffer to uext *)
    for ly = 0 to mysub - 1 do
      let offsetbuf = ly*nvars in
      let offsetue = (ly+1)*dsizex2 in
      RealArray.blitn ~src:bufleft ~spos:offsetbuf
                      ~dst:uext ~dpos:offsetue
                      nvars
    done
  end;

  (* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) *)
  if isubx <> npex-1 then begin
    let bufright = (Mpi.wait_receive request.(3) : RealArray.t) in
    (* Copy the buffer to uext *)
    for ly = 0 to mysub-1 do
      let offsetbuf = ly*nvars in
      let offsetue = (ly+2)*dsizex2 - nvars in
      RealArray.blitn ~src:bufright ~spos:offsetbuf ~dst:uext ~dpos:offsetue nvars
    done
  end

(* fucomm routine. This routine performs all inter-processor
   communication of data in u needed to calculate f.         *)

let fucomm data _ ((udata : RealArray.t),_,_) =
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

let flocal data t ((udata : RealArray.t),_,_) ((dudata : RealArray.t),_,_) =
  (* Get subgrid indices, data sizes, extended work array uext *)
  let isubx    = data.isubx
  and isuby    = data.isuby
  and nvmxsub  = data.nvmxsub
  and nvmxsub2 = data.nvmxsub2
  and uext     = data.uext
  in
  (* Copy local segment of u vector into the working extended array uext *)
  for ly = 0 to mysub-1 do
    RealArray.blitn ~src:udata ~spos:(ly*nvmxsub)
                    ~dst:uext  ~dpos:((ly + 1)*nvmxsub2 + nvars)
                    nvmxsub;
  done;

  (* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of u to uext *)

  (* If isuby = 0, copy x-line 2 of u to uext *)
  if isuby = 0 then RealArray.blitn ~src:udata ~spos:nvmxsub
                                    ~dst:uext  ~dpos:nvars
                                    nvmxsub;

  (* If isuby = NPEY-1, copy x-line MYSUB-1 of u to uext *)
  if isuby = npey-1
    then RealArray.blitn ~src:udata ~spos:((mysub-2)*nvmxsub)
                         ~dst:uext ~dpos:((mysub+1)*nvmxsub2 + nvars)
                         nvmxsub;

  (* If isubx = 0, copy y-line 2 of u to uext *)
  if isubx = 0 then
    for ly = 0 to mysub-1 do
      RealArray.blitn ~src:udata ~spos:(ly*nvmxsub + nvars)
                      ~dst:uext  ~dpos:((ly+1)*nvmxsub2)
                      nvars
    done;

  (* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext *)
  if isubx = npex-1 then
    for ly = 0 to mysub-1 do
      RealArray.blitn ~src:udata ~spos:((ly+1)*nvmxsub - 2*nvars)
                      ~dst:uext  ~dpos:((ly+2)*nvmxsub2 - nvars)
                      nvars
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

(* f routine.  Evaluate f(t,y).  First call fucomm to do communication of
   subgrid boundary data into uext.  Then calculate f by a call to flocal. *)

let f data t u du =
  (* Call fucomm to do inter-processor communication *)
  fucomm data t u;
  (* Call fcalc to calculate all right-hand sides *)
  flocal data t u du

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
      eprintf "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n\n"
                                                              npes (npex*npey);
    exit 1
  end;

  (* Set local length *)
  let local_N = nvars*mxsub*mysub in

  (* Allocate and load user data block *)
  let data = init_user_data my_pe comm in

  (* Allocate u, and set initial values and tolerances *)
  let u = Nvector_parallel.make local_N neq comm 0.0 in
  set_initial_profiles data u;
  let abstol = atol
  and reltol = rtol
  in
  (* Call ARKodeCreate to create the solver memory *)
  let mudq   = nvars * mxsub in
  let mldq   = mudq in
  let mukeep = nvars in
  let mlkeep = mukeep in
  let lsolver = ARKStep.Spils.(spgmr u) in
  let arkode_mem = ARKStep.(
    init
      (implicit (f data)
                 ~lsolver:Spils.(solver lsolver (BBD.prec_left
                                      BBD.({ mudq;   mldq; mukeep; mlkeep})
                            (flocal data))))
      (SStolerances (reltol, abstol))
      t0
      u
  ) in
  ARKStep.set_max_num_steps arkode_mem 10000;
  if ge670 then ARKStep.set_nonlin_conv_coef arkode_mem 0.01;

  (* Print heading *)
  if my_pe = 0 then print_intro npes mudq mldq mukeep mlkeep;

  let solve_problem jpre =
    (* On second run, re-initialize u, the integrator, ARKBBDPRE,
       and ARKSPGMR *)
    if jpre = ARKStep.Spils.PrecRight then begin
      set_initial_profiles data u;
      ARKStep.reinit arkode_mem t0 u;
      BBD.reinit arkode_mem mudq mldq;
      ARKStep.Spils.(set_prec_type lsolver PrecRight);

      if my_pe = 0 then begin
        printf "\n\n-------------------------------------------------------";
        printf "------------\n"
      end
    end;

    if my_pe = 0 then
      printf "\n\nPreconditioner type is:  jpre = %s\n\n"
             (if jpre = ARKStep.Spils.PrecLeft
              then (if lt600 then "PREC_LEFT" else "SUN_PREC_LEFT")
              else (if lt600 then "PREC_RIGHT" else "SUN_PREC_RIGHT"));

    (* In loop over output points, call ARKode, print results, test for error *)
    let tout = ref twohr in
    for _ = 1 to nout do
      let t, _ = ARKStep.evolve_normal arkode_mem !tout u in
      print_output arkode_mem my_pe comm u t;
      tout := !tout +. twohr
    done;

    (* Print final statistics *)
    if my_pe = 0 then print_final_stats arkode_mem
  in
  List.iter solve_problem ARKStep.Spils.([PrecLeft; PrecRight])

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
  for _ = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()
