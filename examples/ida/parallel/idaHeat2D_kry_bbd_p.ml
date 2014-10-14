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
 * Example problem for IDA: 2D heat equation, parallel, GMRES,
 * IDABBDPRE.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver IDASpgmr and BBD
 * preconditioning.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
 * PDE is treated with central differences on a uniform MX x MY
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = MX * MY. Here MX = MY = 10.
 *
 * The system is actually implemented on submeshes, processor by
 * processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
 * processors.
 *
 * The system is solved with IDA using the Krylov linear solver
 * IDASPGMR in conjunction with the preconditioner module IDABBDPRE.
 * The preconditioner uses a tridiagonal approximation
 * (half-bandwidths = 1). The constraints u >= 0 are posed for all
 * components. Local error testing on the boundary values is
 * suppressed. Output is taken at t = 0, .01, .02, .04, ..., 10.24.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module RealArray2 = Sundials.RealArray2
module LintArray = Sundials.LintArray

let fprintf = Printf.fprintf
let printf = Printf.printf

let vconst = Nvector_parallel.DataOps.n_vconst
let vscale = Nvector_parallel.DataOps.n_vscale
let vmaxnorm = Nvector_parallel.Ops.n_vmaxnorm

let slice = Bigarray.Array1.sub

let blit buf buf_offset dst dst_offset len =
  for i = 0 to len-1 do
    dst.{dst_offset + i} <- buf.{buf_offset + i}
  done

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_string (RealArray.create 0) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_string (RealArray.create 1) []) 0
  - header_and_empty_array_size

let bytes x = header_and_empty_array_size + x * float_cell_size

(* Drop the first i elements of a RealArray.t *)
let real_array_drop i a = slice a i (RealArray.length a - i)



let zero =  0.0
let one =   1.0
let two =   2.0

let nout =         11             (* Number of output times *)

let npex =         2              (* No. PEs in x direction of PE array *)
let npey =         2              (* No. PEs in y direction of PE array *)
                                    (* Total no. PEs = npex*my *)
let mxsub =        5              (* No. x points per subgrid *)
let mysub =        5              (* No. y points per subgrid *)

let mx =           (npex*mxsub)   (* mx = number of x mesh points *)
let my =           (npey*mysub)   (* my = number of y mesh points *)
                                    (* Spatial mesh is mx by my *)

type user_data =
  {
    thispe : int;
    mx : int;
    my : int;
    ixsub : int;
    jysub : int;
    npex : int;
    npey : int;
    mxsub : int;
    mysub : int;
    n_local : int;
    dx : float;
    dy : float;
    coeffx : float;
    coeffy : float;
    coeffxy : float;
    uext : RealArray.t;               (* size = (mxsub+2)*(mysub+2) *)
    comm : Mpi.communicator;
  }

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 *)

(*
 * Routine to send boundary data to neighboring PEs.
 *)

let bsend comm thispe ixsub jysub dsizex dsizey uarray =
  let bufleft = RealArray.create mysub
  and bufright = RealArray.create mysub
  in

  (* If jysub > 0, send data from bottom x-line of u. *)

  if jysub <> 0 then
    Mpi.send (slice uarray 0 dsizex) (thispe-npex) 0 comm;

  (* If jysub < npey-1, send data from top x-line of u. *)

  if jysub <> npey-1 then begin
    let offsetu = (mysub-1)*dsizex in
    Mpi.send (slice uarray offsetu dsizex) (thispe+npex) 0 comm
  end;

  (* If ixsub > 0, send data from left y-line of u (via bufleft). *)

  if ixsub <> 0 then begin
    for ly = 0 to mysub-1 do
      let offsetu = ly*dsizex in
      bufleft.{ly} <- uarray.{offsetu}
    done;
    Mpi.send (slice bufleft 0 dsizey) (thispe-1) 0 comm
  end;

  (* If ixsub < npex-1, send data from right y-line of u (via bufright). *)

  if ixsub <> npex-1 then begin
    for ly = 0 to mysub-1 do
      let offsetu = ly*mxsub + (mxsub-1) in
      bufright.{ly} <- uarray.{offsetu}
    done;
    Mpi.send (slice bufright 0 dsizey) (thispe+1) 0 comm
  end

(*
 * Routine to start receiving boundary data from neighboring PEs.
 * Notes:
 *   1) buffer should be able to hold 2*mysub realtype entries, should be
 *      passed to both the BRecvPost and BRecvWait functions, and should not
 *      be manipulated between the two calls.
 *   2) requests returned from brecvpost will have 4 entries, and should be
 *      passed in both a call to brecvwait.
 *)

let brecvpost comm thispe ixsub jysub dsizex dsizey =
  (* If jysub > 0, receive data for bottom x-line of uext. *)
  let r0 = if jysub <> 0
           then Mpi.ireceive (bytes dsizex) (thispe-npex) 0 comm
           else Mpi.null_request
  in

  (* If jysub < npey-1, receive data for top x-line of uext. *)
  let r1 = if jysub <> npey-1
           then Mpi.ireceive (bytes dsizex) (thispe+npex) 0 comm
           else Mpi.null_request
  in

  (* If ixsub > 0, receive data for left y-line of uext (via bufleft). *)
  let r2 = if ixsub <> 0
           then Mpi.ireceive (bytes dsizey) (thispe-1) 0 comm
           else Mpi.null_request
  in

  (* If ixsub < npex-1, receive data for right y-line of uext (via bufright). *)
  let r3 = if ixsub <> npex-1
    then Mpi.ireceive (bytes dsizey) (thispe+1) 0 comm
    else Mpi.null_request
  in

  [|r0; r1; r2; r3|]

(*
 * Routine to finish receiving boundary data from neighboring PEs.
 * Notes:
 *   1) buffer should be able to hold 2*mysub realtype entries, should be
 *      passed to both the BRecvPost and BRecvWait functions, and should not
 *      be manipulated between the two calls.
 *   2) request should have four entries, and should be passed in both
 *      calls also.
 *)

let brecvwait requests ixsub jysub dsizex uext =
  let dsizex2 = dsizex + 2 in

  (* If jysub > 0, receive data for bottom x-line of uext. *)
  if jysub <> 0 then begin
    let buf = (Mpi.wait_receive requests.(0) : RealArray.t) in
    blit buf 0 uext 1 dsizex
  end;

  (* If jysub < npey-1, receive data for top x-line of uext. *)
  if jysub <> npey-1 then begin
    let offsetue = 1 + (mysub+1)*(mxsub+2) in
    let buf = (Mpi.wait_receive requests.(1) : RealArray.t) in
    blit buf 0 uext offsetue dsizex
  end;

  (* If ixsub > 0, receive data for left y-line of uext (via bufleft). *)
  if ixsub <> 0 then begin
    let bufleft = (Mpi.wait_receive requests.(2) : RealArray.t) in
    (* Copy the buffer to uext. *)
    for ly = 0 to mysub-1 do
      let offsetue = (ly+1)*dsizex2 in
      uext.{offsetue} <- bufleft.{ly}
    done
  end;

  (* If ixsub < npex-1, receive data for right y-line of uext (via bufright). *)
  if ixsub <> npex-1 then begin
    let bufright = (Mpi.wait_receive requests.(3) : RealArray.t) in
    (* Copy the buffer to uext *)
    for ly = 0 to mysub-1 do
      let offsetue = (ly+2)*dsizex2 - 1 in
      uext.{offsetue} <- bufright.{ly}
    done
  end

(*
 * reslocal routine.  Compute res = F(t, uu, up).  This routine assumes
 * that all inter-processor communication of data needed to calculate F
 *  has already been done, and that this data is in the work array uext.
 *)

let reslocal data tres uu up res =
  (* Get subgrid indices, array sizes, extended work array uext. *)

  let uext = data.uext in
  let (uuv,_,_) = uu in
  let (upv,_,_) = up in
  let (resv,_,_) = res in
  let ixsub = data.ixsub and jysub = data.jysub in
  let mxsub = data.mxsub and mxsub2 = data.mxsub + 2 in
  let mysub = data.mysub and npex = data.npex and npey = data.npey in

  (* Initialize all elements of res to uu. This sets the boundary
     elements simply without indexing hassles. *)

  vscale one uu res;

  (* Copy local segment of u vector into the working extended array uext.
     This completes uext prior to the computation of the res vector.     *)

  let offsetu = ref 0 in
  let offsetue = ref (mxsub2 + 1) in
  for ly = 0 to mysub-1 do
    for lx = 0 to mxsub-1 do
      uext.{!offsetue+lx} <- uuv.{!offsetu+lx}
    done;
    offsetu := !offsetu + mxsub;
    offsetue := !offsetue + mxsub2
  done;

  (* Set loop limits for the interior of the local subgrid. *)

  let ixbegin = if ixsub = 0 then 1 else 0 in
  let ixend = if ixsub = npex-1 then mxsub-2 else mxsub-1 in
  let jybegin = if jysub = 0 then 1 else 0 in
  let jyend = if jysub = npey-1 then mysub-2 else mysub-1 in

  (* Loop over all grid points in local subgrid. *)

  for ly = jybegin to jyend do
    for lx = ixbegin to ixend do
      let locu = lx + ly*mxsub in
      let locue = (lx+1) + (ly+1)*mxsub2 in
      let termx = data.coeffx *.(uext.{locue-1}      +. uext.{locue+1}) in
      let termy = data.coeffy *.(uext.{locue-mxsub2} +. uext.{locue+mxsub2}) in
      let termctr = data.coeffxy*.uext.{locue} in
      resv.{locu} <- upv.{locu} -. (termx +. termy -. termctr)
   done
  done

(*
 * rescomm routine.  This routine performs all inter-processor
 * communication of data in u needed to calculate G.
 *)

let rescomm data tt uu up =
  let uarray,_,_ = uu in

  (* Get comm, thispe, subgrid indices, data sizes, extended array uext. *)
  let comm = data.comm and  thispe = data.thispe in
  let ixsub = data.ixsub and   jysub = data.jysub in
  let mxsub = data.mxsub and   mysub = data.mysub in
  let uext = data.uext in

  (* Start receiving boundary data from neighboring PEs. *)
  let requests = brecvpost comm thispe ixsub jysub mxsub mysub in

  (* Send data from boundary of local grid to neighboring PEs. *)
  bsend comm thispe ixsub jysub mxsub mysub uarray;

  (* Finish receiving boundary data from neighboring PEs. *)
  brecvwait requests ixsub jysub mxsub uext


(*
 * heatres: heat equation system residual function
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 *
 * This parallel implementation uses several supporting routines.
 * First a call is made to rescomm to do communication of subgrid boundary
 * data into array uext.  Then reslocal is called to compute the residual
 * on individual processors and their corresponding domains.  The routines
 * BSend, BRecvPost, and BREcvWait handle interprocessor communication
 * of uu required to calculate the residual.
 *)

let heatres data tres uu up res =
  let _Nlocal = data.n_local in

  (* Call rescomm to do inter-processor communication. *)
  rescomm data tres uu up;

  (* Call reslocal to calculate res. *)
  reslocal data tres uu up res


(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(*
 * InitUserData initializes the user's data block data.
 *)

let init_user_data thispe comm =
  let dx = one/.(float_of_int mx-.one) in (* Assumes a [0,1] interval in x. *)
  let dy = one/.(float_of_int my-.one) in (* Assumes a [0,1] interval in y. *)
  let coeffx = one/.(dx *. dx) in
  let coeffy = one/.(dy *. dy) in
  let coeffxy = two/.(dx *. dx) +. two/.(dy *. dy)  in
  let jysub = thispe/npex in
  let ixsub = thispe - jysub * npex in
  let n_local = mxsub*mysub in
  {
    thispe = thispe;
    mx = mx;
    my = my;
    ixsub = ixsub;
    jysub = jysub;
    npex = npex;
    npey = npey;
    mxsub = mxsub;
    mysub = mysub;
    n_local = n_local;
    dx = dx;
    dy = dy;
    coeffx = coeffx;
    coeffy = coeffy;
    coeffxy = coeffxy;
    uext = RealArray.create ((mxsub+2)*(mysub+2));
    comm = comm;
  }


(*
 * SetInitialProfile sets the initial values for the problem.
 *)

let set_initial_profile data uu up id res =
  (* Initialize uu. *)

  let udata,_,_ = uu in
  let iddata,_,_ = id in

  (* Set mesh spacings and subgrid indices for this PE. *)
  let ixsub = data.ixsub in
  let jysub = data.jysub in

  (* Set beginning and ending locations in the global array corresponding
     to the portion of that array assigned to this processor. *)
  let ixbegin = mxsub*ixsub in
  let ixend = mxsub*(ixsub+1) - 1 in
  let jybegin = mysub*jysub in
  let jyend = mysub*(jysub+1) - 1 in

  (* Loop over the local array, computing the initial profile value.
     The global indices are (i,j) and the local indices are (iloc,jloc).
     Also set the id vector to zero for boundary points, one otherwise. *)

  vconst one id;
  for jloc = 0 to jyend-jybegin do
    let j = jybegin + jloc in
    let yfact = data.dy*.float_of_int j in
    let offset = jloc*mxsub in
    for iloc = 0 to ixend-ixbegin do
      let i = ixbegin + iloc in
      let xfact = data.dx *. float_of_int i in
      let loc = offset + iloc in
      udata.{loc} <- 16.0 *. xfact *. (one -. xfact) *. yfact *. (one -. yfact);
      if i = 0 || i = mx-1 || j = 0 || j = my-1 then iddata.{loc} <- zero;
    done
  done;

  (* Initialize up. *)

  vconst zero up;    (* Initially set up = 0. *)

  (* heatres sets res to negative of ODE RHS values at interior points. *)
  heatres data zero uu up res;

  (* Copy -res into up to get correct initial up values. *)
  vscale (-.one) res up

(*
 * Print first lines of output (problem description)
 * and table heading
 *)

let print_header _Neq rtol atol =
    printf "\nidaHeat2D_kry_bbd_p: Heat equation, parallel example problem for IDA\n";
    printf "                Discretized heat equation on 2D unit square.\n";
    printf "                Zero boundary conditions,";
    printf " polynomial initial conditions.\n";
    printf "                Mesh dimensions: %d x %d" mx my;
    printf "        Total system size: %d\n\n" _Neq;

    printf "Subgrid dimensions: %d x %d" mxsub mysub;
    printf "         Processor array: %d x %d\n" npex npey;
    printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
    printf "Constraints set to force all solution components >= 0. \n";
    printf "SUPPRESSALG = TRUE to suppress local error testing on";
    printf " all boundary components. \n";
    printf "Linear solver: IDASPGMR.    ";
    printf "Preconditioner: IDABBDPRE - Banded-block-diagonal.\n"

(*
 * Print case and table header
 *)

let print_case case_number mudq mukeep =
  printf "\n\nCase %1d. \n" case_number;
  printf "   Difference quotient half-bandwidths = %d" mudq;
  printf "   Retained matrix half-bandwidths = %d \n" mukeep;

  (* Print output table heading and initial line of table. *)
  printf "\n   Output Summary (umax = max-norm of solution) \n\n";
  printf "  time     umax       k  nst  nni  nli   nre nreLS nge     h      npe nps\n";
  printf " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n"

(*
 * Print integrator statistics and max-norm of solution
 *)

let print_output id mem t uu =
  let umax = vmaxnorm uu in

  if id = 0 then begin

    let kused = Ida.get_last_order mem in
    let nst = Ida.get_num_steps mem in
    let nni = Ida.get_num_nonlin_solv_iters mem in
    let nre = Ida.get_num_res_evals mem in
    let hused = Ida.get_last_step mem in
    let nli = Ida.Spils.get_num_lin_iters mem in
    let nreLS = Ida.Spils.get_num_res_evals mem in
    let nge = Ida_bbd.get_num_gfn_evals mem in
    let npe = Ida.Spils.get_num_prec_evals mem in
    let nps = Ida.Spils.get_num_prec_solves mem in
    printf " %5.2f %13.5e  %d  %3d  %3d  %3d  %4d %4d %4d %9.2e  %3d %3d\n"
           t umax kused nst nni nli nre nreLS nge hused npe nps

  end

(*
 * Print some final integrator statistics
 *)

let print_final_stats mem =

  let netf = Ida.get_num_err_test_fails mem in
  let ncfn = Ida.get_num_nonlin_solv_conv_fails mem in
  let ncfl = Ida.Spils.get_num_conv_fails mem in

  printf "\nError test failures            = %d\n" netf;
  printf "Nonlinear convergence failures = %d\n" ncfn;
  printf "Linear convergence failures    = %d\n" ncfl

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* Get processor number and total number of pe's. *)

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

  (* Set local length local_N and global length Neq. *)

  let local_N = mxsub*mysub in
  let neq = mx * my in

  (* Allocate N-vectors. *)

  let uu = Nvector_parallel.make local_N neq comm 0. in

  let up = Nvector_parallel.make local_N neq comm 0. in

  let res = Nvector_parallel.make local_N neq comm 0. in

  let constraints = Nvector_parallel.make local_N neq comm 0. in

  let id = Nvector_parallel.make local_N neq comm 0. in

  (* Allocate and initialize the data structure. *)

  let data = init_user_data thispe comm in

  (* Initialize the uu, up, id, and constraints profiles. *)

  set_initial_profile data (Nvector.unwrap uu) (Nvector.unwrap up)
    (Nvector.unwrap id) (Nvector.unwrap res);
  Nvector_parallel.Ops.n_vconst one constraints;

  let t0 = zero and t1 = 0.01 in

  (* Scalar relative and absolute tolerance. *)

  let rtol = zero in
  let atol = 1.0e-3 in

  (* Call IDACreate and IDAMalloc to initialize solution *)

  (*
   * -----------------------------
   * Case 1 -- mldq = mudq = mxsub
   * -----------------------------
   *)
  (* Call IDASpgmr to specify the linear solver. *)
  (* Call IDABBDPrecInit to initialize BBD preconditioner. *)

  let mudq = mxsub in
  let mldq = mxsub in
  let mukeep = 1 in
  let mlkeep = 1 in

  let linsolv =
    Ida.Spils.spgmr ~maxl:0
      (Ida_bbd.prec_left ~dqrely:zero
         {
           Ida_bbd.mudq = mudq;
           Ida_bbd.mldq = mldq;
           Ida_bbd.mukeep = mukeep;
           Ida_bbd.mlkeep = mlkeep;
         }
         (reslocal data))
  in
  let mem =
    Ida.init linsolv (Ida.SStolerances (rtol,atol))
      (heatres data) t0 uu up
  in
  Ida.set_var_types mem id;
  Ida.set_suppress_alg mem true;
  Ida.set_constraints mem constraints;

  (* Print problem description *)

  if thispe = 0  then print_header neq rtol atol;

  (* Print output heading (on processor 0 only) and initial solution. *)
  if thispe = 0 then print_case 1 mudq mukeep;

  (* Loop over tout, call IDASolve, print output. *)
  let tout = ref t1 in
  for iout = 1 to nout do
    let (tret, _) = Ida.solve_normal mem !tout uu up in

    print_output thispe mem tret uu;

    tout := !tout *. two
  done;

  (* Print final statistics *)
  if thispe = 0 then print_final_stats mem;

  (*
   * -----------------------------
   * Case 2 -- mldq = mudq = 1
   * -----------------------------
   *)

  let mudq = 1 in
  let mldq = 1 in

  (* Re-initialize the uu and up profiles. *)
  set_initial_profile data (Nvector.unwrap uu) (Nvector.unwrap up)
    (Nvector.unwrap id) (Nvector.unwrap res);

  (* Call IDAReInit to re-initialize IDA. *)
  Ida.reinit mem t0 uu up;

  (* Call IDABBDPrecReInit to re-initialize BBD preconditioner. *)
  Ida_bbd.reinit mem mudq mldq ~dqrely:zero;

  (* Print output heading (on processor 0 only). *)
  if thispe = 0 then print_case 2 mudq mukeep;

  (* Loop over tout, call IDASolve, print output. *)
  let tout = ref t1 in
  for iout = 1 to nout do
    let (tret,_) = Ida.solve_normal mem !tout uu up in

    print_output thispe mem tret uu;

    tout := !tout *. two
  done;

  (* Print final statistics *)
  if thispe = 0 then print_final_stats mem

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
