(* {{{
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * OCaml port: T. Bourke, Oct. 2021
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem (based on cvDiurnal_kry_p.c):
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
 * Each species is stored in its own parallel nvector, and are
 * combined together using the MPIManyVector module.
 *
 * The solution is done with the BDF/GMRES method (i.e. using the
 * SUNLinSol_SPGMR linear solver) and the block-diagonal part of the
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
 * Execution: mpirun -np N cvDiurnal_kry_mpimanyvec
 * with N = NPEX*NPEY (see constants below).
 * -----------------------------------------------------------------
 *}}} *)

open Sundials
module DM = Matrix.ArrayDense

let printf = Printf.printf
let eprintf = Printf.eprintf
let unwrap = Nvector_serial.unwrap

(* helpful macros *)

let sqr x = x *. x

(* Problem Constants *)

let nvars      = 2            (* number of species         *)
let kh         = 4.0e-6       (* horizontal diffusivity Kh *)
let vel        = 0.001        (* advection velocity V      *)
let kv0        = 1.0e-8       (* coefficient in Kv(y)      *)
let q1         = 1.63e-16     (* coefficients q1, q2, c3   *)
let q2         = 4.66e-16
let c3g        = 3.7e16
let a3         = 22.62        (* coefficient in expression for q3(t) *)
let a4         = 7.601        (* coefficient in expression for q4(t) *)
let c1_scale   = 1.0e6        (* coefficients in initial profiles    *)
let c2_scale   = 1.0e12

let t0         = 0.0          (* initial time *)
let nout       = 12           (* number of output times *)
let twohr      = 7200.0       (* number of seconds in two hours  *)
let halfday    = 4.32e4       (* number of seconds in a half day *)
let pi         = 3.1415926535898  (* pi *)

let xmin       = 0.0          (* grid boundaries in x  *)
let xmax       = 20.0
let ymin       = 30.0         (* grid boundaries in y  *)
let ymax       = 50.0

let npex       = 2            (* no. PEs in x direction of PE array *)
let npey       = 2            (* no. PEs in y direction of PE array *)
                              (* Total no. PEs = NPEX*NPEY *)
let mxsub      = 5            (* no. x points per subgrid *)
let mysub      = 5            (* no. y points per subgrid *)

let mx         = (npex*mxsub) (* MX = number of x mesh points *)
let my         = (npey*mysub) (* MY = number of y mesh points *)
                              (* Spatial mesh is MX by MY *)

(* CVodeInit Constants *)

let rtol       = 1.0e-5        (* scalar relative tolerance *)
let floor      = 100.0         (* value of C1 or C2 at which tolerances *)
                               (* change from relative to absolute      *)
let atol       = (rtol*.floor) (* scalar absolute tolerance *)

(* User-defined matrix accessor macro: IJth *)

(* IJth is defined in order to write code which indexes into dense
   matrices with a (row,column) pair, where 1 <= row,column <= NVARS.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NVARS. The small matrix routines in sundials_dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. *)

let ijth a i j = DM.get a (j-1) (i-1)
let set_ijth a i j v = DM.set a (j-1) (i-1) v

(* Type : UserData
   contains problem constants, preconditioner blocks, pivot arrays,
   grid constants, and processor indices, as well as data needed
   for the preconditiner *)

type userdata = {
  mutable q4   : float;
  om   : float;
  dx   : float;
  dy   : float;
  hdco : float;
  haco : float;
  vdco : float;

  c1ext : RealArray.t; (* (MXSUB+2)*(MYSUB+2) *)
  c2ext : RealArray.t; (* (MXSUB+2)*(MYSUB+2) *)

  sendbufferE : RealArray.t; (* NVARS*MYSUB *)
  sendbufferW : RealArray.t; (* NVARS*MYSUB *)
  sendbufferN : RealArray.t; (* NVARS*MXSUB *)
  sendbufferS : RealArray.t; (* NVARS*MXSUB *)

  my_pe    : int;
  isubx    : int;
  isuby    : int;
  nvmxsub2 : int;

  comm : Mpi.communicator;
  request : Mpi.request array; (* 8 *)

  (* For preconditioner *)
  p     : DM.t array array; (* [MXSUB][MYSUB] *)
  jbd   : DM.t array array; (* [MXSUB][MXYUB] *)
  pivot : LintArray.t array array; (* [MXSUB][MYSUB] *)
}

(*********************** Private Helper Functions ************************)

(* Load constants in data *)

let init_user_data my_pe comm =
  let dx = (xmax -. xmin) /. (float (mx - 1)) in
  let dy = (ymax -. ymin) /. (float (my - 1)) in
  let isuby = my_pe / npex in
  {
    q4 = 0.0;

    (* Set problem constants *)
    om = pi /. halfday;
    dx;
    dy;
    hdco = kh /. sqr dx;
    haco = vel /. (2.0 *. dx);
    vdco = (1.0 /. sqr dy) *. kv0;

    c1ext = RealArray.make ((mxsub+2)*(mysub+2)) 0.0;
    c2ext = RealArray.make ((mxsub+2)*(mysub+2)) 0.0;

    sendbufferE = RealArray.make (nvars*mysub) 0.0;
    sendbufferW = RealArray.make (nvars*mysub) 0.0;
    sendbufferN = RealArray.make (nvars*mxsub) 0.0;
    sendbufferS = RealArray.make (nvars*mxsub) 0.0;

    (* Set machine-related constants *)
    comm;
    my_pe;
    request = Array.make 8 Mpi.null_request;

    (* isubx and isuby are the PE grid indices corresponding to my_pe *)
    isuby;
    isubx = my_pe - isuby * npex;

    (* Set the sizes of a boundary x-line in u and uext *)
    nvmxsub2 = mxsub + 2;

    (* Preconditioner-related fields *)
    p     = Array.init mxsub
              (fun _ -> Array.init mysub
                  (fun _ -> DM.create nvars nvars));
    jbd   = Array.init mxsub
              (fun _ -> Array.init mysub
                  (fun _ -> DM.create nvars nvars));
    pivot = Array.init mxsub
              (fun _ -> Array.init mysub
                  (fun _ -> LintArray.create nvars));
  }

(* Set initial conditions in u *)

let set_initial_profiles u { isubx; isuby; dx; dy; _ } =
  (* Set pointer to data array in vector u *)
  let udata, _, _ = Nvector.unwrap u in
  let c1data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 0) in
  let c2data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 1) in

  (* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. *)
  let offset = ref 0 in
  let xmid = 0.5 *. (xmin +. xmax) in
  let ymid = 0.5 *. (ymin +. ymax) in
  for ly = 0 to mysub - 1 do
    let jy = float (ly + isuby*mysub) in
    let y = ymin +. jy *. dy in
    let cy = sqr (0.1 *. (y -. ymid)) in
    let cy = 1.0 -. cy +. 0.5 *. sqr cy in
    for lx = 0 to mxsub - 1 do
      let jx = float (lx + isubx * mxsub) in
      let x = xmin +. jx *. dx in
      let cx = sqr (0.1 *. (x -. xmid)) in
      let cx = 1.0 -. cx +. 0.5 *. sqr cx in
      c1data.{!offset} <- c1_scale *. cx *. cy;
      c2data.{!offset} <- c2_scale *. cx *. cy;
      incr offset
    done
  done

(* Print current t, step count, order, stepsize, and sampled c1,c2 values *)

let print_output cvode_mem my_pe comm u t =
  let npelast = npex*npey - 1 in
  let udata, _, _ = Nvector.unwrap u in
  let c1data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 0) in
  let c2data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 1) in
  let tempu = Array.make 2 0.0 in

  (* Send c1,c2 at top right mesh point to PE 0 *)
  if my_pe = npelast then begin
    Array.set tempu 0 (c1data.{mxsub * mysub - 1});
    Array.set tempu 1 (c2data.{mxsub * mysub - 1});
    if npelast != 0 then Mpi.send_float_array tempu 0 0 comm
  end;

  (* On PE 0, receive c1,c2 at top right, then print performance data
     and sampled solution values *)
  if my_pe = 0 then begin
    if npelast != 0 then Mpi.receive_float_array tempu npelast 0 comm;
    let nst = Cvode.get_num_steps cvode_mem in
    let qu = Cvode.get_last_order cvode_mem in
    let hu = Cvode.get_last_step cvode_mem in
    printf "t = %.2e   no. steps = %d   order = %d   stepsize = %.2e\n" t nst qu hu;
    printf "At bottom left:  c1, c2 = %12.3e %12.3e \n" c1data.{0} c2data.{0};
    printf "At top right:    c1, c2 = %12.3e %12.3e \n\n" tempu.(0) tempu.(1);
  end

(* Print final statistics contained in iopt *)

let print_final_stats cvode_mem =
  let lenrw, leniw = Cvode.get_work_space cvode_mem in
  let nst = Cvode.get_num_steps cvode_mem in
  let nfe = Cvode.get_num_rhs_evals cvode_mem in
  let nsetups = Cvode.get_num_lin_solv_setups cvode_mem in
  let netf = Cvode.get_num_err_test_fails cvode_mem in
  let nni = Cvode.get_num_nonlin_solv_iters cvode_mem in
  let ncfn = Cvode.get_num_nonlin_solv_conv_fails cvode_mem in

  let lenrwLS, leniwLS = Cvode.Spils.get_work_space cvode_mem in
  let nli = Cvode.Spils.get_num_lin_iters cvode_mem in
  let npe = Cvode.Spils.get_num_prec_evals cvode_mem in
  let nps = Cvode.Spils.get_num_prec_solves cvode_mem in
  let ncfl = Cvode.Spils.get_num_lin_conv_fails cvode_mem in
  let nfeLS = Cvode.Spils.get_num_lin_rhs_evals cvode_mem in

  printf "\nFinal Statistics: \n\n";
  printf "lenrw   = %5d     leniw   = %5d\n"   lenrw leniw;
  printf "lenrwls = %5d     leniwls = %5d\n"   lenrwLS leniwLS;
  printf "nst     = %5d\n"                     nst;
  printf "nfe     = %5d     nfels   = %5d\n"   nfe nfeLS;
  printf "nni     = %5d     nli     = %5d\n"   nni nli;
  printf "nsetups = %5d     netf    = %5d\n"   nsetups netf;
  printf "npe     = %5d     nps     = %5d\n"   npe nps;
  printf "ncfn    = %5d     ncfl    = %5d\n\n" ncfn ncfl

(* Routine to start receiving boundary data from neighboring PEs. *)

let header_and_empty_array_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 0) []) 0
let float_cell_size =
  Marshal.total_size (Marshal.to_bytes (RealArray.create 1) []) 0
  - header_and_empty_array_size

let bytes x = header_and_empty_array_size + x * float_cell_size

let brecvpost { my_pe; isubx; isuby; comm; request; _ } =
  (* If isuby > 0, receive data for bottom x-line *)
  if isuby > 0 then
    request.(0) <- Mpi.ireceive (bytes (nvars * mxsub)) (my_pe - npex) 0 comm;
    (* RecvBufferS *)

  (* If isuby < NPEY-1, receive data for top x-line *)
  if isuby < npey - 1 then
    request.(1) <- Mpi.ireceive (bytes (nvars * mxsub)) (my_pe + npex) 1 comm;
    (* RecvBufferN *)

  (* If isubx > 0, receive data for left y-line *)
  if isubx > 0 then
    request.(2) <- Mpi.ireceive (bytes (nvars * mysub)) (my_pe - 1) 2 comm;
    (* RecvBufferW *)

  (* If isubx < NPEX-1, receive data for right y-line *)
  if isubx < npex - 1 then
    request.(3) <- Mpi.ireceive (bytes (nvars * mysub)) (my_pe + 1) 3 comm
    (* RecvBufferE *)

(* Routine to send boundary data to neighboring PEs *)

let bsend { my_pe; isubx; isuby; comm; request;
            sendbufferS; sendbufferN; sendbufferE; sendbufferW; _ }
          c1data c2data =
  (* If isuby > 0, send data from bottom x-line of c1 and c2 *)
  if isuby > 0 then begin
    for lx = 0 to mxsub - 1 do
      sendbufferS.{lx} <- c1data.{lx}
    done;
    for lx = 0 to mxsub - 1 do
      sendbufferS.{mxsub + lx} <- c2data.{lx}
    done;
    request.(4) <- Mpi.isend sendbufferS (my_pe - npex) 1 comm
  end;

  (* If isuby < NPEY-1, send data from top x-line of c1 and c2 *)
  if isuby < npey - 1 then begin
    let offsetu = (mysub - 1) * mxsub in
    for lx = 0 to mxsub - 1 do
      sendbufferN.{lx} <- c1data.{offsetu + lx}
    done;
    for lx = 0 to mxsub - 1 do
      sendbufferN.{mxsub + lx} <- c2data.{offsetu + lx}
    done;
    request.(5) <- Mpi.isend sendbufferN (my_pe + npex) 0 comm
  end;

  (* If isubx > 0, send data from left y-line of c1 and c2 *)
  if isubx > 0 then begin
    for ly = 0 to mysub - 1 do
      sendbufferW.{ly} <- c1data.{ly * mxsub}
    done;
    for ly = 0 to mxsub - 1 do
      sendbufferW.{mysub + ly} <- c2data.{ly * mxsub}
    done;
    request.(6) <- Mpi.isend sendbufferW (my_pe - 1) 3 comm
  end;

  (* If isubx < NPEX-1, send data from right y-line of c1 and c2 *)
  if isubx < npex - 1 then begin
    for ly = 0 to mysub - 1 do
      sendbufferE.{ly} <- c1data.{(ly + 1) * mxsub - 1}
    done;
    for ly = 0 to mxsub - 1 do
      sendbufferE.{mysub + ly} <- c2data.{(ly + 1) * mxsub - 1}
    done;
    request.(7) <- Mpi.isend sendbufferE (my_pe + 1) 2 comm
  end

(* Routine to finish receiving boundary data from neighboring PEs. *)

let brecvwait { isuby; isubx; c1ext; c2ext; request; _ } =
  (* If isuby > 0, wait on communication for bottom x-line *)
  if isuby > 0 then begin
    let recvbuffer0 = Mpi.wait_receive request.(0) in
    Mpi.wait request.(4);

    (* Copy the receive buffer to c1ext and c2ext *)
    for lx = 0 to mxsub - 1 do
      c1ext.{1+lx} <- recvbuffer0.{lx}
    done;
    for lx = 0 to mxsub - 1 do
      c2ext.{1+lx} <- recvbuffer0.{mxsub + lx}
    done
  end;

  (* If isuby < NPEY-1, wait on communication for top x-line *)
  if isuby < npey - 1 then begin
    let recvbuffer0 = Mpi.wait_receive request.(1) in
    Mpi.wait request.(5);

    (* Copy the receive buffer to c1ext and c2ext *)
    for lx = 0 to mxsub - 1 do
      c1ext.{(mysub+1)*(mxsub+2)+1+lx} <- recvbuffer0.{lx}
    done;
    for lx = 0 to mxsub - 1 do
      c2ext.{(mysub+1)*(mxsub+2)+1+lx} <- recvbuffer0.{mxsub + lx}
    done
  end;

  (* If isubx > 0, wait on communication for left y-line *)
  if isubx > 0 then begin
    let recvbuffer0 = Mpi.wait_receive request.(2) in
    Mpi.wait request.(6);

    (* Copy the receive buffer to c1ext and c2ext *)
    for ly = 0 to mysub - 1 do
      c1ext.{(ly+1)*(mxsub+2)} <- recvbuffer0.{ly}
    done;
    for ly = 0 to mysub - 1 do
      c2ext.{(ly+1)*(mxsub+2)} <- recvbuffer0.{mysub + ly}
    done
  end;

  (* If isubx < NPEX-1, wait on communication for right y-line *)
  if isubx < npex - 1 then begin
    let recvbuffer0 = Mpi.wait_receive request.(3) in
    Mpi.wait request.(7);

    (* Copy the receive buffer to c1ext and c2ext *)
    for ly = 0 to mysub - 1 do
      c1ext.{(ly+2)*(mxsub+2)-1} <- recvbuffer0.{ly}
    done;
    for ly = 0 to mysub - 1 do
      c2ext.{(ly+2)*(mxsub+2)-1} <- recvbuffer0.{mysub + ly}
    done
  end

(* PrepareExt routine.
   This routine performs all communication between processors of data needed to calculate f.
   It then copies the data from u to the extended work arrays c1ext and c2ext.
   Then to facilitate homogeneous Neumann boundary conditions, this copies data from the
   first interior mesh lines of c1,c2 to c1ext,c2ext for boundary PEs. *)

let prepare_ext (udata, _, _) ({ c1ext; c2ext; isubx; isuby; _ } as data) =
  (* Access data arrays from u, and extended work arrays c1ext and c2ext *)
  let c1data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 0) in
  let c2data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 1) in

  (* Start receiving boundary data from neighboring PEs *)
  brecvpost data;

  (* Send data from boundary of local grid to neighboring PEs *)
  bsend data c1data c2data;

  (* Copy local segments of c1 and c2 vectors into the working
     extended arrays c1ext and c2ext *)
  for ly = 0 to mysub - 1 do
    for lx = 0 to mxsub - 1 do
      c1ext.{(mxsub+2)*(1+ly)+1+lx} <- c1data.{ly*mxsub+lx}
    done
  done;
  for ly = 0 to mysub - 1 do
    for lx = 0 to mxsub - 1 do
      c2ext.{(mxsub+2)*(1+ly)+1+lx} <- c2data.{ly*mxsub+lx}
    done
  done;

  (* Copy data from the first interior mesh line of c1,c2 to
     c1ext, c2ext for boundary PEs *)

  (* If isuby = 0, copy x-line 2 of c1,c2 to c1ext,c2ext *)
  if isuby = 0 then begin
    for lx = 0 to mxsub - 1 do
      c1ext.{1+lx} <- c1data.{mxsub+lx}
    done;
    for lx = 0 to mxsub - 1 do
      c2ext.{1+lx} <- c2data.{mxsub+lx}
    done
  end;

  (* If isuby = NPEY-1, copy x-line MYSUB-1 of c1,c2 to c1ext,c2ext *)
  if isuby = npey - 1 then begin
    for lx = 0 to mxsub - 1 do
      c1ext.{(mysub+1)*(mxsub+2)+1+lx} <- c1data.{(mysub-2)*mxsub+lx}
    done;
    for lx = 0 to mxsub - 1 do
      c2ext.{(mysub+1)*(mxsub+2)+1+lx} <- c2data.{(mysub-2)*mxsub+lx}
    done
  end;

  (* If isubx = 0, copy y-line 2 of u to uext *)
  if isubx = 0 then begin
    for ly = 0 to mysub - 1 do
      c1ext.{(ly+1)*(mxsub+2)} <- c1data.{ly*mxsub+1}
    done;
    for ly = 0 to mysub - 1 do
      c2ext.{(ly+1)*(mxsub+2)} <- c2data.{ly*mxsub+1}
    done
  end;

  (* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext *)
  if isubx = npex - 1 then begin
    for ly = 0 to mysub - 1 do
      c1ext.{(ly+2)*(mxsub+2)-1} <- c1data.{(ly+1)*mxsub-2};
    done;
    for ly = 0 to mysub - 1 do
      c2ext.{(ly+2)*(mxsub+2)-1} <- c2data.{(ly+1)*mxsub-2};
    done
  end;

  (* Finish receiving boundary data from neighboring PEs *)
  brecvwait data

let fcalc t (udot, _, _)
            ({ isuby; hdco; haco; dy; om; vdco; c1ext; c2ext; _ } as data)
  =
  (* Access output arrays *)
  let c1dot, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udot 0) in
  let c2dot, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udot 1) in

  (* Set diurnal rate coefficients as functions of t, and save q4 in
  data block for use by preconditioner evaluation routine *)
  let s = sin (om *. t) in
  let q3, q4coef =
    if s > 0.0
    then exp (-. a3 /. s), exp(-. a4 /. s)
    else 0.0, 0.0
  in
  data.q4 <- q4coef;

  (* Loop over all grid points in local subgrid *)
  for ly = 0 to mysub - 1 do
    (* get global index of y grid coordinate *)
    let jy = ly + isuby * mysub in

    (* Set vertical diffusion coefficients at jy +- 1/2 *)
    let ydn = ymin +. (float jy -. 0.5) *. dy in
    let yup = ydn +. dy in
    let cydn = vdco *. exp(0.2 *. ydn) in
    let cyup = vdco *. exp(0.2 *. yup) in

    (* Loop over x-direction *)
    for lx = 0 to mxsub - 1 do
      (* Extract c1, c2 over 5-point stencil *)
      let offset = (lx+1) + (ly+1) * (mxsub+2) in
      let c1 = c1ext.{offset} in
      let c2 = c2ext.{offset} in
      let c1dn = c1ext.{offset-(mxsub+2)} in
      let c2dn = c2ext.{offset-(mxsub+2)} in
      let c1up = c1ext.{offset+(mxsub+2)} in
      let c2up = c2ext.{offset+(mxsub+2)} in
      let c1lt = c1ext.{offset-1} in
      let c2lt = c2ext.{offset-1} in
      let c1rt = c1ext.{offset+1} in
      let c2rt = c2ext.{offset+1} in

      (* Set kinetic rate terms *)
      let qq1 = q1 *. c1 *. c3g in
      let qq2 = q2 *. c1 *. c2 in
      let qq3 = q3 *. c3g in
      let qq4 = q4coef *. c2 in
      let rkin1 = -. qq1 -. qq2 +. 2.0 *. qq3 +. qq4 in
      let rkin2 = qq1 -. qq2 -. qq4 in

      (* Set vertical diffusion terms *)
      let vertd1 = cyup *. (c1up -. c1) -. cydn *. (c1 -. c1dn) in
      let vertd2 = cyup *. (c2up -. c2) -. cydn *. (c2 -. c2dn) in

      (* Set horizontal diffusion and advection terms *)
      let hord1 = hdco *. (c1rt -. 2.0 *. c1 +. c1lt) in
      let hord2 = hdco *. (c2rt -. 2.0 *. c2 +. c2lt) in
      let horad1 = haco *. (c1rt -. c1lt) in
      let horad2 = haco *. (c2rt -. c2lt) in

      (* Load all terms into dudata *)
      c1dot.{lx+ly*mxsub} <- vertd1 +. hord1 +. horad1 +. rkin1;
      c2dot.{lx+ly*mxsub} <- vertd2 +. hord2 +. horad2 +. rkin2
    done
  done

(***************** Functions Called by the Solver *************************)

(* f routine.  Evaluate f(t,y).  First call PrepareExt to do communication of
   subgrid boundary data and to copy data from u into uext.
   Then calculate f by a call to fcalc. *)

let f user_data t u udot =
  (* Call PrepareExt to set up work arrays for calculation of RHS
     (includes inter-processor communication) *)
  prepare_ext u user_data;

  (* Call fcalc to calculate all right-hand sides *)
  fcalc t udot user_data

(* Preconditioner setup routine. Generate and preprocess P. *)
let precond user_data Cvode.{ jac_y = (udata, _, _); _ } jok gamma =
  let { isuby; dy; vdco; hdco; q4; p; jbd; pivot; _ } = user_data in
  (* Make local copies of pointers in user_data, pointer to u's data,
     and PE index pair *)
  let c1data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 0) in
  let c2data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get udata 1) in

  (* jok = SUNTRUE: Copy Jbd to P, and update jcurPtr *)
  if jok then
    for ly = 0 to mysub - 1 do
      for lx = 0 to mxsub - 1 do
        DM.blit ~src:jbd.(lx).(ly) ~dst:p.(lx).(ly)
      done
    done
  (* jok = SUNFALSE: Generate Jbd from scratch and copy to P, and update jcurPtr *)
  else
    (* Compute 2x2 diagonal Jacobian blocks (using q4 values
       computed on the last f call).  Load into P. *)
    for ly = 0 to mysub - 1 do
      let jy = float (ly + isuby * mysub) in
      let ydn = ymin +. (jy -. 0.5) *. dy in
      let yup = ydn +. dy in
      let cydn = vdco *. exp (0.2 *. ydn) in
      let cyup = vdco *. exp (0.2 *. yup) in
      let diag = -. (cydn +. cyup +. 2.0 *. hdco) in
      for lx = 0 to mxsub - 1 do
        let c1 = c1data.{lx + ly * mxsub} in
        let c2 = c2data.{lx + ly * mxsub} in
        let j = jbd.(lx).(ly) in
        let a = p.(lx).(ly) in
        set_ijth j 1 1 ((-. q1 *. c3g -. q2 *. c2) +. diag);
        set_ijth j 1 2 (-. q2 *. c1 +. q4);
        set_ijth j 2 1 (q1 *. c3g -. q2 *. c2);
        set_ijth j 2 2 ((-. q2 *. c1 -. q4) +. diag);
        DM.blit ~src:j ~dst:a
      done
    done;

  (* Scale by -gamma *)
  for lx = 0 to mxsub - 1 do
    for ly=0 to mysub - 1 do
      DM.scale (-. gamma) p.(lx).(ly);
    done
  done;

  (* Add identity matrix and do LU decompositions on blocks in place *)
  for lx = 0 to mxsub - 1 do
    for ly = 0 to mysub - 1 do
      DM.add_identity p.(lx).(ly);
      DM.getrf p.(lx).(ly) pivot.(lx).(ly)
    done
  done;

  not jok

(* Preconditioner solve routine *)

let psolve { p; pivot; _ }
           _
           Cvode.Spils.{ rhs = r; _ }
           ((zdata, _, _) as z)
  =
  (* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z.
     First copy vector r to z. *)
  Nvector_mpimany.DataOps.scale 1.0 r z;
  let z1data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get zdata 0) in
  let z2data, _, _ = Nvector_parallel.Any.unwrap (ROArray.get zdata 1) in
  for lx = 0 to mxsub - 1 do
    for ly = 0 to mysub - 1 do
      let v = RealArray.of_list [ z1data.{lx + ly*mxsub};
                                  z2data.{lx + ly*mxsub} ] in
      DM.getrs p.(lx).(ly) pivot.(lx).(ly) v
    done
  done

(***************************** Main Program ******************************)

let main () =
  (* Set per-species problem size, neq *)
  let neq = mx*my in

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

  (* Allocate and load user data block; allocate preconditioner block *)
  let data = init_user_data my_pe comm in

  (* Set local length *)
  let local_N = mxsub*mysub in

  (* Allocate c[0], c[1], u, and set initial values and tolerances *)
  let c0 = Nvector_parallel.Any.make local_N neq comm 0.0 in
  let c1 = Nvector_parallel.Any.make local_N neq comm 0.0 in
  let u = Nvector_mpimany.wrap (ROArray.of_list [c0; c1]) in

  set_initial_profiles u data;
  let abstol = atol in
  let reltol = rtol in

  (* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula *)
  (* Set the pointer to user-defined data *)
  (* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. *)
  (* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances *)
  (* Create SPGMR solver structure with left preconditioning
     and the default Krylov dimension maxl *)
  (* Attach SPGMR solver structure to CVode interface *)
  (* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data *)
  let cvode_mem =
    Cvode.(init BDF
      (SStolerances (reltol, abstol))
      ~lsolver:Spils.(solver (spgmr u)
                             (prec_left ~setup:(precond data) (psolve data)))
      (f data) t0 u)
  in

  if my_pe = 0 then
    printf "\n2-species diurnal advection-diffusion problem\n\n";

  print_output cvode_mem my_pe comm u t0;

  (* In loop over output points, call CVode, print results, test for error *)
  let tout = ref twohr in
  for _ = 1 to nout do
    let t, _ = Cvode.solve_normal cvode_mem !tout u in
    print_output cvode_mem my_pe comm u t;
    tout := !tout +. twohr
  done;

  (* Print final statistics *)
  if my_pe = 0 then print_final_stats cvode_mem

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

