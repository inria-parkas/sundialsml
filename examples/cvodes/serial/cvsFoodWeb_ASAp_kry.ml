(*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2011/11/23 23:53:02 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2014.
 * -----------------------------------------------------------------
 * This program solves a stiff ODE system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions. The dependent variable vector
 * is the following:
 *
 *        1   2        ns
 *  c = (c , c , ..., c  )
 *
 * and the PDEs are as follows:
 *
 *    i               i      i
 *  dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)
 *                    xx     yy        i
 *
 * where
 *
 *                 i          ns         j
 *  f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
 *   i                       j=1
 *
 * The number of species is ns = 2*np, with the first np being prey
 * and the last np being predators. The coefficients a(i,j), b(i),
 * d(i) are:
 *
 *  a(i,i) = -a  (all i)
 *  a(i,j) = -g  (i <= np, j > np)
 *  a(i,j) =  e  (i > np, j <= np)
 *  b(i) =  b*(1 + alpha*x*y)  (i <= np)
 *  b(i) = -b*(1 + alpha*x*y)  (i > np)
 *  d(i) = Dprey  (i <= np)
 *  d(i) = Dpred  (i > np)
 *
 * The spatial domain is the unit square. The final time is 10.
 * The boundary conditions are: normal derivative = 0.
 * A polynomial in x and y is used to set the initial conditions.
 *
 * The PDEs are discretized by central differencing on an MX by
 * MY mesh. The resulting ODE system is stiff.
 *
 * The ODE system is solved by CVODES using Newton iteration and
 * the CVSPGMR linear solver (scaled preconditioned GMRES).
 *
 * The preconditioner matrix used is the product of two matrices:
 * (1) A matrix, only defined implicitly, based on a fixed number
 * of Gauss-Seidel iterations using the diffusion terms only.
 * (2) A block-diagonal matrix based on the partial derivatives
 * of the interaction terms f only, using block-grouping (computing
 * only a subset of the ns by ns blocks).
 *
 * Additionally, CVODES can integrate backwards in time the
 * the semi-discrete form of the adjoint PDE:
 *   d(lambda)/dt = - D^T ( lambda_xx + lambda_yy )
 *                  - F_c^T lambda
 * with homogeneous Neumann boundary conditions and final conditions
 *   lambda(x,y,t=t_final) = - g_c^T(t_final)
 * whose solution at t = 0 represents the sensitivity of
 *   int_x int _y g(t_final,c) dx dy dt
 * with respect to the initial conditions of the original problem.
 *
 * In this example,
 *   g(t,c) = c(ISPEC), with ISPEC defined below.
 *
 * -----------------------------------------------------------------
 * Reference:  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage
 * Matrix Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31
 * (1989), pp. 40-91.  Also available as Lawrence Livermore National
 * Laboratory Report UCRL-95088, Rev. 1, June 1987.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module LintArray = Sundials.LintArray
module Adj = Cvodes.Adjoint
module Densemat = Dls.ArrayDenseMatrix
open Bigarray
let unvec = Sundials.unvec

let nvwrmsnorm = Nvector_serial.DataOps.n_vwrmsnorm
let nvlinearsum = Nvector_serial.DataOps.n_vlinearsum

let printf = Printf.printf
let sqr x = x *. x

let zero  = 0.0
let one   = 1.0
let two   = 2.0

(* Problem Specification Constants *)

let aa    = one       (* aa = a *)
let ee    = 1.0e4     (* ee = e *)
let gg    = 0.5e-6    (* gg = g *)
let bb    = one       (* bb = b *)
let dprey = one    
let dpred = 0.5
let alph  = one
let np    = 3
let ns    = (2*np)

let (+>+) (arr : RealArray.t) off = Array1.sub arr off ns

(* Method Constants *)

let mx    = 20
let my    = 20
let mxns  = mx*ns
let ax    = one
let ay    = one
let dx    = ax/.float (mx-1)
let dy    = ay/.float (my-1)
let mp    = ns
let mq    = mx*my
let mxmp  = mx*mp
let ngx   = 2
let ngy   = 2
let ngrp  = ngx*ngy
let itmax = 5

(* CVodeInit Constants *)

let neq   = (ns*mx*my)
let t0    = 0.0
let rtol  = 1.0e-5
let atol  = 1.0e-5

(* Output Constants *)

let tout  = 10.0

(* Note: The value for species i at mesh point (j,k) is stored in *)
(* component number (i-1) + j*NS + k*NS*MX of an N_Vector,        *)
(* where 1 <= i <= NS, 0 <= j < MX, 0 <= k < MY.                  *)

(* Structure for user data *)

type web_data = {
    p         : Densemat.t array;
    pivot     : LintArray.t array;

    ns        : int;
    mxns      : int;

    mp        : int;
    mq        : int;
    mx        : int;
    my        : int;
    ngrp      : int;
    ngx       : int;
    ngy       : int;
    mxmp      : int;

    jgx       : int array;
    jgy       : int array;
    jigx      : int array;
    jigy      : int array;
    jxr       : int array;
    jyr       : int array;

    acoef     : float array array;
    bcoef     : float array;
    diff      : float array;

    cox       : float array;
    coy       : float array;

    dx        : float;
    dy        : float;
    srur      : float;

    fsave     : RealArray.t;
    fbsave    : RealArray.t;

    rewt      : Nvector_serial.t;

    mutable cvode_mem  : Cvode.serial_session option;
    mutable cvode_memb : Adj.serial_bsession option;
  }

(* Adjoint calculation constants *)
(* g = int_x int_y c(ISPEC) dy dx at t = Tfinal *)

let nsteps = 80  (* check points every NSTEPS steps *)
let ispec  =  6  (* species # in objective *)

(*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 *)

(* Small Vector Kernels *)

let v_sum_prods n (u : RealArray.t) uoffs p (q : RealArray.t) qoffs v (w : RealArray.t) woffs =
  for i = 0 to n - 1 do
    u.{uoffs + i} <- p.(i) *. q.{qoffs + i} +. v.(i) *. w.{woffs + i}
  done

let v_inc_by_prod n (u : RealArray.t) uoffs v (w : RealArray.t) woffs =
  for i = 0 to n - 1 do
    u.{uoffs + i} <- u.{uoffs + i} +. v.(i) *. w.{woffs + i}
  done

let v_prod n (u : RealArray.t) uoffs v (w : RealArray.t) woffs =
  for i = 0 to n - 1 do
    u.{uoffs + i} <- v.(i) *. w.{woffs + i}
  done

let v_zero n u offs =
  for i = 0 to n - 1 do
    u.{offs + i} <- zero
  done

(*
 * This routine sets arrays jg, jig, and jr describing
 * a uniform partition of (0,1,2,...,m-1) into ng groups.
 * The arrays set are:
 *   jg    = length ng+1 array of group boundaries.
 *           Group ig has indices j = jg[ig],...,jg[ig+1]-1.
 *   jig   = length m array of group indices vs node index.
 *           Node index j is in group jig[j].
 *   jr    = length ng array of indices representing the groups.
 *           The index for group ig is j = jr[ig].
 *)

let set_groups m ng jg jig jr =
  let mper = m / ng in (* does integer division *)

  for ig = 0 to ng - 1 do
    jg.(ig) <- ig * mper
  done;
  jg.(ng) <- m;

  let ngm1 = ng - 1 in
  let len1 = ngm1 * mper in

  for j = 0 to len1 - 1 do
    jig.(j) <- j / mper
  done;

  for j = len1 to m - 1 do
    jig.(j) <- ngm1
  done;

  for ig = 0 to ngm1 - 1 do
    jr.(ig) <- ((2 * ig + 1) * mper - 1) / 2;
  done;
  jr.(ngm1) <- (ngm1 * mper + m - 1) / 2

(* Allocate space for user data structure *)

let alloc_user_data () =
  let r =
    {
      p          = Array.init ngrp (fun _ -> Densemat.create ns ns);
      pivot      = Array.init ngrp (fun _ -> LintArray.create ns);

      ns         = ns;
      mxns       = mxns;

      mp         = mp;
      mq         = mq;
      mx         = mx;
      my         = my;
      ngrp       = ngrp;
      ngx        = ngx;
      ngy        = ngy;
      mxmp       = mxmp;

      jgx        = Array.make (ngx + 1) 0;
      jgy        = Array.make (ngy + 1) 0;
      jigx       = Array.make mx 0;
      jigy       = Array.make my 0;
      jxr        = Array.make ngx 0;
      jyr        = Array.make ngy 0;

      acoef      = Array.make_matrix ns ns 0.0;
      bcoef      = Array.make ns 0.0;
      diff       = Array.make ns 0.0;

      cox        = Array.make ns 0.0;
      coy        = Array.make ns 0.0;

      dx         = dx;
      dy         = dy;
      srur       = sqrt Sundials.unit_roundoff;

      fsave      = RealArray.create neq;
      fbsave     = RealArray.create neq;

      rewt       = Nvector_serial.make neq 0.0;

      cvode_mem  = None;
      cvode_memb = None;
    }
  in
  let acoef = r.acoef
  and bcoef = r.bcoef
  and diff  = r.diff
  and cox   = r.cox
  and coy   = r.coy
  in
  for j = 0 to np - 1 do
    for i = 0 to np - 1 do
      acoef.(np + i).(j) <- ee;
      acoef.(i).(np + j) <- -. gg
    done;
    acoef.(j).(j)           <- -. aa;
    acoef.(np + j).(np + j) <- -. aa;
    bcoef.(j)               <- bb;
    bcoef.(np + j)          <- -. bb;
    diff.(j)                <- dprey;
    diff.(np + j)           <- dpred
  done;

  for i = 0 to ns - 1 do
    cox.(i) <- diff.(i) /. sqr dx;
    coy.(i) <- diff.(i) /. sqr dy
  done;

  r

(* Initialize user data structure *)

let init_user_data wdata =
  set_groups mx ngx wdata.jgx wdata.jigx wdata.jxr;
  set_groups my ngy wdata.jgy wdata.jigy wdata.jyr

(* This routine computes and loads the vector of initial values. *)

let cinit wdata (cdata : RealArray.t) =
  let ns   = wdata.ns
  and mxns = wdata.mxns
  and dx   = wdata.dx
  and dy   = wdata.dy
  in
  let x_factor = 4.0 /. sqr ax
  and y_factor = 4.0 /. sqr ay
  in

  for jy = 0 to my - 1 do
    let y     = float jy *. dy in
    let argy  = sqr (y_factor *. y *. (ay -. y)) in
    let iyoff = mxns * jy in
    for jx = 0 to mx - 1 do
      let x = float jx *. dx in
      let argx = sqr (x_factor *. x *. (ax -. x)) in
      let ioff = iyoff + ns * jx in
      for i = 1 to ns do
        let ici = ioff + i - 1 in
        cdata.{ici} <- 10.0 +. float i *. argx *. argy
      done
    done
  done

(* This function computes and loads the final values for the adjoint variables *)

let cb_init wdata (cdata : RealArray.t) is =
  let ns   = wdata.ns
  and mxns = wdata.mxns
  in
  let gu = RealArray.make ns zero in
  gu.{ispec-1} <- one;

  for jy = 0 to my-1 do
    let iyoff = mxns*jy in
    for jx = 0 to mx-1 do
      let ioff = iyoff + ns*jx in
      for i = 1 to ns do
        let ici = ioff + i-1 in
        cdata.{ici} <- gu.{i-1}
      done
    done
  done

(*
 * This routine computes the interaction rates for the species
 * c_1, ... ,c_ns (stored in c[0],...,c[ns-1]), at one spatial point 
 * and at time t.
 *)

let web_rates wdata x y ((c : RealArray.t), c_off)
                        ((rate : RealArray.t), rate_off) =
  let acoef = wdata.acoef
  and bcoef = wdata.bcoef
  in
  for i = rate_off to rate_off + ns - 1 do
    rate.{i} <- zero
  done;
  for j = 0 to ns - 1 do
    let c = c.{c_off + j} in
    for i = 0 to ns - 1 do
      rate.{rate_off + i} <- rate.{rate_off + i} +. c *. acoef.(i).(j)
    done
  done;

  let fac = one +. alph *. x *. y in
  for i = 0 to ns - 1 do
    rate.{rate_off + i} <- c.{c_off + i} *. (bcoef.(i) *. fac
                              +. rate.{rate_off + i})
  done

(* This routine computes the interaction rates for the backward problem *)

let web_rates_b wdata x y ((c : RealArray.t), c_off)
                          ((cB : RealArray.t), cB_off)
                          ((rate : RealArray.t), rate_off)
                          ((rateB : RealArray.t), rateB_off) =
  let ns = wdata.ns
  and acoef = wdata.acoef
  and bcoef = wdata.bcoef in

  let fac = one +. alph *. x *. y in

  for i = 0 to ns - 1 do
    rate.{rate_off + i} <- bcoef.(i) *. fac
  done;
  
  for j = 0 to ns - 1 do
    for i = 0 to ns - 1 do
      rate.{rate_off + i} <-
        rate.{rate_off + i} +. acoef.(i).(j) *. c.{c_off + j}
    done
  done;

  for i = 0 to ns - 1 do
    rateB.{rateB_off + i} <- cB.{cB_off + i} *. rate.{rate_off + i};
    rate.{rate_off + i} <- c.{c_off + i} *. rate.{rate_off + i}
  done;

  for j = 0 to ns - 1 do
    for i = 0 to ns - 1 do
      rateB.{rateB_off + i} <- rateB.{rateB_off + i}
                         +. acoef.(j).(i) *. c.{c_off + j} *. cB.{cB_off + j}
    done
  done

(*
 * This routine computes one block of the interaction terms of the
 * system, namely block (jx,jy), for use in preconditioning.
 * Here jx and jy count from 0.
 *)

let fblock wdata t cdata jx jy cdotdata =
  let iblok = jx + jy * wdata.mx
  and y = float jy *. wdata.dy
  and x = float jx *. wdata.dx
  in
  let ic = wdata.ns * iblok in
  web_rates wdata x y (cdata, ic) (cdotdata, 0)

(*
 * This routine performs ITMAX=5 Gauss-Seidel iterations to compute an
 * approximation to (P-inverse)*z, where P = I - gamma*Jd, and
 * Jd represents the diffusion contributions to the Jacobian.
 * The answer is stored in z on return, and x is a temporary vector.
 * The dimensions below assume a global constant NS >= ns.
 * Some inner loops of length ns are implemented with the small
 * vector kernels v_sum_prods, v_prod, v_inc_by_prod.
 *)

let gs_iter wdata gamma zd xd =
  let ns = wdata.ns
  and mx = wdata.mx
  and my = wdata.my
  and mxns = wdata.mxns
  and cox = wdata.cox
  and coy = wdata.coy
  in

  let beta  = Array.make ns 0.0
  and beta2 = Array.make ns 0.0
  and cof1  = Array.make ns 0.0
  and gam   = Array.make ns 0.0
  and gam2  = Array.make ns 0.0
  in

  (* Write matrix as P = D - L - U.
     Load local arrays beta, beta2, gam, gam2, and cof1. *)
  
  for i = 0 to ns - 1 do
    let temp = one /. (one +. two *. gamma *. (cox.(i) +. coy.(i))) in
    beta.(i)  <- gamma *. cox.(i) *. temp;
    beta2.(i) <- two *. beta.(i);
    gam.(i)   <- gamma *. coy.(i) *. temp;
    gam2.(i)  <- two *. gam.(i);
    cof1.(i)  <- temp
  done;
  
  (* Begin iteration loop.
     Load vector x with (D-inverse)*z for first iteration. *)
  for jy = 0 to my - 1 do
    let iyoff = mxns * jy in
    for jx = 0 to mx - 1 do
      let ic = iyoff + ns*jx in
      v_prod ns xd ic cof1 zd ic (* x[ic+i] = cof1[i]z[ic+i] *)
    done
  done;
  Array1.fill zd zero;
  
  (* Looping point for iterations. *)
  
  for iter = 1 to itmax do
    
    (* Calculate (D-inverse)*U*x if not the first iteration. *)
    
    if (iter > 1) then
      for jy = 0 to my - 1 do
        let iyoff = mxns * jy in
        for jx = 0 to mx - 1 do (* order of loops matters *)
          let ic = iyoff + ns * jx
          and x_loc = if jx = 0 then 0 else if jx = mx - 1 then 2 else 1
          and y_loc = if jy = 0 then 0 else if jy = my - 1 then 2 else 1
          in

          match (3 * y_loc + x_loc) with
          | 0 ->
            (* jx == 0, jy == 0 *)
            (* x[ic+i] = beta2[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] *)
            v_sum_prods ns xd ic beta2 xd (ic + ns) gam2 xd (ic + mxns)

          | 1 ->
            (* 1 <= jx <= mx-2, jy == 0 *)
            (* x[ic+i] = beta[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] *)
            v_sum_prods ns xd ic beta xd (ic + ns) gam2 xd (ic + mxns)

          | 2 ->
            (* jx == mx-1, jy == 0 *)
            (* x[ic+i] = gam2[i]x[ic+mxns+i] *)
            v_prod ns xd ic gam2 xd (ic + mxns)

          | 3 ->
            (* jx == 0, 1 <= jy <= my-2 *)
            (* x[ic+i] = beta2[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] *)
            v_sum_prods ns xd ic beta2 xd (ic + ns) gam xd (ic + mxns)

          | 4 ->
            (* 1 <= jx <= mx-2, 1 <= jy <= my-2 *)
            (* x[ic+i] = beta[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] *)
            v_sum_prods ns xd ic beta xd (ic + ns) gam xd (ic + mxns)

          | 5 ->
            (* jx == mx-1, 1 <= jy <= my-2 *)
            (* x[ic+i] = gam[i]x[ic+mxns+i] *)
            v_prod ns xd ic gam xd (ic + mxns)

          | 6 ->
            (* jx == 0, jy == my-1 *)
            (* x[ic+i] = beta2[i]x[ic+ns+i] *)
            v_prod ns xd ic beta2 xd (ic + ns)

          | 7 ->
            (* 1 <= jx <= mx-2, jy == my-1 *)
            (* x[ic+i] = beta[i]x[ic+ns+i] *)
            v_prod ns xd ic beta xd (ic + ns)

          | 8 ->
            (* jx == mx-1, jy == my-1 *)
            (* x[ic+i] = 0.0 *)
            v_zero ns xd ic

          | _ -> assert false
        done
      done;  (* end if (iter > 1) *)
    
    (* Overwrite x with [(I - (D-inverse)*L)-inverse]*x. *)
    
    for jy = 0 to my - 1 do
      let iyoff = mxns * jy in
      for jx = 0 to mx - 1 do (* order of loops matters *)
        let ic = iyoff + ns * jx
        and x_loc = if jx = 0 then 0 else if jx = mx - 1 then 2 else 1
        and y_loc = if jy = 0 then 0 else if jy = my - 1 then 2 else 1
        in
        match (3 * y_loc + x_loc) with
        | 0 ->
          (* jx == 0, jy == 0 *)
            ()

        | 1 ->
          (* 1 <= jx <= mx-2, jy == 0 *)
          (* x[ic+i] += beta[i]x[ic-ns+i] *)
          v_inc_by_prod ns xd ic beta xd (ic - ns)

        | 2 ->
          (* jx == mx-1, jy == 0 *)
          (* x[ic+i] += beta2[i]x[ic-ns+i] *)
          v_inc_by_prod ns xd ic beta2 xd (ic - ns)

        | 3 ->
          (* jx == 0, 1 <= jy <= my-2 *)
          (* x[ic+i] += gam[i]x[ic-mxns+i] *)
          v_inc_by_prod ns xd ic gam xd (ic - mxns)

        | 4 ->
          (* 1 <= jx <= mx-2, 1 <= jy <= my-2 *)
          (* x[ic+i] += beta[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] *)
          v_inc_by_prod ns xd ic beta xd (ic - ns);
          v_inc_by_prod ns xd ic gam xd (ic - mxns)

        | 5 ->
          (* jx == mx-1, 1 <= jy <= my-2 *)
          (* x[ic+i] += beta2[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] *)
          v_inc_by_prod ns xd ic beta2 xd (ic - ns);
          v_inc_by_prod ns xd ic gam xd (ic - mxns)

        | 6 ->
          (* jx == 0, jy == my-1 *)
          (* x[ic+i] += gam2[i]x[ic-mxns+i] *)
          v_inc_by_prod ns xd ic gam2 xd (ic - mxns)

        | 7 ->
          (* 1 <= jx <= mx-2, jy == my-1 *)
          (* x[ic+i] += beta[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] *)
          v_inc_by_prod ns xd ic beta xd (ic - ns);
          v_inc_by_prod ns xd ic gam2 xd (ic - mxns)

        | 8 ->
          (* jx == mx-1, jy == my-1 *)
          (* x[ic+i] += beta2[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] *)
          v_inc_by_prod ns xd ic beta2 xd (ic - ns);
          v_inc_by_prod ns xd ic gam2 xd (ic - mxns)

        | _ -> assert false
      done
    done;
    
    (* Add increment x to z : z <- z+x *)
    nvlinearsum one zd one xd zd
  done

(* Print maximum sensitivity of G for each species *)

let print_output wdata cdata ns mxns =
  let x = ref zero
  and y = ref zero in

  for i=1 to ns do
    let cmax = ref zero in
    for jy=my-1 downto 0 do
      for jx=0 to mx - 1 do
        let cij = cdata.{(i-1) + jx*ns + jy*mxns} in
        if abs_float(cij) > !cmax then begin
          cmax := cij;
          x := float jx *. wdata.dx;
          y := float jy *. wdata.dy
        end
      done
    done;

    printf "\nMaximum sensitivity with respect to I.C. of species %d\n" i;
    printf "  mu max = %e\n" !cmax;
    printf "at\n";
    printf "  x = %e\n  y = %e\n" !x !y
  done

(* Compute double space integral *)

let double_intgr (cdata : RealArray.t) i wdata =
  let ns   = wdata.ns
  and mx   = wdata.mx
  and my   = wdata.my
  and mxns = wdata.mxns
  and dx   = wdata.dx
  and dy   = wdata.dy in

  let jy = 0 in
  let intgr_x = ref cdata.{(i-1) + jy*mxns} in
  for jx = 1 to mx-2 do
    intgr_x := !intgr_x +. two *. cdata.{(i-1) + jx*ns + jy*mxns}
  done;
  intgr_x := !intgr_x +. cdata.{(i-1) + (mx-1)*ns + jy*mxns};
  intgr_x := !intgr_x *. 0.5 *. dx;
  
  let intgr_xy = ref !intgr_x in
  for jy = 1 to my-2 do
    intgr_x := cdata.{(i-1)+jy*mxns};
    for jx = 1 to mx-2 do
      intgr_x := !intgr_x +. two *. cdata.{(i-1) + jx*ns + jy*mxns}
    done;
    intgr_x := !intgr_x +. cdata.{(i-1) + (mx-1)*ns + jy*mxns};
    intgr_x := !intgr_x *. 0.5 *. dx;
    intgr_xy := !intgr_xy +. two *. !intgr_x
  done;
  
  let jy = my-1 in
  intgr_x := cdata.{(i-1) + jy*mxns};
  for jx = 1 to mx-2 do
    intgr_x := !intgr_x +. two *. cdata.{(i-1) + jx*ns + jy*mxns}
  done;
  intgr_x := !intgr_x +. cdata.{(i-1) + (mx-1)*ns + jy*mxns};
  intgr_x := !intgr_x *. 0.5 *. dx;
  
  intgr_xy := !intgr_xy +. !intgr_x;
  intgr_xy := !intgr_xy *. 0.5 *. dy;

  !intgr_xy

(*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 *)

(*
 * This routine computes the right-hand side of the ODE system and
 * returns it in cdot. The interaction rates are computed by calls to WebRates,
 * and these are saved in fsave for use in preconditioning.
 *)

let f wdata t cdata (cdotdata : RealArray.t) =
  let ns    = wdata.ns
  and fsave = wdata.fsave
  and cox   = wdata.cox
  and coy   = wdata.coy
  and mxns  = wdata.mxns
  and dx    = wdata.dx
  and dy    = wdata.dy
  in
   
  for jy = 0 to my - 1 do
    let y = float jy *. dy in
    let iyoff = mxns*jy in
    let idyu = if jy = my - 1 then - mxns else mxns in
    let idyl = if jy = 0      then - mxns else mxns in
    for jx = 0 to mx - 1 do
      let x = float jx *. dx in
      let ic = iyoff + ns * jx in
      (* Get interaction rates at one point (x,y). *)
      web_rates wdata x y (cdata, ic) (fsave, ic);
      let idxu = if jx = mx - 1 then -ns else ns in
      let idxl = if jx = 0      then -ns else ns in
      for i = 1 to ns do
        let ici = ic + i - 1 in
        (* Do differencing in y. *)
        let dcyli = cdata.{ici} -. cdata.{ici - idyl} in
        let dcyui = cdata.{ici + idyu} -. cdata.{ici} in
        (* Do differencing in x. *)
        let dcxli = cdata.{ici} -. cdata.{ici - idxl} in
        let dcxui = cdata.{ici + idxu} -. cdata.{ici} in
        (* Collect terms and load cdot elements. *)
        cdotdata.{ici} <- coy.(i - 1) *. (dcyui -. dcyli)
                          +. cox.(i - 1) *. (dcxui -. dcxli)
                          +. fsave.{ici}
      done
    done
  done

(*
 * This routine generates the block-diagonal part of the Jacobian
 * corresponding to the interaction rates, multiplies by -gamma, adds
 * the identity matrix, and calls denseGETRF to do the LU decomposition of
 * each diagonal block. The computation of the diagonal blocks uses
 * the preset block and grouping information. One block per group is
 * computed. The Jacobian elements are generated by difference
 * quotients using calls to the routine fblock.
 *
 * This routine can be regarded as a prototype for the general case
 * of a block-diagonal preconditioner. The blocks are of size mp, and
 * there are ngrp=ngx*ngy blocks computed in the block-grouping scheme.
 *)
 
let precond wdata jacarg jok gamma =
  let { Cvode.jac_t   = t;
        Cvode.jac_y   = cdata;
        Cvode.jac_fy  = fc;
        Cvode.jac_tmp = (vtemp1, _, _)
      } = jacarg
  in
  let f1 = vtemp1 in
  let cvode_mem =
    match wdata.cvode_mem with
    | Some c -> c | None -> assert false
  and rewtdata  = unvec wdata.rewt
  in
  Cvode.get_err_weights cvode_mem wdata.rewt;

  let uround = Sundials.unit_roundoff
  and p      = wdata.p
  and pivot  = wdata.pivot
  and jxr    = wdata.jxr
  and jyr    = wdata.jyr
  and mp     = wdata.mp
  and srur   = wdata.srur
  and ngx    = wdata.ngx
  and ngy    = wdata.ngy
  and mxmp   = wdata.mxmp
  and fsave  = wdata.fsave
  in
  (* Make mp calls to fblock to approximate each diagonal block of Jacobian.
     Here, fsave contains the base value of the rate vector and 
     r0 is a minimum increment factor for the difference quotient. *)

  let fac = nvwrmsnorm fc rewtdata in
  let r0 = 1000.0 *. abs_float gamma *. uround *. float neq *. fac in
  let r0 = if r0 = zero then one else r0 in
  
  for igy = 0 to ngy - 1 do
    let jy = jyr.(igy) in
    let if00 = jy * mxmp in
    for igx = 0 to ngx - 1 do
      let jx  = jxr.(igx) in
      let if0 = if00 + jx * mp in
      let ig  = igx + igy * ngx in
      (* Generate ig-th diagonal block *)
      let pdata = Sundials.RealArray2.unwrap p.(ig) in
      for j = 0 to mp - 1 do
        (* Generate the jth column as a difference quotient *)
        let jj = if0 + j in
        let save = cdata.{jj} in
        let r = max (srur *. abs_float save) (r0 /. rewtdata.{jj}) in
        cdata.{jj} <- cdata.{jj} +. r;
        fblock wdata t cdata jx jy f1;
        let fac = -. gamma /. r in
        for i = 0 to mp - 1 do
          pdata.{j, i} <- (f1.{i} -. fsave.{if0 + i}) *. fac
        done;
        cdata.{jj} <- save
      done
    done
  done;
  
  (* Add identity matrix and do LU decompositions on blocks. *)
  let f ig p_ig =
    Densemat.add_identity p_ig;
    Densemat.getrf p_ig pivot.(ig)
  in
  Array.iteri f p;
  true

(*
 * This routine applies two inverse preconditioner matrices
 * to the vector r, using the interaction-only block-diagonal Jacobian
 * with block-grouping, denoted Jr, and Gauss-Seidel applied to the
 * diffusion contribution to the Jacobian, denoted Jd.
 * It first calls GSIter for a Gauss-Seidel approximation to
 * ((I - gamma*Jd)-inverse)*r, and stores the result in z.
 * Then it computes ((I - gamma*Jr)-inverse)*z, using LU factors of the
 * blocks in P, and pivot information in pivot, and returns the result in z.
 *)

let psolve wdata =
  let cache = RealArray.create ns in
  fun jac_arg solve_arg z ->
  let { Cvode.jac_tmp = vtemp; } = jac_arg
  and { Cvode.Spils.rhs = r; Cvode.Spils.gamma = gamma } = solve_arg
  in
  Array1.blit r z;

  (* call GSIter for Gauss-Seidel iterations *)
  gs_iter wdata gamma z vtemp;
  
  (* Do backsolves for inverse of block-diagonal preconditioner factor *)
  let p     = wdata.p
  and pivot = wdata.pivot
  and mx    = wdata.mx
  and my    = wdata.my
  and ngx   = wdata.ngx
  and mp    = wdata.mp
  and jigx  = wdata.jigx
  and jigy  = wdata.jigy
  in
  
  let iv = ref 0 in
  for jy = 0 to my - 1 do
    let igy = jigy.(jy) in
    for jx = 0 to mx - 1 do
      let igx = jigx.(jx) in
      let ig = igx + igy * ngx in

      (* faster to cache and copy in/out than to Bigarray.Array1.sub... *)
      for i=0 to ns - 1 do
        cache.{i} <- z.{!iv + i}
      done;

      Densemat.getrs p.(ig) pivot.(ig) cache;

      for i=0 to ns - 1 do
        z.{!iv + i} <- cache.{i}
      done;
      iv := !iv + mp
    done
  done

(*
 * This routine computes the right-hand side of the adjoint ODE system and
 * returns it in cBdot. The interaction rates are computed by calls to WebRates,
 * and these are saved in fsave for use in preconditioning. The adjoint 
 * interaction rates are computed by calls to WebRatesB.
 *)

let fB wdata t (cdata : RealArray.t) (cBdata : RealArray.t) (cBdotdata : RealArray.t) =
  let mxns   = wdata.mxns
  and ns     = wdata.ns
  and fsave  = wdata.fsave
  and fBsave = wdata.fbsave
  and cox    = wdata.cox
  and coy    = wdata.coy
  and dx     = wdata.dx
  and dy     = wdata.dy in

  for jy = 0 to my - 1 do
    let y = float jy *. dy
    and iyoff = mxns*jy
    and idyu  = if jy = my-1 then -mxns else mxns
    and idyl  = if jy = 0 then -mxns else mxns
    in
    for jx = 0 to mx - 1 do
      let x  = float jx *. dx
      and ic = iyoff + ns*jx in
      (* Get interaction rates at one point (x,y). *)
      web_rates_b wdata x y (cdata, ic) (cBdata, ic) (fsave, ic) (fBsave, ic);
      let idxu = if jx = mx-1 then -ns else ns
      and idxl = if jx = 0 then -ns else ns
      in
      for i = 1 to ns do
        let ici = ic + i-1 in
        (* Do differencing in y. *)
        let dcyli = cBdata.{ici} -. cBdata.{ici - idyl} in
        let dcyui = cBdata.{ici + idyu} -. cBdata.{ici} in
        (* Do differencing in x. *)
        let dcxli = cBdata.{ici} -. cBdata.{ici - idxl} in
        let dcxui = cBdata.{ici + idxu} -. cBdata.{ici} in
        (* Collect terms and load cdot elements. *)
        cBdotdata.{ici} <- -. coy.(i-1)*.(dcyui -. dcyli) 
                           -. cox.(i-1)*.(dcxui -. dcxli)
                           -. fBsave.{ici}
      done
    done
  done

(* Preconditioner setup function for the backward problem *)

let precondb wdata jacarg jok gamma =
  let { Adj.jac_t   = t;
        Adj.jac_y   = cdata;
        Adj.jac_yb  = cBdata;
        Adj.jac_fyb = fcBdata;
        Adj.jac_tmp = (vtemp1, _, _)
      } = jacarg
  in
  let f1 = vtemp1 in
  let cvode_mem =
    match wdata.cvode_memb with
    | Some c -> c | None -> assert false
  and rewtdata = unvec wdata.rewt
  in
  Adj.get_err_weights cvode_mem wdata.rewt;

  let uround = Sundials.unit_roundoff
  and p      = wdata.p
  and pivot  = wdata.pivot
  and jxr    = wdata.jxr
  and jyr    = wdata.jyr
  and mp     = wdata.mp
  and srur   = wdata.srur
  and ngx    = wdata.ngx
  and ngy    = wdata.ngy
  and mxmp   = wdata.mxmp
  and fsave  = wdata.fsave
  in
  (* Make mp calls to fblock to approximate each diagonal block of Jacobian.
     Here, fsave contains the base value of the rate vector and 
     r0 is a minimum increment factor for the difference quotient. *)

  let fac = nvwrmsnorm fcBdata rewtdata in
  let r0 = 1000.0 *. abs_float gamma *. uround *. float neq *. fac in
  let r0 = if r0 = zero then one else r0 in

  for igy = 0 to ngy - 1 do
    let jy = jyr.(igy) in
    let if00 = jy * mxmp in
    for igx = 0 to ngx - 1 do
      let jx  = jxr.(igx) in
      let if0 = if00 + jx * mp in
      let ig  = igx + igy * ngx in
      (* Generate ig-th diagonal block *)
      let pdata = Sundials.RealArray2.unwrap p.(ig) in
      for j = 0 to mp - 1 do
        (* Generate the jth column as a difference quotient *)
        let jj = if0 + j in
        let save = cdata.{jj} in
        let r = max (srur *. abs_float save) (r0 /. rewtdata.{jj}) in
        cdata.{jj} <- cdata.{jj} +. r;
        fblock wdata t cdata jx jy f1;
        let fac = gamma /. r in
        for i = 0 to mp - 1 do
          pdata.{i, j} <- (f1.{i} -. fsave.{if0 + i}) *. fac
        done;
        cdata.{jj} <- save
      done
    done
  done;

  (* Add identity matrix and do LU decompositions on blocks. *)
  let f ig p_ig =
    Densemat.add_identity p_ig;
    Densemat.getrf p_ig pivot.(ig)
  in
  Array.iteri f p;
  true

(* Preconditioner solve function for the backward problem *)

let psolveb wdata =
  let cache = RealArray.create ns in
  fun jac_arg solve_arg (z : RealArray.t) ->
  let { Adj.jac_tmp = vtemp; } = jac_arg
  and { Adj.Spils.rvec = r; Adj.Spils.gamma = gamma } = solve_arg
  in
  Array1.blit r z;

  (* call GSIter for Gauss-Seidel iterations
     (same routine but with gamma=-gamma) *)
  gs_iter wdata (-.gamma) z vtemp;
  
  (* Do backsolves for inverse of block-diagonal preconditioner factor *)
  let p     = wdata.p
  and pivot = wdata.pivot
  and mx    = wdata.mx
  and my    = wdata.my
  and ngx   = wdata.ngx
  and mp    = wdata.mp
  and jigx  = wdata.jigx
  and jigy  = wdata.jigy
  in
  
  let iv = ref 0 in
  for jy = 0 to my - 1 do
    let igy = jigy.(jy) in
    for jx = 0 to mx - 1 do
      let igx = jigx.(jx) in
      let ig = igx + igy * ngx in

      (* faster to cache and copy in/out than to Bigarray.Array1.sub... *)
      for i=0 to ns - 1 do
        cache.{i} <- z.{!iv + i}
      done;

      Densemat.getrs p.(ig) pivot.(ig) cache;

      for i=0 to ns - 1 do
        z.{!iv + i} <- cache.{i}
      done;
      iv := !iv + mp
    done
  done

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  let abstol,  reltol = atol, rtol
  and abstolb, reltolb = atol, rtol
  in
  (* Allocate and initialize user data *)
  let wdata = alloc_user_data () in
  init_user_data wdata;

  (* Set-up forward problem *)
  (* Initializations *)
  let c = Nvector_serial.make neq 0.0 in
  cinit wdata (unvec c);

  (* Call CVodeCreate/CVodeInit for forward run *)
  (* Call CVSpgmr for forward run *)
  printf "\nCreate and allocate CVODES memory for forward run\n";

  let cvode_mem =
    Cvode.init
        Cvode.BDF
        (Cvode.Newton
            (Cvode.Spils.spgmr None Spils.PrecLeft
                { Cvode.Spils.prec_setup_fn = Some (precond wdata);
                  Cvode.Spils.prec_solve_fn = Some (psolve wdata);
                  Cvode.Spils.jac_times_vec_fn = None }))
        (Cvode.SStolerances (reltol, abstol))
        (f wdata) ~t0:t0 c
  in
  wdata.cvode_mem <- Some cvode_mem; (* Used in Precond *)

  (* Set-up adjoint calculations *)
  printf "\nAllocate global memory\n";
  Adj.init cvode_mem nsteps Adj.IHermite;

  (* Perform forward run *)
  printf "\nForward integration\n";
  let t, ncheck, _ = Adj.forward_normal cvode_mem tout c in

  printf "\nncheck = %d\n"  ncheck;

  printf "\n   g = int_x int_y c%d(Tfinal,x,y) dx dy = %f \n\n"
         ispec (double_intgr (unvec c) ispec wdata);

  (* Set-up backward problem *)

  (* Allocate cB *)
  (* Initialize cB = 0 *)
  let cB = Nvector_serial.make neq zero in
  cb_init wdata (unvec cB) ispec;

  (* Create and allocate CVODES memory for backward run *)
  (* Call CVSpgmr *)
  printf "\nCreate and allocate CVODES memory for backward run\n";

  let cvode_memb =
    Adj.init_backward
      cvode_mem
      Cvode.BDF
      (Adj.Newton (Adj.Spils.spgmr None Spils.PrecLeft
          { Adj.Spils.prec_setup_fn = Some (precondb wdata);
            Adj.Spils.prec_solve_fn = Some (psolveb wdata);
            Adj.Spils.jac_times_vec_fn = None }))
      (Adj.SStolerances (reltolb, abstolb))
      (Adj.Basic (fB wdata))
      tout
      cB
  in
  wdata.cvode_memb <- Some cvode_memb;

  (* Perform backward integration *)

  printf "\nBackward integration\n";
  Adj.backward_normal cvode_mem t0;
  let _ = Adj.get cvode_memb cB in
  print_output wdata (unvec cB) ns mxns


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
