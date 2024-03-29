(*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2009/02/17 02:48:46 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Feb 2011.
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
 * 10 x 10 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODE, with the BDF/GMRES
 * method (i.e. using the CVSPGMR linear solver) and the
 * block-diagonal part of the Newton matrix as a left
 * preconditioner. A copy of the block-diagonal part of the
 * Jacobian is saved and conditionally reused within the Precond
 * routine.
 * -----------------------------------------------------------------
 *)

open Sundials

module Direct = Matrix.ArrayDense
let unwrap = Nvector.unwrap

let printf = Printf.printf

(* Problem Constants *)

let zero  = 0.0
let one   = 1.0
let two   = 2.0

let num_species   = 2           (* number of species         *)
let kh            = 4.0e-6      (* horizontal diffusivity Kh *)
let vel           = 0.001       (* advection velocity V      *)
let kv0           = 1.0e-8      (* coefficient in Kv(y)      *)
let q1            = 1.63e-16    (* coefficients q1, q2, c3   *)
let q2            = 4.66e-16
let c3            = 3.7e16
let a3            = 22.62       (* coefficient in expression for q3(t) *)
let a4            = 7.601       (* coefficient in expression for q4(t) *)
let c1_scale      = 1.0e6       (* coefficients in initial profiles    *)
let c2_scale      = 1.0e12

let t0            = zero        (* initial time *)
let nout          = 12          (* number of output times *)
let twohr         = 7200.0      (* number of seconds in two hours  *)
let halfday       = 4.32e4      (* number of seconds in a half day *)
let pi            = 3.1415926535898 (* pi *)

let xmin          = zero        (* grid boundaries in x  *)
let xmax          = 20.0
let ymin          = 30.0        (* grid boundaries in y  *)
let ymax          = 50.0
let xmid          = 10.0        (* grid midpoints in x,y *)
let ymid          = 40.0

let mx            = 10          (* mx = number of x mesh points *)
let my            = 10          (* my = number of y mesh points *)
let nsmx          = 20          (* nsmx = num_species*mx *)
let mm            = (mx * my)   (* mm = mx*my *)

(* CVodeInit Constants *)

let rtol     = 1.0e-5           (* scalar relative tolerance *)
let floor    = 100.0            (* value of C1 or C2 at which tolerances *)
                                (* change from relative to absolute      *)
let atol     = (rtol *. floor)  (* scalar absolute tolerance *)
let neq      = (num_species * mm) (* neq = number of equations *)

(* User-defined vector and matrix accessor macros: IJKth, IJth *)

(* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.

   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MY-1. The vdata array is obtained via
   the macro call vdata = NV_DATA_S(v), where v is an N_Vector.
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in sundials_dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. *)

let ijkth (v : RealArray.t) i j k
  = v.{i - 1 + j * num_species + k * nsmx}
let set_ijkth (v : RealArray.t) i j k e
  = v.{i - 1 + j * num_species + k * nsmx} <- e

let ijth v i j       = Direct.get v (i - 1) (j - 1)
let set_ijth v i j e = Direct.set v (i - 1) (j - 1) e

(* Type : UserData
   contains preconditioner blocks, pivot arrays, and problem constants *)

type user_data = {
        p          : Direct.t array array;
        jbd        : Direct.t array array;
        pivot      : LintArray.t array array;
        mutable q4 : float;
        om         : float;
        dx         : float;
        dy         : float;
        hdco       : float;
        haco       : float;
        vdco       : float;
    }

(* Private Helper Functions *)

let sqr x = x *. x

(* Allocate memory for data structure of type UserData *)

let alloc_user_data () =
  let new_dmat _ = Direct.create num_species num_species in
  let new_int1 _  = LintArray.create num_species in
  let new_y_arr elinit _ = Array.init my elinit in
  let new_xy_arr elinit  = Array.init mx (new_y_arr elinit) in
  {
    p     = new_xy_arr new_dmat;
    jbd   = new_xy_arr new_dmat;
    pivot = new_xy_arr new_int1;
    q4    = 0.0;
    om    = 0.0;
    dx    = 0.0;
    dy    = 0.0;
    hdco  = 0.0;
    haco  = 0.0;
    vdco  = 0.0;
  }

(* Load problem constants in data *)

let init_user_data data =
    let om = pi /. halfday
    and dx = (xmax -. xmin) /. float (mx - 1)
    and dy = (ymax -. ymin) /. float (my - 1)
    in
    { data with
        om = om;
        dx = dx;
        dy = dy;
        hdco = kh /. sqr(dx);
        haco = vel /. (two *. dx);
        vdco = (one /. sqr(dy)) *. kv0;
    }

(* Set initial conditions in u *)

let set_initial_profiles udata dx dy =
  for jy = 0 to my - 1 do
    let y = ymin +. float(jy) *. dy in
    let cy = sqr (0.1 *. (y -. ymid)) in
    let cy = one -. cy +. 0.5 *. sqr(cy) in

    for jx = 0 to mx - 1 do
      let x = xmin +. float jx *. dx in
      let cx = sqr(0.1 *. (x -. xmid)) in
      let cx = one -. cx +. 0.5 *. sqr(cx) in

      set_ijkth udata 1 jx jy (c1_scale *. cx *. cy);
      set_ijkth udata 2 jx jy (c2_scale *. cx *. cy)
    done
  done

(* Print current t, step count, order, stepsize, and sampled c1,c2 values *)

let print_output s udata t =
  let mxh = mx / 2 - 1
  and myh = my / 2 - 1
  and mx1 = mx - 1
  and my1 = my - 1
  in
  let nst = Cvode.get_num_steps s
  and qu  = Cvode.get_last_order s
  and hu  = Cvode.get_last_step s
  in
  printf "t = %.2e   no. steps = %d   order = %d   stepsize = %.2e\n"
         t nst qu hu;
  printf "c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n"
         (ijkth udata 1 0 0)
         (ijkth udata 1 mxh myh)
         (ijkth udata 1 mx1 my1);
  printf "c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n\n"
         (ijkth udata 2 0 0)
         (ijkth udata 2 mxh myh)
         (ijkth udata 2 mx1 my1)

(* Get and print final statistics *)

let print_final_stats s =
  let open Cvode in
  let lenrw, leniw = get_work_space s
  and nst          = get_num_steps s
  and nfe          = get_num_rhs_evals s
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
  printf "\nFinal Statistics.. \n\n";
  printf "lenrw   = %5d     leniw   = %5d\n"   lenrw leniw;
  printf "lenrwLS = %5d     leniwLS = %5d\n"   lenrwLS leniwLS;
  printf "nst     = %5d\n"                      nst;
  printf "nfe     = %5d     nfeLS   = %5d\n"   nfe nfeLS;
  printf "nni     = %5d     nli     = %5d\n"   nni nli;
  printf "nsetups = %5d     netf    = %5d\n"   nsetups netf;
  printf "npe     = %5d     nps     = %5d\n"   npe nps;
  printf "ncfn    = %5d     ncfl    = %5d\n\n" ncfn ncfl

(* Functions Called by the Solver *)

(* f routine. Compute RHS function f(t,u). *)

let f data t (udata : RealArray.t) (dudata : RealArray.t) =
  (* Set diurnal rate coefficients. *)
  let s = sin (data.om *. t) in
  let q3 = if s > zero then exp(-. a3 /.s) else zero in
  data.q4 <- (if s > zero then exp(-. a4 /. s) else zero);

  (* Make local copies of problem variables, for efficiency. *)
  let q4coef  = data.q4
  and dely    = data.dy
  and verdco  = data.vdco
  and hordco  = data.hdco
  and horaco  = data.haco
  in

  (* Loop over all grid points. *)
  for jy = 0 to my - 1 do

    (* Set vertical diffusion coefficients at jy +- 1/2 *)
    let ydn = ymin +. (float jy -. 0.5) *. dely in
    let yup = ydn +. dely in
    let cydn = verdco *. exp(0.2 *. ydn) in
    let cyup = verdco *. exp(0.2 *. yup) in
    let idn = if jy == 0      then  1 else -1 in
    let iup = if jy == my - 1 then -1 else  1 in

    for jx = 0 to mx - 1 do

      (* Extract c1 and c2, and set kinetic rate terms. *)
      let c1 = udata.{0 + jx * num_species + jy * nsmx}
      and c2 = udata.{1 + jx * num_species + jy * nsmx}
      in
      let qq1 = q1 *. c1 *. c3;
      and qq2 = q2 *. c1 *. c2;
      and qq3 = q3 *. c3;
      and qq4 = q4coef *. c2;
      in
      let rkin1 = -. qq1 -. qq2 +. two *. qq3 +. qq4
      and rkin2 = qq1 -. qq2 -. qq4
      in

      (* Set vertical diffusion terms. *)
      let c1dn = udata.{0 + jx * num_species + (jy + idn) * nsmx}
      and c2dn = udata.{1 + jx * num_species + (jy + idn) * nsmx}
      and c1up = udata.{0 + jx * num_species + (jy + iup) * nsmx}
      and c2up = udata.{1 + jx * num_species + (jy + iup) * nsmx}
      in
      let vertd1 = cyup *. (c1up -. c1) -. cydn *. (c1 -. c1dn)
      and vertd2 = cyup *. (c2up -. c2) -. cydn *. (c2 -. c2dn)
      in

      (* Set horizontal diffusion and advection terms. *)
      let ileft  = if jx = 0      then  1 else -1
      and iright = if jx = mx - 1 then -1 else 1
      in
      let c1lt = udata.{0 + (jx + ileft) * num_species + jy * nsmx}
      and c2lt = udata.{1 + (jx + ileft) * num_species + jy * nsmx}
      and c1rt = udata.{0 + (jx + iright) * num_species + jy * nsmx}
      and c2rt = udata.{1 + (jx + iright) * num_species + jy * nsmx}
      in
      let hord1 = hordco *. (c1rt -. two *. c1 +. c1lt)
      and hord2 = hordco *. (c2rt -. two *. c2 +. c2lt)
      and horad1 = horaco *. (c1rt -. c1lt)
      and horad2 = horaco *. (c2rt -. c2lt)
      in

      (* Load all terms into udot. *)
      dudata.{0 + jx * num_species + jy * nsmx}
                                        <- vertd1 +. hord1 +. horad1 +. rkin1;
      dudata.{1 + jx * num_species + jy * nsmx}
                                        <- vertd2 +. hord2 +. horad2 +. rkin2
    done
  done

(* Jacobian-times-vector routine. *)

let jtv data jac_arg (vdata : RealArray.t) (jvdata : RealArray.t) =
  let open Cvode in
  let { jac_t = t; jac_y = (udata : RealArray.t); } = jac_arg in

  (* Set diurnal rate coefficients. *)
  let s = sin(data.om *. t) in
  data.q4 <- (if s > zero then exp(-. a4 /. s) else zero);

  (* Make local copies of problem variables, for efficiency. *)
  let q4coef = data.q4
  and dely   = data.dy
  and verdco = data.vdco
  and hordco = data.hdco
  and horaco = data.haco
  in
  (* Loop over all grid points. *)
  for jy = 0 to my - 1 do

    (* Set vertical diffusion coefficients at jy +- 1/2 *)
    let ydn = ymin +. (float jy -. 0.5) *. dely in
    let yup = ydn +. dely in

    let cydn = verdco *. exp(0.2 *. ydn)
    and cyup = verdco *. exp(0.2 *. yup)
    in
    let idn = (if jy = 0 then 1 else -1)
    and iup = (if jy = my - 1 then -1 else 1)
    in

    for jx = 0 to mx - 1 do

      let jv1 = zero
      and jv2 = zero
      in

      (* Extract c1 and c2 at the current location and at neighbors *)

      let c1 = udata.{0 + jx * num_species + jy * nsmx}
      and c2 = udata.{1 + jx * num_species + jy * nsmx}
      in

      let v1 = vdata.{0 + jx * num_species + jy * nsmx}
      and v2 = vdata.{1 + jx * num_species + jy * nsmx}
      in

      let v1dn = vdata.{0 + jx * num_species + (jy + idn) * nsmx}
      and v2dn = vdata.{1 + jx * num_species + (jy + idn) * nsmx}
      and v1up = vdata.{0 + jx * num_species + (jy + iup) * nsmx}
      and v2up = vdata.{1 + jx * num_species + (jy + iup) * nsmx}
      in

      let ileft  = (if jx = 0      then  1 else -1)
      and iright = (if jx = mx - 1 then -1 else  1)
      in

      let v1lt = vdata.{0 + (jx + ileft) * num_species + jy * nsmx}
      and v2lt = vdata.{1 + (jx + ileft) * num_species + jy * nsmx}
      and v1rt = vdata.{0 + (jx + iright) * num_species + jy * nsmx}
      and v2rt = vdata.{1 + (jx + iright) * num_species + jy * nsmx}
      in

      (* Set kinetic rate terms. *)
      let jv1 = jv1 +. (-.(q1 *. c3 +. q2 *. c2) *. v1
                        +. (q4coef -. q2 *. c1) *. v2)
      and jv2 = jv2 +. ((q1 *. c3 -. q2 *. c2) *. v1
                        -. (q4coef +. q2 *. c1) *. v2)
      in

      (* Set vertical diffusion terms. *)
      let jv1 = jv1 +. (-. (cyup +. cydn) *. v1 +. cyup *. v1up +. cydn *. v1dn)
      and jv2 = jv2 +. (-. (cyup +. cydn) *. v2 +. cyup *. v2up +. cydn *. v2dn)
      in

      (* Set horizontal diffusion and advection terms. *)
      let jv1 = jv1 +. (hordco *. (v1rt -. two *. v1 +. v1lt))
      and jv2 = jv2 +. (hordco *. (v2rt -. two *. v2 +. v2lt))
      in
      let jv1 = jv1 +. (horaco *. (v1rt -. v1lt))
      and jv2 = jv2 +. (horaco *. (v2rt -. v2lt))
      in
      (* Load two components of J*v *)
      jvdata.{0 + jx * num_species + jy * nsmx} <- jv1;
      jvdata.{1 + jx * num_species + jy * nsmx} <- jv2
    done

  done

(* Preconditioner setup routine. Generate and preprocess P. *)

let precond data jacarg jok gamma =
  let open Cvode in
  let { jac_y   = (udata : RealArray.t); jac_tmp = (); _ } = jacarg in

  (* Make local copies of pointers in user_data, and of pointer to u's data *)
  let p     = data.p
  and jbd   = data.jbd
  and pivot = data.pivot
  in

  let r =
    if jok then begin
      (* jok = TRUE: Copy Jbd to P *)
      for jy = 0 to my - 1 do
        for jx = 0 to mx - 1 do
          Direct.blit ~src:jbd.(jx).(jy) ~dst:p.(jx).(jy)
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
      for jy = 0 to my - 1 do
        let ydn = ymin +. (float jy -. 0.5) *. dely in
        let yup = ydn +. dely in
        let cydn = verdco *. exp(0.2 *. ydn)
        and cyup = verdco *. exp(0.2 *. yup)
        in
        let diag = -. (cydn +. cyup +. two *. hordco) in

        for jx = 0 to mx - 1 do
          let c1 = udata.{0 + jx * num_species + jy * nsmx}
          and c2 = udata.{1 + jx * num_species + jy * nsmx}
          and j = jbd.(jx).(jy)
          and a = p.(jx).(jy)
          in
          set_ijth j 1 1 ((-. q1 *. c3 -. q2 *. c2) +. diag);
          set_ijth j 1 2 (-. q2 *. c1 +. q4coef);
          set_ijth j 2 1 (q1 *. c3 -. q2 *. c2);
          set_ijth j 2 2 ((-. q2 *. c1 -. q4coef) +. diag);
          Direct.blit ~src:j ~dst:a
        done
      done;
      true
    end
  in

  (* Scale by -gamma *)
  for jy = 0 to my - 1 do
    for jx = 0 to mx - 1 do
      Direct.scale (-. gamma) p.(jx).(jy)
    done
  done;

  (* Add identity matrix and do LU decompositions on blocks in place. *)
  for jx = 0 to mx - 1 do
    for jy = 0 to my - 1 do
      Direct.add_identity p.(jx).(jy);
      Direct.getrf p.(jx).(jy) pivot.(jx).(jy)
    done
  done;
  r

(* Preconditioner solve routine *)

let psolve data _ solve_arg (zdata : RealArray.t) =
  let open Cvode.Spils in
  let { rhs = (r : RealArray.t); _ } = solve_arg in

  (* Extract the P and pivot arrays from user_data. *)
  let p = data.p
  and pivot = data.pivot
  in

  Bigarray.Array1.blit r zdata;

  (* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. *)
  for jx = 0 to mx - 1 do
    for jy = 0 to my - 1 do
      let off = jx * num_species + jy * nsmx in
      Direct.getrs' p.(jx).(jy) pivot.(jx).(jy) zdata off
    done
  done

(*
 *-------------------------------
 * Main Program
 *-------------------------------
 *)

let main () =

  (* Allocate memory, and set problem data, initial values, tolerances *)
  let u = Nvector_serial.make neq 0.0 in
  let data = init_user_data (alloc_user_data ()) in
  set_initial_profiles (unwrap u) data.dx data.dy;

  let abstol = atol
  and reltol = rtol
  in

  (* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration *)
  (* Set the pointer to user-defined data *)
  (* Call CVSpgmr to specify the linear solver CVSPGMR
   * with left preconditioning and the maximum Krylov dimension maxl *)
  (* set the Jacobian-times-vector function *)
  (* Set the preconditioner solve and setup functions *)
  let lsolver = Cvode.Spils.spgmr u in
  let cvode_mem = Cvode.(
    init BDF
      ~lsolver:Spils.(solver lsolver ~jac_times_vec:(None, jtv data)
                       (prec_left ~setup:(precond data) (psolve data)))
      (SStolerances (reltol, abstol))
      (f data) t0 u
  ) in

  (* In loop over output points, call CVode, print results, test for error *)

  printf " \n2-species diurnal advection-diffusion problem\n\n";
  let tout = ref twohr in
  for _ = 1 to nout do
    let t, _ = Cvode.solve_normal cvode_mem !tout u in
    print_output cvode_mem (unwrap u) t;
    tout := !tout +. twohr
  done;

  print_final_stats cvode_mem

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
