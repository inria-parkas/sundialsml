(*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2009/02/17 02:48:46 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen and Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2014.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(z) = Kv0*exp(z/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= z <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * 10 x 10 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODES, with the BDF/GMRES method
 * (i.e. using the CVSPGMR linear solver) and the block-diagonal
 * part of the Newton matrix as a left preconditioner. A copy of
 * the block-diagonal part of the Jacobian is saved and
 * conditionally reused within the Precond routine.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters q1 and q2.
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on FULL or PARTIAL,
 * respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % cvsDiurnal_FSA_kry -nosensi
 * If sensitivities are to be computed:
 *    % cvsDiurnal_FSA_kry -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module RealArray2 = Sundials.RealArray2
module LintArray = Sundials.LintArray
module Direct = Dls.ArrayDenseMatrix
module Sens = Cvodes.Sensitivity
open Bigarray
let unvec = Nvector.unwrap
let unwrap = RealArray2.unwrap

let printf = Printf.printf

(* Problem Constants *)

let zero  = 0.0
let one   = 1.0
let two   = 2.0

let num_species = 2                (* number of species         *)
let c1_scale    = 1.0e6            (* coefficients in initial profiles    *)
let c2_scale    = 1.0e12

let t0          = zero             (* initial time *)
let nout        = 12               (* number of output times *)
let twohr       = 7200.0           (* number of seconds in two hours  *)
let halfday     = 4.32e4           (* number of seconds in a half day *)
let pi          = 3.1415926535898  (* pi *) 

let xmin        = zero             (* grid boundaries in x  *)
let xmax        = 20.0           
let zmin        = 30.0             (* grid boundaries in z  *)
let zmax        = 50.0
let xmid        = 10.0             (* grid midpoints in x,z *)          
let zmid        = 40.0

let mx          = 15               (* mx = number of x mesh points *)
let mz          = 15               (* my = number of z mesh points *)
let nsmx        = num_species * mx (* nsmx = num_species*mx *)
let mm          = mx * mz          (* mm = mx*mz *)

let rtol        = 1.0e-5           (* scalar relative tolerance *)
let floor       = 100.0            (* value of C1 or C2 at which tolerances *)
                                     (* change from relative to absolute      *)
let atol        = rtol *. floor    (* scalar absolute tolerance *)
let neq         = num_species * mm (* neq = number of equations *)

(* Sensitivity Constants *)
let np = 8
let ns = 2

(* User-defined vector and matrix accessor macros: IJKth, IJth *)

(* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MZ-1. The vdata array is obtained via
   the macro call vdata = NV_DATA_S(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. *)

let ijkth (v : RealArray.t) i j k       = v.{i - 1 + j * num_species + k * nsmx}
let set_ijkth (v : RealArray.t) i j k e = v.{i - 1 + j * num_species + k * nsmx} <- e

(* Type : UserData 
   contains preconditioner blocks, pivot arrays, and problem constants *)

type user_data = {
        params     : RealArray.t;
        p          : Direct.t array array;
        jbd        : Direct.t array array;
        pivot      : LintArray.t array array;
        mutable q4 : float;
        om         : float;
        dx         : float;
        dz         : float;
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
  let new_z_arr elinit _ = Array.init mz elinit in
  let new_xz_arr elinit  = Array.init mx (new_z_arr elinit) in
  {
    params = RealArray.make np 0.0;
    p      = new_xz_arr new_dmat;
    jbd    = new_xz_arr new_dmat;
    pivot  = new_xz_arr new_int1;
    q4     = 0.0;
    om     = 0.0;
    dx     = 0.0;
    dz     = 0.0;
    hdco   = 0.0;
    haco   = 0.0;
    vdco   = 0.0;
  }

(* Load problem constants in data *)

let init_user_data data =
  let q1            = 1.63e-16 in   (* coefficients q1, q2, c3   *) 
  let q2            = 4.66e-16 in
  let c3            = 3.7e16   in
  let a3            = 22.62    in   (* coefficient in expression for q3(t) *)
  let a4            = 7.601    in   (* coefficient in expression for q4(t) *)
  let kh            = 4.0e-6   in   (* horizontal diffusivity Kh *)
  let vel           = 0.001    in   (* advection velocity V      *)
  let kv0           = 1.0e-8   in   (* coefficient in Kv(y)      *)

  data.params.{0} <- q1;
  data.params.{1} <- q2;
  data.params.{2} <- c3;
  data.params.{3} <- a3;
  data.params.{4} <- a4;
  data.params.{5} <- kh;
  data.params.{6} <- vel;
  data.params.{7} <- kv0;

  let dx = (xmax -. xmin) /. float (mx - 1)
  and dz = (zmax -. zmin) /. float (mz - 1)
  in
  { data with
      om = pi /. halfday;
      dx = dx;
      dz = dz;
      hdco = kh /. sqr(dx);
      haco = vel /. (two *. dx);
      vdco = (one /. sqr(dz)) *. kv0;
  }

(* Set initial conditions in u *)

let set_initial_profiles y dx dz =
  for jz = 0 to mz - 1 do
    let z = zmin +. float(jz) *. dz in
    let cz = sqr (0.1 *. (z -. zmid)) in
    let cz = one -. cz +. 0.5 *. sqr(cz) in

    for jx = 0 to mx - 1 do
      let x = xmin +. float jx *. dx in
      let cx = sqr(0.1 *. (x -. xmid)) in
      let cx = one -. cx +. 0.5 *. sqr(cx) in

      set_ijkth y 1 jx jz (c1_scale *. cx *. cz);
      set_ijkth y 2 jx jz (c2_scale *. cx *. cz)
    done
  done

(* Print current t, step count, order, stepsize, and sampled c1,c2 values *)

let print_output s t y =
  let nst = Cvode.get_num_steps s
  and qu  = Cvode.get_last_order s
  and hu  = Cvode.get_last_step s
  and ydata = unvec y
  in
  printf "%8.3e %2d  %8.3e %5d\n"  t qu hu nst;
  printf "                                Solution       ";
  printf "%12.4e %12.4e \n" (ijkth ydata 1 0 0) (ijkth ydata 1 (mx-1) (mz-1)); 
  printf "                                               ";
  printf "%12.4e %12.4e \n" (ijkth ydata 2 0 0) (ijkth ydata 2 (mx-1) (mz-1))

(* Print sampled sensitivities *)

let print_output_s uS =
  let sdata = unvec uS.(0) in
  printf "                                ----------------------------------------\n"; 
  printf "                                Sensitivity 1  ";
  printf "%12.4e %12.4e \n" (ijkth sdata 1 0 0) (ijkth sdata 1 (mx-1) (mz-1)); 
  printf "                                               ";
  printf "%12.4e %12.4e \n" (ijkth sdata 2 0 0) (ijkth sdata 2 (mx-1) (mz-1));

  let sdata = unvec uS.(1) in
  printf "                                ----------------------------------------\n"; 
  printf "                                Sensitivity 2  ";
  printf "%12.4e %12.4e \n" (ijkth sdata 1 0 0) (ijkth sdata 1 (mx-1) (mz-1)); 
  printf "                                               ";
  printf "%12.4e %12.4e \n" (ijkth sdata 2 0 0) (ijkth sdata 2 (mx-1) (mz-1))

(* Get and print final statistics *)

let print_final_stats s sensi =
  let nst     = Cvode.get_num_steps s
  and nfe     = Cvode.get_num_rhs_evals s
  and nsetups = Cvode.get_num_lin_solv_setups s
  and netf    = Cvode.get_num_err_test_fails s
  and nni     = Cvode.get_num_nonlin_solv_iters s
  and ncfn    = Cvode.get_num_nonlin_solv_conv_fails s
  in
  let nli   = Cvode.Spils.get_num_lin_iters s
  and ncfl  = Cvode.Spils.get_num_conv_fails s
  and npe   = Cvode.Spils.get_num_prec_evals s
  and nps   = Cvode.Spils.get_num_prec_solves s
  in
  printf "\nFinal Statistics\n\n";
  printf "nst     = %5d\n\n" nst;
  printf "nfe     = %5d\n"   nfe;
  printf "netf    = %5d    nsetups  = %5d\n" netf nsetups;
  printf "nni     = %5d    ncfn     = %5d\n" nni ncfn;

  if sensi then begin
    let nfSe     = Sens.get_num_rhs_evals s
    and nfeS     = Sens.get_num_rhs_evals_sens s
    and nsetupsS = Sens.get_num_lin_solv_setups s
    and netfS    = Sens.get_num_err_test_fails s
    and nniS     = Sens.get_num_nonlin_solv_iters s
    and ncfnS    = Sens.get_num_nonlin_solv_conv_fails s in
    printf "\n";
    printf "nfSe    = %5d    nfeS     = %5d\n" nfSe nfeS;
    printf "netfs   = %5d    nsetupsS = %5d\n" netfS nsetupsS;
    printf "nniS    = %5d    ncfnS    = %5d\n" nniS ncfnS;
  end;
  printf "\n";
  printf "nli     = %5d    ncfl     = %5d\n" nli ncfl;
  printf "npe     = %5d    nps      = %5d\n" npe nps

(* f routine. Compute f(t,y). *)

let f data t (ydata : RealArray.t) (ydot : RealArray.t) =
  (* Load problem coefficients and parameters *)
  let q1  = data.params.{0} 
  and q2  = data.params.{1}
  and c3  = data.params.{2}
  and a3  = data.params.{3}
  and a4  = data.params.{4}
  in
  (* Set diurnal rate coefficients. *)
  let s = sin (data.om *. t) in
  let q3 = if s > zero then exp(-. a3 /.s) else zero in
  data.q4 <- (if s > zero then exp(-. a4 /. s) else zero);

  (* Make local copies of problem variables, for efficiency. *)
  let q4coef  = data.q4
  and delz    = data.dz
  and verdco  = data.vdco
  and hordco  = data.hdco
  and horaco  = data.haco
  in

  (* Loop over all grid points. *)
  for jz = 0 to mz - 1 do

    (* Set vertical diffusion coefficients at jz +- 1/2 *)
    let zdn = zmin +. (float jz -. 0.5) *. delz in
    let zup = zdn +. delz in
    let czdn = verdco *. exp(0.2 *. zdn) in
    let czup = verdco *. exp(0.2 *. zup) in
    let idn = if jz == 0      then  1 else -1 in
    let iup = if jz == mz - 1 then -1 else  1 in

    for jx = 0 to mx - 1 do

      (* Extract c1 and c2, and set kinetic rate terms. *)
      let c1 = ydata.{jx * num_species + jz * nsmx}
      and c2 = ydata.{1 + jx * num_species + jz * nsmx}
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
      let c1dn = ydata.{0 + jx * num_species + (jz + idn) * nsmx}
      and c2dn = ydata.{1 + jx * num_species + (jz + idn) * nsmx}
      and c1up = ydata.{0 + jx * num_species + (jz + iup) * nsmx}
      and c2up = ydata.{1 + jx * num_species + (jz + iup) * nsmx}
      in
      let vertd1 = czup *. (c1up -. c1) -. czdn *. (c1 -. c1dn)
      and vertd2 = czup *. (c2up -. c2) -. czdn *. (c2 -. c2dn)
      in

      (* Set horizontal diffusion and advection terms. *)
      let ileft  = if jx = 0      then  1 else -1
      and iright = if jx = mx - 1 then -1 else 1
      in
      let c1lt = ydata.{0 + (jx + ileft ) * num_species + jz * nsmx}
      and c2lt = ydata.{1 + (jx + ileft ) * num_species + jz * nsmx}
      and c1rt = ydata.{0 + (jx + iright) * num_species + jz * nsmx}
      and c2rt = ydata.{1 + (jx + iright) * num_species + jz * nsmx}
      in
      let hord1 = hordco *. (c1rt -. two *. c1 +. c1lt)
      and hord2 = hordco *. (c2rt -. two *. c2 +. c2lt)
      and horad1 = horaco *. (c1rt -. c1lt)
      and horad2 = horaco *. (c2rt -. c2lt)
      in

      (* Load all terms into ydot. *)
      ydot.{jx * num_species + jz * nsmx} <- vertd1 +. hord1 +. horad1 +. rkin1;
      ydot.{1 + jx * num_species + jz * nsmx} <- vertd2 +. hord2 +. horad2 +. rkin2
    done
  done

(* Preconditioner setup routine. Generate and preprocess P. *)

let precond data jacarg jok gamma =
  let { Cvode.jac_t   = tn;
        Cvode.jac_y   = (ydata : RealArray.t);
        Cvode.jac_fy  = fydata;
        Cvode.jac_tmp = (vtemp1, vtemp2, vtemp)
      } = jacarg
  in
  (* Load problem coefficients and parameters *)
  let q1  = data.params.{0} 
  and q2  = data.params.{1}
  and c3  = data.params.{2} in

  (* Make local copies of pointers in user_data, and of pointer to u's data *)
  let p     = data.p
  and jbd   = data.jbd
  and pivot = data.pivot
  in

  let r =
    if jok then begin
      (* jok = TRUE: Copy Jbd to P *)
      for jz = 0 to mz - 1 do
        for jx = 0 to mx - 1 do
          Direct.blit jbd.(jx).(jz) p.(jx).(jz)
        done
      done;
      false
    end
    else begin
      (* jok = FALSE: Generate Jbd from scratch and copy to P *)
      (* Make local copies of problem variables, for efficiency. *)
      let q4coef = data.q4
      and delz   = data.dz
      and verdco = data.vdco
      and hordco = data.hdco
      in
      
      (* Compute 2x2 diagonal Jacobian blocks (using q4 values 
         computed on the last f call).  Load into P. *)
      for jz = 0 to mz - 1 do
        let zdn = zmin +. (float jz -. 0.5) *. delz in
        let zup = zdn +. delz in
        let czdn = verdco *. exp(0.2 *. zdn)
        and czup = verdco *. exp(0.2 *. zup)
        in
        let diag = -. (czdn +. czup +. two *. hordco) in

        for jx = 0 to mx - 1 do
          let c1 = ydata.{0 + jx * num_species + jz * nsmx}
          and c2 = ydata.{1 + jx * num_species + jz * nsmx}
          and j = unwrap jbd.(jx).(jz)
          and a = unwrap p.(jx).(jz)
          in
          j.{0, 0} <- (-. q1 *. c3 -. q2 *. c2) +. diag;
          j.{1, 0} <- -. q2 *. c1 +. q4coef;
          j.{0, 1} <- q1 *. c3 -. q2 *. c2;
          j.{1, 1} <- (-. q2 *. c1 -. q4coef) +. diag;
          Array2.blit j a
        done
      done;
      true
    end
  in

  (* Scale by -gamma *)
  for jz = 0 to mz - 1 do
    for jx = 0 to mx - 1 do
      Direct.scale (-. gamma) p.(jx).(jz)
    done
  done;
  
  (* Add identity matrix and do LU decompositions on blocks in place. *)
  for jx = 0 to mx - 1 do
    for jz = 0 to mz - 1 do
      Direct.add_identity p.(jx).(jz);
      Direct.getrf p.(jx).(jz) pivot.(jx).(jz)
    done
  done;
  r

(* Preconditioner solve routine *)

let psolve data jac_arg solve_arg (zdata : RealArray.t) =
  let { Cvode.Spils.rhs = (r : RealArray.t);
        Cvode.Spils.gamma = gamma;
        Cvode.Spils.delta = delta;
        Cvode.Spils.left = lr } = solve_arg
  in

  (* Extract the P and pivot arrays from user_data. *)
  let p = data.p
  and pivot = data.pivot
  in

  Bigarray.Array1.blit r zdata;
  
  (* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. *)
  for jx = 0 to mx - 1 do
    for jz = 0 to mz - 1 do
      (* faster to cache and copy in/out than to Bigarray.Array1.sub... *)
      let off = jx * num_species + jz * nsmx in
      Direct.getrs' p.(jx).(jz) pivot.(jx).(jz) zdata off
    done
  done

(* Process and verify arguments. *)

let wrong_args name =
  printf "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n" name;
  printf "         sensi_meth = sim stg or stg1\n";
  printf "         err_con    = t or f\n";
  exit 0

let process_args () =
  let argv = Sys.argv in
  let argc = Array.length argv in
  if argc < 2 then wrong_args argv.(0);
  
  let sensi =
    if argv.(1) = "-nosensi" then false
    else if argv.(1) = "-sensi" then true
    else wrong_args argv.(0)
  in

  if not sensi then (None, false)
  else begin
    if argc <> 4 then wrong_args argv.(0);
    let sensi_meth =
      if argv.(2) = "sim" then Sens.Simultaneous
      else if argv.(2) = "stg" then Sens.Staggered
      else if argv.(2) = "stg1" then Sens.Staggered1
      else wrong_args argv.(0)
    in
    let err_con =
      if argv.(3) = "t" then true
      else if argv.(3) = "f" then false
      else wrong_args argv.(0)
    in
    (Some sensi_meth, err_con)
  end

(*
 *-------------------------------
 * Main Program
 *-------------------------------
 *)

let main () =
  (* Process arguments *)
  let sensi, err_con = process_args () in

  (* Problem parameters and initial states *)
  let y = Nvector_serial.make neq 0.0 in
  let data = init_user_data (alloc_user_data ()) in
  set_initial_profiles (unvec y) data.dx data.dz;

  (* Tolerances *)
  let abstol = atol
  and reltol = rtol
  in

  (* Create CVODES object *)
  let cvode_mem =
    Cvode.init Cvode.BDF
      (Cvode.Newton
          (Cvode.Spils.spgmr
             (Cvode.Spils.prec_left ~setup:(precond data) (psolve data))))
      (Cvode.SStolerances (reltol, abstol))
      (f data) t0 y
  in
  Cvode.set_max_num_steps cvode_mem 2000;
  printf "\n2-species diurnal advection-diffusion problem\n";

  (* Forward sensitivity analysis *)
  let print_sensi =
    match sensi with
    | None -> (printf "Sensitivity: NO "; (fun _ -> ()))
    | Some sensi_meth -> begin
        let plist = Array.init ns (fun i -> i) in
        let pbar = RealArray.create ns in
        RealArray.mapi (fun is _ -> data.params.{plist.(is)}) pbar;

        let uS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0) in

        Sens.init cvode_mem
                         Sens.EEtolerances
                         sensi_meth
                         ~sens_params:{ Sens.pvals = Some data.params;
                                        Sens.pbar = Some pbar;
                                        Sens.plist = Some plist; }
                         (Sens.OneByOne None)
                         uS;
        Sens.set_err_con cvode_mem err_con;
        Sens.set_dq_method cvode_mem Sens.DQCentered 0.0;

        printf "Sensitivity: YES ";
        (match sensi_meth with
         | Sens.Simultaneous -> printf "( SIMULTANEOUS +"
         | Sens.Staggered    -> printf "( STAGGERED +"
         | Sens.Staggered1   -> printf "( STAGGERED1 +");
        if err_con then printf " FULL ERROR CONTROL )"
                   else printf " PARTIAL ERROR CONTROL )";

        (fun s -> (ignore (Sens.get s uS); print_output_s uS))
      end
  in
  (* In loop over output points, call CVode, print results, test for error *)

  printf "\n\n";
  printf "========================================================================\n";
  printf "     T     Q       H      NST                    Bottom left  Top right \n";
  printf "========================================================================\n";

  let tout = ref twohr in
  for iout = 1 to nout do
    let (t, flag) = Cvode.solve_normal cvode_mem !tout y in
    print_output cvode_mem t y;
    print_sensi cvode_mem;
    tout := !tout +. twohr;
    printf "------------------------------------------------------------------------\n"
  done;

  (* Print final statistics *)
  print_final_stats cvode_mem (sensi <> None)

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
