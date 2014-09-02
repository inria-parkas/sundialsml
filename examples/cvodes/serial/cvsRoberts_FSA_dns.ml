(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2010/12/01 22:58:00 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Jun 2014.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES. The problem is from chemical
 * kinetics, and consists of the following three rate equations:
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions y1 = 1.0, y2 = y3 = 0. The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODES dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters p1, p2, and p3.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on TRUE or FALSE,
 * respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % cvsRoberts_FSA_dns -nosensi
 * If sensitivities are to be computed:
 *    % cvsRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 *)

module Sens = Cvodes.Sensitivity
module RealArray = Sundials.RealArray
module Densemat = Dls.DenseMatrix
let unvec = Sundials.unvec

let printf = Printf.printf

(* Accessor macros *)

(* i-th vector component i=1..NEQ *)
let ith (v : RealArray.t) i = v.{i - 1}
let set_ith (v : RealArray.t) i e = v.{i - 1} <- e

(* (i,j)-th matrix component i,j=1..NEQ *)
let ijth v i j       = Densemat.get v (i - 1) (j - 1)
let set_ijth v i j e = Densemat.set v (i - 1) (j - 1) e

(* Problem Constants *)

let neq   = 3     (* number of equations  *)
let y1    = 1.0   (* initial y components *)
let y2    = 0.0
let y3    = 0.0
let rtol  = 1e-4  (* scalar relative tolerance *)
let atol1 = 1e-8  (* vector absolute tolerance components *)
let atol2 = 1e-14
let atol3 = 1e-6
let t0    = 0.0   (* initial time *)
let t1    = 0.4   (* first output time *)
let tmult = 10.0  (* output time factor *)
let nout  = 12    (* number of output times *)

let np    = 3     (* number of problem parameters *)
let ns    = 3     (* number of sensitivities computed *)

let zero  = 0.0

(* Type : UserData *)

type user_data = { p : float array }

(* Prototypes of functions by CVODES *)

(* f routine. Compute f(t,y). *)

let f data t y ydot =
  let y1 = ith y 1
  and y2 = ith y 2
  and y3 = ith y 3
  and p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2) in

  let yd1 = -.p1*.y1 +. p2*.y2*.y3 in
  let yd3 = p3*.y2*.y2 in
  set_ith ydot 1 yd1;
  set_ith ydot 3 yd3; 
  set_ith ydot 2 (-.yd1 -. yd3)

(* Jacobian routine. Compute J(t,y). *)

let jac data { Cvode.jac_y = y } jmat =
  let y2 = ith y 2
  and y3 = ith y 3
  and p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  in
  set_ijth jmat 1 1 (-.p1);
  set_ijth jmat 1 2 (p2*.y3);
  set_ijth jmat 1 3 (p2*.y2);
  set_ijth jmat 2 1 ( p1);
  set_ijth jmat 2 2 (-.p2*.y3-.2.0*.p3*.y2);
  set_ijth jmat 2 3 (-.p2*.y2);
  set_ijth jmat 3 2 (2.0*.p3*.y2)
 
(* fS routine. Compute sensitivity r.h.s. *)

let fS data t y ydot iS yS ySdot tmp1 tmp2 =
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  and y1 = ith y 1
  and y2 = ith y 2
  and y3 = ith y 3
  and s1 = ith yS 1
  and s2 = ith yS 2
  and s3 = ith yS 3
  in
  let sd1 = -.p1*.s1 +. p2*.y3*.s2 +. p2*.y2*.s3 in
  let sd3 = 2.0*.p3*.y2*.s2 in
  let sd2 = -.sd1-.sd3
  in
  let sd1, sd2, sd3 =
    (match iS with
     | 0 -> (sd1 -. y1, sd2 +. y1, sd3)
     | 1 -> (sd1 +. y2*.y3, sd2 -. y2*.y3, sd3)
     | 2 -> (sd1, sd2 -. y2*.y2, sd3 +. y2*.y2)
     | _ -> assert false);
  in
  set_ith ySdot 1 sd1;
  set_ith ySdot 2 sd2;
  set_ith ySdot 3 sd3

(* EwtSet function. Computes the error weights at the current solution. *)

let ewt data y w =
  let rtol = rtol;
  and atol = Array.of_list [ atol1; atol2; atol3 ]
  in
  for i = 1 to 3 do
    let yy = ith y i in
    let ww = rtol *. (abs_float yy) +. atol.(i-1) in
    if ww <= 0.0 then raise Sundials.NonPositiveEwt;
    set_ith w i (1.0/.ww)
  done

(* Process and verify arguments to cvsfwddenx. *)

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

(* Print current t, step count, order, stepsize, and solution. *)

let print_output s t udata =
  let nst = Cvode.get_num_steps s
  and qu  = Cvode.get_last_order s
  and hu  = Cvode.get_last_step s
  in
  printf "%8.3e %2d  %8.3e %5d\n" t qu hu nst;
  printf "                  Solution       ";
  printf "%12.4e %12.4e %12.4e \n" udata.{0} udata.{1} udata.{2}

(* Print sensitivities. *)

let print_output_s uS =
  let sdata = unvec uS.(0) in
  printf "                  Sensitivity 1  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2};
  let sdata = unvec uS.(1) in
  printf "                  Sensitivity 2  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2};
  let sdata = unvec uS.(2) in
  printf "                  Sensitivity 3  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2}

(* Print some final statistics from the CVODES memory. *)

let print_final_stats s sensi =
  let nst     = Cvode.get_num_steps s
  and nfe     = Cvode.get_num_rhs_evals s
  and nsetups = Cvode.get_num_lin_solv_setups s
  and netf    = Cvode.get_num_err_test_fails s
  and nni     = Cvode.get_num_nonlin_solv_iters s
  and ncfn    = Cvode.get_num_nonlin_solv_conv_fails s
  in
  printf "\nFinal Statistics\n\n";
  printf "nst     = %5d\n\n" nst;
  printf "nfe     = %5d\n" nfe;
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
    printf "nniS    = %5d    ncfnS    = %5d\n" nniS ncfnS
  end;

  let nje   = Cvode.Dls.get_num_jac_evals s
  and nfeLS = Cvode.Dls.get_num_rhs_evals s
  in
  printf "\n";
  printf "nje    = %5d    nfeLS     = %5d\n"  nje  nfeLS

(*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 *)

let main () =
  (* Process arguments *)
  let sensi, err_con = process_args () in

  (* User data structure *)
  let data = { p = Array.of_list [ 0.04; 1.0e4; 3.0e7 ] } in

  (* Initial conditions *)
  let ydata = RealArray.of_list [y1; y2; y3] in
  let y = Nvector_serial.wrap ydata in

  (* Create CVODES object *)
  let cvode_mem =
    Cvode.init Cvode.BDF (Cvode.Newton (Cvode.Dls.dense (Some (jac data))))
      (Cvode.WFtolerances (ewt data)) (f data) ~t0:t0 y
  in

  printf "\n3-species chemical kinetics problem\n";

  (* Sensitivity-related settings *)
  let print_sensi =
    match sensi with
    | None -> (printf "Sensitivity: NO "; (fun _ -> ()))
    | Some sensi_meth -> begin
        let pbar = RealArray.of_array data.p in

        let yS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0) in

        Sens.init cvode_mem
                         Sens.EEtolerances
                         sensi_meth
                         { Sens.pvals = None;
                           Sens.pbar = Some pbar;
                           Sens.plist = None; }
                         (Sens.OneByOne (Some (fS data)))
                         yS;
        Sens.set_err_con cvode_mem err_con;

        printf "Sensitivity: YES ";
        (match sensi_meth with
         | Sens.Simultaneous -> printf "( SIMULTANEOUS +"
         | Sens.Staggered    -> printf "( STAGGERED +"
         | Sens.Staggered1   -> printf "( STAGGERED1 +");
        printf (if err_con then " FULL ERROR CONTROL )"
                           else " PARTIAL ERROR CONTROL )");

        (fun s -> (ignore (Sens.get s yS); print_output_s yS))
      end
  in
  (* In loop over output points, call CVode, print results, test for error *)
  
  printf "\n\n";
  printf "===========================================";
  printf "============================\n";
  printf "     T     Q       H      NST           y1";
  printf "           y2           y3    \n";
  printf "===========================================";
  printf "============================\n";

  let tout = ref t1 in
  for iout = 1 to nout do
    let t, _ = Cvode.solve_normal cvode_mem !tout y in
    print_output cvode_mem t ydata;
    print_sensi cvode_mem;
    printf "-----------------------------------------";
    printf "------------------------------\n";
    tout := !tout *. tmult
  done;

  (* Print final statistics *)
  print_final_stats cvode_mem (sensi<>None)

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
