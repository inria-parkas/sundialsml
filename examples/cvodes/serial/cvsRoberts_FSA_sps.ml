(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2010/12/01 22:58:00 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, and
 *                Radu Serban @ LLNL
 * Programmer(s): Ting Yan @ SMU
 *      Based on cvsRoberts_FSA_dns.c and modified to use SUPERLU_MT
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Dec 2016.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES for Forward Sensitivity
 * Analysis. The problem is from chemical kinetics, and consists
 * of the following three rate equations:
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions y1 = 1.0, y2 = y3 = 0. The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODES SUPERLU_MT linear solver, and a
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
 *    % cvsRoberts_FSA_sps -nosensi
 * If sensitivities are to be computed:
 *    % cvsRoberts_FSA_sps -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 *)

open Sundials

module Sens = Cvodes.Sensitivity
let unwrap = Nvector.unwrap

let printf = Printf.printf

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

let f data _ (y : RealArray.t) (ydot : RealArray.t) =
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2) in

  let yd1 = -.p1*.y.{0} +. p2*.y.{1}*.y.{2} in
  let yd3 = p3*.y.{1}*.y.{1} in
  ydot.{0} <- yd1;
  ydot.{1} <- (-.yd1 -. yd3);
  ydot.{2} <- yd3

(* Jacobian routine. Compute J(t,y). *)

let jac_lt620 data {Cvode.jac_y = (y : RealArray.t)} smat =
  let set_col = Matrix.Sparse.set_col smat in
  let set = Matrix.Sparse.set smat in
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  in
  Matrix.Sparse.set_to_zero smat;

  set_col 0 0;
  set_col 1 3;
  set_col 2 6;
  set_col 3 9;

  set 0 0 (-.p1);
  set 1 1  p1;
  set 2 2  0.00;

  set 3 0 (p2 *. y.{2});
  set 4 1 (-.p2 *. y.{2} -. 2.0 *. p3 *. y.{1});
  set 5 2 (2.0 *. y.{1});

  set 6 0 (p2 *. y.{1});
  set 7 1 (-.p2 *. y.{1});
  set 8 2 0.00

let jac data {Cvode.jac_y = (y : RealArray.t)} smat =
  let set_col = Matrix.Sparse.set_col smat in
  let set = Matrix.Sparse.set smat in
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  in
  Matrix.Sparse.set_to_zero smat;

  (* first column entries start at data[0], two entries (rows 0 and 1) *)
  set_col 0 0;
  set 0 0 (-.p1);
  set 1 1  p1;

  (* second column entries start at data[2], three entries (rows 0, 1, and 2) *)
  set_col 1 2;
  set 2 0 (p2 *. y.{2});
  set 3 1 (-.p2 *. y.{2} -. 2.0 *. p3 *. y.{1});
  set 4 2 (2.0 *. p3 *. y.{1});

  (* third column entries start at data[5], two entries (rows 0 and 1) *)
  set_col 2 5;

  set 5 0 (p2 *. y.{1});
  set 6 1 (-.p2 *. y.{1});

  (* number of non-zeros *)
  set_col 3 7

(* fS routine. Compute sensitivity r.h.s. *)

let fS : user_data -> RealArray.t Sens.sensrhsfn1 =
  fun data iS { Sens.y = y } yS ySdot ->
  let p1 = data.p.(0)
  and p2 = data.p.(1)
  and p3 = data.p.(2)
  and s1 = yS.{0}
  and s2 = yS.{1}
  and s3 = yS.{2}
  in
  let sd1 = -.p1*.s1 +. p2*.y.{2}*.s2 +. p2*.y.{1}*.s3 in
  let sd3 = 2.0*.p3*.y.{1}*.s2 in
  let sd2 = -.sd1-.sd3
  in
  let sd1, sd2, sd3 =
    (match iS with
     | 0 -> (sd1 -. y.{0}, sd2 +. y.{0}, sd3)
     | 1 -> (sd1 +. y.{1}*.y.{2}, sd2 -. y.{1}*.y.{2}, sd3)
     | 2 -> (sd1, sd2 -. y.{1}*.y.{1}, sd3 +. y.{1}*.y.{1})
     | _ -> assert false);
  in
  ySdot.{0} <- sd1;
  ySdot.{1} <- sd2;
  ySdot.{2} <- sd3

(* EwtSet function. Computes the error weights at the current solution. *)

let atol = Array.of_list [ atol1; atol2; atol3 ]
let ewt _ (y : RealArray.t) (w : RealArray.t) =
  let rtol = rtol;
  in
  for i = 0 to 2 do
    let ww = rtol *. (abs_float y.{i}) +. atol.(i) in
    if ww <= 0.0 then raise NonPositiveEwt;
    w.{i} <- (1.0/.ww)
  done

(* Process and verify arguments to cvsfwddenx. *)

let wrong_args name =
  printf "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n" name;
  printf "         sensi_meth = sim, stg, or stg1\n";
  printf "         err_con    = t or f\n";
  exit 0

let process_args () =
  let argv = Sys.argv in
  let argc = Array.length argv in
  let name = argv.(0) in
  if argc < 2 then wrong_args name;

  let sensi =
    if argv.(1) = "-nosensi" then false
    else if argv.(1) = "-sensi" then true
    else wrong_args name
  in

  if not sensi then (None, false)
  else begin
    if argc <> 4 then wrong_args name;
    let sensi_meth =
      if argv.(2) = "sim" then Sens.Simultaneous None
      else if argv.(2) = "stg" then Sens.Staggered None
      else if argv.(2) = "stg1" then Sens.Staggered1 None
      else wrong_args name
    in
    let err_con =
      if argv.(3) = "t" then true
      else if argv.(3) = "f" then false
      else wrong_args name
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
  print_string "                  Solution       ";
  printf "%12.4e %12.4e %12.4e \n" udata.{0} udata.{1} udata.{2}

(* Print sensitivities. *)

let print_output_s uS =
  let sdata = unwrap uS.(0) in
  print_string "                  Sensitivity 1  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2};
  let sdata = unwrap uS.(1) in
  print_string "                  Sensitivity 2  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2};
  let sdata = unwrap uS.(2) in
  print_string "                  Sensitivity 3  ";
  printf "%12.4e %12.4e %12.4e \n"  sdata.{0} sdata.{1} sdata.{2}

(* Print some final statistics from the CVODES memory. *)
(* For high NUM_REPS, the cost of OCaml printf becomes important! *)

let print_string_5d s i =
  print_string s;
  if i < 10 then print_string "    "
  else if i < 100 then print_string "   "
  else if i < 1000 then print_string "  "
  else if i < 10000 then print_string " ";
  print_int i

let print_final_stats s sensi =
  let open Cvode in
  let nst     = get_num_steps s
  and nfe     = get_num_rhs_evals s
  and nsetups = get_num_lin_solv_setups s
  and netf    = get_num_err_test_fails s
  and nni     = get_num_nonlin_solv_iters s
  and nnf     = get_num_nonlin_solv_conv_fails s
  and ncfn    = get_num_step_solve_fails s
  and nje     = Dls.get_num_jac_evals s
  in
  if Sundials_impl.Version.lt620 then begin
    print_string "\nFinal Statistics\n\n";
    print_string_5d "nst     = " nst;
    print_string_5d "\n\nnfe     = " nfe;
    print_string_5d "\nnetf    = " netf;
    print_string_5d "    nsetups  = " nsetups;
    print_string_5d "\nnni     = " nni;
    print_string_5d "    ncfn     = " nnf;
    print_newline ()
  end else begin
    printf "\nFinal Statistics:\n";
    printf "nst = %-6d nfe = %-6d nsetups = %-6d nje = %d\n"
           nst nfe nsetups nje;
    printf "nni = %-6d nnf = %-6d netf = %-6d    ncfn = %-6d\n\n"
           nni nnf netf ncfn
  end;

  if sensi then begin
    let nfSe     = Sens.get_num_rhs_evals s
    and nfeS     = Sens.get_num_rhs_evals_sens s
    and nsetupsS = Sens.get_num_lin_solv_setups s
    and netfS    = Sens.get_num_err_test_fails s
    and nniS     = Sens.get_num_nonlin_solv_iters s
    and nnfS     = Sens.get_num_nonlin_solv_conv_fails s
    and ncfnS    = Sens.get_num_step_solve_fails s
    in
    if Sundials_impl.Version.lt620 then begin
      print_string_5d "\nnfSe    = " nfSe;
      print_string_5d "    nfeS     = " nfeS;
      print_string_5d "\nnetfs   = " netfS;
      print_string_5d "    nsetupsS = " nsetupsS;
      print_string_5d "\nnniS    = " nniS;
      print_string_5d "    ncfnS    = " nnfS;
      print_newline ()
    end else begin
      printf "nfSe = %-6d nfeS = %-6d nsetupsS = %-6d\n"
             nfSe nfeS nsetupsS;
      printf "nniS = %-6d nnfS = %-6d netfS = %-6d ncfnS = %-6d\n\n"
             nniS nnfS netfS ncfnS
    end
  end;

  if Sundials_impl.Version.lt620 then begin
    print_string_5d "\nnje    = " nje;
    print_newline ()
  end

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
  let nthreads = 1 in
  let nnz = neq * neq in
  let m = Matrix.sparse_csc ~nnz neq in
  let jac = (if Sundials_impl.Version.lt620 then jac_lt620 else jac) data in
  let cvode_mem =
    Cvode.(init BDF
                (WFtolerances (ewt data))
                ~lsolver:Dls.(solver ~jac (superlumt ~nthreads y m))
                (f data) t0 y)
  in

  if Sundials_impl.Version.lt620
  then print_string "\n3-species chemical kinetics problem\n"
  else print_string " \n3-species kinetics problem\n";

  (* Sensitivity-related settings *)
  let print_sensi =
    match sensi with
    | None -> (print_string "Sensitivity: NO "; (fun _ -> ()))
    | Some sensi_meth -> begin
        let pbar = RealArray.of_array data.p in

        let yS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0) in

        Sens.(init cvode_mem
                         EEtolerances
                         sensi_meth
                         ~sens_params:{ pvals = None;
                                        pbar = Some pbar;
                                        plist = None; }
                         (OneByOne (Some (fS data)))
                         yS);
        Sens.set_err_con cvode_mem err_con;

        print_string "Sensitivity: YES ";
        (match sensi_meth with
         | Sens.Simultaneous _ -> print_string "( SIMULTANEOUS +"
         | Sens.Staggered _    -> print_string "( STAGGERED +"
         | Sens.Staggered1 _   -> print_string "( STAGGERED1 +");
        print_string (if err_con then " FULL ERROR CONTROL )"
                                 else " PARTIAL ERROR CONTROL )");

        (fun s -> (ignore (Sens.get s yS); print_output_s yS))
      end
  in
  (* In loop over output points, call CVode, print results, test for error *)

  print_string "\n\n\
                ===========================================\
                ============================\n\
               \     T     Q       H      NST           y1\
               \           y2           y3    \n\
                ===========================================\
                ============================\n";

  let tout = ref t1 in
  for _ = 1 to nout do
    let t, _ = Cvode.solve_normal cvode_mem !tout y in
    print_output cvode_mem t ydata;
    print_sensi cvode_mem;
    print_string "-----------------------------------------\
                  ------------------------------\n";
    tout := !tout *. tmult
  done;

  (* Print final statistics *)
  print_final_stats cvode_mem (sensi<>None)

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
