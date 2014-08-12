(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/04/15 16:37:37 $
 * -----------------------------------------------------------------
 * Programmer(s): Cosmin Petra and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Jul 2014.
 * -----------------------------------------------------------------
 * Example problem:
 *
 * This simple example problem for IDA, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:
 *
 *      dy1/dt = -p1*y1 + p2*y2*y3
 *      dy2/dt = p1*y1 - p2*y2*y3 - p3*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7
 *
 * Optionally, IDAS can compute sensitivities with respect to the
 * problem parameters p1, p2, and p3.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of two sensitivity methods (SIMULTANEOUS and STAGGERED can be
 * used and sensitivities may be included in the error test or not 
 *(error control set on TRUE or FALSE, respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % idasRoberts_FSA_dns -nosensi
 * If sensitivities are to be computed:
 *    % idasRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray
module Quad = Idas.Quadrature
module Sens = Idas.Sensitivity
module QuadSens = Idas.Sensitivity.Quadrature
module VarType = Ida.VarType

let printf = Printf.printf

(* As a slight deviation from the sundials/C code, we allow an extra
   argument to repeat the test, used to check that garbage collection
   works properly.  argv is updated to remove that extra argument so
   the rest of the test can exactly mirror the C code.  *)
let argv = ref Sys.argv

(* Problem Constants *)
let neq   = 3             (* number of equations  *)
let t0    = 0.0           (* initial time *)
let t1    = 0.4           (* first output time *)
let tmult = 10.0          (* output time factor *)
let nout  = 12            (* number of output times *)

let np    = 3             (* number of problem parameters *)
let ns    = 3             (* number of sensitivities computed *)

type user_data =
  {
    p : Ida.real_array;         (* problem parameters (size 3) *)
    coef : float;
  }

let wrong_args name =
  printf "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n" name;
  printf "         sensi_meth = sim or stg\n";
  printf "         err_con    = t or f\n";
  exit 0

(* Process and verify arguments to idasfwddenx. *)
let process_args () =
  let argv = !argv in
  let argc = Array.length argv in

  if argc < 2 then wrong_args argv.(0);

  let sensi =
    match argv.(1) with
    | "-nosensi" -> false
    | "-sensi" -> true
    | _ -> wrong_args argv.(0)
  in

  if sensi then
    begin
      if argc <> 4 then wrong_args argv.(0);
      let sensi_meth =
        match argv.(2) with
        | "sim" -> Sens.Simultaneous
        | "stg" -> Sens.Staggered
        | _ -> wrong_args argv.(0)
      and err_con =
        match argv.(3) with
        | "t" -> true
        | "f" -> false
        | _ -> wrong_args argv.(0)
      in (Some sensi_meth, err_con)
    end
  else (None, false)

let res data t yy yp resval =
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}
  and y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}
  and yp1 = yp.{0}
  and yp2 = yp.{1}
  in

  resval.{0} <- yp1 +. p1*.y1 -. p2*.y2*.y3;
  resval.{1} <- yp2 -. p1*.y1 +. p2*.y2*.y3 +. p3*.y2*.y2;
  resval.{2} <- y1 +. y2 +. y3 -. 1.0

let resS : user_data -> Idas.real_array Sens.sensresfn =
  fun data t yy yp resval yyS ypS resvalS tmp1 tmp2 tmp3 ->
  let p1 = data.p.{0}
  and p2 = data.p.{1}
  and p3 = data.p.{2}
  and y1 = yy.{0}
  and y2 = yy.{1}
  and y3 = yy.{2}
  in

  for is = 0 to ns-1 do
    let s1 = yyS.(is).{0};
    and s2 = yyS.(is).{1};
    and s3 = yyS.(is).{2};
    and sd1 = ypS.(is).{0};
    and sd2 = ypS.(is).{1};
    in
    let rs1 = ref (sd1 +. p1*.s1 -. p2*.y3*.s2 -. p2*.y2*.s3);
    and rs2 = ref (sd2 -. p1*.s1 +. p2*.y3*.s2 +. p2*.y2*.s3 +. 2.0*.p3*.y2*.s2);
    and rs3 = ref (s1 +. s2 +. s3)
    in
    (match is with
     | 0 -> rs1 := !rs1 +. y1;
            rs2 := !rs2 -. y1
     | 1 -> rs1 := !rs1 -. y2*.y3;
            rs2 := !rs2 +. y2*.y3;
     | 2 -> rs2 := !rs2 +. y2*.y2
     | _ -> assert false);
    resvalS.(is).{0} <- !rs1;
    resvalS.(is).{1} <- !rs2;
    resvalS.(is).{2} <- !rs3
  done

let rhsQ data t y yp ypQ =
  ypQ.{0} <- y.{2};
  ypQ.{1} <- data.coef *. (y.{0}*.y.{0}+.
                           y.{1}*.y.{1}+.
                           y.{2}*.y.{2})

let print_ic y yp =
  let data = y in
  printf "\n\nConsistent IC:\n";
  printf "\ty = ";
  printf "%12.4e %12.4e %12.4e \n" data.{0} data.{1} data.{2};

  let data = yp in
  printf "\typ= ";
  printf "%12.4e %12.4e %12.4e \n" data.{0} data.{1} data.{2}

let print_sens_ic y yp yS ypS =
  let sdata = Sundials.unvec yS.(0) in

  printf "                  Sensitivity 1  ";

  printf "\n\ts1 = ";
  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2};

  let sdata = Sundials.unvec ypS.(0) in
  printf "\ts1'= ";
  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2};

  printf "                  Sensitivity 2  ";

  let sdata = Sundials.unvec yS.(1) in
  printf "\n\ts2 = ";
  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2};
  let sdata = Sundials.unvec ypS.(1) in
  printf "\ts2'= ";
  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2};


  printf "                  Sensitivity 3  ";
  let sdata = Sundials.unvec yS.(2) in
  printf "\n\ts3 = ";
  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2};

  let sdata = Sundials.unvec ypS.(2) in
  printf "\ts3'= ";
  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2}

let print_output ida_mem t u =
  let udata = u in
  let nst = Ida.get_num_steps ida_mem
  and qu = Ida.get_last_order ida_mem
  and hu = Ida.get_last_step ida_mem
  in
  printf "%8.3e %2d  %8.3e %5d\n" t qu hu nst;

  printf "                  Solution       ";

  printf "%12.4e %12.4e %12.4e \n" udata.{0} udata.{1} udata.{2}

let print_sens_output uS =
  let sdata = Sundials.unvec uS.(0) in
  printf "                  Sensitivity 1  ";

  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2};

  let sdata = Sundials.unvec uS.(1) in
  printf "                  Sensitivity 2  ";

  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2};

  let sdata = Sundials.unvec uS.(2) in
  printf "                  Sensitivity 3  ";

  printf "%12.4e %12.4e %12.4e \n" sdata.{0} sdata.{1} sdata.{2}

let print_final_stats ida_mem sensi =
  let nst = Ida.get_num_steps ida_mem
  and nfe = Ida.get_num_res_evals ida_mem
  and nsetups = Ida.get_num_lin_solv_setups ida_mem
  and netf = Ida.get_num_err_test_fails ida_mem
  and nni = Ida.get_num_nonlin_solv_iters ida_mem
  and ncfn = Ida.get_num_nonlin_solv_conv_fails ida_mem
  in

  let sens_stats =
    if sensi then
      let nfSe = Sens.get_num_res_evals ida_mem
      and nfeS = Sens.get_num_res_evals_sens ida_mem
      and nsetupsS = Sens.get_num_lin_solv_setups ida_mem
      and netfS = Sens.get_num_err_test_fails ida_mem
      and nniS = Sens.get_num_nonlin_solv_iters ida_mem
      and ncfnS = Sens.get_num_nonlin_solv_conv_fails ida_mem
      in lazy (nfSe, nfeS, nsetupsS, netfS, nniS, ncfnS)
    else
      lazy (failwith "bug in C code transcribed to OCaml")
  in

  let nje = Ida.Dls.get_num_jac_evals ida_mem
  and nfeLS = Ida.Dls.get_num_res_evals ida_mem
  in

  printf "\nFinal Statistics\n\n";
  printf "nst     = %5d\n\n" nst;
  printf "nfe     = %5d\n"   nfe;
  printf "netf    = %5d    nsetups  = %5d\n" netf nsetups;
  printf "nni     = %5d    ncfn     = %5d\n" nni ncfn;

  if sensi then
    begin
      let (nfSe, nfeS, nsetupsS, netfS, nniS, ncfnS) = Lazy.force sens_stats in
      printf "\n";
      printf "nfSe    = %5d    nfeS     = %5d\n" nfSe nfeS;
      printf "netfs   = %5d    nsetupsS = %5d\n" netfS nsetupsS;
      printf "nniS    = %5d    ncfnS    = %5d\n" nniS ncfnS
    end;

  printf "\n";
  printf "nje    = %5d    nfeLS     = %5d\n" nje nfeLS


let main () =
  let sensi, err_con = process_args () in
  let data = { p = RealArray.of_array [|0.040; 1.0e4; 3.0e7|];
               coef = 0.5 } in
  let y  = RealArray.of_array [|1.0; 0.0; 0.0|] in

  (* These initial conditions are NOT consistent. See Ida.calc_ic below. *)
  let yp = RealArray.of_array [|0.1; 0.0; 0.0|] in

  (* Wrap y and yp in nvectors.  Operations performed on the wrapped
     representation affect the originals y and yp.  *)
  let wy  = Nvector_serial.wrap y
  and wyp = Nvector_serial.wrap yp in

  let reltol = 1.0e-6
  and abstol = RealArray.of_array [|1.0e-8; 1.0e-14; 1.0e-6|] in
  let tol = Ida.SVtolerances (reltol, Nvector_serial.wrap abstol) in

  let ida_mem = Ida.init (Ida.Dls.dense None) tol (res data) ~t0:t0 wy wyp in

  printf "\n3-species chemical kinetics problem\n";

  (* Sensitivity-related settings *)
  (* with_yS f either calls f with yS and ypS as arguments (if sensi
     is non-None) or is a no-op (if sensi = None).  *)
  let with_yS =
    match sensi with
    | None -> printf "Sensitivity: NO "; (fun f -> ())
    | Some sensi_meth ->
      let pbar = RealArray.clone data.p in
      let yS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0 (* FIXME: clone y then zero out *)) in
      let ypS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0) in
      (*
        * Only non-zero sensitivity I.C. are ypS[0]: 
        * - Ith(ypS[0],1) = -ONE;
        * - Ith(ypS[0],2) =  ONE;
        *
        * They are not set. IDACalcIC also computes consistent IC for sensitivities.
        *)
      let params = { Sens.pvals = Some data.p;
                     Sens.pbar = Some pbar;
                     Sens.plist = None }
      in
      Sens.init ida_mem Sens.EEtolerances sensi_meth params
        (Some (resS data)) yS ypS;
      Sens.set_err_con ida_mem err_con;

      printf "Sensitivity: YES ";
      if sensi_meth = Sens.Simultaneous then
        printf "( SIMULTANEOUS +"
      else
        printf "( STAGGERED +";
      if err_con then
        printf " FULL ERROR CONTROL )"
      else
        printf " PARTIAL ERROR CONTROL )";
      (fun f -> f yS ypS)
  in
  (*----------------------------------------------------------
   *               Q U A D R A T U R E S
   * ---------------------------------------------------------*)
  let yQ = RealArray.of_array [|0.; 0.|] in
  let wyQ = Nvector_serial.wrap yQ in
  Quad.init ida_mem (rhsQ data) wyQ;

  let yQS = Array.init ns (fun _ -> Nvector_serial.make neq 0.0 (* FIXME: clone y then zero out *)) in

  QuadSens.init ida_mem None yQS;

  (* Call IDACalcIC to compute consistent initial conditions. If sensitivity is
     enabled, this function also try to find consistent IC for the sensitivities. *)
  let id = Nvector_serial.wrap (RealArray.of_array [|VarType.differential;
                                                     VarType.differential;
                                                     VarType.algebraic|])
  in
  if sensi = None
  then (Ida.calc_ic_ya_yd' ~y:wy ~y':wyp ida_mem id t1;
        print_ic y yp)
  else with_yS (fun yS ypS ->
      Sens.calc_ic_ya_yd' ~y:wy ~y':wyp ~ys:yS ~y's:ypS ida_mem id t1;
      print_ic y yp;
      print_sens_ic y yp yS ypS);

  (* In loop over output points, call IDA, print results, test for error *)

  printf "\n\n";
  printf "===========================================";
  printf "============================\n";
  printf "     T     Q       H      NST           y1";
  printf "           y2           y3    \n";
  printf "===========================================";
  printf "============================\n";

  let tout = ref t1 in
  for iout = 1 to nout do
    let t, res = Ida.solve_normal ida_mem !tout wy wyp in
    print_output ida_mem t y;
    with_yS (fun yS ypS ->
        let _ = Sens.get ida_mem yS in
        print_sens_output yS);
    printf "-----------------------------------------";
    printf "------------------------------\n";
    tout := !tout *. tmult
  done;

  printf "\nQuadrature:\n";
  let _ = Quad.get ida_mem wyQ in
  printf "G:      %10.4e\n" yQ.{0};

  with_yS (fun yS ypS ->
      let t = QuadSens.get ida_mem yQS in
      printf "\nSensitivities at t=%g:\n" t;
      printf "dG/dp1: %11.4e\n" (Sundials.unvec yQS.(0)).{0};
      printf "dG/dp1: %11.4e\n" (Sundials.unvec yQS.(1)).{0};
      printf "dG/dp1: %11.4e\n" (Sundials.unvec yQS.(2)).{0});

  (* Print final statistics *)
  print_final_stats ida_mem (sensi <> None)


(* Check if the last argument is a repetition count.  *)
let reps =
  let n = Array.length !argv in
  try
    if n >= 2 then
      let reps = int_of_string !argv.(n-1) in
      argv := Array.sub !argv 0 (n-1);
      reps
    else 1
  with _ -> 1
let _ = for i = 1 to reps do main () done
let _ = Gc.full_major ()