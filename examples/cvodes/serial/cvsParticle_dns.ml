(* ----------------------------------------------------------------------------- {{{
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Oct 2021.
 * OCaml port: Timothy Bourke, Inria, Dec 2022.
 * -----------------------------------------------------------------------------
 * Based on an example from Jean-Luc Fattebert @ ORNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves the equation for a particle moving conterclockwise with
 * velocity alpha on the unit circle in the xy-plane. The ODE system is given by
 *
 *   x' = -alpha * y
 *   y' =  alpha * x
 *
 * where x and y are subject to the constraint
 *
 *   x^2 + y^2 - 1 = 0
 *
 * with initial condition x = 1 and y = 0 at t = 0. The system has the analytic
 * solution
 *
 *  x(t) = cos(alpha * t)
 *  y(t) = sin(alpha * t)
 *
 * For a description of the command line options for this example run the
 * program with the --help flag.
 * --------------------------------------------------------------------------- }}} *)

open Sundials
let printf = Printf.printf
let fprintf = Printf.fprintf

(* Problem Constants *)
let pi   = 3.141592653589793238462643383279502884197169
let zero = 0.0
let one  = 1.0
let two  = 2.0

(* User-defined data structure *)
type user_data =
{
  alpha   : float; (* particle velocity *)

  orbits  : int;   (* number of orbits *)
  torbit  : float; (* orbit time       *)

  rtol    : float; (* integration tolerances *)
  atol    : float;

  proj    : bool;  (* enable/disable solution projection *)
  projerr : bool;  (* enable/disable error projection *)

  tstop   : bool;  (* use tstop mode *)
  nout    : int;   (* number of outputs per orbit *)
}

(* -----------------------------------------------------------------------------
 * Functions provided to CVODES
 * ---------------------------------------------------------------------------*)

(* Compute the right-hand side function, y' = f(t,y) *)
let f { alpha; _ } _ (ydata : RealArray.t) (fdata : RealArray.t) =
  fdata.{0} <- -.alpha *. ydata.{1};
  fdata.{1} <-   alpha *. ydata.{0}

(* Compute the Jacobian of the right-hand side function, J(t,y) = df/dy *)
let jac { alpha; _ } _ jacmat =
  let jdata = Matrix.Dense.unwrap jacmat in
  jdata.{0, 0} <-    0.0;
  jdata.{0, 1} <- -. alpha;
  jdata.{1, 0} <-    alpha;
  jdata.{1, 1} <-    0.0

(* Project the solution onto the constraint manifold *)
let proj _ (ydata : RealArray.t) (cdata : RealArray.t)
           _ (edata : RealArray.t option) =
  let  x = ydata.{0} in
  let  y = ydata.{1} in

  (* project onto the unit circle *)
  let r = sqrt (x *. x +. y *. y) in

  let xp = x /. r in
  let yp = y /. r in

  (* correction to the unprojected solution *)
  cdata.{0} <- xp -. x;
  cdata.{1} <- yp -. y;

  (* project the error *)
  match edata with
  | None -> ()
  | Some edata -> begin
      let errxp =    edata.{0} *. yp *. yp -. edata.{1} *. xp *. yp in
      let erryp = -. edata.{0} *. xp *. yp +. edata.{1} *. xp *. xp in
      edata.{0} <- errxp;
      edata.{1} <- erryp
    end

(* -----------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*)

(* Print command line options *)
let input_help =
  "\nCommand line options:\n\
     --alpha <vel>      : particle velocity\n\
     --orbits <orbits>  : number of orbits to perform\n\
     --rtol <rtol>      : relative tolerance\n\
     --atol <atol>      : absoltue tolerance\n\
     --proj <1 or 0>    : enable (1) / disable (0) projection\n\
     --projerr <1 or 0> : enable (1) / disable (0) error projection\n\
     --nout <nout>      : outputs per period\n\
     --tstop            : stop at output time (do not interpolate)\n"

let init_user_data () =
  let alpha = 1.0 in
  let data = ref {
    alpha;

    orbits  = 100;
    torbit  = (2.0 *. pi) /. alpha;

    rtol    = 1.0e-4;
    atol    = 1.0e-9;

    proj    = true;
    projerr = false;

    tstop   = false;
    nout    = 0;
  } in
  let args = Arg.[
    "--alpha", Float (fun alpha ->
        data := { !data with alpha; torbit = (2.0 *. pi) /. alpha }) , "";

    "--orbits", Int (fun orbits -> data := { !data with orbits }) , "";

    "--rtol", Float (fun rtol -> data := { !data with rtol }) , "";
    "--atol", Float (fun atol -> data := { !data with atol }) , "";

    "--proj", Int (fun proj ->
                     data := { !data with proj = proj <> 0 }) , "";
    "--projerr", Int (fun projerr ->
                        data := { !data with projerr = projerr <> 0 }) , "";

    "--nout", Int (fun nout -> data := { !data with nout }) , "";

    "--tstop", Unit (fun () -> data := { !data with tstop = true }) , "";
  ] in
  let anon_fn s = fprintf stderr "ERROR: Invalid input %s\n%s" s input_help in
  Arg.parse args anon_fn input_help;
  if !data.proj then !data else { !data with projerr = false }

let print_user_data { alpha; orbits; rtol; atol; proj; projerr; nout; tstop } =
  printf "\nParticle traveling on the unit circle example\n";
  printf "---------------------------------------------\n";
  printf "alpha      = %0.4e\n" alpha;
  printf "num orbits = %d\n" orbits;
  printf "---------------------------------------------\n";
  printf "rtol       = %g\n" rtol;
  printf "atol       = %g\n" atol;
  printf "proj sol   = %d\n" (if proj then 1 else 0);
  printf "proj err   = %d\n" (if projerr then 1 else 0);
  printf "nout       = %d\n" nout;
  printf "tstop      = %d\n" (if tstop then 1 else 0);
  printf "---------------------------------------------\n"

(* Compute the analytical solution *)
let compute_solution { alpha; _ } t y =
  let ydata = Nvector.unwrap y in
  ydata.{0} <- cos (alpha *. t);
  ydata.{1} <- sin (alpha *. t)

(* Compute the error in the solution and constraint *)
let compute_error udata t y e =
  let ydata = Nvector.unwrap y in
  (* solution error *)
  compute_solution udata t e;
  Nvector_serial.Ops.linearsum 1.0 y (-1.0) e e;
  (* constraint error *)
  ydata.{0} *. ydata.{0} +. ydata.{1} *. ydata.{1} -. 1.0

(* Output the solution to the screen or disk *)
let write_output t y e ec fids =
  let ydata = Nvector.unwrap y in
  let edata = Nvector.unwrap e in

  match fids with
  | None ->
    (* output solution and error to screen *)
    printf "%0.4e %14.6e %14.6e %14.6e %14.6e %14.6e\n"
      t ydata.{0} ydata.{1} edata.{0} edata.{1} ec
  | Some (yfid, efid) -> begin
    (* output solution to disk *)
    fprintf yfid "%24.16e %24.16e %24.16e\n" t ydata.{0} ydata.{1};

    (* output error to disk *)
    fprintf efid "%24.16e %24.16e %24.16e %24.16e\n" t edata.{0} edata.{1} ec
  end

(* Print final statistics *)
let print_stats cvode_mem =
  let open Cvode in
  let nst     = get_num_steps cvode_mem
  and nfe     = get_num_rhs_evals cvode_mem
  and nsetups = get_num_lin_solv_setups cvode_mem
  and netf    = get_num_err_test_fails cvode_mem
  and nni     = get_num_nonlin_solv_iters cvode_mem
  and ncfn    = get_num_nonlin_solv_conv_fails cvode_mem
  and nje     = Dls.get_num_jac_evals cvode_mem
  in
  printf "\nIntegration Statistics:\n";

  printf "Number of steps taken = %-6d\n" nst;
  printf "Number of function evaluations = %-6d\n" nfe;

  printf "Number of linear solver setups = %-6d\n" nsetups;
  printf "Number of Jacobian evaluations = %-6d\n" nje;

  printf "Number of nonlinear solver iterations = %-6d\n" nni;
  printf "Number of convergence failures = %-6d\n" ncfn;
  printf "Number of error test failures = %-6d\n" netf

(* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*)

let main () =
  (* Allocate and initialize user data structure *)
  let udata = init_user_data () in

  (* Create serial vector to store the solution *)
  (* Set initial contion *)
  let y = Nvector_serial.wrap (RealArray.of_array [| 1.0; 0.0 |]) in

  (* Create serial vector to store the solution error *)
  let e = Nvector.clone y in

  (* Set initial error *)
  Nvector_serial.Ops.const 0.0 e;

  (* Create dense SUNMatrix for use in linear solves *)
  (* Create CVODES memory *)
  (* Initialize CVODES *)
  (* Create dense SUNLinearSolver object *)
  (* Attach the matrix and linear solver to CVODES *)
  (* Set a user-supplied Jacobian function *)
  (* Attach user-defined data structure to CVODES *)
  (* Set integration tolerances *)
  let a = Matrix.dense 2 in
  let projfn = if udata.proj then Some proj else None in
  let cvode_mem =
    Cvode.(init BDF ~lsolver:Dls.(solver ~jac:(jac udata) (dense y a))
                    (SStolerances (udata.rtol, udata.atol))
                    (f udata) ?projfn 0.0 y)
  in
  (* Set a user-supplied projection function *)
  if udata.proj then Cvode.set_proj_err_est cvode_mem udata.projerr;

  (* Set max steps between outputs *)
  Cvode.set_max_num_steps cvode_mem 100000;

  (* Output problem setup *)
  print_user_data udata;

  (* Output initial condition *)
  printf "\n     t            x              y";
  printf "             err x          err y       err constr\n";
  write_output 0.0 y e 0.0 None;

  let fids =
    if udata.nout > 0 then
      Some (open_out "cvsParticle_solution.txt",
            open_out "cvsParticle_error.txt")
    else None
  in

  (* Integrate in time and periodically output the solution and error *)
  let totalout, dtout =
    if udata.nout > 0
    then (udata.orbits * udata.nout, udata.torbit /. float udata.nout)
    else (1, udata.torbit *. float udata.orbits)
  in
  let tout = ref dtout in
  let tfinal = ref 0.0 in

  for out = 0 to totalout - 1 do
    (* Stop at output time (do not interpolate output) *)
    if udata.tstop || udata.nout = 0 then
      Cvode.set_stop_time cvode_mem !tout;

    (* Advance in time *)
    let t, _ = Cvode.solve_normal cvode_mem !tout y in
    tfinal := t;

    (* Output solution and error *)
    if udata.nout > 0 then
      let ec = compute_error udata t y e in
      write_output t y e ec fids;

    (* Update output time *)
    tout := if out < totalout - 1
            then !tout +. dtout
            else udata.torbit *. float udata.orbits
  done;

  (* Close output files *)
  (match fids with None -> ()
   | Some (yfid, efid) -> close_out yfid; close_out efid);

  (* Output final solution and error to screen *)
  let ec = compute_error udata !tfinal y e in
  write_output !tfinal y e ec None;

  (* Print some final statistics *)
  print_stats cvode_mem

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

