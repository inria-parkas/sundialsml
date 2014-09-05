(*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2010/12/01 23:02:23 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Aug 2014.
 * -----------------------------------------------------------------
 * This simple example problem for IDA, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:

 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.
 *
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01.
 *
 * The problem is solved with IDA using IDADENSE for the linear
 * solver, with a user-supplied Jacobian. Output is printed at
 * t = .4, 4, 40, ..., 4e10.
 * -----------------------------------------------------------------
 *)
module RealArray = Sundials.RealArray
module Roots = Sundials.Roots
module Alt = Ida.Alternate
module DM = Dls.ArrayDenseMatrix
module LintArray = Sundials.LintArray

let printf = Printf.printf

(* Test the Alt module.  This is a re-implementation of the Dense
   direct linear solver in OCaml.  Original C code is found in
     src/ida/ida_direct_impl.h
     src/ida/ida_dense.c
   of the Sundials source tree.  As a simplification, this solver
   requires a user-supplied Jacobian function, with a
   simplified type.
 *)
module AltDense = struct
  type dense_jacfn =
    float                               (* tt *)
    -> float                            (* cj *)
    -> RealArray.t                      (* y *)
    -> RealArray.t                      (* y' *)
    -> RealArray.t                      (* res *)
    -> DM.t                             (* J *)
    -> RealArray.t Ida.triple_tmp
    -> unit

  (* See IDADlsMem in ida_direct_impl.h *)
  type idadls_mem =
    {
      jj : DM.t;
      pivots : LintArray.t;
    }

  let make (jacfn : dense_jacfn) =
    let nje = ref 0 in

    let linit mem s = (nje := 0) in

    let lsetup mem s yp y'p rrp tmps =
      nje := !nje + 1;
      (* Zero out jj; call Jacobian routine jac; return if it failed. *)
      DM.set_to_zero mem.jj;
      let tn = Ida.get_current_time s in
      let cj = Alt.get_cj s in
      jacfn tn cj yp y'p rrp mem.jj tmps;
      (* Do LU factorization of jj; return success or fail flag. *)
      try DM.getrf mem.jj mem.pivots
      with _ -> raise (Sundials.RecoverableFailure false)
    in

    let lsolve mem s b weight ycur y'cur rescur =
      DM.getrs mem.jj mem.pivots b;

      (* Scale the correction to account for change in cj. *)
      let cjratio = Alt.get_cjratio s in
      Nvector_serial.DataOps.n_vscale (2.0/.(1.0 +. cjratio)) b b
    in

    let solver =
      Alt.make_solver (fun s nv nv' ->
          let n = RealArray.length (Sundials.unvec nv) in
          let mem = { jj = DM.create n n;
                      pivots = LintArray.create n;
                    }
          in
          {
            Alt.linit = Some (linit mem);
            Alt.lsetup = Some (lsetup mem);
            Alt.lsolve = lsolve mem;
          })
    in
    (solver, fun () -> (0, !nje))
end

(* Auxiliary indexing functions *)
(* Translates 1-based indexing into 0-based indexing, just like corresponding
 * macros do in the original C implementation of this example.  *)
let ijth a (i,j) = DM.get a (i-1) (j-1)
and set_ijth a (i,j) x = DM.set a (i-1) (j-1) x
and ith (v : RealArray.t) i = v.{i-1}
and set_ith (v : RealArray.t) i x = v.{i-1} <- x

(* Problem Constants *)

let neq    = 3        (* number of equations  *)
let y1     = 1.0      (* initial y components *)
let y2     = 0.0
let y3     = 0.0
let rtol   = 1.0e-4   (* scalar relative tolerance            *)
let atol1  = 1.0e-8   (* vector absolute tolerance components *)
let atol2  = 1.0e-14
let atol3  = 1.0e-6
let tmult  = 10.0     (* output time factor     *)
let nout   = 12       (* number of output times *)
let nroots = 2        (* number of root functions *)

let print_header rtol avtol yy =
  let open Printf in
  printf "\nidaRoberts_dns: Robertson kinetics DAE serial example problem for IDA\n";
  printf "         Three equation chemical kinetics problem.\n\n";
  printf "Linear solver: IDADENSE, with user-supplied Jacobian.\n";

  printf "Tolerance parameters:  rtol = %g   atol = %g %g %g \n"
    rtol avtol.{0} avtol.{1} avtol.{2};
  printf "Initial conditions y0 = (%g %g %g)\n"
    yy.{0} yy.{1} yy.{2};

  printf "Constraints and id not used.\n\n";
  printf "-----------------------------------------------------------------------\n";
  printf "  t             y1           y2           y3";
  printf "      | nst  k      h\n";
  printf "-----------------------------------------------------------------------\n";

and print_output ida t y =
  let kused = Ida.get_last_order ida
  and nst = Ida.get_num_steps ida
  and hused = Ida.get_last_step ida
  in
  printf "%10.4e %12.4e %12.4e %12.4e | %3d  %1d %12.4e\n" 
    t y.{0} y.{1} y.{2} nst kused hused

and print_final_stats ida nreLS nje =
  let nst = Ida.get_num_steps ida
  and nre = Ida.get_num_res_evals ida
  and nni = Ida.get_num_nonlin_solv_iters ida
  and netf = Ida.get_num_err_test_fails ida
  and ncfn = Ida.get_num_nonlin_solv_conv_fails ida
  and nge = Ida.get_num_g_evals ida
  in
  printf "\nFinal Run Statistics: \n\n";
  printf "Number of steps                    = %d\n" nst;
  printf "Number of residual evaluations     = %d\n" (nre+nreLS);
  printf "Number of Jacobian evaluations     = %d\n" nje;
  printf "Number of nonlinear iterations     = %d\n" nni;
  printf "Number of error test failures      = %d\n" netf;
  printf "Number of nonlinear conv. failures = %d\n" ncfn;
  printf "Number of root fn. evaluations     = %d\n" nge;

and print_root_info root_f1 root_f2 =
  (* For printing root_events.  Normally, string_of_root_event makes the output
   * easier to interpret, but we print them as int here in order to get the same
   * output as the C version of this example that comes with sundials.  *)
  let int_of_root_event = function
    | Roots.NoRoot -> 0
    | Roots.Rising -> 1
    | Roots.Falling -> -1
  in
  printf "    rootsfound[] = %3d %3d\n"
    (int_of_root_event root_f1)
    (int_of_root_event root_f2)

let resrob tres (y : RealArray.t) (yp : RealArray.t) (rr : RealArray.t) =
  rr.{0} <- -.0.04*.y.{0} +. 1.0e4*.y.{1}*.y.{2};
  rr.{1} <- -.rr.{0} -. 3.0e7*.y.{1}*.y.{1} -. yp.{1};
  rr.{0} <-  rr.{0} -. yp.{0};
  rr.{2} <-  y.{0} +. y.{1} +. y.{2} -. 1.0

and jacrob : AltDense.dense_jacfn = fun tt cj y y' res jj tmps ->
  set_ijth jj (1,1) (-. 0.04 -. cj);
  set_ijth jj (2,1) (0.04);
  set_ijth jj (3,1) (1.);
  set_ijth jj (1,2) (1.0e4*.y.{2});
  set_ijth jj (2,2) (-. 1.0e4*.y.{2} -. 6.0e7*.y.{1} -. cj);
  set_ijth jj (3,2) (1.);
  set_ijth jj (1,3) (1.0e4*.y.{1});
  set_ijth jj (2,3) (-.1.0e4*.y.{1});
  set_ijth jj (3,3) (1.)

and grob t (y : RealArray.t) y' (gout : RealArray.t) =
  let y1 = y.{0}
  and y3 = y.{2}
  in
  gout.{0} <- y1 -. 0.0001;
  gout.{1} <- y3 -. 0.01

let main () =
  (* Create and initialize y, y', and absolute tolerance vectors.  For
     larger vectors, you might want to use RealArray.create instead of
     RealArray.of_array to avoid making large temporary OCaml
     arrays.  *)
  let y = RealArray.of_array [|1.; 0.; 0.|]
  and y' = RealArray.of_array [|-0.04; 0.04; 0.|]
  and rtol = 1.0e-4
  and avtol = RealArray.of_array [|1.0e-8; 1.0e-14; 1.0e-6|] in
  (* Integration limits *)
  let t0 = 0.0
  and tout1 = 0.4
  in

  (* Wrap y and y' in nvectors.  Operations performed on the wrapped
     representation affect the originals y and y'.  *)
  let wy = Nvector_serial.wrap y
  and wy' = Nvector_serial.wrap y'
  in

  (* Print header information.  *)
  print_header rtol avtol y;

  (* Call IDACreate, IDAInit, and IDARootInit to initialize IDA memory with
   * a 2-component root function and the dense direct linear solver.  *)
  let altdense, get_stats = AltDense.make jacrob in
  let ida_mem =
    Ida.init altdense
             (Ida.SVtolerances (rtol, Nvector_serial.wrap avtol))
             resrob ~roots:(nroots, grob) ~t0:t0 wy wy'
  in
  (* In loop, call IDASolve, print results, and test for error.  Break out of
   * loop when NOUT preset output times have been reached. *)
  let iout = ref 0
  and tout = ref tout1
  in

  let roots = Roots.create nroots in
  let r = Roots.get roots in

  while (!iout <> nout) do
    let (t, flag) = Ida.solve_normal ida_mem !tout wy wy' in
    print_output ida_mem t y;
    match flag with
    | Sundials.RootsFound ->
        Ida.get_root_info ida_mem roots;
        print_root_info (r 0) (r 1)

    | Sundials.Continue ->
        iout := !iout + 1;
        tout := !tout *. tmult

    | Sundials.StopTimeReached ->
        iout := nout
  done;

  let nre, nje = get_stats () in
  print_final_stats ida_mem nre nje

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
