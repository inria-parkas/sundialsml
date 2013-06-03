(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2008/12/29 22:21:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Ocaml port: Jun Inoue, INRIA, May 2013.
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 *)
module Ida = Ida_serial;;
module Carray = Ida.Carray
module Roots = Ida.Roots
module Dls = Ida.Dls

let printf = Printf.printf

(* Auxiliary indexing functions *)
(* Translates 1-based indexing into 0-based indexing, just like corresponding
 * macros do in the original C implementation of this example.  *)
let ijth a (i,j) = Ida.Densematrix.get a (i-1, j-1)
and set_ijth a (i,j) x = Ida.Densematrix.set a (i-1, j-1) x
and ith v i = v.{i-1}
and set_ith v i x = v.{i-1} <- x


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
and print_final_stats ida =
  let nst = Ida.get_num_steps ida
  and nre = Ida.get_num_res_evals ida
  and nje = Ida.Dls.get_num_jac_evals ida
  and nni = Ida.get_num_nonlin_solv_iters ida
  and netf = Ida.get_num_err_test_fails ida
  and ncfn = Ida.get_num_nonlin_solv_conv_fails ida
  and nreLS = Ida.Dls.get_num_res_evals ida
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
    (int_of_root_event root_f2);
;;

let resrob tres y yp rr =
  rr.{0} <- -.0.04*.y.{0} +. 1.0e4*.y.{1}*.y.{2};
  rr.{1} <- -.rr.{0} -. 3.0e7*.y.{1}*.y.{1} -. yp.{1};
  rr.{0} <-  rr.{0} -. yp.{0};
  rr.{2} <-  y.{0} +. y.{1} +. y.{2} -. 1.0

and jacrob params jj =
  match params with
    { Ida.jac_t=tt; Ida.jac_coef=cj; Ida.jac_y=y; Ida.jac_res=resvec } ->
      set_ijth jj (1,1) (-. 0.04 -. cj);
      set_ijth jj (2,1) (0.04);
      set_ijth jj (3,1) (1.);
      set_ijth jj (1,2) (1.0e4*.y.{2});
      set_ijth jj (2,2) (-. 1.0e4*.y.{2} -. 6.0e7*.y.{1} -. cj);
      set_ijth jj (3,2) (1.);
      set_ijth jj (1,3) (1.0e4*.y.{1});
      set_ijth jj (2,3) (-.1.0e4*.y.{1});
      set_ijth jj (3,3) (1.)
and grob t y y' gout =
  let y1 = y.{0}
  and y3 = y.{2}
  in
  gout.{0} <- y1 -. 0.0001;
  gout.{1} <- y3 -. 0.01
;;

let main () =
  (* Create and initialize y, y', and absolute tolerance vectors.  For larger
   * vectors, you might want to use Carray.create instead of Carray.of_array to
   * avoid making large temporary OCaml arrays.  *)
  let y = Carray.of_array [|1.; 0.; 0.|]
  and y' = Carray.of_array [|-0.04; 0.04; 0.|]
  and rtol = 1.0e-4
  and avtol = Carray.of_array [|1.0e-8; 1.0e-14; 1.0e-6|] in
  (* Integration limits *)
  let t0 = 0.0
  and tout1 = 0.4
  in

  (* Print header information.  *)
  print_header rtol avtol y;

  (* Call IDACreate, IDAInit, and IDARootInit to initialize IDA memory with
   * a 2-component root function and the dense direct linear solver.  *)
  let ida_mem = Ida.init_at_time Ida.Dense resrob (nroots, grob) t0 y y' in

  (* Call IDASVtolerances to set tolerances *)
  Ida.sv_tolerances ida_mem rtol avtol;

  (* Call IDADense and set up the linear solver. *)
  Ida.Dls.set_dense_jac_fn ida_mem jacrob;

  (* In loop, call IDASolve, print results, and test for error.  Break out of
   * loop when NOUT preset output times have been reached. *)
  let iout = ref 0
  and tout = ref tout1
  in

  let roots = Roots.create nroots in
  let r = Roots.get' roots in

  while (!iout <> nout) do
    let (t, flag) = Ida.solve_normal ida_mem !tout y y' in
    print_output ida_mem t y;
    match flag with
    | Ida.RootsFound ->
        Ida.get_root_info ida_mem roots;
        print_root_info (r 0) (r 1)

    | Ida.Continue ->
        iout := !iout + 1;
        tout := !tout *. tmult

    | Ida.StopTimeReached ->
        iout := nout
  done;

  print_final_stats ida_mem

let _ = main ()
let _ = Gc.compact ()
