(* ---------------------------------------------------------------- {{{
 * Programmer(s): Daniel R. Reynolds @ SMU
 * ----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Nov 2021.
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
 *
 *    [u]' = [ G  e ] [(-1+u^2-r)/(2u)] + [      r'(t)/(2u)        ]
 *    [v]    [ e -1 ] [(-2+v^2-s)/(2v)]   [ s'(t)/(2*sqrt(2+s(t))) ]
 *         = [ fs(t,u,v) ]
 *           [ ff(t,u,v) ]
 *
 * where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.
 *
 * This problem has analytical solution given by
 *    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
 *
 * We use the parameters:
 *   e = 0.5 (fast/slow coupling strength) [default]
 *   G = -1e2 (stiffness at slow time scale) [default]
 *   w = 100  (time-scale separation factor) [default]
 *   hs = 0.01 (slow step size) [default]
 *
 * The stiffness of the slow time scale is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to a
 * multirate method that is implicit at the slow time scale.
 *
 * We select the MRI method to use based on an additional input,
 * solve_type; with options (slow type-order/fast type-order):
 * 0. exp-3/exp-3 (standard MIS) [default]
 * 1. none/exp-3 (no slow, explicit fast)
 * 2. none/dirk-3 (no slow, dirk fast)
 * 3. exp-3/none (explicit slow, no fast)
 * 4. dirk-2/none (dirk slow, no fast) -- solve-decoupled
 * 5. exp-4/exp-4 (MRI-GARK-ERK45a / ERK-4-4)
 * 6. exp-4/exp-3 (MRI-GARK-ERK45a / ERK-3-3)
 * 7. dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled
 *
 * The program should be run with arguments in the following order:
 *   $ a.out solve_type h G w e
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out solve_type h G w
 *   $ a.out solve_type h G
 *   $ a.out solve_type h
 *   $ a.out solve_type
 *   $ a.out
 * are acceptable.  We require:
 *   * 0 <= solve_type <= 7
 *   * 0 < h < 1/|G|
 *   * G < 0.0
 *   * w >= 1.0
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 * ---------------------------------------------------------------- }}}*)

open Sundials
module ARKStep = Arkode.ARKStep
module MRIStep = Arkode.MRIStep

module DM = Matrix.Dense

let printf = Printf.printf
let fprintf = Printf.fprintf

(* ------------------------------
 * Private helper functions
 * ------------------------------*)

let r t = 0.5 *. cos t
let s (rpar : RealArray.t) t = cos (rpar.{1} *. t)
let rdot t = -0.5 *. sin t
let sdot (rpar : RealArray.t) t = -. rpar.{1} *. sin(rpar.{1} *. t)
let utrue (rpar : RealArray.t) t = sqrt (1.0 +. r t)
let vtrue (rpar : RealArray.t) t = sqrt (2.0 +. s rpar t)

let ytrue (rpar : RealArray.t) t (y : RealArray.t) =
  y.{0} <- utrue rpar t;
  y.{1} <- vtrue rpar t

(* ------------------------------
 * Functions called by the solver
 * ------------------------------*)

(* ff routine to compute the fast portion of the ODE RHS. *)
let ff (rpar : RealArray.t) t (y : RealArray.t) (ydot : RealArray.t) =
  let e = rpar.{2} in
  let u = y.{0} in
  let v = y.{1} in
  (* fill in the RHS function:
     [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
     [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] *)
  let tmp1 = (-1.0 +. u *. u -. r t) /. (2.0 *. u) in
  let tmp2 = (-2.0 +. v *. v -. s rpar t) /. (2.0 *. v) in
  ydot.{0} <- 0.0;
  ydot.{1} <- e *. tmp1 -. tmp2 +. sdot rpar t /. (2.0 *. vtrue rpar t)

(* fs routine to compute the slow portion of the ODE RHS. *)
let fs (rpar : RealArray.t) t (y : RealArray.t) (ydot : RealArray.t) =
  let g = rpar.{0} in
  let e = rpar.{2} in
  let u = y.{0} in
  let v = y.{1} in
  (* fill in the RHS function:
     [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
     [0 0] [(-2+v^2-s(t))/(2*v)]    [      0      ] *)
  let tmp1 = (-1.0 +. u *. u -. r t) /. (2.0 *. u) in
  let tmp2 = (-2.0 +. v *. v -. s rpar t) /. (2.0 *. v) in
  ydot.{0} <- g *. tmp1 +. e *. tmp2 +. rdot t /. (2.0 *. u);
  ydot.{1} <- 0.0

let fn (rpar : RealArray.t) t (y : RealArray.t) (ydot : RealArray.t) =
  let g = rpar.{0} in
  let e = rpar.{2} in
  let u = y.{0} in
  let v = y.{1} in
  (* fill in the RHS function:
     [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
     [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] *)
  let tmp1 = (-1.0 +. u *. u -. r t) /. (2.0 *. u) in
  let tmp2 = (-2.0 +. v *. v -. s rpar t) /. (2.0 *. v) in
  ydot.{0} <- g *. tmp1 +. e *. tmp2 +. rdot t /. (2.0 *. u);
  ydot.{1} <- e *. tmp1 -. tmp2 +. sdot rpar t /. (2.0 *. vtrue rpar t)

let f0 (rpar : RealArray.t) t (y : RealArray.t) (ydot : RealArray.t) =
  RealArray.fill ydot 0.0

let js (rpar : RealArray.t) MRIStep.{ jac_t = t;
                                      jac_y = (y : RealArray.t);
                                      jac_fy = (ydot : RealArray.t); _ }
                            jmat =
  let g = rpar.{0} in
  let e = rpar.{2} in
  let u = y.{0} in
  let v = y.{1} in
  (* fill in the Jacobian:
     [G/2 + (G*(1+r(t))+rdot(t))/(2*u^2)   e/2+e*(2+s(t))/(2*v^2)]
     [                 0                             0           ] *)
  DM.set jmat 0 0 (g /. 2.0 +. (g *. (1.0 +. r t) +. rdot t) /. (2.0 *. u *. u));
  DM.set jmat 0 1 (e /. 2.0 +. e *. (2.0 +. s rpar t) /. (2.0 *. v *. v));
  DM.set jmat 1 0 0.0;
  DM.set jmat 1 1 0.0

let jn (rpar : RealArray.t) MRIStep.{ jac_t = t;
                                      jac_y = (y : RealArray.t);
                                      jac_fy = (ydot : RealArray.t); _ }
                            jmat =
  let g = rpar.{0} in
  let e = rpar.{2} in
  let u = y.{0} in
  let v = y.{1} in
  (* fill in the Jacobian:
     [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)     e/2+e*(2+s(t))/(2*v^2)]
     [e/2+e*(1+r(t))/(2*u^2)                -1/2 - (2+s(t))/(2*v^2)] *)
  DM.set jmat 0 0 (g /. 2.0 +. (g *. (1.0 +. r t) -. rdot t) /. (2.0 *. u*. u));
  DM.set jmat 0 1 (e /. 2.0 +. e *.(2.0 +. s rpar t) /. (2.0*.v*.v));
  DM.set jmat 1 0 (e /. 2.0 +. e *.(1.0 +. r t) /. (2.0*.u*.u));
  DM.set jmat 1 1 (-1.0 /. 2.0 -. (2.0 +. s rpar t) /. (2.0*.v*.v))

(* Main Program *)
let main () =
  (* general problem parameters *)
  let t0         = 0.0 in     (* initial time *)
  let tf         = 5.0 in     (* final time *)
  let dTout      = 0.1 in     (* time between outputs *)
  let neq        = 2 in       (* number of dependent vars. *)
  let nt         = Int.of_float (ceil (tf /. dTout)) in (* number of output times *)
  let default_reltol =    0.01 in
  let default_abstol = 1e-11 in

  (*
   * Initialization
   *)

  (* Retrieve the command-line options: solve_type h G w e *)
  let argc = Array.length Sys.argv in
  (* problem configuration type *)
  let solve_type = if argc > 1 then int_of_string Sys.argv.(1) else 0 in
  (* slow step size *)
  let hs = if argc > 2 then float_of_string Sys.argv.(2) else  0.01 in
  (* stiffness at slow time scale *)
  let g = if argc > 3 then float_of_string Sys.argv.(3) else -100.0 in
  (* time-scale separation factor *)
  let w = if argc > 4 then float_of_string Sys.argv.(4) else  100.0 in
  (* fast/slow coupling strength *)
  let e = if argc > 5 then float_of_string Sys.argv.(5) else    0.5 in

  (* Check arguments for validity *)
  (*   0 <= solve_type <= 7      *)
  (*   G < 0.0                   *)
  (*   h > 0                     *)
  (*   h < 1/|G| (explicit slow) *)
  (*   w >= 1.0                  *)
  if solve_type < 0 || solve_type > 7
  then failwith "ERROR: solve_type be an integer in [0,7] \n";
  if g >= 0.0
  then failwith "ERROR: G must be a negative real number\n";

  let implicit_slow = solve_type = 4 || solve_type = 7 in

  if hs <= 0.0 then failwith "ERROR: hs must be in positive\n";
  if hs > 1.0 /. abs_float g && not implicit_slow
  then failwith "ERROR: hs must be in (0, 1/|G|)\n";

  if w < 1.0 then failwith "ERROR: w must be >= 1.0\n";

  let rpar = RealArray.of_list [ g; w; e ] in
  let hf = hs /. w in

  (* Initial problem output (and set implicit solver tolerances as needed) *)
  printf "\nMultirate nonlinear Kvaerno-Prothero-Robinson test problem:\n";
  printf "    time domain:  (%g,%g]\n" t0 tf;
  printf "    hs = %g\n" hs;
  printf "    hf = %g\n" hf;
  printf "    G = %g\n" g;
  printf "    w = %g\n" w;
  printf "    e = %g\n" e;
  let reltol, abstol =
    match solve_type with
    | 0 ->
      printf "    solver: exp-3/exp-3 (standard MIS)\n\n";
      default_reltol, default_abstol

    | 1 ->
      printf "    solver: none/exp-3 (no slow, explicit fast)\n\n";
      default_reltol, default_abstol

    | 2 ->
      let reltol = max (hs *. hs *. hs) 1e-10 in
      let abstol = 1e-11 in
      printf "    solver: none/dirk-3 (no slow, dirk fast)\n\n";
      printf "    reltol = %.2e,  abstol = %.2e\n" reltol abstol;
      reltol, abstol

    | 3 ->
      printf "    solver: exp-3/none (explicit slow, no fast)\n";
      default_reltol, default_abstol

    | 4 ->
      let reltol = max (hs *. hs) 1e-10 in
      let abstol = 1e-11 in
      printf "    solver: dirk-2/none (dirk slow, no fast)\n";
      printf "    reltol = %.2e,  abstol = %.2e\n" reltol abstol;
      reltol, abstol

    | 5 ->
      printf "    solver: exp-4/exp-4 (MRI-GARK-ERK45a / ERK-4-4)\n\n";
      default_reltol, default_abstol

    | 6 ->
      printf "    solver: exp-4/exp-3 (MRI-GARK-ERK45a / ERK-3-3)\n\n";
      default_reltol, default_abstol

    | 7 ->
      let reltol = max (hs *. hs *. hs) 1e-10 in
      let abstol = 1e-11 in
      printf "    solver: dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled\n";
      printf "    reltol = %.2e,  abstol = %.2e\n" reltol abstol;
      reltol, abstol

    | _ -> assert false;
  in

  (* Create and initialize serial vector for the solution *)
  let y = Nvector_serial.make neq 0.0 in
  let ydata = Nvector.unwrap y in
  ytrue rpar t0 ydata;

  (*
   * Create the fast integrator and set options
   *)

  (* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, and the
     initial dependent variable vector y. *)
  let inner_arkode_mem =
    match solve_type with
    | 0 | 6 | 7 ->  (* erk-3-3 fast solver *)
        let inner_arkode_mem =
          ARKStep.(init (explicit (ff rpar)) (SStolerances (reltol, abstol)) t0 y)
        in
        let explicit_table = Arkode.ButcherTable.{
            method_order = 3;
            embedding_order = 2;
            stages = 3;
            stage_values = RealArray2.of_lists [
                [ 0.0; 0.5; -1.0 ];
                [ 0.0; 0.0;  2.0 ];
                [ 0.0; 0.0;  0.0 ];
              ];
            stage_times = RealArray.of_list [ 0.0; 0.5; 1.0 ];
            coefficients = RealArray.of_list [
                             1.0 /. 6.0;  2.0 /. 3.0; 1.0 /. 6.0 ];
            bembed = Some (RealArray.of_list [ 0.0; 1.0; 0.0 ]);
          }
        in
        ARKStep.set_tables inner_arkode_mem ~explicit_table ();
        inner_arkode_mem

    | 1 ->  (* erk-3-3 fast solver (full problem) *)
        let inner_arkode_mem =
          ARKStep.(init (explicit (fn rpar)) (SStolerances (reltol, abstol)) t0 y)
        in
        let explicit_table = Arkode.ButcherTable.{
            method_order = 3;
            embedding_order = 2;
            stages = 3;
            stage_values = RealArray2.of_lists [
                [ 0.0; 0.5; -1.0 ];
                [ 0.0; 0.0;  2.0 ];
                [ 0.0; 0.0;  0.0 ];
              ];
            stage_times = RealArray.of_list [ 0.0; 0.5; 1.0 ];
            coefficients = RealArray.of_list [
                             1.0 /. 6.0;  2.0 /. 3.0; 1.0 /. 6.0 ];
            bembed = Some (RealArray.of_list [ 0.0; 1.0; 0.0 ]);
          }
        in
        ARKStep.set_tables inner_arkode_mem ~explicit_table ();
        inner_arkode_mem

    | 5 ->  (* erk-4-4 fast solver *)
        let inner_arkode_mem =
          ARKStep.(init (explicit (ff rpar)) default_tolerances t0 y)
        in
        let explicit_table = Arkode.ButcherTable.{
            method_order = 4;
            embedding_order = 0;
            stages = 4;
            stage_values = RealArray2.of_lists [
                [ 0.0; 0.5;  0.0; 0.0 ];
                [ 0.0; 0.0;  0.5; 0.0 ];
                [ 0.0; 0.0;  0.0; 1.0 ];
                [ 0.0; 0.0;  0.0; 0.0 ];
              ];
            stage_times = RealArray.of_list [ 0.0; 0.5; 0.5; 1.0 ];
            coefficients = RealArray.of_list [
                              1.0 /. 6.0;  1.0 /. 3.0; 1.0 /. 3.0; 1.0 /. 6.0 ];
            bembed = None;
          }
        in
        ARKStep.set_tables inner_arkode_mem ~explicit_table ();
        inner_arkode_mem

    | 2 ->  (* esdirk-3-3 fast solver (full problem) *)
        let a_f = Matrix.dense neq in
        let ls_f = LinearSolver.Direct.dense y a_f in
        let inner_arkode_mem =
          ARKStep.(init (implicit
                           ~lsolver:(Dls.solver ~jac:(jn rpar) ls_f)
                           (fn rpar))
                        (SStolerances (reltol, abstol))
                        t0 y)
        in
        let beta  = (sqrt 3.0) /. 6.0 +. 0.5 in
        let gamma = (-1.0 /. 8.0) *. (sqrt 3.0 +. 1.0) in
        let implicit_table = Arkode.ButcherTable.{
            method_order = 3;
            embedding_order = 0;
            stages = 3;
            stage_values = RealArray2.of_lists [
                [ 0.0 ; 4.0 *. gamma +. 2.0 *. beta        ; 0.5 -. beta -. gamma];
                [ 0.0 ; 1.0 -. 4.0 *. gamma -. 2.0 *. beta ; gamma ];
                [ 0.0 ; 0.0                                ; beta ];
              ];
            stage_times = RealArray.of_list [ 0.0; 1.0; 0.5 ];
            coefficients = RealArray.of_list [
                             1.0 /. 6.0;  1.0 /. 6.0; 2.0 /. 3.0 ];
            bembed = None;
          }
        in
        ARKStep.set_tables inner_arkode_mem ~implicit_table ();
        inner_arkode_mem

    | 3 | 4 -> (* no fast dynamics ('evolve' explicitly w/ erk-3-3) *)
        let inner_arkode_mem =
          ARKStep.(init (explicit (f0 rpar)) (SStolerances (reltol, abstol)) t0 y)
        in
        let explicit_table = Arkode.ButcherTable.{
            method_order = 3;
            embedding_order = 2;
            stages = 3;
            stage_values = RealArray2.of_lists [
                [ 0.0; 0.5; -1.0 ];
                [ 0.0; 0.0;  2.0 ];
                [ 0.0; 0.0;  0.0 ];
              ];
            stage_times = RealArray.of_list [ 0.0; 0.5; 1.0 ];
            coefficients = RealArray.of_list [
                             1.0 /. 6.0;  2.0 /. 3.0; 1.0 /. 6.0 ];
            bembed = Some (RealArray.of_list [ 0.0; 1.0; 0.0 ]);
          }
        in
        ARKStep.set_tables inner_arkode_mem ~explicit_table ();
        inner_arkode_mem

    | _ -> assert false
  in
  (* Set the fast step size *)
  ARKStep.set_fixed_step inner_arkode_mem (Some hf);

  (*
   * Create the slow integrator and set options
   *)

  (* Initialize the slow integrator. Specify the slow right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, the
     initial dependent variable vector y, and the fast integrator. *)
  let arkode_mem = match solve_type with
    | 0 ->  (* KW3 slow solver *)
        let arkode_mem =
          MRIStep.(init (InnerStepper.from_arkstep inner_arkode_mem)
                        default_tolerances
                        (fs rpar)
                        ~slowstep:hs
                        t0 y)
        in
        MRIStep.set_table_num arkode_mem Arkode.ButcherTable.Knoth_Wolke_3_3;
        arkode_mem

    | 3 ->  (* KW3 slow solver (full problem) *)
        let arkode_mem =
          MRIStep.(init (InnerStepper.from_arkstep inner_arkode_mem)
                        default_tolerances
                        (fn rpar)
                        ~slowstep:hs
                        t0 y)
        in
        MRIStep.set_table_num arkode_mem Arkode.ButcherTable.Knoth_Wolke_3_3;
        arkode_mem

    | 5 | 6 -> (* MRI-GARK-ERK45a slow solver *)
        let arkode_mem =
          MRIStep.(init (InnerStepper.from_arkstep inner_arkode_mem)
                        default_tolerances
                        (fs rpar)
                        ~slowstep:hs
                        t0 y)
        in
        MRIStep.(set_coupling arkode_mem Coupling.(load_table GARK_ERK45a));
        arkode_mem

    | 1 | 2 ->  (* no slow dynamics (use ERK-2-2) *)
        let arkode_mem =
          MRIStep.(init (InnerStepper.from_arkstep inner_arkode_mem)
                        default_tolerances
                        (f0 rpar)
                        ~slowstep:hs
                        t0 y)
        in
        let bt = Arkode.ButcherTable.{
            method_order = 2;
            embedding_order = 0;
            stages = 2;
            stage_values = RealArray2.of_arrays [|
                [| 0.0; 2.0 /. 3.0 |];
                [| 0.0;    0.0     |];
              |];
            stage_times = RealArray.of_list [ 0.0; 2.0 /. 3.0 ];
            coefficients = RealArray.of_list [ 0.25; 0.75 ];
            bembed = None;
          }
        in
        MRIStep.set_table arkode_mem 2 bt;
        arkode_mem

    | 4 ->  (* dirk-2 (trapezoidal), solve-decoupled slow solver *)
        let a_s = Matrix.dense neq in
        let ls_s = LinearSolver.Direct.dense y a_s in
        let arkode_mem =
          MRIStep.(init (InnerStepper.from_arkstep inner_arkode_mem)
                        (SStolerances (reltol, abstol))
                         ~lsolver:(Dls.solver ~jac:(jn rpar) ls_s)
                        (fn rpar)
                        ~slowstep:hs
                        t0 y)
        in
        MRIStep.(set_coupling arkode_mem Coupling.(load_table GARK_IRK21a));
        arkode_mem

    | 7 ->  (* MRI-GARK-ESDIRK34a, solve-decoupled slow solver *)
        let a_s = Matrix.dense neq in
        let ls_s = LinearSolver.Direct.dense y a_s in
        let arkode_mem =
          MRIStep.(init (InnerStepper.from_arkstep inner_arkode_mem)
                        (SStolerances (reltol, abstol))
                         ~lsolver:(Dls.solver ~jac:(js rpar) ls_s)
                        (fs rpar)
                        ~slowstep:hs
                        t0 y)
        in
        MRIStep.(set_coupling arkode_mem Coupling.(load_table GARK_ESDIRK34a));
        arkode_mem

    | _ -> assert false
  in

  (*
   * Integrate ODE
   *)

  (* Open output stream for results, output comment line *)
  let ufid = open_out "ark_kpr_mri_solution.txt" in
  fprintf ufid "# t u v uerr verr\n";

  (* output initial condition to disk *)
  fprintf ufid " %.16e %.16e %.16e %.16e %.16e\n"
          t0 ydata.{0} ydata.{1}
          (abs_float (ydata.{0} -. utrue rpar t0))
          (abs_float (ydata.{1} -. vtrue rpar t0));

  (* Main time-stepping loop: calls MRIStepEvolve to perform the
     integration, then prints results. Stops when the final time
     has been reached *)
  let uerr = ref 0.0 in
  let verr = ref 0.0 in
  let uerrtot = ref 0.0 in
  let verrtot = ref 0.0 in
  let errtot  = ref 0.0 in
  printf "        t           u           v       uerr      verr\n";
  printf "   ------------------------------------------------------\n";
  printf "  %10.6f  %10.6f  %10.6f  %.2e  %.2e\n"
         t0 ydata.{0} ydata.{1} !uerr !verr;

  let rec loop tout iout =
    if iout >= nt then ()
    else begin
      (* call integrator *)
      let t, _ = MRIStep.solve_normal arkode_mem tout y in

      (* access/print solution and error *)
      uerr := abs_float (ydata.{0} -. utrue rpar t);
      verr := abs_float (ydata.{1} -. vtrue rpar t);
      printf "  %10.6f  %10.6f  %10.6f  %.2e  %.2e\n"
             t ydata.{0} ydata.{1} !uerr !verr;
      fprintf ufid " %.16e %.16e %.16e %.16e %.16e\n"
              t ydata.{0} ydata.{1} !uerr !verr;
      uerrtot := !uerrtot +. !uerr *. !uerr;
      verrtot := !verrtot +. !verr *. !verr;
      errtot := !errtot +. !uerr *. !uerr +. !verr *. !verr;

      (* successful solve: update time *)
      loop (min tf (tout +. dTout)) (iout + 1)
    end;
  in
  loop (t0 +. dTout) 0;
  let uerrtot = sqrt (!uerrtot /. float nt) in
  let verrtot = sqrt (!verrtot /. float nt) in
  let errtot = sqrt (!errtot /. float nt /. 2.0) in
  printf "   ------------------------------------------------------\n";
  close_out ufid;

  (*
   * Finalize
   *)

  (* Get some slow integrator statistics *)
  let nsts = MRIStep.get_num_steps arkode_mem in
  let nfs = MRIStep.get_num_rhs_evals arkode_mem in

  (* Get some fast integrator statistics *)
  let nstf = ARKStep.get_num_steps inner_arkode_mem in
  let nff, _ = ARKStep.get_num_rhs_evals inner_arkode_mem in

  (* Print some final statistics *)
  printf "\nFinal Solver Statistics:\n";
  printf "   Steps: nsts = %d, nstf = %d\n" nsts nstf;
  printf "   u error = %.3e, v error = %.3e, total error = %.3e\n"
         uerrtot verrtot errtot;
  printf "   Total RHS evals:  Fs = %d,  Ff = %d\n" nfs nff;

  (* Get/print slow integrator decoupled implicit solver statistics *)
  if solve_type = 4 || solve_type = 7 then begin
    let nnis, nncs = MRIStep.get_nonlin_solv_stats arkode_mem in
    let njes = MRIStep.Dls.get_num_jac_evals arkode_mem in
    printf "   Slow Newton iters = %d\n" nnis;
    printf "   Slow Newton conv fails = %d\n" nncs;
    printf "   Slow Jacobian evals = %d\n" njes
  end;

  (* Get/print fast integrator implicit solver statistics *)
  if solve_type = 2 then begin
    let nnif, nncf = ARKStep.get_nonlin_solv_stats inner_arkode_mem in
    let njef = ARKStep.Dls.get_num_jac_evals inner_arkode_mem in
    printf "   Fast Newton iters = %d\n" nnif;
    printf "   Fast Newton conv fails = %d\n" nncf;
    printf "   Fast Jacobian evals = %d\n" njef
  end

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
  for i = 1 to reps do
    main ();
    if gc_each_rep then Gc.compact ()
  done;
  if gc_at_end then Gc.compact ()

