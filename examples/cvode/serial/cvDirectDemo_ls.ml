(*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007/10/25 20:03:29 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Feb 2011.
 * -----------------------------------------------------------------
 * Demonstration program for CVODE - direct linear solvers.
 * Two separate problems are solved using both the CV_ADAMS and CV_BDF
 * linear multistep methods in combination with CV_FUNCTIONAL and
 * CV_NEWTON iterations:
 *
 * Problem 1: Van der Pol oscillator
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.
 * This second-order ODE is converted to a first-order system by
 * defining y0 = x and y1 = xdot.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) dense, user-supplied, (2) dense, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * Problem 2: ydot = A * y, where A is a banded lower triangular
 * matrix derived from 2-D advection PDE.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) band, user-supplied, (2) band, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * For each problem, in the series of eight runs, CVodeInit is
 * called only once, for the first run, whereas CVodeReInit is
 * called for each of the remaining seven runs.
 *
 * Notes: This program demonstrates the usage of the sequential
 * macros NV_Ith_S, NV_DATA_S, DENSE_ELEM, BAND_COL, and
 * BAND_COL_ELEM. The NV_Ith_S macro is used to reference the
 * components of an N_Vector. It works for any size N=NEQ, but
 * due to efficiency concerns it should only by used when the
 * problem size is small. The Problem 1 right hand side and
 * Jacobian functions f1 and Jac1 both use NV_Ith_S. The NV_DATA_S
 * macro gives the user access to the memory used for the component
 * storage of an N_Vector. In the sequential case, the user may
 * assume that this is one contiguous array of reals. The NV_DATA_S
 * macro gives a more efficient means (than the NV_Ith_S macro) to
 * access the components of an N_Vector and should be used when the
 * problem size is large. The Problem 2 right hand side function f2
 * uses the NV_DATA_S macro. The DENSE_ELEM macro used in Jac1
 * gives access to an element of a dense matrix of type DlsMat.
 * It should be used only when the problem size is small (the size
 * of a DlsMat is NEQ x NEQ) due to efficiency concerns. For
 * larger problem sizes, the macro DENSE_COL can be used in order
 * to work directly with a column of a DlsMat. The BAND_COL and
 * BAND_COL_ELEM allow efficient columnwise access to the elements
 * of a band matrix of type DlsMat. These macros are used in the
 * Jac2 function.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module Roots = Sundials.Roots
module Densematrix = Dls.DenseMatrix
module Bandmatrix = Dls.BandMatrix
let unwrap = Nvector.unwrap

let printf = Printf.printf

let bandset bm i j v = Dls.BandMatrix.set bm i j v

(* Shared Problem Constants *)

let atol  = 1.0e-6
let rtol  = 0.0

let zero    = 0.0
let one     = 1.0
let two     = 2.0
let thirty  = 30.0

(* Problem #1 Constants *)

let p1_neq         = 2
let p1_eta         = 3.0
let p1_nout        = 4
let p1_t0          = 0.0
let p1_t1          = 1.39283880203
let p1_dtout       = 2.214773875
let p1_tol_factor  = 1.0e4

(* Problem #2 Constants *)

let p2_meshx       = 5
let p2_meshy       = 5
let p2_neq         = p2_meshx * p2_meshy
let p2_alph1       = 1.0
let p2_alph2       = 1.0
let p2_nout        = 5
let p2_ml          = 5
let p2_mu          = 0
let p2_t0          = 0.0
let p2_t1          = 0.01
let p2_tout_mult   = 10.0
let p2_tol_factor  = 1.0e3

(* Linear Solver Options *)

type miter =
  | Func
  | Dense_User
  | Dense_DQ
  | Diag
  | Band_User
  | Band_DQ

(* Private Helper Functions *)

let sqr x = x *. x

(* Functions Called by the Solver *)

let jac1 { Cvode.jac_y = y } j =
  let y0 = y.{0} in
  let y1 = y.{1} in
  (* previously calls to DENSE_ELEM: *)
  let set = Dls.DenseMatrix.set j in
  set 0 1 one;
  set 1 0 (-. two *. p1_eta *. y0 *. y1 -. one);
  set 1 1 (p1_eta *. (one -. sqr y0))

let jac2 {Cvode.mupper=mu; Cvode.mlower=ml} arg jac =
  (*
     The components of f(t,y) which depend on y    are
                                               i,j
     f    , f      , and f      : 
      i,j    i+1,j        i,j+1

     f    = -2 y    + alpha1 * y      + alpha2 * y
      i,j       i,j             i-1,j             i,j-1

     f      = -2 y      + alpha1 * y    + alpha2 * y
      i+1,j       i+1,j             i,j             i+1,j-1

     f      = -2 y      + alpha1 * y        + alpha2 * y
      i,j+1       i,j+1             i-1,j+1             i,j
  *)
  for j = 0 to p2_meshy - 1 do
    for i = 0 to p2_meshx - 1 do
      let k = i + j * p2_meshx in
      bandset jac k k (-. two);
      if (i != p2_meshx - 1) then bandset jac (k + 1) k p2_alph1;
      if (j != p2_meshy - 1) then bandset jac (k + p2_meshx) k p2_alph2
    done
  done

let prepare_next_run cvode_mem lmm miter mu ml t y =
  printf "\n\n-------------------------------------------------------------";
  
  printf "\n\nLinear Multistep Method : ";
  (match lmm with
   | Cvode.Adams -> printf "ADAMS\n"
   | Cvode.BDF -> printf "BDF\n");
  
  printf "Iteration               : ";

  (* Func comes first and is set via init; for all other cases, reinit.  *)
  if miter = Func
  then printf "FUNCTIONAL\n"
  else begin
    printf "NEWTON\n";
    printf "Linear Solver           : ";

    match miter with
    | Dense_User -> begin
          printf "Dense, User-Supplied Jacobian\n";
          Cvode.(reinit cvode_mem t y
                        ~iter_type:(Newton (Dls.dense ~jac:jac1 ())))
        end

    | Dense_DQ -> begin
          printf("Dense, Difference Quotient Jacobian\n");
          Cvode.(reinit cvode_mem t y
                        ~iter_type:(Newton (Dls.dense ())))
        end

    | Diag -> begin
          printf("Diagonal Jacobian\n");
          Cvode.(reinit cvode_mem t y ~iter_type:(Newton Diag.solver))
        end

    | Band_User -> begin
          printf("Band, User-Supplied Jacobian\n");
          Cvode.(reinit cvode_mem t y
            ~iter_type:
              (Newton (Dls.band { mupper = mu; mlower = ml } ~jac:jac2)))
        end

    | Band_DQ -> begin
          printf("Band, Difference Quotient Jacobian\n");
          Cvode.(reinit cvode_mem t y
            ~iter_type:(Newton (Dls.band { mupper = mu; mlower = ml })))
        end

    | Func -> assert false
  end

let print_final_stats cvode_mem miter ero =
  let open Cvode in
  let (lenrw, leniw) = get_work_space cvode_mem
  and nst            = get_num_steps cvode_mem
  and nfe            = get_num_rhs_evals cvode_mem
  and nsetups        = get_num_lin_solv_setups cvode_mem
  and netf           = get_num_err_test_fails cvode_mem
  and nni            = get_num_nonlin_solv_iters cvode_mem
  and ncfn           = get_num_nonlin_solv_conv_fails cvode_mem
  in
  printf "\n Final statistics for this run:\n\n";
  printf " CVode real workspace length              = %4d \n"   lenrw;
  printf " CVode integer workspace length           = %4d \n"   leniw;
  printf " Number of steps                          = %4d \n"   nst;
  printf " Number of f-s                            = %4d \n"   nfe;
  printf " Number of setups                         = %4d \n"   nsetups;
  printf " Number of nonlinear iterations           = %4d \n"   nni;
  printf " Number of nonlinear convergence failures = %4d \n"   ncfn;
  printf " Number of error test failures            = %4d \n\n" netf;
  
  if miter != Func then begin
    let nje, nfeLS, (lenrwLS, leniwLS) =
      match miter with
      | Dense_User | Dense_DQ ->
          Dls.get_num_jac_evals cvode_mem,
          Dls.get_num_rhs_evals cvode_mem,
          Dls.get_work_space    cvode_mem

      | Band_User | Band_DQ ->
          Dls.get_num_jac_evals cvode_mem,
          Dls.get_num_rhs_evals cvode_mem,
          Dls.get_work_space    cvode_mem

      | Diag ->
          nsetups,
          Diag.get_num_rhs_evals cvode_mem,
          Diag.get_work_space cvode_mem

      | Func -> assert false
    in
    printf " Linear solver real workspace length      = %4d \n" lenrwLS;
    printf " Linear solver integer workspace length   = %4d \n" leniwLS;
    printf " Number of Jacobian evaluations           = %4d  \n" nje;
    printf " Number of f evals. in linear solver      = %4d \n\n" nfeLS
  end;

  printf " Error overrun = %.3f \n" ero


(* Implementation *)

let print_intro1 () =
  printf "Demonstration program for CVODE package - direct linear solvers\n";
  printf "\n\n";
  printf "Problem 1: Van der Pol oscillator\n";
  printf " xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0\n";
  printf " neq = %d,  reltol = %.2g,  abstol = %.2g" p1_neq rtol atol

let print_header1 () =
  printf "\n     t           x              xdot         qu     hu \n"

let print_output1 t y0 y1 qu hu =
  printf "%10.5f    %12.5e   %12.5e   %2d    %6.4e\n" t y0 y1 qu hu

let print_err_output tol_factor =
  printf "\n\n Error exceeds %g * tolerance \n\n" tol_factor

let f1 t (y : RealArray.t) (ydot : RealArray.t) =
  let y0 = y.{0} in
  let y1 = y.{1} in
  ydot.{0} <- y1;
  ydot.{1} <- (one -. sqr y0) *. p1_eta *. y1 -. y0

let snd_true (x, _) = (x, true)

let problem1 () =
  let nerr = ref 0 in
  let y = Nvector_serial.make p1_neq 0.0 in
  let ydata = unwrap y in
  let init_y () = (ydata.{0} <- two; ydata.{1} <- zero) in
  print_intro1 ();

  let run cvode_mem lmm miter =
    begin
      let ero = ref zero in

      (* Func comes first and is set via init; for all other cases, reinit.  *)
      if miter <> Func then init_y ();

      prepare_next_run cvode_mem lmm miter 0 0 p1_t0 y;

      print_header1 ();

      let tout = ref p1_t1 in
      for iout = 1 to p1_nout do
        let (t, success) =
          try
            snd_true (Cvode.solve_normal cvode_mem !tout y)
          with _ -> (incr nerr; (!tout, false))
        in

        let qu =
          try
            Cvode.get_last_order cvode_mem;
          with _ -> (incr nerr; 0)
        and hu =
          try
            Cvode.get_last_step cvode_mem
          with _ -> (incr nerr; 0.0)
        in

        print_output1 t ydata.{0} ydata.{1} qu hu;

        if success && (iout mod 2 = 0) then begin
          let er = abs_float ydata.{0} /. atol in
          ero := max er !ero;
          if er > p1_tol_factor then
            (incr nerr; print_err_output p1_tol_factor)
        end;

        tout := !tout +. p1_dtout
      done;
      print_final_stats cvode_mem miter !ero
    end
  in

  let run_tests lmm =
    init_y ();
    let cvode_mem = Cvode.(init lmm Functional
                                (SStolerances (rtol, atol)) f1 p1_t0 y)
    in
    List.iter (run cvode_mem lmm) [ Func; Dense_User; Dense_DQ; Diag]
  in

  run_tests Cvode.Adams;
  run_tests Cvode.BDF;

  !nerr

let print_intro2 () =
  printf "\n\n-------------------------------------------------------------";
  printf "\n-------------------------------------------------------------";
  printf "\n\nProblem 2: ydot = A * y, where A is a banded lower\n";
  printf "triangular matrix derived from 2-D advection PDE\n\n";
  printf " neq = %d, ml = %d, mu = %d\n" p2_neq p2_ml p2_mu;
  printf " itol = %s, reltol = %.2g, abstol = %.2g" "CV_SS" rtol atol;
  match Sundials.sundials_version with
  | 2,5,_ -> printf "\n      t        max.err      qu     hu \n"
  | _ -> ()

let print_header2 () =
  printf "\n      t        max.err      qu     hu \n"

let print_output2 t erm qu hu =
  printf "%10.3f  %12.4e   %2d   %12.4e\n" t erm qu hu

let f2 t (ydata : RealArray.t) (dydata : RealArray.t) =
  (*
     Excluding boundaries, 

     ydot    = f    = -2 y    + alpha1 * y      + alpha2 * y
         i,j    i,j       i,j             i-1,j             i,j-1
  *)
  for j = 0 to p2_meshy - 1 do
    for i = 0 to p2_meshx - 1 do
      let k = i + j * p2_meshx in 
      let d = -. two *. ydata.{k} in
      let d = d +. (if (i == 0) then 0.0 else p2_alph1 *. ydata.{k - 1}) in
      let d = d +. (if (j == 0) then 0.0 else p2_alph2 *. ydata.{k - p2_meshx})
      in
      dydata.{k} <- d
    done
  done

let max_error (ydata : RealArray.t) t =
  if t = zero then zero
  else
    let ex = if (t <= thirty) then exp(-. two *. t) else zero in
    
    let max_error = ref zero in

    let jfact_inv = ref one in
    for j = 0 to p2_meshy - 1 do

      let ifact_inv = ref one in
      for i = 0 to p2_meshx - 1 do
        let k  = i + j * p2_meshx in
        let yt = t ** float (i + j) *. ex *. !ifact_inv *. !jfact_inv in
        let er = abs_float (ydata.{k} -. yt) in
        max_error := max er !max_error;
        ifact_inv := !ifact_inv /. float (i + 1)
      done;

      jfact_inv := !jfact_inv /. float (j + 1)
    done;
    !max_error

let problem2 () =
  let nerr = ref 0 in
  let y = Nvector_serial.make p2_neq 0.0 in
  let ydata = unwrap y in
  let init_y () = (RealArray.fill ydata zero; ydata.{0} <- one) in
  print_intro2 ();

  let run cvode_mem lmm miter =
    begin
      let ero = ref zero in

      (* Func comes first and is set via init; for all other cases, reinit.  *)
      if miter <> Func then init_y ();
      prepare_next_run cvode_mem lmm miter p2_mu p2_ml p2_t0 y;
      print_header2 ();

      let tout = ref p2_t1 in
      for iout = 1 to p2_nout do
        let (t, success) =
          try
            snd_true (Cvode.solve_normal cvode_mem !tout y)
          with _ -> (incr nerr; (!tout, false))
        in

        let erm = max_error ydata t in
        let qu =
          try
            Cvode.get_last_order cvode_mem;
          with _ -> (incr nerr; 0)
        and hu =
          try
            Cvode.get_last_step cvode_mem
          with _ -> (incr nerr; 0.0)
        in

        print_output2 t erm qu hu;

        if success then begin
          let er = erm /. atol in
          ero := max er !ero;
          if er > p2_tol_factor then
            (incr nerr; print_err_output p2_tol_factor)
        end;

        tout := !tout *. p2_tout_mult
      done;
      print_final_stats cvode_mem miter !ero;
    end
  in

  let run_tests lmm =
    init_y ();
    let cvode_mem = Cvode.(init lmm Functional
                                (SStolerances (rtol, atol)) f2 p2_t0 y)
    in
    List.iter (run cvode_mem lmm) [ Func; Diag; Band_User; Band_DQ]
  in

  run_tests Cvode.Adams;
  run_tests Cvode.BDF;

  !nerr

let print_err_info nerr =
  printf "\n\n-------------------------------------------------------------";
  printf "\n-------------------------------------------------------------";
  printf "\n\n Number of errors encountered = %d \n" nerr

let main () =
  let nerr1 = problem1 () in
  let nerr2 = problem2 () in
  print_err_info (nerr1 + nerr2)

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
