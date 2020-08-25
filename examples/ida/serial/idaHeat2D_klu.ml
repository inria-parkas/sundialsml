(*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2009/09/30 23:25:59 $
 * -----------------------------------------------------------------
 * Programmer(s): Chris Nguyen @ LLNL
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, serial, sparse.
 * Based on idaHeat2D_bnd.c and idaRoberts_klu.c
 * -----------------------------------------------------------------
 * OCaml port: Jun Inoue, Inria, Aug 2014.
 *             Timothy Bourke, Inria, Jan 2016.
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, serial, banded.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the band solver IDABand, and IDACalcIC.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform mgrid x mgrid
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = mgrid^2. Here mgrid = 10.
 *
 * The system is solved with IDA using the banded linear system
 * solver and default difference-quotient Jacobian.
 * For purposes of illustration,
 * IDACalcIC is called to compute correct values at the boundary,
 * given incorrect values as input initial guesses. The constraints
 * u >= 0 are posed for all components. Output is taken at
 * t = 0, .01, .02, .04, ..., 10.24. (Output at t = 0 is for
 * IDACalcIC cost statistics only.)
 * -----------------------------------------------------------------
 *)

open Sundials

module Dls = Ida.Dls

let printf = Printf.printf
let vmax_norm = Nvector_serial.DataOps.n_vmaxnorm
let vscale = Nvector_serial.DataOps.n_vscale

(* Problem Constants *)
let nout = 11
let mgrid = 10                          (* mesh grid size *)
let neq = mgrid * mgrid
let bval = 0.0
let total = 4*mgrid+8*(mgrid-2)+(mgrid-4)*(mgrid+4*(mgrid-2))
            (* total num of nonzero elements *)

type user_data = { mm : int; dx : float; coeff : float }

(* Jacobian matrix setup for mgrid=3  *)
let jac_heat3 { Ida.jac_y = (yval : RealArray.t);
                Ida.jac_coef = cj }
              jacmat =
  let set_col = Matrix.Sparse.set_col jacmat in
  let set = Matrix.Sparse.set jacmat in
  let dx = 1.0/.float(mgrid - 1) in
  let beta = 4.0/.(dx*.dx) +. cj in

  (* initialize Jacobian matrix *)
  Matrix.Sparse.set_to_zero jacmat;

  (* set up number of elements in each column *)
  set_col 0 0;
  set_col 1 1;
  set_col 2 3;
  set_col 3 4;
  set_col 4 6;
  set_col 5 7;
  set_col 6 9;
  set_col 7 10;
  set_col 8 12;
  set_col 9 13;

  (* set up data and row values stored *)
  set 0 0 1.0;
  set 1 1 1.0;
  set 2 4 (-1.0/.(dx*.dx));
  set 3 2 1.0;
  set 4 3 1.0;
  set 5 4 (-1.0/.(dx*.dx));
  set 6 4 (beta);
  set 7 4 (-1.0/.(dx*.dx));
  set 8 5 1.0;
  set 9 6 1.0;
  set 10 4 (-1.0/.(dx*.dx));
  set 11 7 1.0;
  set 12 8 1.0

(* Jacobian matrix setup for mgrid>=4  *)
let jac_heat { Ida.jac_y = (yval : RealArray.t);
               Ida.jac_coef = cj }
              jacmat =
  let get_col = Matrix.Sparse.get_col jacmat in
  let set_col = Matrix.Sparse.set_col jacmat in
  let set_data = Matrix.Sparse.set_data jacmat in
  let set_rowval = Matrix.Sparse.set_rowval jacmat in
  let dx = 1.0/.float(mgrid - 1) in
  let beta = 4.0/.(dx*.dx) +. cj in

  let loop ~from:start ~upto:limit ?inc:(increment=1) f =
    let rec go i = if i < limit then (f i; go (i + increment)) else () in
    go start
  in

  (* initialize Jacobian matrix *)
  Matrix.Sparse.set_to_zero jacmat;

  (* ---- set up number of elements in each column ---- *)

  (**** first column block ****)
  set_col (0) (0);
  set_col (1) (1);
  (* count by twos in the middle  *)
  for i=2 to mgrid-1 do
    set_col i ((get_col (i-1))+2)
  done;
  set_col mgrid (2*mgrid-2);

  (**** second column block ****)
  set_col (mgrid+1) (2*mgrid);
  set_col (mgrid+2) (2*mgrid+3);
  (* count by fours in the middle *)
  for i=0 to mgrid-5 do
    set_col (mgrid+3+i) ((get_col (mgrid+3+i-1))+4)
  done;
  set_col (2*mgrid-1) (2*mgrid+4*(mgrid-2)-2);
  set_col (2*mgrid)   (2*mgrid+4*(mgrid-2));

  (**** repeated (mgrid-4 times) middle column blocks ****)
  for i=0 to mgrid-5 do
    let repeat = i * mgrid in (* shift that accounts for accumulated
                                 number of columns *)

    set_col (2*mgrid+1+repeat)   ((get_col (2*mgrid+1+repeat-1))+2);
    set_col (2*mgrid+1+repeat+1) ((get_col (2*mgrid+1+repeat))+4);

    (* count by fives in the middle *)
    for j=0 to mgrid-5 do
      set_col (2*mgrid+1+repeat+2+j) (get_col (2*mgrid+1+repeat+1+j) + 5)
    done;

    set_col (2*mgrid+1+repeat+(mgrid-4)+2)
            ((get_col (2*mgrid+1+repeat+(mgrid-4)+1))+4);
    set_col (2*mgrid+1+repeat+(mgrid-4)+3)
            ((get_col (2*mgrid+1+repeat+(mgrid-4)+2))+2)
  done;

  (**** last-1 column block ****)
  set_col (mgrid*mgrid-2*mgrid+1) (total-2*mgrid-4*(mgrid-2)+2);
  set_col (mgrid*mgrid-2*mgrid+2) (total-2*mgrid-4*(mgrid-2)+5);
  (* count by fours in the middle *)
  for i=0 to mgrid-5 do
    set_col (mgrid*mgrid-2*mgrid+3+i) ((get_col (mgrid*mgrid-2*mgrid+3+i-1))+4)
  done;
  set_col (mgrid*mgrid-mgrid-1) (total-2*mgrid);
  set_col (mgrid*mgrid-mgrid) (total-2*mgrid+2);

  (**** last column block ****)
  set_col (mgrid*mgrid-mgrid+1) (total-mgrid-(mgrid-2)+1);
  (* count by twos in the middle *)
  for i=0 to mgrid-3 do
    set_col (mgrid*mgrid-mgrid+2+i) ((get_col (mgrid*mgrid-mgrid+2+i-1))+2)
  done;
  set_col (mgrid*mgrid-1) (total-1);
  set_col (mgrid*mgrid)   (total);

  (* ---- set up data stored ---- *)

  (**** first column block ****)
  set_data 0 1.0;
  (* alternating pattern in data, separate loop for each pattern  *)
  loop ~from:1 ~upto:(mgrid+(mgrid-2)) ~inc:2
    (fun i -> set_data i 1.0);
  loop ~from:2 ~upto:(mgrid+(mgrid-2)-1) ~inc:2
    (fun i -> set_data i (-1.0/.(dx*.dx)));

  (**** second column block ****)
  set_data (mgrid+mgrid-2) (1.0);
  set_data (mgrid+mgrid-1) (-1.0/.(dx*.dx));
  set_data (mgrid+mgrid) beta;
  set_data (mgrid+mgrid+1) (-1.0/.(dx*.dx));
  set_data (mgrid+mgrid+2) (-1.0/.(dx*.dx));
  (* middle data elements *)
  for i=0 to mgrid-4-1 do
    set_data (mgrid+mgrid+3+4*i) (-1.0/.(dx*.dx));
    set_data (mgrid+mgrid+4+4*i) (beta);
    set_data (mgrid+mgrid+5+4*i) (-1.0/.(dx*.dx));
    set_data (mgrid+mgrid+6+4*i) (-1.0/.(dx*.dx))
  done;
  set_data (2*mgrid+4*(mgrid-2)-5) (-1.0/.(dx*.dx));
  set_data (2*mgrid+4*(mgrid-2)-4) (beta);
  set_data (2*mgrid+4*(mgrid-2)-3) (-1.0/.(dx*.dx));
  set_data (2*mgrid+4*(mgrid-2)-2) (-1.0/.(dx*.dx));
  set_data (2*mgrid+4*(mgrid-2)-1) (1.0);

  (**** repeated (mgrid-4 times) middle column blocks ****)
  for i=0 to mgrid-4-1 do
    (* shift that accounts for accumulated columns and elements *)
    let repeat = i * (mgrid+4*(mgrid-2)) in

    set_data (2*mgrid+4*(mgrid-2)+repeat)   (1.0);
    set_data (2*mgrid+4*(mgrid-2)+repeat+1) (-1.0/.(dx*.dx));

    set_data (2*mgrid+4*(mgrid-2)+repeat+2) (-1.0/.(dx*.dx));
    set_data (2*mgrid+4*(mgrid-2)+repeat+3) (beta);
    set_data (2*mgrid+4*(mgrid-2)+repeat+4) (-1.0/.(dx*.dx));
    set_data (2*mgrid+4*(mgrid-2)+repeat+5) (-1.0/.(dx*.dx));

    (* 5 in 5*j chosen since there are 5 elements in each column *)
    (* this column loops mgrid-4 times within the outer loop *)
    for j=0 to mgrid-4-1 do
      set_data (2*mgrid+4*(mgrid-2)+repeat+6+5*j)  (-1.0/.(dx*.dx));
      set_data (2*mgrid+4*(mgrid-2)+repeat+7+5*j)  (-1.0/.(dx*.dx));
      set_data (2*mgrid+4*(mgrid-2)+repeat+8+5*j)  beta;
      set_data (2*mgrid+4*(mgrid-2)+repeat+9+5*j)  (-1.0/.(dx*.dx));
      set_data (2*mgrid+4*(mgrid-2)+repeat+10+5*j) (-1.0/.(dx*.dx));
    done;

    set_data (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+6) (-1.0/.(dx*.dx));
    set_data (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+7) (-1.0/.(dx*.dx));
    set_data (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+8) (beta);
    set_data (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+9) (-1.0/.(dx*.dx));

    set_data (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+10) (-1.0/.(dx*.dx));
    set_data (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+11) (1.0)
  done;

  (**** last-1 column block ****)
  set_data (total-6*(mgrid-2)-4) (1.0);
  set_data (total-6*(mgrid-2)-3) (-1.0/.(dx*.dx));
  set_data (total-6*(mgrid-2)-2) (-1.0/.(dx*.dx));
  set_data (total-6*(mgrid-2)-1) (beta);
  set_data (total-6*(mgrid-2)  ) (-1.0/.(dx*.dx));

  (* middle data elements *)
  for i=0 to mgrid-4-1 do
    set_data (total-6*(mgrid-2)+1+4*i) (-1.0/.(dx*.dx));
    set_data (total-6*(mgrid-2)+2+4*i) (-1.0/.(dx*.dx));
    set_data (total-6*(mgrid-2)+3+4*i) (beta);
    set_data (total-6*(mgrid-2)+4+4*i) (-1.0/.(dx*.dx))
  done;
  set_data (total-2*(mgrid-2)-7) (-1.0/.(dx*.dx));
  set_data (total-2*(mgrid-2)-6) (-1.0/.(dx*.dx));
  set_data (total-2*(mgrid-2)-5) (beta);
  set_data (total-2*(mgrid-2)-4) (-1.0/.(dx*.dx));
  set_data (total-2*(mgrid-2)-3) (1.0);

  (**** last column block ****)
  set_data (total-2*(mgrid-2)-2) (1.0);
  (* alternating pattern in data, separate loop for each pattern  *)
  loop ~from:(total-2*(mgrid-2)-1) ~upto:(total-2) ~inc:2
    (fun i -> set_data (i) (-1.0/.(dx*.dx)));
  loop ~from:(total-2*(mgrid-2))   ~upto:(total-1) ~inc:2
    (fun i -> set_data (i) 1.0);
  set_data (total-1) (1.0);

  (* ---- row values ---- *)

  (**** first block ****)
  set_rowval 0 0;
  (* alternating pattern in data, separate loop for each pattern *)
  loop ~from:1 ~upto:(mgrid+(mgrid-2)) ~inc:2
    (fun i -> set_rowval (i) ((i+1)/2));
  loop ~from:2 ~upto:(mgrid+(mgrid-2)-1) ~inc:2
    (fun i -> set_rowval (i) (i/2+mgrid)); (* i+1 unnecessary here *)

  (**** second column block ****)
  set_rowval (mgrid+mgrid-2) (mgrid);
  set_rowval (mgrid+mgrid-1) (mgrid+1);
  set_rowval (mgrid+mgrid)   (mgrid+1);
  set_rowval (mgrid+mgrid+1) (mgrid+2);
  set_rowval (mgrid+mgrid+2) (2*mgrid+1);
  (* middle row values *)
  for i=0 to mgrid-4 do
    set_rowval (mgrid+mgrid+3+4*i) (mgrid+1+i);
    set_rowval (mgrid+mgrid+4+4*i) (mgrid+2+i);
    set_rowval (mgrid+mgrid+5+4*i) (mgrid+3+i);
    set_rowval (mgrid+mgrid+6+4*i) (2*mgrid+2+i)
  done;
  set_rowval (2*mgrid+4*(mgrid-2)-5) (mgrid+(mgrid-2)-1);
  set_rowval (2*mgrid+4*(mgrid-2)-4) (mgrid+(mgrid-2)); (* starting from here,
                                                           add two diag
                                                           patterns *)
  set_rowval (2*mgrid+4*(mgrid-2)-3) (2*mgrid+(mgrid-2));
  set_rowval (2*mgrid+4*(mgrid-2)-2) (mgrid+(mgrid-2));
  set_rowval (2*mgrid+4*(mgrid-2)-1) (mgrid+(mgrid-2)+1);

  (**** repeated (mgrid-4 times) middle column blocks ****)
  for i=0 to mgrid-4-1 do
    (* shift that accounts for accumulated columns and elements *)
    let repeat = i*(mgrid+4*(mgrid-2)) in
    set_rowval (2*mgrid+4*(mgrid-2)+repeat)
               (mgrid+(mgrid-2)+2+mgrid*i);
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+1)
               (mgrid+(mgrid-2)+2+mgrid*i+1);

    set_rowval (2*mgrid+4*(mgrid-2)+repeat+2)
               (mgrid+(mgrid-2)+2+mgrid*i+1-mgrid);
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+3)
               (mgrid+(mgrid-2)+2+mgrid*i+1);
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+4)
               (mgrid+(mgrid-2)+2+mgrid*i+2); (* *this *)
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+5)
               (mgrid+(mgrid-2)+2+mgrid*i+1+mgrid);

    (* 5 in 5*j chosen since there are 5 elements in each column *)
    (* column repeats mgrid-4 times within the outer loop *)
    for j=0 to mgrid-4-1 do
      set_rowval (2*mgrid+4*(mgrid-2)+repeat+6+5*j)
                 (mgrid+(mgrid-2)+2+mgrid*i+1-mgrid+1+j);
      set_rowval (2*mgrid+4*(mgrid-2)+repeat+7+5*j)
                 (mgrid+(mgrid-2)+2+mgrid*i+1+j);
      set_rowval (2*mgrid+4*(mgrid-2)+repeat+8+5*j)
                 (mgrid+(mgrid-2)+2+mgrid*i+2+j);
      set_rowval (2*mgrid+4*(mgrid-2)+repeat+9+5*j)
                 (mgrid+(mgrid-2)+2+mgrid*i+2+1+j);
      set_rowval (2*mgrid+4*(mgrid-2)+repeat+10+5*j)
                 (mgrid+(mgrid-2)+2+mgrid*i+1+mgrid+1+j)
    done;

    set_rowval (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+6)
               (mgrid+(mgrid-2)+2+mgrid*i-2);
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+7)
               (mgrid+(mgrid-2)+2+mgrid*i-2+mgrid-1);
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+8)
               (mgrid+(mgrid-2)+2+mgrid*i-2+mgrid); (* *this+mgrid *)
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+9)
               (mgrid+(mgrid-2)+2+mgrid*i-2+2*mgrid);

    set_rowval (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+10)
               (mgrid+(mgrid-2)+2+mgrid*i-2+mgrid);
    set_rowval (2*mgrid+4*(mgrid-2)+repeat+(mgrid-4)*5+11)
               (mgrid+(mgrid-2)+2+mgrid*i-2+mgrid+1)
  done;

  (**** last-1 column block ****)
  set_rowval (total-6*(mgrid-2)-4)
             (mgrid*mgrid-1-2*(mgrid-1)-1);
  set_rowval (total-6*(mgrid-2)-3)
             (mgrid*mgrid-1-2*(mgrid-1)); (* starting with this as base *)
  set_rowval (total-6*(mgrid-2)-2)
             (mgrid*mgrid-1-2*(mgrid-1)-mgrid);
  set_rowval (total-6*(mgrid-2)-1)
             (mgrid*mgrid-1-2*(mgrid-1));
  set_rowval (total-6*(mgrid-2)  )
             (mgrid*mgrid-1-2*(mgrid-1)+1);
  (* middle row values *)
  for i=0 to mgrid-4-1 do
    set_rowval (total-6*(mgrid-2)+1+4*i)
               (mgrid*mgrid-1-2*(mgrid-1)-mgrid+1+i);
    set_rowval (total-6*(mgrid-2)+2+4*i)
               (mgrid*mgrid-1-2*(mgrid-1)+i);
    set_rowval (total-6*(mgrid-2)+3+4*i)
               (mgrid*mgrid-1-2*(mgrid-1)+1+i); (*copied above*)
    set_rowval (total-6*(mgrid-2)+4+4*i)
               (mgrid*mgrid-1-2*(mgrid-1)+2+i);
  done;
  set_rowval (total-2*(mgrid-2)-7) (mgrid*mgrid-2*mgrid-2);
  set_rowval (total-2*(mgrid-2)-6) (mgrid*mgrid-mgrid-3);
  set_rowval (total-2*(mgrid-2)-5) (mgrid*mgrid-mgrid-2);
  set_rowval (total-2*(mgrid-2)-4) (mgrid*mgrid-mgrid-2);
  set_rowval (total-2*(mgrid-2)-3) (mgrid*mgrid-mgrid-1);

  (* last column block *)
  set_rowval (total-2*(mgrid-2)-2) (mgrid*mgrid-mgrid);
  (* alternating pattern in data, separate loop for each pattern  *)
  for i=0 to mgrid-2-1 do
    set_rowval (total-2*(mgrid-2)-1+2*i) (mgrid*mgrid-2*mgrid+1+i);
    set_rowval (total-2*(mgrid-2)  +2*i) (mgrid*mgrid-mgrid+1+i);
  done;
  set_rowval (total-1) (mgrid*mgrid-1)

  (*; PrintSparseMat(JacMat)*)

(*
 * heatres: heat equation system residual function
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 *)
let heatres t (u : RealArray.t) (u' : RealArray.t) resval data =
  let mm = data.mm
  and coeff = data.coeff
  in
  (* Initialize resval to u, to take care of boundary equations. *)
  RealArray.fill resval 0.0;

  (* Loop over interior points; set res = u' - (central difference). *)
  for j = 1 to mm-2 do
    let offset = mm*j in
    for i = 1 to mm-2 do
      let loc = offset + i in
      resval.{loc} <-
        u'.{loc} -.
        coeff *.
        (u.{loc-1} +. u.{loc+1} +. u.{loc-mm} +. u.{loc+mm} -. 4. *. u.{loc})
    done
  done

let set_initial_profile data u u' id res =
  (* Initialize id to differential. *)
  RealArray.fill id Ida.VarId.differential;

  let mm = data.mm in
  let mm1 = mm - 1 in

  (* Initialize u on all grid points. *)
  for j = 0 to mm-1 do
    let yfact = data.dx *. float_of_int j
    and offset = mm*j in
    for i = 0 to mm-1 do
      let xfact = data.dx *. float_of_int i
      and loc = offset + i in
      u.{loc} <- 16.0 *. xfact *. (1. -. xfact) *. yfact *. (1. -. yfact)
    done
  done;

  (* Initialize u' vector to 0. *)
  RealArray.fill u' 0.;

  (* heatres sets res to negative of ODE RHS values at interior points. *)
  heatres 0. u u' res data;

  (* Copy -res into u' to get correct interior initial u' values. *)
  vscale (-1.) res u';

  (* Finally, set values of u, u', and id at boundary points. *)
  for j = 0 to mm-1 do
    let offset = mm*j in
    for i = 0 to mm-1 do
      let loc = offset + i in
      if j = 0 || j = mm1 || i = 0 || i = mm1
      then (u.{loc} <- bval;
            u'.{loc} <- 0.;
            id.{loc} <- Ida.VarId.algebraic)
    done
  done

let idaklu =
  match Config.sundials_version with 2,_,_ -> "IDAKLU" | _ -> "KLU"

let print_header rtol atol =
  printf "\nidaHeat2D_klu: Heat equation, serial example problem for IDA\n";
  printf "          Discretized heat equation on 2D unit square.\n";
  printf "          Zero boundary conditions,";
  printf " polynomial initial conditions.\n";
  printf "          Mesh dimensions: %d x %d" mgrid mgrid;
  printf "        Total system size: %d\n\n" neq;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Constraints set to force all solution components >= 0. \n";
  printf "Linear solver: %s, sparse direct solver \n" idaklu;
  printf "       difference quotient Jacobian\n";
  printf "IDACalcIC called with input boundary values = %g \n" bval;
  (* Print output table heading and initial line of table.  *)
  printf "\n   Output Summary (umax = max-norm of solution) \n\n";
  printf "  time       umax     k  nst  nni  nje   nre     h       \n";
  printf " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n"

let print_output mem t u =
  let umax = vmax_norm u in
  let open Ida in
  let kused = get_last_order mem
  and nst   = get_num_steps mem
  and nni   = get_num_nonlin_solv_iters mem
  and nre   = get_num_res_evals mem
  and hused = get_last_step mem
  and nje   = Dls.get_num_jac_evals mem
  in
  printf " %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %9.2e \n"
         t umax kused nst nni nje nre hused

let main () =
  (* Create vectors uu, up, res, constraints, id.  *)
  let u = RealArray.create neq
  and u' = RealArray.create neq     (* du/dt *)
  and res = RealArray.create neq
  and constraints = RealArray.create neq
  and id = RealArray.create neq in  (* differentiate between algebraic and
                                       differential *)

  (* Create and load problem data block.  *)
  let data =
    let mm = mgrid in
    let dx = 1. /. (float_of_int mgrid -. 1.) in
    let coeff = 1. /. ( dx *. dx ) in
    { mm=mm; dx=dx; coeff=coeff }
  in

  (* Initialize u, u', id.  *)
  set_initial_profile data u u' id res;

  (* Set constraints to all NonNegative for nonnegative solution values.  *)
  RealArray.fill constraints Constraint.geq_zero;

  (* Set remaining input parameters.  *)
  let t0 = 0.
  and t1 = 0.01
  and rtol = 0.
  and atol = 1.0e-8
  in

  (* Wrap u and u' in nvectors.  Operations performed on the wrapped
     representation affect the originals u and u'.  *)
  let wu = Nvector_serial.wrap u
  and wu' = Nvector_serial.wrap u'
  in

  (* Call IDACreate and IDAMalloc to initialize solution, and call IDABand to
     specify the linear solver.  *)
  let nnz = neq * neq in
  let jacfn =
    if mgrid >= 4 then jac_heat
    else if mgrid = 3 then jac_heat3
    else failwith "mgrid size is too small to run."
  in
  let m = Matrix.sparse_csc ~nnz neq in
  let mem =
    Ida.(init (SStolerances (rtol, atol))
              ~lsolver:Dls.(solver ~jac:jacfn (klu wu m))
              (fun t u u' r -> heatres t u u' r data)
              t0 wu wu')
  in
  Ida.set_constraints mem (Nvector_serial.wrap constraints);

  (* Call IDASetId and IDACalcIC to correct the initial values.  *)
  Ida.calc_ic_ya_yd' mem ~varid:(Nvector_serial.wrap id) t1;

  (* Print output heading. *)
  print_header rtol atol;

  print_output mem t0 u;

  (* Loop over output times, call IDASolve, and print results. *)
  let tout = ref t1 in
  for iout = 1 to nout do
    let (tret, flag) = Ida.solve_normal mem !tout wu wu' in
    print_output mem tret u;
    tout := 2. *. !tout
  done;

  (* Print remaining counters. *)
  let netf = Ida.get_num_err_test_fails mem
  and ncfn = Ida.get_num_nonlin_solv_conv_fails mem
  in
  printf "\n netf = %d,   ncfn = %d \n" netf ncfn


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
