(*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/14 22:15:31 $
 * -----------------------------------------------------------------
 * Programmer(s): Lukas Jager and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * OCaml port: Timothy Bourke, Inria, Aug 2014.
 * -----------------------------------------------------------------
 * Parallel Krylov adjoint sensitivity example problem.
 * -----------------------------------------------------------------
 *)

module RealArray = Sundials.RealArray
module Quad = Cvodes.Quadrature
module Adj = Cvodes.Adjoint
module QuadAdj = Cvodes.Adjoint.Quadrature
module Nvector = Nvector_parallel
module Bbd = Cvode_bbd
module Adjbbd = Cvodes_bbd
open Bigarray

let printf = Printf.printf
let eprintf = Printf.eprintf
let fprintf = Printf.fprintf
let sqr x = x ** 2.0
let unwrap = Nvector_parallel.unwrap

let n_vdotprod = Nvector.DataOps.n_vdotprod
let n_vscale = Nvector.DataOps.n_vscale

(*
 *------------------------------------------------------------------
 * Constants
 *------------------------------------------------------------------
 *)

let dim = 2 (* = 3 (* USE3D *) *)

(* Domain definition *)

let xmin = 0.0
let xmax = 20.0
let mx =   20    (* no. of divisions in x dir. *)
let npx =  2     (* no. of procs. in x dir.    *)

let ymin = 0.0
let ymax = 20.0
let my =   40    (* no. of divisions in y dir. *)
let npy =  2     (* no. of procs. in y dir.    *)

let zmin = 0.0
let zmax = 20.0
let mz =   20    (* no. of divisions in z dir. *)
let npz =  1     (* no. of procs. in z dir.    *)

(* Parameters for source Gaussians *)

let g1_ampl =   1.0
let g1_sigma =  1.7 
let g1_x =      4.0
let g1_y =      8.0
let g1_z =      8.0

let g2_ampl =   0.8
let g2_sigma =  3.0
let g2_x =      16.0
let g2_y =      12.0
let g2_z =      12.0

let g_min =     1.0e-5

(* Diffusion coeff., max. velocity, domain width in y dir. *)

let diff_coef = 1.0
let v_max =     1.0
let l =         (ymax-.ymin)/.2.0
let v_coeff =   v_max/.(ymax-.ymin)/.2.0/.(ymax-.ymin)/.2.0

(* Initial and final times *)

let ti =    0.0
let tf =    10.0

(* Integration tolerances *)

let rtol =    1.0e-8 (* states *)
let atol =    1.0e-6

let rtol_q =  1.0e-8 (* forward quadrature *)
let atol_q =  1.0e-6

let rtol_b =  1.0e-8 (* adjoint variables *)
let atol_b =  1.0e-6

let rtol_qb = 1.0e-8 (* backward quadratures *)
let atol_qb = 1.0e-6

(* Steps between check points *)

let steps = 200

let zero = 0.0
let one =  1.0
let two =  2.0

(*
 *------------------------------------------------------------------
 * Macros
 *------------------------------------------------------------------
 *)

let rec enum i n = if i = n then [] else i::enum (i + 1) n 
let dims = enum 0 dim

let takedim xs =
  let rec take n_xs =
    match n_xs with
      (0, _) -> []
    | (n, x::xs) -> x::take (n - 1, xs)
    | _ -> assert false
  in take (dim, xs)

let g1 = Array.of_list (takedim [ g1_x; g1_y; g1_z ])
let g2 = Array.of_list (takedim [ g2_x; g2_y; g2_z ])

(* IJth:     (i[0],i[1],i[2])-th vector component                       *)
(* IJth_ext: (i[0],i[1],i[2])-th vector component in the extended array *)

let ijth =
  if dim = 3 then
    fun l_m y i -> y.{i.(0)+(l_m.(0)*(i.(1)+i.(2)*l_m.(1)))}
  else
    fun l_m y i -> y.{i.(0)+i.(1)*l_m.(0)}

let set_ijth =
  if dim = 3 then
    fun l_m y i v -> y.{i.(0)+(l_m.(0)*(i.(1)+i.(2)*l_m.(1)))} <- v
  else
    fun l_m y i v -> y.{i.(0)+i.(1)*l_m.(0)} <- v

let add_to_ijth =
  if dim = 3 then
    fun l_m y i v ->
      let idx = i.(0)+(l_m.(0)*(i.(1)+i.(2)*l_m.(1))) in
      y.{idx} <- y.{idx} +. v
  else
    fun l_m y i v ->
      let idx = i.(0)+i.(1)*l_m.(0) in
      y.{idx} <- y.{idx} +. v

let ijth_ext =
  if dim = 3 then
    fun l_m y i -> y.{(i.(0)+1)+((l_m.(0)+2)*((i.(1)+1)+(i.(2)+1)*(l_m.(1)+2)))}
  else
    fun l_m y i -> y.{(i.(0)+1) + (i.(1)+1) * (l_m.(0)+2)}

let set_ijth_ext =
  if dim = 3 then
    fun l_m y i v ->
      y.{(i.(0)+1)+((l_m.(0)+2)*((i.(1)+1)+(i.(2)+1)*(l_m.(1)+2)))} <- v
  else
    fun l_m y i v -> y.{(i.(0)+1) + (i.(1)+1) * (l_m.(0)+2)} <- v

(*
 *------------------------------------------------------------------
 * Type definition: ProblemData 
 *------------------------------------------------------------------
 *)

type problem_data = {
    (* Domain *)
    xmin      : float array; (* "left" boundaries *)  
    xmax      : float array; (* "right" boundaries *)
    m         : int array;   (* number of grid points *)
    dx        : float array; (* grid spacing *)
    dOmega    : float;       (* differential volume *)

    (* Parallel stuff *)
    comm : Mpi.communicator; (* MPI communicator *)
    myId      : int;         (* process id *) 
    npes      : int;         (* total number of processes *)
    num_procs : int array;   (* number of processes in each direction *)
    nbr_left  : int array;   (* MPI ID of "left" neighbor *)
    nbr_right : int array;   (* MPI ID of "right" neighbor *)
    m_start   : int array;   (* "left" index in the global domain *)
    l_m       : int array;   (* number of local grid points *) 

    y_ext            : RealArray.t; (* extended data array *)
    mutable buf_send : RealArray.t; (* Send buffer *)
    mutable buf_size : int;         (* Buffer size *)

    (* Source *)
    p         : Nvector.t;   (* Source parameters *) 
  }

(*
 *------------------------------------------------------------------
 * Private functions
 *------------------------------------------------------------------
 *)

let set_source d =
  let l_m  = d.l_m in
  let m_start = d.m_start in
  let xmin = d.xmin in
  let dx = d.dx in
  let pdata = unwrap d.p in

  let ii = Array.make dim 0 in

  let rec f p1 p2 dims =
    match dims with
      [] -> let g = p1 +. p2 in
            set_ijth l_m pdata ii (if g < g_min then zero else g)
    | d::ds ->
        for i = 0 to l_m.(d) - 1 do
          ii.(d) <- i;
          let x = xmin.(d) +. float (m_start.(d) + i) *. dx.(d) in
          let p1 = p1 *. exp (-. sqr (g1.(d) -. x) /. sqr g1_sigma) in
          let p2 = p2 *. exp (-. sqr (g2.(d) -. x) /. sqr g2_sigma) in
          f p1 p2 ds
        done
  in
  f g1_ampl g2_ampl dims

(*
 *------------------------------------------------------------------
 * SetData:
 * Allocate space for the ProblemData structure.
 * Set fields in the ProblemData structure.
 * Return local and global problem dimensions.
 *
 * SetSource:
 * Instantiates the source parameters for a combination of two
 * Gaussian sources.
 *------------------------------------------------------------------
 *)

let array_of_third (x, y, z) = Array.of_list z
let array_of_second (x, y) = Array.of_list y

let set_data comm npes myId =
  (* Set domain boundaries *)
  let xmin = Array.of_list (takedim [ xmin; ymin; zmin ]) in
  let xmax = Array.of_list (takedim [ xmax; ymax; zmax ]) in
  let m    = Array.of_list (takedim [ mx + 1; my + 1; mz + 1 ]) in

  (* Calculate grid spacing and differential volume *)
  let dx = Array.init dim (fun i -> (xmax.(i) -. xmin.(i)) /. float (m.(i) - 1)) in
  let dOmega = Array.fold_left (fun a dx -> a *. dx) one dx in
  
  (* Set partitioning *)
  let num_procs = Array.of_list (takedim [ npx; npy; npz ]) in
  let n         = Array.copy num_procs in
  let nd        = Array.init dim (fun i -> m.(i) / n.(i)) in

  (* Compute the neighbors *)
  let nbr_left  = array_of_third (
    List.fold_left (function (dv, off, r) -> fun i ->
        dv / n.(i), off * n.(i),
        r @ [if (dv mod n.(i)) = 0 then myId else myId - off])
      (myId, 1, []) dims)
  in
  let nbr_right = array_of_third (
    List.fold_left (function (dv, off, r) -> fun i ->
        dv / n.(i), off * n.(i),
        r @ [if (dv mod n.(i)) = n.(i)-1 then myId else myId + off])
      (myId, 1, []) dims)
  in
 
  (* Compute the local subdomains 
     m_start: left border in global index space 
     l_m:     length of the subdomain *)
  let m_start = array_of_second (
    List.fold_left (function (dv, r) -> fun i ->
        dv / n.(i), r @ [ (dv mod n.(i)) * nd.(i) ])
      (myId, []) dims)
  in
  let l_m = Array.of_list (List.map (fun i -> 
        if nbr_right.(i) = myId then m.(i) - m_start.(i) else nd.(i)) dims)
  in

  (* Allocate memory for the y_ext array 
     (local solution + data from neighbors) *)
  let yext_size = List.fold_left (fun s i -> s * (l_m.(i) + 2)) 1 dims in

  (* Allocate space for the source parameters *)
  let neq = List.fold_left (fun s i -> s * m.(i)) 1 dims in
  let l_neq = List.fold_left (fun s i -> s * l_m.(i)) 1 dims in

  let data = {
    xmin   = xmin;
    xmax   = xmax;
    m      = m;
    dx     = dx;
    dOmega = dOmega;

    (* Set MPI communicator, id, and total number of processes *)

    comm = comm;
    myId = myId;
    npes = npes;

    num_procs = num_procs;
    nbr_left  = nbr_left;
    nbr_right = nbr_right;
    m_start   = m_start;
    l_m       = l_m;

    y_ext = RealArray.create yext_size;

    (* Initialize Buffer field.
       Size of buffer is checked when needed *)
    buf_send = RealArray.create 0;
    buf_size = 0;

    p = Nvector.make l_neq neq comm zero;
  } in

  (* Initialize the parameters for a source with Gaussian profile *)
  set_source data;
  data, neq, l_neq

(*
 *------------------------------------------------------------------
 * f_comm: 
 * Function for inter-process communication
 * Used both for the forward and backward phase.
 *------------------------------------------------------------------
 *)

let f_comm d t (ydata, _, _) =
  let comm = d.comm in
  let id = d.myId in
  
  (* extract data from domain*)
  let n = Array.copy d.num_procs in
  let l_m = Array.copy d.l_m in
  let yextdata = d.y_ext in
  
  (* Calculate required buffer size *)
  let size = (Array.fold_left (fun s l -> s * l) 1 l_m)
             / (Array.fold_left min max_int l_m)
  in

  (* Adjust buffer size if necessary *)
  if d.buf_size < size then begin
    d.buf_send <- RealArray.create size;
    d.buf_size <- size
  end;
  let buf_send = d.buf_send in
  
  (* Compute the communication pattern; who sends first? *)
  (* if proc_cond==1 , process sends first in this dimension *)
  let proc_cond = Array.make dim 0 in
  proc_cond.(0) <- (id mod n.(0)) mod 2;
  proc_cond.(1) <- ((id/n.(0)) mod n.(1)) mod 2;
  if dim = 3 then
    proc_cond.(2) <- (id/n.(0)/n.(1)) mod 2;

  (* Compute the actual communication pattern *)
  (* nbr[dim][0] is first proc to communicate with in dimension dim *)
  (* nbr[dim][1] the second one *)
  let nbr = Array.init dim (fun i ->
    Array.init 2 (fun j ->
      if j = proc_cond.(i) then d.nbr_left.(i) else d.nbr_right.(i)))
  in

  (* Communication: loop over dimension and direction (left/right) *)
  let i = Array.make dim 0 in
  for d = 0 to dim - 1 do
    for dir = 0 to 1 do
      (* If subdomain at boundary, no communication in this direction *)
      if id <> nbr.(d).(dir) then begin
        (* Compute the index of the boundary (right or left) *)
        i.(d) <- if (dir lxor proc_cond.(d) <> 0) then l_m.(d)-1 else 0;

        (* Loop over all other dimensions and copy data into buf_send *)
        let l = Array.init (dim - 1) (fun i -> (d + 1 + i) mod dim) in

        let c = ref 0 in
        let rec fill_send_buf () =
          for k=0 to l_m.(l.(0)) - 1 do
            i.(l.(0)) <- k;
            buf_send.{!c} <- ijth l_m ydata i;
            incr c
          done
        in
        if dim = 3 then
          for j=0 to l_m.(l.(1)) - 1 do
            i.(l.(1)) <- j;
            fill_send_buf ()
          done
        else fill_send_buf ();
          
        let buf_recv =
          if proc_cond.(d) = 1 then begin
            (* Send buf_send and receive into buf_recv *)
            Mpi.send buf_send nbr.(d).(dir) 0 comm;
            (Mpi.receive nbr.(d).(dir) 0 comm : RealArray.t)
          end else begin
            (* Receive into buf_recv and send buf_send*)
            let r = (Mpi.receive nbr.(d).(dir) 0 comm : RealArray.t) in
            Mpi.send buf_send nbr.(d).(dir) 0 comm;
            r
          end
        in

        (* Compute the index of the boundary (right or left) in yextdata *)
        i.(d) <- if dir lxor proc_cond.(d) <> 0 then l_m.(d) else -1;

        (* Loop over all other dimensions and copy data into yextdata *)
        let c = ref 0 in
        let rec empty_recv_buf () =
          for k=0 to l_m.(l.(0)) - 1 do
            i.(l.(0)) <- k;
            set_ijth_ext l_m yextdata i buf_recv.{!c};
            incr c
          done
        in
        if dim = 3 then
          for j=0 to l_m.(l.(1)) - 1 do
            i.(l.(1)) <- j;
            empty_recv_buf ()
          done
        else empty_recv_buf ()
      end
    done (* end loop over direction *)
  done (* end loop over dimension *) 

(*
 *------------------------------------------------------------------
 * Load_yext: 
 * copies data from src (y or yB) into y_ext, which already contains
 * data from neighboring processes.
 *------------------------------------------------------------------
 *)

let load_yext data src =
  let i = Array.make dim 0 in
  let l_m = data.l_m in
  (* copy local segment *)
  let rec copy d =
    if d < 0 then
      set_ijth_ext l_m data.y_ext i (ijth l_m src i)
    else
      for j = 0 to l_m.(d) - 1 do
        i.(d) <- j;
        copy (d - 1)
      done
  in copy (dim - 1)

(*
 *------------------------------------------------------------------
 * PrintHeader:
 * Print first lins of output (problem description)
 *------------------------------------------------------------------
 *)

let print_header () =
    printf "\nParallel Krylov adjoint sensitivity analysis example\n";
    printf "%1dD Advection diffusion PDE with homogeneous Neumann B.C.\n" dim;
    printf "Computes gradient of G = int_t_Omega ( c_i^2 ) dt dOmega\n";
    printf "with respect to the source values at each grid point.\n\n";

    printf "Domain:\n";

    printf "   %f < x < %f   mx = %d  npe_x = %d \n" xmin xmax mx npx;
    printf "   %f < y < %f   my = %d  npe_y = %d \n" ymin ymax my npy;
    if dim = 3 then
      printf "   %f < z < %f   mz = %d  npe_z = %d \n" zmin zmax mz npz;
    printf "\n"

(*
 *------------------------------------------------------------------
 * PrintFinalStats:
 * Print final statistics contained in cvode_mem
 *------------------------------------------------------------------
 *)

let print_final_stats s =
  let lenrw, leniw = Cvode.get_work_space s
  and nst          = Cvode.get_num_steps s
  and nfe          = Cvode.get_num_rhs_evals s
  and nsetups      = Cvode.get_num_lin_solv_setups s
  and netf         = Cvode.get_num_err_test_fails s
  and nni          = Cvode.get_num_nonlin_solv_iters s
  and ncfn         = Cvode.get_num_nonlin_solv_conv_fails s
  in
  let lenrwSPGMR, leniwSPGMR = Cvode.Spils.get_work_space s
  and nli      = Cvode.Spils.get_num_lin_iters s
  and npe      = Cvode.Spils.get_num_prec_evals s
  and nps      = Cvode.Spils.get_num_prec_solves s
  and ncfl     = Cvode.Spils.get_num_conv_fails s
  and nfeSPGMR = Cvode.Spils.get_num_rhs_evals s
  in
  printf "\nFinal Statistics.. \n\n";
  printf "lenrw   = %6d     leniw = %6d\n"   lenrw leniw;
  printf "llrw    = %6d     lliw  = %6d\n"   lenrwSPGMR leniwSPGMR;
  printf "nst     = %6d\n"                   nst;
  printf "nfe     = %6d     nfel  = %6d\n"   nfe nfeSPGMR;
  printf "nni     = %6d     nli   = %6d\n"   nni nli;
  printf "nsetups = %6d     netf  = %6d\n"   nsetups netf;
  printf "npe     = %6d     nps   = %6d\n"   npe nps;
  printf "ncfn    = %6d     ncfl  = %6d\n\n" ncfn ncfl

let print_final_statsB s =
  let lenrw, leniw = Adj.get_work_space s
  and nst          = Adj.get_num_steps s
  and nfe          = Adj.get_num_rhs_evals s
  and nsetups      = Adj.get_num_lin_solv_setups s
  and netf         = Adj.get_num_err_test_fails s
  and nni          = Adj.get_num_nonlin_solv_iters s
  and ncfn         = Adj.get_num_nonlin_solv_conv_fails s
  in
  let lenrwSPGMR, leniwSPGMR = Adj.Spils.get_work_space s
  and nli      = Adj.Spils.get_num_lin_iters s
  and npe      = Adj.Spils.get_num_prec_evals s
  and nps      = Adj.Spils.get_num_prec_solves s
  and ncfl     = Adj.Spils.get_num_conv_fails s
  and nfeSPGMR = Adj.Spils.get_num_rhs_evals s
  in
  printf "\nFinal Statistics.. \n\n";
  printf "lenrw   = %6d     leniw = %6d\n"   lenrw leniw;
  printf "llrw    = %6d     lliw  = %6d\n"   lenrwSPGMR leniwSPGMR;
  printf "nst     = %6d\n"                   nst;
  printf "nfe     = %6d     nfel  = %6d\n"   nfe nfeSPGMR;
  printf "nni     = %6d     nli   = %6d\n"   nni nli;
  printf "nsetups = %6d     netf  = %6d\n"   nsetups netf;
  printf "npe     = %6d     nps   = %6d\n"   npe nps;
  printf "ncfn    = %6d     ncfl  = %6d\n\n" ncfn ncfl

(*
 *------------------------------------------------------------------
 * OutputGradient:
 * Generate matlab m files for visualization
 * One file gradXXXX.m from each process + a driver grad.m
 *------------------------------------------------------------------
 *)

let output_gradient data myId qB =
  let filename = Printf.sprintf "grad%03d.m" myId in
  let fid = open_out filename in

  let l_m     = data.l_m in
  let m_start = data.m_start in
  let xmin    = data.xmin in
  let dx      = data.dx in

  let qBdata = unwrap qB in
  let pdata = unwrap data.p in

  (* Write matlab files with solutions from each process *)
  let i = Array.make dim 0 in
  let x = Array.make dim zero in
  for j=0 to l_m.(0) - 1 do
    i.(0) <- j;
    x.(0) <- xmin.(0) +. float (m_start.(0)+i.(0)) *. dx.(0);
    for k=0 to l_m.(1) - 1 do
      i.(1) <- k;
      x.(1) <- xmin.(1) +. float (m_start.(1)+i.(1)) *. dx.(1);
      if dim = 3 then
        for j=0 to l_m.(2) - 1 do
          i.(2) <- j;
          x.(2) <- xmin.(2) +. float (m_start.(2)+i.(2)) *. dx.(2);
          let g = ijth l_m qBdata i in
          let p = ijth l_m pdata i in
          fprintf fid "x%d(%d,1) = %e; \n" myId (i.(0)+1) x.(0);
          fprintf fid "y%d(%d,1) = %e; \n" myId (i.(1)+1) x.(1);
          fprintf fid "z%d(%d,1) = %e; \n" myId (i.(2)+1) x.(2);
          fprintf fid "p%d(%d,%d,%d) = %e; \n" myId (i.(1)+1) (i.(0)+1) (i.(2)+1) p;
          fprintf fid "g%d(%d,%d,%d) = %e; \n" myId (i.(1)+1) (i.(0)+1) (i.(2)+1) g
        done
      else begin
        let g = ijth l_m qBdata i in
        let p = ijth l_m pdata i in
        fprintf fid "x%d(%d,1) = %e; \n" myId (i.(0)+1) x.(0);
        fprintf fid "y%d(%d,1) = %e; \n" myId (i.(1)+1) x.(1);
        fprintf fid "p%d(%d,%d) = %e; \n" myId (i.(1)+1) (i.(0)+1) p;
        fprintf fid "g%d(%d,%d) = %e; \n" myId (i.(1)+1) (i.(0)+1) g
      end
    done
  done;
  close_out fid;

  (* Write matlab driver *)
  if myId == 0 then begin
    let fid = open_out "grad.m" in

    if dim = 3 then begin
      fprintf fid "clear;\nfigure;\nhold on\n";
      fprintf fid "trans = 0.7;\n";
      fprintf fid "ecol  = 'none';\n";
      fprintf fid "xp=[%f %f];\n" g1_x g2_x;
      fprintf fid "yp=[%f %f];\n" g1_y g2_y;
      fprintf fid "zp=[%f %f];\n" g1_z g2_z;
      fprintf fid "ns = length(xp)*length(yp)*length(zp);\n";

      for ip=0 to data.npes - 1 do
        fprintf fid "\ngrad%03d;\n" ip;
        fprintf fid "[X Y Z]=meshgrid(x%d,y%d,z%d);\n" ip ip ip;
        fprintf fid "s%d=slice(X,Y,Z,g%d,xp,yp,zp);\n" ip ip;
        fprintf fid "for i = 1:ns\n";
        fprintf fid "  set(s%d(i),'FaceAlpha',trans);\n" ip;
        fprintf fid "  set(s%d(i),'EdgeColor',ecol);\n" ip;
        fprintf fid "end\n"
      done;
      
      fprintf fid "view(3)\n";
      fprintf fid "\nshading interp\naxis equal\n"
    end else begin
      fprintf fid "clear;\nfigure;\n";
      fprintf fid "trans = 0.7;\n";
      fprintf fid "ecol  = 'none';\n";

      for ip=0 to data.npes - 1 do
        fprintf fid "\ngrad%03d;\n" ip;

        fprintf fid "\nsubplot(1,2,1)\n";
        fprintf fid "s=surf(x%d,y%d,g%d);\n" ip ip ip;
        fprintf fid "set(s,'FaceAlpha',trans);\n";
        fprintf fid "set(s,'EdgeColor',ecol);\n";
        fprintf fid "hold on\n";
        fprintf fid "axis tight\n";
        fprintf fid "box on\n";
        
        fprintf fid "\nsubplot(1,2,2)\n";
        fprintf fid "s=surf(x%d,y%d,p%d);\n" ip ip ip;
        fprintf fid "set(s,'CData',g%d);\n" ip;
        fprintf fid "set(s,'FaceAlpha',trans);\n";
        fprintf fid "set(s,'EdgeColor',ecol);\n";
        fprintf fid "hold on\n";
        fprintf fid "axis tight\n";
        fprintf fid "box on\n"
      done
    end;
    close_out fid
  end

(*
 *------------------------------------------------------------------
 * Interface functions to CVODES
 *------------------------------------------------------------------
 *)

(*
 *------------------------------------------------------------------
 * f and f_local:
 * Forward phase ODE right-hand side
 *------------------------------------------------------------------
 *)

let f_local data t (ydata, _, _) (dydata, _, _) =
  (* Extract stuff from data structure *)
  let id = data.myId in
  let xmin = data.xmin in
  let l_m = data.l_m in
  let m_start = data.m_start in
  let dx = data.dx in
  let nbr_left = data.nbr_left in
  let nbr_right = data.nbr_right in

  (* Get pointers to vector data *)
  let pdata = unwrap data.p in

  (* Copy local segment of y to y_ext *)
  load_yext data ydata;
  let ydata = data.y_ext in

  (* Velocity components in x1 and x2 directions (Poiseuille profile) *)
  let v = Array.make dim zero in
  let i = Array.make dim 0 in
  let x = Array.make dim zero in

  (* Local domain is [xmin+(m_start+1)*dx, xmin+(m_start+1+l_m-1)*dx] *)
  let f () =
    for j=0 to l_m.(1)-1 do
      i.(1) <- j;
      x.(1) <- xmin.(1) +. float (m_start.(1)+i.(1))*.dx.(1);

      (* Velocity component in x0 direction (Poiseuille profile) *)
      let x1 = x.(1) -. xmin.(1) -. l in
      v.(0) <- v_coeff *. (l +. x1) *. (l -. x1);

      for k=0 to l_m.(0)-1 do
        i.(0) <- k;
        x.(0) <- xmin.(0) +. float (m_start.(0)+i.(0))*.dx.(0);

        let c = ijth_ext l_m ydata i in               

        (* Source term*)
        set_ijth l_m dydata i (ijth l_m pdata i);

        let fid dim =
          i.(dim)  <- i.(dim) + 1;
          let cr = ijth_ext l_m ydata i in
          i.(dim)  <- i.(dim) - 2;
          let cl = ijth_ext l_m ydata i in
          i.(dim)  <- i.(dim) + 1;

          (* Boundary conditions for the state variables *)
          let cr, cl =
            if (i.(dim) = l_m.(dim)-1) && (nbr_right.(dim) = id) then cl, cl
            else if i.(dim) = 0 && nbr_left.(dim) = id then cr, cr
            else cr, cl
          in
          let adv  = v.(dim) *. (cr-.cl) /. (two*.dx.(dim)) in
          let diff = diff_coef *. (cr-.two*.c+.cl) /. sqr dx.(dim) in

          add_to_ijth l_m dydata i (diff -. adv)
        in
        List.iter fid dims
      done
    done
  in
  if dim = 3 then
    for j=0 to l_m.(2) - 1 do
      i.(2) <- j;
      x.(2) <- xmin.(2) +. float (m_start.(2)+i.(2))*.dx.(2);
      f ()
    done
  else f ()

let f data t y ydot =
  (* Do all inter-processor communication *)
  f_comm data t y;

  (* Compute right-hand side locally *)
  f_local data t y ydot

(*
 *------------------------------------------------------------------
 * fQ:
 * Right-hand side of quadrature equations on forward integration.
 * The only quadrature on this phase computes the local contribution
 * to the function G.
 *------------------------------------------------------------------
 *)

let fQ data t y (dqdata, _, _) =
  dqdata.{0} <- (n_vdotprod y y) *. 0.5 *. data.dOmega

(*
 *------------------------------------------------------------------
 * fB and fB_local:
 * Backward phase ODE right-hand side (the discretized adjoint PDE)
 *------------------------------------------------------------------
 *)

let fB_local data t (ydata, _, _) (yBdata, _, _) (dyBdata, _, _) =
  (* Extract stuff from data structure *)
  let id = data.myId in
  let xmin = data.xmin in
  let l_m = data.l_m in
  let m_start = data.m_start in
  let dx = data.dx in
  let nbr_left = data.nbr_left in
  let nbr_right = data.nbr_right in

  (* Copy local segment of yB to y_ext *)
  load_yext data yBdata;
  let yBdata = data.y_ext in

  (* Velocity components in x1 and x2 directions (Poiseuille profile) *)
  let v = Array.make dim zero in
  let i = Array.make dim 0 in
  let x = Array.make dim zero in
 
  (* local domain is [xmin+(m_start)*dx, xmin+(m_start+l_m-1)*dx] *)
  let f () =
    for j=0 to l_m.(1)-1 do
      i.(1) <- j;
      x.(1) <- xmin.(1) +. float (m_start.(1)+i.(1))*.dx.(1);

      (* Velocity component in x0 direction (Poiseuille profile) *)
      let x1 = x.(1) -. xmin.(1) -. l in
      v.(0) <- v_coeff *. (l +. x1) *. (l -. x1);

      for k=0 to l_m.(0)-1 do
        i.(0) <- k;
        x.(0) <- xmin.(0) +. float (m_start.(0)+i.(0))*.dx.(0);

        let c = ijth_ext l_m yBdata i in               

        (* Source term for adjoint PDE *)
        set_ijth l_m dyBdata i (-. ijth l_m ydata i);

        let fid dim =
          i.(dim)  <- i.(dim) + 1;
          let cr = ijth_ext l_m yBdata i in
          i.(dim)  <- i.(dim) - 2;
          let cl = ijth_ext l_m yBdata i in
          i.(dim)  <- i.(dim) + 1;

          (* Boundary conditions for the adjoint variables *)
          let cr, cl =
            if (i.(dim) = l_m.(dim)-1) && (nbr_right.(dim) = id) then
              cl-.(two*.dx.(dim)*.v.(dim)/.diff_coef)*.c, cl
            else if i.(dim) = 0 && nbr_left.(dim) = id then
              cr, cr+.(two*.dx.(dim)*.v.(dim)/.diff_coef)*.c
            else cr, cl
          in
          let adv  = v.(dim) *. (cr-.cl) /. (two*.dx.(dim)) in
          let diff = diff_coef *. (cr-.two*.c+.cl) /. sqr dx.(dim) in

          add_to_ijth l_m dyBdata i (-. (diff +. adv))
        in
        List.iter fid dims
      done
    done
  in
  if dim = 3 then
    for j=0 to l_m.(2) - 1 do
      i.(2) <- j;
      x.(2) <- xmin.(2) +. float (m_start.(2)+i.(2))*.dx.(2);
      f ()
    done
  else f ()

let fB data t y yB yBdot =
  (* Do all inter-processor communication *)
  f_comm data t yB;

  (* Compute right-hand side locally *)
  fB_local data t y yB yBdot

(*
 *------------------------------------------------------------------
 * fQB:
 * Right-hand side of quadrature equations on backward integration
 * The i-th component of the gradient is nothing but int_t yB_i dt
 *------------------------------------------------------------------
 *)

let fQB dataB t y yB qBdot = n_vscale (-.dataB.dOmega) yB qBdot

(*
 *------------------------------------------------------------------
 * Main program
 *------------------------------------------------------------------
 *)

let main () =
  (* Initialize MPI and set Ids *)
  let comm   = Mpi.comm_world in
  let myId   = Mpi.comm_rank comm in

  (* Check number of processes *)
  let npes_needed = npx * npy * (if dim = 3 then npz else 1) in
  let npes = Mpi.comm_size comm in

  if npes_needed <> npes then begin
    if myId = 0 then
      eprintf "I need %d processes but I only got %d\n" npes_needed npes;
    exit 1
  end;

  (* Test if matlab output is requested *)
  let output = Array.length Sys.argv > 1 in

  (* Allocate and set problem data structure *)
  let d, neq, l_neq = set_data comm npes myId in

  if myId = 0 then print_header ();

  (*-------------------------- 
    Forward integration phase
    --------------------------*)

  (* Allocate space for y and set it with the I.C. *)
  let y = Nvector.make l_neq neq comm zero in
  
  (* Allocate and initialize qB (local contribution to cost) *)
  let q = Nvector.make 1 npes comm zero in

  (* Attach preconditioner and linear solver modules *)
  let spgmr = Cvode.Spils.spgmr
                (Bbd.prec_left
                   { Bbd.mudq = d.l_m.(0) + 1;
                     Bbd.mldq = d.l_m.(0) + 1;
                     Bbd.mukeep = 2;
                     Bbd.mlkeep = 2; }
                   { Bbd.local_fn = f_local d;
                     Bbd.comm_fn = None; })
  in
  (* Create CVODES object, attach user data, and allocate space *)
  let abstol, reltol = atol, rtol in
  let cvode_mem = Cvode.init Cvode.BDF (Cvode.Newton spgmr)
                    (Cvode.SStolerances (reltol, abstol))
                    (f d) ti y
  in
  
  (* Initialize quadrature calculations *)
  let abstolQ = atol_q in
  let reltolQ = rtol_q in

  Quad.init cvode_mem (fQ d) q;
  Quad.set_tolerances cvode_mem (Quad.SStolerances (reltolQ, abstolQ));

  (* Allocate space for the adjoint calculation *)
  Adj.init cvode_mem steps Adj.IHermite;

  (* Integrate forward in time while storing check points *)
  if  myId = 0 then printf "Begin forward integration... ";
  ignore (Adj.forward_normal cvode_mem tf y);
  if  myId = 0 then printf "done. ";

  (* Extract quadratures *)
  ignore (Quad.get cvode_mem q);
  let qdata = unwrap q in
  let g = Mpi.allreduce_float qdata.{0} Mpi.Float_sum comm in
  if myId = 0 then printf "  G = %e\n" g;

  (* Print statistics for forward run *)
  if myId = 0 then print_final_stats cvode_mem;

  (*-------------------------- 
    Backward integration phase
    --------------------------*)
 
  (* Allocate and initialize yB *)
  let yB = Nvector.make l_neq neq comm zero in

  (* Allocate and initialize qB (gradient) *)
  let qB = Nvector.make l_neq neq comm zero in

  (* Attach preconditioner and linear solver modules *)
  let bspgmr = Adj.Spils.spgmr
                (Adjbbd.prec_left
                   { Adjbbd.mudq = d.l_m.(0) + 1;
                     Adjbbd.mldq = d.l_m.(0) + 1;
                     Adjbbd.mukeep = 2;
                     Adjbbd.mlkeep = 2; }
                   { Adjbbd.local_fn = fB_local d;
                     Adjbbd.comm_fn = None; })
  in

  (* Create and allocate backward CVODE memory *)
  let abstolB = atol_b in
  let reltolB = rtol_b in
  let cvode_memB =
    Adj.init_backward cvode_mem
      Cvode.BDF
      (Adj.Newton bspgmr)
      (Adj.SStolerances (reltolB, abstolB))
      (Adj.NoSens (fB d))
      tf yB
  in

  (* Initialize quadrature calculations *)
  let abstolQB = atol_qb in
  let reltolQB = rtol_qb in
  QuadAdj.init cvode_memB (QuadAdj.NoSens (fQB d)) qB;
  QuadAdj.set_tolerances cvode_memB (QuadAdj.SStolerances (reltolQB, abstolQB));

  (* Integrate backwards *)
  if myId = 0 then printf "Begin backward integration... ";
  Adj.backward_normal cvode_mem ti;
  if myId = 0 then printf "done.\n";
  
  (* Extract solution *)
  ignore (Adj.get cvode_memB yB);

  (* Extract quadratures *)
  ignore (QuadAdj.get cvode_memB qB);

  (* Print statistics for backward run *)
  if myId = 0 then print_final_statsB cvode_memB;

  (* Process 0 collects the gradient components and prints them *)
  if output then begin
    output_gradient d myId qB;
    if myId = 0 then printf "Wrote matlab file 'grad.m'.\n"
  end

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
