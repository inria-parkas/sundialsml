open Sundials
module ARKStep = Arkode.ARKStep

(* 1. Define right-hand-side functions. *)
let fe t y yd = yd.{0} <- y.{1}
let fi t y yd = yd.{1} <- -9.81

(* 2. Optionally define a root function. *)
let g t y gout = gout.{0} <- 1.0 -. y.{0}

(* 3. Set vector of initial values.
      The length of this vector determines the problem size. *)
let yd = RealArray.of_list [ 10.0; 0.0 ]
let y = Nvector_serial.wrap yd

(* 4. Create and initialize a solver session.
      This will determine whether the problem is purely explicit, purely
      implicit, or both. Initialize a specific linear solver if there is an
      implicit component, and root-finding if necessary. It is also possible to
      specify a mass matrix solver. *)
let m = Matrix.dense 2
let s = ARKStep.(
  init
    (imex ~lsolver:Dls.(solver (dense y m))
          ~linearity:(Linear true)
          ~fi
          fe)
    (SStolerances (1e-4, 1e-9))
    ~roots:(1, g)
    0.0
    y);;

(* 5. Set optional inputs, e.g.,
      call [set_*] functions to change solver parameters. *)
ARKStep.set_stop_time s 10.0;;
ARKStep.set_all_root_directions s RootDirs.Increasing;;

(* 6. Advance the solution in time,
      by repeatedly calling [solve_normal] or [solve_one_step]. *)
let rec go (t, r) =
  Printf.printf "% .10e\t% .10e\t% .10e\n" t yd.{0} yd.{1};
  match r with
  | ARKStep.Success -> go (ARKStep.solve_normal s (t +. 0.5) y)
  | ARKStep.RootsFound -> begin
        yd.{1} <- -0.8 *. yd.{1};
        ARKStep.reinit s t y;
        go (t, ARKStep.Success)
      end
  | ARKStep.StopTimeReached -> ();;

Printf.printf "time\ty\ty'\n";;
go (0.0, ARKStep.Success);;

(* 7. Get optional outputs,
      call the [get_*] functions to examine solver statistics. *)
let ns = ARKStep.get_num_steps s
