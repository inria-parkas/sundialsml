open Sundials
module MRIStep = Arkode.MRIStep
module ARKStep = Arkode.ARKStep

(* 1. Define right-hand-side functions. *)
let fslow t y yd = yd.{0} <- y.{1};
                   yd.{1} <- 0.0
let ffast t y yd = yd.{0} <- 0.0;
                   yd.{1} <- -9.81

(* 2. Optionally define a root function. *)
let g t y gout = gout.{0} <- 1.0 -. y.{0}

(* 3. Set vector of initial values.
      The length of this vector determines the problem size. *)
let y = RealArray.of_list [ 10.0; 0.0 ]
let ynv = Nvector_serial.wrap y

(* 4. Create an ARKStep solver for the fast (inner) integration *)
let sf = ARKStep.(init (explicit ffast) default_tolerances 0.0 ynv)

(* 5. Create and initialize a solver session. *)
let slowstep = 0.01
let s = MRIStep.(init sf default_tolerances fslow ~roots:(1, g) ~slowstep 0.0 ynv);;

(* 6. Set optional inputs, e.g.,
      call [set_*] functions to change solver parameters. *)
MRIStep.set_stop_time s 10.0;;
MRIStep.set_all_root_directions s RootDirs.Increasing;;

(* 7. Advance the solution in time,
      by repeatedly calling [solve_normal] or [solve_one_step]. *)
let rec go (t, r) =
  Printf.printf "% .10e\t% .10e\t% .10e\n" t y.{0} y.{1};
  match r with
  | MRIStep.Success -> go (MRIStep.solve_normal s (t +. 0.5) ynv)
  | MRIStep.RootsFound -> begin
        y.{1} <- -0.8 *. y.{1};
        MRIStep.reinit s t ynv;
        ARKStep.reinit sf t ynv;
        go (t, MRIStep.Success)
      end
  | MRIStep.StopTimeReached -> ();;

Printf.printf "time\ty\ty'\n";;
go (0.0, MRIStep.Success);;

(* 8. Get optional outputs,
      call the [get_*] functions to examine solver statistics. *)
let ns = MRIStep.get_num_steps s
