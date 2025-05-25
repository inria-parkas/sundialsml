open Sundials
module SPRKStep = Arkode.SPRKStep

(* 1. Define right-hand-side functions. *)
let f1 _t y yd  = yd.{0} <- y.{1}
let f2 _t y yd = yd.{1} <- -2.0 *. y.{0}

(* 2. Optionally define a root function. *)
let g _t y gout = gout.{0} <- 1.0 -. y.{0}

(* 3. Set vector of initial values.
      The length of this vector determines the problem size. *)
let y = RealArray.of_list [ 10.0 *. cos(0.0); -20.0 *. sin(0.0) ]
let ynv = Nvector_serial.wrap y

(* 4. Create and initialize a solver session. *)
let slowstep = 0.01
let s = SPRKStep.(init ~step:0.01 ~order:4 ~f1 ~f2 ~roots:(1, g) 0.0 ynv);;

(* 5. Set optional inputs, e.g.,
      call [set_*] functions to change solver parameters. *)
SPRKStep.set_stop_time s 100.0;;
SPRKStep.set_use_compensated_sums s true;;
SPRKStep.set_all_root_directions s RootDirs.Increasing;;

(* 6. Advance the solution in time,
      by repeatedly calling [evolve_normal] or [evolve_one_step]. *)
let rec go (t, r) =
  Printf.printf "% .10e\t% .10e\t% .10e\n" t y.{0} y.{1};
  match r with
  | SPRKStep.Success -> go (SPRKStep.evolve_normal s (t +. 0.5) ynv)
  | SPRKStep.RootsFound -> begin
        y.{1} <- -0.8 *. y.{1};
        SPRKStep.reinit s t ynv;
        go (t, SPRKStep.Success)
      end
  | SPRKStep.StopTimeReached -> ();;

Printf.printf "time\ty\ty'\n";;
go (0.0, SPRKStep.Success);;

(* 7. Get optional outputs,
      call the [get_*] functions to examine solver statistics. *)
let ns = SPRKStep.get_num_steps s
