(* 1. Define a right-hand-side function. *)
let f t y yd = yd.{0} <- y.{1}; yd.{1} <- -9.81

(* 2. Optionally define a root function. *)
let g t y gout = gout.{0} <- 1.0 -. y.{0}

(* 3. Set vector of initial values.
      The length of this vector determines the problem size. *)
let yd = Sundials.RealArray.of_list [ 10.0; 0.0 ]
let y = Nvector_serial.wrap yd

(* 4. Create and initialize a solver session.
      This will initialize a specific linear solver and the root-finding
      mechanism, if necessary. *)
let s = Cvode.init Cvode.Adams Cvode.Functional
                   (Cvode.SStolerances (1e-4, 1e-8))
                   f ~roots:(1, g) 0.0 y;;

(* 5. Set optional inputs, e.g.,
      call [set_*] functions to change solver parameters. *)
Cvode.set_stop_time s 10.0;;
Cvode.set_all_root_directions s Sundials.RootDirs.Increasing;;

(* 6. Advance the solution in time,
      by repeatedly calling [solve_normal] or [solve_one_step]. *)
let rec go (t, r) =
  Printf.printf "% .10e\t% .10e\t% .10e\n" t yd.{0} yd.{1};
  match r with
  | Cvode.Success -> go (Cvode.solve_normal s (t +. 0.5) y)
  | Cvode.RootsFound -> begin
        yd.{1} <- -0.8 *. yd.{1};
        Cvode.reinit s t y;
        go (t, Cvode.Success)
      end
  | Cvode.StopTimeReached -> ();;

Printf.printf "time\ty\ty'\n";;
go (0.0, Cvode.Success);;

(* 7. Get optional outputs,
      call the [get_*] functions to examine solver statistics. *)
let ns = Cvode.get_num_steps s
