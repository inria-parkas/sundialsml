(* 1. Define right-hand-side functions. *)
let f_e t y yd = yd.{0} <- y.{1}
let f_i t y yd = yd.{1} <- -9.81

(* 2. Optionally define a root function. *)
let g t y gout = gout.{0} <- 1.0 -. y.{0}

(* 3. Set vector of initial values.
      The length of this vector determines the problem size. *)
let yd = Sundials.RealArray.of_list [ 10.0; 0.0 ]
let y = Nvector_serial.wrap yd

(* 4. Create and initialize a solver session.
      This will determine whether the problem is purely explicit, purely
      implicit, or both. Initialize a specific linear solver if there is an
      implicit component, and root-finding if necessary. It is also possible to
      specify a mass matrix solver. *)
let m = Matrix.dense 2
let s = Arkode.(
  init
    (ImEx { explicit = f_e;
            implicit = (f_i, Newton Arkode.Dls.(solver (dense y m) m),
                        Linear true); })
    (SStolerances (1e-4, 1e-9))
    ~roots:(1, g)
    0.0
    y);;

(* 5. Set optional inputs, e.g.,
      call [set_*] functions to change solver parameters. *)
Arkode.set_stop_time s 10.0;;
Arkode.set_all_root_directions s Sundials.RootDirs.Increasing;;

(* 6. Advance the solution in time,
      by repeatedly calling [solve_normal] or [solve_one_step]. *)
let rec go (t, r) =
  Printf.printf "% .10e\t% .10e\t% .10e\n" t yd.{0} yd.{1};
  match r with
  | Arkode.Success -> go (Arkode.solve_normal s (t +. 0.5) y)
  | Arkode.RootsFound -> begin
        yd.{1} <- -0.8 *. yd.{1};
        Arkode.reinit s t y;
        go (t, Arkode.Success)
      end
  | Arkode.StopTimeReached -> ();;

Printf.printf "time\ty\ty'\n";;
go (0.0, Arkode.Success);;

(* 7. Get optional outputs,
      call the [get_*] functions to examine solver statistics. *)
let ns = Arkode.get_num_steps s
