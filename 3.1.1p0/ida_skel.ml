open Sundials

(* 1. Define a residual function. *)
let resf t v v' r =
  r.{0} <- v.{2}  -. v'.{0};
  r.{1} <- v.{3}  -. v'.{1};
  r.{2} <- v'.{2} -. v.{4} *. v.{0};
  r.{3} <- v'.{3} -. v.{4} *. v.{1} +. 9.81;
  r.{4} <- v.{2}*.v.{2} +. v'.{2}*.v.{0} +. v.{3}*.v'.{1} +. v'.{3}*.v.{1}

(* 2. Optionally define a root function. *)
let g t v v' gout = gout.{0} <- v.{0} -. v.{1} *. 0.5774

(* 3. Set vector of initial values.
      The length of this vector determines the problem size. *)
let vd = RealArray.of_list [ 0.9848; 0.1736; 0.0; 0.0; 0.0 ]
let v  = Nvector_serial.wrap vd
let v' = Nvector_serial.make 5 0.0

(* 4. Create and initialize a solver session.
      This will initialize a specific linear solver and the root-finding
      mechanism, if necessary. *)
let m = Matrix.dense 5
let s = Ida.(init Dls.(solver (dense v m))
                  (SStolerances (1e-9, 1e-9))
                  resf ~roots:(1, g) 0.0 v v');;

(* 5. Set optional inputs, e.g.,
      call [set_*] functions to change solver parameters. *)
Ida.set_stop_time s 10.0;
Ida.set_all_root_directions s RootDirs.Decreasing;;

(* 6. Correct initial values *)
let vids = Nvector_serial.make 5 Ida.VarId.differential;;
(Nvector_serial.unwrap vids).{4} <- Ida.VarId.algebraic;

Ida.set_suppress_alg s ~varid:vids true;
Ida.calc_ic_ya_yd' s ~y:v ~y':v' 0.1;;

(* 7. Advance the solution in time,
      by repeatedly calling [solve_normal] or [solve_one_step]. *)
let rec go (t, r) =
  Printf.printf "% .10e\t% .10e\t% .10e\t% .10e\t% .10e\t% .10e\n"
                t vd.{0} vd.{1} vd.{2} vd.{3} vd.{4};
  match r with
  | Ida.Success -> go (Ida.solve_normal s (t +. 0.05) v v')
  | Ida.RootsFound -> begin
         vd.{2} <- -. 0.5 *. vd.{2};
         vd.{3} <- -. 0.5 *. vd.{3};
         Ida.reinit s t v v';
         Ida.calc_ic_ya_yd' ~y:v ~y':v' s ~varid:vids (t +. 0.05);
          go (t, Ida.Success)
      end
  | Ida.StopTimeReached -> ();;

Printf.printf "time\tx\ty\tx'\ty'\tp\n";;
go (0.0, Ida.Success);;

(* 8. Get optional outputs,
      call the [get_*] functions to examine solver statistics. *)
let ns = Ida.get_num_steps s
