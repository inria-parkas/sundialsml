
let ypos_i = 0
let yvel_i = 1
let yacc_i = 2
let xpos_i = 3

let under_i = 0

let gravity = -9.81 (* m/s/s *)
let t_delta = ref 0.005  (* s *)
let x_vel   = 0.8   (* m/s *)
let x_limit = 14.0  (* m *)

let f t y yd =
  yd.{xpos_i} <- x_vel;
  yd.{ypos_i} <- y.{yvel_i};
  yd.{yvel_i} <- y.{yacc_i};
  yd.{yacc_i} <- 0.0

let ground =
  [|
    4.0; (* 0 <= x < 2 *)
    3.0; (* 2 <= x < 4 *)
    2.0; (* 4 <= x < 6 *)
    1.0; (* 6 <= x < 8 *)
    0.0; (* 8 <= x *)
  |]
let ground_limits = [| 2.0; 5.0; 7.0; 11.0 |]
let ground_maxidx = Array.length ground - 1

let lookup_limit x =
  let rec f idx =
    if idx = ground_maxidx || x < ground_limits.(idx) then idx
    else f (idx + 1)
  in f 0

let g t y gout =
  let idx = lookup_limit y.{xpos_i} in
  gout.{under_i} <- y.{ypos_i} -. ground.(idx)

let y = Sundials.RealArray.make 4
let y_nvec = Nvector_serial.wrap y
let _ = y.{xpos_i} <- 0.0;
        y.{ypos_i} <- 10.0;
        y.{yvel_i} <- 0.0;
        y.{yacc_i} <- gravity

let n_roots = 1
let rootdata = Sundials.Roots.create n_roots
let ball_event s t y =
  Cvode.get_root_info s rootdata;

  if (Sundials.Roots.detected rootdata under_i && y.{yvel_i} <= 0.0) then
    (print_endline "hit ground!";
     y.{yvel_i} <- (-0.8 *. y.{yvel_i});
     Cvode.reinit s t y_nvec)

let s = Cvode.init Cvode.Adams Cvode.Functional Cvode.default_tolerances
                   f ~roots:(n_roots, g) y_nvec

let trace = ref false
let log = ref false
let show = ref true
let delay = ref true

let args = [
    ("-trace", Arg.Set trace, "Show a trace of ball positions.");
    ("-d", Arg.Set_float t_delta, "Set the default time step.");
    ("-noshow", Arg.Clear show, "Disable the graphical display.");
    ("-nodelay", Arg.Clear delay, "No delays between frames.");
    ("-log",
     Arg.Unit (fun () -> log := true; show := false; delay := false),
     "Log state variables to stdout (implies -noshow and -nodelay).");
]

let _ =
  Arg.parse args (fun _ -> ()) "ball: simulate a ball bouncing down steps using sundials";
  if !show then Showball.start !trace !t_delta (ground, ground_limits);
  if !log then Sundials.RealArray.print_with_time 0.0 y;
  let t = ref !t_delta in
  while (y.{xpos_i} < x_limit) do
    let (t', result) = Cvode.solve_normal s !t y_nvec in
        if (result = Sundials.RootsFound) then ball_event s t' y;

        if !log then Sundials.RealArray.print_with_time t' y;
        if !show then Showball.show (y.{xpos_i}, y.{ypos_i});

        t := t' +. !t_delta
  done;
  if !show then Showball.stop ()

