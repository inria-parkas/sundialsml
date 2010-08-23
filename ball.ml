
let ypos_i = 0
let yvel_i = 1
let yacc_i = 2
let xpos_i = 3

let under_i = 0

let gravity = -9.81 (* m/s/s *)
let t_delta = 0.05  (* s *)
let x_vel   = 0.8   (* m/s *)
let x_limit = 15.0  (* m *)

let real_time_delay () =
  Unix.sleep 1

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

let y = Cvode_serial.create 4
let _ = y.{xpos_i} <- 0.0;
        y.{ypos_i} <- 10.0;
        y.{yvel_i} <- 0.0;
        y.{yacc_i} <- gravity

let rootdata = Cvode_serial.int_array 1
let ball_event s t y =
  Cvode_serial.get_roots s rootdata;

  if (rootdata.{under_i} != 0l && y.{yvel_i} <= 0.0) then
    (print_endline "hit ground!";
     y.{yvel_i} <- (-0.8 *. y.{yvel_i});
     Cvode_serial.reinit s (t +. t_delta (* XXX NB XXX *)) y)

let s = Cvode_serial.init f (1, g) y

let _ =
  Showball.start true (ground, ground_limits);
  Cvode_serial.print_results 0.0 y;
  let t = ref t_delta in
  while (y.{xpos_i} < x_limit) do
    let (t', roots) = Cvode_serial.advance s !t y in
        Cvode_serial.print_results t' y;
        Showball.show (y.{xpos_i}, y.{ypos_i});
        real_time_delay ();
        if (roots) then
          (ball_event s t' y;
           t := t' +. t_delta +. t_delta) (* XXX NB XXX *)
        else if (t' >= !t) then
          t := !t +. t_delta
  done;
  Unix.sleep 30;
  Showball.stop ()

let _ = Cvode_serial.free s

