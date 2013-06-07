
let scale = ref 1.0
let ball_radius = 4

let obstacle = ref [||]
let pivot = ref (400, 480)

let xc x = truncate (x *. !scale) + fst !pivot
let yc y = truncate (y *. !scale) + snd !pivot

let show_objects () =
  Graphics.set_color Graphics.blue;
  Graphics.fill_poly !obstacle;
  Graphics.set_color Graphics.black;
  Graphics.fill_circle (fst !pivot) (snd !pivot) (ball_radius / 2)

let leave_trace = ref false
let refresh_delay = ref (1.0 /. 60.)

let width = 800
let height = 600

(** The first argument specifies how many meters of the model the window height
    will show.  For example, a value of 2 means a 1-meter tall object is
    rendered to fill half the window height.  *)
let start scale_val timestep trace (piv_x,piv_y) (obs_x,obs_y) =
  Graphics.open_graph "";
  (*Graphics.resize_window width height;*)
  Graphics.auto_synchronize false;
  Graphics.clear_graph ();

  (* Actual size of window might be different. *)
  let width = Graphics.size_x ()
  and height = Graphics.size_y () in

  (* x[m] * height[px/h] / scale_val[m/h] = x'[px] *)
  scale := float_of_int height /. scale_val;

  let piv_x = piv_x *. float_of_int width
  and piv_y = piv_y *. float_of_int height in
  pivot := (int_of_float piv_x, int_of_float piv_y);
  (* Find the farthest point on screen from the pivot in the (obs_x,obs_y)
     direction.  *)
  let swap (x,y) = (y,x)
  and id p = p in
  let (kmin_from_x, kmax_from_x) =
    (if obs_x > 0. then id else swap)
      (-.piv_x /. obs_x, (float_of_int width -. piv_x) /. obs_x)
  and (kmin_from_y, kmax_from_y) =
    (if obs_y > 0. then id else swap)
      (-.piv_y /. obs_y, (float_of_int width -. piv_y) /. obs_y)
  in
  let kmin = max kmin_from_x kmin_from_y
  and kmax = min kmax_from_x kmax_from_y
  in
  let (x1,y1) = (int_of_float (piv_x +. obs_x *. kmin),
                 int_of_float (piv_y +. obs_y *. kmin))
  and (x2,y2) = (int_of_float (piv_x +. obs_x *. kmax),
                 int_of_float (piv_y +. obs_y *. kmax))
  in
  obstacle :=
    if y1 > y2 then [|(0,y1); (x1,y1); (x2,y2); (0,y2)|]
    else            [|(0,y2); (x1,y2); (x2,y1); (0,y1)|];

  show_objects ();
  Graphics.synchronize ();
  leave_trace := trace;
  refresh_delay := timestep

(* A trick to sleep with subsecond precision, suggested at
     http://caml.inria.fr/mantis/print_bug_page.php?bug_id=4023
   The let rec works around a defect when coupled with Graphics, noted in
     http://caml.inria.fr/pub/ml-archives/ocaml-beginners/2002/03/5984aeef10678485ba179f25a269845c.en.html
   but this workaround is also slightly defective.  Fantastic.
 *)
let rec minisleep (sec : float) =
  try ignore (Unix.select [] [] [] sec)
  with Unix.Unix_error (Unix.EINTR, _, _) -> minisleep sec

let last_x = ref 0.0
let last_y = ref 0.0

(* x,y should measure the displacement from the pivot in model-meters.  *)
let show (x, y) =
  Graphics.set_color Graphics.background;
  (* Erase the previous string.  *)
  Graphics.moveto (fst !pivot) (snd !pivot);
  Graphics.lineto (xc !last_x) (yc !last_y);
  if (not !leave_trace) then begin
    Graphics.fill_circle (xc !last_x) (yc !last_y) ball_radius;
    show_objects ()
  end;
  Graphics.set_color Graphics.black;
  Graphics.moveto (fst !pivot) (snd !pivot);
  Graphics.lineto (xc x) (yc y);
  Graphics.set_color Graphics.red;
  Graphics.fill_circle (xc x) (yc y) ball_radius;
  Graphics.synchronize ();
  minisleep !refresh_delay;
  last_x := x;
  last_y := y

let stop () =
  Graphics.set_color Graphics.white;
  Graphics.moveto 2 2;
  Graphics.draw_string "Simulation finished.  Press a key to exit.";
  Graphics.synchronize ();
  let _ = Graphics.wait_next_event [Graphics.Key_pressed] in
  Graphics.close_graph ()

