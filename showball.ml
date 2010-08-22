
let scale = 50.0
let ball_radius = 3

let x_off = 50
let y_off = 0

let xc x = truncate (x *. scale) + x_off
let yc y = truncate (y *. scale) + y_off

let show_floors (height, extent) =
  let maxidx = min (Array.length height) (Array.length extent) in
  let rec f min idx =
    if (idx < maxidx) then begin
      let h = yc height.(idx) in
      let max = xc extent.(idx) in
      Graphics.moveto min h;
      Graphics.lineto max h;
      f max (idx + 1)
    end
  in
  f (xc 0.0) 0

let leave_trace = ref false

let start trace floors =
  Graphics.open_graph "";
  Graphics.resize_window 800 600;
  Graphics.clear_graph ();
  Unix.sleep 2; (* TODO: why is this necessary? *)
  show_floors floors;
  leave_trace := trace

let last_x = ref 0.0
let last_y = ref 0.0

let show (x, y) =
  if (not !leave_trace) then begin
    Graphics.set_color Graphics.background;
    Graphics.fill_circle (xc !last_x) (yc !last_y) ball_radius
  end;
  Graphics.set_color Graphics.red;
  Graphics.fill_circle (xc x) (yc y) ball_radius;
  last_x := x;
  last_y := y

let stop () =
  Graphics.close_graph ()

