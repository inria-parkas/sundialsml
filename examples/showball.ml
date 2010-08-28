
let scale = 90.0
let ball_radius = 4

let x_off = 50
let y_off = 0

let xc x = truncate (x *. scale) + x_off
let yc y = truncate (y *. scale) + y_off

let floors = ref ([| 0.0 |], [| 0.0 |])

let show_floors () =
  let (heights, extents) = !floors in
  let maxidx = min (Array.length heights) (Array.length extents) in
  let rec f idx =
    if (idx < maxidx) then
      (Graphics.lineto (Graphics.current_x ()) (yc heights.(idx));
       Graphics.lineto (xc extents.(idx)) (yc heights.(idx));
       f (idx + 1))
  in
  Graphics.set_color Graphics.blue;
  Graphics.moveto (xc 0.0) (yc heights.(0));
  f 0;
  Graphics.lineto (Graphics.current_x ()) (yc 0.0)

let leave_trace = ref false

let start trace floor_details =
  Graphics.open_graph "";
  Graphics.resize_window 800 600;
  (* Graphics.auto_synchronize false; *)
  Graphics.clear_graph ();
  Unix.sleep 1; (* TODO: why is this necessary? *)
  floors := floor_details;
  show_floors ();
  leave_trace := trace

let last_x = ref 0.0
let last_y = ref 0.0

let show (x, y) =
  if (not !leave_trace) then begin
    Graphics.set_color Graphics.background;
    Graphics.fill_circle (xc !last_x) (yc !last_y) ball_radius;
    show_floors ()
  end;
  Graphics.set_color Graphics.red;
  Graphics.fill_circle (xc x) (yc y) ball_radius;
  (* Graphics.synchronize (); *)
  last_x := x;
  last_y := y

let stop () =
  ignore (Graphics.read_key ());
  Graphics.close_graph ()

(*
let print_stats s =
  let in_stats = integrator_stats s
  (* and ls_stats = TODO *)
  (* and jac_evals = TODO CVDlsGetNumJacEvals *)
  (* and root_evals = TODO CVodeGetNumGEvals *)
  and printf = Printf.printf
  in
    printf("\nFinal Statistics:\n");
    printf "nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n"
      in_stats.steps
      in_stats.rhs_evals
      in_stats.linear_solver_setups
      ls_stats.rhs_evals
      jac_evals;
    printf "nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n\n"
      ls_stats.iterations
      ls_stats.convergence_failures
      in_stats.error_test_failures
      root_evals
*) 

