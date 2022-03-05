(*
  A double pendulum from Wikipedia:
  https://en.wikipedia.org/wiki/Double_pendulum
 *)

let pi = Float.pi
let pi2 = 2.0 *. Float.pi
let g = 9.81

let width = 400
let height = 400

type pendulum = {
  mutable thetai : float;
  mutable length : float;
  mutable mass   : float;
}

let p1 = {
  thetai = pi2 /. 8.0;
  length = 1.0;
  mass   = 1.0;
}

let p2 = {
  thetai = 0.0;
  length = 1.0;
  mass   = 1.0;
}

let color1 = Graphics.rgb 30 80 210   (* blue *)
let color2 = Graphics.rgb 15 200 20   (* green *)
let length_factor = 0.45 *. float height /. 2.0

let show_graphics () =
  Graphics.open_graph "";
  Graphics.set_window_title "Double Pendulum";
  Graphics.resize_window width height

let draw_pendulum theta1 theta2 =
  let x1, y1 = width / 2, height / 2 in
  let x2, y2 = x1 + int_of_float (p1.length *. length_factor *. sin theta1),
               y1 - int_of_float (p1.length *. length_factor *. cos theta1)
  in
  let x3, y3 = x2 + int_of_float (p2.length *. length_factor *. sin theta2),
               y2 - int_of_float (p2.length *. length_factor *. cos theta2)
  in
  Graphics.clear_graph ();
  Graphics.set_line_width 5;
  Graphics.moveto x1 y1;
  Graphics.set_color color1;
  Graphics.draw_circle x1 y1 1;
  Graphics.lineto x2 y2;
  Graphics.set_color color2;
  Graphics.lineto x3 y3;
  Graphics.draw_circle x3 y3 1

let draw_axis x y =
  Graphics.moveto (x - 20) y;
  Graphics.lineto (x + 20) y;
  Graphics.moveto x y;
  Graphics.lineto x (y - 30)

let draw_pendulum_with_parameters theta1 theta2 =
  draw_pendulum theta1 theta2;
  let x1, y1 = width / 2, height / 2 in
  let x2, y2 = x1 + int_of_float (p1.length *. length_factor *. sin theta1),
               y1 - int_of_float (p1.length *. length_factor *. cos theta1)
  in
  (* draw parameters *)
  Graphics.set_line_width 1;
  Graphics.(set_color black);
  draw_axis x1 y1;
  draw_axis x2 y2

let theta1, theta2, ptheta1, ptheta2 = 0, 1, 2, 3

module F = struct
  let ( + ) = Float.add
  let ( - ) = Float.sub
  let ( * ) = Float.mul
  let ( / ) = Float.div
end

let sqr x = x *. x

let f _t s sd =
  let open! F in
  sd.{theta1} <-
    (6. / p1.mass * sqr p1.length)
    *
    ((2. * s.{ptheta1} - 3. * cos (s.{theta1} - s.{theta2}) * s.{ptheta2})
     / (16. - 9. * sqr (cos (s.{theta1} - s.{theta2}))));
  sd.{theta2} <-
    (6. / p2.mass * sqr p2.length)
    *
    ((8. * s.{ptheta2} - 3. * cos (s.{theta1} - s.{theta2}) * s.{ptheta1})
     / (16. - 9. * sqr (cos (s.{theta1} - s.{theta2}))));
  sd.{ptheta1} <-
    -0.5 * p1.mass * sqr p1.length
    * (sd.{theta1} * sd.{theta2} * sin (s.{theta1} - s.{theta2})
       + 3. * (g / p1.length) * sin s.{theta1});
  sd.{ptheta2} <-
    -0.5 * p2.mass * sqr p2.length
    * (-. sd.{theta1} * sd.{theta2} * sin (s.{theta1} - s.{theta2})
       + (g / p2.length) * sin s.{theta2})

let initial_sleep = ref 1

let simulate () =
  let ya = Sundials.RealArray.of_list [ p1.thetai; p2.thetai; 0.; 0. ] in
  let y = Nvector_serial.wrap ya in

  let print_state t =
    Printf.printf "% 10.7f\t% 10.7f\t% 10.7f\t% 10.7f\t% 10.7f\n%!"
      t ya.{theta1} ya.{theta2} ya.{ptheta1} ya.{ptheta2};
  in

  (* advance solver step-by-step *)
  let rec go init s (t, r) =
    if not init then print_state t;
    draw_pendulum ya.{theta1} ya.{theta2};
    Unix.sleepf 0.05;
    match r with
    | Cvode.Success -> go false s (Cvode.solve_normal s (t +. 0.05) y)
    | Cvode.RootsFound -> assert false
    | Cvode.StopTimeReached -> ()
  in

  (* create a session with the cvode numeric solver *)
  let s = Cvode.(init Adams default_tolerances f 0.0 y) in
  Cvode.set_stop_time s 10.0;

  (* show initial state with pause *)
  Printf.printf "  t               theta1          theta2          ptheta1         ptheta2\n";
  print_state 0.0;
  show_graphics ();
  draw_pendulum_with_parameters p1.thetai p2.thetai;
  Unix.sleep !initial_sleep;

  go true s (0.0, Cvode.Success)

(* command-line driver *)

let set_init p v = (p.thetai <- pi2 *. v)
let set_length p v = (p.length <- v)
let set_mass p v = (p.mass <- v)

let _ =
  Arg.parse [
    ("-i1", Arg.Float (set_init p1), "initial angle for upper pendulum (2π * i1)");
    ("-i2", Arg.Float (set_init p2), "initial angle for lower pendulum (2π * i2)");
    ("-l1", Arg.Float (set_length p1), "length of upper pendulum (0, 1.0]");
    ("-l2", Arg.Float (set_length p2), "length of lower pendulum (0, 1.0]");
    ("-m1", Arg.Float (set_mass p1), "mass of upper pendulum");
    ("-m2", Arg.Float (set_mass p2), "mass of upper pendulum");
    ("-p",  Arg.Set_int initial_sleep, "initial sleep time");
  ] (fun _ -> ()) "dpendulum";
  simulate ()

