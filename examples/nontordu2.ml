
module Cvode = Cvode_serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

(*
 * Example 'non-tordu' of Albert and Benoît
 *
 * der(ydot) = 0 init -1 reset
 *                       |  1 every up(x)
 *                       | -1 every up(-x)
 * der(y) = ydot init 0
 *
 * der(xdot) = 0 init -1 reset
 *                       | -1 every up(y)
 *                       |  1 every up(-y)
 *                       |  1 every up(z)
 * der(x) = xdot init 0
 *
 * der(z) = 1 init -1
 *
 *)

(*
 * Note this example can be compiled using the `synchronous semantics' (where up
 * values are calculated internally rather than in the discrete solver), because
 * there are no cyclic dependencies.
 *)

let max_sim_time = 20.0
let max_step_size = 0.1

(* index elements of v and der *)
let x    = 0
and y    = 1
and z    = 2
and xdot = 3
and ydot = 4
and n_eq = 5

(* index elements of up and up_e *)
and zc_y  = 0       (* up(y)  *)
and zc_my = 1       (* up(-y) *)
and zc_z  = 2       (* up(z)  *)
and zc_x  = 3       (* up(x)  *)
and zc_mx = 4       (* up(-x) *)

and n_zc = 5

let f init      (* boolean: true => initialization *)
      up_arr    (* array of booleans: zero-crossings, value of up() *)
      v         (* array of floats: continuous state values *)
      der       (* array of floats: continuous state derivatives *)
      up_e =    (* array of floats: value of expressions inside up() *)
  begin
    if init then
      begin    (* initialization: calculate v *)
        v.{y}    <-  (0.0);        (* y: init 0 *)
        v.{x}    <-  (0.0);        (* x: init 0 *)
        v.{z}    <- (-1.0);        (* z: init -1 *)
        v.{xdot} <- (-1.0);        (* xdot: init -1 *)
        v.{ydot} <- (-1.0)         (* ydot: init -1 *)
      end
    else
    if Roots.exists up_arr
    then begin (* discrete mode: using up, calculate v *)
      let up = Roots.get up_arr in
      
      v.{ydot} <- (if up(zc_x) then 1.0           (*  1 every up(x)  *)
                   else if up(zc_mx) then -1.0    (* -1 every up(-x) *)
                   else v.{ydot});                (* unchanged *)

      v.{xdot} <- (if up(zc_y) then -1.0          (* -1 every up(y)  *)
                   else if up(zc_my) then 1.0     (*  1 every up(-y) *)
                   else if up(zc_z) then 1.0      (*  1 every up(z)  *)
                   else v.{xdot})                 (* unchanged *)
    end
    else begin (* continuous mode: using v, calculate der *)
      der.{y} <- v.{ydot};      (* der(y) = ydot *)
      der.{x} <- v.{xdot};      (* der(x) = xdot *)
      der.{z} <- 1.0;           (* der(z) = 1 *)
      der.{xdot} <- 0.0;        (* der(xdot) = 0 *)
      der.{ydot} <- 0.0         (* der(ydot) = 0 *)
    end
  end;
  begin        (* discrete and continuous: calculate up_e *)
    up_e.{zc_y}  <- v.{y};      (* up(y)  *)
    up_e.{zc_my} <- (-.v.{y});  (* up(-y) *)
    up_e.{zc_z}  <- v.{z};      (* up(z)  *)
    up_e.{zc_x}  <- v.{x};      (* up(x)  *)
    up_e.{zc_mx} <- (-.v.{x})   (* up(-x) *)
  end;
  true

let _ = Arg.parse (Solvelucy.args n_eq) (fun _ -> ())
        "nontordu2: chattering on velocity"

let _ =
  Solvelucy.enable_logging ();
  print_endline "";
  print_endline "C: result of continuous solver";
  print_endline "D: result of discrete solver";
  print_endline "";
  print_endline "    up(y)";
  print_endline "   /  up(-y)";
  print_endline "   | /  up(z)";
  print_endline "   | | /  up(x)";
  print_endline "   | | | /  up(-x)";
  print_endline "   | | | | /";
  print_endline "R: 0 0 0 0 0";
  print_endline ""

let _ = print_endline "        time\t      x\t\t      y\t\t      z\t\t     xdot\t     ydot"
let _ = Solvelucy.run_delta
          (Some max_sim_time) f (fun t -> t +. max_step_size) n_eq n_zc

