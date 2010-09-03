
module Cvode = Cvode_serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

(*
 * Example 'non-tordu' of Albert and Benoît
 *
 * der(y) = 0 init -1 reset
 *                    |  1 every up(x)
 *                    | -1 every up(-x)
 *
 * der(x) = 0 init -1 reset
 *                    | -1 every up(y)
 *                    |  1 every up(-y)
 *                    |  1 every up(z)
 *
 * der(z) = 1 init -1
 *
 *)

(*
 * Note this example cannot be compiled using the `synchronous semantics' (where
 * up values are calculated internally rather than in the discrete solver),
 * because there are cyclic dependencies:
 *     v.{y} depends on up(zc_x)
 *     which depends on v.{x}
 *     which depends on up(zc_y)
 *     which depends on v.{y}
 *
 *)

let max_sim_time = 5.0
let max_step_size = 0.1

(* index elements of v and der *)
let x = 0
and y = 1
and z = 2
and n_eq = 3

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
        v.{y} <- (-1.0);        (* y: init -1 *)
        v.{x} <- (-1.0);        (* x: init -1 *)
        v.{z} <- (-1.0);        (* z: init -1 *)
      end
    else
    if Roots.exists up_arr
    then begin (* discrete mode: using up, calculate v *)
      let up = Roots.get up_arr in

      v.{y} <- (if up(zc_x) then 1.0           (*  1 every up(x)  *)
                else if up(zc_mx) then -1.0    (* -1 every up(-x) *)
                else v.{y});                   (* unchanged *)

      v.{x} <- (if up(zc_y) then -1.0          (* -1 every up(y)  *)
                else if up(zc_my) then 1.0     (*  1 every up(-y) *)
                else if up(zc_z) then 1.0      (*  1 every up(z)  *)
                else v.{x})                    (* unchanged *)
    end
    else begin (* continuous mode: using v, calculate der *)
      der.{y} <- 0.0;           (* der(y) = 0 *)
      der.{x} <- 0.0;           (* der(x) = 0 *)
      der.{z} <- 1.0            (* der(z) = 1 *)
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

and zc_y  = 0       (* up(y)  *)
and zc_my = 1       (* up(-y) *)
and zc_z  = 2       (* up(z)  *)
and zc_x  = 3       (* up(x)  *)
and zc_mx = 4       (* up(-x) *)
let _ =
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

let _ = print_endline "        time\t      x\t\t      y\t\t      z"
let _ = Solvelucy.run_delta
          (Some max_sim_time) f (fun t -> t +. max_step_size) n_eq n_zc

