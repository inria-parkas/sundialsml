
module Cvode = Cvode_serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

(*
 * Albert and Benoît's 'reasonable' sliding mode control example
 *
 * der(x) = 0 init - sgn(y0) reset
 *                           | -1 every up(y)
 *                           |  1 every up(-y)
 * der(y) = x init y0
 *
 *)

let y0 = ref (-1.0)
let multiple_discrete = ref (true)

let sgn x = if x < 0.0 then -1.0
            else if x > 0.0 then 1.0
            else 0.0

(* index elements of v and der *)
let states = [| "x"; "y" |]
let n_eq = Array.length states
let x = 0
and y = 1

(* index elements of up and up_e *)
let roots = [| "up(y)"; "up(-y)" |]
let n_zc = Array.length roots
and zc_y  = 0       (* up(y)  *)
and zc_my = 1       (* up(-y) *)


let f init      (* boolean: true => initialization *)
      up_arr    (* array of booleans: zero-crossings, value of up() *)
      v         (* array of floats: continuous state values *)
      der       (* array of floats: continuous state derivatives *)
      up_e =    (* array of floats: value of expressions inside up() *)
  begin
    if init then
      begin    (* initialization: calculate v *)
        v.{x} <- -. (sgn !y0);  (* x: init - sgn(y0) *)
        v.{y} <- !y0            (* y: init y0 *)
      end
    else
    if Roots.exists up_arr
    then begin (* discrete mode: using up, calculate v *)
      let up = Roots.get up_arr in

      v.{x} <- (if up(zc_y) then -1.0       (* -1 every up(y)  *)
                else if up(zc_my) then 1.0  (*  1 every up(-y) *)
                else v.{x})                 (* unchanged *)
    end
    else begin (* continuous mode: using v, calculate der *)
      der.{x} <- 0.0;           (* der(x) = 0 *)
      der.{y} <- v.{x}          (* der(y) = x *)
    end
  end;
  begin        (* discrete and continuous: calculate up_e *)
    up_e.{zc_y}  <- v.{y};      (* up(y)  *)
    up_e.{zc_my} <- (-.v.{y})   (* up(-y) *)
  end;
  true

let args =
  [
    ("-y0", Solvelucy.set_float_delta y0, "initial value of y");
    ("-single", Arg.Clear multiple_discrete,
     "restrict to 1 discrete step between continuous steps");
  ]

let _ = Arg.parse (args @ Solvelucy.args n_eq) (fun _ -> ())
        "nontordu3: non-standard chattering"

let _ =
  Cvode.extra_time_precision := true;

  if !multiple_discrete
  then print_endline "! allow multiple discrete steps: (C+D+C+)*\n\n"
  else print_endline "! single discrete step (C+DC+)*";

  Solvelucy.enable_logging ();
  Solvelucy.run !multiple_discrete f None states roots

