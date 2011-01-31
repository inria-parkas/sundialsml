
module Cvode = Cvode.Serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

(*
 * Example 3 from the HSCC submission
 *
 * der(x) = 0 init 0 reset
 *                   | last x + 1 every up(y)
 *                   | last x + 2 every up(z)
 * der(y) = 0 init -1 reset 1 every up(z)
 * der(z) = 1 init -1
 *
 *)

let multiple_discrete = ref (true)

(* index elements of v and der *)
let states = [| "x"; "y"; "z" |]
let n_eq = Array.length states
let x = 0
and y = 1
and z = 2

(* index elements of up and up_e *)
let roots = [| "up(y)"; "up(z)" |]
let n_zc = Array.length roots
and zc_y  = 0       (* up(y)  *)
and zc_z = 1        (* up(z) *)

let f init      (* boolean: true => initialization *)
      up_arr    (* array of booleans: zero-crossings, value of up() *)
      v         (* array of floats: continuous state values *)
      der       (* array of floats: continuous state derivatives *)
      up_e =    (* array of floats: value of expressions inside up() *)
  begin
    if init then
      begin    (* initialization: calculate v *)
        v.{x} <- 0.0;
        v.{y} <- -1.0;
        v.{z} <- -1.0
      end
    else
    if Roots.exists up_arr
    then begin (* discrete mode: using up, calculate v *)
      let up = Roots.get up_arr in

      v.{x} <- (if up(zc_y) then v.{x} +. 1.0
                else if up(zc_z) then v.{x} +. 2.0
                else v.{x});

      v.{y} <- (if up(zc_z) then 1.0 else v.{y})

    end
    else begin (* continuous mode: using v, calculate der *)
      der.{x} <- 0.0;
      der.{y} <- 0.0;
      der.{z} <- 1.0
    end
  end;
  begin        (* discrete and continuous: calculate up_e *)
    up_e.{zc_y} <- v.{y};
    up_e.{zc_z} <- v.{z}
  end;
  true

let args =
  [
    ("-single", Arg.Clear multiple_discrete,
     "restrict to 1 discrete step between continuous steps");
  ]

let _ = Arg.parse (args @ Solvelucy.args n_eq) (fun _ -> ())
        "example3: double increment"

let _ =
  if !multiple_discrete
  then print_endline "! allow multiple discrete steps: (C+D+C+)*\n\n"
  else print_endline "! single discrete step (C+DC+)*";

  Solvelucy.enable_logging ();
  Solvelucy.run !multiple_discrete f None states roots

