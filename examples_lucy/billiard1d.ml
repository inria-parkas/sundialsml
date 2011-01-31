
module Cvode = Cvode.Serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

(*
 * Benoît's 'billiard à 1 dimension'
 *
 * assume:
 *   d1 <= d2 <= 0
 *
 * der(x1) = v1 init d1
 * der(x2) = v2 init d2
 * der(v1) = 0 init w1 reset
 *                     | last(v2) when up(x1 - x2)
 * der(v2) = 0 init w2 reset
 *                     | last(v1) when up(x1 - x2)
 *                     | -last(v2) when up(x2)
 *)

(* initial values *)
let d1 = ref (-5.0)
and w1 = ref ( 2.0)
and d2 = ref (-3.0)
and w2 = ref ( 1.0)

(* index elements of v and der *)
let states = [| "x1"; "x2"; "v1"; "v2" |]
let n_eq = Array.length states
and x1 = 0
and x2 = 1
and v1 = 2
and v2 = 3

(* index elements of up and up_e *)
let roots = [| "up(x1 - x2)"; "up(x2)" |]
let n_zc = Array.length roots
and zc_x1minx2  = 0 (* up(x1 - x2)  *)
and zc_x2 = 1       (* up(x2) *)

let f init      (* boolean: true => initialization *)
      up_arr    (* array of booleans: zero-crossings, value of up() *)
      v         (* array of floats: continuous state values *)
      der       (* array of floats: continuous state derivatives *)
      up_e =    (* array of floats: value of expressions inside up() *)
  begin
    if init then
      begin    (* initialization: calculate v *)
        v.{x1} <- !d1;
        v.{x2} <- !d2;
        v.{v1} <- !w1;
        v.{v2} <- !w2
      end
    else
    if Roots.exists up_arr
    then begin (* discrete mode: using up, calculate v *)
      let up = Roots.get up_arr in
      let last_v1 = v.{v1} in

      v.{v1} <- (if up(zc_x1minx2) then (* last *) v.{v2} else v.{v1});
      v.{v2} <- (if up(zc_x1minx2) then last_v1
                 else if up(zc_x2) then -. v.{v2}
                 else v.{v2})
    end
    else begin (* continuous mode: using v, calculate der *)
      der.{x1} <- v.{v1};
      der.{x2} <- v.{v2};
      der.{v1} <- 0.0;
      der.{v2} <- 0.0
    end
  end;
  begin        (* discrete and continuous: calculate up_e *)
    up_e.{zc_x1minx2} <- v.{x1} -. v.{x2};
    up_e.{zc_x2} <- v.{x2}
  end;
  true

let args =
  [
    ("-d1", Solvelucy.set_float_delta d1, "initial position of ball 1");
    ("-w1", Solvelucy.set_float_delta w1, "initial velocity of ball 1");
    ("-d2", Solvelucy.set_float_delta d2, "initial position of ball 2");
    ("-w2", Solvelucy.set_float_delta w2, "initial velocity of ball 2");
  ]

let _ = Arg.parse (args @ Solvelucy.args n_eq) (fun _ -> ())
        "billiard1d: 1-dimensional billiard balls"

let _ =
  Solvelucy.enable_logging ();
  (* Solvelucy.enable_zeroc_logging (); *)
  Solvelucy.run_delta f None states roots

