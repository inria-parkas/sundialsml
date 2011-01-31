
module Cvode = Cvode.Serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

(* Simple example of cascaded zero-crossings
 *
 * der(x) = rx init x0 reset 0 every up(y)
 * der(y) = ry init y0 reset 0 every up(z)
 * der(z) = 1 init z0 
 *)

(* initial values *)
let x0 = ref (-1.0)
and y0 = ref (-1.0)
and z0 = ref (-1.0)
let rx = ref (0.0)
and ry = ref (1.0)

(* index elements of v and der *)
let states = [| "x"; "y"; "z" |]
let n_eq = Array.length states
and x = 0
and y = 1
and z = 2

(* index elements of up and up_e *)
let roots = [| "up(y)"; "up(z)" |]
let n_zc = Array.length roots
and zc_y = 0
and zc_z = 1

let f init      (* boolean: true => initialization *)
      up_arr    (* array of booleans: zero-crossings, value of up() *)
      v         (* array of floats: continuous state values *)
      der       (* array of floats: continuous state derivatives *)
      up_e =    (* array of floats: value of expressions inside up() *)
  begin
    if init then
      begin    (* initialization: calculate v *)
        v.{x} <- !x0;
        v.{y} <- !y0;
        v.{z} <- !z0
      end
    else
    if Roots.exists up_arr
    then begin (* discrete mode: using up, calculate v *)
      let up = Roots.get up_arr in
      v.{x} <- (if up(zc_y) then !rx else v.{x});
      v.{y} <- (if up(zc_z) then !ry else v.{y});
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
    ("-x0", Solvelucy.set_float_delta x0, "initial value of x");
    ("-y0", Solvelucy.set_float_delta y0, "initial value of y");
    ("-z0", Solvelucy.set_float_delta z0, "initial value of z");
    ("-rx", Solvelucy.set_float_delta rx, "reset value for x");
    ("-ry", Solvelucy.set_float_delta ry, "reset value for y (try 0.0)");
  ]

let _ = Arg.parse (args @ Solvelucy.args n_eq) (fun _ -> ())
        "cascade: simple zero-crossing cascade"

let _ =
  Solvelucy.max_sim_time := Some 2.0;
  Solvelucy.enable_logging ();
  (* Solvelucy.enable_zeroc_logging (); *)
  Solvelucy.run_delta f None states roots

