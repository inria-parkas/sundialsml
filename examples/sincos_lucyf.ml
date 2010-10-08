
module Cvode = Cvode_serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

let states = [| "y.{0}"; "y.{1}"; "y.{2}" |]
let n_eq = Array.length states
let roots = [| "rout.{0}"; "rout.{1}" |]
let n_zc = Array.length roots

let f init rin y der rout =
  if init then
    begin (* calculate: y *)
      y.{0} <- 0.0;
      y.{1} <- 0.0;
      y.{2} <- 0.0;
      true
    end
  else if Roots.get rin 0 || Roots.get rin 1 then
    begin (* using: rin, calculate: y, rout *)
      true
    end
  else begin (* using: y, calculate: der, rout *)
      der.{0} <- cos(y.{1});
      der.{1} <- 1.0;
      der.{2} <- 0.0;

      rout.{0} <- y.{0};
      rout.{1} <- y.{1};
      true
    end

let _ = Solvelucy.max_sim_time := Some 10.0;
        Arg.parse (Solvelucy.args n_eq) (fun _ -> ())
        "sincos_lucyf: simple sinusoidal output"

let _ =
  Solvelucy.enable_logging ();
  Solvelucy.run_delta f None states roots

