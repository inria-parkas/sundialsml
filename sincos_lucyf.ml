
module Cvode = Cvode_serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

open Solvelucy

let f mode rin y der rout =
  match mode with
  | Init -> begin (* calculate: y *)
        print_endline "init"; (* XXX *)

        y.{0} <- 0.0;
        y.{1} <- 0.0;
        y.{2} <- 0.0;
        true
      end

  | Discrete -> begin (* using: rin, calculate: y, rout *)
        print_endline "discrete"; (* XXX *)
        true
      end

  | Continuous -> begin (* using: y, calculate: der, rout *)
        print_endline "continuous"; (* XXX *)

        der.{0} <- cos(y.{1});
        der.{1} <- 1.0;
        der.{2} <- 0.0;

        rout.{0} <- y.{0};
        rout.{1} <- y.{1};
        true
      end

let _ = sundialify (Some 10.0) f (fun t -> t +. 0.1) 3 2

