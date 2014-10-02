
open Sundials
open Spils

let nvec = Nvector_serial.wrap

let printf = Format.printf
let fprintf = Format.fprintf

let print_vec out m =
  let nr = RealArray.length m in
  fprintf out "@[<h>";
  for i = 0 to nr - 1 do
    fprintf out "@ % e" m.{i}
  done;
  fprintf out "@]"

(* example:
     2*x0   - x1   + x2 = -1
       x0 + 2*x1   - x2 =  6
       x0   - x1 + 2*x2 = -3
 *)
let atimes x z =
  z.{0} <- 2.0 *. x.{0} -.        x.{1} +.        x.{2};
  z.{1} <-        x.{0} +. 2.0 *. x.{1} -.        x.{2};
  z.{2} <-        x.{0} -.        x.{1} +. 2.0 *. x.{2}

let main () =
  let b = RealArray.of_array [| -1.0; 6.0; -3.0 |] in

  (* SPGMR *)

  let s = SPGMR.make 3 (nvec b) in

  let x = RealArray.of_array [|  0.0; 0.0;  0.0 |] in
  let solved, res_norm, nli, nps =
      SPGMR.solve s (nvec x) (nvec b) Spils.PrecNone Spils.ModifiedGS
                  1.0e-4 0 None None atimes None in

  printf "SPGMR solution: x=@\n%a@\n" print_vec x;
  printf "  (solved=%B res_norm=%e nli=%d nps=%d)@\n" solved res_norm nli nps;

  (* SPBCG *)

  let s = SPBCG.make 3 (nvec b) in

  let x = RealArray.of_array [|  0.0; 0.0;  0.0 |] in
  let solved, res_norm, nli, nps =
      SPBCG.solve s (nvec x) (nvec b) Spils.PrecNone
                  1.0e-4 None None atimes None in

  printf "SPBCG solution: x=@\n%a@\n" print_vec x;
  printf "  (solved=%B res_norm=%e nli=%d nps=%d)@\n" solved res_norm nli nps;

  (* SPTFQMR *)

  let s = SPTFQMR.make 3 (nvec b) in

  let x = RealArray.of_array [|  0.0; 0.0;  0.0 |] in
  let solved, res_norm, nli, nps =
      SPTFQMR.solve s (nvec x) (nvec b) Spils.PrecNone
                    1.0e-4 None None atimes None in

  printf "SPTFQMR solution: x=@\n%a@\n" print_vec x;
  printf "  (solved=%B res_norm=%e nli=%d nps=%d)@\n" solved res_norm nli nps;;

main ();;
Gc.compact ();;

