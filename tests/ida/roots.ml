module Ida = Ida_serial
module Carray = Ida.Carray

(* A perennial bug in the binding was to confuse number of roots with number of
   equations.  Trying it with a large root and a small y vector, or vice versa,
   will expose this kind of bug in the form of index-out-of-bounds errors.  *)

let test nvars nroots =
  let rootsfn t y y' roots =
    for i = 0 to nroots - 1 do
      roots.{i} <- 1.
    done
  and resfn t y y' res =
    for i = 0 to nvars - 1 do
      res.{i} <- y.{i} -. t
    done
  in
  let y0 = Carray.init nvars 0.
  and y0' = Carray.init nvars 1.
  in
  let ida = Ida.init_at_time (Ida.Spgmr 0) resfn (nroots, rootsfn) 0. y0 y0' in
  for t = 1 to 10 do
    ignore (Ida.solve_normal ida (float_of_int t) y0 y0')
  done

let _ = test 1000000 1                  (* many variables, few roots *)
let _ = test 1 1000000                  (* few variables, many roots *)
