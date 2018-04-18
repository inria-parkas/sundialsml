
open Sundials
open Bigarray

module M = Matrix.Dense
let printf = Format.printf
let fprintf = Format.fprintf

let print_mat out m =
  let (nr, nc) = M.size m in
  fprintf out "@[<v>";
  for i = 0 to nr - 1 do
    fprintf out "@[<h>";
    for j = 0 to nc - 1 do
      fprintf out "@ % e" (M.get m i j)
    done;
    fprintf out "@]@ "
  done;
  fprintf out "@]"

let print_vec out m =
  let nr = Array1.dim m in
  fprintf out "@[<h>";
  for i = 0 to nr - 1 do
    fprintf out "@ % e" m.{i}
  done;
  fprintf out "@]"

let print_p out p =
  let nc = Array1.dim p in
  fprintf out "@[<h>";
  for i = 0 to nc - 1 do
    fprintf out "@ % d" p.{i}
  done;
  fprintf out "@]"

let nrows, ncols = 5, 5
let a = M.create nrows ncols;;

M.set a 0 0 ( 1.0);
M.set a 0 1 ( 2.0);
M.set a 0 2 ( 0.0);
M.set a 0 3 ( 0.0);
M.set a 0 4 ( 0.0);

M.set a 1 0 ( 2.0);
M.set a 1 1 ( 2.0);
M.set a 1 2 ( 3.0);
M.set a 1 3 ( 0.0);
M.set a 1 4 ( 0.0);

M.set a 2 0 ( 0.0);
M.set a 2 1 ( 3.0);
M.set a 2 2 ( 3.0);
M.set a 2 3 ( 4.0);
M.set a 2 4 ( 0.0);

M.set a 3 0 ( 0.0);
M.set a 3 1 ( 0.0);
M.set a 3 2 ( 4.0);
M.set a 3 3 ( 4.0);
M.set a 3 4 ( 5.0);

M.set a 4 0 ( 0.0);
M.set a 4 1 ( 0.0);
M.set a 4 2 ( 0.0);
M.set a 4 3 ( 5.0);
M.set a 4 4 ( 5.0);

printf "initially: a=@\n%a@\n" print_mat a;;

let a_ra2 = Sundials.RealArray2.wrap (M.unwrap a) in

let p = LintArray.create nrows in
Array1.fill p 0;
Matrix.ArrayDense.getrf a_ra2 p;
printf "getrf: a=@\n%a@\n" print_mat a;
printf "       p=@\n%a@\n@\n" print_p p;

let s = RealArray.create nrows in
s.{0} <-  5.0;
s.{1} <- 15.0;
s.{2} <- 31.0;
s.{3} <- 53.0;
s.{4} <- 45.0;
Matrix.ArrayDense.getrs a_ra2 p s;
printf "getrs: s=@\n%a@\n" print_vec s;;

