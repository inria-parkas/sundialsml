
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

let nrows, ncols = 3, 3;;

let main () =
  let a = M.create nrows ncols in
  let zero = M.make nrows ncols 0.0 in

  M.set a 0 0 ( 1.0);
  M.set a 0 1 ( 2.0);
  M.set a 0 2 ( 3.0);

  M.set a 1 0 ( 2.0);
  M.set a 1 1 (-4.0);
  M.set a 1 2 ( 6.0);

  M.set a 2 0 ( 3.0);
  M.set a 2 1 (-9.0);
  M.set a 2 2 (-3.0);

  printf "initially: a=@\n%a@\n" print_mat a;

  (try
    let x = RealArray.of_array [| 1.0; 2.0; 3.0 |] in
    let y = RealArray.create nrows in
    M.matvec a x y;
    printf "matvec: y=@\n%a@\n\n" print_vec y
  with Config.NotImplementedBySundialsVersion -> ());

  let b = M.create nrows ncols in
  M.blit ~src:a ~dst:b;

  M.scale_add 2.0 b zero;
  printf "scale copy x2: b=@\n%a@\n" print_mat b;

  M.scale_addi 1.0 b;
  printf "add identity: b=@\n%a@\n" print_mat b;

  let a_ra2 = Sundials.RealArray2.wrap (M.unwrap a) in

  let p = LintArray.create nrows in
  Array1.fill p 0;
  Matrix.ArrayDense.getrf a_ra2 p;
  printf "getrf: a=@\n%a@\n" print_mat a;
  printf "       p=@\n%a@\n@\n" print_p p;

  let s = RealArray.create nrows in
  s.{0} <-  5.0;
  s.{1} <- 18.0;
  s.{2} <-  6.0;
  Matrix.ArrayDense.getrs a_ra2 p s;
  printf "getrs: s=@\n%a@\n" print_vec s;
  ();;

main ();;
Gc.compact ();;

