
open Sundials
open Bigarray

module M = Dls.ArrayBandMatrix
let printf = Format.printf
let fprintf = Format.fprintf

let print_mat_data out m =
  let d = RealArray2.unwrap m in
  let nc = Array2.dim1 d in
  let nr = Array2.dim2 d in
  fprintf out "@[<v>";
  for j = 0 to nc - 1 do
    fprintf out "@[<h>";
    for i = 0 to nr - 1 do
      fprintf out "@ % e" (d.{j, i})
    done;
    fprintf out "@]@ "
  done;
  fprintf out "@]"

let print_mat mu ml smu out m =
  let (nr, nc) = RealArray2.size m in
  fprintf out "@[<v>";
  for i = 0 to nr - 1 do
    fprintf out "@[<h>";
    for j = 0 to nr - 1 do (* square *)
      if (i > j + ml) || (j > i + mu)
      then fprintf out "       --     "
      else fprintf out "@ % e" (M.get m smu i j)
    done;
    fprintf out "@]@ "
  done;
  fprintf out "@]"

let print_factored_mat mu ml smu out m =
  let (nr, nc) = RealArray2.size m in
  fprintf out "@[<v>";
  for i = 0 to nr - 1 do
    fprintf out "@[<h>";
    for j = 0 to nr - 1 do (* square *)
      if (j > i + mu) && (j <= i + smu)
      then fprintf out "@ (% e)" (M.get m smu i j)
      else if (i > j + ml) || (j  > i + mu)
      then fprintf out "        --      "
      else fprintf out "@  % e " (M.get m smu i j)
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

let n =5
and mu = 1
and ml = 1;;
let smu = min (n - 1) (mu + ml);;

let a = M.make n smu ml;;
Array2.fill (RealArray2.unwrap a) 0.0;;

M.set a smu 0 0 ( 1.0);
M.set a smu 0 1 ( 2.0);

M.set a smu 1 0 ( 2.0);
M.set a smu 1 1 ( 2.0);
M.set a smu 1 2 ( 3.0);

M.set a smu 2 1 ( 3.0);
M.set a smu 2 2 ( 3.0);
M.set a smu 2 3 ( 4.0);

M.set a smu 3 2 ( 4.0);
M.set a smu 3 3 ( 4.0);
M.set a smu 3 4 ( 5.0);

M.set a smu 4 3 ( 5.0);
M.set a smu 4 4 ( 5.0);

printf "initially: a.data=@\n%a@\n" print_mat_data a;;

printf "initially: a=@\n%a@\n" (print_mat mu ml smu) a;;

let b = RealArray2.copy a;;

M.scale 2.0 b mu ml smu;
printf "scale copy x2: b=@\n%a@\n" (print_mat mu ml smu) b;;

M.add_identity b smu;
printf "add identity: b=@\n%a@\n" (print_mat mu ml smu) b;;

let p = LintArray.make 5;;
Array1.fill p 0;
M.gbtrf a mu ml smu p;
printf "getrf: a=@\n%a@\n" (print_factored_mat mu ml smu) a;
printf "       p=@\n%a@\n@\n" print_p p;;

let s = RealArray.make n;;
s.{0} <-  5.0;
s.{1} <- 15.0;
s.{2} <- 31.0;
s.{3} <- 53.0;
s.{4} <- 45.0;
M.gbtrs a smu ml p s;
printf "getrs: s=@\n%a@\n" print_vec s;;

