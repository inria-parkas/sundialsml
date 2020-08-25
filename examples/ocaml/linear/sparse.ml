
open Sundials
open Bigarray

module S = Matrix.Sparse
module D = Matrix.Dense
let printf = Format.printf
let fprintf = Format.fprintf

let print_mat out mat =
  let m, n = S.size mat in
  let nnz, _ = S.dims mat in
  fprintf out "@[<v>matrix (M=%d, N=%d, NNZ=%d):@\n" m n nnz;
  for i = 0 to n - 1 do
    fprintf out "  col %d:" i;
    for j = S.get_col mat i to (S.get_col mat (i + 1)) - 1 do
      let r, v = S.get mat j in
      printf " (%d: % .02e)" r v
    done;
    fprintf out "@\n"
  done;
  fprintf out "@]"

let print_vec out m =
  let nr = Array1.dim m in
  fprintf out "@[<h>";
  for i = 0 to nr - 1 do
    fprintf out "@ % e" m.{i}
  done;
  fprintf out "@]"

let nrows, ncols, nzeros = 3, 3, 5;;

let main () =
  let da = D.create nrows ncols in
  let b = S.(make CSC nrows ncols 2) in
  let zero = S.(make CSC nrows ncols 0) in
  let x = Sundials.RealArray.of_list [2.0; 3.0; 4.0] in
  let y = Sundials.RealArray.create nrows in

  D.set da 0 0 ( 0.0);
  D.set da 0 1 ( 2.0);
  D.set da 0 2 ( 0.0);

  D.set da 1 0 ( 0.0);
  D.set da 1 1 (-4.0);
  D.set da 1 2 ( 0.0);

  D.set da 2 0 ( 0.0);
  D.set da 2 1 (-9.0);
  D.set da 2 2 (-3.0);

  printf "initially da=@\n%a@\n%!" D.pp da;

  let a = S.(from_dense CSC 0.0 da) in
  printf "initially a=@\n%a@\n%a@\n%!" S.pp a print_mat a;

  S.scale_addi 1.0 a;
  printf "a + 1=@\n%a@\n" print_mat a;

  S.blit ~src:a ~dst:b;
  S.scale_add 2.0 b zero;
  printf "scale copy x2: b=@\n%a@\n" print_mat b;

  S.set_to_zero b;
  S.resize ~nnz:9 b;
  printf "set to zero (NNZ=9): b=@\n%a@\n" print_mat b;

  S.set_col b 0 0;
  S.set_col b 1 0;
  S.set_col b 2 0;
  S.set_col b 3 1;

  S.set b 0 1 7.0;

  printf "b with element [r=1, c=2] set to 7.0:@\n%a@\n" print_mat b;

  S.scale_add 1.0 a b;
  printf "a = a + b: a=@\n%a@\n" print_mat a;

  S.matvec a x y;
  printf "y = A*x: x=@\n%a@\ny=@\n%a@\n" print_vec x print_vec y;

  ();;

main ();;
Gc.compact ();;

