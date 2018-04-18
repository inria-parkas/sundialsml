
open Sundials
open Bigarray

(* Custom implementation using OCaml arrays of arrays *)

module Custom = struct

  let m_size a = (Array.length a, Array.length a.(0))

  let m_copy a b =
    let m, n = m_size a in
    Array.iteri (fun i c1 -> Array.blit c1 0 b.(i) 0 n) a

  let m_clone a =
    let m, n = m_size a in
    let b = Array.make_matrix m n 0.0 in
    m_copy a b;
    b

  let m_zero a =
    let m, n = m_size a in
    for i = 0 to m - 1 do
      for j = 0 to m - 1 do
        a.(i).(j) <- 0.0
      done
    done

  let m_scale_add c a b =
    let m, n = m_size a in
    for i = 0 to m - 1 do
      for j = 0 to m - 1 do
        a.(i).(j) <- c *. a.(i).(j) +. b.(i).(j)
      done
    done

  let m_scale_addi c a =
    let m, n = m_size a in
    for i = 0 to m - 1 do
      for j = 0 to m - 1 do
        a.(i).(j) <- c *. a.(i).(j) +. (if i = j then 1.0 else 0.0)
      done
    done

  let m_matvec a x y =
    let m, n = m_size a in
    for i = 0 to m - 1 do
      y.{i} <- 0.0;
      for j = 0 to n - 1 do
        y.{i} <- y.{i} +. a.(i).(j) *. x.{j}
      done
    done

  let m_space a =
    let m, n = m_size a in
    (m * n * 2, m + n)

  let wrap = Matrix.wrap_custom
    {
      Matrix.m_clone      = m_clone;
      Matrix.m_zero       = m_zero;
      Matrix.m_copy       = m_copy;
      Matrix.m_scale_add  = m_scale_add;
      Matrix.m_scale_addi = m_scale_addi;
      Matrix.m_matvec     = m_matvec;
      Matrix.m_space      = m_space;
    }

  let make m n x = wrap (Array.make_matrix m n x)
  let create m n = make m n 0.0

  let get a i j = (Matrix.unwrap a).(i).(j)
  let set a i j x = (Matrix.unwrap a).(i).(j) <- x
  let size a = m_size (Matrix.unwrap a)

  let dumb_realarray2_adapter f ma =
    let a = Matrix.unwrap ma in
    let m, n = m_size a in
    let ra = Sundials.RealArray2.create m n in
    for i = 0 to m - 1 do
      for j = 0 to n - 1 do
        Sundials.RealArray2.set ra i j a.(i).(j)
      done
    done;
    f ra;
    for i = 0 to m - 1 do
      for j = 0 to n - 1 do
        a.(i).(j) <- Sundials.RealArray2.get ra i j
      done
    done

  let getrf a p =
    dumb_realarray2_adapter (fun a -> Matrix.ArrayDense.getrf a p) a

  let getrs a p s =
    dumb_realarray2_adapter (fun a -> Matrix.ArrayDense.getrs a p s) a

end

module M = Custom
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
    Matrix.matvec a (Nvector_serial.wrap x) (Nvector_serial.wrap y);
    printf "matvec: y=@\n%a@\n\n" print_vec y
  with NotImplementedBySundialsVersion -> ());

  let b = M.create nrows ncols in
  Matrix.blit a b;

  Matrix.scale_add 2.0 b zero;
  printf "scale copy x2: b=@\n%a@\n" print_mat b;

  Matrix.scale_addi 1.0 b;
  printf "add identity: b=@\n%a@\n" print_mat b;

  let p = LintArray.create nrows in
  Array1.fill p 0;
  M.getrf a p;
  printf "getrf: a=@\n%a@\n" print_mat a;
  printf "       p=@\n%a@\n@\n" print_p p;

  let s = RealArray.create nrows in
  s.{0} <-  5.0;
  s.{1} <- 18.0;
  s.{2} <-  6.0;
  M.getrs a p s;
  printf "getrs: s=@\n%a@\n" print_vec s;
  ();;

main ();;
Gc.compact ();;

