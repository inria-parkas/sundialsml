(* Basic quickcheck infrastructure.  Random generation and shrinking of basic
   data types, as well as a lazy stream implementation with a usable
   interface for shrinking.  *)
module Stream =
struct
  type 'a t = Cons of 'a * 'a t Lazy.t | Nil
  let empty = Nil
  let null = function
    | Nil -> true
    | _ -> false
  let singleton x = Cons (x, lazy Nil)
  let rec of_list = function
    | [] -> Nil
    | x::xs -> Cons (x, lazy (of_list xs))
  let rec map f = function
    | Cons (x, xs) -> Cons (f x, lazy (map f (Lazy.force xs)))
    | Nil -> Nil
  let rec iter f = function
    | Cons (x, xs) -> f x; iter f (Lazy.force xs)
    | Nil -> ()
  let to_list xs =
    let rec go acc = function
      | Nil -> List.rev acc
      | Cons (x, xs) -> go (x::acc) (Lazy.force xs)
    in go [] xs

  let append xs ys =
    let rec go = function
      | Nil -> ys
      | Cons (x,xs) -> Cons (x, lazy (go (Lazy.force xs)))
    in go xs
  let cons x xs = Cons (x, lazy xs)

  let concat xss =
    let rec go = function
      | Nil -> Nil
      | Cons (xs, xss) -> flatten xs xss
    and flatten xs xss =
      match xs with
      | Nil -> go (Lazy.force xss)
      | Cons (x,xs) -> Cons (x, lazy (flatten (Lazy.force xs) xss))
    in go xss

  let rec take n xs =
    match n, xs with
    | 0, _ -> Nil
    | n, Nil -> Nil
    | n, Cons (x,xs) -> Cons (x, lazy (take (n-1) (Lazy.force xs)))
  let rec take_while p = function
    | Cons (x, xs) when p x -> Cons (x, lazy (take_while p (Lazy.force xs)))
    | _ -> Nil
  let rec drop_while p = function
    | Cons (x, xs) when not (p x) ->
      Cons (x, lazy (take_while p (Lazy.force xs)))
    | _ -> Nil

  let rec generate f =
    match f () with
    | Some x -> Cons (x, lazy (generate f))
    | None -> Nil

  let rec find p = function
    | Nil -> None
    | Cons (x, xs) -> if p x then Some x else find p (Lazy.force xs)

  let rec find_some f = function
    | Nil -> None
    | Cons (x, xs) ->
      match f x with
      | None -> find_some f (Lazy.force xs)
      | Some x -> Some x

  let rec iterate f x = Cons (x, lazy (iterate f (f x)))
  let iterate_on x f = iterate f x

  let rec repeat f = Cons (f (), lazy (repeat f))
  let repeat_n n f = take n (repeat f)

  let guard b x  = if b then x else Nil
  let guard1 b x = guard b (singleton x)

  let rec filter p = function
    | Nil -> Nil
    | Cons (x, lazy xs) when p x -> Cons (x, lazy (filter p xs))
    | Cons (_, lazy xs) -> filter p xs

  let rec enum istart iend =
    if istart <= iend
    then Cons (istart, lazy (enum (istart + 1) iend))
    else Nil

  let length xs =
    let rec go i = function
      | Nil -> i
      | Cons (_, lazy xs) -> go (i+1) xs
    in go 0 xs
end

let size = ref 5

(* Generation *)
let gen_nat () = Random.int !size
let gen_pos_int () = gen_nat () + 1
let gen_int () = Random.int (!size * 2) - !size

let gen_float () =
  let max_abs = 2. *. float_of_int !size in
  Random.float max_abs -. max_abs /. 2.

(* This generator produces floating point values that are multiples of a fixed
   floating point value df.  With a sufficiently large df, discretized values
   makes it easy to detect problems like values being "too close" -- detecting,
   and properly mimicking, the behavior of the sundials solver would otherwise
   require writing code that depends on the details of its internals.  *)
let discrete_unit = ref (1. /. 2.**10.)
let gen_discrete_float () =
  let df = !discrete_unit in
  floor (gen_float () /. df +. 0.5) *. df

let gen_choice choices =
  choices.(Random.int (Array.length choices))

let enum istart iend =
  let rec go acc i =
    if istart <= i then go (i::acc) (i-1)
    else acc
  in go [] iend

let gen_list g = List.map (fun _ -> g ()) (enum 1 (gen_nat ()))

let gen_array ?(size=gen_pos_int ()) gen_elem =
  let v = Array.make size (gen_elem ()) in
  for i = 1 to size-1 do
    v.(i) <- gen_elem ()
  done;
  v

(* Shrinking.  Some of the algorithms are adapted from Haskell's QuickCheck.
   Keep in mind that the shrink-and-test cycle requires all shrinkers to ensure
   the outputs are strictly simpler than the input in some sense.  *)
let is_nan (x : float) = not (x = x)

let (@@) = Stream.append

let shrink_int n =
  let less_complex k = abs n > abs k in
  Stream.guard1 (n < -n) (- n)
  @@ Stream.filter less_complex
     (Stream.of_list [0;1;2]
      @@ Stream.repeat_n (gen_nat ()) (fun () ->
        int_of_float (floor (float_of_int n *. (Random.float 2. -. 1.)))))
let shrink_nat n = Stream.map abs (shrink_int n)

let shrink_float f =
  let less_complex g = abs_float f > abs_float g in
  Stream.guard1 (f < 0.) (-. f)
  @@ Stream.filter less_complex
     (Stream.of_list [0.;1.;floor f]
      (* Shrink magnitude.  *)
      @@ Stream.guard (f = floor f) (Stream.map float_of_int
                                       (shrink_int (int_of_float f)))
      (* Shrink number of significant figures.  *)
      @@ Stream.map (fun pos -> floor (f *. pos) /. pos)
         (Stream.take_while ((>=) f) (Stream.iterate (( *.) 2.) 1.)))
let shrink_discrete_float f =
  Stream.map (fun k -> float_of_int k *. !discrete_unit)
    (shrink_int (int_of_float (f /. !discrete_unit)))

let rec shrink_list shrink_elem = function
  | [] -> Stream.of_list []
  | x::xs ->
    Stream.singleton xs
    @@ Stream.map (fun xs -> x::xs) (shrink_list shrink_elem xs)
    @@ Stream.map (fun x -> x::xs) (shrink_elem x)

let shrink_fixed_size_list shrink_elem = function
  | [] -> Stream.of_list []
  | x::xs ->
    Stream.map (fun xs -> x::xs) (shrink_list shrink_elem xs)
    @@ Stream.map (fun x -> x::xs) (shrink_elem x)

let shrink_array shrink_elem a =
  Stream.map Array.of_list (shrink_list shrink_elem (Array.to_list a))

let shrink_bigarray1 ?(shrink_size=true) shrink_elem a =
  let open Bigarray in
  let n = Array1.dim a in
  let create = Array1.create (Array1.kind a) (Array1.layout a) in
  let copy () =
    let b = create (n - 1) in
    for j = 0 to n-1 do
      b.{j} <- a.{j}
    done;
    b
  and drop i =
    let b = create (n - 1) in
    for j = 0 to i-1 do
      b.{j} <- a.{j}
    done;
    for j = i+1 to n-1 do
      b.{j} <- a.{j}
    done;
    b
  in
  let shrink_at i =
    Stream.map (fun ai -> let a = copy () in
                          a.{i} <- ai; a)
      (shrink_elem a.{i})
  in
  Stream.guard shrink_size (Stream.map drop
                              (Stream.enum 0 (Bigarray.Array1.dim a)))
  @@ Stream.concat (Stream.map shrink_at (Stream.enum 0 (n-1)))

(* Helper functions for printing data.  *)
let show_sequence show_elem xs = String.concat "; " (List.map show_elem xs)
let show_list show_elem xs =
  "[ " ^ show_sequence show_elem xs ^ " ]"
let show_array show_elem xs =
  "[| " ^ show_sequence show_elem (Array.to_list xs) ^ " |]"
let list_of_bigarray1 xs =
  let a = ref [] in
  for i = Bigarray.Array1.dim xs - 1 downto 0 do
    a := xs.{i}::!a
  done;
  !a
let show_bigarray1 show_elem xs =
  "[|| " ^ show_sequence show_elem (list_of_bigarray1 xs) ^ " ||]"
let dump_bigarray1 dump_elem xs =
  "(Carray.of_array [|"
  ^ show_sequence dump_elem (list_of_bigarray1 xs) ^ "|])"

(* Types and functions for modeling IDA.  *)
module IdaModel =
struct
  module Carray = Sundials.Carray
  type cmd = SolveNormal of float
  type result = Unit | Int of int | Float of float
                | Any
                | Aggr of result list
                | Carray of Carray.t    (* NB: always copy the array! *)
                | SolverResult of Ida.solver_result
                | Exn of exn
  type resfn_type = ResFnLinear of Carray.t

  let carray x = Carray (Carray.of_carray x)

  (* Whole-test results *)
  type failure_type = ResultMismatch of int * result * result
                      | TestCodeDied of result list
                      | TestCodeOverrun
  type test_result = OK | Failed of failure_type

  let cmp_eps = ref 1e-5

  let show_carray = show_bigarray1 string_of_float
  let dump_carray = dump_bigarray1 string_of_float

  let show_cmd = function
    | SolveNormal f ->
      if f >= 0. then "SolveNormal " ^ string_of_float f
      else "SolveNormal (" ^ string_of_float f ^ ")"
  let show_cmds cmds = show_list show_cmd cmds

  let show_solver = function
    | Ida.Dense -> "Dense"
    | Ida.Band range -> Printf.sprintf "Band { mupper=%d; mlower=%d }"
                                       range.Ida.mupper range.Ida.mlower
    | Ida.Sptfqmr _ | Ida.Spbcg _ | Ida.Spgmr _
    | Ida.LapackBand _ | Ida.LapackDense _ ->
      raise (Failure "linear solver not implemented")

  let dump_solver solver = "Ida." ^ show_solver solver

  let rec show_result = function
    | Any -> "_"
    | Unit -> "()"
    | Int i -> string_of_int i
    | Float f -> string_of_float f
    | Carray ca -> show_carray ca
    | SolverResult Ida.Continue -> "Continue"
    | SolverResult Ida.RootsFound -> "RootsFound"
    | SolverResult Ida.StopTimeReached -> "StopTimeReached"
    | Aggr rs -> "Aggr [" ^ String.concat "; " (List.map show_result rs) ^ "]"
    | Exn exn -> "exception " ^ Printexc.to_string exn

  let show_resfn_type ?(dump=false) = function
    | ResFnLinear slope when dump -> "ResFnLinear " ^ dump_carray slope
    | ResFnLinear slope -> "ResFnLinear " ^ show_carray slope

  let dump_resfn_type = show_resfn_type ~dump:true
      
  let rec results_equal r1 r2 =
    match r1, r2 with
    | Any, _ -> true
    | _, Any -> true
    | Unit, Unit -> true
    | Int i1, Int i2 -> i1 = i2
    | Float f1, Float f2 -> abs_float (f1 -. f2) < !cmp_eps
    | Aggr l1, Aggr l2 -> result_lists_equal l1 l2
    | Carray v1, Carray v2 -> carrays_equal v1 v2
    | SolverResult r1, SolverResult r2 -> r1 = r2
    | Exn e1, Exn e2 -> exns_equal e1 e2
    | _, _ -> false
  and result_lists_equal r1 r2 =
    match r1, r2 with
    | [], [] -> true
    | x::xs, y::ys -> results_equal x y && result_lists_equal xs ys
    | _, _ -> false
  and carrays_equal v1 v2 =
    let n = Carray.length v1 in
    let rec go i =
      if i < n then abs_float (v1.{i} -. v2.{i}) < !cmp_eps && go (i+1)
      else true
    in
    n = Carray.length v2 && go 0
  and exns_equal e1 e2 =
  (* Compare only the tags *)
    match e1, e2 with
    | Failure _, Failure _ -> true
    | Invalid_argument _, Invalid_argument _ -> true
    | _, _ -> e1 = e2

  let is_exn = function
    | Exn _ -> true
    | _ -> false
  let not_exn x = not (is_exn x)
  let shrink_carray = shrink_bigarray1
end
