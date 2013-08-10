open Pprint

let (@@) = Fstream.append

let size = ref 5

let with_size s f =
  let old_size = !size in
  try size := s;
      let ret = f () in
      size := old_size;
      ret
  with exn ->
    size := old_size;
    f ()

let gen_nat () = Random.int (max 1 !size)
let gen_pos () = gen_nat () + 1
let gen_neg () = - gen_pos ()
let gen_int () = let size = max !size 1 in
                 Random.int (size * 2) - size

let shrink_int n =
  (* abs k < abs n, but taking care of overflow.  *)
  let less_complex k =
    match k >= 0, n >= 0 with
    | true, true   -> k < n
    | false, false -> k > n
    | false, true  -> k+n > 0
    | true, false  -> k+n < 0
  in
  Fstream.guard1 (n < -n) (- n)
  @@ Fstream.filter less_complex
     (Fstream.of_list [0;1;2]
      @@ Fstream.filter (fun x -> x <> 0 && x <> 1 && x <> 2)
         (Fstream.take_while less_complex
            (Fstream.map (fun higher_bits -> n - higher_bits)
               (Fstream.iterate (fun x -> x / 2) n))))

let shrink_nat n = Fstream.map abs (shrink_int n)
let shrink_pos n = Fstream.map ((+) 1) (shrink_nat (n-1))
let shrink_neg n = Fstream.map (fun x -> -x) (shrink_pos (-n))

let gen_choice choices = choices.(Random.int (Array.length choices))

let shrink_choice choices c =
  Fstream.take_while (fun c' -> c' <> c)
    (Fstream.map (Array.get choices)
       (Fstream.enum 0 (Array.length choices - 1)))


let gen_shrink_choice choices =
  (fun () -> gen_choice choices), shrink_choice choices

(* Find the least index i s.t. x < a.(i) in a non-strictly increasing array, a.
   Assumes that an index satisfying this constraint exists.  Note ties are
   broken by taking the minimum.  *)
let search_gt a x =
  let rec go l u =
    (* Invariant: l < answer <= u *)
    if u-l <= 1 then u
    else let m = (l+u)/2 in
         if x < a.(m) then go l m
         else go m u
  in go (-1) (Array.length a - 1)

(* Fills dest with the prefix-sum of src, and returns the total.  src and dest
   can be physically equal.  The elements are checked to be non-negative.  *)
let prefix_sum_weights ?context weight_of src dest =
  let invalid_arg =
    match context with
    | None -> (fun s -> invalid_arg ("prefix_sum_nat: " ^ s))
    | Some f -> (fun s -> invalid_arg (f ^ ": " ^ s))
  and s = ref 0 in
  for i = 0 to Array.length src - 1 do
    let w = weight_of src.(i) in
    if w < 0 then invalid_arg "negative weight";
    s := !s + w;
    (* Not a perfect way to detect overflow, but better than doing nothing.  *)
    if !s < 0 then invalid_arg "weights overflow";
    dest.(i) <- !s
  done;
  !s

(* Like prefix_sum_weights but doesn't check for negative weights or overflow,
   and treats src.(i) as 0 if valid.(i) x returns false.  *)
let unsafe_cond_prefix_sum_weights ?context weight_of valid x src dest =
  let s = ref 0 in
  for i = 0 to Array.length src - 1 do
    let src_i = src.(i) in
    if valid src_i x then s := !s + weight_of src_i;
    dest.(i) <- !s
  done;
  !s

(* Like unsafe_cond_prefix_sum but with an implicit src that is filled with 1
   everywhere.  *)
let cond_collect_valid ?context valid x src dest =
  let num_valid = ref 0 in
  for i = 0 to Array.length dest - 1 do
    if valid src.(i) x then
      (dest.(!num_valid) <- i;
       num_valid := !num_valid + 1)
  done;
  !num_valid

let ceil_lg x =
  let rec go ret = function
    | 0 -> ret
    | n -> go (ret + 1) (n / 2)
  in go 0 x

let gen_weighted_choice choices =
  let cdf = Array.make (Array.length choices) 0 in
  let total =
    prefix_sum_weights ~context:"gen_weighted_choice" fst choices cdf
  in
  fun () -> snd (choices.(search_gt cdf (Random.int total)))

let gen_cond_choice ?show_input choices =
  let n = Array.length choices in
  let perm_index = Array.make (Array.length choices) 0 in
  let context = "gen_cond_choice" in
  let no_choice x =
    let msg = context ^ ": all choices masked out"; in
    match show_input with
    | None -> invalid_arg msg
    | Some show -> invalid_arg (msg ^ " for " ^ show x)
  in
  let slow_path x =
    let num_valid = cond_collect_valid ~context fst x choices perm_index in
    if num_valid = 0 then no_choice x;
    snd choices.(perm_index.(Random.int num_valid))
  in
  (* The fast path assumes the number of invalidated entries are small.
     Provided the majority (> 1/2) of elements are valid, we only need to fall
     back on the slow path with probability at most 1/1024 for n < 1024, and
     1/n^2 for n >= 1024.  *)
  let rec fast_path tries x =
    if tries = 0 then slow_path x
    else
      let i = Random.int n in
      if fst choices.(i) x then snd choices.(i)
      else fast_path (tries - 1) x
  in
  fast_path (max 10 (2 * ceil_lg n))

let gen_weighted_cond_choice ?show_input choices =
  let context = "gen_weighted_cond_choice" in
  let weight_of (_,w,_) = w
  and predicate_of (p,_,_) = p
  and candidate_of (_,_,c) = c
  in
  let n = Array.length choices in
  let cdf = Array.make (Array.length choices) 0 in
  let masked_cdf = Array.copy cdf in
  let total = prefix_sum_weights ~context weight_of choices cdf in
  let no_choice x =
    let msg = context ^ ": all choices masked out"; in
    match show_input with
    | None -> invalid_arg msg
    | Some show -> invalid_arg (msg ^ " for " ^ show x)
  in
  let slow_path x =
    let total =
      unsafe_cond_prefix_sum_weights ~context
        weight_of predicate_of x choices masked_cdf
    in
    if total = 0 then no_choice x;
    candidate_of choices.(search_gt masked_cdf (Random.int total))
  in
  (* The fast path assumes the number of invalidated entries are small.
     Provided the majority (> 1/2) of elements are valid, we only need to fall
     back on the slow path with probability at most 1/1024 for n < 1024, and
     1/n^2 for n >= 1024.  *)
  let rec fast_path tries x =
    if tries = 0 then slow_path x
    else
      let (valid, _, cand) = choices.(search_gt cdf (Random.int total)) in
      if valid x then cand
      else fast_path (tries - 1) x
  in
  fast_path (max 10 (2 * ceil_lg n))

let enum istart iend =
  let rec go acc i =
    if istart <= i then go (i::acc) (i-1)
    else acc
  in go [] iend

let gen_list g () = List.map (fun _ -> g ()) (enum 1 (gen_nat ()))

let rec shrink_list shrink_elem ?(shrink_size=true) = function
  | [] -> Fstream.of_list []
  | x::xs ->
    Fstream.guard1 shrink_size xs
    @@ Fstream.map (fun xs -> x::xs) (shrink_list shrink_elem xs)
    @@ Fstream.map (fun x -> x::xs) (shrink_elem x)

let gen_1pass_list gen seed ?(size=gen_nat ()) () =
  (* In haskell notation,
     let (seeds_tl, ys) = unzip $ map gen seeds
         seeds = seed:seeds_tl
     in take size ys
   *)
  let rec seeds_tl_and_ys = lazy (Fstream.unzip (Fstream.map gen seeds))
  and seeds = lazy (Fstream.Cons (seed, fst (Lazy.force seeds_tl_and_ys))) in
  Fstream.to_list
    (Fstream.take size (snd (Lazy.force seeds_tl_and_ys)))

let fixup_list fixup seed =
  let rec go seed acc = function
    | [] -> List.rev acc
    | x::xs -> let (seed, x) = fixup seed x in
               go seed (x::acc) xs
  in go seed []


let shrink_1pass_list shrink fixup seed xs =
  let rec go seed = function
    | [] -> Fstream.nil
    | x::xs ->
      Fstream.cons (fixup_list fixup seed xs)  (* drop x *)
        (Fstream.map                           (* keep x *)
           (fun xs -> x::xs)
           (go (fst (fixup seed x)) xs)
         @@ Fstream.map                        (* shrink x *)
             (fun (seed, x) -> x::fixup_list fixup seed xs)
             (shrink seed x))
  in
  go seed xs

let array_like_drop_elem make length get set a i =
  let n = length a in
  assert (i < n);
  let a' = make (n-1) (get a 0) in
  for j = 1 to i-1 do
    set a' j (get a j)
  done;
  for j = i+1 to n-1 do
    set a' (j-1) (get a j)
  done;
  a'

let gen_array_like make set gen_elem ?(size=gen_pos ()) () =
  let a = make size (gen_elem ()) in
  for i = 1 to size-1 do
    set a i (gen_elem ())
  done;
  a

let shorten_array_like make length get set a =
  Fstream.map
    (fun i -> (i, array_like_drop_elem make length get set a i))
    (Fstream.enum 0 (length a - 1))

let shrink_array_like_elem make length get set shrink_elem a =
  let copy a =
    let n = length a in
    if n = 0 then a
    else
      let a' = make n (get a 0) in
      for i = 1 to n-1 do
        set a' i (get a i)
      done;
      a'
  in
  let shrink_one a i =
    Fstream.map (fun x -> let a = copy a in set a i x; a)
      (shrink_elem (get a i))
  in
  Fstream.concat (Fstream.map (shrink_one a) (Fstream.enum 0 (length a - 1)))

let shrink_array_like make length get set shrink_elem ?(shrink_size=true) a =
  Fstream.guard shrink_size
    (Fstream.map snd (shorten_array_like make length get set a))
  @@ shrink_array_like_elem make length get set shrink_elem a

let shorten_shrink_array_like make length get set shrink_elem a =
  shorten_array_like make length get set a
  @@ Fstream.map (fun x -> (-1, x))
       (shrink_array_like_elem make length get set shrink_elem a)

let gen_array gen_elem =
  gen_array_like Array.make Array.set gen_elem

let shrink_array shrink_elem =
  shrink_array_like Array.make Array.length Array.get Array.set shrink_elem

let shorten_array shrink_elem =
  shorten_array_like Array.make Array.length Array.get Array.set shrink_elem

let shorten_shrink_array a =
  shorten_shrink_array_like Array.make Array.length Array.get Array.set a

let array_drop_elem a i =
  array_like_drop_elem Array.make Array.length Array.get Array.set a i

let gen_bigarray1 kind layout gen_elem =
  let make n x =
    let a = Bigarray.Array1.create kind layout n in
    Bigarray.Array1.fill a x;
    a
  in
  gen_array_like make Bigarray.Array1.set gen_elem

let bigarray1_drop_elem a i =
  let make n x =
    let a =
      Bigarray.Array1.create
        (Bigarray.Array1.kind a)
        (Bigarray.Array1.layout a)
        n
    in
    Bigarray.Array1.fill a x;
    a
  in
  array_like_drop_elem make Bigarray.Array1.dim
    Bigarray.Array1.get Bigarray.Array1.set a i

let shrink_bigarray1 kind layout shrink_elem =
  let make n x =
    let a = Bigarray.Array1.create kind layout n in
    Bigarray.Array1.fill a x;
    a
  in
  shrink_array_like
    make
    Bigarray.Array1.dim
    Bigarray.Array1.get
    Bigarray.Array1.set
    shrink_elem

let uniq_list ls =
  let seen = Hashtbl.create 10 in
  let rec go acc = function
    | [] -> List.rev acc
    | x::xs when Hashtbl.mem seen x -> go acc xs
    | x::xs -> Hashtbl.add seen x (); go (x::acc) xs
  in go [] ls

let uniq_array a =
  let seen = Hashtbl.create 10 in
  let b = Array.copy a in
  let bsize = ref 0 in
  for i = 0 to Array.length a - 1 do
    if not (Hashtbl.mem seen a.(i)) then
      (b.(!bsize) <- a.(i);
       bsize := !bsize + 1;
       Hashtbl.add seen a.(i) ())
  done;
  Array.sub b 0 !bsize

let gen_pair gen_x gen_y () = (gen_x (), gen_y ())

let shrink_pair shrink_x shrink_y (x,y) =
  Fstream.map (fun x -> (x,y)) (shrink_x x)
  @@ Fstream.map (fun y -> (x,y)) (shrink_y y)

let shrink_fst shrink_x (x,y) = Fstream.map (fun x -> (x,y)) (shrink_x x)
let shrink_snd shrink_y (x,y) = Fstream.map (fun y -> (x,y)) (shrink_y y)

let no_shrink _ = Fstream.nil


type ('a,'b) property = 'a -> 'b test_result
and  'reason test_result = OK | Falsified of 'reason | Failed of exn

let isOK = function
  | OK -> true
  | _ -> false

let boolean_prop prop x =
  if prop x then OK
  else Falsified ()

let null_formatter = Format.make_formatter (fun _ _ _ -> ()) (fun _ -> ())


let test_in_sandbox prop x =
  try prop x
  with exn -> Failed exn

let minimize ?pp_input ?(pp_formatter=Format.err_formatter) shrink prop x res =
  let trace, pp_input =
    match pp_input with
    | Some s -> true, s
    | None -> false, (fun _ -> failwith "internal error")
  in
  if trace then Format.pp_print_char pp_formatter '\n';
  let rec go ct x reason =
    let failure x =
      let res = test_in_sandbox prop x in
      if trace then
        (Format.fprintf pp_formatter "@[<2>Trying:@\n";
         pp_input pp_formatter x;
         Format.fprintf pp_formatter "@]@\n -> %s@."
           (if isOK res then "triggers bug"
            else "not a counterexample"));
      if res = OK then None
      else Some (x, res)
    in
    match Fstream.find_some failure (shrink x) with
    | None -> (ct, x, reason)
    | Some (x, reason) -> go (ct+1) x reason
  in go 0 x res

let test_case_number = ref 0

let quickcheck gen shrink ?pp_input ?(pp_formatter=Format.err_formatter)
    prop max_tests =
  let old_size = !size in
  let minimize x res =
    if shrink == no_shrink then (x, res)
    else
      begin
        Printf.fprintf stderr "Shrinking...";
        flush stderr;
        let (ct, x, res) =
          minimize ?pp_input ~pp_formatter shrink prop x res
        in
        Printf.fprintf stderr "%d shrinks.\n" ct;
        flush stderr;
        (x, res)
      end
  in
  let trace, pp_input =
    match pp_input with
    | Some pp_input -> true, pp_input
    | None -> false, (fun _ -> failwith "internal error")
  in
  let rec test num_passed =
    let gen () =
      try size := num_passed; gen ()
      with exc ->
        Printf.fprintf stderr
          "Error: the generator failed (after %d tests)\n%s\n"
          num_passed (Printexc.to_string exc);
        size := old_size;
        raise exc
    in
    let check_for_bug x =
      if trace then
        (Format.fprintf pp_formatter "@[<2>Test Case %d:@\n" (num_passed+1);
         pp_input pp_formatter x;
         Format.fprintf pp_formatter "@.");
      let res = test_in_sandbox prop x in
      match res with
      | OK -> res
      | Falsified _ ->
        Printf.fprintf stderr "Failed! (after %d test(s))\n" num_passed;
        flush stderr;
        res
      | Failed exc -> 
        Printf.fprintf stderr
          "Failed! Exception raised (after %d test(s)):\n%s\n"
          num_passed
          (Printexc.to_string exc);
        flush stderr;
        res
    in
    test_case_number := num_passed + 1;
    if num_passed < max_tests then
      let x = gen () in
      match check_for_bug x with
      | OK -> if (not trace) then print_char '*';
              flush stdout;
              test (num_passed + 1)
      | res -> Some (minimize x res)
    else
      (Printf.printf "\n+++ OK, passed %d tests." max_tests;
       None)
  in
  let ret = test 0 in
  size := old_size;
  ret

