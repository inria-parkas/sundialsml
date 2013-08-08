open Pprint

let (@@) = Fstream.append

(* Controls the size of the generated test case.  *)
let size = ref 5

(* Generation & shrinking *)
let gen_nat () = Random.int (max 1 !size)
let gen_pos () = gen_nat () + 1
let gen_int () = let size = max !size 1 in
                 Random.int (size * 2) - size

(* Rearrange as gen -> shrink -> gen -> shrink.  Make gen_* for each type of
   time value.  Algorithm taken from Haskell's quickcheck.  *)
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

(** Randomly choose an element from an array.  The array elements may be
    generators themselves, but in that case don't forget to delay the choice
    until you need it, e.g:
    {[let gen_int_nat () = gen_choice [|gen_int; gen_nat|]]} ()
    not
    {[let gen_int_nat = gen_choice [|gen_int;gen_nat|]]}  *)
let gen_choice choices = choices.(Random.int (Array.length choices))

(** Shrink an element chosen from an array of candidates, which must be
    comparable by [=].  Candidates listed earlier are considered smaller.  *)
let shrink_choice choices c =
  Fstream.take_while (fun c' -> c' <> c)
    (Fstream.map (Array.get choices)
       (Fstream.enum 0 (Array.length choices - 1)))

(** Make a generator and a shrinker for an array of choices for some
    first-order data type.  Note the generator (the [fst] component of
    [gen_shrink_choice choices]) is not [gen_choice choices] but rather
    [fun () -> gen_choice choices].  *)
let gen_shrink_choice choices =
  (fun () -> gen_choice choices), shrink_choice choices

let enum istart iend =
  let rec go acc i =
    if istart <= i then go (i::acc) (i-1)
    else acc
  in go [] iend

let rec repeat_apply f n x =
  if n = 0 then x
  else if n < 0 then invalid_arg "repeat_apply: negative exponent"
  else repeat_apply f (n-1) (f x)


let gen_list g () = List.map (fun _ -> g ()) (enum 1 (gen_nat ()))
let rec shrink_list shrink_elem = function
  | [] -> Fstream.of_list []
  | x::xs ->
    Fstream.cons xs
      (Fstream.map (fun xs -> x::xs) (shrink_list shrink_elem xs)
       @@ Fstream.map (fun x -> x::xs) (shrink_elem x))

(** Generate a list of values that satisfy an invariant which can be checked by
    scanning the list once, in order.

    [gen_1pass_list gen seed] returns {[[y1, y2, ..., yn]]} where
    [(seed1, y1) = gen seed,
     (seed2, y2) = gen seed1,
     (seed3, y3) = gen seed2,
     ...].

    [gen x] should produce a value taking into account some information [x]
    about previously generated elements, and return that value along with [x]
    updated with information about the new value.  [seed] is the initial value
    of [x].

    For example, [gen_1pass_list (fun x -> let x = gen_nat () + x in (x,x)) 0]
    generates non-strictly increasing lists of natural numbers.

 *)
let gen_1pass_list gen seed () =
  (* In haskell notation,
     let (seeds_tl, ys) = unzip $ map gen seeds
         seeds = seed:seeds_tl
     in take (gen_nat ()) ys
   *)
  let rec seeds_tl_and_ys = lazy (Fstream.unzip (Fstream.map gen seeds))
  and seeds = lazy (Fstream.Cons (seed, fst (Lazy.force seeds_tl_and_ys))) in
  Fstream.to_list
    (Fstream.take (gen_nat ()) (snd (Lazy.force seeds_tl_and_ys)))

(** This function fixes up a list given a function to fix up an element.
   [fixup] should be as specified in {!shrink_1pass_list}.  *)
let fixup_list fixup seed =
  let rec go seed acc = function
    | [] -> List.rev acc
    | x::xs -> let (seed, x) = fixup seed x in
               go seed (x::acc) xs
  in go seed []

(** Shrink a list of values while maintaining an invariant that can be checked
    by scanning the list once, in order.

    {[shrink_1pass_list shrink fixup seed xs]} assumes [shrink] and [fixup] are
    purely functional.  It must be the case that [xs] is one of the lists that
    can be produced by [gen_1pass_list gen seed] using the same [seed] and some
    impure function [gen]; the shrunk lists will also be such lists.

    [shrink] is a function that shrinks one element of the list.  [shrink s x]
    should produce some or all of the possible return values of [gen s] whose
    [snd]'s are "smaller" than [x].

    [fixup] is used when an element of the list is shrunk.  Its job is to
    update all subsequent elements and restore the invariant if it is broken.
    [fixup s x] should produce one of the possible return values of [gen s].
    The [snd] of the return value need not be "smaller" than [x], but it should
    be equal to [x] whenever possible (i.e. it should return [x] as-is if it
    already satisfies the invariant).  Unlike [gen], [fixup] can (and probably
    should) be purely functional.

    Currently, [shrink_1pass_list] enumerates lists produced by the following
    procedure: either drop an element of the list or replace it by a smaller
    value produced by [shrink]; then pass [fixup] through the whole list to
    restore any invariants broken by shrinking.

 *)
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

(** Drop an element of an array-like structure.  *)
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

(** Generate an array-like data structure.  *)
let gen_array_like make set gen_elem ?(size=gen_pos ()) () =
  let a = make size (gen_elem ()) in
  for i = 1 to size-1 do
    set a i (gen_elem ())
  done;
  a

(** Shrink an array-like data structure by dropping elements.  Each array will
    be returned with the index of the dropped element.  This function doesn't
    match the usual signature ['a -> 'a Fstream.t] for shrinkers, so it's not
    given a [shrink_] prefix.  *)
let shorten_array_like make length get set a =
  Fstream.map
    (fun i -> (i, array_like_drop_elem make length get set a i))
    (Fstream.enum 0 (length a - 1))

(** Shrink an array-like data structure by shrinking its elements.  The
    returned arrays will have the same length as the input array.  *)
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

(** Like [shrink_array_like ~shrink_size:true], but returns the index of the
    element that was dropped.  For arrays obtained by shrinking an element
    without shortening the array, the index will be [-1].  *)
let shorten_shrink_array_like make length get set shrink_elem a =
  shorten_array_like make length get set a
  @@ Fstream.map (fun x -> (-1, x))
       (shrink_array_like_elem make length get set shrink_elem a)

let gen_array gen_elem =
  gen_array_like Array.make Array.set gen_elem

let shrink_array shrink_elem =
  shrink_array_like Array.make Array.length Array.get Array.set shrink_elem

let shorten_array a =
  shorten_array_like Array.make Array.length Array.get Array.set a

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

(* Remove duplicates without reordering.  *)
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

let gen_pair gen_x gen_y = (gen_x (), gen_y ())

let shrink_pair shrink_x shrink_y (x,y) =
  Fstream.map (fun x -> (x,y)) (shrink_x x)
  @@ Fstream.map (fun y -> (x,y)) (shrink_y y)

let shrink_fst shrink_x (x,y) = Fstream.map (fun x -> (x,y)) (shrink_x x)
let shrink_snd shrink_y (x,y) = Fstream.map (fun y -> (x,y)) (shrink_y y)

let shrink_fixed_size_list shrink_elem = function
  | [] -> Fstream.of_list []
  | x::xs ->
    Fstream.map (fun xs -> x::xs) (shrink_list shrink_elem xs)
    @@ Fstream.map (fun x -> x::xs) (shrink_elem x)


(** A property is a function from some type ['a] to a ['b test_result].  The
    test result [OK] means the property holds for that data, [Falsified foo]
    means the data falsifies the property, where [foo] is a user-defined
    description of why/how it was falsified, and [Failed] means the property
    raised an exception.  Properties should not return [Failed] but rather let
    exceptions propagate; the quickcheck framework catches them and conerts
    them to [Failed].
 *)
type ('a,'b) property = 'a -> 'b test_result
and  'reason test_result = OK | Falsified of 'reason | Failed of exn

let isOK = function
  | OK -> true
  | _ -> false

(** Convert a function of type ['a -> bool] to a property.  *)
let boolean_prop prop x =
  if prop x then OK
  else Falsified ()

(** A formatter that performs no output.  Useful for disabling output from
    {!quickcheck}, {!minimize}, etc.  *)
let null_formatter = Format.make_formatter (fun _ _ _ -> ()) (fun _ -> ())

(** A shrinker that always fails to produce a shrunk value.  Used as a stub
    when you can't supply a shrinking function.  *)
let no_shrink _ = Fstream.nil

let test_in_sandbox prop x =
  try prop x
  with exn -> Failed exn

(** [minimize shrink prop x reason] minimizes a counterexample [x] of property
    [prop].  [reason] is the return value of [prop x] and is returned when no
    counterexample smaller than [x] is found.  It's up to the caller to ensure
    [reason <> OK]; this function just assumes that's the case.

    If the optional argument [pp_input] is supplied, it is used to print each
    shrunk test case before trying it, along with additional output.  The
    optional argument [pp_formatter] tells where to direct this output.

    Returns (<number of shrinks performed>, <shrunk data>, <prop result>)
 *)
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

(** Returns an input that fails the property with the return value of the
    property, or lack thereof.  If the optional argument [pp_input] is
    specified, dumps each test case before trying it.  [pp_formatter] is where
    this output goes; all other outputs are sent directly to stdout and
    stdin.  *)
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

