open Pprint

(** Alias of [Fstream.append].  *)
val (@@) : 'a Fstream.t -> 'a Fstream.t -> 'a Fstream.t

(** Controls the size of the generated test case.  *)
val size : int ref

(** [with_size s f] temporarily sets the [size] to [s], calls [f], restores
    [size], and returns.  *)
val with_size : int -> (unit -> 'a) -> 'a

(** Generate a non-negative integer.  Suitable for deciding the size of arrays,
    lists, etc. to generate.  *)
val gen_nat : unit -> int

(** Generate a positive integer.  Suitable for deciding the size of arrays,
    lists, etc. to generate.  *)
val gen_pos : unit -> int

(** Generate a negative integer.  The generated integer is usually small
    enough not to trigger overflows, which may or may not be desirable.  *)
val gen_neg : unit -> int

(** Generate an arbitrary integer.  The generated integer is usually small
    enough not to trigger overflows, which may or may not be desirable.  *)
val gen_int : unit -> int

(** Shrink an int.  *)
val shrink_int : int -> int Fstream.t

(** Shrink a non-negative integer to non-negative integers.  *)
val shrink_nat : int -> int Fstream.t

(** Shrink a positive integer to positive integers.  *)
val shrink_pos : int -> int Fstream.t

(** Shrink a negative integer to negative integers.  *)
val shrink_neg : int -> int Fstream.t

(** Randomly choose an element from an array of candidates.

    The candidates may be generators, but in that case the choice should be
    delayed until it's needed, e.g:
    {[let gen_int_nat () = gen_choice [|gen_int; gen_nat|] ()]}
    not
    {[let gen_int_nat = gen_choice [|gen_int;gen_nat|]]}
 *)
val gen_choice : 'a array -> 'a

(** Weighted-probability random selection.  Given an array
    {[[|(w0,c0); (w1,c1); ... |]]}, this function returns the candidate [c0]
    with probability proportional to [w0], the candidate [c1] with probability
    proportional to [w1], and so on.
*)
val gen_weighted_choice : (int * 'a) array -> unit -> 'a

(** Uniform-probability random selection with conditional masking.  This
    function selects one element of an array of candidates, but only one for
    which the condition is right.  Given {[[|(p0,c0); (p1,c1); ... |]]} and
    some [x], this function returns one of the [ci]'s such that the
    corresponding predicate [pi x] returns [true].  The predicates are
    recomputed for each [x].  The [pi] must be pure functions.

    This function has an O(n) setup cost when it receives the array, and
    thereafter for every [x] it runs in O(log n) most of the time, provided a
    majority of the [pi]'s return [true].  More precisely, it tries an O(1)
    opportunistic algorithm O(log n) times, each of which can fail with
    probability k/n where k is the number of [pi]'s that return [false].  If
    all tries fail, it falls back to an O(n) algorithm.  The probability of
    falling back on the O(n) algorithm is (k/n)^(2 lg(n)) which, if k <= n/2,
    is at most 1/n^2.

    It is an error if all of the [pi]'s return [false].  In this case an
    [Invalid_argument] is thrown.  If the optional argument [show_input] is
    specified, it is used to get a string representation of the [x] that
    triggered the error, which is included in the error message.
*)
val gen_cond_choice :
  ?show_input:('a -> string)
  -> (('a -> bool) * 'b) array
  -> 'a
  -> 'b

(** Weighted-probability random selection with conditional masking.  This
    function selects one element of an array of candidates according to their
    weights, but only one for which the condition is right.

    Given {[[|(p0,w0,c0); (p1,w1,c1); ... |]]} and some [x], this function
    returns one of the [ci]'s such that the corresponding predicate [pi x]
    returns [true].  [c0] is selected with probability proportional to [w0]
    provided [p0 x] is [true], [c1] is selected with probability proportional
    to [w1] provided [p1 x] is [true], and so on.  The [pi]'s are recomputed
    for each [x] (although the function tries to avoid actually evaluating the
    recompuations as much as it can).  The [pi] must be pure functions.

    This function has an O(n) setup cost when it receives the array, and
    thereafter for every [x] it runs in O((log n)^2) most of the time, provided
    a majority of the [pi]'s return [true].  More precisely, it tries an
    O(log n) opportunistic algorithm O(log n) times, each of which can fail
    with probability k/n where k is the number of [pi]'s that return [false].
    If all tries fail, it falls back to an O(n) algorithm.  The probability of
    falling back on the O(n) algorithm is (k/n)^(2 lg(n)) which, if k <= n/2,
    is at most 1/n^2.

    It is an error if all of the [pi]'s return [false].  In this case an
    [Invalid_argument] is thrown.  If the optional argument [show_input] is
    specified, it is used to get a string representation of the [x] that
    triggered the error, which is included in the error message.
*)
val gen_weighted_cond_choice :
  ?show_input:('a -> string)
  -> (('a -> bool) * int * 'b) array
  -> 'a
  -> 'b

(** Shrink an element chosen from an array of candidates (a la [gen_choice]),
    which must be comparable by [=].  Candidates listed earlier are considered
    smaller.  *)
val shrink_choice : 'a array -> 'a -> 'a Fstream.t

(** Make a generator and a shrinker for an array of choices for some
    first-order data type.  Note the generator (the [fst] component of
    [gen_shrink_choice choices]) is not [gen_choice choices] but rather
    [fun () -> gen_choice choices].  *)
val gen_shrink_choice : 'a array -> ((unit -> 'a) * ('a -> 'a Fstream.t))

(** Utility function.  [enum a b] returns the list [a; a+1; ...; b].  *)
val enum : int -> int -> int list

(** Make a generator for ['a list] given a generator for ['a].  *)
val gen_list : (unit -> 'a) -> (unit -> 'a list)

(** Make a shrinker for ['a list] given a shrinker for ['a].  Some of the
    shrunk lists may have shorter lengths than the input list iff [shrink_size]
    is [true] (which is the default behavior).  *)
val shrink_list :
  ('a -> 'a Fstream.t)
  -> ?shrink_size:bool
  -> ('a list -> 'a list Fstream.t)

(** Make a generator of lists which satisfy an invariant that can be checked by
    scanning the list once.

    [gen_1pass_list gen seed ()] returns {[[y1, y2, ..., yn]]} where
    [(seed1, y1) = gen seed,
     (seed2, y2) = gen seed1,
     (seed3, y3) = gen seed2,
     ...].

    [gen x] should return [(x',v)] where [v] is some randomly selected value
    taking into account some information [x] about previously generated
    elements, and [x'] is just [x] updated with information about [v].  [seed]
    is the initial value of [x].

    For example, [gen_1pass_list (fun x -> let x = gen_nat () + x in (x,x)) 0]
    generates non-strictly increasing lists of natural numbers.

 *)
val gen_1pass_list : ('a -> 'a * 'b) -> 'a -> ?size:int -> unit -> 'b list

(** Make a shrinker to go with a generator create by [gen_1pass_list].  This
    shrinks a list while maintaining an invariant that can be checked by
    scanning the list once.

    {[shrink_1pass_list shrink fixup seed xs]} assumes [shrink] and [fixup] are
    purely functional.  It must be the case that [xs] is one of the lists that
    can be produced by [gen_1pass_list gen seed] using the same [seed] and some
    impure function [gen]; the shrunk lists will also be such lists.

    [shrink] is a function that shrinks one element of the list.  [shrink s x]
    should produce some or all of the possible return values of [gen s] whose
    [snd]'s are smaller than [x].

    [fixup] is used when an element of the list is shrunk.  Its job is to
    update subsequent elements to restore the invariant if it is broken.
    [fixup s x] should produce one of the possible return values of [gen s].
    The [snd] of the return value need not be "smaller" than [x], but it should
    be equal to [x] whenever possible (i.e. it should return [x] as-is if it
    already satisfies the invariant).  Unlike [gen], [fixup] can (and probably
    should) be purely functional.

    Currently, [shrink_1pass_list] enumerates lists produced by the following
    procedure: either drop an element of the list or replace it by a smaller
    value produced by [shrink]; then pass [fixup] through the whole list to
    restore any invariants broken by the shrinking.

 *)
val shrink_1pass_list :
  ('a -> 'b -> ('a * 'b) Fstream.t)
  -> ('a -> 'b -> 'a * 'b)
  -> 'a
  -> 'b list
  -> 'b list Fstream.t

(** This is just Haskell's [mapAccumL] that extracts the list and throws away
    the final state.

    If [f] is a function such that [f ctx x = (ctx',x')] where [ctx] is some
    contextual information, [x] is an input, [x'] is [x] modified to be
    allowable under that context, and [ctx'] is [ctx] updated with information
    about [x'], then [fixup_list f ctx0 xs] modifies all elements of [xs] by
    applying [f] with the context [ctx0] updated with elements that came
    before.

    This function fixes up a list to satisfy a constraint, in the sense that
    {!shrink_1pass_list} fixes up lists to restore some invariant.  *)
val fixup_list : ('a -> 'b -> 'a * 'b) -> 'a -> 'b list -> 'b list

(** [array_like_drop_elem make length get set a i] creates a new array that
    contains all the elements of [a] except the one at index [i].  So the new
    array will be on element shorter.  [make], [length], [get], and [set]
    should behave like the functions of those names from the [Array] module.

    This is useful for controlling shrinking of arrays.  *)
val array_like_drop_elem :
  (int -> 'a -> 'array)
  -> ('array -> int)
  -> ('array -> int -> 'a)
  -> ('array -> int -> 'a -> unit)
  -> 'array
  -> int
  -> 'array

(** [gen_array_like make set gen_elem ()] generates an array-like data
    structure.  The [make] and [set] should be to the array-like data structure
    what [Array.make] and [Array.set] are to regular arrays.  *)
val gen_array_like :
  (int -> 'a -> 'array)
  -> ('array -> int -> 'a -> unit)
  -> (unit -> 'a)
  -> ?size:int
  -> unit
  -> 'array

(** [shorten_array_like make length get set shrink_elem a] shrinks an
    array-like data structure by dropping elements.  Each array will be
    returned with the index of the dropped element.  This function doesn't
    match the usual signature ['a -> 'a Fstream.t] for shrinkers, so it's not
    given a [shrink_] prefix.  [make], [length], [get], and [set] should be to
    the array-like data structure what functions of the same names in the
    [Array] module are to regular arrays.  *)
val shorten_array_like :
  (int -> 'a -> 'array)
  -> ('array -> int)
  -> ('array -> int -> 'a)
  -> ('array -> int -> 'a -> unit)
  -> 'array
  -> (int * 'array) Fstream.t

(** [shrink_array_like_elem make legnth get set shrink_elem a] shrinks an
    array-like data structure by shrinking its elements.  The returned arrays
    will have the same length as the input array.  [make], [length], [get], and
    [set] should be to the array-like data structure what functions of the same
    names in the [Array] module are to regular arrays.  *)
val shrink_array_like_elem :
  (int -> 'a -> 'array)
  -> ('array -> int)
  -> ('array -> int -> 'a)
  -> ('array -> int -> 'a -> unit)
  -> ('a -> 'a Fstream.t)
  -> 'array
  -> 'array Fstream.t

(** [shrink_array_like make legnth get set shrink_elem a] shrinks an array-like
    data structure.  The length is shrunk if the optional argument
    [shrink_size] is set to [true] (which is the default behavior), but the
    length is kept the same if [shrink_size] is [false].  [make], [length],
    [get], and [set] should be to the array-like data structure what functions
    of the same names in the [Array] module are to regular arrays.  *)
val shrink_array_like :
  (int -> 'a -> 'array)
  -> ('array -> int)
  -> ('array -> int -> 'a)
  -> ('array -> int -> 'a -> unit)
  -> ('a -> 'a Fstream.t)
  -> ?shrink_size:bool
  -> 'array
  -> 'array Fstream.t

(** Like [shrink_array_like ~shrink_size:true], but returns the index of the
    element that was dropped.  For arrays obtained by shrinking an element
    without shortening the array, the index will be [-1].  *)
val shorten_shrink_array_like :
  (int -> 'a -> 'array)
  -> ('array -> int)
  -> ('array -> int -> 'a)
  -> ('array -> int -> 'a -> unit)
  -> ('a -> 'a Fstream.t)
  -> 'array
  -> (int * 'array) Fstream.t

(** {!gen_array_like} instantiated for ordinary arrays.  *)
val gen_array : (unit -> 'a) -> (?size:int -> unit -> 'a array)

(** {!shrink_array_like} instantiated for ordinary arrays.  *)
val shrink_array :
  ('a -> 'a Fstream.t)
  -> (?shrink_size:bool -> 'a array -> 'a array Fstream.t)

(** {!array_like_drop_elem} instantiated to ordinary arrays.  *)
val array_drop_elem : 'a array -> int -> 'a array

(** {!shorten_array_like} instantiated for ordinary arrays.  *)
val shorten_array : 'a array -> (int * 'a array) Fstream.t

(** {!shorten_shrink_array_like} instantiated for ordinary arrays.  *)
val shorten_shrink_array :
  ('a -> 'a Fstream.t)
  -> 'a array
  -> (int * 'a array) Fstream.t

(** {!gen_array_like} instantiated for 1-dimensional bigarrays.  *)
val gen_bigarray1 :
  ('a, 'b) Bigarray.kind
  -> 'c Bigarray.layout
  -> (unit -> 'a)
  -> (?size:int -> unit -> ('a,'b,'c) Bigarray.Array1.t)

(** {!shrink_array_like} instantiated for 1-dimensional bigarrays.  *)
val shrink_bigarray1 :
  ('a, 'b) Bigarray.kind
  -> 'c Bigarray.layout
  -> ('a -> 'a Fstream.t)
  -> (?shrink_size:bool
      -> ('a,'b,'c) Bigarray.Array1.t
      -> ('a,'b,'c) Bigarray.Array1.t Fstream.t)

(** {!array_like_drop_elem} instantiated to 1-dimensional bigarrays.  *)
val bigarray1_drop_elem :
  ('a,'b,'c) Bigarray.Array1.t
  -> int
  -> ('a,'b,'c) Bigarray.Array1.t

(** Remove duplicates from a list without reordering.  *)
val uniq_list : 'a list -> 'a list

(** Remove duplicates from an array without reordering.  *)
val uniq_array : 'a array -> 'a array

(** Make a generator for a pair given a pair of generators.  *)
val gen_pair : (unit -> 'a) -> (unit -> 'b) -> (unit -> 'a * 'b)

(** Make a shrinker for a pair given a pair of shrinkers.  *)
val shrink_pair :
  ('a -> 'a Fstream.t)
  -> ('b -> 'b Fstream.t)
  -> ('a * 'b -> ('a * 'b) Fstream.t)

(** Make a function that shrinks a pair by applying the given shrinker to the
   [fst].  *)
val shrink_fst : ('a -> 'a Fstream.t) -> ('a * 'b -> ('a * 'b) Fstream.t)

(** Make a function that shrinks a pair by applying the given shrinker to the
   [snd].  *)
val shrink_snd : ('b -> 'b Fstream.t) -> ('a * 'b -> ('a * 'b) Fstream.t)


(** A shrinker that always fails to produce a shrunk value.  Used as a stub
    for an unimplemented shrinking function.  *)
val no_shrink : 'a -> 'a Fstream.t

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

val pp_property : 'a pp -> 'b pp -> ('a,'b) property pp
val pp_test_result : 'a pp -> 'a test_result pp

val isOK : 'a test_result -> bool

(** Convert a function of type ['a -> bool] to a property.  *)
val boolean_prop : ('a -> bool) -> ('a, unit) property

(** A formatter that performs no output.  Useful for disabling output from
    {!quickcheck}, {!minimize}, etc.  *)
val null_formatter : Format.formatter

(** [minimize shrink prop x reason] minimizes a counterexample [x] of property
    [prop].  [reason] is the return value of [prop x] and is returned when no
    counterexample smaller than [x] is found.  It's up to the caller to ensure
    [reason <> OK]; this function just assumes that's the case.

    If the optional argument [pp_input] is supplied, it is used to print each
    shrunk test case before trying it, along with additional output.  The
    optional argument [pp_formatter] tells where to direct this output.

    Returns (<number of shrinks performed>, <shrunk data>, <prop result>)
 *)
val minimize :
  ?pp_input : 'a pp
  -> ?pp_formatter : Format.formatter
  -> ('a -> 'a Fstream.t)
  -> ('a,'b) property
  -> 'a
  -> 'b test_result
  -> int * 'a * 'b test_result

(** During a callback from {!quickcheck} or {!minimize}, this variable contains
    the test case number of the current test.  *)
val test_case_number : int ref

(** Returns an input that fails the property along with the corresponding
    [test_result], or lack thereof.  If the optional argument [pp_input] is
    specified, dumps each test case before trying it.  [pp_formatter] is where
    this output goes; all other outputs are sent directly to stdout and stderr.
    All callbacks can refer to the variable [test_case_number] to get a serial
    number for the test case it is called upon.  *)
val quickcheck :
  (unit -> 'a)
  -> ('a -> 'a Fstream.t)
  -> ?pp_input:'a pp
  -> ?pp_formatter:Format.formatter
  -> ('a,'b) property
  -> int
  -> ('a * 'b test_result) option


