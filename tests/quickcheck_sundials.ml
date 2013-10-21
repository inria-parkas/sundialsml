(* Generators, shrinkers, and pretty-printers for sundials types that are
   shared between IDA and CVODE.  Pretty-printers derived in pprint_sundials.ml
   are wrapped or overridden to fine-tune the output.  *)
open Pprint
open Quickcheck
open Sundials
include Pprint_sundials

(* Generators and shrinkers for different kinds of time values.  Time values
   can occur in several places in a script:

   - the query time of SolveNormal

   - the time of the zero of a root function

   - the stop time

   We have to make sure that these types of events do not happen simultaneously
   because otherwise the order in which IDA detects them is dictated by
   floating point error and is unpredictable.

   Therefore, we:

   - discretize the time values, i.e. make them a multiple of a fixed value
     called discrete_unit

   - offset each type of time value by a sub-discrete_unit value, so that
     e.g. the stop time is always distinct from a query time modulo
     discrete_unit.

   Be careful, because this means that the set of valid times is not
   necessarily closed under addition!  For instance, a valid stop time + a
   valid stop time is not a valid stop time.  This is not a big deal, with one
   exception: we need to do a fair amount of arithmetic on query times, so we
   specially require query_time_offs to be zero.

 *)
let discrete_unit = 1.
(* NB: lots of code depend on the fact that query_time_offs is zero! *)
let query_time_offs = 0.
let root_time_offs = discrete_unit /. 2.
let stop_time_offs = discrete_unit /. 4.

(* This must be smaller than any of the offsets.  *)
let time_epsilon = discrete_unit /. 8.

type sign = Positive | NonNegative | ArbitrarySign
(* Returns (gen, shrink), where shrink's first argument gives a *non-strict*
   lower bound on the allowable time values.  *)
let discrete_float_type ?(sign=NonNegative) offset =
  let gen_rank, shrink_rank =
    match sign with
    | Positive -> gen_pos, shrink_pos
    | NonNegative -> gen_nat, shrink_nat
    | ArbitrarySign -> gen_int, shrink_int
  in
  ((fun t0 () -> float_of_int (gen_rank ()) *. discrete_unit +. offset +. t0),
   (fun t0 f ->
     Fstream.map
       (fun x -> float_of_int x *. discrete_unit +. offset +. t0)
       (shrink_rank (int_of_float ((f -. t0 -. offset) /. discrete_unit)))))

let gen_t0, shrink_t0 = discrete_float_type ~sign:ArbitrarySign 0.
let gen_t0, shrink_t0 = gen_t0 0., shrink_t0 0.

let gen_query_time, shrink_query_time = discrete_float_type query_time_offs
let gen_root_time, shrink_root_time = discrete_float_type root_time_offs
let gen_stop_time, shrink_stop_time = discrete_float_type stop_time_offs

(* Similar arrangement for non-time values.  *)
let gen_discrete_float, shrink_discrete_float = discrete_float_type 0.


let gen_carray gen_elem =
  gen_bigarray1 Bigarray.float64 Bigarray.c_layout gen_elem

let shrink_carray ?(shrink_size=true) shrink_elem =
  shrink_bigarray1 Bigarray.float64 Bigarray.c_layout ~shrink_size shrink_elem

let carray_drop_elem : Carray.t -> int -> Carray.t = bigarray1_drop_elem

let pp_carray, dump_carray, show_carray, display_carray,
  print_carray, prerr_carray =
  printers_of_pp (fun ?(prec=0) fmt xs ->
    if !read_write_invariance
    then pp_parens (prec >= Prec.app)
          (pp_array_like Carray.length Bigarray.Array1.get
           "Carray.of_array [|@[<hov>" "@]|]" pp_float) fmt xs
    else pp_array_like Carray.length Bigarray.Array1.get
           "[<@[<hov>" "@]>]" pp_float fmt xs)


let gen_root_event, shrink_root_event =
  gen_shrink_choice [| Roots.NoRoot; Roots.Rising; Roots.Falling |]

let gen_root_info = gen_array_like Roots.make Roots.set gen_root_event

let pp_root_info, dump_root_info, show_root_info, display_root_info,
  print_root_info, prerr_root_info
    =
  printers_of_pp (fun ?(prec=0) fmt xs ->
    let get a i =
    (* If Roots is used improperly or if there's a bug in the binding, a root
       info array can contain garbage that doesn't correspond to any of NoRoot,
       Rising, or Falling.  Such values trigger a Failure in Roots.get.  *)
      try show_root_event (Roots.get a i)
      with Failure _ -> "<garbage>"
    in
    if !read_write_invariance
    then pp_parens (prec >= Prec.app)
           (pp_array_like Roots.length get "Roots.of_array [|@[<hov>" "@]|]"
              pp_string_noquote)
           fmt xs
    else pp_array_like Roots.length get "[<" ">]" pp_string_noquote fmt xs)

let gen_root_direction, shrink_root_direction =
  gen_shrink_choice [| RootDirs.Increasing; RootDirs.Decreasing;
                       RootDirs.IncreasingOrDecreasing |]

let gen_root_dirs =
  gen_array_like RootDirs.make RootDirs.set gen_root_direction

let shrink_root_dirs =
  shrink_array_like RootDirs.make RootDirs.length RootDirs.get RootDirs.set
    shrink_root_direction

let pp_root_dirs, dump_root_dirs, show_root_dirs, display_root_dirs
  , print_root_dirs, prerr_root_dirs =
  printers_of_pp
    (fun ?(prec=0) fmt xs ->
      let get a i =
        (* If RootDirs is used improperly or if there's a bug in the binding, a
           root info array can contain garbage that doesn't correspond to any
           of Increasing, Decreasing, or IncreasingOrDecreasing.  Such values
           trigger a Failure in RootDirs.get.  *)
        try RootDirs.string_of_root_direction (RootDirs.get a i)
        with Failure _ -> "<garbage>"
      in
      if !read_write_invariance
      then pp_parens (prec >= Prec.app)
             (pp_array_like RootDirs.length get "RootDirs.of_array [|" "|]"
                pp_string_noquote)
             fmt xs
      else pp_array_like RootDirs.length get "[<" ">]"
             pp_string_noquote fmt xs)
