(** Helper functions for pretty-printing OCaml data.  *)

(** When set, the functions will try to generate output that can be fed
    directly back into the OCaml toplevel.  If not set, many functions will try
    to make the output concise rather than read-write invariant.  For example,
    floats that happen to have integral values may be printed without a
    dot.  *)
val read_write_invariance : bool ref

(** Locally enable read-write invariance.  This can be nested.  *)
val with_read_write_invariance : (unit -> 'a) -> 'a

(** Naming convention:

     - [pp] pretty-prints to a Format.formatter.  This is the most primitive
            type of pretty-printer.
     - [show] returns a string (it creates and uses its own Buffer.t)

  [pp] and [show] must respect the [read_write_invariance] flag.  [pp] is the
  most "basic" type of pretty-printer in that all higher-order pretty-printers
  like [pp_list] or [show_list] take [pp] as arguments.

     - [dump] is [pp] wrapped in {!with_read_write_invariance}
     - [display] is [show] wrapped in {!with_read_write_invariance}

     - [ppout] prints to [stdout]; this does not go through [std_formatter].
     - [pperr] prints to [stderr]; this does not go through [err_formatter].

  The optional argument ~prec is the ambient precedence, used to determine if
  the output needs to be enclosed in parentheses.  For example, the precedence
  is [Prec.app] when the pretty-printing an object that appears as an argument
  to a contructor or function application.  Common precedences are defined in
  the submodule [Prec], but since OCaml has no infix constructors, it's almost
  always sufficient to use [0] (i.e. immediately inside parens) or [Prec.app]
  (i.e. in operand position).

  Note that [ppout] is not the same as [pp Format.std_formatter].  [ppout]
  sends its output directly to [stdout], so it cannot be mixed with calls
  to [Format.std_formatter].  You have to flush [Format.std_formatter] before
  calling [ppout], just like you have to flush that formatter before calling
  [print_string], [print_int], and the like.

  Note also that [show_*] is different from [string_of_*].  The standard
  library's [string_of_*] functions generally strive for read-write invariance
  and print everything in one line, whereas [show_*] tries to be more succinct
  and does not shy away from inserting line breaks (so it can be poorly laid
  out if it's printed from the middle of a line).

 *)
type 'a pp = ?prec:int -> Format.formatter -> 'a -> unit
and  'a show = ?prec:int -> 'a -> string
and  'a dump = ?prec:int -> Format.formatter -> 'a -> unit
and  'a display = ?prec:int -> 'a -> string
and  'a ppout = ?prec:int -> 'a -> unit
and  'a pperr = ?prec:int -> 'a -> unit

module Prec :
sig
  val app : int
  val plus : int
  val minus : int
  val times : int
  val div : int
end

val show_of_pp : 'a pp -> 'a show
val display_of_pp : 'a pp -> 'a display
val dump_of_pp : 'a pp -> 'a dump
val ppout_of_pp : 'a pp -> 'a ppout
val pperr_of_pp : 'a pp -> 'a pperr

val pp_of_show : 'a show -> 'a pp
val dump_of_show : 'a show -> 'a dump
val display_of_show : 'a show -> 'a display
val ppout_of_show : 'a show -> 'a ppout
val pperr_of_show : 'a show -> 'a pperr


(** Generate all the pretty-printer interfaces given an implementation of
    show.  *)
val printers_of_show :
  'a show -> 'a pp * 'a dump * 'a show * 'a display * 'a ppout * 'a pperr

(** Generate all the pretty-printer interfaces given an implementation of
    pp.  *)
val printers_of_pp :
  'a pp -> 'a pp * 'a dump * 'a show * 'a display * 'a ppout * 'a pperr

(*
  Removed because it doesn't work due to value restriction.

  (** Like printers_of_pp, but works on a function of type ['a pp -> 'a pp].  *)
  val ho_printers_of_pp :
  ('a pp -> 'b pp) ->
  ('a pp -> 'b pp)
  * ('a pp -> 'b dump)
  * ('a pp -> 'b show)
  * ('a pp -> 'b display)
  * ('a pp -> 'b ppout)
  * ('a pp -> 'b pperr)

  (** Like printers_of_show, but works on a function of type
  ['a show -> 'a show].  Note that due to technical limitations [ppout] and
  [pperr] functions take [pp] as inputs, not [ppout] or [pperr]. *)
  val ho_printers_of_show :
  ('a show -> 'b show) ->
  ('a pp -> 'b pp)
  * ('a pp -> 'b dump)
  * ('a pp -> 'b show)
  * ('a pp -> 'b display)
  * ('a pp -> 'b ppout)
  * ('a pp -> 'b pperr)
 *)

val show_int : int show
val dump_int : int dump
val pp_int : int pp
val display_int : int display
val ppout_int : int ppout
val pperr_int : int pperr

val show_unit : unit show
val dump_unit : unit dump
val pp_unit : unit pp
val display_unit : unit display
val ppout_unit : unit ppout
val pperr_unit : unit pperr

val show_bool : bool show
val dump_bool : bool dump
val pp_bool : bool pp
val display_bool : bool display
val ppout_bool : bool ppout
val pperr_bool : bool pperr

val show_char : char show
val dump_char : char dump
val pp_char : char pp
val display_char : char display
val ppout_char : char ppout
val pperr_char : char pperr

val show_string : string show
val dump_string : string dump
val pp_string : string pp
val display_string : string display
val ppout_string : string ppout
val pperr_string : string pperr

(** Returns the string as-is, without adding quotations or escapes.  *)
val show_string_noquote : string show
(** Prints the string as-is, without adding quotations or escapes.  *)
val dump_string_noquote : string dump
(** Prints the string as-is, without adding quotations or escapes.  *)
val pp_string_noquote : string pp
(** Returns the string as-is, without adding quotations or escapes.  *)
val display_string_noquote : string display
(** Prints the string as-is, without adding quotations or escapes.  *)
val ppout_string_noquote : string ppout
(** Prints the string as-is, without adding quotations or escapes.  *)
val pperr_string_noquote : string pperr

(** Returns a string containing the char as-is, without adding quotations or
    escapes.  *)
val show_char_noquote : char show
(** Prints the char as-is, without adding quotations or escapes.  *)
val dump_char_noquote : char dump
(** Prints the char as-is, without adding quotations or escapes.  *)
val pp_char_noquote : char pp
(** Returns a string containing the char as-is, without adding quotations or
    escapes.  *)
val display_char_noquote : char display
(** Prints the char as-is, without adding quotations or escapes.  *)
val ppout_char_noquote : char ppout
(** Prints the char as-is, without adding quotations or escapes.  *)
val pperr_char_noquote : char pperr

val show_float : float show
val dump_float : float dump
val pp_float : float pp
val display_float : float display
val ppout_float : float ppout
val pperr_float : float pperr

(** Make a pretty-printer that prints a fixed string, regardless of its
    input.  [~paren_prec], if set, denotes the minimum precedence at which
    parens are output around the fixed string.  *)
val pp_const : ?paren_prec:int -> string -> 'a pp
val show_const : ?paren_prec:int -> string -> 'a show
val dump_const : ?paren_prec:int -> string -> 'a dump
val display_const : ?paren_prec:int -> string -> 'a display
val ppout_const : ?paren_prec:int -> string -> 'a ppout
val pperr_const : ?paren_prec:int -> string -> 'a pperr

(** Just returns "<fun>".  *)
val show_fun : ('a -> 'b) show
(** Just returns "<fun>".  *)
val dump_fun : ('a -> 'b) dump
(** Just prints "<fun>".  *)
val pp_fun : ('a -> 'b) pp
(** Just prints "<fun>".  *)
val display_fun : ('a -> 'b) display
(** Just prints "<fun>".  *)
val ppout_fun : ('a -> 'b) ppout
(** Just prints "<fun>".  *)
val pperr_fun : ('a -> 'b) pperr

(** [pp_enclose left right cond fmt pp] outputs a pair of delimiters [left] and
    [right], respectively before and after [pp], if [cond] is [true].
    Otherwise just runs [pp].  *)
val pp_enclose :
  (unit, Format.formatter, unit) format
  -> (unit, Format.formatter, unit) format
  -> bool
  -> (Format.formatter -> 'a -> unit)
  -> Format.formatter
  -> 'a
  -> unit

(** [pp_parens cond fmt pp] outputs a pair of parentheses around [pp], provided
    [cond] is true.  Otherwise just runs [pp].  *)
val pp_parens :
  bool
  -> (Format.formatter -> 'a -> unit)
  -> Format.formatter
  -> 'a
  -> unit

(** [pp_brackets cond fmt pp] outputs a pair of sqaure-brackets around [pp],
    provided [cond] is true.  Otherwise just runs [pp].  *)
val pp_brackets :
  bool
  -> (Format.formatter -> 'a -> unit)
  -> Format.formatter
  -> 'a
  -> unit

(** [pp_braces cond fmt pp] outputs a pair of braces around [pp], provided
    [cond] is true.  Otherwise just runs [pp].  *)
val pp_braces :
  bool
  -> (Format.formatter -> 'a -> unit)
  -> Format.formatter
  -> 'a
  -> unit

(** [pp_array_brackets cond fmt pp] outputs a pair of array brackets [| and |]
    around [pp], provided [cond] is true.  Otherwise just runs [pp].  *)
val pp_array_brackets :
  bool
  -> (Format.formatter -> 'a -> unit)
  -> Format.formatter
  -> 'a
  -> unit

(** [pp_double_quotes cond fmt pp] outputs a pair of double-quotes around [pp],
    provided [cond] is true.  Otherwise just runs [pp].  *)
val pp_double_quotes :
  bool
  -> (Format.formatter -> 'a -> unit)
  -> Format.formatter
  -> 'a
  -> unit

(** [pp_quotes cond fmt pp] outputs a pair of single-quotes around [pp],
    provided [cond] is true.  Otherwise just runs [pp].  *)
val pp_quotes :
  bool
  -> (Format.formatter -> 'a -> unit)
  -> Format.formatter
  -> 'a
  -> unit

(** [show_enclose left right cond str] returns [left ^ str ^ right] if [cond]
    is [true], otherwise returns [str].  *)
val show_enclose : string -> string -> bool -> string -> string

(** [show_parens cond str] adds a pair of parentheses around [str] if [cond] is
    [true]; it otherwise returns [str] without change.  *)
val show_parens : bool -> string -> string

(** [show_brackets cond str] adds a pair of square-brackets around [str] if
    [cond] is [true]; it otherwise returns [str] without change.  *)
val show_brackets : bool -> string -> string

(** [show_braces cond str] adds a pair of braces around [str] if [cond] is
    [true]; it otherwise returns [str] without change.  *)
val show_braces : bool -> string -> string

(** [show_double_quotes cond str] adds a pair of double-quotes around [str] if
    [cond] is [true]; it otherwise returns [str] without change.  *)
val show_double_quotes : bool -> string -> string

(** [show_quotes cond str] adds a pair of single-quotes around [str] if [cond]
    is [true]; it otherwise returns [str] without change.  *)
val show_quotes : bool -> string -> string

(** [pp_seq opening delim closing fmt pps] runs a sequence of pretty-printers
    with delimiters.

    - [opening] and [closing] are the opening and closing delimiters to output
      at the beginning and end.  They should open and close a box, usually
      hv box if elements are expected to have irregular lengths or hov box
      if elements are expected to be rather uniform.
    - [delim] is inserted along with a breakable space after every element
      except for the last one.
    - [pps] is a stream of pretty-printing functions (of type
      [(Format.formatter -> unit) Fstream.t]) to run.

 *)
val pp_seq :
  (unit, Format.formatter, unit) format
  -> (unit, Format.formatter, unit) format
  -> (unit, Format.formatter, unit) format
  -> Format.formatter
  -> (Format.formatter -> unit) Fstream.t
  -> unit

(** [pp_fixarg pp x] is just convenient shorthand for [fun fmt -> pp fmt x].
    Used with {!pp_seq}.  *)
val pp_fixarg : 'a pp -> 'a -> (Format.formatter -> unit)

(** Uesd to indicate whether a certain part of the output is optional or
    required.  Optional parts are generally omitted unless
    [read_write_invariance] is set.  *)
type 'a pp_needed = Required of 'a | Optional of 'a

(** [pp_record fields] generates a [pp] implementation for a record type.  The
    argument [fields] is a list of the form ["foo", pp_foo] where [foo] is the
    name of a field of the record and [pp_foo] extracts that field and
    pretty-prints it.  *)
val pp_record : (string * 'a pp) pp_needed array -> 'a pp

(** Pretty-prints a list in a hov-box, i.e. output as many elements as fits on
    each line, just as you would fill a line with words in ordinary English
    writing.  *)
val pp_list : 'a pp -> 'a list pp
val dump_list : 'a pp -> 'a list dump
val show_list : 'a pp -> 'a list show
val display_list : 'a pp -> 'a list display
val ppout_list : 'a pp -> 'a list ppout
val pperr_list : 'a pp -> 'a list pperr

(** Like {!pp_list}, but the element pretty-printer gets the index.  *)
val pp_list_ix : (int -> 'a pp) -> 'a list pp
val dump_list_ix : (int -> 'a pp) -> 'a list dump
val show_list_ix : (int -> 'a pp) -> 'a list show
val display_list_ix : (int -> 'a pp) -> 'a list display
val ppout_list_ix : (int -> 'a pp) -> 'a list ppout
val pperr_list_ix : (int -> 'a pp) -> 'a list pperr

(** Like {!pp_list_ix}, but uses an hv-box, i.e. either print everything in one
    line if that's possible, or else print each element in a separate line.  *)
val pp_hvlist_ix : (int -> 'a pp) -> 'a list pp
val dump_hvlist_ix : (int -> 'a pp) -> 'a list dump
val show_hvlist_ix : (int -> 'a pp) -> 'a list show
val display_hvlist_ix : (int -> 'a pp) -> 'a list display
val ppout_hvlist_ix : (int -> 'a pp) -> 'a list ppout
val pperr_hvlist_ix : (int -> 'a pp) -> 'a list pperr

(** Like {!pp_list}, but uses an hv-box, i.e. either print everything in one
    line if that's possible, or else print each element in a separate line.  *)
val pp_hvlist : 'a pp -> 'a list pp
val dump_hvlist : 'a pp -> 'a list dump
val show_hvlist : 'a pp -> 'a list show
val display_hvlist : 'a pp -> 'a list display
val ppout_hvlist : 'a pp -> 'a list ppout
val pperr_hvlist : 'a pp -> 'a list pperr

(** Like {!pp_list_ix}, but uses a v-box, i.e. exactly one element is listed
    per line.  *)
val pp_vlist_ix : (int -> 'a pp) -> 'a list pp
val dump_vlist_ix : (int -> 'a pp) -> 'a list dump
val show_vlist_ix : (int -> 'a pp) -> 'a list show
val display_vlist_ix : (int -> 'a pp) -> 'a list display
val ppout_vlist_ix : (int -> 'a pp) -> 'a list ppout
val pperr_vlist_ix : (int -> 'a pp) -> 'a list pperr

(** Like {!pp_list}, but uses a v-box, i.e. exactly one element is listed
    per line.  *)
val pp_vlist : 'a pp -> 'a list pp
val dump_vlist : 'a pp -> 'a list dump
val show_vlist : 'a pp -> 'a list show
val display_vlist : 'a pp -> 'a list display
val ppout_vlist : 'a pp -> 'a list ppout
val pperr_vlist : 'a pp -> 'a list pperr

(** {!pp_seq} specialized to array-like data types.  Should be invoked like
    [pp_array_like length get opening closing] where [length], [get] should
    behave like the functions in the module [Array].  *)
val pp_array_like :
  ('a -> int)
  -> ('a -> int -> 'b)
  -> (unit, Format.formatter, unit) format
  -> (unit, Format.formatter, unit) format
  -> 'b pp -> 'a pp

val pp_array : 'a pp -> 'a array pp
val dump_array : 'a pp -> 'a array dump
val show_array : 'a pp -> 'a array show
val display_array : 'a pp -> 'a array display
val ppout_array : 'a pp -> 'a array ppout
val pperr_array : 'a pp -> 'a array pperr

open Bigarray

val pp_bigarray_c_layout : c_layout pp
val pp_bigarray_fortran_layout : fortran_layout pp

(** [pp_bigarray_float32_elt] is a placeholder and always fails when called.
    It is only used as arguments to {!pp_bigarray1} so as to make the interface
    consistent.  The same goes for all other [_elt] functions.  *)
val pp_bigarray_float32_elt : float32_elt pp
val pp_bigarray_float64_elt : float64_elt pp
val pp_bigarray_complex32_elt : complex32_elt pp
val pp_bigarray_complex64_elt : complex64_elt pp
val pp_bigarray_int8_signed_elt : int8_signed_elt pp
val pp_bigarray_int8_unsigned_elt : int8_unsigned_elt pp
val pp_bigarray_int16_signed_elt : int16_signed_elt pp
val pp_bigarray_int16_unsigned_elt : int16_unsigned_elt pp
val pp_bigarray_int_elt : int_elt pp
val pp_bigarray_int32_elt : int32_elt pp
val pp_bigarray_int32_elt : int32_elt pp
val pp_bigarray_int64_elt : int64_elt pp
val pp_bigarray_nativeint_elt : nativeint_elt pp
val pp_bigarray_char_elt : int8_unsigned_elt pp



val pp_bigarray1 :
  'a pp -> 'b pp -> 'c pp -> ('a,'b,'c) Array1.t pp
val dump_bigarray1 :
  'a pp -> 'b pp -> 'c pp -> ('a,'b,'c) Array1.t dump
val show_bigarray1 :
  'a pp -> 'b pp -> 'c pp -> ('a,'b,'c) Array1.t show
val display_bigarray1 :
  'a pp -> 'b pp -> 'c pp -> ('a,'b,'c) Array1.t display
val ppout_bigarray1 :
  'a pp -> 'b pp -> 'c pp -> ('a,'b,'c) Array1.t ppout
val pperr_bigarray1 :
  'a pp -> 'b pp -> 'c pp -> ('a,'b,'c) Array1.t pperr

(** [pp_tuple fields] generates a [pp] implementation for a tuple type.  The
    argument [fields] is a list of functions of type ['a pp] that extracts a
    field of the tuple and pretty-prints it.  *)
val pp_tuple : 'a pp list -> 'a pp
val show_tuple : 'a pp list -> 'a show
val dump_tuple : 'a pp list -> 'a dump
val display_tuple : 'a pp list -> 'a display
val ppout_tuple : 'a pp list -> 'a ppout
val pperr_tuple : 'a pp list -> 'a pperr

val pp_pair : 'a pp -> 'b pp -> ('a * 'b) pp
val show_pair : 'a pp -> 'b pp -> ('a * 'b) show
val dump_pair : 'a pp -> 'b pp -> ('a * 'b) dump
val display_pair : 'a pp -> 'b pp -> ('a * 'b) display
val ppout_pair : 'a pp -> 'b pp -> ('a * 'b) ppout
val pperr_pair : 'a pp -> 'b pp -> ('a * 'b) pperr

val pp_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) pp
val show_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) show
val dump_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) dump
val display_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) display
val ppout_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) ppout
val pperr_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) pperr

(** [pp_precompose f pp] creates a pretty-printer that pre-processes the input
    with the function [f] before passing it to [pp].  For example,
    [pp_precompose fst pp] pretty-prints only the [fst] of its argument.  *)
val pp_precompose : ('a -> 'b) -> 'b pp -> 'a pp
val show_precompose : ('a -> 'b) -> 'b pp -> 'a show
val dump_precompose : ('a -> 'b) -> 'b pp -> 'a dump
val display_precompose : ('a -> 'b) -> 'b pp -> 'a display
val ppout_precompose : ('a -> 'b) -> 'b pp -> 'a ppout
val pperr_precompose : ('a -> 'b) -> 'b pp -> 'a pperr

val pp_fst : 'a pp -> ('a * 'b) pp
val show_fst : 'a pp -> ('a * 'b) show
val dump_fst : 'a pp -> ('a * 'b) dump
val display_fst : 'a pp -> ('a * 'b) display
val ppout_fst : 'a pp -> ('a * 'b) ppout
val pperr_fst : 'a pp -> ('a * 'b) pperr

val pp_snd : 'b pp -> ('a * 'b) pp
val show_snd : 'b pp -> ('a * 'b) show
val dump_snd : 'b pp -> ('a * 'b) dump
val display_snd : 'b pp -> ('a * 'b) display
val ppout_snd : 'b pp -> ('a * 'b) ppout
val pperr_snd : 'b pp -> ('a * 'b) pperr

(** Pretty-print a constructor name and an argument.  For example,
    [pp_ctor (Required "Some") pp_int (-1)] prints "Some (-1)".  Mult-argument
    constructors must be destructured on the caller's side.  For some types
    like [type 'a point = Point of 'a * 'a], it may make sense to do, e.g.
    [match p with
     Point (x,y) -> pp_ctor (Optional "Point") (pp_pair pp_int pp_int) (x,y)]
    with does not output the constructor [Point] unless [read_write_invariance]
    is set.
  *)
val pp_ctor : string pp_needed -> 'a pp -> 'a pp

val pp_option : 'a pp -> 'a option pp
val show_option : 'a pp -> 'a option show
val dump_option : 'a pp -> 'a option dump
val display_option : 'a pp -> 'a option display
val ppout_option : 'a pp -> 'a option ppout
val pperr_option : 'a pp -> 'a option pperr

(** Pretty-prints a lazy value.  If the computation has already been forced,
    this pretty-printer prints the forced value.  Otherwise, this just prints
    "<lazy>" without forcing computation, just like the top level.  *)
val pp_lazy_t : 'a pp -> 'a lazy_t pp
val show_lazy_t : 'a pp -> 'a lazy_t show
val dump_lazy_t : 'a pp -> 'a lazy_t dump
val display_lazy_t : 'a pp -> 'a lazy_t display
val ppout_lazy_t : 'a pp -> 'a lazy_t ppout
val pperr_lazy_t : 'a pp -> 'a lazy_t pperr

(** Alias for {!pp_lazy_t}.  *)
val pp_lazy : 'a pp -> 'a lazy_t pp
val show_lazy : 'a pp -> 'a lazy_t show
val dump_lazy : 'a pp -> 'a lazy_t dump
val display_lazy : 'a pp -> 'a lazy_t display
val ppout_lazy : 'a pp -> 'a lazy_t ppout
val pperr_lazy : 'a pp -> 'a lazy_t pperr

val show_exn : exn show
val dump_exn : exn dump
val pp_exn : exn pp
val display_exn : exn display
val ppout_exn : exn ppout
val pperr_exn : exn pperr
