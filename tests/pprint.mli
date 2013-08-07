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

     - [pp] pretty-prints to a Format.formatter
     - [show] returns a string (it creates and uses its own Buffer.t)

  [pp] and [show] must respect the [read_write_invariance] flag.  [pp] is the
  most "basic" type of pretty-printer in that all higher-order pretty-printers
  like [pp_list] or [show_list] take [pp] as arguments.

     - [dump] is [pp] wrapped in {!with_read_write_invariance}
     - [display] is [show] wrapped in {!with_read_write_invariance}

     - [print] is [pp] with the formatter fixed to [std_formatter]
     - [prerr] is [pp] with the formatter fixed to [err_formatter]

  Note [print_float], [print_string], and [print_char] are already defined in
  [Pervasives] and their behaviors do not match the description above.  For
  example, [Pervasives.print_string] never produces double-quotation marks
  whereas this module's [print_string] does.  To avoid confusion, the version
  of these function defined in this module are hidden inside a submodule named
  [PrettyPrim] by default.

 *)
type 'a pp = Format.formatter -> 'a -> unit
and  'a show = 'a -> string
and  'a dump = 'a pp
and  'a display = 'a show
and  'a print = 'a -> unit
and  'a prerr = 'a -> unit

val show_of_pp : 'a pp -> 'a show
val display_of_pp : 'a pp -> 'a display
val dump_of_pp : 'a pp -> 'a dump
val print_of_pp : 'a pp -> 'a print
val prerr_of_pp : 'a pp -> 'a prerr

val pp_of_show : 'a show -> 'a pp
val dump_of_show : 'a show -> 'a dump
val display_of_show : 'a show -> 'a display
val print_of_show : 'a show -> 'a print
val prerr_of_show : 'a show -> 'a prerr

(** Generate all the pretty-printer interfaces given an implementation of
    show.  *)
val printers_of_show :
  'a show -> 'a pp * 'a dump * 'a show * 'a display * 'a print * 'a prerr

(** Generate all the pretty-printer interfaces given an implementation of
    pp.  *)
val printers_of_pp :
  'a pp -> 'a pp * 'a dump * 'a show * 'a display * 'a print * 'a prerr

(*
  Removed because it doesn't work due to value restriction.

  (** Like printers_of_pp, but works on a function of type ['a pp -> 'a pp].  *)
  val ho_printers_of_pp :
  ('a pp -> 'b pp) ->
  ('a pp -> 'b pp)
  * ('a pp -> 'b dump)
  * ('a pp -> 'b show)
  * ('a pp -> 'b display)
  * ('a pp -> 'b print)
  * ('a pp -> 'b prerr)

  (** Like printers_of_show, but works on a function of type
  ['a show -> 'a show].  Note that due to technical limitations [print] and
  [prerr] functions take [pp] as inputs, not [print] or [prerr]. *)
  val ho_printers_of_show :
  ('a show -> 'b show) ->
  ('a pp -> 'b pp)
  * ('a pp -> 'b dump)
  * ('a pp -> 'b show)
  * ('a pp -> 'b display)
  * ('a pp -> 'b print)
  * ('a pp -> 'b prerr)
 *)

val show_int : int show
val dump_int : int dump
val pp_int : int pp
val display_int : int display

val show_bool : bool show
val dump_bool : bool dump
val pp_bool : bool pp
val display_bool : bool display

val show_char : char show
val dump_char : char dump
val pp_char : char pp
val display_char : char display
(* print_char, prerr_char exported through PP *)

val show_string : string show
val dump_string : string dump
val pp_string : string pp
val display_string : string display
(* print_string, prerr_string exported through PP *)

(** Returns the string as-is, without adding quotations or escapes.  *)
val show_string_verbatim : string show
(** Prints the string as-is, without adding quotations or escapes.  *)
val dump_string_verbatim : string dump
(** Prints the string as-is, without adding quotations or escapes.  *)
val pp_string_verbatim : string pp
(** Returns the string as-is, without adding quotations or escapes.  *)
val display_string_verbatim : string display
(** Prints the string as-is, without adding quotations or escapes.  *)
val print_string_verbatim : string print
(** Prints the string as-is, without adding quotations or escapes.  *)
val prerr_string_verbatim : string prerr

(** Returns a string containing the char as-is, without adding quotations or
    escapes.  *)
val show_char_verbatim : char show
(** Prints the char as-is, without adding quotations or escapes.  *)
val dump_char_verbatim : char dump
(** Prints the char as-is, without adding quotations or escapes.  *)
val pp_char_verbatim : char pp
(** Returns a string containing the char as-is, without adding quotations or
    escapes.  *)
val display_char_verbatim : char display
(** Prints the char as-is, without adding quotations or escapes.  *)
val print_char_verbatim : char print
(** Prints the char as-is, without adding quotations or escapes.  *)
val prerr_char_verbatim : char prerr

val show_float : float show
val dump_float : float dump
val pp_float : float pp
val display_float : float display
(* print_float, prerr_float exported through PP *)

(** Just returns "<fun>".  *)
val show_fun : ('a -> 'b) show
(** Just returns "<fun>".  *)
val dump_fun : ('a -> 'b) dump
(** Just prints "<fun>".  *)
val pp_fun : ('a -> 'b) pp
(** Just prints "<fun>".  *)
val display_fun : ('a -> 'b) display
(** Just prints "<fun>".  *)
val print_fun : ('a -> 'b) print
(** Just prints "<fun>".  *)
val prerr_fun : ('a -> 'b) prerr

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

(** [pp_seq opening delim closing pps] runs a sequence of pretty-printers with
    delimiters.

    - [opening] and [closing] are the opening and closing delimiters to output
    at the beginning and end.
    - [delim] is inserted along with a breakable space after every element
    except for the last one.
    - [pps] is a stream of pretty-printing functions (of type
    [(Format.formatter -> unit) Fstream.t]) to run.

 *)
val pp_seq :
  string
  -> string
  -> string
  -> Format.formatter
  -> (Format.formatter -> unit) Fstream.t
  -> unit

val pp_list : 'a pp -> 'a list pp
val dump_list : 'a pp -> 'a list dump
val show_list : 'a pp -> 'a list show
val display_list : 'a pp -> 'a list display
val print_list : 'a pp -> 'a list print
val prerr_list : 'a pp -> 'a list prerr

(** {!pp_seq} specialized to array-like data types.  Should be invoked like
    [pp_array_like length get opening closing] where [length], [get] should
    behave like the functions in the module [Array].  *)
val pp_array_like : ('a -> int) -> ('a -> int -> 'b) -> string -> string
  -> 'b pp -> 'a pp

val pp_array : 'a pp -> 'a array pp
val dump_array : 'a pp -> 'a array dump
val show_array : 'a pp -> 'a array show
val display_array : 'a pp -> 'a array display
val print_array : 'a pp -> 'a array print
val prerr_array : 'a pp -> 'a array prerr

open Bigarray
val pp_bigarray1 :
  ('a,'b) kind -> 'c layout -> 'a pp -> ('a,'b,'c) Array1.t pp
val dump_bigarray1 :
  ('a,'b) kind -> 'c layout -> 'a pp -> ('a,'b,'c) Array1.t dump
val show_bigarray1 :
  ('a,'b) kind -> 'c layout -> 'a pp -> ('a,'b,'c) Array1.t show
val display_bigarray1 :
  ('a,'b) kind -> 'c layout -> 'a pp -> ('a,'b,'c) Array1.t display
val print_bigarray1 :
  ('a,'b) kind -> 'c layout -> 'a pp -> ('a,'b,'c) Array1.t print
val prerr_bigarray1 :
  ('a,'b) kind -> 'c layout -> 'a pp -> ('a,'b,'c) Array1.t prerr

(** [pp_record fields] generates a [pp] implementation for a record type.  The
    argument [fields] is a list of the form ["foo", pp_foo] where [foo] is the
    name of a field of the record and [pp_foo] extracts that field and
    pretty-prints it.  *)
val pp_record : (string * 'a pp) list -> 'a pp

(** [pp_tuple fields] generates a [pp] implementation for a tuple type.  The
    argument [fields] is a list of functions of type ['a pp] that extracts a
    field of the tuple and pretty-prints it.  *)
val pp_tuple : 'a pp list -> 'a pp
val show_tuple : 'a pp list -> 'a show
val dump_tuple : 'a pp list -> 'a dump
val display_tuple : 'a pp list -> 'a display
val print_tuple : 'a pp list -> 'a print
val prerr_tuple : 'a pp list -> 'a prerr

val pp_pair : 'a pp -> 'b pp -> ('a * 'b) pp
val show_pair : 'a pp -> 'b pp -> ('a * 'b) show
val dump_pair : 'a pp -> 'b pp -> ('a * 'b) dump
val display_pair : 'a pp -> 'b pp -> ('a * 'b) display
val print_pair : 'a pp -> 'b pp -> ('a * 'b) print
val prerr_pair : 'a pp -> 'b pp -> ('a * 'b) prerr

val pp_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) pp
val show_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) show
val dump_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) dump
val display_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) display
val print_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) print
val prerr_triple : 'a pp -> 'b pp -> 'c pp -> ('a * 'b * 'c) prerr

module PrettyPrim : sig
  val print_string : string print
  val prerr_string : string prerr
  val print_float : float print
  val prerr_float : float prerr
  val print_char : char print
  val prerr_char : char prerr
end
