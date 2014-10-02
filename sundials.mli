(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(** Generic definitions for the Sundials suite.

 @version VERSION()
 @author Timothy Bourke (Inria)
 @author Jun Inoue (Inria)
 @author Marc Pouzet (LIENS)
 *)

(** {2 Constants} *)

(** [true] iff this binding is compiled with BLAS/LAPACK support.  *)
val blas_lapack_supported : bool

(** The BIG_REAL constant.
    @cvode <node5#s:types> Data Types
 *)
val big_real : float

(** The UNIT_ROUNDOFF constant.
    @cvode <node5#s:types> Data Types
 *)
val unit_roundoff : float

(** {2 Exceptions} *)

(** This exception may be thrown inside most callback functions to indicate a
    recoverable failure. For callbacks that return a boolean, the accompanying
    value is used as the result, otherwise it is ignored. Throwing any other
    kind of exception normally indicates an unrecoverable failure. *)
exception RecoverableFailure

(** This exception may be thrown from error-weight functions to
    indicate that an error weight became non-positive.  See
    [WFtolerances] in {!Cvode.tolerance} or {!Ida.tolerance}. *)
exception NonPositiveEwt

(** {2 Nvectors} *)

(** The type representing an nvector with underlying data of type ['data]. The
    additional type argument, ['kind], is typically either [serial], [parallel],
    or [custom]. It is needed because some linear solvers make additional
    assumptions about the underlying vector representation.

    @cvode <node7#s:nvector> N_Vector *)
type ('data, 'kind) nvector

(** [unvec nv] returns the data underlying the nvector [nv]. *)
val unvec : ('data, 'kind) nvector -> 'data

(** {2 Arrays} *)

module RealArray :
  sig
    (** Utilities for manipulating 1-dimensional bigarrays of floats,
        suitable for use with sundials.  Most of the functions in this
        submodule mimic those of {{:OCAML_DOC_ROOT(Array)} (Array)},
        with the exceptions of {!create} and {!sub}, which are modeled
        after {{:OCAML_DOC_ROOT(Bigarray.Array1)} (Bigarray)}.  *)

    (** A {{:OCAML_DOC_ROOT(Bigarray.Array1)} (Bigarray)} vector of floats. *)
    type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

    (** [make n x] returns an array with [n] elements each initialized to x. *)
    val make : int -> float -> t

    (** [create n] returns an array with [n] elements. *)
    val create : int -> t

    (** [init n f] returns an array with [n] elements, with element [i]
        initialized to [f i]. *)
    val init : int -> (int -> float) -> t

    (** Creates an array by copying the contents of an {{:OCAML_DOC_ROOT(Array)}
        Array}. *)
    val of_array : float array -> t

    (** Creates an array by copying the contents of a
        {{:OCAML_DOC_ROOT(List)} List}. *)
    val of_list : float list -> t

    (** Copies into a new {{:OCAML_DOC_ROOT(Array)} Array}. *)
    val to_array : t -> float array

    (** Copies into an existing {{:OCAML_DOC_ROOT(Array)} Array}. *)
    val into_array : t -> float array -> unit

    (** Copies into a {{:OCAML_DOC_ROOT(List)} List}. *)
    val to_list : t -> float list

    (** Create a new array with the same contents as an existing one. *)
    val copy : t -> t

    (** Extract a sub-array of the given real array.  This function
        imitates that of Bigarray rather than Array, i.e. it doesn't
        copy anything.  *)
    val sub : t -> int -> int -> t

    (** [blit src isrc dst idst len] copies [len] elements of [src] at
        offset [isrc] to [dst] at offset [idst], i.e. it writes
        [src.{isrc}], [src.{isrc+1}], ... [src.{isrc+len-1}]
        to
        [dst.{idst}], [dst.{idst+1}], ... [dst.{idst+len-1}].

        @raise Invalid_argument "RealArray.blit" if [isrc], [idst], and
        [len] don't specify valid subarrays of [src] and [dst].
      *)
    val blit : t -> int -> t -> int -> int -> unit

    (** An alias of [Bigarray.Array1.blit].  *)
    val blit_all : t -> t -> unit

    (** [fill a c] sets all elements of the array [a] to the constant [c]. *)
    val fill : t -> float -> unit

    (** Returns the length of an array *)
    val length : t -> int

    (** [print_with_time t a] prints a line containing the current time (see
        {!print_time}) followed by a tab-delimited list of the values of [a],
        and then a newline. See also {!extra_precision}. *)
    val print_with_time : float -> t -> unit

    (** [fold_left f b a] returns [f (f (f b a.{0}) a.{1}) ...]. *)
    val fold_left : ('a -> float -> 'a) -> 'a -> t -> 'a

    (** [fold_right f b a] returns [f (f (f b a.{0}) a.{1}) ...]. *)
    val fold_right : (float -> 'a -> 'a) -> t -> 'a -> 'a

    (** [iter f a] applies [f] to the values of each element in [a]. *)
    val iter : (float -> unit) -> t -> unit

    (** [iteri f a] applies [f] to the indexes and values of each element
        in [a]. *)
    val iteri : (int -> float -> unit) -> t -> unit

    (** [map f a] applies [f] to the value of each element in [a] and
        stores the result back into the same element. *)
    val map : (float -> float) -> t -> unit

    (** [map f a] applies [f] to the index and value of each element
        in [a] and stores the result back into the same element. *)
    val mapi : (int -> float -> float) -> t -> unit
  end

module RealArray2 :
  sig
    (** A {{:OCAML_DOC_ROOT(Bigarray.Array2)} (Bigarray)} 2D vector of floats. *)
    type data = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

    (** [make m n] returns an [m] by [n] array. *)
    val make_data : int -> int -> data

    (** The underlying data is stored as a two-dimensional
       {{:OCAML_DOC_ROOT(Bigarray.Array2)}Bigarray} of floats.
       Note that, in the underlying array, Sundials stores columns
       in the first dimension and rows in the second. So, the value
       at row [i] and column [j] in an array [m] is [m.{j}.{i}].
     *)
    type t

    (** [make nr nc v] creates an [nr] by [nc] wrapped array with all elements
        set to [v]. *)
    val make : int -> int -> float -> t

    (** [create nr nc] creates an [nr] by [nc] wrapped array. *)
    val create : int -> int -> t

    (** [get a i j] gives the value of the [(i, j)]th element of [a]. *)
    val get : t -> int -> int -> float

    (** [col a j] slices the [j]th column of [a].  The slice shares
        storage with [a]. *)
    val col : t -> int -> RealArray.t

    (** [set a i j v] sets the value of the [(i, j)]th element of [a] to [v]. *)
    val set : t -> int -> int -> float -> unit

    (** [nr, nc = size a] gives the number of rows, [nr], and the number of
        columns, [nc], in [a] *)
    val size : t -> int * int

    (** [copy a] creates a copy of [a] and its underlying {!data} array. *)
    val copy : t -> t

    (** [copyinto a b] copies the contents of [a] into [b]. Both arrays
        must have the same dimensions. *)
    val copyinto : t -> t -> unit

    (** Creates a new array from an existing {!data} array; changes to either
        array affect the other (i.e., they share the same underlying storage). *)
    val wrap : data -> t

    (** Returns the underlying {!data} array; changes to either array affect the
        other (i.e., they share the same underlying storage). *)
    val unwrap : t -> data
  end

module LintArray :
  sig
    (** A {{:OCAML_DOC_ROOT(Bigarray.Array1)} (Bigarray)} vector of integers. *)
    type t = (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t

    (** [make n v] returns an new array with [n] elements each set to [v]. *)
    val make  : int -> int -> t

    (** [create n] returns an uninitialized new array with [n] elements. *)
    val create  : int -> t
  end

(** {2 Arrays of roots (zero-crossings)} *)

(** Utility functions for arrays of roots (zero-crossings). *)
module Roots :
  sig
    type t
    type val_array = RealArray.t

    type root_event =
      | NoRoot      (** No root (0)       *)
      | Rising      (** Rising root (1)   *)
      | Falling     (** Falling root (-1) *)

    (** [create n] returns an array with [n] elements, each set to NoRoot. *)
    val create : int -> t

    (** [make n x] returns an array with [n] elements, each set to [x]. *)
    val make : int -> root_event -> t

    (** Returns the length of an array *)
    val length : t -> int

    (** [detected r i] returns [true] if the value of the [i]th element of [r]
        is either Rising or Falling. *)
    val detected : t -> int -> bool

    (** [rising r i] returns [true] if the value of the [i]th element of [r] is
        Rising. *)
    val rising : t -> int -> bool

    (** [falling r i] returns [true] if the value of the [i]th element of [r] is
        Falling. *)
    val falling : t -> int -> bool

    (** [get r i] returns the value of the [i]th element of [r]. *)
    val get : t -> int -> root_event

    (** [set r i v] sets the value of the [i]th element of [r]. *)
    val set : t -> int -> root_event -> unit

    (** [copy r] creates a new array with the contents as [r]. *)
    val copy : t -> t

    (** [set_noroot r i] sets the value of the [i]th element of [r] to
        NoRoot.  *)
    val set_noroot : t -> int -> unit

    (** [set_rising r i] sets the value of the [i]th element of [r] to
        Rising.  *)
    val set_rising : t -> int -> unit

    (** [set_falling r i] sets the value of the [i]th element of [r] to
        Falling. *)
    val set_falling : t -> int -> unit

    (** Returns 0 for NoRoot, 1 for Rising, and -1 for Falling. *)
    val int_of_root_event : root_event -> int

    (** Resets all elements to NoRoot. *)
    val reset : t -> unit

    (** [print r] prints a line containing a tab-delimited list of the values of
        [r] (by their constructor names), and then a newline. *)
    val print : t -> unit

    (** [string_of_root_event r] returns the name of the data constructor [r]
        of type [root_event] as a string. *)
    val string_of_root_event : root_event -> string

    (** Returns [true] if any elements are equal to Rising or Falling. *)
    val exists : t -> bool

    (** [iter f r] applies [f] to the values of each element in
        [r]. *)
    val iter : (root_event -> unit) -> t -> unit

    (** [iteri f r] applies [f] to the indexes and values of each element
        in [r]. *)
    val iteri : (int -> root_event -> unit) -> t -> unit

    (** Makes a [Roots.t] from a list of root events.  *)
    val of_list : root_event list -> t

    (** Copies the contents of an {{:OCAML_DOC_ROOT(Array)} Array} into an
        opaque array of type [Roots.t].  *)
    val of_array : root_event array -> t

    (** Copies the contents of an opaque array of type [Roots.t] into an
        {{:OCAML_DOC_ROOT(Array)} Array}.  *)
    val to_array : t -> root_event array

    (** Copies the contents of a [Roots.t] into a list.  *)
    val to_list : t -> root_event list

    (** [fill_all a x] sets the values of [a] to [x] everywhere. *)
    val fill_all : t -> root_event -> unit

    (** [fill a i len x] sets the values of [a] from [i] through [i+len-1] to
        [x]. *)
    val fill : t -> int -> int -> root_event -> unit
  end

(** Utility functions for arrays of directions to detect on root functions
    (increasing/decreasing/either). *)
module RootDirs :
  sig
    type t

    type root_direction =
      | Increasing                      (** Monitor rising zero-crossings *)
      | Decreasing                      (** Monitor falling zero-crossings *)
      | IncreasingOrDecreasing          (** Monitor all zero-crossings *)

    (** [string_of_root_direction d] returns d as a human-readable string.  *)
    val string_of_root_direction : root_direction -> string

    (** [make n] returns an array with [n] elements, each set to the specified
        value. *)
    val make : int -> root_direction -> t

    (** [create n] returns an array with [n] elements, each set to
        IncreasingOrDecreasing. *)
    val create : int -> t

    (** [copy_n n a] returns a fresh array with [n] elements, initialized from
        the contents of a.  If [n > Array.length a], then the extra space is
        initialized to IncreasingOrDecreasing.  *)
    val copy_n : int -> root_direction array -> t

    (** Returns the length of an array *)
    val length : t -> int

    (** [get r i] returns the value of the [i]th element of [r]. *)
    val get : t -> int -> root_direction

    (** [set r i v] sets the value of the [i]th element of [r]. *)
    val set : t -> int -> root_direction -> unit

    (** [fill a i len x] sets the values of [a] from [i] through [i+len-1] to
        [x]. *)
    val fill : t -> int -> int -> root_direction -> unit

    (** [fill_all a x] sets the values of [a] to [x] everywhere. *)
    val fill_all : t -> root_direction -> unit

    (** [blit a oa b ob len] copies the values of [a] at indices
        [oa, oa+1, ..., oa+len-1] to [b] at indices
        [ob, ob+1, ..., ob+len-1]. *)
    val blit : t -> int -> t -> int -> int -> unit

    (** [blit_all a b] copies the values of [a] to [b].  If
        [length a > length b], then [b] is filled with a prefix of [a].
        If [length a < length b], then only a prefix of [b] is modified.  *)
    val blit_all : t -> t -> unit

    (** [init n f] creates an array of length [n] and sets it to [f i] for each
        index [i]. *)
    val init : int -> (int -> root_direction) -> t

    (** Makes a [RootDirs.t] from a list of root events.  *)
    val of_list : root_direction list -> t

    (** Copies the contents of an {{:OCAML_DOC_ROOT(Array)} Array} into an
        opaque array of type [RootDirs.t].  *)
    val of_array : root_direction array -> t

    (** Copies the contents of an opaque array of type [RootDirs.t] into an
        {{:OCAML_DOC_ROOT(Array)} Array}.  *)
    val to_array : t -> root_direction array

    (** Copies the contents of a [RootDirs.t] into a list.  *)
    val to_list : t -> root_direction list
  end

(** {2 Solver results and error reporting} *)

(**
 Possible values returned when a CVODE/IDA solver step function succeeds.
 Failures are indicated by exceptions.

 @cvode <node5#sss:cvode> CVode
 @ida <node5#sss:ida> IDASolve
 *)
type solver_result =
  | Continue            (** CV_SUCCESS / IDA_SUCCESS *)
  | RootsFound          (** CV_ROOT_RETURN / IDA_ROOT_RETURN *)
  | StopTimeReached     (** CV_TSTOP_RETURN / IDA_TSTOP_RETURN *)

(**
 Type of values passed to a registered error handler function.

 @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
 @ida <node5#sss:optin_main> IDASetErrHandlerFn
 *)
type error_details = {
    error_code : int;
    module_name : string;               (** IDA, CVODE, CVSPGMR, etc. *)
    function_name : string;
    error_message : string;
  }

(** {2 Miscellaneous utility functions} *)

(** [print_time (s1, s2) t] prints [t] with [s1] on the left and [s2] on the
    right.  *)
val print_time : string * string -> float -> unit

(** Controls the precision of {!print_time} and {!RealArray.print_with_time}.
 
    If [true] the format [%.15e] is used, otherwise [%e]
    (the default) is used. *)
val extra_precision : bool ref

(** [format_float fmt f] formats [f] according to the format string [fmt],
    using the low-level [caml_format_float] function. *)
val format_float : string -> float -> string

(** Equivalent to [format_float "%a"].
  
    [floata f] returns the bit-level representation of [f] in
    hexadecimal as a string. *)
val floata : float -> string

