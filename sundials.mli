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

(** Generic definitions, arrays, and utility functions.

 @version VERSION()
 @author Timothy Bourke (Inria)
 @author Jun Inoue (Inria)
 @author Marc Pouzet (LIENS)
 *)

(** {2 Constants} *)

(** Indicates whether the interface was compiled with BLAS/LAPACK support. *)
val blas_lapack_supported : bool

(** The largest value representable as a real.

    @cvode <node5#s:types> Data Types
 *)
val big_real : float

(** The smallest value representable as a real.

    @cvode <node5#s:types> Data Types
 *)
val small_real : float

(** The difference bewteen 1.0 and the minimum real greater than 1.0.

    @cvode <node5#s:types> Data Types
 *)
val unit_roundoff : float

(** {2 Exceptions} *)

(** Indicates a recoverable failure within a callback function.
    Any other exception normally indicates an unrecoverable failure. *)
exception RecoverableFailure

(** Raised by error-weight functions on non-positive error weights. See
    {{!Cvode.tolerance}[Cvode.WFtolerances]} or
    {{!Ida.tolerance}[Ida.WFtolerances]}. *)
exception NonPositiveEwt

(** {2 Arrays} *)

(** Vectors of floats (one-dimensional bigarrays). *)
module RealArray :
  sig
    (** A {{:OCAML_DOC_ROOT(Bigarray.Array1)} Bigarray} of floats. *)
    type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

    (** [make n x] returns an array with [n] elements each set to [x]. *)
    val make : int -> float -> t

    (** [create n] returns an uninitialized array with [n] elements. *)
    val create : int -> t

    (** [init n f] returns an array with [n] elements, with element [i]
        set to [f i]. *)
    val init : int -> (int -> float) -> t

    (** Creates an array by copying the contents of a {{:OCAML_DOC_ROOT(Array)}
        float array}. *)
    val of_array : float array -> t

    (** Creates an array by copying the contents of a
        {{:OCAML_DOC_ROOT(List)} [float list]}. *)
    val of_list : float list -> t

    (** Copies into a new {{:OCAML_DOC_ROOT(Array)} [float array]}. *)
    val to_array : t -> float array

    (** Copies into an existing {{:OCAML_DOC_ROOT(Array)} [float array]}. *)
    val into_array : t -> float array -> unit

    (** Copies into a {{:OCAML_DOC_ROOT(List)} [float list]}. *)
    val to_list : t -> float list

    (** Create a new array with the same contents as an existing one. *)
    val copy : t -> t

    (** Access a sub-array of the given array without copying. *)
    val sub : t -> int -> int -> t

    (** [blit_some src isrc dst idst len] copies [len] elements of [src] at
        offset [isrc] to [dst] at offset [idst].

        @raise Invalid_argument "RealArray.blit_some" if [isrc], [idst], and
        [len] do not specify valid subarrays of [src] and [dst]. *)
    val blit_some : t -> int -> t -> int -> int -> unit

    (** Copy the first array into the second one.
        See {{:OCAML_DOC_ROOT(Bigarray.Genarray#VALblit)}
        [Bigarray.Genarray.blit]} for more details. *)
    val blit : t -> t -> unit

    (** [fill a c] sets all elements of [a] to the constant [c]. *)
    val fill : t -> float -> unit

    (** Returns the length of an array. *)
    val length : t -> int

    (** [fold_left f b a] returns [f (f (f b a.{0}) a.{1}) ...)]. *)
    val fold_left : ('a -> float -> 'a) -> 'a -> t -> 'a

    (** [fold_right f b a] returns [(f ... (f a.{n-2} (f a.{n-1} b)))]. *)
    val fold_right : (float -> 'a -> 'a) -> t -> 'a -> 'a

    (** [iter f a] successively applies [f] to the elements of [a]. *)
    val iter : (float -> unit) -> t -> unit

    (** [iteri f a] successively applies [f] to the indexes and values
        of [a]. *)
    val iteri : (int -> float -> unit) -> t -> unit

    (** [map f a] replaces each element [a.{i}] with [f a.{i}]. *)
    val map : (float -> float) -> t -> unit

    (** [map f a] replaces each element [a.{i}] with [f i a.{i}]. *)
    val mapi : (int -> float -> float) -> t -> unit
  end

(** Matrices of floats (wrappers around two-dimensional bigarrays). *)
module RealArray2 :
  sig
    (** A two-dimensional matrix. The underlying data can be accessed as
        a {{:OCAML_DOC_ROOT(Bigarray.Array2)}Bigarray} via {!unwrap},
        but note that the first index specifies the column. *)
    type t

    (** An alias for the underlying
        {{:OCAML_DOC_ROOT(Bigarray.Array2)}Bigarray}. *)
    type data = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

    (** [make nr nc v] returns an array with [nr] rows and [nc] columns, and
        with elements set to [v]. *)
    val make : int -> int -> float -> t

    (** [create nr nc] returns an uninitialized array with [nr] rows and [nc]
        columns. *)
    val create : int -> int -> t

    (** [get a i j] returns the value at row [i] and column [j] of [a]. *)
    val get : t -> int -> int -> float

    (** [col a j] returns the [j]th column of [a]. The slice shares storage
        with [a]. *)
    val col : t -> int -> RealArray.t

    (** [set a i j v] sets the value at row [i] and column [j] of [a] to [v]. *)
    val set : t -> int -> int -> float -> unit

    (** [nr, nc = size a] returns the numbers of rows [nr] and columns [nc]
        of [a] *)
    val size : t -> int * int

    (** Creates a new array with the same contents as an existing one. *)
    val copy : t -> t

    (** Copy the first array into the second one.
        See {{:OCAML_DOC_ROOT(Bigarray.Genarray#VALblit)}
        [Bigarray.Genarray.blit]} for more details. *)
    val blit : t -> t -> unit

    (** [make m n] returns an uninitialized [m] by [n] array. *)
    val make_data : int -> int -> data

    (** Creates a new matrix from an existing {!data} array. Changes to one
        affect the other since they share the same underlying storage. *)
    val wrap : data -> t

    (** Returns the {!data} array behind a matrix. Changes to one affect the
        other since they share the same underlying storage. Note that the
        array is accessed column-first, that is,
        [get a i j = (unwrap a).{j, i}]. *)
    val unwrap : t -> data
  end

(** Vectors of integers (one-dimensional bigarrays). *)
module LintArray :
  sig
    (** A {{:OCAML_DOC_ROOT(Bigarray.Array1)} Bigarray} of integers. *)
    type t = (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t

    (** [make n x] returns an array with [n] elements each set to [v]. *)
    val make  : int -> int -> t

    (** [create n] returns an uninitialized array with [n] elements. *)
    val create  : int -> t
  end

(** {2 Arrays of roots (zero-crossings)} *)

(** Vectors of root (zero-crossing) statuses. *)
module Roots :
  sig
    type t

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

(** Vectors of root (zero-crossing) directions. *)
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

    (** [blit_some a oa b ob len] copies the values of [a] at indices
        [oa, oa+1, ..., oa+len-1] to [b] at indices
        [ob, ob+1, ..., ob+len-1]. *)
    val blit_some : t -> int -> t -> int -> int -> unit

    (** [blit a b] copies the values of [a] to [b].  If
        [length a > length b], then [b] is filled with a prefix of [a].
        If [length a < length b], then only a prefix of [b] is modified.  *)
    val blit : t -> t -> unit

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

(** Result of a successful {{!Cvode.step_normal}CVODE} or
    {{!Ida.step_normal}IDA} step. Failures are indicated by exceptions.
 
 @cvode <node5#sss:cvode> CVode
 @ida <node5#sss:idasolve> IDASolve *)
type solver_result =
  | Continue            (** CV_SUCCESS / IDA_SUCCESS *)
  | RootsFound          (** CV_ROOT_RETURN / IDA_ROOT_RETURN *)
  | StopTimeReached     (** CV_TSTOP_RETURN / IDA_TSTOP_RETURN *)

(** Information passed to registered error handler functions.
    See {!Cvode.set_err_handler_fn}, {!Ida.set_err_handler_fn}, and
    {!Kinsol.set_err_handler_fn}.

 @cvode <node5#ss:ehFn> CVodeErrHandlerFn
 @ida <node5#ss:ehFn> IDAErrHandlerFn
 @kinsol <node5#ss:ehFn> KINErrHandlerFn *)
type error_details = {
    error_code : int;
    module_name : string;        (** IDA, CVODE, CVSPGMR, etc. *)
    function_name : string;
    error_message : string;
  }

(** {2 Miscellaneous utility functions} *)

(** [format_float fmt f] formats [f] according to the format string [fmt].
    It uses the low-level [caml_format_float] function. *)
val format_float : string -> float -> string

(** Returns the bit-level representation of a float in hexadecimal as a string.
    Equivalent to [format_float "%a"]. *)
val floata : float -> string

