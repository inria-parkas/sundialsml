(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Generic definitions, arrays, and utility functions.

 @version VERSION()
 @author Timothy Bourke (Inria/ENS)
 @author Jun Inoue (Inria/ENS)
 @author Marc Pouzet (UPMC/ENS/Inria)
 *)

(** Installation specific constants. *)
module Config = Sundials_Config

(** Index values for sparse matrices. *)
module Index = Sundials_Index

(** {2:exceptions Exceptions} *)

(** Indicates a recoverable failure within a callback function.
    Any other exception normally indicates an unrecoverable failure. *)
exception RecoverableFailure

(** Raised by error-weight functions on non-positive error weights. See
    {{!Cvode.tolerance}[Cvode.WFtolerances]} or
    {{!Ida.tolerance}[Ida.WFtolerances]}. *)
exception NonPositiveEwt

(** {2:arrays Arrays} *)

(** Vectors of floats (one-dimensional bigarrays). *)
module RealArray = Sundials_RealArray

(** Matrices of floats (wrappers around two-dimensional bigarrays with
    extra internal information required by Sundials). *)
module RealArray2 = Sundials_RealArray2

(** Vectors of integers (one-dimensional bigarrays). *)
module LintArray = Sundials_LintArray

(** {2:roots Arrays of roots (zero-crossings)} *)

(** Vectors of root (zero-crossing) statuses. *)
module Roots : sig (* {{{ *)
  (** Arrays that communicate the occurrence of zero-crossings. The
      underlying representation is hidden to isolate compatability
      issues related to integers. *)
  type t

  (** Values indicating the status of root functions.
      @cvode <node5#sss:optout_root> CVodeGetRootInfo
      @ida <node5#sss:optout_root> IdaGetRootInfo *)
  type r =
    | NoRoot      (** No root was found on the corresponding function (0). *)
    | Rising      (** The corresponding root function is increasing (1). *)
    | Falling     (** The corresponding root function is decreasing (-1). *)

  (** [create n] returns an array with [n] elements each set to
      {{!r}NoRoot}. *)
  val create : int -> t

  (** [make n x] returns an array with [n] elements each set to [x]. *)
  val make : int -> r -> t

  (** [init n f] returns an array with [n] elements, with element [i] set
      to [f i]. *)
  val init : int -> (int -> r) -> t

  (** Returns the length of an array. *)
  val length : t -> int

  (** Pretty-print a root array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
  val pp : Format.formatter -> t -> unit

  (** Pretty-print a root array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module.
      The defaults are: [start="\["], [stop="\]"], [sep="; "], and
      [item] prints '_' for {{!r}NoRoot}, 'R' for {{!r}Rising}, and 'F' for
      {{!r}Falling}. *)
  val ppi : ?start:string -> ?stop:string -> ?sep:string
            -> ?item:(Format.formatter -> int -> r -> unit)
            -> unit
            -> Format.formatter -> t -> unit

  (** Returns [true] only if the specified element is either {{!r}Rising} or
      {{!r}Falling}. *)
  val detected : t -> int -> bool

  (** Returns [true] only if the specified element is {{!r}Rising}. *)
  val rising : t -> int -> bool

  (** Returns [true] only if the specified element is {{!r}Falling}. *)
  val falling : t -> int -> bool

  (** [get r i] returns the [i]th element of [r]. *)
  val get : t -> int -> r

  (** [set r i v] sets the [i]th element of [r] to [v]. *)
  val set : t -> int -> r -> unit

  (** [set_noroot r i] sets the [i]th element of [r] to {{!r}NoRoot}. *)
  val set_noroot : t -> int -> unit

  (** [set_rising r i] sets the [i]th element of [r] to {{!r}Rising}. *)
  val set_rising : t -> int -> unit

  (** [set_falling r i] sets the [i]th element of [r] to {{!r}Falling}. *)
  val set_falling : t -> int -> unit

  (** [fill a x] sets all elements in [a] to [x]. *)
  val fill : t -> r -> unit

  (** Creates a new array with the same contents as an existing one. *)
  val copy : t -> t

  (** Returns [0] for {{!r}NoRoot}, [1] for {{!r}Rising}, and [-1] for
      {{!r}Falling}. *)
  val int_of_root : r -> int

  (** Resets all elements to {{!r}NoRoot}. *)
  val reset : t -> unit

  (** [true] if any elements are equal to {{!r}Rising} or {{!r}Falling}. *)
  val exists : t -> bool

  (** [iter f r] successively applies [f] to each element in [r]. *)
  val iter : (r -> unit) -> t -> unit

  (** [iteri f r] successively applies [f] to the indexes and
      elements of [r]. *)
  val iteri : (int -> r -> unit) -> t -> unit

  (** Creates an array by copying the contents of a [r list]. *)
  val of_list : r list -> t

  (** Copies into a list. *)
  val to_list : t -> r list

  (** Creates a new value from the contents of an
      {{:OCAML_DOC_ROOT(Array.html)} array}. *)
  val of_array : r array -> t

  (** Creates a new array from the contents of a given value. *)
  val to_array : t -> r array
end (* }}} *)

(** Vectors of root (zero-crossing) directions. *)
module RootDirs : sig (* {{{ *)
  type t
  (** Arrays that communicate which zero-crossings are sought. The
      underlying representation is hidden to isolate compatability
      issues related to integers. *)

  (** Values indicating which types of roots are sought.

      @cvode <node5#sss:optin_root> CVodeSetRootDirection
      @ida <node5#sss:optin_root> IdaSetRootDirection *)
  type d =
    | Increasing              (** Only look for rising zero-crossings. *)
    | Decreasing              (** Only look for falling zero-crossings. *)
    | IncreasingOrDecreasing  (** Look for any zero-crossing. *)

  (** [make n x] returns an array with [n] elements each set to [x]. *)
  val make : int -> d -> t

  (** [create n] returns an array with [n] elements each set to
      {{!d}IncreasingOrDecreasing}. *)
  val create : int -> t

  (** [init n f] returns an array with [n] elements, with element [i] set
      to [f i]. *)
  val init : int -> (int -> d) -> t

  (** Pretty-print a root direction array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
  val pp : Format.formatter -> t -> unit

  (** Pretty-print a root direction array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module.
      The defaults are: [start="\["], [stop="\]"], [sep="; "], and
      [item] prints 'R' for {{!d}Increasing}, 'F' for {{!d}Decreasing}, and
      'E' (either) for {{!d}IncreasingOrDecreasing}. *)
  val ppi : ?start:string -> ?stop:string -> ?sep:string
            -> ?item:(Format.formatter -> int -> d -> unit)
            -> unit
            -> Format.formatter -> t -> unit

  (** [copy n a] returns an array with [n] elements, initialized from
      the contents of a. If [n > Array.length a] then the extra space is
      initialized to {{!d}IncreasingOrDecreasing}. *)
  val copy : int -> d array -> t

  (** Returns the length of an array *)
  val length : t -> int

  (** [get r i] returns the [i]th element of [r]. *)
  val get : t -> int -> d

  (** [set r i v] sets the [i]th element of [r] to [v]. *)
  val set : t -> int -> d -> unit

  (** [fill_all a x] sets the values of [a] to [x] everywhere. *)
  val fill : t -> d -> unit

  (** [blit_some src isrc dst idst len] copies [len] elements of [src] at
      offset [isrc] to [dst] at offset [idst].

      @raise Invalid_argument "RootDirs.blit_some" if [isrc], [idst], and
      [len] do not specify valid subarrays of [src] and [dst]. *)
  val blit_some : t -> int -> t -> int -> int -> unit

  (** Copy the first array into the second one.
      See {{:OCAML_DOC_ROOT(Bigarray.Genarray.html#VALblit)}
      [Bigarray.Genarray.blit]} for more details. *)
  val blit : t -> t -> unit

  (** Creates an array by copying the contents of a [d list]. *)
  val of_list : d list -> t

  (** Copies into a list. *)
  val to_list : t -> d list

  (** Creates a new value from the contents of an
      {{:OCAML_DOC_ROOT(Array.html)} array}. *)
  val of_array : d array -> t

  (** Creates a new array from the contents of a given value. *)
  val to_array : t -> d array
end (* }}} *)

(** {2:constraints Constraints} *)

(** Symbolic names for variable constraints. These names describe
    the constants passed to {!Ida.set_constraints} and
    {!Kinsol.set_constraints}.

    @ida <node5#sss:optin_main> IDASetConstraints
    @kinsol <node5#ss:optin_main> KINSetConstraints *)
module Constraint : sig (* {{{ *)
  (** The constant [0.0]. *)
  val unconstrained : float

  (** The constant [1.0]. *)
  val geq_zero : float

  (** The constant [-1.0]. *)
  val leq_zero : float

  (** The constant [2.0]. *)
  val gt_zero : float

  (** The constant [-2.0]. *)
  val lt_zero : float

  (** For pattern-matching on constraints. See {!of_float}. *)
  type t =
  | Unconstrained   (** [true] *)
  | GeqZero         (** [>= 0] *)
  | LeqZero         (** [<= 0] *)
  | GtZero          (** [> 0]  *)
  | LtZero          (** [< 0]  *)

  (** Map constraint values to floating-point constants. *)
  val to_float : t -> float

  (** Map floating-point constants to constraint values.

      @raise Invalid_argument The given value is not a legal constraint. *)
  val of_float : float -> t
end (* }}} *)

(** {2:matlin Matrices and Linear Solvers} *)

module Matrix = Sundials_Matrix

module LinearSolver = Sundials_LinearSolver

(** {2:results Solver results and error reporting} *)

(** Files for error and diagnostic information. File values are passed
    to functions like {!Cvode.set_error_file} and {!Kinsol.set_info_file} to
    log solver errors and diagnostics. *)
module Logfile = Sundials_Logfile

module Util : sig (* {{{ *)

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

  (** {2:misc Miscellaneous utility functions} *)

  (** [format_float fmt f] formats [f] according to the format string [fmt].
      It uses the low-level [caml_format_float] function. *)
  val format_float : string -> float -> string

  (** Returns the bit-level representation of a float in hexadecimal as a string.
      Equivalent to [format_float "%a"]. *)
  val floata : float -> string

end (* }}} *)

