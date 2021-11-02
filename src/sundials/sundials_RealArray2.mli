(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

(** A two-dimensional matrix. The underlying data can be accessed as
    a {{:OCAML_DOC_ROOT(Bigarray.Array2.html)}Bigarray} via {!unwrap},
    but note that the first index specifies the column. *)
type t

(** An alias for two-dimensional
    {{:OCAML_DOC_ROOT(Bigarray.Array2.html)}Bigarray}s of floating-point
    numbers. *)
type data = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

(** [make nr nc v] returns an array with [nr] rows and [nc] columns, and
    with elements set to [v]. *)
val make : int -> int -> float -> t

(** [create nr nc] returns an uninitialized array with [nr] rows and [nc]
    columns. *)
val create : int -> int -> t

(** [of_lists xxs] constructs an array from lists of rows. The inner lists
    must all have the same length. *)
val of_lists : (float list) list -> t

(** [of_lists xxs] constructs an array from an array of rows. The inner arrays
    must all have the same length. *)
val of_arrays : (float array) array -> t

(** An array with no elements. *)
val empty : t

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

(** Pretty-print an array using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

(** Pretty-print an array using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module.
    The defaults are: [start="\["], [stop="\]"], [rowsep=";"],
    [indent=4], [sep=" "], and
    [item=fun f r c->Format.fprintf f "(%2d,%2d)=% -15e" r c] (see
    {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}).
    The [indent] argument specifies the indent for wrapped rows. *)
val ppi : ?start:string -> ?rowstart:string
          -> ?stop:string -> ?rowstop:string
          -> ?sep:string -> ?rowsep:string
          -> ?item:(Format.formatter -> int -> int -> float -> unit)
          -> unit
          -> Format.formatter -> t -> unit

(** Creates a new array with the same contents as an existing one. *)
val copy : t -> t

(** Copy the first array into the second one.
    See {{:OCAML_DOC_ROOT(Bigarray.Genarray.html#VALblit)}
    [Bigarray.Genarray.blit]} for more details. *)
val blit : src:t -> dst:t -> unit

(** [fill a c] sets all elements of [a] to the constant [c]. *)
val fill : t -> float -> unit

(** [make m n] returns an uninitialized [m] by [n] array. *)
val make_data : int -> int -> data

(** Creates a new matrix from an existing {!data} array. Changes to
    one affect the other since they share the same underlying storage. *)
val wrap : data -> t

(** Returns the {!data} array behind a matrix. Changes to one affect
    the other since they share the same underlying storage. Note that the
    array is accessed column-first, that is,
    [get a i j = (unwrap a).{j, i}]. *)
val unwrap : t -> data

