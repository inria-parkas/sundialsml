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

(** A {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} Bigarray} of floats. *)
type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

(** [make n x] returns an array with [n] elements each set to [x]. *)
val make : int -> float -> t

(** [create n] returns an uninitialized array with [n] elements. *)
val create : int -> t

(** An array with no elements. *)
val empty : t

(** [init n f] returns an array with [n] elements, with element [i]
    set to [f i]. *)
val init : int -> (int -> float) -> t

(** Pretty-print an array using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

(** Pretty-print an array using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module.
    The defaults are: [start="\["], [stop="\]"], [sep=" "], and
    [item=fun f->Format.fprintf f "%2d=% -15e"] (see
  {{:OCAML_DOC_ROOT(Format.html#VALfprintf)} fprintf}). *)
val ppi : ?start:string -> ?stop:string -> ?sep:string
          -> ?item:(Format.formatter -> int -> float -> unit)
          -> unit
          -> Format.formatter -> t -> unit

(** Creates an array by copying the contents of a
    {{:OCAML_DOC_ROOT(Array.html)} [float array]}. *)
val of_array : float array -> t

(** Creates an array by copying the contents of a
    {{:OCAML_DOC_ROOT(List.html)} [float list]}. *)
val of_list : float list -> t

(** Copies into a new {{:OCAML_DOC_ROOT(Array.html)} [float array]}. *)
val to_array : t -> float array

(** Copies into an existing
    {{:OCAML_DOC_ROOT(Array.html)} [float array]}. *)
val into_array : t -> float array -> unit

(** Copies into a {{:OCAML_DOC_ROOT(List.html)} [float list]}. *)
val to_list : t -> float list

(** Creates a new array with the same contents as an existing one. *)
val copy : t -> t

(** Access a sub-array of the given array without copying. *)
val sub : t -> int -> int -> t

(** [blit_some src isrc dst idst len] copies [len] elements of [src] at
    offset [isrc] to [dst] at offset [idst].

    @raise Invalid_argument "RealArray.blit_some" if [isrc], [idst], and
    [len] do not specify valid subarrays of [src] and [dst]. *)
val blit_some : t -> int -> t -> int -> int -> unit

(** Copy the first array into the second one.
    See {{:OCAML_DOC_ROOT(Bigarray.Genarray.html#VALblit)}
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

