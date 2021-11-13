(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2021 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Read-only polymorphic arrays.

    An immutable version of the [Array] module from the OCaml standard
    library. Used when the values in an array are shared with the underlying C
    library and must not be removed. Holding the values in the array stops
    them from being garbage collected until the array is no longer required.
    This, in turn, prevents the corresponding C values from being freed while
    the array is still reachable. *)

(** An immutable array. *)
type 'a t

(** Copy a mutable array into an immutable one. *)
val from_array : 'a array -> 'a t

(** Copy an immutable array into a mutable one. *)
val to_array : 'a t -> 'a array

(** Return the length of the array. *)
val length : 'a t -> int

(** Return the element at the given index. *)
val get : 'a t -> int -> 'a

(** Create an immutable array of the given length and apply a function to
    set each element. *)
val init : int -> (int -> 'a) -> 'a t

(** Create a new immutable array by joining two existing ones. *)
val append : 'a t -> 'a t -> 'a t

(** Create a new immutable array by concatenating a list of existing ones. *)
val concat : 'a t list -> 'a t

(** Create a new immutable array from a slice of an existing one. *)
val sub : 'a t -> int -> int -> 'a t

(** Copy an immutable array. *)
val copy : 'a t -> 'a t

(** Copy the elements of an immutable array into a list. *)
val to_list : 'a t -> 'a list

(** Copy the elements of a list into a new immutable array. *)
val of_list : 'a list -> 'a t

(** Call a function successively on elements of an immutable array
    starting at index [0]. *)
val iter : ('a -> unit) -> 'a t -> unit

(** Call a function successively on index values and elements of an
    immutable array starting at index [0]. *)
val iteri : (int -> 'a -> unit) -> 'a t -> unit

(** Create a new immutable array by mapping a function across the elements
    of an existing one. *)
val map : ('a -> 'b) -> 'a t -> 'b t

(** Create a new immutable array by mapping a function across the elements,
    and their indexes, of an existing one. *)
val mapi : (int -> 'a -> 'b) -> 'a t -> 'b t

(** Fold a function across the elements of an immutable array from the
    [0]th element upward. *)
val fold_left : ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

(** Fold a function across the elements of an immutable array from the
    last element downward. *)
val fold_right : ('b -> 'a -> 'a) -> 'b t -> 'a -> 'a

(** Call a function successively on paired elements from two immutable arrays
    starting at index [0].

    @raise Invalid_Argument if the arrays are not of the same size. *)
val iter2 : ('a -> 'b -> unit) -> 'a t -> 'b t -> unit

(** Call a function successively on paired elements, and their indexes, from
    two immutable arrays starting at index [0].

    @raise Invalid_Argument if the arrays are not of the same size. *)
val iteri2 : (int -> 'a -> 'b -> unit) -> 'a t -> 'b t -> unit

(** Create a new immutable array by mapping a function across the elements
    of two existing ones.

    @raise Invalid_Argument if the arrays are not of the same size. *)
val map2 : ('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t

(** Fold a function across the elements of two immutable arrays from the
    [0]th elements upward.

    @raise Invalid_Argument if the arrays are not of the same size. *)
val fold_left2 : ('a -> 'b -> 'c -> 'a) -> 'a -> 'b t -> 'c t -> 'a

(** Returns true only if a predicate is true for all elements of an
    immutable array. *)
val for_all : ('a -> bool) -> 'a t -> bool

(** Returns true only if a predicate is true for all paired elements of two
    immutable arrays.

    @raise Invalid_Argument if the arrays are not of the same size. *)
val for_all2 : ('a -> 'b -> bool) -> 'a t -> 'b t -> bool

(** Returns true only if a predicate is true for at least one element of an
    immutable array. *)
val exists : ('a -> bool) -> 'a t -> bool

(** Returns true only if the immutable array contains an element that is
    structurally equal to the given one. *)
val mem : 'a -> 'a t -> bool

(** Returns true only if the immutable array contains an element that is
    (physically) identical to the given one. *)
val memq : 'a -> 'a t -> bool

(** Call a function successively on triples from three immutable arrays
    starting at index [0].

    @raise Invalid_Argument if the arrays are not of the same size. *)
val iter3 : ('a -> 'b -> 'c -> unit) -> 'a t -> 'b t -> 'c t -> unit

(** Call a function successively on triples, and their indexes, from three
    immutable arrays starting at index [0].

    @raise Invalid_Argument if the arrays are not of the same size. *)
val iteri3 : (int -> 'a -> 'b -> 'c -> unit) -> 'a t -> 'b t -> 'c t -> unit

(** Create a new immutable array by mapping a function across the elements
    of three existing ones.

    @raise Invalid_Argument if the arrays are not of the same size. *)
val map3 : ('a -> 'b -> 'c -> 'd) -> 'a t -> 'b t -> 'c t -> 'd t

(** Fold a function across the elements of three immutable arrays from the
    [0]th elements upward.

    @raise Invalid_Argument if the arrays are not of the same size. *)
val fold_left3 : ('a -> 'b -> 'c -> 'd -> 'a) -> 'a -> 'b t -> 'c t -> 'd t -> 'a

