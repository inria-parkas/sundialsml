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

(** Immutable polymorphic arrays.

    An immutable version of the [Array] module from the OCaml standard
    library. Used when the values in an array are shared with the underlying C
    library and must not be removed. Holding the values in the array stops
    them from being garbage collected until the array is no longer required.
    This, in turn, prevents the corresponding C values from being freed while
    the array is still reachable. *)

type 'a t

val from_array : 'a array -> 'a t

val to_array : 'a t -> 'a array

val length : 'a t -> int

val get : 'a t -> int -> 'a

val init : int -> (int -> 'a) -> 'a t

val append : 'a t -> 'a t -> 'a t

val concat : 'a t list -> 'a t

val sub : 'a t -> int -> int -> 'a t

val copy : 'a t -> 'a t

val to_list : 'a t -> 'a list

val of_list : 'a list -> 'a t

val iter : ('a -> unit) -> 'a t -> unit

val iteri : (int -> 'a -> unit) -> 'a t -> unit

val map : ('a -> 'b) -> 'a t -> 'b t

val mapi : (int -> 'a -> 'b) -> 'a t -> 'b t

val fold_left : ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

val fold_right : ('b -> 'a -> 'a) -> 'b t -> 'a -> 'a

val iter2 : ('a -> 'b -> unit) -> 'a t -> 'b t -> unit

val iteri2 : (int -> 'a -> 'b -> unit) -> 'a t -> 'b t -> unit

val map2 : ('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t

val fold_left2 : ('a -> 'b -> 'c -> 'a) -> 'a -> 'b t -> 'c t -> 'a

val for_all : ('a -> bool) -> 'a t -> bool

val for_all2 : ('a -> 'b -> bool) -> 'a t -> 'b t -> bool

val exists : ('a -> bool) -> 'a t -> bool

val mem : 'a -> 'a t -> bool

val memq : 'a -> 'a t -> bool

val iter3 : ('a -> 'b -> 'c -> unit) -> 'a t -> 'b t -> 'c t -> unit

val iteri3 : (int -> 'a -> 'b -> 'c -> unit) -> 'a t -> 'b t -> 'c t -> unit

val map3 : ('a -> 'b -> 'c -> 'd) -> 'a t -> 'b t -> 'c t -> 'd t

val fold_left3 : ('a -> 'b -> 'c -> 'd -> 'a) -> 'a -> 'b t -> 'c t -> 'd t -> 'a

