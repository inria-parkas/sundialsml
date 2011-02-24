(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*       Timothy Bourke (INRIA Rennes) and Marc Pouzet (LIENS)         *)
(*                                                                     *)
(*  Copyright 2011 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

(** Generic definitions for the Sundials suite

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 *)

val extra_time_precision : bool ref
val print_time : string * string -> float -> unit

val format_float : string -> float -> string
val floata : float -> string

val big_real : float
val unit_roundoff : float

type real_array =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

val make_real_array : int -> real_array

type real_array2 =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

val make_real_array2 : int -> int -> real_array2

type int_array =
  (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t

val make_int_array  : int -> int_array

module Carray :
  sig
    type t = real_array

    val empty : t
    val create : int -> t
    val of_array : float array -> t
    val fill : t -> float -> unit
    val length : t -> int

    val print_with_time : float -> t -> unit
    val print_with_time' : float -> t -> unit

    val app : (float -> unit) -> t -> unit
    val appi : (int -> float -> unit) -> t -> unit

    val map : (float -> float) -> t -> unit
    val mapi : (int -> float -> float) -> t -> unit

    val clamp : float -> t -> unit
  end

module Roots :
  sig
    type t = int_array
    type val_array = Carray.t

    val empty : t
    val create : int -> t
    val length : t -> int

    val print : t -> unit
    val print' : t -> unit

    val get : t -> int -> bool
    val get' : t -> int -> int

    val set : t -> int -> bool -> unit

    val reset : t -> unit
    val exists : t -> bool

    val app : (bool -> unit) -> t -> unit
    val appi : (int -> bool -> unit) -> t -> unit
  end

