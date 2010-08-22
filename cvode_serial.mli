(* Aug 2010, Timothy Bourke (INRIA) *)

val kind : (float, Bigarray.float64_elt) Bigarray.kind
val layout : Bigarray.c_layout Bigarray.layout
type c_array =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
type int_array =
  (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t

type val_array = c_array
type der_array = c_array
type root_array = c_array

val create : int -> c_array
val of_array : float array -> c_array 
val int_array : int -> int_array 

val print_results : float -> c_array -> unit
val print_roots : int_array -> unit

type session

val no_roots : (int * (float -> val_array -> root_array -> int))

exception RecoverableFailure
val init :
    (float -> val_array -> der_array -> unit)
    -> (int * (float -> val_array -> root_array -> unit))
    -> val_array
    -> session
val set_tolerances : session -> float -> c_array -> unit
val get_roots : session -> int_array -> unit

val advance : session -> float -> val_array -> float * bool
val step : session -> float -> val_array -> float * bool
val free : session -> unit

