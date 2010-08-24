(* Aug 2010, Timothy Bourke (INRIA) *)

val kind : (float, Bigarray.float64_elt) Bigarray.kind
val layout : Bigarray.c_layout Bigarray.layout
type c_array =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
val empty : c_array

type val_array = c_array
type der_array = c_array
type rootval_array = c_array

val create : int -> c_array
val of_array : float array -> c_array 
val fill : c_array -> float -> unit

val print_results : float -> c_array -> unit

module Roots :
  sig
    type t
    val empty : t
    val create : int -> t
    val print : t -> unit
    val get : t -> int -> bool
    val set : t -> int -> bool -> unit
    val length : t -> int 
    val reset : t -> unit
  end

type lmm =
| Adams
| BDF

type bandrange = {mupper : int; mlower : int}
type sprange = { pretype : int; maxl: int }

type linear_solver =
| Dense
| LapackDense
| Band of bandrange
| LapackBand of bandrange
| Diag
| Spgmr of sprange
| Spbcg of sprange
| Sptfqmr of sprange

type iter =
| Newton of linear_solver
| Functional

(* Solver exceptions *)
exception IllInput
exception TooClose
exception TooMuchWork
exception TooMuchAccuracy
exception ErrFailure
exception ConvergenceFailure
exception LinearInitFailure
exception LinearSetupFailure
exception LinearSolveFailure
exception RhsFuncErr
exception FirstRhsFuncFailure
exception RepeatedRhsFuncErr
exception UnrecoverableRhsFuncErr
exception RootFuncFailure

(* TODO: add optional output functions *)
(* TODO: add optional input functions *)

type session

val no_roots : (int * (float -> val_array -> Roots.t -> int))

exception RecoverableFailure
val init :
    lmm
    -> iter
    -> (float -> val_array -> der_array -> unit)
    -> (int * (float -> val_array -> rootval_array -> unit))
    -> val_array
    -> session

val reinit : session -> float -> val_array -> unit
val set_tolerances : session -> float -> c_array -> unit
val get_roots : session -> Roots.t -> unit

val advance : session -> float -> val_array -> float * bool
val step : session -> float -> val_array -> float * bool
val free : session -> unit

