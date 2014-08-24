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

(** Nvectors for standard arrays and one-dimensional bigarrays.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
  @see <OCAML_DOC_ROOT(Array)> Array
  @see <OCAML_DOC_ROOT(Bigarray.Array1)> Bigarray.Array1
 *)

(** An abstract set of functions for working manipulating nvectors
    where the underlying data structure is an array of [float]s.  *)
module type ARRAY_NVECTOR =
  sig
    include Nvector.NVECTOR_OPS

    (** The set of nvector operations on an array. *)
    val array_nvec_ops  : t Nvector_custom.nvector_ops

    (** [make n x] creates an nvector containing an array
        of [n] elements, each of which is equal to [x]. *)
    val make            : int -> float -> t Nvector_custom.t

    (** Lifts an array to an nvector. *)
    val wrap            : t -> t Nvector_custom.t
  end

module MakeOps : functor (A : sig
      type data
      val get       : data -> int -> float
      val set       : data -> int -> float -> unit
      val fill      : data -> float -> unit
      val make      : int -> float -> data
      val clone     : data -> data
      val length    : data -> int
      val fold_left : ('a -> float -> 'a) -> 'a -> data -> 'a
    end) -> ARRAY_NVECTOR with type t = A.data

(** Nvector on {{:OCAML_DOC_ROOT(Array)} Array}s of [float]s. *)
include ARRAY_NVECTOR with type t = float array
 
