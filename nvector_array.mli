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

(** A custom nvector based on {!float array}s.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

(** An abstract set of functions for working manipulating nvectors
    where the underlying data structure is an array of [float]s.  *)
module type ARRAY_NVECTOR =
  sig
    type data
    type kind = Nvector_custom.kind

    (** The set of nvector operations on an array. *)
    val array_nvec_ops  : data Nvector_custom.nvector_ops

    (** [make n x] creates an nvector containing an array
        of [n] elements, each of which is equal to [x]. *)
    val make            : int -> float -> data Nvector_custom.t

    (** Lifts an array to an nvector. *)
    val wrap            : data -> data Nvector_custom.t

    (** Returns the array underlying an nvector. *)
    val unwrap          : data Nvector_custom.t -> data

    module Ops : Nvector.NVECTOR_OPS with type t = data Nvector_custom.t
    module DataOps : Nvector.NVECTOR_OPS with type t = data
  end

(** Produce a custom nvector from basic operations on
    an underlying array. *)
module MakeOps : functor (A : sig
      type data
      val get       : data -> int -> float
      val set       : data -> int -> float -> unit
      val fill      : data -> float -> unit
      val make      : int -> float -> data
      val clone     : data -> data
      val length    : data -> int
    end) -> ARRAY_NVECTOR with type data = A.data

(** Nvector on {{:OCAML_DOC_ROOT(Array)} Array}s of [float]s. *)
include ARRAY_NVECTOR with type data = float array
 
