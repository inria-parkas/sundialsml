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

(** A custom nvector based on float arrays.

  @version VERSION()
  @author Timothy Bourke (Inria/ENS)
  @author Jun Inoue (Inria/ENS)
  @author Marc Pouzet (UPMC/ENS/Inria)
 *)

(** An abstract set of functions for working manipulating nvectors
    where the underlying data structure is an array of [float]s.  *)
module type ARRAY_NVECTOR =
  sig
    (** Type of the underlying array. *)
    type data

    (** Array nvectors are custom nvectors. *)
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

    (** Standard operations over array nvectors. *)
    module Ops : Nvector.NVECTOR_OPS with type t = data Nvector_custom.t

    (** Standard operations over the underlying array. *)
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

(** Nvector on {{:OCAML_DOC_ROOT(Array.html)} Array}s of [float]s. *)
include ARRAY_NVECTOR with type data = float array

