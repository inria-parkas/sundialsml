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

(** Interface to Serial nvectors that are manipulated by the standard operations
    provided by Sundials.
 
    @cvode <node7#ss:nvec_ser> NVECTOR_SERIAL *)

type kind
type t = (Sundials.RealArray.t, kind) Sundials.nvector
type data = Sundials.RealArray.t

(** [make n iv] creates a new serial nvector with [n] elements, each initialized
    to [iv]. *)
val make : int -> float -> t

(** [wrap a] creates a new serial nvector over the elements of [a]. *)
val wrap : Sundials.RealArray.t -> t

(** Nvector operations on serial nvectors implemented in OCaml. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations on {!Sundials.RealArray}s implemented in OCaml. *)
module DataOps : Nvector.NVECTOR_OPS with type t = Sundials.RealArray.t

