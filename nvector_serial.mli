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

(** The standard serial nvectors of Sundials.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @cvode <node7#ss:nvec_ser> NVECTOR_SERIAL *)

(** Serial nvectors are based on
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarrays} of floats. *)
type data = Sundials.RealArray.t

(** Represents the internal layout of a serial nvector. *)
type kind

(** The type of serial nvectors. *)
type t = (data, kind) Nvector.t

(** [make n iv] creates a new serial nvector with [n] elements, each initialized
    to [iv]. *)
val make : int -> float -> t

(** [wrap a] creates a new serial nvector over the elements of [a]. *)
val wrap : Sundials.RealArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> Sundials.RealArray.t

(** Nvector operations on serial nvectors implemented in OCaml. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations on {!Sundials.RealArray}s implemented in OCaml. *)
module DataOps : Nvector.NVECTOR_OPS with type t = Sundials.RealArray.t

