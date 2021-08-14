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

(** The standard many-vector nvectors of Sundials.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @nocvode <node> NVECTOR_MANYVECTOR
    @since 5.0.0 *)

open Sundials

(** The data in underlying nvectors is exposed as an array of wrapped values
    paired with the sum of their lengths. *)
type data = Nvector.any ROArray.t * int

(** Represents the internal layout of a many-vector nvector. *)
type kind

(** The type of many-vector nvectors. *)
type t = (data, kind) Nvector.t

(** A generic data wrapper for {!data}. *)
type Nvector.gdata += Many of data

(** Creates a many-vector nvector from an array of generic nvectors.

    @since 5.0.0 *)
val wrap : Nvector.any ROArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> data

(** Returns the sum of the lengths of the component nvectors. *)
val length : t -> int

(** Returns the number of subectors in the array. *)
val num_subvectors : t -> int

(** Underlying nvector operations on many-vector nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

(** Nvector operations implemented in OCaml on many-vector payloads. *)
module DataOps : Nvector.NVECTOR_OPS with type t = data

(** {2:genvec Generic nvector interface}

    Create many-vector nvectors using the generic nvector interface where the
    payload is wrapped with the {{!Nvector.gdata}Many} constructor. *)
module Any : sig

  (** Creates a generic nvector from an array of generic nvectors.

      @since 5.0.0 *)
  val wrap : Nvector.any ROArray.t -> Nvector.any

end

