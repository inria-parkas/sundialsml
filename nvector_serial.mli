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

(** [make n iv] creates a new serial nvector with [n] elements, each initialized
    to [iv]. *)
val make : int -> float -> t

(** [wrap a] creates a new serial nvector over the elements of [a]. *)
val wrap : Sundials.RealArray.t -> t

module Ops : sig
  val n_vlinearsum    : float -> t -> float -> t -> t -> unit
  val n_vconst        : float -> t -> unit
  val n_vprod         : t -> t -> t -> unit
  val n_vdiv          : t -> t -> t -> unit
  val n_vscale        : float -> t -> t -> unit
  val n_vabs          : t -> t -> unit
  val n_vinv          : t -> t -> unit
  val n_vaddconst     : t -> float -> t -> unit
  val n_vdotprod      : t -> t -> float
  val n_vmaxnorm      : t -> float
  val n_vwrmsnorm     : t -> t -> float
  val n_vwrmsnormmask : t -> t -> t -> float
  val n_vmin          : t -> float
  val n_vwl2norm      : t -> t -> float
  val n_vl1norm       : t -> float
  val n_vcompare      : float -> t -> t -> unit
  val n_vinvtest      : t -> t -> bool
  val n_vconstrmask   : t -> t -> t -> bool
  val n_vminquotient  : t -> t -> float
end

