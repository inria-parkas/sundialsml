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

type cnvec
type ('data, 'kind) nvector = 'data * cnvec * ('data -> bool)
and ('data, 'kind) t = ('data, 'kind) nvector

let unwrap (payload, _, _) = payload

exception IncompatibleNvector

let check (_, _, checkfn) = (function (payload, _, _) ->
                              if not (checkfn payload) then
                                raise IncompatibleNvector)

module type NVECTOR_OPS =
  sig
    type t

    val n_vclone        : t -> t
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
    val n_vmin          : t -> float
    val n_vcompare      : float -> t -> t -> unit
    val n_vinvtest      : t -> t -> bool

    val n_vwl2norm      : t -> t -> float
    val n_vl1norm       : t -> float
    val n_vwrmsnormmask : t -> t -> t -> float
    val n_vconstrmask   : t -> t -> t -> bool
    val n_vminquotient  : t -> t -> float
  end

module type NVECTOR =
  sig
    type kind
    type data
    type t = (data, kind) nvector
    val wrap : data -> t
    module Ops : NVECTOR_OPS with type t = t
    module DataOps : NVECTOR_OPS with type t = data
  end

