(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*       Timothy Bourke (INRIA Rennes) and Marc Pouzet (LIENS)         *)
(*                                                                     *)
(*  Copyright 2011 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

(** NVectors for standard arrays and one-dimensional bigarrays

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 @see <OCAML_DOC_ROOT(Array)> Array
 @see <OCAML_DOC_ROOT(Bigarray.Array1)> Bigarray.Array1
 *)

module type ARRAY_NVECTOR =
  sig
    type t

    val array_nvec_ops  : t Nvector.Mutable.nvector_ops
    val make            : int -> float -> t Nvector.nvector
    val wrap            : t -> t Nvector.nvector
    val data            : t Nvector.nvector -> t
  end

include ARRAY_NVECTOR with type t = float array

module Bigarray :
  ARRAY_NVECTOR
  with
    type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

