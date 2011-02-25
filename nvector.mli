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

(***********************************************************************)
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Define abstract NVectors using imperative or functional operations

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 @see <CVODE_DOC_ROOT(node7#s:nvector)> CVode:6
 *)

type 'a nvector

module Mutable :
  sig
    type 'a nvector_ops = {
      nvclone           : 'a -> 'a;
      nvdestroy         : ('a -> unit) option;
      nvspace           : ('a -> int * int) option;
      nvlinearsum       : float -> 'a -> float -> 'a -> 'a -> unit;
      nvconst           : float -> 'a -> unit;
      nvprod            : 'a -> 'a -> 'a -> unit;
      nvdiv             : 'a -> 'a -> 'a -> unit;
      nvscale           : float -> 'a -> 'a -> unit;
      nvabs             : 'a -> 'a -> unit;
      nvinv             : 'a -> 'a -> unit;
      nvaddconst        : 'a -> float -> 'a -> unit;
      nvmaxnorm         : 'a -> float;
      nvwrmsnorm        : 'a -> 'a -> float;
      nvmin             : 'a -> float;
      nvdotprod         : ('a -> 'a -> float) option;
      nvcompare         : (float -> 'a -> 'a -> unit) option;
      nvinvtest         : ('a -> 'a -> bool) option;
      nvwl2norm         : ('a -> 'a -> float) option;
      nvl1norm          : ('a -> float) option;
      nvwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
      nvconstrmask      : ('a -> 'a -> 'a -> bool) option;
      nvminquotient     : ('a -> 'a -> float) option;
    }

    val make_nvector    : 'a nvector_ops -> 'a -> 'a nvector
    val nvector_data    : 'a nvector -> 'a

    val add_tracing     : string -> 'a nvector_ops -> 'a nvector_ops
  end

module Immutable :
  sig
    type 'a nvector_ops = {
      nvclone           : 'a -> 'a;
      nvdestroy         : ('a -> unit) option;
      nvspace           : ('a -> int * int) option;
      nvlinearsum       : float -> 'a -> float -> 'a -> 'a;
      nvconst           : float -> 'a;
      nvprod            : 'a -> 'a -> 'a;
      nvdiv             : 'a -> 'a -> 'a;
      nvscale           : float -> 'a -> 'a;
      nvabs             : 'a -> 'a;
      nvinv             : 'a -> 'a;
      nvaddconst        : 'a -> float -> 'a;
      nvmaxnorm         : 'a -> float;
      nvwrmsnorm        : 'a -> 'a -> float;
      nvmin             : 'a -> float;
      nvdotprod         : ('a -> 'a -> float) option;
      nvcompare         : (float -> 'a -> 'a) option;
      nvinvtest         : ('a -> 'a -> bool) option;
      nvwl2norm         : ('a -> 'a -> float) option;
      nvl1norm          : ('a -> float) option;
      nvwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
      nvconstrmask      : ('a -> 'a -> 'a -> bool) option;
      nvminquotient     : ('a -> 'a -> float) option;
    }

    val from_immutable  : 'a nvector_ops -> 'a ref Mutable.nvector_ops
    val make_nvector    : 'a nvector_ops -> 'a ref -> 'a ref nvector
    val nvector_data    : 'a ref nvector -> 'a
  end

