(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(* Types shared between Dls and Sls.  *)

(* This module's purpose is just to define types shared between Dls and Sls.
   The implementations of the DenseMatrix.t and BandMatrix.t types must be
   transparent to both modules, but opaque outside the sundialsml library.

   To satisfy this requirement, we define the type in a separate module
   - this module - that exports everything, then omit Dls_impl.cmi during
   installation. This way, the compiled library doesn't divulge
   implementation details.
 *)

module DenseTypes = struct
  type data = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

  (* Must correspond with dls_ml.h:dls_densematrix_index *)
  type t = {
    payload : data;
    dlsmat  : Obj.t;
    mutable valid : bool;
  }
end

module BandTypes = struct
  type data = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

  (* Must correspond with dls_ml.h:dls_bandmatrix_index *)
  type t = {
    payload : data;
    dlsmat  : Obj.t;
    ismu    : int;
    mutable valid : bool;
  }
end

