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

(* Types underlying Sls.  *)

(* This module defines the Sls matrix type which is manipulated (abstractly)
   via the Sls module. It is declared here so that it can also be used
   (concretely) from the <solver>_impl modules (in callbacks from the linear
   solver data type). To ensure that this type will be opaque outside of
   Sundials/ML, we simply do not install the Sls_impl.cmi file. *)

type int_array =
  (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t

type real_array = Sundials.RealArray.t

type slsmat

type sformat =
  | CSC_MAT
  | CSR_MAT

type 'f t = {
  idxptrs : int_array;
  idxvals : int_array;
  data    : real_array;
  slsmat  : slsmat;
  sformat : sformat;
  mutable valid : bool;
}

let invalidate v = v.valid <- false

