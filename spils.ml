(***********************************************************************)
(*                                                                     *)
(*               OCaml interface to (serial) Sundials                  *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

type gramschmidt_type =
  | ModifiedGS
        (** Modified Gram-Schmidt orthogonalization (MODIFIED_GS) *)
  | ClassicalGS
        (** Classical Gram Schmidt orthogonalization (CLASSICAL_GS) *)

type preconditioning_type =
  | PrecNone
  | PrecLeft
  | PrecRight
  | PrecBoth

exception MemoryRequestFailure

exception ConvFailure
exception QRfactFailure
exception PSolveFailure of bool
exception ATimesFailure of bool
exception PSetFailure of bool
exception GSFailure
exception QRSolFailure

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("spils_MemoryRequestFailure",      MemoryRequestFailure);
    ("spils_ConvFailure",               ConvFailure);
    ("spils_QRfactFailure",             QRfactFailure);
    ("spils_PSolveFailure",             PSolveFailure false);
    ("spils_ATimesFailure",             ATimesFailure false);
    ("spils_PSetFailure",               PSetFailure false);
    ("spils_GSFailure",                 GSFailure);
    ("spils_QRSolFailure",              QRSolFailure);
  ]

external qr_fact : int
                   -> Sundials.RealArray2.t
                   -> Sundials.RealArray.t
                   -> bool
                   -> int
    = "c_spils_qr_fact"
 
external qr_sol : int
                 -> Sundials.RealArray2.t
                 -> Sundials.RealArray.t
                 -> Sundials.RealArray.t
                 -> int
    = "c_spils_qr_sol"

