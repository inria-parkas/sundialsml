(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

exception LinearSolverInUse = Lsolver_impl.LinearSolverInUse
exception UnrecoverableFailure of bool
exception MatrixNotSquare
exception MatrixVectorMismatch
exception InsufficientStorageUpperBandwidth
exception ATimesFailure of bool
exception PSetFailure of bool
exception PSolveFailure of bool
exception GSFailure
exception QRSolFailure
exception ResReduced
exception ConvFailure
exception QRfactFailure
exception LUfactFailure
exception PackageFailure of bool
exception IllegalPrecType
exception InternalFailure of (string * int)

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "ml_lsolver_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       lsolver_exn_index.  *)
    [|UnrecoverableFailure false;
      MatrixNotSquare;
      MatrixVectorMismatch;
      InsufficientStorageUpperBandwidth;
      Invalid_argument ""; (* Standard OCaml exception *)
      ATimesFailure false;
      PSetFailure false;
      PSolveFailure false;
      GSFailure;
      QRSolFailure;
      ResReduced;
      ConvFailure;
      QRfactFailure;
      LUfactFailure;
      PackageFailure false;
      IllegalPrecType;
      InternalFailure ("", 0);
    |]

