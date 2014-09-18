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

type gramschmidt_type =
  | ModifiedGS
        (** Modified Gram-Schmidt orthogonalization (MODIFIED_GS) *)
  | ClassicalGS
        (** Classical Gram Schmidt orthogonalization (CLASSICAL_GS) *)

type preconditioning_type =
  | PrecTypeNone
  | PrecTypeLeft
  | PrecTypeRight
  | PrecTypeBoth

exception MemoryRequestFailure

exception ConvFailure
exception QRfactFailure
exception PSolveFailure of bool
exception ATimesFailure of bool
exception PSetFailure of bool
exception GSFailure
exception QRSolFailure

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

type ('a, 'k) nvector = ('a, 'k) Sundials.nvector

type 'a atimes = 'a -> 'a -> unit

type 'a psolve = 'a -> 'a -> bool -> unit

external modified_gs : (('a, 'k) nvector) array
                       -> Sundials.RealArray2.t
                       -> int
                       -> int
                       -> float
    = "c_spils_modified_gs"

external classical_gs' : (('a, 'k) nvector) array
                         * Sundials.RealArray2.t
                         * int
                         * int
                         * ('a, 'k) nvector
                         * Sundials.RealArray.t
                         -> float
    = "c_spils_classical_gs"

let classical_gs v h k p temp s = classical_gs' (v, h, k, p, temp, s)

module SPGMR =
  struct
    type 'a t

    external make  : int -> ('a, 'k) nvector -> 'a t
        = "c_spils_spgmr_make"

    external solve' : 'a t                            (*  0 *)
                      * ('a, 'k) nvector                    (*  1 *)
                      * ('a, 'k) nvector                    (*  2 *)
                      * preconditioning_type    (*  3 *)
                      * gramschmidt_type        (*  4 *)
                      * float                         (*  5 *)
                      * int                           (*  6 *)
                      * (('a, 'k) nvector) option           (*  7 *)
                      * (('a, 'k) nvector) option           (*  8 *)
                      * 'a atimes                     (*  9 *)
                      * ('a psolve) option            (* 10 *)
                      -> bool * float * int * int
        = "c_spils_spgmr_solve"

    let solve s x b pretype gstype delta max_restarts s1 s2 atimes psolve
        = solve' (s, x, b, pretype, gstype, delta, max_restarts, s1, s2,
                  atimes, psolve)

  end

module SPBCG =
  struct
    type 'a t

    external make  : int -> ('a, 'k) nvector -> 'a t
        = "c_spils_spbcg_make"

    external solve' : 'a t
                      * ('a, 'k) nvector
                      * ('a, 'k) nvector
                      * preconditioning_type
                      * float
                      * (('a, 'k) nvector) option
                      * (('a, 'k) nvector) option
                      * 'a atimes
                      * ('a psolve) option
                      -> bool * float * int * int
        = "c_spils_spbcg_solve"
    let solve s x b pretype delta sx sb atimes psolve =
      solve' (s, x, b, pretype, delta, sx, sb, atimes, psolve)
 end

module SPTFQMR =
  struct
    
    type 'a t

    external make  : int -> ('a, 'k) nvector -> 'a t
        = "c_spils_sptfqmr_make"

    external solve' : 'a t
                      * ('a, 'k) nvector
                      * ('a, 'k) nvector
                      * preconditioning_type
                      * float
                      * (('a, 'k) nvector) option
                      * (('a, 'k) nvector) option
                      * 'a atimes
                      * ('a psolve) option
                      -> bool * float * int * int
        = "c_spils_sptfqmr_solve"
    let solve s x b pretype delta sx sb atimes psolve =
      solve' (s, x, b, pretype, delta, sx, sb, atimes, psolve)

 end

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "c_spils_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       spils_exn_index.  *)
    [|MemoryRequestFailure;
      ConvFailure;
      QRfactFailure;
      PSolveFailure false;
      ATimesFailure false;
      PSetFailure false;
      GSFailure;
      QRSolFailure;
    |]
