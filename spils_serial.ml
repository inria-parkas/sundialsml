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

include Spils (* Vital for registering Spils exceptions. *)

type real_array = Sundials.RealArray.t

type atimes = real_array -> real_array -> unit

type psolve = real_array -> real_array -> bool -> unit

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("c_ba_spils_RecoverableFailure", Sundials.RecoverableFailure);
  ]

external modified_gs : real_array array
                       -> Sundials.RealArray2.t
                       -> int
                       -> int
                       -> float
    = "c_ba_spils_modified_gs"

external classical_gs' : real_array array
                         * Sundials.RealArray2.t
                         * int
                         * int
                         * real_array
                         * real_array
                         -> float
    = "c_ba_spils_classical_gs"

let classical_gs v h k p temp s = classical_gs' (v, h, k, p, temp, s)

module SPGMR =
  struct
    type t

    external make  : int -> real_array -> t
        = "c_ba_spils_spgmr_make"

    external solve' : t                               (*  0 *)
                      * real_array                    (*  1 *)
                      * real_array                    (*  2 *)
                      * Spils.preconditioning_type    (*  3 *)
                      * Spils.gramschmidt_type        (*  4 *)
                      * float                         (*  5 *)
                      * int                           (*  6 *)
                      * real_array option             (*  7 *)
                      * real_array option             (*  8 *)
                      * atimes                        (*  9 *)
                      * psolve option                 (* 10 *)
                      -> bool * float * int * int
        = "c_ba_spils_spgmr_solve"

    let solve s x b pretype gstype delta max_restarts s1 s2 atimes psolve
        = solve' (s, x, b, pretype, gstype, delta, max_restarts, s1, s2,
                  atimes, psolve)

  end

module SPBCG =
  struct
    type t

    external make  : int -> real_array -> t
        = "c_ba_spils_spbcg_make"

    external solve' : t                               (* 0 *)
                      * real_array                    (* 1 *)
                      * real_array                    (* 2 *)
                      * Spils.preconditioning_type    (* 3 *)
                      * float                         (* 4 *)
                      * real_array option             (* 5 *)
                      * real_array option             (* 6 *)
                      * atimes                        (* 7 *)
                      * psolve option                 (* 8 *)
                      -> bool * float * int * int
        = "c_ba_spils_spbcg_solve"
    let solve s x b pretype delta sx sb atimes psolve =
      solve' (s, x, b, pretype, delta, sx, sb, atimes, psolve)
 end

module SPTFQMR =
  struct
    type t

    external make  : int -> real_array -> t
        = "c_ba_spils_sptfqmr_make"

    external solve' : t                               (* 0 *)
                      * real_array                    (* 1 *)
                      * real_array                    (* 2 *)
                      * Spils.preconditioning_type    (* 3 *)
                      * float                         (* 4 *)
                      * real_array option             (* 5 *)
                      * real_array option             (* 6 *)
                      * atimes                        (* 7 *)
                      * psolve option                 (* 8 *)
                      -> bool * float * int * int
        = "c_ba_spils_sptfqmr_solve"
    let solve s x b pretype delta sx sb atimes psolve =
      solve' (s, x, b, pretype, delta, sx, sb, atimes, psolve)
 end

