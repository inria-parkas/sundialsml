(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

external qr_fact : Sundials.RealArray2.t
                   -> Sundials.RealArray.t
                   -> bool
                   -> unit
    = "c_spils_qr_fact"
 
external qr_sol : Sundials.RealArray2.t
                 -> Sundials.RealArray.t
                 -> Sundials.RealArray.t
                 -> unit
    = "c_spils_qr_sol"

external modified_gs : (('a, 'k) Nvector.t) array
                       -> Sundials.RealArray2.t
                       -> int
                       -> int
                       -> float
    = "c_spils_modified_gs"

external classical_gs' : (('a, 'k) Nvector.t) array
                         * Sundials.RealArray2.t
                         * int
                         * int
                         * ('a, 'k) Nvector.t
                         * Sundials.RealArray.t
                         -> float
    = "c_spils_classical_gs"

let classical_gs v h k p temp s = classical_gs' (v, h, k, p, temp, s)

