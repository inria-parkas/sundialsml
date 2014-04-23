(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

type 'a nvector = 'a Nvector.nvector

type ('a, 'adata) atimes = 'adata -> 'a nvector -> 'a nvector -> int

type ('a, 'pdata) psolve = 'pdata -> 'a nvector -> 'a nvector -> bool -> int

external modified_gs' : 'a nvector
                        * Sundials.Realarray2.t
                        * int
                        * int
                        * Sundials.real_array
                        -> unit
    = "c_nvec_modified_gs"

let modified_gs v h k p new_vk_norm = modified_gs' (v, h, k, p, new_vk_norm)

external classical_gs' : 'a nvector
                         * Sundials.Realarray2.t
                         * int
                         * int
                         * Sundials.real_array
                         * 'a nvector
                         * Sundials.real_array
                         -> unit
    = "c_nvec_classical_gs"

let classical_gs v h k p new_vk_norm temp s = 
      classical_gs' (v, h, k, p, new_vk_norm, temp, s)

module SPGMR =
  struct
    type 'a t

    type 'a nvector = 'a Nvector.nvector

    external make  : int -> 'a nvector -> 'a t
        = "c_nvec_spgmr_make"

    external solve' : 'a t
                      * 'adata
                      * 'a nvector
                      * 'a nvector
                      * Spils.preconditioning_type
                      * Spils.gramschmidt_type 
                      * float
                      * int
                      * 'pdata
                      * ('a nvector) option
                      * ('a nvector) option
                      * ('a, 'adata) atimes
                      * (('a, 'pdata) psolve) option
                      * float array
                      -> bool * int * int
        = "c_nvec_spgmr_solve"

    let solve s adata x b pretype gstype delta max_restarts pdata
              s1 s2 atimes psolve res_norm
        = solve' (s, adata, x, b, pretype, gstype, delta, max_restarts, pdata,
                  s1, s2, atimes, psolve, res_norm)

  end

module SPBCG =
  struct
    type 'a t

    type 'a nvector = 'a Nvector.nvector

    external make  : int -> 'a nvector -> 'a t
        = "c_nvec_spbcg_make"

    external solve' : 'a t
                      * 'adata
                      * 'a nvector
                      * 'a nvector
                      * Spils.preconditioning_type
                      * float
                      * 'pdata
                      * ('a nvector) option
                      * ('a nvector) option
                      * ('a, 'adata) atimes
                      * (('a, 'pdata) psolve) option
                      * float array
                      -> bool * int * int
        = "c_nvec_spbcg_solve"
    let solve s adata x b pretype delta pdata sx sb atimes psolve res_norm =
      solve' (s, adata, x, b, pretype, delta, pdata, sx, sb, atimes,
              psolve, res_norm)
 end

module SPTFQMR =
  struct
    
    type 'a t

    type 'a nvector = 'a Nvector.nvector

    external make  : int -> 'a nvector -> 'a t
        = "c_nvec_sptfqmr_make"

    external solve' : 'a t
                      * 'adata
                      * 'a nvector
                      * 'a nvector
                      * Spils.preconditioning_type
                      * float
                      * 'pdata
                      * ('a nvector) option
                      * ('a nvector) option
                      * ('a, 'adata) atimes
                      * (('a, 'pdata) psolve) option
                      * float array
                      -> bool * int * int
        = "c_nvec_sptfqmr_solve"
    let solve s adata x b pretype delta pdata sx sb atimes psolve res_norm =
      solve' (s, adata, x, b, pretype, delta, pdata, sx, sb, atimes,
              psolve, res_norm)

 end

