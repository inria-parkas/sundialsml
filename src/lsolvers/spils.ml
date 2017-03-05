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

exception ZeroDiagonalElement of int
exception ConvFailure
exception QRfactFailure
exception PSolveFailure of bool
exception ATimesFailure of bool
exception PSetFailure of bool
exception GSFailure
exception QRSolFailure

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

type 'a atimes = 'a -> 'a -> unit

type 'a psolve = 'a -> 'a -> bool -> unit

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

module SPGMR =
  struct
    type ('d, 'k) memrec
    type ('d, 'k) t = ('d, 'k) memrec * (('d, 'k) Nvector.t -> unit)

    external make' : int -> ('d, 'k) Nvector.t -> ('d, 'k) memrec
        = "c_spils_spgmr_make"

    let make lmax temp = (make' lmax temp, Nvector.check temp)

    external solve' : ('d, 'k) memrec               (*  0 *)
                      * ('d, 'k) Nvector.t          (*  1 *)
                      * ('d, 'k) Nvector.t          (*  2 *)
                      * preconditioning_type        (*  3 *)
                      * gramschmidt_type            (*  4 *)
                      * float                       (*  5 *)
                      * int                         (*  6 *)
                      * (('d, 'k) Nvector.t) option (*  7 *)
                      * (('d, 'k) Nvector.t) option (*  8 *)
                      * 'd atimes                   (*  9 *)
                      * ('d psolve) option          (* 10 *)
                      -> bool * float * int * int
        = "c_spils_spgmr_solve"

    let solve (s, checkvec) ~x ~b ~delta ?maxr:(mr=0) ?s1 ?s2
                ?psolve atimes pretype gstype
        = if Sundials_config.safe then begin
            checkvec x;
            checkvec b;
            (match s1 with None -> () | Some v -> checkvec v);
            (match s2 with None -> () | Some v -> checkvec v)
          end;
          solve' (s, x, b, pretype, gstype, delta, mr, s1, s2,
                  atimes, psolve)
  end

module SPFGMR =
  struct
    type ('d, 'k) memrec
    type ('d, 'k) t = ('d, 'k) memrec * (('d, 'k) Nvector.t -> unit)

    external make' : int -> ('d, 'k) Nvector.t -> ('d, 'k) memrec
        = "c_spils_spfgmr_make"

    let make lmax temp = (make' lmax temp, Nvector.check temp)

    external solve' : ('d, 'k) memrec               (*  0 *)
                      * ('d, 'k) Nvector.t          (*  1 *)
                      * ('d, 'k) Nvector.t          (*  2 *)
                      * preconditioning_type        (*  3 *)
                      * gramschmidt_type            (*  4 *)
                      * float                       (*  5 *)
                      * int                         (*  6 *)
                      * (('d, 'k) Nvector.t) option (*  7 *)
                      * (('d, 'k) Nvector.t) option (*  8 *)
                      * 'd atimes                   (*  9 *)
                      * ('d psolve) option          (* 10 *)
                      * int                         (* 11 *)
                      -> bool * float * int * int
        = "c_spils_spfgmr_solve"

    let solve (s, checkvec) ~x ~b ~delta ?maxr:(mr=0)
                ?maxi:(mi=max_int) ?s1 ?s2 ?psolve atimes pretype gstype
        = if Sundials_config.safe then begin
            checkvec x;
            checkvec b;
            (match s1 with None -> () | Some v -> checkvec v);
            (match s2 with None -> () | Some v -> checkvec v)
          end;
          solve' (s, x, b, pretype, gstype, delta, mr, s1, s2,
                  atimes, psolve, mi)
  end

module SPBCG =
  struct
    type ('d, 'k) memrec
    type ('d, 'k) t = ('d, 'k) memrec * (('d, 'k) Nvector.t -> unit)

    external make'  : int -> ('d, 'k) Nvector.t -> ('d, 'k) memrec
        = "c_spils_spbcg_make"

    let make lmax temp = (make' lmax temp, Nvector.check temp)

    external solve' : ('d, 'k) memrec
                      * ('d, 'k) Nvector.t
                      * ('d, 'k) Nvector.t
                      * preconditioning_type
                      * float
                      * (('d, 'k) Nvector.t) option
                      * (('d, 'k) Nvector.t) option
                      * 'd atimes
                      * ('d psolve) option
                      -> bool * float * int * int
        = "c_spils_spbcg_solve"

    let solve (s, checkvec) ~x ~b ~delta ?sx ?sb ?psolve atimes pretype
        = if Sundials_config.safe then begin
            checkvec x;
            checkvec b;
            (match sx with None -> () | Some v -> checkvec v);
            (match sb with None -> () | Some v -> checkvec v)
          end;
          solve' (s, x, b, pretype, delta, sx, sb, atimes, psolve)
 end

module SPTFQMR =
  struct
    type ('d, 'k) memrec
    type ('d, 'k) t = ('d, 'k) memrec * (('d, 'k) Nvector.t -> unit)

    external make'  : int -> ('d, 'k) Nvector.t -> ('d, 'k) memrec
        = "c_spils_sptfqmr_make"

    let make lmax temp = (make' lmax temp, Nvector.check temp)

    external solve' : ('d, 'k) memrec
                      * ('d, 'k) Nvector.t
                      * ('d, 'k) Nvector.t
                      * preconditioning_type
                      * float
                      * (('d, 'k) Nvector.t) option
                      * (('d, 'k) Nvector.t) option
                      * 'd atimes
                      * ('d psolve) option
                      -> bool * float * int * int
        = "c_spils_sptfqmr_solve"

    let solve (s, checkvec) ~x ~b ~delta ?sx ?sb ?psolve atimes pretype
        = if Sundials_config.safe then begin
            checkvec x;
            checkvec b;
            (match sx with None -> () | Some v -> checkvec v);
            (match sb with None -> () | Some v -> checkvec v)
          end;
          solve' (s, x, b, pretype, delta, sx, sb, atimes, psolve)
 end

module PCG =
  struct
    type ('d, 'k) memrec
    type ('d, 'k) t = ('d, 'k) memrec * (('d, 'k) Nvector.t -> unit)

    external make'  : int -> ('d, 'k) Nvector.t -> ('d, 'k) memrec
        = "c_spils_pcg_make"

    let make lmax temp = (make' lmax temp, Nvector.check temp)

    external solve' : ('d, 'k) memrec               (*  0 *)
                      * ('d, 'k) Nvector.t          (*  1 *)
                      * ('d, 'k) Nvector.t          (*  2 *)
                      * ('d, 'k) Nvector.t          (*  3 *)
                      * preconditioning_type        (*  4 *)
                      * float                       (*  5 *)
                      * 'd atimes                   (*  6 *)
                      * ('d psolve) option          (*  7 *)
                      -> bool * float * int * int
        = "c_spils_pcg_solve"

    let solve (s, checkvec) ~x ~b ~delta ~w ?psolve atimes pretype
        = if Sundials_config.safe then begin
            checkvec x;
            checkvec b;
            checkvec w
          end;
          solve' (s, x, b, w, pretype, delta, atimes, psolve)
  end

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "c_spils_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       spils_exn_index.  *)
    [|ZeroDiagonalElement 0;
      ConvFailure;
      QRfactFailure;
      PSolveFailure false;
      ATimesFailure false;
      PSetFailure false;
      GSFailure;
      QRSolFailure;
    |]

