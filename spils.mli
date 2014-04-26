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

(***********************************************************************)
(* The documentation text is adapted from comments in the Sundials     *)
(* header files by Scott D. Cohen, Alan C. Hindmarsh and Radu Serban   *)
(* at the Center for Applied Scientific Computing, Lawrence Livermore  *)
(* National Laboratory.                                                *)
(***********************************************************************)

(**
  Scaled Preconditioned Iterative Linear Solvers (SPILS) routines.
  @cvode <node9#s:spils>  The SPILS Modules
 *)

(**
 Constants representing the types of Gram-Schmidt orthogonalization possible for
 the SPGMR linear solver ({!Spils_nvector.SPGMR.solve}/{!Spils_serial.SPGMR.solve}).
 @cvode <node9.html#ss:spgmr> ModifiedGS/ClassicalGS
 *)
type gramschmidt_type =
  | ModifiedGS
        (** Modified Gram-Schmidt orthogonalization (MODIFIED_GS) *)
  | ClassicalGS
        (** Classical Gram Schmidt orthogonalization (CLASSICAL_GS) *)

(**
 Type of preconditioning for Krylov solvers.
 @cvode <node3#s:preconditioning> Preconditioning
 @cvode <node5#sss:lin_solv_init> CVSpgmr/CVSpbcg/CVSptfqrm
 *)
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

(**
  [r = qr_fact n h q newjob] performs a QR factorization of the Hessenberg
  matrix [h], where
  - [n] is the problem size,
  - [h] is the [n+1] by [n] Hessenberg matrix (stored row-wise) to be factored,
  - [q] is an array of length [2*n] containing the Givens rotations    computed
    by this function. A Givens rotation has the form [| c -s; s c |]. The
    components of the Givens rotations are stored in [q] as
    [(c, s, c, s, ..., c, s)], and,
  - if [newjob] is true then a new QR factorization is performed, otherwise it
    is assumed that the first [n-1] columns of [h] have already been factored,
    and only the last column needs to be updated.

  The result, [r], is 0 if successful. If a zero is encountered on the diagonal
  of the triangular factor [R], then QRfact returns the equation number of the
  zero entry, where the equations are numbered from 1, not 0. If {!qr_sol} is
  subsequently called in this situation, it will return an error because it
  could not divide by the zero diagonal entry.                             
 *)
val qr_fact : int
              -> Sundials.Realarray2.t
              -> Sundials.real_array
              -> bool
              -> int
 
(**
  [r = qr_sol n h q b] solves the linear least squares problem
  [min (b - h*x, b - h*x), x in R^n] where
  - [n] is the problem size,
  - [h] is computed by {!qr_fact} containing the upper triangular factor
    [R] of the original Hessenberg matrix,
  - [q] is the array computed by {!qr_fact} containing the Givens rotations used
    to factor [h], and,
  - [b] is the [n+1]-vector which, on successful return, will contain the
    solution [x] of the least squares problem.

  The result, [r], is 0 if successful. Otherwise, a zero was encountered on the
  diagonal of the triangular factor [R]. In this case, QRsol returns the
  equation number (numbered from 1, not 0) of the zero entry.
 *)
val qr_sol : int
             -> Sundials.Realarray2.t
             -> Sundials.real_array
             -> Sundials.real_array
             -> int

