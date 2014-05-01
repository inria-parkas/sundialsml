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

(** Abstract nvector interface to the SPILS routines.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

(** The type of vectors passed to the solver. *)
type 'a nvector = 'a Nvector.nvector

(**
  The type of a function [f v z] that calculates [z = A v] using an internal
  representation of [A]. The vector [v] must not be changed. Results are stored
  in [z]. The {!Sundials.RecoverableFailure} exception can be raised to indicate
  a recoverable failure. Any other exception indicates an unrecoverable failure.
 *)
type 'a atimes = 'a -> 'a -> unit

(**
  The type of a fucntion [f r z lr] that solves the preconditioner equation
  [P z = r] for the vector [z]. If [lr] is true then [P] should be taken as the
  left preconditioner and otherwise as the right preconditioner. The
  {!Sundials.RecoverableFailure} exception can be raised to indicate a
  recoverable failure. Any other exception indicates an unrecoverable failure.
 *)
type 'a psolve = 'a -> 'a -> bool -> unit

(**
  [new_vk_norm = modified_gs v h k p new_vk_norm] performs a modified
  Gram-Schmidt orthogonalization  of [v[k]] against the [p] unit vectors at
  [v.{k-1}], [v.{k-2}], ..., [v.{k-p}]. Its arguments are:
  - [v] an array of [k + 1] vectors assumed to have an L2-norm of 1,
  - [h] is the output [k] by [k] Hessenberg matrix of inner products.
    This matrix must be allocated row-wise so that the [(i,j)]th entry is
    [h.{i}.{j}]. The inner products [(v.{i}, v.{k})], [i=i0], [i0+1], ...,
    [k-1], are stored at [h.{i}.{k-1}] where [i0=MAX(0,k-p)],
  - [k] is the index of the vector in [v] that needs to be
    orthogonalized against previous vectors in [v],
  - [p] is the number of previous vectors in [v] against     
    which [v.{k}] is to be orthogonalized, and,
  The returned value, [new_vk_norm], is the Euclidean norm of the orthogonalized
  vector [v.{k}].
                                                                
  If [(k-p) < 0], then [modified_gs] uses [p=k]. The orthogonalized [v.{k}] is
  not normalized and is stored over the old [v.{k}]. Once the orthogonalization
  has been performed, the Euclidean norm of [v.{k}] is stored in [new_vk_norm].                           
 *)
val modified_gs : ('a nvector) array
                 -> Sundials.Realarray2.t
                 -> int
                 -> int
                 -> float

(**
  [new_vk_norm = classical_gs v h k p new_vk_norm temp s] performs a classical
  Gram-Schmidt orthogonalization  of [v[k]] against the [p] unit vectors at
  [v.{k-1}], [v.{k-2}], ..., [v.{k-p}]. Its arguments are:
  - [v] an array of [k + 1] vectors assumed to have an L2-norm of 1,
  - [h] is the output [k] by [k] Hessenberg matrix of inner products.
    This matrix must be allocated row-wise so that the [(i,j)]th entry is
    [h.{i}.{j}]. The inner products [(v.{i}, v.{k})], [i=i0], [i0+1], ...,
    [k-1], are stored at [h.{i}.{k-1}] where [i0=MAX(0,k-p)],
  - [k] is the index of the vector in [v] that needs to be
    orthogonalized against previous vectors in [v],
  - [p] is the number of previous vectors in [v] against     
    which [v.{k}] is to be orthogonalized,
  - [temp] is used as a workspace, and,
  - [s] is another workspace.
  The returned value, [new_vk_norm], is the Euclidean norm of the orthogonalized
  vector [v.{k}].

  If [(k-p) < 0], then [modifiedGS] uses [p=k]. The orthogonalized [v.{k}] is
  not normalized and is stored over the old [v.{k}]. Once the orthogonalization
  has been performed, the Euclidean norm of [v.{k}] is stored in [new_vk_norm].                           
 *)
val classical_gs : ('a nvector) array
                  -> Sundials.Realarray2.t
                  -> int
                  -> int
                  -> 'a nvector
                  -> Sundials.real_array
                  -> float

(** {3 Scaled Preconditioned GMRES Method }
    @cvode <node9#ss:spgmr> The SPGMR Module *)

(** Implementation of the Scaled Preconditioned Generalized Minimum Residual
   (GMRES) method. *)
module SPGMR :
  sig
    
    (**
     This type represents a solver instance, based on {!Nvector.nvector}s,
     returned from a call to {!make}.
    *)
    type 'a t

    (**
     [make lmax temp] returns a solver session, where [lmax] is the maximum
     Krylov subspace dimension that the solver will be permitted to use, and
     [temp] indirectly specifies the problem size.

     @cvode <node9#ss:spgmr> SpgmrMalloc
     @raise MemoryRequestFailure Memory could not be allocated.
     *)
    val make  : int -> 'a nvector -> 'a t

    (**
     [solved, res_norm, nli, nps = solve s x b pretype gstype delta max_restarts
     s1 s2 atimes psolve res_norm] solves the linear system [Ax = b] using the
     SPGMR iterative method where
      - [s] is a solver session (allocated with {!make}),
      - [x] is the initial guess upon entry, and the solution on return,
      - [b] is the right-hand side vector,
      - [pretype] is the type of preconditioning to use,
      - [gstype] is the type of Gram-Schmidt orthogonalization to use,
      - [delta] is the tolerance on the L2 norm of the scaled, preconditioned
        residual which will satisfy [|| s1 P1_inv (b - Ax) ||_2 <= delta] if
        [solved] is true,
      - [max_restarts] is the maximum number of allowed restarts,
      - [s1] are the optional positive scale factors for [P1-inv b] where
        [P1] is the left preconditioner,
      - [s2] are the optional positive scale factors for [P2 x] where [P2]
        is the right preconditioner,
      - [atimes] multiplies the coefficients [A] by a given vector,
      - [psolve] optionally solves the preconditioner system, and,

      The returned value, [solved], indicates whether the system converged or
      whether it only managed to reduce the norm of the residual, [res_norm] is
      the L2 norm of the scaled preconditioned residual, [|| s1 P1_inv (b - Ax)
      ||_2], [nli] indicates the number of linear iterations performed, and
      [nps] indicates the number of calls made to [psolve].

      Repeated calls can be made to [solve] with varying input arguments, but a
      new session must be created with [make] if the problem size or the maximum
      Krylov dimension changes.

      @cvode <node9#ss:spgmr> SpgmrSolve
      @raise ConvFailure Failed to converge
      @raise QRfactFailure QRfact found singular matrix
      @raise PSolveFailure psolve failed (recoverable or not)
      @raise ATimesFailure atimes failed (recoverable or not)
      @raise PSetFailure pset failed (recoverable or not)
      @raise GSFailure Gram-Schmidt routine failed.
      @raise QRSolFailure QRsol found singular R.
     *)
    val solve : 'a t
                -> 'a nvector
                -> 'a nvector
                -> Spils.preconditioning_type
                -> Spils.gramschmidt_type 
                -> float
                -> int
                -> ('a nvector) option
                -> ('a nvector) option
                -> 'a atimes
                -> ('a psolve) option
                -> bool * float * int * int
  end

(** {3 Scaled Preconditioned Bi-CGStab Method }
    @cvode <node9#ss:spgmr> The SPBCG Module *)

(** Implementation of the Scaled Preconditioned Biconjugate Gradient Stabilized
   (Bi-CGStab) method. *)
module SPBCG :
  sig
    
    (**
     This type represents a solver instance, based on {!Nvector.nvector}s,
     returned from a call to {!make}.
    *)
    type 'a t

    (**
     [make lmax temp] returns a solver session, where [lmax] is the maximum
     Krylov subspace dimension that the solver will be permitted to use, and
     [temp] indirectly specifies the problem size.

     @cvode <node9#ss:spbcg> SpbcgMalloc
     @raise MemoryRequestFailure Memory could not be allocated.
     *)
    val make  : int -> 'a nvector -> 'a t

    (**
     [solved, res_norm, nli, nps = solve s x b pretype delta sx sb atimes
     psolve] solves the linear system [Ax = b] using the SPBCG iterative method
     where
      - [s] is a solver session (allocated with {!make}),
      - [x] is the initial guess upon entry, and the solution on return,
      - [b] is the right-hand side vector,
      - [pretype] is the type of preconditioning to use,
      - [delta] is the tolerance on the L2 norm of the scaled, preconditioned
        residual which will satisfy [||sb*P1_inv*(b-Ax)||_L2 <= delta] if
        [solved] is true,
      - [sx] are the optional positive scaling factors for [x],
      - [sb] are the optional positive scaling factors for [b],
      - [atimes] multiplies the coefficients [A] by a given vector,
      - [psolve] optionally solves the preconditioner system, and,

      The returned value, [solved], indicates whether the system converged or
      whether it only managed to reduce the norm of the residual, [res_norm] is
      used for returning the L2 norm of the scaled preconditioned residual,
      [||sb*P1_inv*(b-Ax)||_L2], [nli] indicates the number of linear iterations
      performed, and [nps] indicates the number of calls made to [psolve].

      Repeated calls can be made to [solve] with varying input arguments, but a
      new session must be created with [make] if the problem size or the maximum
      Krylov dimension changes.

      @cvode <node9#ss:spbcg> SpbcgSolve
      @raise ConvFailure Failed to converge
      @raise PSolveFailure psolve failed (recoverable or not)
      @raise ATimesFailure atimes failed (recoverable or not)
      @raise PSetFailure pset failed (recoverable or not)
     *)
    val solve : 'a t
                -> 'a nvector
                -> 'a nvector
                -> Spils.preconditioning_type
                -> float
                -> ('a nvector) option
                -> ('a nvector) option
                -> 'a atimes
                -> ('a psolve) option
                -> bool * float * int * int

 end

(** {3 Scaled Preconditioned TFQMR Method }
    @cvode <node9#ss:sptfqmr> The SPTFQMR Module *)

(** Implementation of the Scaled Preconditioned Transpose-Free Quasi-Minimal
    Residual (SPTFQMR) method *)
module SPTFQMR :
  sig
    
    (**
     This type represents a solver instance, based on {!Nvector.nvector}s,
     returned from a call to {!make}.
    *)
    type 'a t

    (**
     [make lmax temp] returns a solver session, where [lmax] is the maximum
     Krylov subspace dimension that the solver will be permitted to use, and
     [temp] indirectly specifies the problem size.

     @cvode <node9#ss:sptfqmr> SptfqmrMalloc
     @raise MemoryRequestFailure Memory could not be allocated.
     *)
    val make  : int -> 'a nvector -> 'a t

    (**
     [solved, res_norm, nli, nps = solve s x b pretype delta sx sb atimes
     psolve] solves the linear system [Ax = b] using the SPTFQMR iterative
     method where
      - [s] is a solver session (allocated with {!make}),
      - [x] is the initial guess upon entry, and the solution on return,
      - [b] is the right-hand side vector,
      - [pretype] is the type of preconditioning to use,
      - [delta] is the tolerance on the L2 norm of the scaled, preconditioned
        residual which will satisfy [||sb*P1_inv*(b-Ax)||_L2 <= delta] if
        [solved] is true,
      - [sx] are the optional positive scaling factors for [x],
      - [sb] are the optional positive scaling factors for [b],
      - [atimes] multiplies the coefficients [A] by a given vector,
      - [psolve] optionally solves the preconditioner system, and,

      The value returned by this function, [solved], indicates whether the
      system converged or whether it only managed to reduce the norm of the
      residual, [res_norm] is used for returning the L2 norm of the scaled
      preconditioned residual, [||sb*P1_inv*(b-Ax)||_L2], [nli] indicates the
      number of linear iterations performed, and [nps] indicates the number of
      calls made to [psolve].

      Repeated calls can be made to [solve] with varying input arguments, but a
      new session must be created with [make] if the problem size or the maximum
      Krylov dimension changes.

      @cvode <node9#ss:sptfqmr> SptfqmrSolve
      @raise ConvFailure Failed to converge
      @raise PSolveFailure psolve failed (recoverable or not)
      @raise ATimesFailure atimes failed (recoverable or not)
      @raise PSetFailure pset failed (recoverable or not)
     *)
    val solve : 'a t
                -> 'a nvector
                -> 'a nvector
                -> Spils.preconditioning_type
                -> float
                -> ('a nvector) option
                -> ('a nvector) option
                -> 'a atimes
                -> ('a psolve) option
                -> bool * float * int * int

 end

