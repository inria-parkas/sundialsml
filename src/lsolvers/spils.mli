(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(***********************************************************************)
(* The documentation text is adapted from comments in the Sundials     *)
(* header files by Scott D. Cohen, Alan C. Hindmarsh and Radu Serban   *)
(* at the Center for Applied Scientific Computing, Lawrence Livermore  *)
(* National Laboratory.                                                *)
(***********************************************************************)

(** Scaled Preconditioned Iterative Linear Solvers routines.

    Global constants and general purpose solver routines.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvode <node9#s:spils>  The SPILS Modules *)

(** Performs a QR factorization of a Hessenberg matrix.
    The call [qr_fact h q newjob], where [h] is the [n+1] by [n]
    Hessenberg matrix (stored row-wise), [q] stores the computed Givens
    rotation, and [newjob=false] indicates that the first [n-1] columns of
    [h] have already been factored. The computed Givens rotation has the form
    {% $\begin{bmatrix} c & -s \\ s & c \end{bmatrix}$ %}. It is stored in
    the [2n] elements of [q] as [[|c; s; c; s; ...; c; s|]].

    @raise Matrix.ZeroDiagonalElement Zero found in matrix diagonal *)
val qr_fact : Sundials.RealArray2.t
              -> Sundials.RealArray.t
              -> bool
              -> unit

(** Solve the linear least squares problem. In
    [qr_sol h q b], [h] and [q] are, respectively, the upper triangular
    factor $R$ of the original Hessenberg matrix and [Q] the Givens
    rotations used to factor it—both computed by {!qr_fact}. The function
    computes the [n+1] elements of [b] to solve $Rx = Qb$.

    @raise Matrix.ZeroDiagonalElement Zero found in matrix diagonal *)
val qr_sol : Sundials.RealArray2.t
             -> Sundials.RealArray.t
             -> Sundials.RealArray.t
             -> unit

(** Performs a modified Gram-Schmidt orthogonalization. In
    [modified_gs v h k p],
  - [v] is an array of at least [k + 1] vectors with an L2-norm of 1,
  - [h] is the output [k] by [k] Hessenberg matrix of inner products,
  - [k] specifies the vector in [v] to be orthogonalized against previous ones,
        and,
  - [p] is the number of previous vectors in [v] to orthogonalize against.

  The vector [v[k]] is orthogonalized against the [p] unit vectors at
  [v.{k-1}], [v.{k-2}], ..., [v.{k-p}].
  The matrix [h] must be allocated row-wise so that the [(i,j)]th entry is
  [h.{i}.{j}].
  The inner products are computed, {% $\mathtt{h.\\{}i\mathtt{, k-1\\}} =
    \mathtt{v.\\{}i\mathtt{\\}} \cdot \mathtt{v.\\{k\\}}$ %}, for
  {% $i=\max(0, \mathtt{k}-\mathtt{p})\ldots \mathtt{k}-1$ %}.
  The orthogonalized [v.{k}] is {b not} normalized and is stored over the old
  [v.{k}]. The function returns the Euclidean norm of the orthogonalized
  vector. *)
val modified_gs : (('d, 'k) Nvector.t) array
                 -> Sundials.RealArray2.t
                 -> int
                 -> int
                 -> float

(** Performs a classical Gram-Schmidt orthogonalization. In
    [classical_gs v h k p temp s],
  - [v] is an array of at least [k + 1] vectors with an L2-norm of 1,
  - [h] is the output [k] by [k] Hessenberg matrix of inner products,
  - [k] specifies the vector in [v] to be orthogonalized against previous ones,
        and,
  - [p] is the number of previous vectors in [v] to orthogonalize against.
  - [temp] and [s] are used as workspaces.

  The vector [v[k]] is orthogonalized against the [p] unit vectors at
  [v.{k-1}], [v.{k-2}], ..., [v.{k-p}].
  The matrix [h] must be allocated row-wise so that the [(i,j)]th entry is
  [h.{i}.{j}].
  The inner products are computed, {% $\mathtt{h.\\{}i\mathtt{, k-1\\}} =
    \mathtt{v.\\{}i\mathtt{\\}} \cdot \mathtt{v.\\{k\\}}$ %}, for
  {% $i=\max(0, \mathtt{k}-\mathtt{p})\ldots \mathtt{k}-1$ %}.
  The orthogonalized [v.{k}] is {b not} normalized and is stored over the old
  [v.{k}]. The function returns the Euclidean norm of the orthogonalized
  vector. *)
val classical_gs : (('d, 'k) Nvector.t) array
                  -> Sundials.RealArray2.t
                  -> int
                  -> int
                  -> ('d, 'k) Nvector.t
                  -> Sundials.RealArray.t
                  -> float

