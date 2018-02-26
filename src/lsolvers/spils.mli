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

(** {2 Types} *)

(** The type of Gram-Schmidt orthogonalization in SPGMR linear solvers.

    @cvode <node9#ss:spgmr> ModifiedGS/ClassicalGS *)
type gramschmidt_type =
  | ModifiedGS   (** Modified Gram-Schmidt orthogonalization
                     {cconst MODIFIED_GS} *)
  | ClassicalGS  (** Classical Gram Schmidt orthogonalization
                     {cconst CLASSICAL_GS} *)

(** The type of preconditioning in Krylov solvers.

    @cvode <node3#s:preconditioning> Preconditioning
    @cvode <node5#sss:lin_solv_init> CVSpgmr/CVSpbcg/CVSptfqrm *)
type preconditioning_type =
  | PrecNone    (** No preconditioning *)
  | PrecLeft    (** {% $(P^{-1}A)x = P^{-1}b$ %} *)
  | PrecRight   (** {% $(AP^{-1})Px = b$ %} *)
  | PrecBoth    (** {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)

(** {2 Exceptions} *)

(** Raised by {!qr_fact} and {!qr_sol} on a zero diagonal element during
    factorization. The argument gives the equation number (from 1). *)
exception ZeroDiagonalElement of int

(** Raised when a solver fails to converge.
    {cconst *_CONVFAIL} *)
exception ConvFailure

(** Raised when QR factorization yields a singular matrix.
    {cconst *_SPTFQMR_QRFACT_FAIL} *)
exception QRfactFailure

(** Raised when a preconditioner solver fails. The argument is [true] for a
    recoverable failure and [false] for an unrecoverable one.
    {cconst *_PSOLVE_FAIL_REC/_UNREC} *)
exception PSolveFailure of bool

(** Raised when an atimes function fails. The argument is [true] for a
    recoverable failure and [false] for an unrecoverable one.
    {cconst *_ATIMES_FAIL_REC/_UNREC} *)
exception ATimesFailure of bool

(** Raised when a preconditioner setup routine fails. The argument is [true]
    for a recoverable failure and [false] for an unrecoverable one.
    {cconst *_PSET_FAIL_REC/_UNREC} *)
exception PSetFailure of bool

(** Raised when a Gram-Schmidt routine fails. {cconst *_GS_FAIL} *)
exception GSFailure

(** Raised QR solution finds a singular result. {cconst *_QRSOL_FAIL} *)
exception QRSolFailure

(** {2 Basic routines} *)

(** Performs a QR factorization of a Hessenberg matrix.
    The call [qr_fact h q newjob], where [h] is the [n+1] by [n]
    Hessenberg matrix (stored row-wise), [q] stores the computed Givens
    rotation, and [newjob=false] indicates that the first [n-1] columns of
    [h] have already been factored. The computed Givens rotation has the form
    {% $\begin{bmatrix} c & -s \\ s & c \end{bmatrix}$ %}. It is stored in
    the [2n] elements of [q] as [[|c; s; c; s; ...; c; s|]].

    @raise ZeroDiagonalElement Zero found in matrix diagonal *)
val qr_fact : Sundials.RealArray2.t
              -> Sundials.RealArray.t
              -> bool
              -> unit

(** Solve the linear least squares problem. In
    [qr_sol h q b], [h] and [q] are, respectively, the upper triangular
    factor $R$ of the original Hessenberg matrix and [Q] the Givens
    rotations used to factor it—both computed by {!qr_fact}. The function
    computes the [n+1] elements of [b] to solve $Rx = Qb$.

    @raise ZeroDiagonalElement Zero found in matrix diagonal *)
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

(** {2 Solvers} *)

(** Functions [f v z] that calculate [z = A v] using an internal
    representation of [A]. Results are stored in [z], [v] must not be changed.
    Raise {!Sundials.RecoverableFailure} to indicate a recoverable failure,
    any other exception indicates an unrecoverable failure. *)
type 'd atimes = 'd -> 'd -> unit

(** Functions [f r z tol lr] that solve the preconditioner equation [P z = r]
    for the vector [z] such that
    $\left\lVert Pz - r \right\rVert_\mathrm{wrms} < \mathit{tol}$.
    If [lr = true] then [P] is used as a left preconditioner and otherwise as
    a right preconditioner.
    Raise {!Sundials.RecoverableFailure} to indicate a recoverable failure,
    any other exception indicates an unrecoverable failure.

    In versions of Sundials prior to 3.0.0, the [tol] parameter is always [0.0]
    and should be ignored.
*)
type 'd psolve = 'd -> 'd -> float -> bool -> unit

(** The Scaled Preconditioned Generalized Minimum Residual (GMRES) method. *)
module SPGMR :
  sig
    
    (** An instance of the SPGMR solver.

        @cvode <node9#ss:spgmr> The SPGMR Module *)
    type ('d, 'k) t

    (** [make lmax temp] returns a solver session. [lmax] is the maximum
        Krylov subspace dimension to use, and [temp] sets the problem size.

        @cvode <node9#ss:spgmr> SpgmrMalloc *)
    val make  : int -> ('d, 'k) Nvector.t -> ('d, 'k) t

    (** Solves the linear system [Ax = b] using the SPGMR iterative method.
        The [atimes] function computes the matrix vector product [Ax], the
        other arguments are described below. The function returns a tuple
        [(solved, res_norm, nli, nps)] where [solved] indicates whether the
        system converged, [res_norm] is the L2 norm of the scaled
        preconditioned residual
        {% $\lVert s_1 P_1^{-1} (b - Ax) \rVert_{L2}$ %}, and,
        [nli] and [nps] count, respectively, linear iterations performed and
        calls to [psolve]. Repeated calls can be made to [solve] with
        varying input arguments, but a new session must be created if either
        the problem size or the maximum Krylov dimension change.

        @cvode <node9#ss:spgmr> SpgmrSolve
        @param x initial guess on entry; result on return
        @param b right-hand side vector
        @param delta tolerance of the L2 norm: [res_norm <= delta]
        @param max_restarts allowed restarts before failure (defaults to 0)
        @param s1 optional positive scale factors for {% $P_1 - b^{-1}$ %},
                  where {% $P_1$ %} is the left preconditioner.
        @param s2 optional positive scale factors for {% $P_2 x$ %},
                  where {% $P_2$ %} is the right preconditioner.
        @param psolve optionally solves the preconditioner system.

        @raise ConvFailure Failed to converge
        @raise QRfactFailure QRfact found a singular matrix
        @raise PSolveFailure The [psolve] function failed
        @raise ATimesFailure The [atimes] function failed
        @raise PSetFailure pset failed
        @raise GSFailure The Gram-Schmidt routine failed.
        @raise QRSolFailure QRsol found a singular [R]. *)
    val solve : ('d, 'k) t
                -> x:('d, 'k) Nvector.t
                -> b:('d, 'k) Nvector.t
                -> delta:float
                -> ?max_restarts:int
                -> ?s1:(('d, 'k) Nvector.t)
                -> ?s2:(('d, 'k) Nvector.t)
                -> ?psolve:('d psolve)
                -> 'd atimes
                -> preconditioning_type
                -> gramschmidt_type 
                -> bool * float * int * int
  end

(** The Scaled Preconditioned Flexible Generalized Minimum Residual (GMRES)
    method. *)
module SPFGMR :
  sig
    
    (** An instance of the SPFGMR solver.

        @nocvode <node9#ss:spfgmr> The SPFGMR Module
        @since 2.6.0 *)
    type ('d, 'k) t

    (** [make lmax temp] returns a solver session. [lmax] is the maximum
        Krylov subspace dimension to use, and [temp] sets the problem size.

        @nocvode <node9#ss:spfgmr> SpfgmrMalloc *)
    val make  : int -> ('d, 'k) Nvector.t -> ('d, 'k) t

    (** Solves the linear system [Ax = b] using the SPFGMR iterative method.
        The [atimes] function computes the matrix vector product [Ax], the
        other arguments are described below. The function returns a tuple
        [(solved, res_norm, nli, nps)] where [solved] indicates whether the
        system converged, [res_norm] is the L2 norm of the scaled
        preconditioned residual
        {% $\lVert s_1 P_1^{-1} (b - Ax) \rVert_{L2}$ %}, and,
        [nli] and [nps] count, respectively, linear iterations performed and
        calls to [psolve]. Repeated calls can be made to [solve] with
        varying input arguments, but a new session must be created if either
        the problem size or the maximum Krylov dimension change.

        @nocvode <node9#ss:spgmr> SpfgmrSolve
        @param x initial guess on entry; result on return
        @param b right-hand side vector
        @param delta tolerance of the L2 norm: [res_norm <= delta]
        @param max_restarts allowed restarts before failure (defaults to 0)
        @param max_iters maximum number of iterations (defaults to [lmax]).
        @param s1 optional positive scale factors for {% $P_1 - b^{-1}$ %},
                  where {% $P_1$ %} is the left preconditioner.
        @param s2 optional positive scale factors for {% $P_2 x$ %},
                  where {% $P_2$ %} is the right preconditioner.
        @param psolve optionally solves the preconditioner system.

        @raise ConvFailure Failed to converge
        @raise QRfactFailure QRfact found a singular matrix
        @raise PSolveFailure The [psolve] function failed
        @raise ATimesFailure The [atimes] function failed
        @raise PSetFailure pset failed
        @raise GSFailure The Gram-Schmidt routine failed.
        @raise QRSolFailure QRsol found a singular [R]. *)
    val solve : ('d, 'k) t
                -> x:('d, 'k) Nvector.t
                -> b:('d, 'k) Nvector.t
                -> delta:float
                -> ?max_restarts:int
                -> ?max_iters:int
                -> ?s1:(('d, 'k) Nvector.t)
                -> ?s2:(('d, 'k) Nvector.t)
                -> ?psolve:('d psolve)
                -> 'd atimes
                -> preconditioning_type
                -> gramschmidt_type 
                -> bool * float * int * int
  end

(** The Scaled Preconditioned Biconjugate Gradient Stabilized (Bi-CGStab)
    method. *)
module SPBCG :
  sig
    
    (** An instance of the SPBCG solver.

        @cvode <node9#ss:spgmr> The SPBCG Module *)
    type ('d, 'k) t

    (** [make lmax temp] returns a solver session. [lmax] is the maximum
        Krylov subspace dimension to use, and [temp] sets the problem size.

        @cvode <node9#ss:spbcg> SpbcgMalloc *)
    val make  : int -> ('d, 'k) Nvector.t -> ('d, 'k) t

    (** Solves the linear system [Ax = b] using the SPBCG iterative method.
        The [atimes] function computes the matrix vector product [Ax], the
        other arguments are described below. The function returns a tuple
        [(solved, res_norm, nli, nps)] where [solved] indicates whether the
        system converged, [res_norm] is the L2 norm of the scaled
        preconditioned residual
        {% $\lVert s_b P_1^{-1} (b - Ax) \rVert_{L2}$ %}, and,
        [nli] and [nps] count, respectively, linear iterations performed and
        calls to [psolve]. Repeated calls can be made to [solve] with
        varying input arguments, but a new session must be created if either
        the problem size or the maximum Krylov dimension change.

        @cvode <node9#ss:spgmr> SpgmrSolve
        @param x initial guess on entry; result on return
        @param b right-hand side vector
        @param delta tolerance of the L2 norm: [res_norm <= delta]
        @param sx optional positive scale factors for [x],
        @param sb optional positive scale factors for [b],
        @param psolve optionally solves the preconditioner system.

        @raise ConvFailure Failed to converge
        @raise PSolveFailure The [psolve] function failed
        @raise ATimesFailure The [atimes] function failed
        @raise PSetFailure pset failed *)
    val solve : ('d, 'k) t
                -> x:('d, 'k) Nvector.t
                -> b:('d, 'k) Nvector.t
                -> delta:float
                -> ?sx:(('d, 'k) Nvector.t)
                -> ?sb:(('d, 'k) Nvector.t)
                -> ?psolve:('d psolve)
                -> 'd atimes
                -> preconditioning_type
                -> bool * float * int * int

 end


(** The Scaled Preconditioned Transpose-Free Quasi-Minimal Residual
    (SPTFQMR) method *)
module SPTFQMR :
  sig
    
    (** An instance of the SPTFQMR solver.

        @cvode <node9#ss:sptfqmr> The SPTFQMR Module *)
    type ('d, 'k) t

    (** [make lmax temp] returns a solver session. [lmax] is the maximum
        Krylov subspace dimension to use, and [temp] sets the problem size.

        @cvode <node9#ss:sptfqmr> SptfqmrMalloc *)
    val make  : int -> ('d, 'k) Nvector.t -> ('d, 'k) t

    (** Solves the linear system [Ax = b] using the SPTFQMR iterative method.
        The [atimes] function computes the matrix vector product [Ax], the
        other arguments are described below. The function returns a tuple
        [(solved, res_norm, nli, nps)] where [solved] indicates whether the
        system converged, [res_norm] is the L2 norm of the scaled
        preconditioned residual
        {% $\lVert s_b P_1^{-1} (b - Ax) \rVert_{L2}$ %}, and,
        [nli] and [nps] count, respectively, linear iterations performed and
        calls to [psolve]. Repeated calls can be made to [solve] with
        varying input arguments, but a new session must be created if either
        the problem size or the maximum Krylov dimension change.

        @cvode <node9#ss:spgmr> SpgmrSolve
        @param x initial guess on entry; result on return
        @param b right-hand side vector
        @param delta tolerance of the L2 norm: [res_norm <= delta]
        @param sx optional positive scale factors for [x],
        @param sb optional positive scale factors for [b],
        @param psolve optionally solves the preconditioner system.

        @raise ConvFailure Failed to converge
        @raise PSolveFailure The [psolve] function failed
        @raise ATimesFailure The [atimes] function failed
        @raise PSetFailure pset failed *)
    val solve : ('d, 'k) t
                -> x:('d, 'k) Nvector.t
                -> b:('d, 'k) Nvector.t
                -> delta:float
                -> ?sx:(('d, 'k) Nvector.t)
                -> ?sb:(('d, 'k) Nvector.t)
                -> ?psolve:('d psolve)
                -> 'd atimes
                -> preconditioning_type
                -> bool * float * int * int

 end

(** The Preconditioned Conjugate-Gradient (PCG) method. *)
module PCG :
  sig
    
    (** An instance of the PCG solver.

        @noarkode <node9#ss:pcg> The PCG Module
        @since 2.6.0 *)
    type ('d, 'k) t

    (** [make lmax temp] returns a solver session. [lmax] is the maximum
        Krylov subspace dimension to use, and [temp] sets the problem size.

        @noarkode <node9#ss:pcg> PcgMalloc *)
    val make  : int -> ('d, 'k) Nvector.t -> ('d, 'k) t

    (** Solves the linear system [Ax = b] using the PCG iterative method.
        The [atimes] function computes the matrix vector product [Ax], the
        other arguments are described below. The function returns a tuple
        [(solved, res_norm, nli, nps)] where [solved] indicates whether the
        system converged, [res_norm] is the L2 norm of the scaled
        preconditioned residual
        {% $\lVert b - Ax \rVert_{L2}$ %}, and,
        [nli] and [nps] count, respectively, linear iterations performed and
        calls to [psolve]. Repeated calls can be made to [solve] with
        varying input arguments, but a new session must be created if either
        the problem size or the maximum Krylov dimension change.

        @arkode <node9#ss:pcg> PcgSolve
        @param x initial guess on entry; result on return
        @param b right-hand side vector
        @param delta tolerance of the L2 norm: [res_norm <= delta]
        @param w used in computing the residual norm for stopping solver
        @param psolve optionally solves the preconditioner system.

        @raise ConvFailure Failed to converge
        @raise PSolveFailure The [psolve] function failed
        @raise ATimesFailure The [atimes] function failed
        @raise PSetFailure pset failed *)
    val solve : ('d, 'k) t
                -> x:('d, 'k) Nvector.t
                -> b:('d, 'k) Nvector.t
                -> delta:float
                -> w:('d, 'k) Nvector.t
                -> ?psolve:('d psolve)
                -> 'd atimes
                -> preconditioning_type
                -> bool * float * int * int
  end

