(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Generic iterative linear Solvers.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNLinearSolver module
    @since 3.0.0 *)

(** An iterative linear solver.
    The type variables specify the {!Nvector.nvector} data (['data]) and
    kind (['kind]), and the iterative method (['iter]).

    @nocvode <node> Description of the SUNLinearSolver module
    @nocvode <node> SUNLinearSolver *)
type ('data, 'kind, 'iter) t
  = ('data, 'kind, 'iter) Lsolver_impl.Iterative.t

(** The type of Gram-Schmidt orthogonalization in iterative linear solvers.

    @nocvode <node> ModifiedGS/ClassicalGS *)
type gramschmidt_type = Lsolver_impl.Iterative.gramschmidt_type =
  | ModifiedGS   (** Modified Gram-Schmidt orthogonalization
                     {cconst MODIFIED_GS} *)
  | ClassicalGS  (** Classical Gram Schmidt orthogonalization
                     {cconst CLASSICAL_GS} *)

(** {3:solvers Solvers} *)

(** Krylov iterative solver using the scaled preconditioned biconjugate
    stabilized (Bi-CGStab) method. The [maxl] arguments gives the maximum
    dimension of the Krylov subspace (defaults to 5). The nvector argument
    is used as a template.

    @nocvode <node> SUNSPBCGS *)
val spbcgs : ?maxl:int -> ('d, 'k) Nvector.t -> ('d, 'k, [`Spbcgs]) t

(** Krylov iterative solver using the scaled preconditioned flexible
    generalized minimum residual (GMRES) method. The [maxl] arguments gives
    the maximum dimension of the Krylov subspace (defaults to 5). The
    nvector argument is used as a template.

    NB: [max_restarts] is ignored in Sundials < 3.0.0.

    @nocvode <node> SUNSPFGMR *)
val spfgmr : ?maxl:int -> ?max_restarts:int
             -> ?gs_type:gramschmidt_type
             -> ('d, 'k) Nvector.t -> ('d, 'k, [`Spfgmr]) t

(** Krylov iterative solver using the scaled preconditioned generalized
    minimum residual (GMRES) method. The [maxl] arguments gives the maximum
    dimension of the Krylov subspace (defaults to 5). The nvector argument
    is used as a template.

    NB: [max_restarts] is ignored in Sundials < 3.0.0.

    @nocvode <node> SUNSPGMR *)
val spgmr : ?maxl:int -> ?max_restarts:int
            -> ?gs_type:gramschmidt_type
            -> ('d, 'k) Nvector.t -> ('d, 'k, [`Spgmr]) t

(** Krylov iterative with the scaled preconditioned transpose-free
    quasi-minimal residual (SPTFQMR) method. The [maxl] arguments gives the
    maximum dimension of the Krylov subspace (defaults to 5). The nvector
    argument is used as a template.

    @nocvode <node> SUNSPTFQMR *)
val sptfqmr : ?maxl:int -> ('d, 'k) Nvector.t -> ('d, 'k, [`Sptfqmr]) t

(** Krylov iterative solver using the preconditioned conjugate gradient
    (PCG) method. The [maxl] arguments gives the maximum dimension of the
    Krylov subspace (defaults to 5). The nvector argument is used as a
    template.

    @nocvode <node> SUNPCG *)
val pcg : ?maxl:int -> ('d, 'k) Nvector.t -> ('d, 'k, [`Pcg]) t

(** Custom iterative linear solvers. *)
module Custom : sig (* {{{ *)

  (** Used to distinguish custom linear solvers *)
  type 'lsolver tag = [`Custom of 'lsolver]

  (** A function [atimesfn v z] computes the action of a matrix on the
      vector [v], storing the result in [z]. *)
  type ('data, 'kind) atimesfn =
       ('data, 'kind) Nvector.t
    -> ('data, 'kind) Nvector.t
    -> unit

  (** Functions that set up any problem data in preparation for calls to
      [psolvefn]. *)
  type psetupfn = unit -> unit

  (** A function [psolvefn r z tol lr] that solves the preconditioner
      equation {% $Pz = r$ %} for the vector [z] such that
      {% $\left\lVert Pz - r \right\rVert_\mathrm{wrms} < \mathit{tol}$ %}.
      If [lr] is [true] then {% $P$ %} should be treated as a left
      preconditioner and otherwise as a right preconditioner. *)
  type ('data, 'kind) psolvefn =
       ('data, 'kind) Nvector.t
    -> ('data, 'kind) Nvector.t
    -> float
    -> bool
    -> unit

  (** The operations required to implement an iterative linear solver.
      Failure should be indicated by raising an exception (preferably
      one of the exceptions in the {!module:Lsolver} package). Raising
      {!exception:Sundials.RecoverableFailure} indicates a generic
      recoverable failure. *)
  type ('data, 'kind, 'lsolver) ops = {

    init : 'lsolver -> unit;
    (** Performs linear solver initalization. *)

    setup : 'lsolver -> unit;
    (** Performs linear solver setup. *)

    solve : 'lsolver
            -> ('data, 'kind) Nvector.t
            -> ('data, 'kind) Nvector.t
            -> float
            -> unit;
    (** The call [solve ls x b tol] should solve the linear system
        {% $Ax = b$ %} to within the weight 2-norm tolerance [tol].
        {% $A$ %} is only available indirectly via the [atimes] function. *)

    set_atimes
      : ('lsolver -> ('data, 'kind) atimesfn -> unit) option;
      (** Provides the linear solver with a problem-specific {!atimesfn}.
          The given function may only be used within [init], [setup],
          and [solve]. *)

    set_preconditioner
      : ('lsolver
         -> psetupfn option
         -> ('data, 'kind) psolvefn option
         -> unit) option;
      (** Provides the linear solver with preconditioner routines.
          The given functions may only be used within [init], [setup],
          and [solve]. *)

    set_scaling_vectors
      : ('lsolver
         -> ('data, 'kind) Nvector.t option
         -> ('data, 'kind) Nvector.t option
         -> unit) option;
    (** Passes the left/right scaling vectors for use in [solve].
        The call [set_scaling_vectors ls s1 s2] provides diagonal matrices
        of scale factors for solving the system
        {% $\tilde{A}\tilde{x} = \tilde{b}$ %}
        where
        {% $\tilde{A} = S_1 P_1^{-1} A P_2^{-1} S_2^{-1}$ %},
        {% $\tilde{b} = S_1 P_1^{-1} b$ %}, and
        {% $\tilde{x} = S_2 P_2 x$ %}.
        A [None] argument indicates an identity scaling matrix. *)

    get_num_iters : ('lsolver -> int) option;
    (** The number of linear iterations performed in the last
       [solve] call. *)

    get_res_norm : ('lsolver -> float) option;
    (** The final residual norm from the last [solve] call. *)

    get_res_id : ('lsolver -> ('data, 'kind) Nvector.t) option;
    (** The preconditioned initial residual vector. This vector may be
        requested if the iterative method computes the preconditioned
        initial residual and returns from [solve] successfully without
        performing any iterations (i.e., either the initial guess or the
        preconditioner is sufficiently accurate). *)

    get_work_space : ('lsolver -> int * int) option;
    (** Return the storage requirements for the linear solver.
        The result [(lrw, liw)] gives the number of words used for
        storing real values ([lrw]) and the number of words used
        for storing integer values ([liw]). * *)
  }

  (** Create an iterative linear solver given a set of operations and an
      internal state.

      NB: This feature is only available for Sundials >= 3.0.0. *)
  val make : ('data, 'kind, 'lsolver) ops
             -> 'lsolver
             -> ('data, 'kind, 'lsolver tag) t

  (** Return the internal state from an custom iterative linear solver. *)
  val unwrap : ('data, 'kind, 'lsolver tag) t -> 'lsolver

end (* }}} *)

(** Low-level routines on arrays. *)
module Algorithms : sig (* {{{ *)

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
      [h] have already been factored. The computed Givens rotation has the
      form {% $\begin{bmatrix} c & -s \\ s & c \end{bmatrix}$ %}. It is
      stored in the [2n] elements of [q] as [[|c; s; c; s; ...; c; s|]].

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
    - [k] specifies the vector in [v] to be orthogonalized against previous
          ones, and,
    - [p] is the number of previous vectors in [v] to orthogonalize against.

    The vector [v[k]] is orthogonalized against the [p] unit vectors at
    [v.{k-1}], [v.{k-2}], ..., [v.{k-p}].
    The matrix [h] must be allocated row-wise so that the [(i,j)]th entry is
    [h.{i}.{j}].
    The inner products are computed, {% $\mathtt{h.\\{}i\mathtt{, k-1\\}} =
      \mathtt{v.\\{}i\mathtt{\\}} \cdot \mathtt{v.\\{k\\}}$ %}, for
    {% $i=\max(0, \mathtt{k}-\mathtt{p})\ldots \mathtt{k}-1$ %}.
    The orthogonalized [v.{k}] is {b not} normalized and is stored over the
    old [v.{k}]. The function returns the Euclidean norm of the orthogonalized
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
    - [k] specifies the vector in [v] to be orthogonalized against previous
          ones, and,
    - [p] is the number of previous vectors in [v] to orthogonalize against.
    - [temp] and [s] are used as workspaces.

    The vector [v[k]] is orthogonalized against the [p] unit vectors at
    [v.{k-1}], [v.{k-2}], ..., [v.{k-p}].
    The matrix [h] must be allocated row-wise so that the [(i,j)]th entry is
    [h.{i}.{j}].
    The inner products are computed, {% $\mathtt{h.\\{}i\mathtt{, k-1\\}} =
      \mathtt{v.\\{}i\mathtt{\\}} \cdot \mathtt{v.\\{k\\}}$ %}, for
    {% $i=\max(0, \mathtt{k}-\mathtt{p})\ldots \mathtt{k}-1$ %}.
    The orthogonalized [v.{k}] is {b not} normalized and is stored over the
    old [v.{k}]. The function returns the Euclidean norm of the orthogonalized
    vector. *)
  val classical_gs : (('d, 'k) Nvector.t) array
                    -> Sundials.RealArray2.t
                    -> int
                    -> int
                    -> ('d, 'k) Nvector.t
                    -> Sundials.RealArray.t
                    -> float

end (* }}} *)

(** {3:parameters Solver parameters} *)

(** Updates the number of linear solver iterations to allow.

    @nocvode <node> SUNSPBCGSSetMaxl
    @nocvode <node> SUNSPTFQMRSetMaxl
    @nocvode <node> SUNPCGSetMaxl *)
val set_maxl : ('d, 'k, [< `Spbcgs|`Sptfqmr|`Pcg]) t -> int -> unit

(** Sets the Gram-Schmidt orthogonalization to use.

    @nocvode <node> SUNSPGMRSetGSType
    @nocvode <node> SUNSPFGMRSetGSType *)
val set_gs_type :
  ('d, 'k, [< `Spfgmr|`Spgmr]) t -> gramschmidt_type -> unit

(** Sets the number of GMRES restarts to allow.

    NB: This feature is not supported by Sundials < 3.0.0.

    @nocvode <node> SUNSPGMRSetMaxRestarts
    @nocvode <node> SUNSPFGMRSetMaxRestarts *)
val set_max_restarts : ('d, 'k, [< `Spfgmr|`Spgmr]) t -> int -> unit

(** The type of preconditioning in Krylov solvers.

    @nocvode <node> Preconditioning *)
type preconditioning_type = Lsolver_impl.Iterative.preconditioning_type =
  | PrecNone    (** No preconditioning *)
  | PrecLeft    (** {% $(P^{-1}A)x = P^{-1}b$ %} *)
  | PrecRight   (** {% $(AP^{-1})Px = b$ %} *)
  | PrecBoth    (** {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)

(** Change the preconditioning direction without modifying callback
    functions.

    Raises {!Lsolver.IllegalPrecType} if the current preconditioner
    is {{!preconditioning_type}PrecNone} and the given argument is not
    (since no callback functions are specified in this case. May raise
    {!Lsolver.IllegalPrecType} if the given type is not allowed by the
    underlying solver.

    @nocvode <node> SUNPCGSetPrecType
    @nocvode <node> SUNSPBCGSSetPrecType
    @nocvode <node> SUNSPFGMRSetPrecType
    @nocvode <node> SUNSPGMRSetPrecType
    @nocvode <node> SUNSPTFQMRSetPrecType *)
val set_prec_type : ('d, 'k, 'f) t -> preconditioning_type -> unit

