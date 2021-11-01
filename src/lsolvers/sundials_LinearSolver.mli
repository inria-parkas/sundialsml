(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2020 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Generic linear solvers.

    Sundials provides a set of functions for instantiating linear solvers from
    two families: {!module:Direct} and {!module:Iterative}. Any instance may
    be associated with at most one solver session.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNLinearSolver module
    @since 3.0.0 *)

open Sundials

(** A linear solver.
    The type variables specify the Jacobian matrix (['matrix]), the
    {!Nvector.nvector} data (['data]) and kind (['kind]), and a
    ['tag] used to identify specific solver features (like the iterative
    method).

    A linear solver of this type must be converted to session-specific
    form by {!Cvode.Dls.solver}, {!Ida.Dls.solver}, {!Cvode.Spils.solver},
    {!Ida.Spils.solver}, etc., before being
    attached to a session via [init] or [reinit].

  @nocvode <node> Description of the SUNLinearSolver module
  @nocvode <node> SUNLinearSolver *)
type ('matrix, 'data, 'kind, 'tag) t
  = ('matrix, 'data, 'kind, 'tag) Sundials_LinearSolver_impl.linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type ('mat, 'kind, 'tag) serial_t
  = ('mat, Nvector_serial.data, [>Nvector_serial.kind] as 'kind, 'tag) t

(** Broadly classifies the operations provided by a linear solver and its
    operating principle. *)
type linear_solver_type =
  | Direct          (** Performs an exact computation based on a matrix. *)
  | Iterative       (** Computes an inexact approximation without a matrix. *)
  | MatrixIterative (** Computes an inexact approximation using a matrix. *)
  | MatrixEmbedded  (** Computes without an explicit input matrix *)

(** The identifier of a linear solver. *)
type linear_solver_id =
  | Band               (* {cconst SUNLINEARSOLVER_BAND} *)
  | Dense              (* {cconst SUNLINEARSOLVER_DENSE} *)
  | Klu                (* {cconst SUNLINEARSOLVER_KLU} *)
  | LapackBand         (* {cconst SUNLINEARSOLVER_LAPACKBAND} *)
  | LapackDense        (* {cconst SUNLINEARSOLVER_LAPACKDENSE} *)
  | Pcg                (* {cconst SUNLINEARSOLVER_PCG} *)
  | Spbcgs             (* {cconst SUNLINEARSOLVER_SPBCGS} *)
  | Spfgmr             (* {cconst SUNLINEARSOLVER_SPFGMR} *)
  | Spgmr              (* {cconst SUNLINEARSOLVER_SPGMR} *)
  | Sptfqmr            (* {cconst SUNLINEARSOLVER_SPTFQMR} *)
  | Superludist        (* {cconst SUNLINEARSOLVER_SUPERLUDIST} *)
  | Superlumt          (* {cconst SUNLINEARSOLVER_SUPERLUMT} *)
  | Custom             (* {cconst SUNLINEARSOLVER_CUSTOM} *)

(** {3:callbacks Callback Routines} *)

(** A function [atimesfn v z] computes the action of the system
    matrix on the vector [v], storing the result in [z]. The matrix is
    represented implicitly by the effect of the function.

    If a problem occurs the function raises {!ATimesFailure}. *)
type 'data atimesfn = 'data -> 'data -> unit

(** Functions that set up any problem data in preparation for calls to
    [psolvefn].

    If a problem occurs the function raises {!PSetFailure}. *)
type psetupfn = unit -> unit

(** A function [psolvefn r z tol lr] that solves the preconditioner
    equation {% $Pz = r$ %} for the vector [z] such that
    {% $\left\lVert Pz - r \right\rVert_\mathrm{wrms} < \mathit{tol}$ %}.
    If [lr] is [true] then {% $P$ %} should be treated as a left
    preconditioner and otherwise as a right preconditioner.

    If a problem occurs the function raises {!PSolveFailure}. *)
type 'data psolvefn = 'data -> 'data -> float -> bool -> unit

(** {2:lsolvers Linear Solver Families} *)

(** Direct Linear Solvers *)
module Direct : sig (* {{{ *)
  (** Definitions in this module are more conveniently accessed
      through session-specific direct linear solver modules like
      {!Cvode.Dls} and {!Ida.Dls}.  For example, [Cvode.Dls.dense] is
      an alias for {!dense}.  *)

  (** {3:solvers Types} *)

  (** Creates a direct linear solver on dense matrices. The nvector and matrix
    argument are used to determine the linear system size and to assess
    compatibility with the linear solver implementation.
    The matrix is used internally after the linear solver is attached to a
    session.

  @nocvode <node> SUNLinSol_Dense *)
  val dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k, [`Dls]) serial_t

  (** Creates a direct linear solver on dense matrices using LAPACK.
      See {!dense}. Only available if
      {{!Sundials_Config.lapack_enabled}Config.lapack_enabled}.

  @nocvode <node> SUNLinSol_LapackDense *)
  val lapack_dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k, [`Dls]) serial_t

  (** Creates a direct linear solver on banded matrices. The nvector and matrix
    argument are used to determine the linear system size and to assess
    compatibility with the linear solver implementation.
    The matrix is used internally after the linear solver is attached to a
    session.

  @nocvode <node> SUNLinSol_Band *)
  val band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k, [`Dls]) serial_t

  (** Creates a direct linear solver on banded matrices using LAPACK.
      See {!band}.
      Only available if {{!Sundials_Config.lapack_enabled}Config.lapack_enabled}.

  @nocvode <node> SUNLinSol_LapackBand *)
  val lapack_band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k, [`Dls]) serial_t

  (** KLU direct linear solver operating on sparse matrices (requires KLU). *)
  module Klu : sig (* {{{ *)

    (** The ordering algorithm used for reducing fill. *)
    type ordering =
      Sundials_LinearSolver_impl.Klu.ordering =
      | Amd      (** Approximate minimum degree permutation. *)
      | ColAmd   (** Column approximate minimum degree permutation. *)
      | Natural  (** Natural ordering. *)

    (** Creates a direct linear solver on sparse matrices using KLU. The
      nvector and matrix argument are used to determine the linear system
      size and to assess compatibility with the linear solver implementation.
      The matrix is used internally after the linear solver is attached to a
      session.

      @raise Config.NotImplementedBySundialsVersion Solver not available.
      @nocvode <node> SUNLinSol_KLU *)
    val make :
      ?ordering:ordering
      -> 'k Nvector_serial.any
      -> ('s, 'k) Matrix.sparse
      -> ('s Matrix.Sparse.t, 'k, [`Dls|`Klu]) serial_t

    (** Reinitializes memory and flags for a new factorization (symbolic and
      numeric) at the next solver setup call. In the call [reinit ls a nnz],
      [a] is the Jacobian matrix, which is reinitialized with the given
      number of non-zeros if [nnz] if given. New symbolic and numeric
      factorizations will be completed at the next solver step.

      @nocvode <node> SUNLinSol_KLUReInit *)
    val reinit : ('s Matrix.Sparse.t, 'k, [>`Klu]) serial_t
                 -> ('s, 'k) Matrix.sparse -> ?nnz:int -> unit -> unit

    (** Sets the ordering algorithm used to minimize fill-in.

      @nocvode <node> SUNLinSol_KLUSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k, [>`Klu]) serial_t
                       -> ordering -> unit

  end (* }}} *)

  (** Creates a direct linear solver on sparse matrices using KLU.
    See {!Klu.make}.

    @raise Config.NotImplementedBySundialsVersion Solver not available.
    @nocvode <node> SUNLinSol_KLU *)
  val klu :
    ?ordering:Klu.ordering
    -> 'k Nvector_serial.any
    -> ('s, 'k) Matrix.sparse
    -> ('s Matrix.Sparse.t, 'k, [`Klu|`Dls]) serial_t

  (** SuperLUMT direct linear solver operating on sparse matrices (requires
    SuperLUMT). *)
  module Superlumt : sig (* {{{ *)

    (** The ordering algorithm used for reducing fill. *)
    type ordering =
      Sundials_LinearSolver_impl.Superlumt.ordering =
      | Natural       (** Natural ordering. *)
      | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
      | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
      | ColAmd        (** Column approximate minimum degree permutation. *)

    (** Creates a direct linear solver on sparse matrices using SuperLUMT. The
      nvector and matrix argument are used to determine the linear system
      size and to assess compatibility with the linear solver implementation.
      The matrix is used internally after the linear solver is attached to a
      session.

      NB: The {{!Sundials_Matrix.Sparse.csr}Matrix.Sparse.csr} format is only
          supported for
          {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0.

      @raise Config.NotImplementedBySundialsVersion Solver not available.
      @nocvode <node> SUNLinSol_SuperLUMT *)
    val make :
      ?ordering:ordering
      -> nthreads:int
      -> 'k Nvector_serial.any
      -> ('s, 'k) Matrix.sparse
      -> ('s Matrix.Sparse.t, 'k, [`Dls|`Slu]) serial_t

    (** Sets the ordering algorithm used to minimize fill-in.

      @nocvode <node> SUNLinSol_SuperLUMTSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k, [>`Slu]) serial_t
                       -> ordering -> unit

  end (* }}} *)

  (** Creates a direct linear solver on sparse matrices using SuperLUMT.
    See {!Superlumt.make}.

    NB: The {{!Sundials_Matrix.Sparse.csr}Matrix.Sparse.csr} format is only
        supported for
        {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0.

    @raise Config.NotImplementedBySundialsVersion Solver not available.
    @nocvode <node> SUNLinSol_SuperLUMT *)
  val superlumt :
    ?ordering:Superlumt.ordering
    -> nthreads:int
    -> 'k Nvector_serial.any
    -> ('s, 'k) Matrix.sparse
    -> ('s Matrix.Sparse.t, 'k, [>`Slu|`Dls]) serial_t

end (* }}} *)

(** Iterative Linear Solvers *)
module Iterative : sig (* {{{ *)

  (** Definitions in this module are more conveniently accessed
      through session-specific iterative linear solver modules like
      {!Cvode.Spils} and {!Ida.Spils}.  For example,
      [Cvode.Spils.spbcgs] is an alias for {!spbcgs}.  *)

  (** {3:solvers Types} *)

  (** The type of Gram-Schmidt orthogonalization in iterative linear solvers.

    @nocvode <node> ModifiedGS/ClassicalGS *)
  type gramschmidt_type =
    Sundials_LinearSolver_impl.Iterative.gramschmidt_type =
    | ModifiedGS   (** Modified Gram-Schmidt orthogonalization
                       {cconst MODIFIED_GS} *)
    | ClassicalGS  (** Classical Gram Schmidt orthogonalization
                       {cconst CLASSICAL_GS} *)

  (** {3:solvers Solvers} *)

  (** Krylov iterative solver using the scaled preconditioned biconjugate
    stabilized (Bi-CGStab) method. The [maxl] arguments gives the maximum
    dimension of the Krylov subspace (defaults to 5). The nvector argument
    is used as a template.

    @nocvode <node> SUNLinSol_SPBCGS *)
  val spbcgs : ?maxl:int -> ('d, 'k) Nvector.t
               -> ('m, 'd, 'k, [`Iter|`Spbcgs]) t

  (** Krylov iterative solver using the scaled preconditioned flexible
    generalized minimum residual (GMRES) method. The [maxl] arguments gives
    the maximum dimension of the Krylov subspace (defaults to 5). The
    nvector argument is used as a template.

    NB: [max_restarts] is ignored by CVODE, CVODES, and ARKODE
        for {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

    @nocvode <node> SUNLinSol_SPFGMR *)
  val spfgmr : ?maxl:int -> ?max_restarts:int
               -> ?gs_type:gramschmidt_type
               -> ('d, 'k) Nvector.t
               -> ('m, 'd, 'k, [`Iter|`Spfgmr]) t

  (** Krylov iterative solver using the scaled preconditioned generalized
    minimum residual (GMRES) method. The [maxl] arguments gives the maximum
    dimension of the Krylov subspace (defaults to 5). The nvector argument
    is used as a template.

    NB: [max_restarts] is ignored by CVODE, CVODES, and ARKODE
        for {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

    @nocvode <node> SUNLinSol_SPGMR *)
  val spgmr : ?maxl:int -> ?max_restarts:int
              -> ?gs_type:gramschmidt_type
              -> ('d, 'k) Nvector.t
              -> ('m, 'd, 'k, [`Iter|`Spgmr]) t

  (** Krylov iterative with the scaled preconditioned transpose-free
    quasi-minimal residual (SPTFQMR) method. The [maxl] arguments gives the
    maximum dimension of the Krylov subspace (defaults to 5). The nvector
    argument is used as a template.

    @nocvode <node> SUNLinSol_SPTFQMR *)
  val sptfqmr : ?maxl:int -> ('d, 'k) Nvector.t
                -> ('m, 'd, 'k, [`Iter|`Sptfqmr]) t

  (** Krylov iterative solver using the preconditioned conjugate gradient
    (PCG) method. The [maxl] arguments gives the maximum dimension of the
    Krylov subspace (defaults to 5). The nvector argument is used as a
    template.

    @nocvode <node> SUNLinSol_PCG *)
  val pcg : ?maxl:int -> ('d, 'k) Nvector.t
            -> ('m, 'd, 'k, [`Iter|`Pcg]) t

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
      The call [qr_fact n h q factored], where [h] is the [n+1] by [n]
      Hessenberg matrix (stored row-wise), [q] stores the computed Givens
      rotation, and [factored=true] indicates that the first [n-1] columns of
      [h] have already been factored. The computed Givens rotation has the
      form {% $\begin{bmatrix} c & -s \\ s & c \end{bmatrix}$ %}. It is
      stored in the [2n] elements of [q] as [[|c; s; c; s; ...; c; s|]].

      @raise Matrix.ZeroDiagonalElement Zero found in matrix diagonal *)
    val qr_fact : int
                  -> RealArray2.t
                  -> RealArray.t
                  -> bool
                  -> unit

    (** Solve the linear least squares problem. In
      [qr_sol n h q b], [h] and [q] are, respectively, the upper triangular
      factor $R$ of the original Hessenberg matrix and [Q] the Givens
      rotations used to factor it—both computed by {!qr_fact}. The function
      computes the [n+1] elements of [b] to solve $Rx = Qb$.

      @raise Matrix.ZeroDiagonalElement Zero found in matrix diagonal *)
    val qr_sol : int
                 -> RealArray2.t
                 -> RealArray.t
                 -> RealArray.t
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
                      -> RealArray2.t
                      -> int
                      -> int
                      -> float

    (** Performs a classical Gram-Schmidt orthogonalization. In
      [classical_gs v h k p s temp],
    - [v] is an array of at least [k + 1] vectors with an L2-norm of 1,
    - [h] is the output [k] by [k] Hessenberg matrix of inner products,
    - [k] specifies the vector in [v] to be orthogonalized against previous
          ones,
    - [p] is the number of previous vectors in [v] to orthogonalize against,
        and,
    - [s] and [temp] are arrays of at least [k + 1] elements used
      as workspaces.

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
                       -> RealArray2.t
                       -> int
                       -> int
                       -> RealArray.t
                       -> (('d, 'k) Nvector.t) array
                       -> float

  end (* }}} *)

  (** {3:parameters Solver parameters} *)

  (** Updates the number of linear solver iterations to allow.

    @nocvode <node> SUNLinSol_SPBCGSSetMaxl
    @nocvode <node> SUNLinSol_SPTFQMRSetMaxl
    @nocvode <node> SUNLinSol_PCGSetMaxl *)
  val set_maxl : ('m, 'd, 'k, [< `Iter|`Spbcgs|`Sptfqmr|`Pcg]) t
                 -> int -> unit

  (** Sets the Gram-Schmidt orthogonalization to use.

    @nocvode <node> SUNLinSol_SPGMRSetGSType
    @nocvode <node> SUNLinSol_SPFGMRSetGSType *)
  val set_gs_type : ('m, 'd, 'k, [< `Iter|`Spfgmr|`Spgmr]) t
                    -> gramschmidt_type -> unit

  (** Sets the number of GMRES restarts to allow.

    NB: This feature is not supported by CVODE, CVODES, and ARKODE
        for {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

    @nocvode <node> SUNLinSol_SPGMRSetMaxRestarts
    @nocvode <node> SUNLinSol_SPFGMRSetMaxRestarts *)
  val set_max_restarts : ('m, 'd, 'k, [< `Iter|`Spfgmr|`Spgmr]) t
                         -> int -> unit

  (** The type of preconditioning in Krylov solvers.

    @nocvode <node> Preconditioning *)
  type preconditioning_type =
    Sundials_LinearSolver_impl.Iterative.preconditioning_type =
    | PrecNone    (** No preconditioning *)
    | PrecLeft    (** {% $(P^{-1}A)x = P^{-1}b$ %} *)
    | PrecRight   (** {% $(AP^{-1})Px = b$ %} *)
    | PrecBoth    (** {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)

  (** Change the preconditioning direction without modifying callback
    functions.

    Raises {!IllegalPrecType} if the current preconditioner
    is {{!preconditioning_type}PrecNone} and the given argument is not
    (since no callback functions are specified in this case. May raise
    {!IllegalPrecType} if the given type is not allowed by the
    underlying solver.

    @nocvode <node> SUNLinSol_PCGSetPrecType
    @nocvode <node> SUNLinSol_SPBCGSSetPrecType
    @nocvode <node> SUNLinSol_SPFGMRSetPrecType
    @nocvode <node> SUNLinSol_SPGMRSetPrecType
    @nocvode <node> SUNLinSol_SPTFQMRSetPrecType *)
  val set_prec_type : ('m, 'd, 'k, [>`Iter]) t
                      -> preconditioning_type -> unit

  (** Sets the output file for informative (non-error) messages. The default
      is to send such messages to stdout.
      The optional argument is a convenience for invoking {!set_print_level}.

      Sundials must be built with {cconst SUNDIALS_BUILD_WITH_MONITORING} to
      use this function.

      @nocvode <node> SUNLinSolSetInfoFile_PCG
      @nocvode <node> SUNLinSolSetInfoFile_SPBCGS
      @nocvode <node> SUNLinSolSetInfoFile_SPFGMR
      @nocvode <node> SUNLinSolSetInfoFile_SPGMR
      @nocvode <node> SUNLinSolSetInfoFile_SPTFQMR
      @since 5.3.0 *)
  val set_info_file
    : ('m, 'd, 'k, [>`Iter]) t -> ?print_level:bool -> Sundials.Logfile.t -> unit

  (** Sets the level of output verbosity. When [false] (the default) no
      information is printed, when [true] the residual norm is printed for
      each linear iteration.

      Sundials must be built with {cconst SUNDIALS_BUILD_WITH_MONITORING} to
      use this function.

      @nocvode <node> SUNLinSolSetPrintLevel_PCG
      @nocvode <node> SUNLinSolSetPrintLevel_SPBCGS
      @nocvode <node> SUNLinSolSetPrintLevel_SPFGMR
      @nocvode <node> SUNLinSolSetPrintLevel_SPGMR
      @nocvode <node> SUNLinSolSetPrintLevel_SPTFQMR
      @since 5.3.0 *)
  val set_print_level : ('m, 'd, 'k, [>`Iter]) t -> bool -> unit

end (* }}} *)

(** Custom linear solvers. *)
module Custom : sig (* {{{ *)

  (** A function [atimesfn v z] computes the action of the system
      matrix on the vector [v], storing the result in [z]. The matrix is
      represented implicitly by the effect of the function.

      If a problem occurs the function raises {!ATimesFailure}. *)
  type ('data, 'kind) atimesfn =
    ('data, 'kind) Nvector.t
    -> ('data, 'kind) Nvector.t
    -> unit

  (** Functions that set up any problem data in preparation for calls to
      [psolvefn].

      If a problem occurs the function raises {!PSetFailure}. *)
  type psetupfn = unit -> unit

  (** A function [psolvefn r z tol lr] that solves the preconditioner
      equation {% $Pz = r$ %} for the vector [z] such that
      {% $\left\lVert Pz - r \right\rVert_\mathrm{wrms} < \mathit{tol}$ %}.
      If [lr] is [true] then {% $P$ %} should be treated as a left
      preconditioner and otherwise as a right preconditioner.

      If a problem occurs the function raises {!PSolveFailure}. *)
  type ('data, 'kind) psolvefn =
    ('data, 'kind) Nvector.t
    -> ('data, 'kind) Nvector.t
    -> float
    -> bool
    -> unit

  (** The operations required to implement an iterative linear solver.
      Failure should be indicated by raising an exception (preferably
      one of the exceptions in this package). Raising
      {!exception:Sundials.RecoverableFailure} indicates a generic
      recoverable failure. *)
  type ('matrix, 'data, 'kind, 'lsolver) ops = {
      solver_type : linear_solver_type;
      (** Broadly classifies the operations provided by a linear solver and
          its operating principle. *)

      solver_id : linear_solver_id;
      (** Identifies the linear solver. This value should normally be set
          to {{!linear_solver_id}Custom}. *)

      init : ('lsolver -> unit) option;
      (** Performs linear solver initalization. *)

      setup : ('lsolver -> 'matrix -> unit) option;
      (** Performs linear solver setup. *)

      solve : 'lsolver -> 'matrix -> 'data -> 'data -> float -> unit;
      (** The call [solve ls a x b tol] should solve the linear system
          {% $Ax = b$ %}.

          Direct solvers can ignore [tol].
          Matrix-free solvers can ignore [a] (relying instead on the [atimes]
          function).
          Iterative solvesr should attempt to solve to within the weighted
          2-norm tolerance, [tol]. *)

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
        : ('lsolver -> 'data option -> 'data option -> unit) option;
      (** Passes the left/right scaling vectors for use in [solve].
          The call [set_scaling_vectors ls s1 s2] provides diagonal matrices
          of scale factors for solving the system
          {% $\tilde{A}\tilde{x} = \tilde{b}$ %}
          where
          {% $\tilde{A} = S_1 P_1^{-1} A P_2^{-1} S_2^{-1}$ %},
          {% $\tilde{b} = S_1 P_1^{-1} b$ %}, and
          {% $\tilde{x} = S_2 P_2 x$ %}.
          A [None] argument indicates an identity scaling matrix. *)

      set_zero_guess : ('lsolver -> bool -> unit) option;
      (** Used to indicate whether the next call to the solve function will
          be maed with a zero initial guess. *)

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

      get_last_flag : ('lsolver -> int) option;
      (** Use to indicate the last error encountered by the linear solver. *)

      get_work_space : ('lsolver -> int * int) option;
      (** Return the storage requirements for the linear solver.
          The result [(lrw, liw)] gives the number of words used for
          storing real values ([lrw]) and the number of words used
          for storing integer values ([liw]). * *)

      (** Called by Iterative.set_prec_type and when a linear solver is
          associated with an integrator. *)
      set_prec_type :
        ('lsolver -> Iterative.preconditioning_type -> unit) option;
    }

  (** Create a linear solver given a set of operations and an internal state.
      The resulting solver is tagged with both [`Dls] and [`Iter], although
      its suitability depends on which operations are defined.
      Use a type cast to filter inappropriate tags.

      The [solver_type] may not be {{!type:linear_solver_type}MatrixEmbedded}.

    NB: This feature is only available for
        {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0. *)
  val make_with_matrix : ('matrix, 'data, 'kind, 'lsolver) ops
     -> 'lsolver
     -> ('matrixkind, 'matrix, 'data, 'kind) Matrix.t
     -> ('matrix, 'data, 'kind, [`Dls|`Iter|`Custom of 'lsolver]) t

  (** Create a linear solver given a set of operations and an internal state.
      The resulting solver is tagged with both [`Dls] and [`Iter], although
      its suitability depends on which operations are defined.
      Use a type cast to filter inappropriate tags.

      The [solver_type] may not be {{!type:linear_solver_type}Direct}.

      This feature is only available for
      {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0.

      The {{!type:linear_solver_type}MatrixEmbedded} solver_type is only
      available for
      {{!Sundials_Config.sundials_version}Config.sundials_version} >= 5.8.0. *)
  val make_without_matrix : (unit, 'data, 'kind, 'lsolver) ops
     -> 'lsolver
     -> (unit, 'data, 'kind, [`Iter|`MatE|`Custom of 'lsolver]) t

  (** Return the internal state from an custom iterative linear solver. *)
  val unwrap : ('m, 'data, 'kind, [>`Custom of 'lsolver]) t -> 'lsolver

  (** {3:dlscustom Wrapper for making custom Direct linear solvers.} *)

  (** The operations required to implement a direct linear solver.
      Failure should be indicated by raising an exception (preferably
      one of the exceptions in this package). Raising
      {!exception:Sundials.RecoverableFailure} indicates a generic
      recoverable failure. *)
  type ('matrix, 'data, 'kind, 'lsolver) dls_ops = {
      init : ('lsolver -> unit) option;
      (** Performs linear solver initalization. *)

      setup : ('lsolver -> 'matrix -> unit) option;
      (** Performs linear solver setup based on an updated matrix. *)

      solve : 'lsolver -> 'matrix -> 'data -> 'data -> float -> unit;
      (** The call [solve ls A x b tols] should solve the linear system
      {% $Ax = b$ %} to the requested tolerance. *)

      space : ('lsolver -> int * int) option;
      (** Return the storage requirements for the linear solver.
      The result [(lrw, liw)] gives the number of words used for
      storing real values ([lrw]) and the number of words used
      for storing integer values ([liw]). *)
    }

  (** Create a direct linear solver given a set of operations and an
      internal state.

      NB: This feature is only available for
          {{!Sundials_Config.sundials_version}Config.sundials_version} >= 3.0.0. *)
  val make_dls : ('matrix, 'data, 'kind, 'lsolver) dls_ops
           -> 'lsolver
           -> ('matrixkind, 'matrix, 'data, 'kind) Matrix.t
           -> ('matrix, 'data, 'kind, [`Dls|`Custom of 'lsolver]) t

end (* }}} *)

(** {2:ops Operations} *)

(** Set the linear solver's problem-specific {!atimesfn}.

    @nocvode <node> SUNLinSolSetATimes *)
val set_atimes : ('m, 'd, 'k, 't) t -> 'd atimesfn -> unit

(** Set the linear solver's preconditioner routines.

    @nocvode <node> SUNLinSolSetPreconditioner *)
val set_preconditioner :
    ('m, 'd, 'k, 't) t
  -> psetupfn
  -> 'd psolvefn
  -> unit

(** Sets the linear solver's left/right scaling vectors for use in {!solve}.
    The call [set_scaling_vectors ls s1 s2] provides diagonal matrices
    of scale factors for solving the system
    {% $\tilde{A}\tilde{x} = \tilde{b}$ %}
    where
    {% $\tilde{A} = S_1 P_1^{-1} A P_2^{-1} S_2^{-1}$ %},
    {% $\tilde{b} = S_1 P_1^{-1} b$ %}, and
    {% $\tilde{x} = S_2 P_2 x$ %}.

    Note that the underlying data structures may be used directly by the
    linear solver, i.e., without a copy.

    @nocvode <node> SUNLinSolSetScalingVectors *)
val set_scaling_vectors :
     ('m, 'd, 'k, 't) t
  -> ('d, 'k) Nvector.t
  -> ('d, 'k) Nvector.t
  -> unit

(** Indicates that the next call to {!solve} will be made with a zero initial
    guess.

    @since 5.8.0
    @nocvode <node> SUNLinSolSetZeroGuess *)
val set_zero_guess :
     ('m, 'd, 'k, 't) t
  -> bool
  -> unit

(** Initializes a linear solver.

    @nocvode <node> SUNLinSolInitialize *)
val init : ('m, 'd, 'k, 't) t -> unit

(** Instruct the linear solver to prepare to solve using an updated
    system matrix.

    @nocvode <node> SUNLinSolSetup *)
val setup : ('m, 'd, 'k, 't) t -> ('a, 'm, 'd, 'k) Matrix.t -> unit

(** Solve a linear system.
    The call [solve ls a x b tol] solves the linear system
    {% $Ax = b$ %}.

    Direct solvers ignore [tol].
    Matrix-free solvers ignore [a], relying instead on the function passed to
    {!set_atimes}.
    Iterative solvesr attempt to respect the weighted 2-norm tolerance,
    [tol].

    @nocvode <node> SUNLinSolSolve *)
val solve :
     ('m, 'd, 'k, 't) t
  -> ('a, 'm, 'd, 'k) Matrix.t
  -> ('d, 'k) Nvector.t
  -> ('d, 'k) Nvector.t
  -> float
  -> unit

(* The number of linear iterations performed in the last {!solve} call.

    @nocvode <node> SUNLinNumIters *)
val get_num_iters : ('m, 'd, 'k, 't) t -> int

(* The final residual norm from the last {!solve} call.

    @nocvode <node> SUNLinResNorm *)
val get_res_norm : ('m, 'd, 'k, 't) t -> float

(* The preconditioned initial residual vector.
   May be called when an iterative method computes the preconditioned initial
   residual and {!solve} succeeds without performing any iterations, i.e., if
   the intial guess or the preconditioner is sufficiently accurate.

   The linear solver may return a reference to its internal array, i.e., it
   may not return a copy.

   @nocvode <node> SUNLinResNorm *)
val get_res_id : ('m, 'd, 'k, 't) t -> 'd

(** Returns the type of the linear solver.

    @nocvode <node> SUNLinSolGetType *)
val get_type : ('m, 'd, 'k, 't) t -> linear_solver_type

(** Returns the identifier of the linear solver.

    @nocvode <node> SUNLinSolGetID
    @since 5.0.0 *)
val get_id : ('m, 'd, 'k, 't) t -> linear_solver_id

(** Returns an indication of the last error encountered by a linear solver.

    @nocvode <node> SUNLinSolGetLastFlag
    @since 5.0.0 *)
val get_last_flag : ('m, 'd, 'k, 't) t -> int

(** The storage requirements of the linear solver.
    The result [(lrw, liw)] gives the number of words used for
    storing real values ([lrw]) and the number of words used
    for storing integer values ([liw]).

    @nocvode <node> SUNLinSpace *)
val get_work_space : ('m, 'd, 'k, 't) t -> int * int

(** {2:exceptions Exceptions} *)

(** Raised on invalid use of linear solver functions. For instance,
    initializing a session with {!Cvode.Diag} and then calling
    {!Cvode.Spils.get_num_lin_iters}, which rather requires a
    linear solver from {!Iterative}. *)
exception InvalidLinearSolver

(** Raised on an unrecoverable failure in a linear solver. The argument is
    [true] for a recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_PACKAGE_FAIL_REC/_UNREC} *)
exception UnrecoverableFailure of bool

(** Raised when creating a linear solver if the given matrix is not square. *)
exception MatrixNotSquare

(** Raised when creating a linear solver if the number of matrix rows and the
    vector length are not equal. *)
exception MatrixVectorMismatch

(** Raised when the storage upper bandwidth ([smu]) of a
    {{!Sundials_Matrix.Band.t}Matrix.Band.t} is
    insufficient for use in a particular linear solver. *)
exception InsufficientStorageUpperBandwidth

(** Raised on an attempt to associate a linear solver instance with more than
    one session. *)
exception LinearSolverInUse

(** Indicates failure of an atimes function. The argument is [true] for a
    recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_ATIMES_FAIL_REC/_UNREC} *)
exception ATimesFailure of bool

(** Indicates failure of a preconditioner setup routine. The argument is
    [true] for a recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_PSET_FAIL_REC/_UNREC} *)
exception PSetFailure of bool

(** Indicates failure of a preconditioner solver. The argument is [true] for a
    recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_PSOLVE_FAIL_REC/_UNREC} *)
exception PSolveFailure of bool

(** Indicates failure of a Gram-Schmidt routine. {cconst SUNLS_GS_FAIL} *)
exception GSFailure

(** Indicates that the QR solution found a singular result.
    {cconst SUNLS_QRSOL_FAIL} *)
exception QRSolFailure

(** An error occurred in a vector operation.
    {cconst SUNLS_VECTOROP_ERR} *)
exception VectorOpError

(** Indicates that the residual is reduced but without convergence to the
    desired tolerance. {cconst SUNLS_RES_REDUCED} *)
exception ResReduced

(** Indicates that a solver failed to converge. {cconst SUNLS_CONV_FAIL} *)
exception ConvFailure

(** Indicates that QR factorization encountered a singular matrix.
    {cconst SUNLS_QRFACT_FAIL} *)
exception QRfactFailure

(** Indicates that LU factorization encountered a singular matrix.
    {cconst SUNLS_LUFACT_FAIL} *)
exception LUfactFailure

(** Indicates failure in an external linear solver package. The argument
    is [true] for a recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_PACKAGE_FAIL_REC/_UNREC} *)
exception PackageFailure of bool

(** Raised by {!Iterative.set_prec_type} if the given type is not allowed. *)
exception IllegalPrecType

(** Indicates that an internal callback, identified by the first argument,
    returned the given unknown error code. *)
exception InternalFailure of (string * int)

(** Setup failed due to a zero diagonal element during LU factorization.
    The argument indicates the column index numbered from one. *)
exception ZeroInDiagonal of int

