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

(** {2:exceptions Exceptions} *)

(** Raised on an unrecoverable failure in a linear solver. The argument is
    [true] for a recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_PACKAGE_FAIL_REC/_UNREC} *)
exception UnrecoverableFailure of bool

(** Raised when creating a linear solver if the given matrix is not square. *)
exception MatrixNotSquare

(** Raised when creating a linear solver if the number of matrix rows and the
    vector length are not equal. *)
exception MatrixVectorMismatch

(** Raised when the storage upper bandwidth ([smu]) of a {!Band.t} is
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

(** {2:iterative Iterative Linear Solvers} *)

module Iterative : sig
  (** Definitions in this module are more conveniently accessed
      through session-specific iterative linear solver modules like
      {!Cvode.Spils} and {!Ida.Spils}.  For example,
      {!Cvode.Spils.spbcgs} is an alias for
      {!Lsolver.Iterative.spbcgs}.  *)

  (** {3:solvers Types} *)

  (** An iterative linear solver.
    The type variables specify the {!Nvector.nvector} data (['data]) and
    kind (['kind]), and the iterative method (['iter]).

    A linear solver of this type must be converted to session-specific
    form by {!Cvode.Spils.solver}, {!Ida.Spils.solver}, etc. before
    being attached to a session.


    @nocvode <node> Description of the SUNLinearSolver module
    @nocvode <node> SUNLinearSolver *)
  type ('data, 'kind, 'iter) linear_solver
    = ('data, 'kind, 'iter) LinearSolver_impl.Iterative.linear_solver

  (** The type of Gram-Schmidt orthogonalization in iterative linear solvers.

    @nocvode <node> ModifiedGS/ClassicalGS *)
  type gramschmidt_type = LinearSolver_impl.Iterative.gramschmidt_type =
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
  val spbcgs : ?maxl:int -> ('d, 'k) Nvector.t
               -> ('d, 'k, [`Spbcgs]) linear_solver

  (** Krylov iterative solver using the scaled preconditioned flexible
    generalized minimum residual (GMRES) method. The [maxl] arguments gives
    the maximum dimension of the Krylov subspace (defaults to 5). The
    nvector argument is used as a template.

    NB: [max_restarts] is ignored in Sundials < 3.0.0.

    @nocvode <node> SUNSPFGMR *)
  val spfgmr : ?maxl:int -> ?max_restarts:int
               -> ?gs_type:gramschmidt_type
               -> ('d, 'k) Nvector.t
               -> ('d, 'k, [`Spfgmr]) linear_solver

  (** Krylov iterative solver using the scaled preconditioned generalized
    minimum residual (GMRES) method. The [maxl] arguments gives the maximum
    dimension of the Krylov subspace (defaults to 5). The nvector argument
    is used as a template.

    NB: [max_restarts] is ignored in Sundials < 3.0.0.

    @nocvode <node> SUNSPGMR *)
  val spgmr : ?maxl:int -> ?max_restarts:int
              -> ?gs_type:gramschmidt_type
              -> ('d, 'k) Nvector.t
              -> ('d, 'k, [`Spgmr]) linear_solver

  (** Krylov iterative with the scaled preconditioned transpose-free
    quasi-minimal residual (SPTFQMR) method. The [maxl] arguments gives the
    maximum dimension of the Krylov subspace (defaults to 5). The nvector
    argument is used as a template.

    @nocvode <node> SUNSPTFQMR *)
  val sptfqmr : ?maxl:int -> ('d, 'k) Nvector.t
                -> ('d, 'k, [`Sptfqmr]) linear_solver

  (** Krylov iterative solver using the preconditioned conjugate gradient
    (PCG) method. The [maxl] arguments gives the maximum dimension of the
    Krylov subspace (defaults to 5). The nvector argument is used as a
    template.

    @nocvode <node> SUNPCG *)
  val pcg : ?maxl:int -> ('d, 'k) Nvector.t
            -> ('d, 'k, [`Pcg]) linear_solver

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
      one of the exceptions in the {!module:LinearSolver} package). Raising
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
               -> ('data, 'kind, 'lsolver tag) linear_solver

    (** Return the internal state from an custom iterative linear solver. *)
    val unwrap : ('data, 'kind, 'lsolver tag) linear_solver -> 'lsolver

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
  val set_maxl : ('d, 'k, [< `Spbcgs|`Sptfqmr|`Pcg]) linear_solver
                 -> int -> unit

  (** Sets the Gram-Schmidt orthogonalization to use.

    @nocvode <node> SUNSPGMRSetGSType
    @nocvode <node> SUNSPFGMRSetGSType *)
  val set_gs_type :
    ('d, 'k, [< `Spfgmr|`Spgmr]) linear_solver -> gramschmidt_type -> unit

  (** Sets the number of GMRES restarts to allow.

    NB: This feature is not supported by Sundials < 3.0.0.

    @nocvode <node> SUNSPGMRSetMaxRestarts
    @nocvode <node> SUNSPFGMRSetMaxRestarts *)
  val set_max_restarts : ('d, 'k, [< `Spfgmr|`Spgmr]) linear_solver
                         -> int -> unit

  (** The type of preconditioning in Krylov solvers.

    @nocvode <node> Preconditioning *)
  type preconditioning_type = LinearSolver_impl.Iterative.preconditioning_type =
    | PrecNone    (** No preconditioning *)
    | PrecLeft    (** {% $(P^{-1}A)x = P^{-1}b$ %} *)
    | PrecRight   (** {% $(AP^{-1})Px = b$ %} *)
    | PrecBoth    (** {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)

  (** Change the preconditioning direction without modifying callback
    functions.

    Raises {!LinearSolver.IllegalPrecType} if the current preconditioner
    is {{!preconditioning_type}PrecNone} and the given argument is not
    (since no callback functions are specified in this case. May raise
    {!LinearSolver.IllegalPrecType} if the given type is not allowed by the
    underlying solver.

    @nocvode <node> SUNPCGSetPrecType
    @nocvode <node> SUNSPBCGSSetPrecType
    @nocvode <node> SUNSPFGMRSetPrecType
    @nocvode <node> SUNSPGMRSetPrecType
    @nocvode <node> SUNSPTFQMRSetPrecType *)
  val set_prec_type : ('d, 'k, 'f) linear_solver -> preconditioning_type -> unit

end

(** {2:direct Direct Linear Solvers} *)

module Direct : sig
  (** Definitions in this module are more conveniently accessed
      through session-specific direct linear solver modules like
      {!Cvode.Dls} and {!Ida.Dls}.  For example, {!Cvode.Dls.dense} is
      an alias for {!Lsolver.Direct.dense}.  *)

  (** {3:solvers Types} *)

  (** Used to identify generic direct solvers. *)
  type tag = [`Basic]

  (** A generic direct linear solver.
    The type variables specify the Jacobian matrix (['matrix]), the
    {!Nvector.nvector} data (['data]) and kind (['kind]), and a
    ['tag] used to identify specific solver features.

    A linear solver of this type must be converted to session-specific
    form by {!Cvode.Dls.solver}, {!Ida.Dls.solver}, etc., before being
    attached to a session via [init] or [reinit].

    @nocvode <node> Description of the SUNLinearSolver module
    @nocvode <node> SUNLinearSolver *)
  type ('matrix, 'data, 'kind, 'tag) linear_solver
    = ('matrix, 'data, 'kind, 'tag) LinearSolver_impl.Direct.linear_solver

  (** Alias for linear solvers that are restricted to serial nvectors. *)
  type ('mat, 'kind, 'tag) serial_linear_solver
    = ('mat, Nvector_serial.data, [>Nvector_serial.kind] as 'kind, 'tag)
        linear_solver

  (** Creates a direct linear solver on dense matrices. The nvector and matrix
    argument are used to determine the linear system size and to assess
    compatibility with the linear solver implementation.
    The matrix is used internally after the linear solver is attached to a
    session.

  @nocvode <node> SUNDenseLinearSolver *)
  val dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k, tag) serial_linear_solver

  (** Creates a direct linear solver on dense matrices using LAPACK.
    See {!dense}. Only available if {!Sundials.lapack_enabled}.

  @nocvode <node> SUNLapackDense *)
  val lapack_dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k, tag) serial_linear_solver

  (** Creates a direct linear solver on banded matrices. The nvector and matrix
    argument are used to determine the linear system size and to assess
    compatibility with the linear solver implementation.
    The matrix is used internally after the linear solver is attached to a
    session.

  @nocvode <node> SUNBandLinearSolver *)
  val band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k, tag) serial_linear_solver

  (** Creates a direct linear solver on banded matrices using LAPACK.
    See {!band}. Only available if {!Sundials.lapack_enabled}.

  @nocvode <node> SUNLapackBand *)
  val lapack_band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k, tag) serial_linear_solver

  (** KLU direct linear solver operating on sparse matrices (requires KLU). *)
  module Klu : sig (* {{{ *)

    (** Used to distinguish KLU direct solvers. *)
    type tag = [`Klu]

    (** The ordering algorithm used for reducing fill. *)
    type ordering = LinearSolver_impl.Klu.ordering =
      | Amd      (** Approximate minimum degree permutation. *)
      | ColAmd   (** Column approximate minimum degree permutation. *)
      | Natural  (** Natural ordering. *)

    (** Creates a direct linear solver on sparse matrices using KLU. The
      nvector and matrix argument are used to determine the linear system
      size and to assess compatibility with the linear solver implementation.
      The matrix is used internally after the linear solver is attached to a
      session.

      @raise Sundials.NotImplementedBySundialsVersion Solver not available.
      @nocvode <node> SUNKLU *)
    val make :
      ?ordering:ordering
      -> 'k Nvector_serial.any
      -> ('s, 'k) Matrix.sparse
      -> ('s Matrix.Sparse.t, 'k, tag) serial_linear_solver

    (** Reinitializes memory and flags for a new factorization (symbolic and
      numeric) at the next solver setup call. In the call [reinit ls a nnz],
      [a] is the Jacobian matrix, which is reinitialized with the given
      number of non-zeros if [nnz] if given. New symbolic and numeric
      factorizations will be completed at the next solver step.

      @nocvode <node> SUNKLUReInit *)
    val reinit : ('s Matrix.Sparse.t, 'k, [>tag]) serial_linear_solver
                 -> ('s, 'k) Matrix.sparse -> ?nnz:int -> unit -> unit

    (** Sets the ordering algorithm used to minimize fill-in.

      @nocvode <node> SUNKLUSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k, [>tag]) serial_linear_solver
                       -> ordering -> unit

  end (* }}} *)

  (** Creates a direct linear solver on sparse matrices using KLU.
    See {!Klu.make}.

    @raise Sundials.NotImplementedBySundialsVersion Solver not available.
    @nocvode <node> SUNKLU *)
  val klu :
    ?ordering:Klu.ordering
    -> 'k Nvector_serial.any
    -> ('s, 'k) Matrix.sparse
    -> ('s Matrix.Sparse.t, 'k, Klu.tag) serial_linear_solver

  (** SuperLUMT direct linear solver operating on sparse matrices (requires
    SuperLUMT). *)
  module Superlumt : sig (* {{{ *)

    (** Used to distinguish SuperLUMT direct solvers. *)
    type tag = [`Superlumt]

    (** The ordering algorithm used for reducing fill. *)
    type ordering = LinearSolver_impl.Superlumt.ordering =
      | Natural       (** Natural ordering. *)
      | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
      | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
      | ColAmd        (** Column approximate minimum degree permutation. *)

    (** Creates a direct linear solver on sparse matrices using SuperLUMT. The
      nvector and matrix argument are used to determine the linear system
      size and to assess compatibility with the linear solver implementation.
      The matrix is used internally after the linear solver is attached to a
      session.

      @raise Sundials.NotImplementedBySundialsVersion Solver not available.
      @nocvode <node> SUNSuperLUMT *)
    val make :
      ?ordering:ordering
      -> nthreads:int
      -> 'k Nvector_serial.any
      -> (Matrix.Sparse.csc, 'k) Matrix.sparse
      -> (Matrix.Sparse.csc Matrix.Sparse.t, 'k, tag) serial_linear_solver

    (** Sets the ordering algorithm used to minimize fill-in.

      @nocvode <node> SUNSuperLUMTSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k, [>tag]) serial_linear_solver
                       -> ordering -> unit

  end (* }}} *)

  (** Creates a direct linear solver on sparse matrices using SuperLUMT.
    See {!Superlumt.make}.

    @raise Sundials.NotImplementedBySundialsVersion Solver not available.
    @nocvode <node> SUNSuperLUMT *)
  val superlumt :
    ?ordering:Superlumt.ordering
    -> nthreads:int
    -> 'k Nvector_serial.any
    -> (Matrix.Sparse.csc, 'k) Matrix.sparse
    -> (Matrix.Sparse.csc Matrix.Sparse.t, 'k, Superlumt.tag)
         serial_linear_solver

  (** Custom direct linear solvers. *)
  module Custom : sig (* {{{ *)

    (** Used to distinguish custom linear solvers *)
    type 'lsolver tag = [`Custom of 'lsolver]

    (** The operations required to implement a direct linear solver.
      Failure should be indicated by raising an exception (preferably
      one of the exceptions in the {!module:LinearSolver} package). Raising
      {!exception:Sundials.RecoverableFailure} indicates a generic
      recoverable failure. *)
    type ('matrix, 'data, 'kind, 'lsolver) ops = {
        init : 'lsolver -> unit;
        (** Performs linear solver initalization. *)

        setup : 'lsolver -> 'matrix -> unit;
        (** Performs linear solver setup based on an updated matrix. *)

        solve : 'lsolver
                -> 'matrix
                -> ('data, 'kind) Nvector.t
                -> ('data, 'kind) Nvector.t
                -> unit;
        (** The call [solve ls A x b] should solve the linear system
        {% $Ax = b$ %}. *)

        get_work_space : ('lsolver -> int * int) option;
        (** Return the storage requirements for the linear solver.
        The result [(lrw, liw)] gives the number of words used for
        storing real values ([lrw]) and the number of words used
        for storing integer values ([liw]). *)
      }

    (** Create a direct linear solver given a set of operations and an
      internal state.

      NB: This feature is only available for Sundials >= 3.0.0. *)
    val make : ('matrix, 'data, 'kind, 'lsolver) ops
               -> 'lsolver
               -> ('matrixkind, 'matrix, 'data, 'kind) Matrix.t
               -> ('matrix, 'data, 'kind, 'lsolver tag) linear_solver

    (** Return the internal state from an custom direct linear solver. *)
    val unwrap : ('matrix, 'data, 'kind, 'lsolver tag) linear_solver -> 'lsolver

  end (* }}} *)
end
