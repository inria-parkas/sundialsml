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

    A set of functions for instantiating linear solvers from two families:
    {!module:Direct} and {!module:Iterative}. Any instance may be associated
    with at most one solver session.

    The generic operations defined herein apply to the linear solver types
    {!Direct.t} and {!Iterative.t}. Solver-specific configuration and status
    operations are defined in the different solver modules and applied to
    session types, for example, {!Cvode.Direct.get_work_space} and
    {!Arkode.Iterative.set_eps_lin}.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNLinearSolver module
    @since 3.0.0 *)

(** Direct Linear Solvers. *)
module Direct : sig (* {{{ *)

  (** Used to identify generic direct solvers. *)
  type tag = [`Basic]

  (** A direct linear solver.
      The type variables specify the Jacobian matrix (['matrix]), the
      {!Nvector.nvector} data (['data]) and kind (['kind]), and a
      ['tag] used to identify specific solver features.

      @nocvode <node> Description of the SUNLinearSolver module
      @nocvode <node> SUNLinearSolver *)
  type ('matrix, 'data, 'kind, 'tag) t
    = ('matrix, 'data, 'kind, 'tag) Lsolver_impl.Direct.t

  (** Alias for linear solvers that are restricted to serial nvectors. *)
  type ('mat, 'kind, 'tag) serial_t
    = ('mat, Nvector_serial.data, [>Nvector_serial.kind] as 'kind, 'tag) t

  (** Creates a direct linear solver on dense matrices. The nvector and matrix
      argument are used to determine the linear system size and to assess
      compatibility with the linear solver implementation.

    @nocvode <node> SUNDenseLinearSolver *)
  val dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k, tag) serial_t

  (** Creates a direct linear solver on dense matrices using LAPACK.
      See {!make}. Only available if {!Sundials.lapack_enabled}.

    @nocvode <node> SUNLapackDense *)
  val lapack_dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k, tag) serial_t

  (** Creates a direct linear solver on banded matrices. The nvector and matrix
      argument are used to determine the linear system size and to assess
      compatibility with the linear solver implementation.

    @nocvode <node> SUNBandLinearSolver *)
  val band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k, tag) serial_t

  (** Creates a direct linear solver on banded matrices using LAPACK.
      See {!make}. Only available if {!Sundials.lapack_enabled}.

    @nocvode <node> SUNLapackBand *)
  val lapack_band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k, tag) serial_t

  (** KLU direct linear solver operating on sparse matrices (requires KLU). *)
  module Klu : sig (* {{{ *)

    (** Used to distinguish KLU direct solvers. *)
    type tag = [`Klu]

    (** The ordering algorithm used for reducing fill. *)
    type ordering = Lsolver_impl.Klu.ordering =
         Amd      (** Approximate minimum degree permutation. *)
       | ColAmd   (** Column approximate minimum degree permutation. *)
       | Natural  (** Natural ordering. *)

    (** Creates a direct linear solver on sparse matrices using KLU. The
        nvector and matrix argument are used to determine the linear system
        size and to assess compatibility with the linear solver implementation.

        @raise Sundials.NotImplementedBySundialsVersion Solver not available.
        @nocvode <node> SUNKLU *)
    val make :
      ?ordering:ordering
      -> 'k Nvector_serial.any
      -> ('s, 'k) Matrix.sparse
      -> ('s Matrix.Sparse.t, 'k, tag) serial_t

    (** Reinitializes memory and flags for a new factorization (symbolic and
        numeric) at the next solver setup call. In the call [reinit ls a nnz],
        [a] is the Jacobian matrix, which is reinitialized with the given
        number of non-zeros if [nnz] if given. New symbolic and numeric
        factorizations will be completed at the next solver step.

        @nocvode <node> SUNKLUReInit *)
    val reinit : ('s Matrix.Sparse.t, 'k, [>tag]) serial_t
      -> ('s, 'k) Matrix.sparse -> ?nnz:int -> unit -> unit

    (** Sets the ordering algorithm used to minimize fill-in.

        @nocvode <node> SUNKLUSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k, [>tag]) serial_t
      -> ordering -> unit

  end (* }}} *)

  (** SuperLUMT direct linear solver operating on sparse matrices (requires
      SuperLUMT). *)
  module Superlumt : sig (* {{{ *)

    (** Used to distinguish SuperLUMT direct solvers. *)
    type tag = [`Superlumt]

    (** The ordering algorithm used for reducing fill. *)
    type ordering = Lsolver_impl.Superlumt.ordering =
         Natural       (** Natural ordering. *)
       | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
       | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
       | ColAmd        (** Column approximate minimum degree permutation. *)

    (** Creates a direct linear solver on sparse matrices using SuperLUMT. The
        nvector and matrix argument are used to determine the linear system
        size and to assess compatibility with the linear solver implementation.

        @raise Sundials.NotImplementedBySundialsVersion Solver not available.
        @nocvode <node> SUNSuperLUMT *)
    val make :
      ?ordering:ordering
      -> nthreads:int
      -> 'k Nvector_serial.any
      -> (Matrix.Sparse.csc, 'k) Matrix.sparse
      -> (Matrix.Sparse.csc Matrix.Sparse.t, 'k, tag) serial_t

    (** Sets the ordering algorithm used to minimize fill-in.

        @nocvode <node> SUNSuperLUMTSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k, [>tag]) serial_t
      -> ordering -> unit

  end (* }}} *)

  (** Custom direct linear solvers. *)
  module Custom : sig (* {{{ *)

    (** Used to distinguish custom linear solvers *)
    type 'lsolver tag = [`Custom of 'lsolver]

    (** The operations required to implement a direct linear solver.
        Failure should be indicated by raising an exception (preferably
        one of the exceptions in the {!module:Lsolver} package). Raising
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
               -> ('matrix, 'data, 'kind, 'lsolver tag) t

    (** Return the internal state from an custom direct linear solver. *)
    val unwrap : ('matrix, 'data, 'kind, 'lsolver tag) t -> 'lsolver

  end (* }}} *)

end (* }}} *)

(** Iterative linear Solvers. *)
module Iterative : sig (* {{{ *)

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

      Raises {!IllegalPrecType} if the current preconditioner
      is {{!preconditioning_type}PrecNone} and the given argument is not
      (since no callback functions are specified in this case. May raise
      {!IllegalPrecType} if the given type is not allowed by the underlying
      solver.

      @nocvode <node> SUNPCGSetPrecType
      @nocvode <node> SUNSPBCGSSetPrecType
      @nocvode <node> SUNSPFGMRSetPrecType
      @nocvode <node> SUNSPGMRSetPrecType
      @nocvode <node> SUNSPTFQMRSetPrecType *)
  val set_prec_type : ('d, 'k, 'f) t -> preconditioning_type -> unit

end (* }}} *)

(** {2:exceptions Exceptions} *)

(** Raised on an unrecoverable failure in a linear solver. The argument is
    [true] for a recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_PACKAGE_FAIL_REC/_UNREC} *)
exception UnrecoverableFailure of bool

(** Raised when creating a linear solver if the given matrix is not square. *)
exception MatrixNotSquare

(** Raised on an attempt to associate a linear solver instance with more than
    one session. *)
exception SolverInUse

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

(** Raised by {!set_prec_type} if the given type is not allowed. *)
exception IllegalPrecType

(** Indicates that an internal callback, identified by the first argument,
    returned the given unknown error code. *)
exception InternalFailure of (string * int)

(** Raised if a zero diagonal element is found during factorization using a
    low-level routine like {!Spils.qr_fact}, {!Spils.qr_sol},
    {!Dls.DenseMatrix.getrf}, or {!Dls.ArrayDenseMatrix.getrf}.
    The argument gives the equation number (from 1). *)
exception ZeroDiagonalElement of int

