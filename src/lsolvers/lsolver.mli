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

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNLinearSolver module
    @since 3.0.0 *)

(* TODO: Write up notes on how the linear solvers work and explain the
   important differences before and after Sundials 3.0.0. *)

(** Direct Linear Solvers. *)
module Direct : sig (* {{{ *)
  (** A direct linear solver.
      The type variables specify the Jacobian matrix (['matrix]) and the
      {!Nvector.nvector} data (['data]) and kind (['kind]).

      @nocvode <node> Description of the SUNLinearSolver module
      @nocvode <node> SUNLinearSolver *)
  type ('matrix, 'data, 'kind) t
    = ('matrix, 'data, 'kind) Lsolver_impl.Direct.t

  (** Alias for linear solvers that are restricted to serial nvectors. *)
  type ('mat, 'kind) serial_t
    = ('mat, Nvector_serial.data, [>Nvector_serial.kind] as 'kind) t

  (** Creates a direct linear solver on dense matrices. The nvector and matrix
      argument are used to determine the linear system size and to assess
      compatibility with the linear solver implementation.

    @nocvode <node> SUNDenseLinearSolver *)
  val dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k) serial_t

  (** Creates a direct linear solver on dense matrices using LAPACK.
      See {!make}. Only available if {!Sundials.lapack_enabled}.

    @nocvode <node> SUNLapackDense *)
  val lapack_dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k) serial_t

  (** Creates a direct linear solver on banded matrices. The nvector and matrix
      argument are used to determine the linear system size and to assess
      compatibility with the linear solver implementation.

    @nocvode <node> SUNBandLinearSolver *)
  val band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k) serial_t

  (** Creates a direct linear solver on banded matrices using LAPACK.
      See {!make}. Only available if {!Sundials.lapack_enabled}.

    @nocvode <node> SUNLapackBand *)
  val lapack_band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k) serial_t

  (** KLU direct linear solver operating on sparse matrices (requires KLU). *)
  module Klu : sig
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
      -> ('s Matrix.Sparse.t, 'k) serial_t

    (** Reinitializes memory and flags for a new factorization (symbolic and
        numeric) at the next solver setup call. In the call [reinit ls a nnz],
        [a] is the Jacobian matrix, which is reinitialized with the given
        number of non-zeros if [nnz] if given. New symbolic and numeric
        factorizations will be completed at the next solver step.

        @nocvode <node> SUNKLUReInit *)
    val reinit : ('s Matrix.Sparse.t, 'k) serial_t
      -> ('s, 'k) Matrix.sparse -> ?nnz:int -> unit -> unit

    (** Sets the ordering algorithm used to minimize fill-in.

        @nocvode <node> SUNKLUSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k) serial_t -> ordering -> unit
  end

  (** SuperLUMT direct linear solver operating on sparse matrices (requires
      SuperLUMT). *)
  module Superlumt : sig
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
      -> (Matrix.csc, 'k) Matrix.sparse
      -> (Matrix.csc Matrix.Sparse.t, 'k) serial_t

    (** Sets the ordering algorithm used to minimize fill-in.

        @nocvode <node> SUNSuperLUMTSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k) serial_t -> ordering -> unit
  end

  (* TODO: Add custom direct lsolver (subset of generic functions). *)
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
  type gramschmidt_type = Spils.gramschmidt_type =
    | ModifiedGS   (** Modified Gram-Schmidt orthogonalization
                       {cconst MODIFIED_GS} *)
    | ClassicalGS  (** Classical Gram Schmidt orthogonalization
                       {cconst CLASSICAL_GS} *)

  (** {3:iter_solvers Solvers} *)

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

  (* TODO: Add custom spils lsolver (subset of generic functions). *)

  (** {3:iter_param Solver parameters} *)

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
  type preconditioning_type = Spils.preconditioning_type =
    | PrecNone    (** No preconditioning *)
    | PrecLeft    (** {% $(P^{-1}A)x = P^{-1}b$ %} *)
    | PrecRight   (** {% $(AP^{-1})Px = b$ %} *)
    | PrecBoth    (** {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)

  (** Change the preconditioning direction without modifying callback
      functions.

      Raises {IllegalPrecType} if the current preconditioner
      is {{!preconditioning_type}PrecNone} and the given argument is not
      (since no callback functions are specified in this case. May raise
      {IllegalPrecType} if the given type is not allowed by the underlying
      solver.

      @nocvode <node> SUNPCGSetPrecType
      @nocvode <node> SUNSPBCGSSetPrecType
      @nocvode <node> SUNSPFGMRSetPrecType
      @nocvode <node> SUNSPGMRSetPrecType
      @nocvode <node> SUNSPTFQMRSetPrecType *)
  val set_prec_type : ('d, 'k, 'f) t -> preconditioning_type -> unit

  (** {2 Exceptions} *)

  (** Raised when an atimes function fails. The argument is [true] for a
      recoverable failure and [false] for an unrecoverable one.
      {cconst SUNLS_ATIMES_FAIL_REC/_UNREC} *)
  exception ATimesFailure of bool

  (** Raised when a preconditioner setup routine fails. The argument is [true]
      for a recoverable failure and [false] for an unrecoverable one.
      {cconst SUNLS_PSET_FAIL_REC/_UNREC} *)
  exception PSetFailure of bool

  (** Raised when a preconditioner solver fails. The argument is [true] for a
      recoverable failure and [false] for an unrecoverable one.
      {cconst SUNLS_PSOLVE_FAIL_REC/_UNREC} *)
  exception PSolveFailure of bool

  (** Raised when a Gram-Schmidt routine fails. {cconst SUNLS_GS_FAIL} *)
  exception GSFailure

  (** Raised QR solution finds a singular result. {cconst SUNLS_QRSOL_FAIL} *)
  exception QRSolFailure

  (** Raised if the residual is reduced but without convergence to the desired
      tolerance. {cconst SUNLS_RES_REDUCED} *)
  exception ResReduced

  (** Raised when a solver fails to converge.
      {cconst SUNLS_CONV_FAIL} *)
  exception ConvFailure

  (** Raised when QR factorization encounters a singular matrix.
      {cconst SUNLS_QRFACT_FAIL} *)
  exception QRfactFailure

  (** Raised when LU factorization encounters a singular matrix.
      {cconst SUNLS_LUFACT_FAIL} *)
  exception LUfactFailure

  (** Raised by {!set_prec_type} if the given type is not allowed. *)
  exception IllegalPrecType

end (* }}} *)

(* TODO: Add custom linear solvers
   - custom direct lsolver (subset of generic functions)
   - custom spils lsolver (subset of generic functions)
   - custom lsolver (four functions to each solver as before)
 *)

(** {2 Exceptions} *)

(** Raised on an unrecoverable failure in a linear solver. The argument is
    [true] for a recoverable failure and [false] for an unrecoverable one.
    {cconst SUNLS_PACKAGE_FAIL_REC/_UNREC} *)
exception UnrecoverableFailure of bool

(** Raised when creating a linear solver if the given matrix is not square. *)
exception MatrixNotSquare

