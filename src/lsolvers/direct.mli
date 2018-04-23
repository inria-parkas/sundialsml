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

(** Direct Linear Solvers.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNLinearSolver module
    @since 3.0.0 *)

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
    The matrix is used internally after the linear solver is attached to a
    session.

  @nocvode <node> SUNDenseLinearSolver *)
val dense :
  'k Nvector_serial.any
  -> 'k Matrix.dense
  -> (Matrix.Dense.t, 'k, tag) serial_t

(** Creates a direct linear solver on dense matrices using LAPACK.
    See {!dense}. Only available if {!Sundials.lapack_enabled}.

  @nocvode <node> SUNLapackDense *)
val lapack_dense :
  'k Nvector_serial.any
  -> 'k Matrix.dense
  -> (Matrix.Dense.t, 'k, tag) serial_t

(** Creates a direct linear solver on banded matrices. The nvector and matrix
    argument are used to determine the linear system size and to assess
    compatibility with the linear solver implementation.
    The matrix is used internally after the linear solver is attached to a
    session.

  @nocvode <node> SUNBandLinearSolver *)
val band :
  'k Nvector_serial.any
  -> 'k Matrix.band
  -> (Matrix.Band.t, 'k, tag) serial_t

(** Creates a direct linear solver on banded matrices using LAPACK.
    See {!band}. Only available if {!Sundials.lapack_enabled}.

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
      The matrix is used internally after the linear solver is attached to a
      session.

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

(** Creates a direct linear solver on sparse matrices using KLU.
    See {!Klu.make}.

    @raise Sundials.NotImplementedBySundialsVersion Solver not available.
    @nocvode <node> SUNKLU *)
val klu :
  ?ordering:Klu.ordering
  -> 'k Nvector_serial.any
  -> ('s, 'k) Matrix.sparse
  -> ('s Matrix.Sparse.t, 'k, Klu.tag) serial_t

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
      The matrix is used internally after the linear solver is attached to a
      session.

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

(** Creates a direct linear solver on sparse matrices using SuperLUMT.
    See {!Superlumt.make}.

    @raise Sundials.NotImplementedBySundialsVersion Solver not available.
    @nocvode <node> SUNSuperLUMT *)
val superlumt :
  ?ordering:Superlumt.ordering
  -> nthreads:int
  -> 'k Nvector_serial.any
  -> (Matrix.Sparse.csc, 'k) Matrix.sparse
  -> (Matrix.Sparse.csc Matrix.Sparse.t, 'k, Superlumt.tag) serial_t

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
             -> ('matrixkind, 'matrix, 'data, 'kind) Matrix.t
             -> ('matrix, 'data, 'kind, 'lsolver tag) t

  (** Return the internal state from an custom direct linear solver. *)
  val unwrap : ('matrix, 'data, 'kind, 'lsolver tag) t -> 'lsolver

end (* }}} *)

