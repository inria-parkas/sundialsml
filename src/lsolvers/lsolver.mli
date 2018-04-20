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

(** Generic linear solvers: shared definitions.

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

