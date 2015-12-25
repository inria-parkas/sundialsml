(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** SuperLU_MT sparse-direct linear solver module for ARKODE
    (requires SuperLU_MT).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @noarkode <node> Linear solver specification functions *)

(** Callback functions that compute sparse approximations to a Jacobian
    matrix. In the call [sparse_jac_fn arg jac], [arg] is a {!jacobian_arg}
    with three work vectors and the computed Jacobian must be stored
    in [jac].

    The callback should load the [(i,j)]th entry of [jac] with
    {% $\partial f_i/\partial y_j$%}, i.e., the partial derivative of the
    [i]th equation with respect to the [j]th variable, evaluated at the
    values of [t] and [y] obtained from [arg]. Only nonzero elements need
    be loaded into [jac].

    Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
    Any other exception is treated as an unrecoverable error.

    {warning Neither the elements of [arg] nor the matrix [jac] should
             be accessed after the function has returned.}

    @noarkode <node5#ss:sjacFn> ARKSlsSparseJacFn *)
type sparse_jac_fn =
  (Sundials.RealArray.t Arkode.triple, Sundials.RealArray.t) Arkode.jacobian_arg
  -> Sls.SparseMatrix.t -> unit

(** A direct linear solver on sparse matrices. In the call,
    [superlumt jfn ~nnz ~nthreads], [jfn] is a callback function that computes
    an approximation to the Jacobian matrix, [nnz] is the maximum number of
    nonzero entries in that matrix, and [nthreads] is the number of threads
    to use when factorizing/solving.

    @noarkode <node5#sss:lin_solv_init> ARKSuperLUMT
    @noarkode <node5#sss:optin_sls> ARKSlsSetSparseJacFn
    @noarkode <node5#ss:sjacFn> ARKSlsSparseJacFn *)
val superlumt
      : sparse_jac_fn -> nnz:int -> nthreads:int -> Arkode.serial_linear_solver

(** The ordering algorithm used for reducing fill. *)
type ordering =
     Natural       (** Natural ordering. *)
   | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
   | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
   | ColAmd        (** Column approximate minimum degree permutation. *)

(** Sets the ordering algorithm used to minimize fill-in.

    @noarkode <node5#ss:sls_optin> ARKSuperLUMTSetOrdering *)
val set_ordering : Arkode.serial_session -> ordering -> unit

(** Returns the number of calls made by a sparse linear solver to the
    Jacobian approximation function.

    @noarkode <node5#sss:optout_sls> ARKSlsGetNumJacEvals *)
val get_num_jac_evals : Arkode.serial_session -> int

module Mass : sig
  (** Callback functions that compute sparse approximations to the mass
      matrix. In the call [sparse_fn t tmp m],
      - [t] is the independent variable,
      - [tmp] is workspace data, and
      - [m] is the output sparse mass matrix.

      The callback should load the compressed-sparse-column matrix [m] with an
      approximation to the mass matrix {% $M(t)$%}.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the elements of [tmp] nor the matrix [m] should
               be accessed after the function has returned.}

      @noarkode <node5#ss:sjacFn> ARKSlsSparseMassFn *)
  type sparse_fn
    = float -> Sundials.RealArray.t Arkode.triple -> Sls.SparseMatrix.t -> unit

  (** A direct linear solver on sparse matrices. In the call,
      [superlumt mfn ~nnz ~nthreads], [mfn] is a callback function that computes
      an approximation to the mass matrix, [nnz] is the maximum number of
      nonzero entries in that matrix, and [nthreads] is the number of threads
      to use when factorizing/solving.

      @noarkode <node5#sss:lin_solv_init> ARKMassSuperLUMT
      @noarkode <node5#ss:smassFn> ARKSlsSparseMassFn *)
  val superlumt
      : sparse_fn -> nnz:int -> nthreads:int -> Arkode.serial_linear_solver

  (** Sets the ordering algorithm used to minimize fill-in.

      @noarkode <node5#ss:sls_optin> ARKMassSuperLUMTSetOrdering *)
  val set_ordering : Arkode.serial_session -> ordering -> unit

  (** Returns the number of calls made by a sparse linear solver to the
      mass matrix approximation function.

      @noarkode <node5#sss:optout_sls> ARKSlsGetNumMassEvals *)
  val get_num_evals : Arkode.serial_session -> int
end

