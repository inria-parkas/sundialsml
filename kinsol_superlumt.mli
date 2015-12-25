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

(** SuperLU_MT sparse-direct linear solver module for KINSOL
    (requires SuperLU_MT).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @nokinsol <node5#sss:kinsuperlumt> The SuperLUMT Solver *)

(** Callback functions that compute sparse approximations to a Jacobian
    matrix. In the call [sparse_jac_fn arg jac], [arg] is a {!jacobian_arg}
    with two work vectors and the computed Jacobian must be stored
    in [jac].

    The callback should load the [(i,j)]th entry of [jac] with
    {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of the
    [i]th equation with respect to the [j]th variable, evaluated at the
    values of [t] and [y] obtained from [arg]. Only nonzero elements need
    be loaded into [jac].

    Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
    Any other exception is treated as an unrecoverable error.

    {warning Neither the elements of [arg] nor the matrix [jac] should
             be accessed after the function has returned.}

    @nokinsol <node5#ss:sjacFn> KINSlsSparseJacFn *)
type sparse_jac_fn =
  (Sundials.RealArray.t Kinsol.double, Sundials.RealArray.t) Kinsol.jacobian_arg
  -> Sls.SparseMatrix.t -> unit

(** A direct linear solver on sparse matrices. In the call,
    [superlumt jfn nnz nthreads], [jfn] is a callback function that computes
    an approximation to the Jacobian matrix, [nnz] is the maximum number of
    nonzero entries in that matrix, and [nthreads] is the number of threads
    to use when factorizing/solving.

    @nokinsol <node5#sss:lin_solv_init> KINSuperLUMT
    @nokinsol <node5#sss:optin_sls> KINSlsSetSparseJacFn
    @nokinsol <node5#ss:sjacFn> KINSlsSparseJacFn *)
val superlumt
      : sparse_jac_fn -> nnz:int -> nthreads:int -> Kinsol.serial_linear_solver

(** The ordering algorithm used for reducing fill. *)
type ordering =
     Natural       (** Natural ordering. *)
   | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
   | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
   | ColAmd        (** Column approximate minimum degree permutation. *)

(** Sets the ordering algorithm used to minimize fill-in.

    @nokinsol <node5#ss:sls_optin> KINSuperLUMTSetOrdering *)
val set_ordering : Kinsol.serial_session -> ordering -> unit

(** Returns the number of calls made by a sparse linear solver to the
    Jacobian approximation function.

    @nokinsol <node5#sss:optout_sls> KINSlsGetNumJacEvals *)
val get_num_jac_evals : Kinsol.serial_session -> int

