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

(** KLU sparse-direct linear solver module for CVODE (requires KLU).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @nocvode <node5#sss:cvklu> The KLU Solver *)

(** Callback functions that compute sparse approximations to a Jacobian
    matrix. In the call [sparse_jac_fn arg jac], [arg] is a {!jacobian_arg}
    with three work vectors and the computed Jacobian must be stored
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

    @nocvode <node5#ss:sjacFn> CVSlsSparseJacFn *)
type sparse_jac_fn =
  (Sundials.RealArray.t Cvode.triple, Sundials.RealArray.t) Cvode.jacobian_arg
  -> Sls.SparseMatrix.t -> unit

(** A direct linear solver on sparse matrices. In the call,
    [klu jfn nnz], [jfn] specifies a callback function that computes an
    approximation to the Jacobian matrix and [nnz] is the maximum number
    of nonzero entries in that matrix.

    @nocvode <node5#sss:lin_solv_init> CVKLU
    @nocvode <node5#sss:optin_sls> CVSlsSetSparseJacFn
    @nocvode <node5#ss:sjacFn> CVSlsSparseJacFn *)
val klu : sparse_jac_fn -> int -> Cvode.serial_linear_solver

(** The ordering algorithm used for reducing fill. *)
type ordering =
     Amd      (** Approximate minimum degree permutation. *)
   | ColAmd   (** Column approximate minimum degree permutation. *)
   | Natural  (** Natural ordering. *)

(** Sets the ordering algorithm used to minimize fill-in.

    @nocvode <node5#ss:sls_optin> CVKLUSetOrdering *)
val set_ordering : Cvode.serial_session -> ordering -> unit

(** Reinitializes the Jacobian matrix memory and flags.
    In the call, [reinit s n nnz realloc], [n] is the number of system state
    variables, and [nnz] is the number of non-zeroes in the Jacobian matrix.
    New symbolic and numeric factorizations will be completed at the next solver
    step. If [realloc] is true, the Jacobian matrix will be reallocated based on
    [nnz].

    @nocvode <node5#ss:sls_optin> CVKLUReInit *)
val reinit : Cvode.serial_session -> int -> int -> bool -> unit

