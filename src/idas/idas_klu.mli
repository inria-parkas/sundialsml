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

(** KLU sparse-direct linear solver module for IDAS (requires KLU).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @noidas <node5#sss:cvklu> The KLU Solver *)

(** Callback functions that compute sparse approximations to a Jacobian
    matrix without forward sensitivites. In the call [sparse_jac_fn arg jac],
    [arg] is a {!Idas.Adjoint.jacobian_arg} with three work vectors and the
    computed Jacobian must be stored in [jac].

    The callback should load the [(i,j)]th entry of [jac] with
    {% $\frac{\partial F_i}{\partial y_j} + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%},
    i.e., the partial derivative of the [i]th equation with respect to
    the [j]th variable, evaluated at the values of [t], [y], and [y']
    obtained from [arg]. Only nonzero elements need be loaded into [jac].

    Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
    Any other exception is treated as an unrecoverable error.

    {warning Neither the elements of [arg] nor the matrix [jac] should
             be accessed after the function has returned.}

    @noidas <node5#ss:sjacFnB> IDASlsSparseJacFnB *)
type sparse_jac_fn_no_sens =
  (Sundials.RealArray.t Ida.triple, Sundials.RealArray.t)
      Idas.Adjoint.jacobian_arg
  -> Sls.SparseMatrix.t -> unit

(** Callback functions that compute sparse approximations to a Jacobian
    matrix with forward sensitivities. In the call [sparse_jac_fn arg s s' jac],
    [arg] is a {!Idas.Adjoint.jacobian_arg} with three work vectors, [s] is an
    array of forward sensitivity vectors, [s'] is an array of forward
    sensitivity derivative vectors and the computed Jacobian must be stored in
    [jac].

    The callback should load the [(i,j)]th entry of [jac] with
    {% $\frac{\partial F_i}{\partial y_j} + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%},
    i.e., the partial derivative of the [i]th equation with respect to
    the [j]th variable, evaluated at the values of [t], [y], and [y']
    obtained from [arg]. Only nonzero elements need be loaded into [jac].

    Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
    Any other exception is treated as an unrecoverable error.

    {warning Neither the elements of [arg] nor the matrix [jac] should
             be accessed after the function has returned.}

    @noidas <node5#ss:sjacFnBS> IDASlsSparseJacFnBS *)
type sparse_jac_fn_with_sens =
  (Sundials.RealArray.t Ida.triple, Sundials.RealArray.t)
      Idas.Adjoint.jacobian_arg
  -> Sundials.RealArray.t array
  -> Sundials.RealArray.t array
  -> Sls.SparseMatrix.t -> unit

(** Callback functions that compute sparse approximations to a Jacobian
    matrix.

    @noidas <node5#ss:sjacFnB> IDASlsSparseJacFnB
    @noidas <node5#ss:sjacFnBS> IDASlsSparseJacFnBS *)
type sparse_jac_fn =
    NoSens of sparse_jac_fn_no_sens
    (** Does not depend on forward sensitivities. *)
  | WithSens of sparse_jac_fn_with_sens
    (** Depends on forward sensitivities. *)

(** A direct linear solver on sparse matrices. In the call,
    [klu jfn nnz], [jfn] is a callback function that computes an
    approximation to the Jacobian matrix and [nnz] is the maximum number
    of nonzero entries in that matrix.

    @noidas <node5#sss:lin_solv_init> IDAKLUB
    @noidas <node5#sss:optin_sls> IDASlsSetSparseJacFnB
    @noidas <node5#sss:optin_sls> IDASlsSetSparseJacFnBS
    @noidas <node5#ss:sjacFnB> IDASlsSparseJacFnB
    @noidas <node5#ss:sjacFnBS> IDASlsSparseJacFnBS *)
val klu : sparse_jac_fn -> int -> 'k Idas.Adjoint.serial_linear_solver

(** The ordering algorithm used for reducing fill. *)
type ordering = Ida_klu.ordering =
     Amd      (** Approximate minimum degree permutation. *)
   | ColAmd   (** Column approximate minimum degree permutation. *)
   | Natural  (** Natural ordering. *)

(** Sets the ordering algorithm used to minimize fill-in.

    @noidas <node5#ss:sls_optin> IDAKLUSetOrdering
    @idas <node7#ss:optional_output_b> IDAGetAdjIdaBmem *)
val set_ordering : 'k Idas.Adjoint.serial_bsession -> ordering -> unit

(** Reinitializes the Jacobian matrix memory and flags.
    In the call, [reinit s n nnz realloc], [n] is the number of system state
    variables, and [nnz] is the number of non-zeroes in the Jacobian matrix.
    New symbolic and numeric factorizations will be completed at the next solver
    step. If [realloc] is true, the Jacobian matrix will be reallocated based on
    [nnz].

    @noidas <node5#ss:sls_optin> IDAKLUReInit
    @idas <node7#ss:optional_output_b> IDAGetAdjIdaBmem *)
val reinit : 'k Idas.Adjoint.serial_bsession -> int -> int -> bool -> unit

(** Returns the number of calls made by a sparse linear solver to the
    Jacobian approximation function.

    @noidas <node5#sss:optout_sls> IDASlsGetNumJacEvals
    @idas <node7#ss:optional_output_b> IDAGetAdjIdaBmem *)
val get_num_jac_evals : 'k Idas.Adjoint.serial_bsession -> int

