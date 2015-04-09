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

(** SuperLU_MT sparse-direct linear solver module for CVODES
    (requires SuperLU_MT).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @nocvodes <node5#sss:cvsuperlumt> The SuperLUMT Solver *)

(** Callback functions that compute sparse approximations to a Jacobian
    matrix without forward sensitivites. In the call [sparse_jac_fn arg jac],
    [arg] is a {!jacobian_arg} with three work vectors and the computed
    Jacobian must be stored in [jac].

    The callback should load the [(i,j)]th entry of [jac] with
    {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of the
    [i]th equation with respect to the [j]th variable, evaluated at the
    values of [t] and [y] obtained from [arg]. Only nonzero elements need
    be loaded into [jac].

    Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
    Any other exception is treated as an unrecoverable error.

    {warning Neither the elements of [arg] nor the matrix [jac] should
             be accessed after the function has returned.}

    @nocvodes <node5#ss:sjacFnB> CVSlsSparseJacFnB *)
type sparse_jac_fn_no_sens =
  (Sundials.RealArray.t Cvode.triple, Sundials.RealArray.t)
      Cvodes.Adjoint.jacobian_arg
  -> Sls.SparseMatrix.t -> unit

(** Callback functions that compute sparse approximations to a Jacobian
    matrix with forward sensitivities. In the call [sparse_jac_fn arg s jac],
    [arg] is a {!jacobian_arg} with three work vectors, [s] is an array
    holding the forward sensitivity vectors, and the computed
    Jacobian must be stored in [jac].

    The callback should load the [(i,j)]th entry of [jac] with
    {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of the
    [i]th equation with respect to the [j]th variable, evaluated at the
    values of [t] and [y] obtained from [arg]. Only nonzero elements need
    be loaded into [jac].

    Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
    Any other exception is treated as an unrecoverable error.

    {warning Neither the elements of [arg] nor the matrix [jac] should
             be accessed after the function has returned.}

    @nocvodes <node5#ss:sjacFnBS> CVSlsSparseJacFnBS *)
type sparse_jac_fn_with_sens =
  (Sundials.RealArray.t Cvode.triple, Sundials.RealArray.t)
      Cvodes.Adjoint.jacobian_arg
  -> Sundials.RealArray.t array -> Sls.SparseMatrix.t -> unit

(** Callback functions that compute sparse approximations to a Jacobian
    matrix.

    @nocvodes <node5#ss:sjacFnB> CVSlsSparseJacFnB
    @nocvodes <node5#ss:sjacFnBS> CVSlsSparseJacFnBS *)
type sparse_jac_fn =
    NoSens of sparse_jac_fn_no_sens
    (** Does not depend on forward sensitivities. *)
  | WithSens of sparse_jac_fn_with_sens
    (** Depends on forward sensitivities. *)

(** A direct linear solver on sparse matrices. In the call,
    [superlumt jfn nnz], [jfn] specifies a callback function that computes an
    approximation to the Jacobian matrix and [nnz] is the maximum number
    of nonzero entries in that matrix.

    @nocvodes <node5#sss:lin_solv_init> CVSuperLUMTB
    @nocvodes <node5#sss:optin_sls> CVSlsSetSparseJacFnB
    @nocvodes <node5#sss:optin_sls> CVSlsSetSparseJacFnBS
    @nocvodes <node5#ss:sjacFnB> CVSlsSparseJacFnB
    @nocvodes <node5#ss:sjacFnBS> CVSlsSparseJacFnBS *)
val superlumt : sparse_jac_fn -> nnz:int -> nthreads:int
                  -> Cvodes.Adjoint.serial_linear_solver

(** The ordering algorithm used for reducing fill. *)
type ordering = Cvode_superlumt.ordering =
     Natural       (** Natural ordering. *)
   | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
   | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
   | ColAmd        (** Column approximate minimum degree permutation. *)

(** Sets the ordering algorithm used to minimize fill-in.

    @nocvodes <node5#ss:sls_optin> CVSuperLUMTSetOrdering
    @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
val set_ordering : Cvodes.Adjoint.serial_bsession -> ordering -> unit

(** Returns the number of calls made by a sparse linear solver to the
    Jacobian approximation function.

    @nocvodes <node5#sss:optout_sls> CVSlsGetNumJacEvals
    @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
val get_num_jac_evals : Cvodes.Adjoint.serial_bsession -> int

