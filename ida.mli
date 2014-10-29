(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(***********************************************************************)
(* Parts of the comment text are taken directly from:                  *)
(*                                                                     *)
(*               User Documentation for IDA v2.7.0                     *)
(*         Alan C. Hindmarsh, Radu Serban, and Aaron Collier           *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Variable-step solution of DAE initial value problems with
    zero-crossing detection.

    This module solves numerically problems of the form
    {% $F(t, y, \dot{y})=0$%}, {% $y(t_0) = y_0$%},
    {% $\dot{y}(t_0)=\dot{y}_0$%}.

    This documented interface is structured as follows.
    {ol
      {- {{:#linear}Linear solvers}}
      {- {{:#tols}Tolerances}}
      {- {{:#init}Initialization}}
      {- {{:#solver}Solving}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}
      {- {{:#exceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS) *)

open Sundials

(** A session with the IDA solver.

    An example session with Ida ({openfile ida_skel.ml}): {[
#include "examples/ocaml/skeletons/ida_skel.ml"
    ]}

    @ida <node5#ss:skeleton_sim> Skeleton of main program
 *)
type ('a, 'k) session = ('a, 'k) Ida_impl.session

(** Alias for sessions based on serial nvectors. *)
type serial_session = (RealArray.t, Nvector_serial.kind) session

(** {2:linear Linear Solvers} *)

(** Linear solvers used by Ida.

    @ida <node5#sss:lin_solv_init> Linear Solver Specification Functions *)
type ('data, 'kind) linear_solver = ('data, 'kind) Ida_impl.linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type serial_linear_solver =
      (Nvector_serial.data, Nvector_serial.kind) linear_solver

(** Workspaces with two temporary vectors. *)
type 'a double = 'a * 'a

(** Workspaces with three temporary vectors. *)
type 'a triple = 'a * 'a * 'a

(** Arguments common to Jacobian callback functions.    
   
    @ida <node5#ss:djacFn> IDADlsDenseJacFn
    @ida <node5#ss:bjacFn> IDADlsBandJacFn
    @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn
    @ida <node5#ss:psolveFn> IDASpilsPrecSolveFn
    @ida <node5#ss:precondFn> IDASpilsPrecSetupFn *)
type ('t, 'a) jacobian_arg =
  {
    jac_t    : float;        (** The independent variable. *)
    jac_y    : 'a;           (** The dependent variable vector. *)
    jac_y'   : 'a;           (** The derivative vector (i.e. dy/dt). *)
    jac_res  : 'a;           (** The current value of the residual vector. *)
    jac_coef : float;        (** The coefficient $c_j$ in
                                 {% $J = \frac{\partial F}{\partial y} + c_j
                                     \frac{\partial F}{\partial\dot{y}}$%}. *)
    jac_tmp  : 't            (** Workspace data. *)
  }

(** The range of nonzero entries in a band matrix.  *)
type bandrange = Ida_impl.bandrange =
  { mupper : int; (** The upper half-bandwidth.  *)
    mlower : int; (** The lower half-bandwidth.  *) }

(** Direct Linear Solvers operating on dense and banded matrices.

    @ida <node5#sss:optin_dls> Direct linear solvers optional input functions
    @ida <node5#sss:optout_dls> Direct linear solvers optional output functions *)
module Dls :
  sig

    (** Callback functions that compute dense approximations to a Jacobian
        matrix. In the call [dense_jac_fn arg jac], [arg] is a {!jacobian_arg}
        with three work vectors and the computed Jacobian must be stored
        in [jac].

        The callback should load the [(i,j)]th entry of [jac] with
        {% $\frac{\partial F_i}{\partial y_j}
            + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%}, i.e., the partial
        derivative of the [i]th equation with respect to the [j]th variable,
        evaluated at the values of [t], [y], and [y'] obtained from [arg].
        Only nonzero elements need be loaded into [jac].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jac] should
                 be accessed after the function has returned.}

        @ida <node5#ss:djacFn> IDADlsDenseJacFn *)
    type dense_jac_fn = (RealArray.t triple, RealArray.t) jacobian_arg
                         -> Dls.DenseMatrix.t -> unit

    (** A direct linear solver on dense matrices. The optional argument
        specifies a callback function for computing an approximation to the
        Jacobian matrix. If this argument is omitted, then a default
        implementation based on difference quotients is used.

        @ida <node5#sss:lin_solv_init> IDADense
        @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
        @ida <node5#ss:djacFn> IDADlsDenseJacFn *)
    val dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

    (** A direct linear solver on dense matrices using LAPACK. See {!Dense}.
        Only available if {!lapack_enabled}.

        @ida <node5#sss:lin_solv_init> IDALapackDense
        @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
        @ida <node5#ss:djacFn> IDADlsDenseJacFn *)
    val lapack_dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

    (** Callback functions that compute banded approximations to
        a Jacobian matrix. In the call [band_jac_fn {mupper; mlower} arg jac],
        - [mupper] is the upper half-bandwidth of the Jacobian,
        - [mlower] is the lower half-bandwidth of the Jacobian,
        - [arg] is a {!jacobian_arg} with three work vectors, and,
        - [jac] is storage for the computed Jacobian.

        The callback should load the [(i,j)]th entry of [jac] with
        {% $\frac{\partial F_i}{\partial y_j}
            + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%}, i.e., the partial
        derivative of the [i]th equation with respect to the [j]th variable,
        evaluated at the values of [t] and [y] obtained from [arg]. Only
        nonzero elements need be loaded into [jac].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jac] should
                 be accessed after the function has returned.}

        @ida <node5#ss:bjacFn> IDADlsBandJacFn *)
    type band_jac_fn = bandrange
                        -> (RealArray.t triple, RealArray.t) jacobian_arg
                        -> Dls.BandMatrix.t -> unit

    (** A direct linear solver on banded matrices. The optional argument
        specifies a callback function for computing an approximation to the
        Jacobian matrix. If this argument is omitted, then a default
        implementation based on difference quotients is used. The other
        argument gives the width of the bandrange.

        @ida <node5#sss:lin_solv_init> IDABand
        @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
        @ida <node5#ss:bjacFn> IDADlsBandJacFn *)
    val band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

    (** A direct linear solver on banded matrices using LAPACK. See {!band}.
        Only available if {!lapack_enabled}.

        @ida <node5#sss:lin_solv_init> IDALapackBand
        @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
        @ida <node5#ss:bjacFn> IDADlsBandJacFn *)
    val lapack_band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by a direct
        linear solver.

        @ida <node5#sss:optout_dls> IDADlsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : serial_session -> int * int


    (** Returns the number of calls made by a direct linear solver to the
        Jacobian approximation function.

        @ida <node5#sss:optout_dls> IDADlsGetNumJacEvals *)
    val get_num_jac_evals : serial_session -> int

    (** Returns the number of calls to the residual callback due to
        the finite difference Jacobian approximation.

        @ida <node5#sss:optout_dls> IDADlsGetNumResEvals *)
    val get_num_res_evals : serial_session -> int

    (** {3:lowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change a Jacobian function. These low-level functions are provided for
        experts who want to avoid resetting internal counters and other
        associated side-effects. *)

    (** Change the dense Jacobian function.

        @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn *)
    val set_dense_jac_fn : serial_session -> dense_jac_fn -> unit

    (** Remove a dense Jacobian function and use the default
        implementation.

        @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn *)
    val clear_dense_jac_fn : serial_session -> unit

    (** Change the band Jacobian function.

        @ida <node5#sss:optin_dls> IDADlsSetBandJacFn *)
    val set_band_jac_fn : serial_session -> band_jac_fn -> unit

    (** Remove a banded Jacobian function and use the default
        implementation.

        @ida <node5#sss:optin_dls> IDADlsSetBandJacFn *)
    val clear_band_jac_fn : serial_session -> unit
  end

(** Scaled Preconditioned Iterative Linear Solvers.

    @ida <node5#sss:optin_spils> Iterative linear solvers optional input functions.
    @ida <node5#sss:optout_spils> Iterative linear solvers optional output functions. *)
module Spils :
  sig
    (** {3:precond Preconditioners} *)

    (** Callback functions that solve a linear system involving a
        preconditioner matrix.
        In the call [prec_solve_fn jac r z delta],
        [jac] is a {!jacobian_arg} with one work vector,
        [r] is the right-hand side vector,
        [z] is computed to solve {% $Pz = r$%},
        and [delta] is the input tolerance.
        $P$ is a preconditioner matrix, which approximates, however crudely,
        the Jacobian matrix {% $\frac{\partial F}{\partial y}
                  + \mathtt{arg.jac\_coef}\frac{\partial F}{\partial\dot{y}}$%}.
        If the solution is found via an iterative method, it must satisfy
        {% $\sqrt{\sum_i (\mathit{Res}_i \cdot \mathit{ewt}_i)^2}
              < \mathtt{delta}$%},
        where {% $\mathit{Res} = r - Pz$%} and {% $\mathit{ewt}$%} comes from
        {!get_err_weights}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac], [r], and [z] should not
                 be accessed after the function has returned.}

        @ida <node5#ss:psolveFn> IDASpilsPrecSolveFn *)
    type 'a prec_solve_fn =
      ('a, 'a) jacobian_arg
      -> 'a
      -> 'a
      -> float
      -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related data
        need by {!prec_solve_fn}. The sole argument is a {!jacobian_arg} with
        three work vectors.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The elements of the argument should not be accessed after the
                 function has returned.}

        @ida <node5#ss:precondFn> IDASpilsPrecSetupFn *)
    type 'a prec_setup_fn = ('a triple, 'a) jacobian_arg -> unit

    (** Callback functions that compute the Jacobian times a vector. In the
        call [jac_times_vec_fn arg v jv], [arg] is a {!jacobian_arg} with two
        work vectors, [v] is the vector multiplying the Jacobian, and [jv] is
        the vector in which to store the
        result—{% $\mathtt{jv} = J\mathtt{v}$%}.
      
        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor [v] or [jv] should be
                 accessed after the function has returned.}

        @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn *)
    type 'a jac_times_vec_fn =
      ('a double, 'a) jacobian_arg
      -> 'a
      -> 'a
      -> unit

    (** Specifies a preconditioner and its callback functions.
        The following functions and those in {!Ida_bbd} construct
        preconditioners.

        The {!prec_solve_fn} is usually mandatory. The {!prec_setup_fn} can be
        omitted if not needed. If the {!jac_times_vec_fn} is omitted, a
        default implementation based on difference quotients is used.

        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
        @ida <node5#ss:psolveFn> IDASpilsPrecSolveFn
        @ida <node5#ss:precondFn> IDASpilsPrecSetupFn
        @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn *)
    type ('a, 'k) preconditioner = ('a, 'k) Ida_impl.SpilsTypes.preconditioner

    (** No preconditioning.  *)
    val prec_none : ('a, 'k) preconditioner

    (** Left preconditioning. {% $Pz = r$%}, where $P$ approximates, perhaps
        crudely, {% $J = \frac{\partial F}{\partial y}
                            + c_j\frac{\partial F}{\partial\dot{y}}$%}. *)
    val prec_left :
      ?setup:'a prec_setup_fn
      -> ?jac_times_vec:'a jac_times_vec_fn
      -> 'a prec_solve_fn
      -> ('a, 'k) preconditioner

    (** {3:lsolvers Solvers} *)

    (** Krylov iterative solver using the scaled preconditioned generalized
        minimum residual (GMRES) method.
        In the call [spgmr ~maxl:maxl ~max_restarts:maxr prec],
        [maxl] is the maximum dimension of the Krylov subspace (defaults to 5),
        [maxr] is the maximum number of restarts (defaults to 5),
        and [prec] is a {!preconditioner}.

        @ida <node5#sss:lin_solv_init> IDASpgmr
        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#sss:optin_spils> IDASpilsSetMaxl
        @ida <node5#sss:optin_spils> IDASpilsSetMaxRestarts *)
    val spgmr : ?maxl:int -> ?max_restarts:int
      -> ('a, 'k) preconditioner -> ('a, 'k) linear_solver

    (** Krylov iterative solver using the scaled preconditioned biconjugate
        stabilized (Bi-CGStab) method.
        In the call [spbcg ~maxl:maxl prec], [maxl] is the maximum dimension of
        the Krylov subspace (defaults to 5), and [prec] is a {!preconditioner}.

        @ida <node5#sss:lin_solv_init> IDASpbcg
        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#sss:optin_spils> IDASpilsSetMaxl *)
    val spbcg : ?maxl:int -> ('a, 'k) preconditioner -> ('a, 'k) linear_solver

    (** Krylov iterative with the scaled preconditioned transpose-free
        quasi-minimal residual (SPTFQMR) method.
        In the call [sptfqmr ~maxl:maxl prec], [maxl] is the maximum dimension
        of the Krylov subspace (defaults to 5), and [prec] is a
        {!preconditioner}.

        @ida <node5#sss:lin_solv_init> IDASptfqmr
        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#sss:optin_spils> IDASpilsSetMaxl *)
    val sptfqmr : ?maxl:int -> ('a, 'k) preconditioner -> ('a, 'k) linear_solver

    (** {3:set Solver parameters} *)

    (** The type of Gram-Schmidt orthogonalization.

        @ida <node9#ss:spgmr> ModifiedGS/ClassicalGS *)
    type gramschmidt_type = Spils.gramschmidt_type =
      | ModifiedGS   (** Modified Gram-Schmidt orthogonalization
                         {cconst MODIFIED_GS} *)
      | ClassicalGS  (** Classical Gram Schmidt orthogonalization
                         {cconst CLASSICAL_GS} *)

    (** Sets the Gram-Schmidt orthogonalization to be used with the
        Spgmr {!linear_solver}.

        @ida <node5#sss:optin_spils> IDASpilsSetGSType *)
    val set_gs_type : ('a, 'k) session -> Spils.gramschmidt_type -> unit

    (** Sets the factor by which the Krylov linear solver's convergence test
        constant is reduced from the Newton iteration test constant.
        This factor must be >= 0; passing 0 specifies the default (0.05).

        @ida <node5#sss:optin_spils> IDASpilsSetEpsLin *)
    val set_eps_lin : ('a, 'k) session -> float -> unit

    (** Resets the maximum Krylov subspace dimension for the Bi-CGStab and
        TFQMR methods. A value <= 0 specifies the default (5.0).

        @ida <node5#sss:optin_spils> IDASpilsSetMaxl *)
    val set_maxl : ('a, 'k) session -> int -> unit

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the spils
        linear solver.

        @ida <node5#sss:optout_spils> IDASpilsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space       : ('a, 'k) session -> int * int

    (** Returns the cumulative number of linear iterations.

        @ida <node5#sss:optout_spils> IDASpilsGetNumLinIters *)
    val get_num_lin_iters    : ('a, 'k) session -> int

    (** Returns the cumulative number of linear convergence failures.

        @ida <node5#sss:optout_spils> IDASpilsGetNumConvFails *)
    val get_num_conv_fails   : ('a, 'k) session -> int

    (** Returns the number of calls to the setup function.

        @ida <node5#sss:optout_spils> IDASpilsGetNumPrecEvals *)
    val get_num_prec_evals   : ('a, 'k) session -> int

    (** Returns the cumulative number of calls to the preconditioner solve
        function.

        @ida <node5#sss:optout_spils> IDASpilsGetNumPrecSolves *)
    val get_num_prec_solves  : ('a, 'k) session -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        function.

        @ida <node5#sss:optout_spils> IDASpilsGetNumJtimesEvals *)
    val get_num_jtimes_evals : ('a, 'k) session -> int

    (** Returns the number of calls to the residual callback for
        finite difference Jacobian-vector product approximation. This counter is
        only updated if the default difference quotient function is used.

        @ida <node5#sss:optout_spils> IDASpilsGetNumResEvals *)
    val get_num_res_evals    : ('a, 'k) session -> int

    (** {3:lowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change preconditioner functions. These low-level functions are provided
        for experts who want to avoid resetting internal counters and other
        associated side-effects. *)

    (** Change the preconditioner functions.

        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#ss:psolveFn> IDASpilsPrecSolveFn
        @ida <node5#ss:precondFn> IDASpilsPrecSetupFn *)
     val set_preconditioner :
       ('a,'k) session
       -> ?setup:'a prec_setup_fn
       -> 'a prec_solve_fn
       -> unit

    (** Change the Jacobian-times-vector function.

        @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
        @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn *)
    val set_jac_times_vec_fn :
      ('a,'k) session
      -> 'a jac_times_vec_fn
      -> unit

    (** Remove a Jacobian-times-vector function and use the default
        implementation.

        @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
        @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn *)
    val clear_jac_times_vec_fn : ('a, 'k) session -> unit
  end

(** Alternate Linear Solvers.

    @ida <node8#s:new_linsolv> Providing Alternate Linear Solver Modules *)
module Alternate :
  sig

    (** Functions that initialize linear solver data, like counters and
        statistics.

        Raising any exception in this function (including
        {!Sundials.RecoverableFailure}) is treated as an unrecoverable error.

        @ida <node8#SECTION00810000000000000000> linit *)
    type ('data, 'kind) linit = ('data, 'kind) session -> unit

    (** Functions that prepare the linear solver for subsequent calls to
        {!callbacks.lsolve}. The call [lsetup s y y' res tmp] has as
        arguments

        - [s], the solver session,
        - [y],  the predicted $y$ vector for the current internal step,
        - [y'], the predicted {% $\dot{y}$%} vector for the current internal
                step,
        - [res], the value of the residual function at [y] and [y'], i.e.
                 {% $F(t_n, y_{\text{pred}}, \dot{y}_{\text{pred}})$%}, and,
        - [tmp], temporary variables for use by the routine.

        This function may raise a {!Sundials.RecoverableFailure} exception to
        indicate that a recoverable error has occurred. Any other exception is
        treated as an unrecoverable error.

        {warning The vectors [y], [y'], [res], and those in [tmp] should not
                 be accessed after the function returns.}

        @ida <node8#SECTION00820000000000000000> lsetup *)
    type ('data, 'kind) lsetup =
      ('data, 'kind) session
      -> 'data
      -> 'data
      -> 'data
      -> 'data triple
      -> unit

    (** Functions that solve the linear equation $Mx = b$.
        $M$ is a preconditioning matrix chosen by the user, and $b$ is the
        right-hand side vector calculated within the function.
        $M$ should approximate {% $J = \frac{\partial F}{\partial y}
            + c_j\frac{\partial F}{\partial \dot{y}}$%}, and $c_j$ is
        available through {!get_cj}.
        The call [lsolve s b weight ycur y'cur rescur] has as arguments:

        - [s], the solver session,
        - [b], for returning the calculated solution,
        - [weight], the error weights,
        - [ycur], the solver's current approximation to $y(t_n)$,
        - [y'cur], the solver's current approximation to {% $\dot{y}(t_n)$%},
                   and,
        - [rescur], a vector containing the current residual value.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        @ida <node8#SECTION00830000000000000000> lsolve
        @ida <node3#e:DAE_Jacobian> IVP solution (Eq. 2.5) *)
    type ('data, 'kind) lsolve =
      ('data, 'kind) session
      -> 'data
      -> 'data
      -> 'data
      -> 'data
      -> 'data
      -> unit

    (** The callbacks needed to implement an alternate linear solver. *)
    type ('data, 'kind) callbacks =
      {
        linit  : ('data, 'kind) linit option;
        lsetup : ('data, 'kind) lsetup option;
        lsolve : ('data, 'kind) lsolve
      }

    (** Creates a linear solver from a function returning a set of
        callbacks. The creation function is passed a session and a vector.
        The latter indicates the problem size and can, for example, be
        cloned. *)
    val make :
          (('data, 'kind) session
            -> ('data, 'kind) Nvector.t
            -> ('data, 'kind) Nvector.t
            -> ('data, 'kind) callbacks)
          -> ('data, 'kind) linear_solver

    (** {3:internals Solver internals} *)

    (** Returns the current [cj] value. *)
    val get_cj : ('data, 'kind) session -> float

    (** Returns the current [cjratio] value. *)
    val get_cjratio : ('data, 'kind) session -> float
  end

(** {2:tols Tolerances} *)

(** Functions that set the multiplicative error weights for use in the weighted
    RMS norm. The call [efun y ewt] takes the dependent variable vector [y] and
    fills the error-weight vector [ewt] with positive values or raises
    {!NonPositiveEwt}. Other exceptions are eventually propagated, but
    should be avoided ([efun] is not allowed to abort the solver). *)
type 'data error_fun = 'data -> 'data -> unit

type ('data, 'kind) tolerance =
  | SStolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVtolerances of float * ('data, 'kind) Nvector.t
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFtolerances of 'data error_fun
    (** Set the multiplicative error weights for the weighted RMS norm. *)

(** A default relative tolerance of 1.0e-4 and absolute tolerance of 1.0e-8. *)
val default_tolerances : ('data, 'kind) tolerance

(** {2:init Initialization} *)

(** The DAE's residual function.  Called by the solver like [f t y y' r],
    where:
    - [t] is the current value of the independent variable,
          i.e., the simulation time.
    - [y] is a vector of dependent-variable values, i.e. y(t).
    - [y'] is the derivative of [y] with respect to [t], i.e. dy/dt.
    - [r] is the output vector to fill in with the value of the system
          residual for the given values of t, y, and y'.
    The residual function should return normally if successful, raise
    {!Sundials.RecoverableFailure} if a recoverable error occurred (e.g. [y] has
    an illegal value), or raise some other exception if a nonrecoverable error
    occurred.  If a recoverable error occurred, the integrator will attempt to
    correct and retry.  If a nonrecoverable error occurred, the integrator will
    halt and propagate the exception to the caller.

    {b NB:} [y], [y'], and [r] must no longer be accessed after [f] has
            returned a result, i.e. if their values are needed outside of
            the function call, then they must be copied to separate physical
            structures.

    @ida <node5#ss:resFn> IDAResFn *)
type 'a resfn = float -> 'a -> 'a -> 'a -> unit

(** A function called by the solver to calculate the values of root
    functions (zero-crossing expressions) which are used to detect
    significant events.  It is passed four arguments [t], [y], [y'],
    and [gout]:
    - [t] is the current value of the independent variable,
          i.e., the simulation time.
    - [y] is a vector of dependent-variable values, i.e. y(t).
    - [y'] is the derivative of [y] with respect to [t], i.e. dy/dt.
    - [gout] is a vector for storing the values of g(t, y, y').
    Note that [t], [y], [y'] are the same as for {!resfn}.  If the
    labeled argument ~roots is omitted, then no root finding is
    performed.  If the root function raises an exception, the
    integrator will halt immediately and propagate the exception to
    the caller.

    {b NB:} [y] and [gout] must no longer be accessed after [g] has returned
            a result, i.e. if their values are needed outside of the function
            call, then they must be copied to separate physical structures.

    See also {!init} and {!reinit}.

    @ida <node5#ss:rootFn> IDARootFn *)
type 'a rootsfn = float -> 'a -> 'a -> RealArray.t -> unit

(** [init linsolv tol f ~roots:(nroots, g) ~t0:t0 y0 y'0] initializes
    the IDA solver to solve the DAE f t y y' = 0 and returns a
    {!session}.

    - [linsolv] is the linear solver to attach to this session,
    - [tol]     specifies the integration tolerances,
    - [f]       is the residual function (see below),
    - [nroots]  specifies the number of root functions (zero-crossings),
    - [g]       calculates the values of the root functions,
    - [t0]      is the initial value of the independent variable t, which
                defaults to 0,
    - [y0]      is a vector of initial values for the dependent-variable vector
                [y].  This vector's size determines the number of equations
                in the session, see {!RealArray.t}, and,
    - [y'0]     is a vector of initial values for [y'], i.e. the derivative
                of [y] with respect to t.  This vector's size must match the
                size of [y0].

    The labeled arguments [roots] and [t0] are both optional and default to
    {!no_roots} (i.e. no root finding is done) and [0.0], respectively.

    This function calls IDACreate, IDAInit, IDARootInit, an
    appropriate linear solver function, and one of IDASStolerances,
    IDASVtolerances, or IDAWFtolerances. It does everything necessary
    to initialize an IDA session; the {!solve_normal} or
    {!solve_one_step} functions can be called directly afterward.

    @ida <node5#sss:idainit>       IDACreate/IDAInit
    @ida <node5#ss:idarootinit>    IDARootInit
    @ida <node5#sss:lin_solv_init> Linear solvers
    @ida <node5#sss:idatolerances> IDASStolerances
    @ida <node5#sss:idatolerances> IDASVtolerances
    @ida <node5#sss:idatolerances> IDAWFtolerances
    @ida <node5#ss:ewtsetFn>       IDAEwtFn *)
val init :
    ('a, 'kind) linear_solver
    -> ('a, 'kind) tolerance
    -> 'a resfn
    -> ?varid:('a, 'kind) Nvector.t
    -> ?roots:(int * 'a rootsfn)
    -> float
    -> ('a, 'kind) Nvector.t
    -> ('a, 'kind) Nvector.t
    -> ('a, 'kind) session

(** This is a convenience value for signalling that there are no
    roots (zero-crossings) to monitor. *)
val no_roots : (int * 'a rootsfn)

(** Return the number of root functions. *)
val nroots : ('a, 'k) session -> int

(** {3 Initial Value Calculation} *)

(** Symbolic names for variable type classifications needed to
    calculate initial values (see {!calc_ic_ya_yd'}) or to suppress
    local error tests on some variables (see {!suppress_alg}).

    Those functions require you to pass in an nvector populated with
    magic constants specifying each variable as algebraic or
    differential.  This module gives symbolic names to those
    constants, for your convenience.

    Note: variable type classification is called "id" in the C
    interface, e.g. [IDASetId].

    @ida <node5#sss:idasetid> IDASetId
    @ida <node5#sss:optin_main> IDASetSuppressAlg *)
module VarId :
  sig
    (** A symbolic name for the magic floating-point constant [0.0]. *)
    val algebraic : float
    (** A symbolic name for the magic floating-point constant [1.0]. *)
    val differential : float

    (** An ADT representation of the magic constants specifying
        variable types, useful for pattern-matching.  *)
    type t =
    | Algebraic    (** Algebraic variable; residual function must not depend
                       on this component's derivative.  Corresponds to
                       numerical value 0.0.  *)
    | Differential (** Differential variable; residual function can depend on
                       this component's derivative.  Corresponds to numerical
                       value 1.0.  *)

    (** Encode an [Algebraic] / [Differential] into the corresponding
        magic floating-point constant.  *)
    val to_float : t -> float

    (** Decode a magic float-point constant into an [Algebraic] /
        [Differential] specification.  Raises [Invalid_argument] if
        the given floating point value is not a legal variable type
        specification.  *)
    val of_float : float -> t

    (** Maps [algebraic -> "Algebraic"] and [differential ->
        "Differential"].  Raises [Invalid_argument] if the given
        floating point value is not a legal variable type
        specification.  *)
    val string_of_float : float -> string

    (** Returns ["Algebraic"] or ["Differential"].  *)
    val string_of_var_type : t -> string
  end

(** Specify which variables are algebraic and which variables are
    differential, needed for {!set_suppress_alg}.  This function must
    not be called if you already called {!calc_ic_ya_yd'}.

    The SUNDIALS manual is not clear about whether it's safe to change the
    variable types after you've already set it.

    [set_var_types] corresponds to [IDASetId] in the C interface, and an alias
    {!set_id} is also available in this binding.  We prefer the more
    descriptive name {!set_var_types}, however.

    @ida <node5#sss:optin_main> IDASetId *)
val set_id : ('a, 'k) session -> ('a,'k) Nvector.t -> unit

(** Indicate whether or not to ignore algebraic variables in the local
    error test.  This is set to [false] by default.  Before you can
    set it to [true], you must specify which variables are algebraic
    through {!calc_ic_ya_yd'} or {!set_var_types}, but not both.

    Exactly one of these functions should be called, exactly once,
    before the first call to {!solve_normal}, {!solve_one_step}, or
    {!calc_ic_ya_yd'}.  Forgetting to do so will cause an
    {!Ida.IllInput} exception.

    Note: {!set_var_types} is the preferred alias to {!set_id}, which
    corresponds to [IDASetId] in the C interface.

    In general, suppressing local error tests for algebraic variables
    is {i discouraged} when solving DAE systems of index 1, whereas it
    is generally {i encouraged} for systems of index 2 or more.  See
    pp. 146-147 of the following reference for more on this issue:

    K. E. Brenan, S. L. Campbell, and L. R. Petzold.  Numerical
    Solution of Initial-Value Problems in Differential-Algebraic
    Equations.  SIAM, Philadelphia, Pa, 1996.

    @ida <node5#sss:optin_main> IDASetId
    @ida <node5#sss:optin_main> IDASetSuppressAlg *)
val set_suppress_alg : ('a, 'k) session
                       -> ?varid:('a, 'k) Nvector.t -> bool -> unit

(** [calc_ic_y ida ~y:yvar tout1] corrects the initial values y0 at
    time t0, using the initial values of the derivatives y'0.  That
    is, if the {i t0,y0,y'0} that were given to {!init} or {!reinit}
    does not satisfy {i f(t0,y0,y'0) = 0}, where {i f} is the residual
    function, then [calc_ic_y] will modify {i y'0} so that this
    equation holds.  If {i f(t0,y0,y'0) = 0} is already true, a call
    to [calc_ic_y] is unnecessary.  [calc_ic_y] must not be called
    after any calls to {!solve_normal} or {!solve_one_step} without a
    {!reinit} in between.

    The optional parameter [~y], if given, will receive the corrected
    {i y} vector.  [tout1] is the first value of {i t} at which a
    solution will be requested (using {!solve_normal} or
    {!solve_one_step}). This value is needed here only to determine
    the direction of integration and rough scale in the independent
    variable {i t}.

    [calc_ic_y] differs from {!calc_ic_ya_yd'} in that,
    {!calc_ic_ya_yd'} computes parts of y and y' using parts of y as
    input, whereas [calc_ic_y] computes all of y using all of y' as
    input.  Here, y means the vector formed by collecting scalar
    variables that appear in the mathematical description of the DAE
    system, and y' is its derivative.  This is not to be confused with
    the labeled argument whose name is [~y]: y and y' are mathematical
    objects whereas [~y] is a programming construct.

    IDA's initial value correction works for certain index-one
    problems including a class of systems of semi-implicit form, and
    uses Newton iteration combined with a linesearch algorithm.  See
    Section 2.1 of the IDA User Guide and the following reference for
    more information:

    P. N. Brown, A. C. Hindmarsh, and L. R. Petzold. Consistent Initial Condition Calculation for Differential-Algebraic Systems. SIAM J. Sci. Comput., 19:1495-1512, 1998.

    @ida <node5#ss:idacalcic> IDACalcIC
    @ida <node5#sss:optout_iccalc> IDAGetConsistentIC *)
val calc_ic_y : ('a, 'k) session -> ?y:('a, 'k) Nvector.t -> float -> unit

(** [calc_ic_ya_yd' ida ~y:yvar ~y':y'var vartypes tout1] corrects the
    initial values y0 and y0' at time t0.  That is, if the {i
    t0,y0,y'0} that were given to {!init} or {!reinit} does not
    satisfy {i f(t0,y0,y'0) = 0}, where {i f} is the residual
    function, then [calc_ic_ya_yd'] will modify parts of {i y0} and {i
    y'0} so that this equation holds.  If {i f(t0,y0,y'0) = 0} is
    already true, a call to [calc_ic_ya_yd'] is unnecessary.
    [calc_ic_ya_yd'] must not be called after any calls to
    {!solve_normal} or {!solve_one_step} without a {!reinit} in
    between.

    The optional parameters [~y] and [~y'], if given, will receive the
    corrected vectors.  [tout1] is the first value of t at which a
    solution will be requested (using {!solve_normal} or
    {!solve_one_step}), and is needed here only to determine the
    direction of integration and rough scale in the independent
    variable t.

    {!calc_ic_y} differs from [calc_ic_ya_yd'] in that,
    [calc_ic_ya_yd'] computes parts of y and y' using parts of y as
    input, whereas {!calc_ic_y} computes all of y using all of y' as
    input.  Here, y means the vector formed by collecting scalar
    variables that appear in the mathematical description of the DAE
    system, and y' means its derivative.  These are not to be confused
    with the labeled arguments whose names are [~y] and [~y']: y and
    y' are mathematical objects whereas [~y] and [~y'] are programming
    constructs.

    The non-optional nvector argument, named [vartypes] at the
    beginning, specifies some components of y as algebraic (i.e. their
    derivatives do not appear in the DAE) and others as differential
    (i.e. their derivatives appear in the DAE).  [calc_ic_ya_yd']
    modifies the algebraic components of y and differential components
    of y', using the differential components of y as input.  So if we
    let Ia be the set of indices at which [vartypes] is [Algebraic]
    and Id be the set of indices at which [vartypes] is
    [Differential], then y and y' are each partitioned into two
    sub-vectors (we use OCaml's array-indexing notation to denote
    subscripting):

      - y  splits into A  = \{ y.(i)  | i in Ia \} and D  = \{ y.(i)  | i in Id \}
      - y' splits into A' = \{ y'.(i) | i in Ia \} and D' = \{ y'.(i) | i in Id \}

    The residual function must be such that it ignores all values in
    A'.  [calc_ic_ya_yd'] computes (i.e. modifies) A and D' while
    treating D as read-only and ignoring A'.

      input:   D
      output:  A, D'
      ignored: A'

    Note: [vartypes] is called "id" in the C interface, e.g. [IDASetId].

    [calc_ic_ya_yd'] sets the variable types that {!set_suppress_alg}
    uses, so you do not need to set it again with {!set_var_types} (or
    its alias {!set_id}) before calling {!set_suppress_alg}.

    Note: the nvector interface gives no way of checking that the [~y]
    and [~y'] vectors have the right sizes.  Passing incorrectly sized
    vectors leads to memory corruption, so beware!  It's a good idea
    to always reuse the nvectors you gave to {!init} or {!reinit}.

    @ida <node5#ss:idacalcic> IDACalcIC
    @ida <node5#sss:optin_main> IDASetId
    @ida <node5#sss:optout_iccalc> IDAGetConsistentIC *)
val calc_ic_ya_yd' :
  ('a, 'k) session
  -> ?y:('a, 'k) Nvector.t
  -> ?y':('a, 'k) Nvector.t
  -> ?varid:('a, 'k) Nvector.t
  -> float
  -> unit

(** [get_num_backtrack_ops ida] gets the number of backtrack operations done in
    the linesearch algorithm in {!calc_ic_ya_yd'} or {!calc_ic_y}.

    @ida <node5#sss:optout_iccalc> IDAGetNumBcktrackOps *)
val get_num_backtrack_ops : ('a, 'k) session -> int

(** {2:solver Solving} *)

(** Possible values returned when a IDA solver step function succeeds.
    Failures are indicated by exceptions.

    @ida <node5#sss:ida> IDASolve *)
type solver_result =
  | Success             (** IDA_SUCCESS *)
  | RootsFound          (** IDA_ROOT_RETURN *)
  | StopTimeReached     (** IDA_TSTOP_RETURN *)

(** [(tret, r) = solve_normal s tout yout y'out] integrates the DAE
    over an interval in t.

   The arguments are:
   - [s] a session with the solver.
   - [tout] the next time at which a computed solution is desired.
   - [yout] a vector to store the computed solution. The same vector that was
            passed to {!init} can be (but does not have to be) reused for this
            argument.
   - [y'out] a vector to store the computed solution's derivative.
             The same vector that was passed to {!init} can be (but
             does not have to be) reused for this argument.

   Two values are returned:
    - [tret] the time reached by the solver, which will be equal to [tout] if
      no errors occur.
    - [r] indicates whether roots were found, or whether an optional stop time,
          set by {!set_stop_time}, was reached; see {!Sundials.solver_result}.

   This routine will throw one of the solver {!exceptions} if an error
   occurs.

   Note: the nvector interface gives no way of checking that the
   [yout] and [y'out] vectors have the right sizes.  Passing
   incorrectly sized vectors leads to memory corruption, so beware!
   It's a good idea to always reuse the nvectors you gave to {!init}
   or {!reinit}.

   @ida <node5#sss:idasolve> IDASolve (IDA_NORMAL) *)
val solve_normal : ('a, 'k) session -> float
                   -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t
                   -> float * solver_result

(** This function is identical to {!solve_normal}, except that it
    returns after one internal solver step.

    @ida <node5#sss:idasolve> IDASolve (IDA_ONE_STEP) *)
val solve_one_step : ('a, 'k) session -> float
                     -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t
                     -> float * solver_result

(**
  [get_dky s t k dky] computes the [k]th derivative of the function y at time
  [t], i.e. d(k)y/dt(k)(t). The function requires that tn - hu <= [t] <=
  tn, where tn denotes the current internal time reached, and hu is the last
  internal step size successfully used by the solver.
  The user may request [k] = 0, 1,..., qu, where qu is the current order.

  This function may only be called after a successful return from either
  {!solve_normal} or {!solve_one_step}.

  Values for the limits may be obtained:
    - tn = {!get_current_time}
    - qu = {!get_last_order}
    - hu = {!get_last_step}

  @ida <node5#sss:optin_root> IDAGetDky
 *)
val get_dky : ('a, 'k) session -> float -> int -> ('a, 'k) Nvector.t -> unit

(** [reinit s ~linsolv:linsolv ~roots:(nroots, g) t0 y0 y'0]
    reinitializes the solver session [s] with a new time [t0] and new
    values for the variables [y0].  There are two optional arguments
    to change the linear solver and the set of root functions.

    The optional argument [linsolv] sets the linear solver.  If
    omitted, the current linear solver will be kept.  If a session is
    created with, say, [Dense (Some f)], and then reinitialized with
    [Dense None], then the linear solver is reset and [f] is removed
    from the session.  The same goes for all other optional callbacks.

    The optional argument [roots] sets the root functions; see {!init}
    for what each component does.  {!Ida.no_roots} may be passed in to
    turn off root finding.  If omitted, the current root functions
    will be kept.

    @ida <node5#sss:cvreinit> IDAReInit *)
val reinit :
  ('a, 'k) session
  -> ?linsolv:('a, 'k) linear_solver
  -> ?roots:(int * 'a rootsfn)
  -> float
  -> ('a, 'k) Nvector.t
  -> ('a, 'k) Nvector.t
  -> unit

(** {2:set Modifying the solver (optional input functions)} *)

(** Set the integration tolerances.

    @ida <node5#sss:idatolerances> IDASStolerances
    @ida <node5#sss:idatolerances> IDASVtolerances
    @ida <node5#sss:idatolerances> IDAWFtolerances
    @ida <node5#ss:ewtsetFn> Error weight function *)
val set_tolerances : ('a, 'k) session -> ('a, 'k) tolerance -> unit

(** [set_error_file s fname trunc] opens the file named [fname] and to
    which all messages from the default error handler are then
    directed.  If the file already exists it is either trunctated
    ([trunc] = [true]) or appended to ([trunc] = [false]).

    The error file is closed if set_error_file is called again, or
    otherwise when the session is garbage collected.

    @ida <node5#sss:optin_main> IDASetErrFile *)
val set_error_file : ('a, 'k) session -> string -> bool -> unit

(** [set_err_handler_fn s efun] specifies a custom function [efun] for
    handling error messages.  The error handler function must not fail
    -- any exceptions raised from it will be captured and discarded.

    @ida <node5#sss:optin_main> IDASetErrHandlerFn
    @ida <node5#ss:ehFn> IDAErrHandlerFn *)
val set_err_handler_fn : ('a, 'k) session -> (Sundials.error_details -> unit)
                         -> unit

(** This function restores the default error handling function. It is
    equivalent to calling IDASetErrHandlerFn with an argument of [NULL].

    @ida <node5#sss:optin_main> IDASetErrHandlerFn *)
val clear_err_handler_fn : ('a, 'k) session -> unit

(** Specifies the maximum order of the linear multistep method.

    @ida <node5#sss:optin_main> IDASetMaxOrd *)
val set_max_ord : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of steps to be taken by the solver in
    its attempt to reach the next output time.

    @ida <node5#sss:optin_main> IDASetMaxNumSteps *)
val set_max_num_steps : ('a, 'k) session -> int -> unit

(** Specifies the initial step size.

    @ida <node5#sss:optin_main> IDASetInitStep *)
val set_init_step : ('a, 'k) session -> float -> unit

(** Specifies an upper bound on the magnitude of the step size.

    @ida <node5#sss:optin_main> IDASetMaxStep *)
val set_max_step : ('a, 'k) session -> float -> unit

(** Specifies the value of the independent variable t past which the
    solution is not to proceed.  The default, if this routine is not
    called, is that no stop time is imposed.

    @ida <node5#sss:optin_main> IDASetStopTime *)
val set_stop_time : ('a, 'k) session -> float -> unit

(** Specifies the maximum number of error test failures permitted in
    attempting one step.

    @ida <node5#sss:optin_main> IDASetMaxErrTestFails *)
val set_max_err_test_fails : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver iterations
    permitted per step.

    @ida <node5#sss:optin_main> IDASetMaxNonlinIters *)
val set_max_nonlin_iters : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver convergence
    failures permitted during one step.

    @ida <node5#sss:optin_main> IDASetMaxConvFails *)
val set_max_conv_fails : ('a, 'k) session -> int -> unit

(** Specifies the safety factor used in the nonlinear convergence test.

    @ida <node5#sss:optin_main> IDASetNonlinConvCoef
    @ida <node3#ss:ivp_sol> IVP Solution *)
val set_nonlin_conv_coef : ('a, 'k) session -> float -> unit

(** Set inequality constraints on variables.  See {!Constraint}.

    @ida <node5#sss:optin_main> IDASetConstraints *)
val set_constraints : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit

(** {2:get Querying the solver (optional output functions)} *)

(** Returns the real and integer workspace sizes.

    @ida <node5#sss:optout_main> IDAGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space          : ('a, 'k) session -> int * int

(** Returns the cumulative number of internal steps taken by the
    solver.

    @ida <node5#sss:optout_main> IDAGetNumSteps *)
val get_num_steps           : ('a, 'k) session -> int

(** Returns the number of calls to the user's right-hand side
    function.

    @ida <node5#sss:optout_main> IDAGetNumResEvals *)
val get_num_res_evals       : ('a, 'k) session -> int

(** Returns the number of calls made to the linear solver's setup
    function.

    @ida <node5#sss:optout_main> IDAGetNumLinSolvSetups *)
val get_num_lin_solv_setups : ('a, 'k) session -> int

(** Returns the number of local error test failures that have
    occurred.

    @ida <node5#sss:optout_main> IDAGetNumErrTestFails *)
val get_num_err_test_fails  : ('a, 'k) session -> int

(** Returns the integration method order used during the last internal
    step.

    @ida <node5#sss:optout_main> IDAGetLastOrder *)
val get_last_order          : ('a, 'k) session -> int

(** Returns the integration method order to be used on the next
    internal step.

    @ida <node5#sss:optout_main> IDAGetCurrentOrder *)
val get_current_order       : ('a, 'k) session -> int

(** Returns the integration step size taken on the last internal step.

    @ida <node5#sss:optout_main> IDAGetLastStep *)
val get_last_step           : ('a, 'k) session -> float

(** Returns the integration step size to be attempted on the next
    internal step.

    @ida <node5#sss:optout_main> IDAGetCurrentStep *)
val get_current_step        : ('a, 'k) session -> float

(** Returns the the value of the integration step size used on the
    first step.

    @ida <node5#sss:optout_main> IDAGetActualInitStep *)
val get_actual_init_step    : ('a, 'k) session -> float

(** Returns the the current internal time reached by the solver.

    @ida <node5#sss:optout_main> IDAGetCurrentTime *)
val get_current_time        : ('a, 'k) session -> float

(* IDAGetNumStabLimOrderReds appears in the sundials 2.5.0 manual on
   p.52 but there's no such function in the implementation.  It's
   probably a leftover from earlier versions or something.

(** Returns the number of order reductions dictated by the BDF
    stability limit detection algorithm.

    @ida <node5#sss:optout_main> IDAGetNumStabLimOrderReds
    @ida <node3#s:bdf_stab> BDF stability limit detection *)
val get_num_stab_lim_order_reds : session -> int
*)

(** Returns a suggested factor by which the user's tolerances should
    be scaled when too much accuracy has been requested for some
    internal step.

    @ida <node5#sss:optout_main> IDAGetTolScaleFactor *)
val get_tol_scale_factor : ('a, 'k) session -> float

(** Returns the solution error weights at the current time.

    @ida <node5#sss:optout_main> IDAGetErrWeights
    @ida <node3#ss:ivp_sol> IVP solution (W_i) *)
val get_err_weights : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit

(** Returns the vector of estimated local errors.

    @ida <node5#sss:optout_main> IDAGetEstLocalErrors *)
val get_est_local_errors : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit

(** Aggregated integrator statistics.
 *
    @ida <node5#sss:optout_main> IDAGetIntegratorStats *)
type integrator_stats = {
    num_steps : int;
    num_res_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

(** Returns the integrator statistics as a group.

    @ida <node5#sss:optout_main> IDAGetIntegratorStats *)
val get_integrator_stats    : ('a, 'k) session -> integrator_stats

(** Convenience function that calls get_integrator_stats and prints
    the results to stdout.

    @ida <node5#sss:optout_main> IDAGetIntegratorStats *)
val print_integrator_stats  : ('a, 'k) session -> unit


(** Returns the number of nonlinear (functional or Newton) iterations
    performed.

    @ida <node5#sss:optout_main> IDAGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('a, 'k) session -> int

(** Returns the number of nonlinear convergence failures that have
    occurred.

    @ida <node5#sss:optout_main> IDAGetNumNonlinSolvConvFails *)
val get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int

(** [nniters, nncfails = get_nonlin_solv_stats s] obtains both the numbers of
    nonlinear iterations performed [nniters] and of nonlinear convergence
    failures that have occurred [nncfails].

    @ida <node5#sss:optout_main> IDAGetNonlinSolvStats *)
val get_nonlin_solv_stats : ('a, 'k) session -> int *int

(** {2:roots Additional root finding functions} *)

(** [set_root_direction s dir] specifies the direction of
    zero-crossings to be located and returned. [dir] may contain one
    entry of type {!Ida.root_direction} for each root function.

    @ida <node5#sss:optin_root> IDASetRootDirection *)
val set_root_direction : ('a, 'k) session -> Sundials.RootDirs.d array -> unit

(** Like {!set_root_direction} but specifies a single direction of
    type {!Ida.root_direction} for all root functions.

  @ida <node5#sss:optin_root> IDASetRootDirection *)
val set_all_root_directions : ('a, 'k) session -> Sundials.RootDirs.d -> unit

(**
  Disables issuing a warning if some root function appears to be identically
  zero at the beginning of the integration.

  @ida <node5#sss:optin_root> IDASetNoInactiveRootWarn *)
val set_no_inactive_root_warn : ('a, 'k) session -> unit

(**
  Fills an array showing which functions were found to have a root.

  @ida <node5#sss:optout_root> IDAGetRootInfo *)
val get_root_info : ('a, 'k) session -> Sundials.Roots.t -> unit

(** Returns the cumulative number of calls made to the user-supplied
    root function g.

    @ida <node5#sss:optout_root> IDAGetNumGEvals *)
val get_num_g_evals : ('a, 'k) session -> int

(** {2:exceptions Exceptions} *)

(** Raised on missing or illegal solver inputs. Also raised if an element
    of the error weight vector becomes zero during time stepping, or the
    linear solver initialization function failed, or a root was found both at
    [t] and very near [t].
 
    @ida <node5#sss:idasolve> IDA_ILL_INPUT *)
exception IllInput

(** The requested time could not be reached in [mxstep] internal steps.
    See {!set_max_num_steps}

    @ida <node5#sss:idasolve> IDA_TOO_MUCH_WORK *)
exception TooMuchWork

(** The requested accuracy could not be satisfied.

    @ida <node5#sss:idasolve> IDA_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** Too many error test failures within a step or at the minimum step size.
    See {!set_max_err_test_fails} and {!set_min_step}.

    @ida <node5#sss:idasolve> IDA_ERR_FAIL *)
exception ErrFailure                

(** Too many convergence test failures within a step or at the minimum step
    size, or Newton convergence failed.
    See {!set_max_conv_fails} and {!set_min_step}.

    @ida <node5#sss:idasolve> IDA_CONV_FAIL *)
exception ConvergenceFailure        

(** Linear solver initialization failed.

    @ida <node5#sss:idasolve> IDA_LINIT_FAIL *)
exception LinearInitFailure         

(** Linear solver setup failed in an unrecoverable manner.

    @ida <node5#sss:idasolve> IDA_LSETUP_FAIL *)
exception LinearSetupFailure        

(** Linear solver solution failed in an unrecoverable manner.

    @ida <node5#sss:idasolve> IDA_LSOLVE_FAIL *)
exception LinearSolveFailure        

(** The residual function failed in an unrecoverable manner.

    @ida <node5#ss:idasolve> IDA_RES_FAIL *)
exception ResFuncFailure

(** The residual function had a recoverable error when first called.

    @ida <node5#ss:idacalcic> IDA_FIRST_RES_FAIL *)
exception FirstResFuncFailure       

(** Too many convergence test failures, or unable to estimate the initial step
    size, due to repeated recoverable errors in the residual function.

    @ida <node5#sss:idasolve> IDA_REP_RES_ERR *)
exception RepeatedResFuncFailure

(** The rootfinding function failed.

    @ida <node5#sss:idasolve> IDA_RTFUNC_FAIL *)
exception RootFuncFailure           

(** No solution satisfying the inequality constraints could be found.
 
    @ida <node5#ss:idacalcic> IDA_CONSTR_FAIL *)
exception ConstraintFailure

(** Linesearch could not find a solution with a step larger than steptol in
    weighted RMS norm.
 
    @ida <node5#ss:idacalcic> IDA_LINESEARCH_FAIL *)
exception LinesearchFailure

(** A recoverable error occurred in a callback but no recovery was possible.
 
    @ida <node5#ss:idacalcic> IDA_NO_RECOVERY *)
exception NoRecovery

(** A component of the error weight vector, either for the input value or a
    corrected value, is zero.
 
    @ida <node5#ss:idacalcic> IDA_BAD_EWT *)
exception BadEwt

(** Raised by {!get_dky} for invalid order values.
 
    @ida <node5#ss:optional_dky> IDAGetDky (IDA_BAD_K) *)
exception BadK

(** Raised by {!get_dky} for invalid time values.

    @ida <node5#ss:optional_dky> IDAGetDky (IDA_BAD_T) *)
exception BadT

(** Raised by {!get_dky} on an invalid derivative vector.

    @ida <node5#ss:optional_dky> IDAGetDky (IDA_BAD_DKY) *)
exception BadDky

(** Variable ids are required but not set. See {!set_id}. *)
exception IdNotSet

