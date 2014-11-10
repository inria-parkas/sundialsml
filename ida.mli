(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
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
      {- {{:#solver}Solution}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}
      {- {{:#exceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** A session with the IDA solver.

    An example session with Ida ({openfile ida_skel.ml}): {[
#include "examples/ocaml/skeletons/ida_skel.ml"
    ]}

    @ida <node5#ss:skeleton_sim> Skeleton of main program
 *)
type ('d, 'k) session = ('d, 'k) Ida_impl.session

(** Alias for sessions based on serial nvectors. *)
type serial_session = (RealArray.t, Nvector_serial.kind) session

(** {2:linear Linear solvers} *)

(** Linear solvers used by Ida.

    @ida <node5#sss:lin_solv_init> Linear Solver Specification Functions *)
type ('data, 'kind) linear_solver = ('data, 'kind) Ida_impl.linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type serial_linear_solver =
      (Nvector_serial.data, Nvector_serial.kind) linear_solver

(** Workspaces with two temporary vectors. *)
type 'd double = 'd * 'd

(** Workspaces with three temporary vectors. *)
type 'd triple = 'd * 'd * 'd

(** Arguments common to Jacobian callback functions.    
   
    @ida <node5#ss:djacFn> IDADlsDenseJacFn
    @ida <node5#ss:bjacFn> IDADlsBandJacFn
    @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn
    @ida <node5#ss:psolveFn> IDASpilsPrecSolveFn
    @ida <node5#ss:precondFn> IDASpilsPrecSetupFn *)
type ('t, 'd) jacobian_arg =
  {
    jac_t    : float;        (** The independent variable. *)
    jac_y    : 'd;           (** The dependent variable vector. *)
    jac_y'   : 'd;           (** The derivative vector (i.e.
                                  {% $\frac{\mathrm{d}y}{\mathrm{d}t}$%}). *)
    jac_res  : 'd;           (** The current value of the residual vector. *)
    jac_coef : float;        (** The coefficient $c_j$ in
                                 {% $J = \frac{\partial F}{\partial y} + c_j \frac{\partial F}{\partial\dot{y}}$%}. *)
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
  sig (* {{{ *)

    (** Callback functions that compute dense approximations to a Jacobian
        matrix. In the call [dense_jac_fn arg jac], [arg] is a {!jacobian_arg}
        with three work vectors and the computed Jacobian must be stored
        in [jac].

        The callback should load the [(i,j)]th entry of [jac] with
        {% $\frac{\partial F_i}{\partial y_j} + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%},
        i.e., the partial derivative of the [i]th equation with respect to
        the [j]th variable, evaluated at the values of [t], [y], and [y']
        obtained from [arg]. Only nonzero elements need be loaded into [jac].

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

    (** A direct linear solver on dense matrices using LAPACK. See {!dense}.
        Only available if {!Sundials.lapack_enabled}.

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
        {% $\frac{\partial F_i}{\partial y_j} + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%},
        i.e., the partial derivative of the [i]th equation with respect to
        the [j]th variable, evaluated at the values of [t] and [y] obtained
        from [arg]. Only nonzero elements need be loaded into [jac].

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
        Only available if {!Sundials.lapack_enabled}.

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
  end (* }}} *)

(** Scaled Preconditioned Iterative Linear Solvers.

    @ida <node5#sss:optin_spils> Iterative linear solvers optional input functions.
    @ida <node5#sss:optout_spils> Iterative linear solvers optional output functions. *)
module Spils :
  sig (* {{{ *)
    (** {3:precond Preconditioners} *)

    (** Callback functions that solve a linear system involving a
        preconditioner matrix.
        In the call [prec_solve_fn jac r z delta],
        [jac] is a {!jacobian_arg} with one work vector,
        [r] is the right-hand side vector,
        [z] is computed to solve {% $Pz = r$%},
        and [delta] is the input tolerance.
        $P$ is a preconditioner matrix, which approximates, however crudely,
        the Jacobian matrix
        {% $\frac{\partial F}{\partial y} + \mathtt{arg.jac\_coef}\frac{\partial F}{\partial\dot{y}}$%}.
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
    type 'd prec_solve_fn =
      ('d, 'd) jacobian_arg
      -> 'd
      -> 'd
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
    type 'd prec_setup_fn = ('d triple, 'd) jacobian_arg -> unit

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
    type 'd jac_times_vec_fn =
      ('d double, 'd) jacobian_arg
      -> 'd
      -> 'd
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
    type ('d, 'k) preconditioner = ('d, 'k) Ida_impl.SpilsTypes.preconditioner

    (** No preconditioning.  *)
    val prec_none : ('d, 'k) preconditioner

    (** Left preconditioning. {% $Pz = r$%}, where $P$ approximates, perhaps
        crudely,
        {% $J = \frac{\partial F}{\partial y} + c_j\frac{\partial F}{\partial\dot{y}}$%}. *)
    val prec_left :
      ?setup:'d prec_setup_fn
      -> ?jac_times_vec:'d jac_times_vec_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

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
      -> ('d, 'k) preconditioner -> ('d, 'k) linear_solver

    (** Krylov iterative solver using the scaled preconditioned biconjugate
        stabilized (Bi-CGStab) method.
        In the call [spbcg ~maxl:maxl prec], [maxl] is the maximum dimension of
        the Krylov subspace (defaults to 5), and [prec] is a {!preconditioner}.

        @ida <node5#sss:lin_solv_init> IDASpbcg
        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#sss:optin_spils> IDASpilsSetMaxl *)
    val spbcg : ?maxl:int -> ('d, 'k) preconditioner -> ('d, 'k) linear_solver

    (** Krylov iterative with the scaled preconditioned transpose-free
        quasi-minimal residual (SPTFQMR) method.
        In the call [sptfqmr ~maxl:maxl prec], [maxl] is the maximum dimension
        of the Krylov subspace (defaults to 5), and [prec] is a
        {!preconditioner}.

        @ida <node5#sss:lin_solv_init> IDASptfqmr
        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#sss:optin_spils> IDASpilsSetMaxl *)
    val sptfqmr : ?maxl:int -> ('d, 'k) preconditioner -> ('d, 'k) linear_solver

    (** {3:set Solver parameters} *)

    (** Sets the Gram-Schmidt orthogonalization to be used with the
        Spgmr {!linear_solver}.

        @ida <node5#sss:optin_spils> IDASpilsSetGSType *)
    val set_gs_type : ('d, 'k) session -> Spils.gramschmidt_type -> unit

    (** Sets the factor by which the Krylov linear solver's convergence test
        constant is reduced from the Newton iteration test constant.
        This factor must be >= 0; passing 0 specifies the default (0.05).

        @ida <node5#sss:optin_spils> IDASpilsSetEpsLin *)
    val set_eps_lin : ('d, 'k) session -> float -> unit

    (** Resets the maximum Krylov subspace dimension for the Bi-CGStab and
        TFQMR methods. A value <= 0 specifies the default (5.0).

        @ida <node5#sss:optin_spils> IDASpilsSetMaxl *)
    val set_maxl : ('d, 'k) session -> int -> unit

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the spils
        linear solver.

        @ida <node5#sss:optout_spils> IDASpilsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space       : ('d, 'k) session -> int * int

    (** Returns the cumulative number of linear iterations.

        @ida <node5#sss:optout_spils> IDASpilsGetNumLinIters *)
    val get_num_lin_iters    : ('d, 'k) session -> int

    (** Returns the cumulative number of linear convergence failures.

        @ida <node5#sss:optout_spils> IDASpilsGetNumConvFails *)
    val get_num_conv_fails   : ('d, 'k) session -> int

    (** Returns the number of calls to the setup function.

        @ida <node5#sss:optout_spils> IDASpilsGetNumPrecEvals *)
    val get_num_prec_evals   : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the preconditioner solve
        function.

        @ida <node5#sss:optout_spils> IDASpilsGetNumPrecSolves *)
    val get_num_prec_solves  : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        function.

        @ida <node5#sss:optout_spils> IDASpilsGetNumJtimesEvals *)
    val get_num_jtimes_evals : ('d, 'k) session -> int

    (** Returns the number of calls to the residual callback for
        finite difference Jacobian-vector product approximation. This counter is
        only updated if the default difference quotient function is used.

        @ida <node5#sss:optout_spils> IDASpilsGetNumResEvals *)
    val get_num_res_evals    : ('d, 'k) session -> int

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
       ('d,'k) session
       -> ?setup:'d prec_setup_fn
       -> 'd prec_solve_fn
       -> unit

    (** Change the Jacobian-times-vector function.

        @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
        @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn *)
    val set_jac_times_vec_fn :
      ('d,'k) session
      -> 'd jac_times_vec_fn
      -> unit

    (** Remove a Jacobian-times-vector function and use the default
        implementation.

        @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
        @ida <node5#ss:jtimesFn> IDASpilsJacTimesVecFn *)
    val clear_jac_times_vec_fn : ('d, 'k) session -> unit
  end (* }}} *)

(** Alternate Linear Solvers.

    @ida <node8#s:new_linsolv> Providing Alternate Linear Solver Modules *)
module Alternate :
  sig (* {{{ *)

    (** Functions that initialize linear solver data, like counters and
        statistics.

        Raising any exception in this function (including
        {!Sundials.RecoverableFailure}) is treated as an unrecoverable error.

        @ida <node8#SECTION00810000000000000000> linit *)
    type ('data, 'kind) linit = ('data, 'kind) session -> unit

    (** Functions that prepare the linear solver for subsequent calls to
        {{!callbacks}lsolve}. The call [lsetup s y y' res tmp] has as
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
        $M$ should approximate {% $J = \frac{\partial F}{\partial y} + c_j\frac{\partial F}{\partial \dot{y}}$%},
        and $c_j$ is available through {!get_cj}.
        The call [lsolve s b weight ycur y'cur rescur] has as arguments:

        - [s], the solver session,
        - [b], for returning the calculated solution,
        - [weight], the error weights,
        - [ycur], the solver's current approximation to $y(t_n)$,
        - [ycur'], the solver's current approximation to {% $\dot{y}(t_n)$%},
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
  end (* }}} *)

(** {2:tols Tolerances} *)

(** Functions that set the multiplicative error weights for use in the weighted
    RMS norm. The call [efun y ewt] takes the dependent variable vector [y] and
    fills the error-weight vector [ewt] with positive values or raises
    {!Sundials.NonPositiveEwt}. Other exceptions are eventually propagated, but
    should be avoided ([efun] is not allowed to abort the solver). *)
type 'data error_fun = 'data -> 'data -> unit

(** Tolerance specifications. *)
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

(** Residual functions that define a DAE problem. They are passed four
 * arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$,
    - [y'], the vector of dependent-variable derivatives, i.e.,
            {% $\dot{y} = \frac{\mathrm{d}y}{\mathrm{d}t}$%}, and,
    - [r] a vector for storing the residual value, {% $F(t, y, \dot{y})$%}.

    Within the function, raising a {!Sundials.RecoverableFailure} exception
    indicates a recoverable error. Any other exception is treated as an
    unrecoverable error.

    {warning [y], [y'], and [r] should not be accessed after the function
             returns.}

    @ida <node5#ss:resFn> IDAResFn *)
type 'd resfn = float -> 'd -> 'd -> 'd -> unit

(** Called by the solver to calculate the values of root functions. These
    ‘zero-crossings’ are used to detect significant events. The function is
    passed four arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$,
    - [y'], the vector of dependent-variable derivatives, i.e.,
            {% $\dot{y} = \frac{\mathrm{d}y}{\mathrm{d}t}$%}, and,
    - [gout], a vector for storing the value of {% $g(t, y, \dot{y})$%}.

    {warning [y], [y'], and [gout] should not be accessed after the function
             has returned.}

    @ida <node5#ss:rootFn> IDARootFn *)
type 'd rootsfn = float -> 'd -> 'd -> RealArray.t -> unit

(** Creates and initializes a session with the solver. The call
    {[init linsolv tol f ~varid:varid ~roots:(nroots, g) t0 y0 y'0]} has
    as arguments:
    - [linsolv], the linear solver to use,
    - [tol],     the integration tolerances,
    - [f],       the DAE residual function,
    - [varid],   optionally classifies variables as algebraic or differential,
    - [nroots],  the number of root functions,
    - [g],       the root function ([(nroots, g)] defaults to {!no_roots}),
    - [t0],      the initial value of the independent variable,
    - [y0],      a vector of initial values that also determines the number
                 of equations, and,
    - [y'0],     the initial values for
                 {% $\dot{y} = \frac{\mathrm{d}y}{\mathrm{d}t}$%}.

    This function does everything necessary to initialize a session, i.e.,
    it makes the calls referenced below. The {!solve_normal} and
    {!solve_one_step} functions may be called directly.

    @ida <node5#sss:idainit>       IDACreate/IDAInit
    @ida <node5#ss:idarootinit>    IDARootInit
    @ida <node5#sss:lin_solv_init> Linear solvers
    @ida <node5#sss:idatolerances> IDASStolerances
    @ida <node5#sss:idatolerances> IDASVtolerances
    @ida <node5#sss:idatolerances> IDAWFtolerances
    @ida <node5#ss:ewtsetFn>       IDAEwtFn
    @ida <node5#sss:idasetid>      IDASetId *)
val init :
    ('d, 'kind) linear_solver
    -> ('d, 'kind) tolerance
    -> 'd resfn
    -> ?varid:('d, 'kind) Nvector.t
    -> ?roots:(int * 'd rootsfn)
    -> float
    -> ('d, 'kind) Nvector.t
    -> ('d, 'kind) Nvector.t
    -> ('d, 'kind) session

(** A convenience value for signalling that there are no roots to monitor. *)
val no_roots : (int * 'd rootsfn)

(** {3:calcic Initial Condition Calculation} *)

(** Symbolic names for constants used when calculating initial values or
    supressing local error tests. See {!calc_ic_ya_yd'} and
    {!set_suppress_alg}.

    @ida <node5#sss:idasetid> IDASetId *)
module VarId :
  sig (* {{{ *)
    (** The constant [0.0]. *)
    val algebraic : float

    (** The constant [1.0]. *)
    val differential : float

    (** For pattern-matching on constraints. See {!of_float}. *)
    type t =
    | Algebraic    (** Residual functions must not depend on the derivatives
                       of algebraic variables. *)
    | Differential (** Residual functions may depend on the derivatives of
                       differential variables. *)

    (** Map id values to floating-point constants. *)
    val to_float : t -> float

    (** Map floating-point constants to id values.
        
        @raise Invalid_argument The given value is not a legal id. *)
    val of_float : float -> t
  end (* }}} *)

(** Class components of the state vector as either algebraic or differential.
    These classifications are required by {!calc_ic_ya_yd'} and
    {!set_suppress_alg}. See also {!VarId}.

    @ida <node5#sss:optin_main> IDASetId *)
val set_id : ('d, 'k) session -> ('d,'k) Nvector.t -> unit

(** Indicates whether or not to ignore algebraic variables in the local error
    test. When ignoring algebraic variables ([true]), a [varid] vector must be
    specified either in the call or by a prior call to {!init} or {!set_id}.
    Suppressing local error tests for algebraic variables is {i discouraged}
    for DAE systems of index 1 and {i encouraged} for systems of index 2 or
    more. 

    @ida <node5#sss:optin_main> IDASetId
    @ida <node5#sss:optin_main> IDASetSuppressAlg *)
val set_suppress_alg : ('d, 'k) session
                       -> ?varid:('d, 'k) Nvector.t -> bool -> unit

(** Computes the initial state vector for certain index-one problems.
    All components of $y$ are computed, using {% $\dot{y}$%}, to satisfy
    the constraint {% $F(t_0, y_0, \dot{y}_0) = 0$%}. If given, the
    [~y] vector is filled with the corrected values. The last argument is
    the first vale of $t$ at which a solution will be requested.
    A {!reinit} is required before calling this function after
    {!solve_normal} or {!solve_one_step}.

    @ida <node5#ss:idacalcic> IDACalcIC (IDA_Y_INIT)
    @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
    @raise IdNotSet Variable ids have not been set (see {!set_id}). *)
val calc_ic_y : ('d, 'k) session -> ?y:('d, 'k) Nvector.t -> float -> unit

(** Computes the algebraic components of the initial state and derivative
    vectors for certain index-one problems.
    The elements of $y$ and {% $\dot{y}$%} marked as algebraic are computed,
    using {% $\dot{y}$%}, to satisfy the constraint
    {% $F(t_0, y_0, \dot{y}_0) = 0$%}.
    The variable ids must be given in [~varid] or by a prior call to {!init} or
    {!set_id}.
    If given, the [~y] and [~y'] vectors are filled with the corrected values.
    The last argument is the first vale of $t$ at which a solution will be
    requested. A {!reinit} is required before calling this function after
    {!solve_normal} or {!solve_one_step}.

    @ida <node5#ss:idacalcic> IDACalcIC (IDA_YA_YDP_INIT)
    @ida <node5#sss:optin_main> IDASetId
    @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
    @raise IdNotSet Variable ids have not been set (see {!set_id}). *)
val calc_ic_ya_yd' :
  ('d, 'k) session
  -> ?y:('d, 'k) Nvector.t
  -> ?y':('d, 'k) Nvector.t
  -> ?varid:('d, 'k) Nvector.t
  -> float
  -> unit

(** Returns the number of backtrack operations during {!calc_ic_ya_yd'} or
    {!calc_ic_y}.

    @ida <node5#sss:optout_iccalc> IDAGetNumBcktrackOps *)
val get_num_backtrack_ops : ('d, 'k) session -> int

(** {2:solver Solution} *)

(** Values returned by the step functions. Failures are indicated by
    exceptions.

    @ida <node5#sss:ida> IDASolve *)
type solver_result =
  | Success             (** The solution was advanced. {cconst IDA_SUCCESS} *)
  | RootsFound          (** A root was found. See {!get_root_info}.
                            {cconst IDA_ROOT_RETURN} *)
  | StopTimeReached     (** The stop time was reached. See {!set_stop_time}.
                            {cconst IDA_TSTOP_RETURN} *)

(** Integrates a DAE system over an interval. The call
    [tret, r = solve_normal s tout yout y'out] has as arguments
    - [s], a solver session,
    - [tout] the next time at which the solution is desired,
    - [yout], a vector to store the computed solution, and,
    - [y'out], a vector to store the computed derivative.

    It returns [tret], the time reached by the solver, which will be equal to
    [tout] if no errors occur, and, [r], a {!solver_result}.

    @ida <node5#sss:idasolve> IDASolve (IDA_NORMAL)
    @raise IllInput Missing or illegal solver inputs.
    @raise TooMuchWork The requested time could not be reached in [mxstep] internal steps.
    @raise TooMuchAccuracy The requested accuracy could not be satisfied.
    @raise ErrFailure Too many error test failures within a step or at the minimum step size.
    @raise ConvergenceFailure Too many convergence test failures within a step or at the minimum step size.
    @raise LinearInitFailure Linear solver initialization failed.
    @raise LinearSetupFailure Linear solver setup failed unrecoverably.
    @raise LinearSolveFailure Linear solver solution failed unrecoverably.
    @raise ConstraintFailure Inequality constraints were violated and recovery is not possible.
    @raise RepeatedResFuncFailure The residual function repeatedly returned a recoverable error but the solver could not recover.
    @raise ResFuncFailure The residual function failed unrecoverably.
    @raise RootFuncFailure Failure in the rootfinding function [g]. *)
val solve_normal : ('d, 'k) session -> float
                   -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t
                   -> float * solver_result

(** Like {!solve_normal} but returns after one internal solver step.

    @ida <node5#sss:idasolve> IDASolve (IDA_ONE_STEP) *)
val solve_one_step : ('d, 'k) session -> float
                     -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t
                     -> float * solver_result

(** Returns the interpolated solution or derivatives.
    [get_dky s dky t k] computes the [k]th derivative of the function at time
    [t], i.e., {% $\frac{d^\mathtt{k}y(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
    and stores it in [dky]. The arguments must satisfy
    {% $t_n - h_u \leq \mathtt{t} \leq t_n$%}—where $t_n$
    denotes {!get_current_time} and $h_u$ denotes {!get_last_step},—
    and {% $0 \leq \mathtt{k} \leq q_u$%}—where $q_u$ denotes
    {!get_last_order}.

    This function may only be called after a successful return from either
    {!solve_normal} or {!solve_one_step}.

    @ida <node5#sss:optin_root> IDAGetDky
    @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
    @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
val get_dky : ('d, 'k) session -> ('d, 'k) Nvector.t -> float -> int -> unit

(** Reinitializes the solver with new parameters and state values. The
    values of the independent variable, i.e., the simulation time, the
    state variables, and the derivatives must be given.
    If the argument [~linsolv] is not given, the current linear solver
    remains unchanged. The argument [~roots] works similarly; pass
    {!no_roots} to disable root finding.

    @ida <node5#sss:cvreinit> IDAReInit *)
val reinit :
  ('d, 'k) session
  -> ?linsolv:('d, 'k) linear_solver
  -> ?roots:(int * 'd rootsfn)
  -> float
  -> ('d, 'k) Nvector.t
  -> ('d, 'k) Nvector.t
  -> unit

(** {2:set Modifying the solver (optional input functions)} *)

(** Set the integration tolerances.

    @ida <node5#sss:idatolerances> IDASStolerances
    @ida <node5#sss:idatolerances> IDASVtolerances
    @ida <node5#sss:idatolerances> IDAWFtolerances
    @ida <node5#ss:ewtsetFn>       IDAEwtFn *)
val set_tolerances : ('d, 'k) session -> ('d, 'k) tolerance -> unit

(** Opens the named file to receive messages from the default error handler.
    If the file already exists it is either truncated ([true]) or extended
    ([false]).
    The file is closed if the function is called again or when the session is
    garbage collected.

    @ida <node5#sss:optin_main> IDASetErrFile *)
val set_error_file : ('d, 'k) session -> string -> bool -> unit

(** Specifies a custom function for handling error messages.
    This function must not fail: any exceptions are trapped and discarded.

    @ida <node5#sss:optin_main> IDASetErrHandlerFn
    @ida <node5#ss:ehFn> IDAErrHandlerFn *)
val set_err_handler_fn : ('d, 'k) session -> (error_details -> unit) -> unit

(** Restores the default error handling function.

    @ida <node5#sss:optin_main> IDASetErrHandlerFn *)
val clear_err_handler_fn : ('d, 'k) session -> unit

(** Specifies the maximum order of the linear multistep method.

    @ida <node5#sss:optin_main> IDASetMaxOrd *)
val set_max_ord : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of steps taken in attempting to reach
    a given output time.

    @ida <node5#sss:optin_main> IDASetMaxNumSteps *)
val set_max_num_steps : ('d, 'k) session -> int -> unit

(** Specifies the initial step size.

    @ida <node5#sss:optin_main> IDASetInitStep *)
val set_init_step : ('d, 'k) session -> float -> unit

(** Specifies an upper bound on the magnitude of the step size.

    @ida <node5#sss:optin_main> IDASetMaxStep *)
val set_max_step : ('d, 'k) session -> float -> unit

(** Limits the value of the independent variable [t] when solving.
    By default no stop time is imposed.

    @ida <node5#sss:optin_main> IDASetStopTime *)
val set_stop_time : ('d, 'k) session -> float -> unit

(** Specifies the maximum number of error test failures permitted in attempting
    one step.

    @ida <node5#sss:optin_main> IDASetMaxErrTestFails *)
val set_max_err_test_fails : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver iterations permitted per
    step.

    @ida <node5#sss:optin_main> IDASetMaxNonlinIters *)
val set_max_nonlin_iters : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver convergence failures
    permitted during one step.

    @ida <node5#sss:optin_main> IDASetMaxConvFails *)
val set_max_conv_fails : ('d, 'k) session -> int -> unit

(** Specifies the safety factor used in the nonlinear convergence test.

    @ida <node5#sss:optin_main> IDASetNonlinConvCoef
    @ida <node3#ss:ivp_sol> IVP Solution *)
val set_nonlin_conv_coef : ('d, 'k) session -> float -> unit

(** Specifies a vector defining inequality constraints for each
    component of the solution vector [u].  See {!Sundials.Constraint}.

    @ida <node5#sss:optin_main> IDASetConstraints *)
val set_constraints : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** {2:get Querying the solver (optional output functions)} *)

(** Returns the sizes of the real and integer workspaces.

    @ida <node5#sss:optout_main> IDAGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space          : ('d, 'k) session -> int * int

(** Returns the cumulative number of internal solver steps.

    @ida <node5#sss:optout_main> IDAGetNumSteps *)
val get_num_steps           : ('d, 'k) session -> int

(** Returns the number of calls to the residual function.

    @ida <node5#sss:optout_main> IDAGetNumResEvals *)
val get_num_res_evals       : ('d, 'k) session -> int

(** Returns the number of calls made to the linear solver's setup function.

    @ida <node5#sss:optout_main> IDAGetNumLinSolvSetups *)
val get_num_lin_solv_setups : ('d, 'k) session -> int

(** Returns the number of local error test failures that have occurred.

    @ida <node5#sss:optout_main> IDAGetNumErrTestFails *)
val get_num_err_test_fails  : ('d, 'k) session -> int

(** Returns the integration method order used during the last internal step.

    @ida <node5#sss:optout_main> IDAGetLastOrder *)
val get_last_order          : ('d, 'k) session -> int

(** Returns the integration method order to be used on the next internal step.

    @ida <node5#sss:optout_main> IDAGetCurrentOrder *)
val get_current_order       : ('d, 'k) session -> int

(** Returns the integration step size taken on the last internal step.

    @ida <node5#sss:optout_main> IDAGetLastStep *)
val get_last_step           : ('d, 'k) session -> float

(** Returns the integration step size to be attempted on the next internal step.

    @ida <node5#sss:optout_main> IDAGetCurrentStep *)
val get_current_step        : ('d, 'k) session -> float

(** Returns the the value of the integration step size used on the first step.

    @ida <node5#sss:optout_main> IDAGetActualInitStep *)
val get_actual_init_step    : ('d, 'k) session -> float

(** Returns the the current internal time reached by the solver.

    @ida <node5#sss:optout_main> IDAGetCurrentTime *)
val get_current_time        : ('d, 'k) session -> float

(** Returns a suggested factor by which the user's tolerances should be scaled
    when too much accuracy has been requested for some internal step.

    @ida <node5#sss:optout_main> IDAGetTolScaleFactor *)
val get_tol_scale_factor : ('d, 'k) session -> float

(** Returns the solution error weights at the current time.

    @ida <node5#sss:optout_main> IDAGetErrWeights
    @ida <node3#ss:ivp_sol> IVP solution (W_i) *)
val get_err_weights : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Returns the vector of estimated local errors.

    @ida <node5#sss:optout_main> IDAGetEstLocalErrors *)
val get_est_local_errors : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Summaries of integrator statistics. *)
type integrator_stats = {
    num_steps : int;
      (** Cumulative number of internal solver steps. *)
    num_res_evals : int;
      (** Number of calls to the residual function. *)
    num_lin_solv_setups : int;
      (** Number of setups calls to the linear solver. *)
    num_err_test_fails : int;
      (** Number of local error test failures. *)
    last_order : int;
      (** Integration method order used in the last internal step. *)
    current_order : int;
      (** Integration method order to be used in the next internal step. *)
    actual_init_step : float;
      (** Integration step sized used on the first step. *)
    last_step : float;
      (** Integration step size of the last internal step. *)
    current_step : float;
      (** Integration step size to attempt on the next internal step. *)
    current_time : float
      (** Current internal time reached by the solver. *)
  }

(** Returns the integrator statistics as a group.

    @ida <node5#sss:optout_main> IDAGetIntegratorStats *)
val get_integrator_stats    : ('d, 'k) session -> integrator_stats

(** Prints the integrator statistics on the given channel.

    @ida <node5#sss:optout_main> IDAGetIntegratorStats *)
val print_integrator_stats  : ('d, 'k) session -> out_channel -> unit


(** Returns the number of nonlinear (functional or Newton) iterations performed.

    @ida <node5#sss:optout_main> IDAGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('d, 'k) session -> int

(** Returns the number of nonlinear convergence failures that have occurred.

    @ida <node5#sss:optout_main> IDAGetNumNonlinSolvConvFails *)
val get_num_nonlin_solv_conv_fails : ('d, 'k) session -> int

(** Returns both the numbers of nonlinear iterations performed [nniters] and
    nonlinear convergence failures [nncfails].

    @ida <node5#sss:optout_main> IDAGetNonlinSolvStats
    @return ([nniters], [nncfails]) *)
val get_nonlin_solv_stats : ('d, 'k) session -> int *int

(** {2:roots Additional root-finding functions} *)

(** [set_root_direction s dir] specifies the direction of zero-crossings to be
    located and returned. [dir] may contain one entry for each root function.

    @ida <node5#sss:optin_root> IDASetRootDirection *)
val set_root_direction : ('d, 'k) session -> RootDirs.d array -> unit

(** Like {!set_root_direction} but specifies a single direction for all root
    functions.

    @ida <node5#sss:optin_root> IDASetRootDirection *)
val set_all_root_directions : ('d, 'k) session -> RootDirs.d -> unit

(** Disables issuing a warning if some root function appears to be identically
    zero at the beginning of the integration.

    @ida <node5#sss:optin_root> IDASetNoInactiveRootWarn *)
val set_no_inactive_root_warn : ('d, 'k) session -> unit

(** Returns the number of root functions. *)
val get_num_roots : ('d, 'k) session -> int

(** Fills an array showing which functions were found to have a root.

    @ida <node5#sss:optout_root> IDAGetRootInfo *)
val get_root_info : ('d, 'k) session -> Roots.t -> unit

(** Returns the cumulative number of calls made to the user-supplied root
    function g.

    @ida <node5#sss:optout_root> IDAGetNumGEvals *)
val get_num_g_evals : ('d, 'k) session -> int

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

(** Too many error test failures within a step. See {!set_max_err_test_fails}.

    @ida <node5#sss:idasolve> IDA_ERR_FAIL *)
exception ErrFailure                

(** Too many convergence test failures within a step,
    or Newton convergence failed. See {!set_max_conv_fails}.

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

(** Variable ids are required but not set. See {!set_id}. *)
exception IdNotSet

