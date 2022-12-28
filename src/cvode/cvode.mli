(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(***********************************************************************)
(* Parts of the comment text are taken directly from:                  *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Variable-step solution of ODE initial value problems with
    zero-crossing detection.

    This module solves numerically problems of the form
    {% $\dot{y} = f(t, y)$%}, {% $y(t_0) = y_0$%}.

    This documented interface is structured as follows.
    {ol
      {- {{:#linear}Linear solvers}}
      {- {{:#tols}Tolerances}}
      {- {{:#solver}Solver initialization and use}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}
      {- {{:#exceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** A session with the CVODE solver.

    An example session with Cvode ({openfile cvode_skel.ml}): {[
#include "../../examples/ocaml/skeletons/cvode_skel.ml"
    ]}

    @cvode <Usage/index.html#a-skeleton-of-the-users-main-program> Skeleton of main program *)
type ('d, 'k) session = ('d, 'k) Cvode_impl.session

(** Alias for sessions based on serial nvectors. *)
type 'k serial_session = (Nvector_serial.data, 'k) session
                         constraint 'k = [>Nvector_serial.kind]

(** {2:linear Linear solvers} *)

(** Linear solvers used by Cvode.

    @cvode <Usage/index.html#linear-solver-interface-functions> Linear Solver Interface Functions *)
type ('data, 'kind) linear_solver = ('data, 'kind) Cvode_impl.linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type 'kind serial_linear_solver =
  (Nvector_serial.data, 'kind) linear_solver
  constraint 'kind = [>Nvector_serial.kind]

(** Workspaces with two temporary vectors. *)
type 'd double = 'd * 'd

(** Workspaces with three temporary vectors. *)
type 'd triple = 'd * 'd * 'd

(** Arguments common to Jacobian callback functions.

    @cvode CVLsJacFn
    @cvode CVLsJacTimesVecFn
    @cvode CVLsPrecSolveFn
    @cvode CVLsPrecSetupFn *)
type ('t, 'd) jacobian_arg = ('t, 'd) Cvode_impl.jacobian_arg =
  {
    jac_t   : float;        (** The independent variable. *)
    jac_y   : 'd;           (** The dependent variable vector. *)
    jac_fy  : 'd;           (** The derivative vector
                            (i.e., {% $\frac{\mathrm{d}y}{\mathrm{d}t}$%}). *)
    jac_tmp : 't            (** Workspace data. *)
  }

(** Right-hand side functions for calculating ODE derivatives. They are passed
    three arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
    - [y'], a vector for storing the value of $f(t, y)$.

    Within the function, raising a {!Sundials.RecoverableFailure} exception
    indicates a recoverable error. Any other exception is treated as an
    unrecoverable error.

    {warning [y] and [y'] should not be accessed after the function
             returns.}

    @cvode CVRhsFn *)
type 'd rhsfn = float -> 'd -> 'd -> unit

(** Diagonal approximation of Jacobians by difference quotients. *)
module Diag : sig (* {{{ *)
  (** A linear solver based on Jacobian approximation by difference
      quotients.

      @cvode CVDiag *)
  val solver : ('data, 'kind) linear_solver

  (** Returns the sizes of the real and integer workspaces used by the
      Diagonal linear solver.

      @cvode CVDiagGetWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space : ('d, 'k) session -> int * int

  (** Returns the number of calls made to the right-hand side
      function due to finite difference Jacobian approximation in the
      Diagonal linear solver.

      @cvode CVDiagGetNumRhsEvals *)
  val get_num_rhs_evals : ('d, 'k) session -> int
end (* }}} *)

(** Direct Linear Solvers operating on dense, banded, and sparse matrices. *)
module Dls : sig (* {{{ *)
  include module type of Sundials_LinearSolver.Direct

  (** Callback functions that compute approximations to a Jacobian
      matrix. In the call [jac arg jm], [arg] is a {!jacobian_arg} with
      three work vectors and the computed Jacobian must be stored in [jm].

      The callback should load the [(i,j)]th entry of [jm] with
      {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of the
      [i]th equation with respect to the [j]th variable, evaluated at the
      values of [t] and [y] obtained from [arg]. Only nonzero elements need
      be loaded into [jm].

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the elements of [arg] nor the matrix [jm] should
               be accessed after the function has returned.}

      @cvode CVLsJacFn *)
  type 'm jac_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg -> 'm -> unit

  (** Function to compute the linear system matrix {% $M = I - \gamma J$ %}
      or an approximation of it. Offers an alternative to evaluating the
      Jacobian of the right-hand-side function.

      In addition to those shared with the Jacobian function, the arguments of
      this function are
      - [m], storage for the computed linear system matrix,
      - [jok], indicates whether the Jacobian-related data needs to be
               updated, and
      - [gamma], the scalar in the formula above.

      The function should return true only if the Jacobian data was
      recomputed.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the Jacobian argument elements nor the matrix [m]
               should be accessed after the function has returned.}

      @cvode CVLsLinSysFn
      @since 5.0.0 *)
  type 'm linsys_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg -> 'm -> bool -> float -> bool

  (** Create a Cvode-specific linear solver from a Jacobian approximation
      function and generic direct linear solver.

      The Jacobian approximation function is optional for dense and banded
      solvers (if not given an internal difference quotient approximation is
      used), but must be provided for other solvers (or [Invalid_argument]
      is raised).

      The [linsys] argument allows to override the standard linear system
      function that calls [jac] to compute {% $M$ %}. This feature is only
      available in Sundials >= 5.0.0.

      @cvode CVodeSetLinearSolver
      @cvode CVodeSetJacFn
      @cvode CVodeSetLinSysFn *)
  val solver :
    ?jac:'m jac_fn ->
    ?linsys:'m linsys_fn ->
    ('m, RealArray.t, 'kind, [>`Dls]) LinearSolver.t ->
    'kind serial_linear_solver

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by a direct
      linear solver.

      @cvode CVodeGetLinWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space : 'kind serial_session -> int * int

  (** Returns the number of calls made by a direct linear solver to the
      Jacobian approximation function.

      @cvode CVodeGetNumJacEvals *)
  val get_num_jac_evals : 'kind serial_session -> int

  (** Returns the number of calls to the right-hand side callback due to
      the finite difference Jacobian approximation.

      @cvode CVodeGetNumLinRhsEvals *)
  val get_num_lin_rhs_evals : 'kind serial_session -> int

end (* }}} *)

(** Scaled Preconditioned Iterative Linear Solvers. *)
module Spils : sig (* {{{ *)
  include module type of Sundials_LinearSolver.Iterative

  (** {3:precond Preconditioners} *)

  (** Arguments passed to the preconditioner solver function.

      @cvode CVLsPrecSolveFn *)
  type 'd prec_solve_arg =
    {
      rhs   : 'd;         (** Right-hand side vector of the linear system. *)
      gamma : float;      (** Scalar $\gamma$ in the Newton
                              matrix given by $M = I - \gamma J$. *)
      delta : float;      (** Input tolerance for iterative methods. *)
      left  : bool;       (** [true] for left preconditioning and
                              [false] for right preconditioning. *)
    }

  (** Callback functions that solve a linear system involving a
      preconditioner matrix. In the call [prec_solve_fn jac arg z],
      [jac] is a {!jacobian_arg} with no work vectors, [arg] is
      a {!prec_solve_arg} that specifies the linear system, and [z] is
      computed to solve {% $P\mathtt{z} = \mathtt{arg.rhs}$%}.
      $P$ is a preconditioner matrix, which approximates, however crudely,
      the Newton matrix {% $M = I - \gamma J$%} where
      {% $J = \frac{\partial f}{\partial y}$%}.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning The elements of [jac], [arg], and [z] should not
               be accessed after the function has returned.}

      @cvode CVLsPrecSolveFn *)
  type 'd prec_solve_fn =
    (unit, 'd) jacobian_arg
    -> 'd prec_solve_arg
    -> 'd
    -> unit

  (** Callback functions that preprocess or evaluate Jacobian-related data
      needed by {!prec_solve_fn}. In the call [prec_setup_fn jac jok gamma],
      [jac] is a {!jacobian_arg} with no work vectors, [jok] indicates
      whether any saved Jacobian-related data can be reused with the current
      value of [gamma], and [gamma] is the scalar $\gamma$ in the Newton
      matrix {% $M = I - \gamma J$%} where $J$ is the Jacobian matrix.
      A function should return [true] if Jacobian-related data was updated
      and [false] if saved data was reused.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning The elements of [jac] should not be accessed after the
               function has returned.}

      @cvode CVodeSetPreconditioner
      @cvode CVLsPrecSetupFn *)
  type 'd prec_setup_fn =
    (unit, 'd) jacobian_arg
    -> bool
    -> float
    -> bool

  (** Specifies a preconditioner, including the type of preconditioning
      (none, left, right, or both) and callback functions.
      The following functions and those in {!Banded} and {!Cvode_bbd}
      construct preconditioners.

      The {!prec_solve_fn} is mandatory. The {!prec_setup_fn} can be
      omitted if not needed.

      @cvode CVodeSetPreconditioner
      @cvode CVLsPrecSetupFn
      @cvode CVLsPrecSolveFn *)
  type ('d,'k) preconditioner = ('d,'k) Cvode_impl.SpilsTypes.preconditioner

  (** No preconditioning.  *)
  val prec_none : ('d, 'k) preconditioner

  (** Left preconditioning. {% $(P^{-1}A)x = P^{-1}b$ %}. *)
  val prec_left :
    ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> ('d, 'k) preconditioner

  (** Right preconditioning. {% $(AP^{-1})Px = b$ %}. *)
  val prec_right :
    ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> ('d, 'k) preconditioner

  (** Left and right preconditioning.
      {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)
  val prec_both :
    ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> ('d, 'k) preconditioner

  (** Banded preconditioners.  *)
  module Banded : sig (* {{{ *)

    (** The range of nonzero entries in a band matrix. *)
    type bandrange =
      { mupper : int; (** The upper half-bandwidth. *)
        mlower : int; (** The lower half-bandwidth. *) }

    (** A band matrix {!preconditioner} based on difference quotients.
        The call [prec_left br] instantiates a left preconditioner which
        generates a banded approximation to the Jacobian with [br.mlower]
        sub-diagonals and [br.mupper] super-diagonals.

        @cvode CVBandPrecInit *)
    val prec_left : bandrange -> (Nvector_serial.data,
                                  [>Nvector_serial.kind]) preconditioner

    (** Like {!prec_left} but preconditions from the right.

        @cvode CVBandPrecInit *)
    val prec_right : bandrange -> (Nvector_serial.data,
                                   [>Nvector_serial.kind]) preconditioner

    (** Like {!prec_left} but preconditions from both sides.

        @cvode CVBandPrecInit *)
    val prec_both : bandrange -> (Nvector_serial.data,
                                  [>Nvector_serial.kind]) preconditioner

    (** {4:stats Banded statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the
        banded preconditioner module.

        @cvode CVBandPrecGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : 'kind serial_session -> int * int

    (** Returns the number of calls to the right-hand side callback for the
        difference banded Jacobian approximation. This counter is only updated
        if the default difference quotient function is used.

        @cvode CVBandPrecGetNumRhsEvals *)
    val get_num_rhs_evals : 'kind serial_session -> int
  end (* }}} *)

  (** {3:lsolvers Solvers} *)

  (** Callback functions that preprocess or evaluate Jacobian-related data
      needed by the jac_times_vec_fn. In the call [jac_times_setup_fn arg],
      [arg] is a {!jacobian_arg} with no work vectors.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning The elements of [arg] should not be accessed after the
               function has returned.}

      @cvode CVLsJacTimesSetupFn *)
  type 'd jac_times_setup_fn = (unit, 'd) jacobian_arg -> unit

  (** Callback functions that compute the Jacobian times a vector. In the
      call [jac_times_vec_fn arg v jv], [arg] is a {!jacobian_arg} with one
      work vector, [v] is the vector multiplying the Jacobian, and [jv] is
      the vector in which to store the
      result—{% $\mathtt{jv} = J\mathtt{v}$%}.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the elements of [arg] nor [v] or [jv] should be
               accessed after the function has returned.}

      @cvode CVLsJacTimesVecFn *)
  type 'd jac_times_vec_fn =
    ('d, 'd) jacobian_arg
    -> 'd (* v *)
    -> 'd (* Jv *)
    -> unit

  (** Create a Cvode-specific linear solver from a generic iterative
      linear solver.

      The [jac_times_rhs] argument specifies an alternative right-hand-side
      function for use in the internal Jacobian-vector product difference
      quotient approximation. It is incorrect to specify both this argument
      and [jac_times_vec].

      NB: a [jac_times_vec] function is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

      NB: a [jac_times_rhs] function is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 5.3.0.

      @cvode CVodeSetLinearSolver
      @cvode CVodeSetJacTimes
      @cvode CVodeSetJacTimesRhsFn *)
  val solver :
    ('m, 'd, 'k, [>`Iter]) LinearSolver.t
    -> ?jac_times_vec:'d jac_times_setup_fn option * 'd jac_times_vec_fn
    -> ?jac_times_rhs:'d rhsfn
    -> ('d, 'k) preconditioner
    -> ('d, 'k) linear_solver

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by the spils
      linear solver.

      @cvode CVodeGetLinWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space       : ('d, 'k) session -> int * int

  (** Returns the cumulative number of linear iterations.

      @cvode CVodeGetNumLinIters *)
  val get_num_lin_iters    : ('d, 'k) session -> int

  (** Returns the cumulative number of linear convergence failures.

      @cvode CVodeGetNumLinConvFails *)
  val get_num_lin_conv_fails: ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the setup function with
      [jok=false].

      @cvode CVodeGetNumPrecEvals *)
  val get_num_prec_evals   : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the preconditioner solve
      function.

      @cvode CVodeGetNumPrecSolves *)
  val get_num_prec_solves  : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the Jacobian-vector
      setup function.

      @cvode CVodeGetNumJTSetupEvals
      @since 3.0.0 *)
  val get_num_jtsetup_evals : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the Jacobian-vector
      function.

      @cvode CVodeGetNumJtimesEvals *)
  val get_num_jtimes_evals : ('d, 'k) session -> int

  (** Returns the number of calls to the right-hand side callback for
      finite difference Jacobian-vector product approximation. This counter is
      only updated if the default difference quotient function is used.

      @cvode CVodeGetNumLinRhsEvals *)
  val get_num_lin_rhs_evals : ('d, 'k) session -> int

  (** {3:lowlevel Low-level solver manipulation}

      The {!init} and {!reinit} functions are the preferred way to set or
      change preconditioner functions. These low-level functions are provided
      for experts who want to avoid resetting internal counters and other
      associated side-effects. *)

  (** Change the preconditioner functions.

      @cvode CVodeSetPreconditioner
      @cvode CVLsPrecSolveFn
      @cvode CVLsPrecSetupFn *)
  val set_preconditioner :
    ('d, 'k) session
    -> ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> unit

  (** Change the Jacobian-times-vector function.

      NB: the [jac_times_setup] argument is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

      @cvode CVodeSetJacTimes
      @cvode CVLsJacTimesSetupFn
      @cvode CVLsJacTimesVecFn *)
  val set_jac_times :
    ('d, 'k) session
    -> ?jac_times_setup:'d jac_times_setup_fn
    -> 'd jac_times_vec_fn
    -> unit

  (** Remove a Jacobian-times-vector function and use the default
      implementation.

      @cvode CVodeSetJacTimes
      @cvode CVodeJacTimesSetupFn
      @cvode CVodeJacTimesVecFn *)
  val clear_jac_times : ('d, 'k) session -> unit

end (* }}} *)

(** Create a CVode-specific linear solver from a generic matrix embedded
    solver.

    @cvode CVodeSetLinearSolver
    @since 5.8.0 *)
val matrix_embedded_solver :
       (unit, 'data, 'kind, [>`MatE]) LinearSolver.t
    -> ('data, 'kind) linear_solver

(** {2:tols Tolerances} *)

(** Functions that set the multiplicative error weights for use in the weighted
    RMS norm. The call [efun y ewt] takes the dependent variable vector [y] and
    fills the error-weight vector [ewt] with positive values or raises
    {!Sundials.NonPositiveEwt}. Other exceptions are eventually propagated, but
    should be avoided ([efun] is not allowed to abort the solver). *)
type 'data error_weight_fun = 'data -> 'data -> unit

(** Tolerance specifications. *)
type ('data, 'kind) tolerance =
  | SStolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVtolerances of float * ('data, 'kind) Nvector.t
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFtolerances of 'data error_weight_fun
    (** Set the multiplicative error weights for the weighted RMS norm. *)

(** A default relative tolerance of 1.0e-4 and absolute tolerance of 1.0e-8. *)
val default_tolerances : ('data, 'kind) tolerance

(** {2:solver Solver initialization and use} *)

(** Choice of linear multistep method.

    @cvode CVodeCreate *)
type lmm =
  | Adams   (** Adams-Moulton formulas (non-stiff systems). *)
  | BDF     (** Backward Differentiation Formulas (stiff systems). *)

(** Called by the solver to calculate the values of root functions. These
    ‘zero-crossings’ are used to detect significant events. The function is
    passed three arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
    - [gout], a vector for storing the value of $g(t, y)$.

    {warning [y] and [gout] should not be accessed after the function has
             returned.}

    @cvode cvRootFn *)
type 'd rootsfn = float -> 'd -> RealArray.t -> unit

(** A function to compute the projection of the solution and, if enabled, the
    error on the constraint manifold.

    Such functions take the following arguments:
    - [t], the independent variable,
    - [ycur], the dependent variable vector,
    - [corr], the correction to the dependent variable vector so that
              {% $y(t) + c$ %} satisifies the constraint equation,
    - [eps], the tolerance to use in the nonlinear stopping test when
             solving the nonlinear contrained least-squares problem, and
    - [err], if error project is enabled (the default), the current error
             estimate on input and updated to the projected error for output.

    (The [err] argument is [None] if error projection is not enabled.)

    Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
    Any other exception is treated as an unrecoverable error. The integrator
    will, in most cases, try to correct and reattempt the step.

    The solve should stop when the WRMS norm of the current iterate update is
    less than [eps]. The projection routine can access the error weight vector
    with {!get_err_weights}.

    @cvode CVProjFn
    @since 5.3.0 *)
type 'd proj_fn = float -> 'd -> 'd -> float -> 'd option -> unit

(** Creates and initializes a session with the solver. The call
    {[init lmm tol ~nlsolver ~nlsrhsfn ~lsolver f ~roots:(nroots, g) ~projfn t0 y0]}
    has as arguments:
    - [lmm],      the linear multistep method (see {!lmm}),
    - [tol],      the integration tolerances,
    - [nlsolver], the solver to use to calculate integration steps,
    - [nlsrhsfn], alternative right-hand-side function to use in
                  nonlinear system function evaluations,
    - [lsolver],  used by [nlsolver]s based on Newton interation,
    - [f],        the ODE right-hand-side function,
    - [nroots],   the number of root functions,
    - [g],        the root function ([(nroots, g)] defaults to {!no_roots}),
    - [projfn],   enables projection onto the constraint manifold using the
                  given function after each time step,
    - [t0],       the initial value of the independent variable, and,
    - [y0],       a vector of initial values that also determines the number
                  of equations.

    This function does everything necessary to initialize a session, i.e.,
    it makes the calls referenced below. The {!solve_normal} and
    {!solve_one_step} functions may be called directly.

    If an [nlsolver] is not specified, then the
    {{!Sundials_NonlinearSolver.Newton}Newton} module is used by default.
    In this case only, [lsolver] defaults to {!Diag.solver} if not otherwise
    specified. Specifying an [nlsolver] that requires a linear solver without
    specifying an [lsolver] results in a {!NonlinearInitFailure} (or
    {!IllInput} for Sundials < 4.0.0) exception on the first call to
    {!solve_normal} or {!solve_one_step}.

    The projection feature is only supported for Sundials >= 5.3.0 and the
    {{!lmm}BDF} method.

    The alternative right-hand-side function for nonlinear system function
    evaluations is only supported for Sundials >= 5.8.0.

    By default, the session is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing an
    optional [context] argument.

    @cvode CVodeCreate
    @cvode CVodeInit
    @cvode CVodeRootInit
    @cvode CVodeSetLinearSolver
    @cvode CVodeSetNonlinearSolver
    @cvode CVodeSStolerances
    @cvode CVodeSVtolerances
    @cvode CVodeWFtolerances
    @cvode CVEwtFn
    @cvode CVodeSetProjFn
    @cvode CVodeSetNlsRhsFn *)
val init :
       ?context:Context.t
    -> lmm
    -> ('data, 'kind) tolerance
    -> ?nlsolver
         : ('data, 'kind, ('data, 'kind) session, [`Nvec])
             Sundials_NonlinearSolver.t
    -> ?nlsrhsfn:'data rhsfn
    -> ?lsolver  : ('data, 'kind) linear_solver
    -> 'data rhsfn
    -> ?roots:(int * 'data rootsfn)
    -> ?projfn:'data proj_fn
    -> float
    -> ('data, 'kind) Nvector.t
    -> ('data, 'kind) session

(** A convenience value for signalling that there are no roots to monitor. *)
val no_roots : (int * 'd rootsfn)

(** Values returned by the step functions. Failures are indicated by
    exceptions.

    @cvode CVode *)
type solver_result =
  | Success             (** The solution was advanced. {cconst CV_SUCCESS} *)
  | RootsFound          (** A root was found. See {!get_root_info}.
                            {cconst CV_ROOT_RETURN} *)
  | StopTimeReached     (** The stop time was reached. See {!set_stop_time}.
                            {cconst CV_TSTOP_RETURN} *)

(** Integrates an ODE system over an interval. The call
    [tret, r = solve_normal s tout yout] has as arguments
    - [s], a solver session,
    - [tout], the next time at which a solution is desired, and,
    - [yout], a vector to store the computed solution.

    It returns [tret], the time reached by the solver, which will be equal to
    [tout] if no errors occur, and, [r], a {!solver_result}.

    @cvode CVode (CV_NORMAL)
    @raise IllInput Missing or illegal solver inputs.
    @raise TooClose The initial and final times are too close to each other and not initial step size was specified.

    @raise TooMuchWork The requested time could not be reached in [mxstep] internal steps.
    @raise TooMuchAccuracy The requested accuracy could not be satisfied.
    @raise ErrFailure Too many error test failures within a step or at the minimum step size.
    @raise ConvergenceFailure Too many convergence test failures within a step or at the minimum step size.
    @raise LinearInitFailure Linear solver initialization failed.
    @raise LinearSetupFailure Linear solver setup failed unrecoverably.
    @raise LinearSolveFailure Linear solver solution failed unrecoverably.
    @raise RhsFuncFailure Unrecoverable failure in the RHS function [f].
    @raise FirstRhsFuncFailure Initial unrecoverable failure in the RHS function [f].
    @raise RepeatedRhsFuncFailure Too many convergence test failures, or unable to estimate the initial step size, due to repeated recoverable errors in the right-hand side function.
    @raise UnrecoverableRhsFuncFailure The right-hand side function had a recoverable error, but no recovery was possible. This error can only occur after an error test failure at order one.
    @raise RootFuncFailure Failure in the rootfinding function [g]. *)
val solve_normal : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                        -> float * solver_result

(** Like {!solve_normal} but returns after one internal solver step.

    @cvode CVode (CV_ONE_STEP) *)
val solve_one_step : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
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

    @cvode CVodeGetDky
    @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
    @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
val get_dky : ('d, 'k) session -> ('d, 'k) Nvector.t -> float -> int -> unit

(** Reinitializes the solver with new parameters and state values. The
    values of the independent variable, i.e., the simulation time, and the
    state variables must be given. If given, [nlsolver] specifies a nonlinear
    solver, [nlsrhsfn] specifies an alternative rhs function for nonlinear
    system function evaluations,[lsolver] specifies a linear solver,
    [roots] specifies a new root finding function, and [rhsfn] specifies
    a new rhs function; all default to unchanged.

    If the new problem does not have a constraint equation, but the old one
    did, then {!set_proj_frequency} must a zero argument to disable
    projection.

    @cvode CVodeReInit
    @cvode CVodeSetLinearSolver
    @cvode CVodeSetNonlinearSolver
    @cvode CVodeSetNlsRhsFn *)
val reinit :
  ('d, 'k) session
  -> ?nlsolver:('d, 'k, ('d, 'k) session, [`Nvec]) Sundials_NonlinearSolver.t
  -> ?nlsrhsfn:'d rhsfn
  -> ?lsolver:('d, 'k) linear_solver
  -> ?roots:(int * 'd rootsfn)
  -> ?rhsfn:'d rhsfn
  -> float
  -> ('d, 'k) Nvector.t
  -> unit

(** {2:set Modifying the solver} *)

(** Sets the integration tolerances.

    @cvode CVodeSStolerances
    @cvode CVodeSVtolerances
    @cvode CVodeWFtolerances
    @cvode CVEwtFn *)
val set_tolerances : ('d, 'k) session -> ('d, 'k) tolerance -> unit

(** Configure the default error handler to write messages to a file.
    By default it writes to Logfile.stderr.

    @cvode CVodeSetErrFile *)
val set_error_file : ('d, 'k) session -> Logfile.t -> unit

(** Specifies a custom function for handling error messages.
    The handler must not fail: any exceptions are trapped and discarded.

    @cvode CVodeSetErrHandlerFn
    @cvode CVErrHandlerFn *)
val set_err_handler_fn : ('d, 'k) session -> (Util.error_details -> unit) -> unit

(** Restores the default error handling function.

    @cvode CVodeSetErrHandlerFn *)
val clear_err_handler_fn : ('d, 'k) session -> unit

(** Specifies a function to be called after the given number of successful
    steps.

    The solver solution may be read by the monitoring function, but it should
    not be changed.

    This function requires that the underlying library was explicitly built
    with support for monitoring (see {!Sundials_Config.monitoring_enabled}).

    @cvode CVodeSetMonitorFn
    @cvode CVodeSetMonitorFrequency
    @raise NotImplementedBySundialsVersion if not provided by the underlying library
    @since 5.3.0 *)
val set_monitor_fn
  : ('d, 'k) session -> int -> (('d, 'k) session -> unit) -> unit

(** Sets the number of successful steps between calls to the monitoring
    function.

    @cvode CVodeSetMonitorFrequency
    @raise NotImplementedBySundialsVersion if not provided by the underlying library
    @since 5.3.0 *)
val set_monitor_frequency : ('d, 'k) session -> int -> unit

(** Turns monitoring off.

    @cvode CVodeSetMonitorFn
    @since 5.3.0 *)
val clear_monitor_fn : ('d, 'k) session -> unit

(** {3:set_main Main solver optional input functions} *)

(** Specifies the maximum order of the linear multistep method.

    @cvode CVodeSetMaxOrd *)
val set_max_ord : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of steps taken in attempting to reach
    a given output time.

    @cvode CVodeSetMaxNumSteps *)
val set_max_num_steps : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of messages warning that [t + h = t] on
    the next internal step.

    @cvode CVodeSetMaxHnilWarns *)
val set_max_hnil_warns : ('d, 'k) session -> int -> unit

(** Indicates whether the BDF stability limit detection algorithm should be
    used.

    @cvode CVodeSetStabLimDet
    @cvode <Mathematics_link.html#cvode-mathematics-stablimit> BDF Stability Limit Detection *)
val set_stab_lim_det : ('d, 'k) session -> bool -> unit

(** Specifies the initial step size.

    @cvode CVodeSetInitStep *)
val set_init_step : ('d, 'k) session -> float -> unit

(** Specifies a lower bound on the magnitude of the step size.

    @cvode CVodeSetMinStep *)
val set_min_step : ('d, 'k) session -> float -> unit

(** Specifies an upper bound on the magnitude of the step size.

    @cvode CVodeSetMaxStep *)
val set_max_step : ('d, 'k) session -> float -> unit

(** Limits the value of the independent variable [t] when solving.
    By default no stop time is imposed.

    @cvode CVodeSetStopTime *)
val set_stop_time : ('d, 'k) session -> float -> unit

(** Specifies the maximum number of error test failures permitted in attempting
    one step.

    @cvode CVodeSetMaxErrTestFails *)
val set_max_err_test_fails : ('d, 'k) session -> int -> unit

(** Specifies a vector defining inequality constraints for each
    component of the solution vector [y].  See {!Sundials.Constraint}.

    @cvode CVodeSetConstraints *)
val set_constraints : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Disables constraint checking.

    @cvode CVodeSetConstraints *)
val clear_constraints : ('d, 'k) session -> unit

(** {3:set_ls Linear solver interface optional input functions} *)

(** Specifies the frequency of calls to the linear solver setup routine.
    Positive values specify the number of time steps between setup calls,
    negative values force recomputation at each Newton step, and zero values
    reset to the default (20).

    @cvode CVodeSetLSetupFrequency
    @since 5.4.0 *)
val set_lsetup_frequency : ('d, 'k) session -> int -> unit

(** Sets the maximum number of time steps to wait before recomputation of
    the Jacobian or recommendation to update the preconditioner.
    If the integer argument is less than or equal to 0, a default value of
    50 is used.

    @cvode CVodeSetJacEvalFrequency
    @since 4.0.0 *)
val set_jac_eval_frequency : ('d, 'k) session -> int -> unit

(** Enables or disables scaling of the linear system solution to account
    for a change in {% $\gamma$ %} in the linear system.
    Linear solution scaling is enabled by default when a matrix-based
    linear solver is attached.

    @since 5.2.0
    @cvode CVodeSetLinearSolutionScaling *)
val set_linear_solution_scaling : ('d, 'k) session -> bool -> unit

(** Sets the factor by which the Krylov linear solver's convergence test
    constant is reduced from the Newton iteration test constant.
    This factor must be >= 0; passing 0 specifies the default (0.05).

    @cvode CVodeSetEpsLin *)
val set_eps_lin : ('d, 'k) session -> float -> unit

(** Sets the factor for converting from the integrator tolerance (WRMS
    norm) to the linear solver tolerance (L2 norm). That is,
    {% $\mathit{tol}_{\mathsf{L2}} =
        \mathit{fact}\cdot\mathit{tol}_{\mathsf{WRMS}}$ %}.
    The given value is used directly if it is greater than zero.
    If it is zero (the default), then the square root of the state
    vector length is used.
    If it is less than zero, then the square root of the dot product of a
    state vector full of ones with itself is used.

    @cvode CVodeSetLSNormFactor
    @since 5.4.0 *)
val set_ls_norm_factor : ('d, 'k) session -> float -> unit

(** {3:set_nls Nonlinear solver interface optional input functions} *)

(** Specifies the maximum number of nonlinear solver iterations permitted per
    step.

    @cvode CVodeSetMaxNonlinIters *)
val set_max_nonlin_iters : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver convergence failures
    permitted during one step.

    @cvode CVodeSetMaxConvFails *)
val set_max_conv_fails : ('d, 'k) session -> int -> unit

(** Specifies the safety factor used in the nonlinear convergence test.

    @cvode CVodeSetNonlinConvCoef *)
val set_nonlin_conv_coef : ('d, 'k) session -> float -> unit

(** {3:set_tsa Time step adaptivity optional input functions} *)

(** Specifies the lower {% $\eta_\mathrm{min\_fx}$ %} and
    upper {% $\eta_\mathrm{max\_fx}$ %} bounds in which the step size will
    remain unchanged. The default values are {% $\eta_\mathrm{min\_fx} = 0$ %}
    and {% $\eta_\mathrm{max\_fx} = 1.5$ %}.
    If the step size change factor satisfies,
    {% $\eta_\mathrm{min\_fx} < \eta < \eta_\mathrm{max\_fx}$ %},
    the current step size is retained.

    If [min < 0] or [min >= 1], the default minimum value is used.
    If [max < 1], the default maximum value is used.

    @cvode CVodeSetEtaFixedStepBounds
    @since 6.2.0 *)
val set_eta_fixed_step_bounds
      : ('d, 'k) session -> min:float -> max:float -> unit -> unit

(** Specifies the maximum step size factor after the first time step
    {% $\eta_\mathrm{max\_fs}$ %}.
    If the argument is [<= 1], the default value of {% $10^4$ %} is used.

    @cvode CVodeSetEtaMaxFirstStep
    @since 6.2.0 *)
val set_eta_max_first_step           : ('d, 'k) session -> float -> unit

(** Specifies the maximum step size factor for steps early in the integration,
    {% $\eta_\mathrm{max\_es}$ %}.
    If the argument is [<= 1], the default value of 10 is used.

    @cvode CVodeSetEtaMaxEarlyStep
    @since 6.2.0 *)
val set_eta_max_early_step           : ('d, 'k) session -> float -> unit

(** Specifies the number of steps to use the early integration maximum step
    size factor. If the argument is [< 0], the default value of 10 is used.

    @cvode CVodeSetNumStepsEtaMaxEarlyStep
    @since 6.2.0 *)
val set_num_steps_eta_max_early_step : ('d, 'k) session -> int -> unit

(** Specifies the minimum step size factor, {% $\eta_\mathrm{min\_gs}$ %}.
    If the argument is [<= 0] or [>= 1], the default value of 1.0 is used.

    @cvode CVodeSetEtaMin
    @since 6.2.0 *)
val set_eta_min                      : ('d, 'k) session -> float -> unit

(** Specifies the maximum step size factor, {% $\eta_\mathrm{max\_gs}$ %}.
    If the argument is [<= 1], the default value of 10 is used.

    @cvode CVodeSetEtaMax
    @since 6.2.0 *)
val set_eta_max                      : ('d, 'k) session -> float -> unit

(** Specifies the minimum step size factor after an error test failure,
    {% $\eta_\mathrm{min\_ef}$ %}.
    If the argument is [<= 0] or [>= 1], the default value of 0.1 is used.

    @cvode CVodeSetEtaMinErrFail
    @since 6.2.0 *)
val set_eta_min_err_fail             : ('d, 'k) session -> float -> unit

(** Specifies the maximum step size factor after multiple error test failures,
    {% $\eta_\mathrm{max\_ef}$ %}.
    If the argument is [<= 0] or [>= 1], the default value of 0.2 is used.

    @cvode CVodeSetEtaMaxErrFail
    @since 6.2.0 *)
val set_eta_max_err_fail             : ('d, 'k) session -> float -> unit

(** Specifies the number of error test failures necessary to enforce the
    maximum step size factor.
    If the argument is [< 0], the default value of 2 is used.
    If the argument is [0], the value set by {!set_eta_max} is used.

    @cvode CVodeSetNumFailsEtaMaxErrFail
    @since 6.2.0 *)
val set_num_fails_eta_max_err_fail   : ('d, 'k) session -> int -> unit

(** Specifies the step size factor after a nonlinear solver failure,
    {% $\eta_\mathrm{cf}$ %}.
    If the argument is [<= 0] or [>= 1], the default value of 0.25 is used.

    @cvode CVodeSetEtaConvFail
    @since 6.2.0 *)
val set_eta_conv_fail                : ('d, 'k) session -> float -> unit

(** {3:set_proj Projection optional input functions} *)

(** Enables or disables projection of the error estimate by the projection
    function.

    @cvode CVodeSetProjErrEst
    @since 5.3.0 *)
val set_proj_err_est : ('d, 'k) session -> bool -> unit

(** Set the frequency with which the projection is performed. The default is
    1, that is, every time step. A value of 0 disables projection and a value
    less than zero restores the default.

    @cvode CVodeSetProjFrequency
    @since 5.3.0 *)
val set_proj_frequency : ('d, 'k) session -> int -> unit

(** Set the maximum number of projection failures in a step attempt before an
    unrecoverable error is returned. The default is 10. A value less than 1
    restores the default.

    @cvode CVodeSetMaxNumProjFails
    @since 5.3.0 *)
val set_max_num_proj_fails : ('d, 'k) session -> int -> unit

(** Set the tolerance for the nonlinear-constrained least-squares problem
    solved by the projection function. The default is 0.1. A value less than
    or equal to zero restores the default.

    @cvode CVodeSetEpsProj
    @since 5.3.0 *)
val set_eps_proj : ('d, 'k) session -> float -> unit

(** Sets the time-step reduction factor to apply on a projection function
    failure. The default is 0.25. A value less than or equal to 1, or greater
    than 1 restores the default.

    @cvode CVodeSetProjFailEta
    @since 5.3.0 *)
val set_proj_fail_eta : ('d, 'k) session -> float -> unit

(** {2:get Querying the solver (optional output functions)} *)

(** Returns the real and integer workspace sizes.

    @cvode CVodeGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space          : ('d, 'k) session -> int * int

(** Returns the cumulative number of internal steps taken by the solver.

    @cvode CVodeGetNumSteps *)
val get_num_steps           : ('d, 'k) session -> int

(** Returns the number of calls to the right-hand side function.

    @cvode CVodeGetNumRhsEvals *)
val get_num_rhs_evals       : ('d, 'k) session -> int

(** Returns the number of calls made to the linear solver's setup function.

    @cvode CVodeGetNumLinSolvSetups *)
val get_num_lin_solv_setups : ('d, 'k) session -> int

(** Returns the number of local error test failures that have occurred.

    @cvode CVodeGetNumErrTestFails *)
val get_num_err_test_fails  : ('d, 'k) session -> int

(** Returns the number of failed steps due to a nonlinear solver failure.

    @cvode CVodeGetNumStepSolveFails
    @since 6.2.0 *)
val get_num_step_solve_fails : ('d, 'k) session -> int

(** Returns the integration method order used during the last internal step.

    @cvode CVodeGetLastOrder *)
val get_last_order          : ('d, 'k) session -> int

(** Returns the integration method order to be used on the next internal step.

    @cvode CVodeGetCurrentOrder *)
val get_current_order       : ('d, 'k) session -> int

(** Returns the current state vector. This vector provides direct access to
    the data within the integrator.

    @cvode CVodeGetCurrentState
    @since 5.0.0 *)
val get_current_state : ('d, 'k) session -> 'd

(** Internal data required to construct the current nonlinear implicit
    system within a nonlinear solver. *)
type 'd nonlin_system_data = {
  tn    : float;
    (** Independent variable value {% $t_n$ %}. *)
  ypred : 'd;
    (** Predicted state vector {% $y_{\mathit{pred}}$ %} at {% $t_n$ %}. This
        data must not be changed. *)
  yn    : 'd;
    (** State vector {% $y^n$ %}. This data may not be current and may
        need to be filled. *)
  fn    : 'd;
    (** The right-hand side function evaluated at the current time and state,
        {% $f(t_n, y^n)$ %}. * This data may not be current and may need to
        be filled. *)
  gamma : float;
    (** Current value of {% $\gamma$ %}. *)
  rl1   : float;
      (** A scaling factor used to compute {% $\tilde{a}_n =
          \mathtt{rl1}\cdot\mathtt{zn1}$ %}. *)
  zn1   : 'd;
      (** A vector used to compute {% $\tilde{a}_n =
          \mathtt{rl1}\cdot\mathtt{zn1}$ %}. *)
}

(** Gives direct access to the internal data required to construct the
    current nonlinear system within a nonlinear solver. This
    function should be called inside the nonlinear system function.
    If the nonlinear solver uses the [lsetup] or [lsolve] functions, then
    the nonlinear solver system function must fill the [zi] and [fi]
    vectors with, respectively, the current state and corresponding
    evaluation of the right-hand-side function:
    {% $y^n = y_{\mathit{pred}} + y_{\mathit{cor}}$ %} and
    {% $f_n = f(t_n, y^n)$ %} where {% $y_{\mathit{cor}}$ %} is the
    first argument of the nonlinear solver system function. Within a custom
    linear solver, then the vectors [yn] and [fn] are only current after
    an evaluation of the nonlinear system function.

    @cvode CVodeGetNonlinearSystemData
    @since 5.4.0 *)
val get_nonlin_system_data : ('d, 'k) session -> 'd nonlin_system_data

(** Computes the current stage state vector using the stored prediction and
    the supplied correction from the nonlinear solver. The call
    [compute_state s ycor yn] computes {% $y^n = y_{\mathit{pred}}
    + y_{\mathit{cor}}$ %}.

    @cvode CVodeComputeState
    @since 5.4.0 *)
val compute_state : ('d, 'k) session
                    -> ('d, 'k) Nvector.t
                    -> ('d, 'k) Nvector.t
                    -> unit

(** Returns the current value of {% $\gamma$ %}.
    This scalar appears in the internal Newton equation,
    {% $M = I - \gamma J$ %}.

    @cvode CVodeGetCurrentGamma
    @since 5.0.0 *)
external get_current_gamma : ('d, 'k) session -> (float [@unboxed])
    = "sunml_cvode_get_current_gamma"
      "sunml_cvode_get_current_gamma_unboxed"

(** Returns the integration step size taken on the last internal step.

    @cvode CVodeGetLastStep *)
val get_last_step           : ('d, 'k) session -> float

(** Returns the integration step size to be attempted on the next internal step.

    @cvode CVodeGetCurrentStep *)
val get_current_step        : ('d, 'k) session -> float

(** Returns the the value of the integration step size used on the first step.

    @cvode CVodeGetActualInitStep *)
val get_actual_init_step    : ('d, 'k) session -> float

(** Returns the the current internal time reached by the solver.

    @cvode CVodeGetCurrentTime *)
val get_current_time        : ('d, 'k) session -> float

(** Returns the number of order reductions dictated by the BDF stability limit
    detection algorithm.

    @cvode CVodeGetNumStabLimOrderReds
    @cvode <Mathematics_link.html#cvode-mathematics-stablimit> BDF Stability Limit Detection *)
val get_num_stab_lim_order_reds : ('d, 'k) session -> int

(** Returns a suggested factor by which the user's tolerances should be scaled
    when too much accuracy has been requested for some internal step.

    @cvode CVodeGetTolScaleFactor *)
val get_tol_scale_factor : ('d, 'k) session -> float

(** Returns the solution error weights at the current time.

    @cvode CVodeGetErrWeights *)
val get_err_weights : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Returns the vector of estimated local errors.

    @cvode CVodeGetEstLocalErrors *)
val get_est_local_errors : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Summaries of integrator statistics. *)
type integrator_stats = {
    num_steps : int;
      (** Cumulative number of internal solver steps. *)
    num_rhs_evals : int;
      (** Number of calls to the right-hand side function. *)
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

    @cvode CVodeGetIntegratorStats *)
val get_integrator_stats    : ('d, 'k) session -> integrator_stats

(** Prints the integrator statistics on the given channel.

    @cvode CVodeGetIntegratorStats *)
val print_integrator_stats  : ('d, 'k) session -> out_channel -> unit

(** Summaries of linear solver statistics. *)
type linear_solver_stats = {
    jac_evals : int;
      (** Number of calls made by a linear solver to the Jacobian
          approximation function. *)
    lin_rhs_evals : int;
      (** Number of calls to the right-hand side callback due to
          the finite difference Jacobian approximation. *)
    lin_iters : int;
      (** The cumulative number of linear iterations. *)
    lin_conv_fails : int;
      (** The cumulative number of linear convergence failures. *)
    prec_evals : int;
      (** The cumulative number of calls to the setup function. *)
    prec_solves : int;
      (** The cumulative number of calls to the solve function. *)
    jtsetup_evals : int;
      (** The cumulative number of calls to the Jacobian-vector
          setup function. *)
    jtimes_evals : int;
      (** The cumulative number of calls to the Jacobian-vector function. *)
  }

(** Returns linear solver statistics as a group.

    @cvode CVodeGetLinSolveStats
    @since 5.3.0 *)
val get_linear_solver_stats : ('d, 'k) session -> linear_solver_stats

(** Returns the cumulative number of nonlinear (functional or Newton)
    iterations.

    @cvode CVodeGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('d, 'k) session -> int

(** Returns the cumulative number of nonlinear convergence failures.

    @cvode CVodeGetNumNonlinSolvConvFails *)
val get_num_nonlin_solv_conv_fails : ('d, 'k) session -> int

(** Returns both the numbers of nonlinear iterations performed [nniters] and
    nonlinear convergence failures [nncfails].

    @cvode CVodeGetNonlinSolvStats
    @return ([nniters], [nncfails]) *)
val get_nonlin_solv_stats : ('d, 'k) session -> int * int

(** Outputs all of the integrator, nonlinear solver, linear solver, and other
    statistics.

    @cvode CVodePrintAllStats
    @since 6.2.0 *)
val print_all_stats
      : ('d, 'k) session -> Logfile.t -> Sundials.output_format -> unit

(** {2:roots Additional root-finding functions} *)

(** [set_root_direction s dir] specifies the direction of zero-crossings to be
    located and returned. [dir] may contain one entry for each root function.

    @cvode CVodeSetRootDirection *)
val set_root_direction : ('d, 'k) session -> RootDirs.d array -> unit

(** Like {!set_root_direction} but specifies a single direction for all root
    functions.

    @cvode CVodeSetRootDirection *)
val set_all_root_directions : ('d, 'k) session -> RootDirs.d -> unit

(** Disables issuing a warning if some root function appears to be identically
    zero at the beginning of the integration.

    @cvode CVodeSetNoInactiveRootWarn *)
val set_no_inactive_root_warn : ('d, 'k) session -> unit

(** Returns the number of root functions. *)
val get_num_roots : ('d, 'k) session -> int

(** Fills an array showing which functions were found to have a root.

    @cvode CVodeGetRootInfo *)
val get_root_info : ('d, 'k) session -> Roots.t -> unit

(** Returns the cumulative number of calls made to the user-supplied root
    function g.

    @cvode CVodeGetNumGEvals *)
val get_num_g_evals : ('d, 'k) session -> int

(** Returns the current total number of projection evaluations.

    @cvode CVodeGetNumProjEvals *)
val get_num_proj_evals : ('d, 'k) session -> int

(** Returns the current total number of projection evaluation failures.
    @cvode CVodeGetNumProjFails *)
val get_num_proj_fails : ('d, 'k) session -> int

(** {2:exceptions Exceptions} *)

(** Raised on missing or illegal solver inputs. Also raised if an element
    of the error weight vector becomes zero during time stepping, or the
    linear solver initialization function failed, or a root was found both at
    [t] and very near [t].

    @cvode <Constants_link.html> CV_ILL_INPUT *)
exception IllInput

(** The initial and final times are too close to each other and an initial step
    size was not specified.

    @cvode <Constants_link.html> CV_TOO_CLOSE *)
exception TooClose

(** The requested time could not be reached in [mxstep] internal steps.
    See {!set_max_num_steps}

    @cvode <Constants_link.html> CV_TOO_MUCH_WORK *)
exception TooMuchWork

(** The requested accuracy could not be satisfied.

    @cvode <Constants_link.html> CV_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** Too many error test failures within a step or at the minimum step size.
    See {!set_max_err_test_fails} and {!set_min_step}.

    @cvode <Constants_link.html> CV_ERR_FAILURE *)
exception ErrFailure

(** Too many convergence test failures within a step or at the minimum step
    size. See {!set_max_conv_fails} and {!set_min_step}.

    @cvode <Constants_link.html> CV_CONV_FAILURE *)
exception ConvergenceFailure

(** Linear solver initialization failed.

    @cvode <Constants_link.html> CV_LINIT_FAIL *)
exception LinearInitFailure

(** Linear solver setup failed in an unrecoverable manner.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.PSetFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @cvode CVodeGetLastLinFlag
    @cvode <Constants_link.html> CV_LSETUP_FAIL *)
exception LinearSetupFailure of exn option

(** Linear solver solution failed in an unrecoverable manner.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.ATimesFailure},
    {!Sundials_LinearSolver.PSolveFailure},
    {!Sundials_LinearSolver.GSFailure},
    {!Sundials_LinearSolver.QRSolFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @cvode CVodeGetLastLinFlag
    @cvode <Constants_link.html> CV_LSOLVE_FAIL *)
exception LinearSolveFailure of exn option

(** The nonlinear solver failed in a general way.

    @cvode <Constants_link.html> CV_NLS_FAIL
    @since 5.0.0 *)
exception NonlinearSolverFailure

(** Nonlinear solver initialization failed.

    @cvode <Constants_link.html> CV_NLS_INIT_FAIL *)
exception NonlinearInitFailure

(** Nonlinear solver setup failed in an unrecoverable manner.

    @cvode <Constants_link.html> CV_NLS_SETUP_FAIL *)
exception NonlinearSetupFailure

(** The right-hand side function failed in an unrecoverable manner.

    @cvode <Constants_link.html> CV_RHSFUNC_FAIL *)
exception RhsFuncFailure

(** The right-hand side function had a recoverable error when first called.

    @cvode <Constants_link.html> CV_FIRST_RHSFUNC_ERR *)
exception FirstRhsFuncFailure

(** Too many convergence test failures, or unable to estimate the initial step
    size, due to repeated recoverable errors in the right-hand side function.

    @cvode <Constants_link.html> CV_REPTD_RHSFUNC_ERR *)
exception RepeatedRhsFuncFailure

(** The right-hand side function had a recoverable error, but no recovery was
    possible. This error can only occur after an error test failure at order
    one.

    @cvode <Constants_link.html> CV_UNREC_RHSFUNC_ERR *)
exception UnrecoverableRhsFuncFailure

(** The rootfinding function failed.

    @cvode <Constants_link.html> CV_RTFUNC_FAIL *)
exception RootFuncFailure

(** No solution satisfying the inequality constraints could be found.

    @cvode <Constants_link.html> CV_CONSTR_FAIL *)
exception ConstraintFailure

(** Raised by {!get_dky} for invalid order values.

    @cvode <Constants_link.html> CVodeGetDky (CV_BAD_K) *)
exception BadK

(** Raised by {!get_dky} for invalid time values.

    @cvode <Constants_link.html> CVodeGetDky (CV_BAD_T) *)
exception BadT

(** A fused vector operation failed.

    @cvode <Constants_link.html> CV_VECTOROP_ERR *)
exception VectorOpErr

(** The projection function failed.

    @cvode <Constants_link.html> CV_PROJFUNC_FAIL
    @since 5.0.0 *)
exception ProjFuncFailure

(** The projection function failed repeatedly.

    @cvode <Constants_link.html> CV_REPTD_PROJFUNC_ERR
    @since 5.0.0 *)
exception RepeatedProjFuncError

(** The project functionality is not enabled. A projection function must be
    given in the call to {!init} and the last call to {!set_proj_frequency}
    must not have set the frequency to zero.

    @cvode <Constants_link.html> CV_PROJ_MEM_NULL
    @since 5.0.0 *)
exception ProjectionNotEnabled

