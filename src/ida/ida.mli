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
      {- {{:#nls}Advanced nonlinear solver functions}}
      {- {{:#exceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** A session with the IDA solver.

    An example session with Ida ({openfile ida_skel.ml}): {[
#include "../../examples/ocaml/skeletons/ida_skel.ml"
    ]}

    @ida <Usage/index.html#a-skeleton-of-the-users-main-program> Skeleton of main program
 *)
type ('d, 'k) session = ('d, 'k) Ida_impl.session

(** Alias for sessions based on serial nvectors. *)
type 'k serial_session = (RealArray.t, 'k) session
                         constraint 'k = [>Nvector_serial.kind]

(** {2:linear Linear solvers} *)

(** Linear solvers used by Ida.

    @ida <Usage/index.html#linear-solver-interface-functions> Linear solver interface functions *)
type ('data, 'kind) linear_solver = ('data, 'kind) Ida_impl.linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type 'kind serial_linear_solver =
  (Nvector_serial.data, 'kind) linear_solver
  constraint 'kind = [>Nvector_serial.kind]

(** Workspaces with two temporary vectors. *)
type 'd double = 'd * 'd

(** Workspaces with three temporary vectors. *)
type 'd triple = 'd * 'd * 'd

(** Arguments common to Jacobian callback functions.

    @ida IDALsJacFn
    @ida IDALsJacTimesVecFn
    @ida IDALsPrecSolveFn
    @ida IDALsPrecSetupFn *)
type ('t, 'd) jacobian_arg = ('t, 'd) Ida_impl.jacobian_arg =
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

    @ida IDAResFn *)
type 'd resfn = float -> 'd -> 'd -> 'd -> unit

(** Direct Linear Solvers operating on dense, banded and sparse matrices. *)
module Dls : sig (* {{{ *)
  include module type of Sundials_LinearSolver.Direct

  (** Callback functions that compute dense approximations to a Jacobian
      matrix. In the call [jac arg jm], [arg] is a {!jacobian_arg}
      with three work vectors and the computed Jacobian must be stored
      in [jm].

      The callback should load the [(i,j)]th entry of [jm] with
      {% $\frac{\partial F_i}{\partial y_j} + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%},
      i.e., the partial derivative of the [i]th equation with respect to
      the [j]th variable, evaluated at the values of [t], [y], and [y']
      obtained from [arg]. Only nonzero elements need be loaded into [jm].

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the elements of [arg] nor the matrix [jm] should
               be accessed after the function has returned.}

      @ida IDALsJacFn *)
  type 'm jac_fn = (RealArray.t triple, RealArray.t) jacobian_arg -> 'm -> unit

  (** Create an Ida-specific linear solver from a Jacobian approximation
      function and a generic direct linear solver.
      The Jacobian approximation function is optional for dense and banded
      solvers (if not given an internal difference quotient approximation is
      used), but must be provided for other solvers (or [Invalid_argument] is
      raised).

      @ida IDASetLinearSolver
      @ida IDASetJacFn *)
  val solver :
    ?jac:'m jac_fn ->
    ('m, RealArray.t, 'kind, [>`Dls]) LinearSolver.t ->
    'kind serial_linear_solver

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by a direct
      linear solver.

      @ida IDAGetLinWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space : 'k serial_session -> int * int


  (** Returns the number of calls made by a direct linear solver to the
      Jacobian approximation function.

      @ida IDAGetNumJacEvals *)
  val get_num_jac_evals : 'k serial_session -> int

  (** Returns the number of calls to the residual callback due to
      the finite difference Jacobian approximation.

      @ida IDAGetNumLinResEvals *)
  val get_num_lin_res_evals : 'k serial_session -> int

end (* }}} *)

(** Scaled Preconditioned Iterative Linear Solvers. *)
module Spils : sig (* {{{ *)
  include module type of Sundials_LinearSolver.Iterative

  (** {3:precond Preconditioners} *)

  (** Callback functions that solve a linear system involving a
      preconditioner matrix.
      In the call [prec_solve_fn jac r z delta],
      [jac] is a {!jacobian_arg} with no work vectors,
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

      @ida IDALsPrecSolveFn *)
  type 'd prec_solve_fn =
    (unit, 'd) jacobian_arg
    -> 'd
    -> 'd
    -> float
    -> unit

  (** Callback functions that preprocess or evaluate Jacobian-related data
      need by {!prec_solve_fn}. The sole argument is a {!jacobian_arg} with
      no work vectors.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning The elements of the argument should not be accessed after the
               function has returned.}

      @ida IDALsPrecSetupFn *)
  type 'd prec_setup_fn = (unit, 'd) jacobian_arg -> unit

  (** Specifies a preconditioner and its callback functions.
      The following functions and those in {!Ida_bbd} construct
      preconditioners.

      The {!prec_solve_fn} is mandatory. The {!prec_setup_fn} can be
      omitted if not needed.

      @ida IDASetPreconditioner
      @ida IDALsPrecSolveFn
      @ida IDALsPrecSetupFn *)
  type ('d, 'k) preconditioner = ('d, 'k) Ida_impl.SpilsTypes.preconditioner

  (** No preconditioning.  *)
  val prec_none : ('d, 'k) preconditioner

  (** Left preconditioning. {% $Pz = r$%}, where $P$ approximates, perhaps
      crudely,
      {% $J = \frac{\partial F}{\partial y} + c_j\frac{\partial F}{\partial\dot{y}}$%}. *)
  val prec_left :
    ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> ('d, 'k) preconditioner

  (** {3:lsolvers Solvers} *)

  (** Callback functions that preprocess or evaluate Jacobian-related data
      needed by the jac_times_vec_fn. In the call [jac_times_setup_fn arg],
      [arg] is a {!jacobian_arg} with no work vectors.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning The elements of [arg] should not be accessed after the
               function has returned.}

      @ida IDALsJacTimesSetupFn *)
  type 'd jac_times_setup_fn = (unit, 'd) jacobian_arg -> unit

  (** Callback functions that compute the Jacobian times a vector. In the
      call [jac_times_vec_fn arg v jv], [arg] is a {!jacobian_arg} with two
      work vectors, [v] is the vector multiplying the Jacobian, and [jv] is
      the vector in which to store the
      result—{% $\mathtt{jv} = J\mathtt{v}$%}.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the elements of [arg] nor [v] or [jv] should be
               accessed after the function has returned.}

      @ida IDALsJacTimesVecFn *)
  type 'd jac_times_vec_fn =
    ('d double, 'd) jacobian_arg
    -> 'd
    -> 'd
    -> unit

  (** Create an Ida-specific linear solver from a generic iterative
      linear solver.

      The [jac_times_res] argument specifies an alternative DAE residual
      function for use in the internal Jacobian-vector product difference
      quotient approximation. It is incorrect to specify both this argument
      and [jac_times_vec].

      NB: a [jac_times_setup_fn] is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

      NB: a [jac_times_res] function is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 5.3.0.

      @ida IDASetLinearSolver
      @ida IDASetJacTimes
      @ida IDASetJacTimesResFn *)
  val solver :
    ('m, 'd, 'k, [>`Iter]) LinearSolver.t
    -> ?jac_times_vec:'d jac_times_setup_fn option * 'd jac_times_vec_fn
    -> ?jac_times_res:'d resfn
    -> ('d, 'k) preconditioner
    -> ('d, 'k) linear_solver

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by the spils
      linear solver.

      @ida IDAGetLinWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space       : ('d, 'k) session -> int * int

  (** Returns the cumulative number of linear iterations.

      @ida IDAGetNumLinIters *)
  val get_num_lin_iters    : ('d, 'k) session -> int

  (** Returns the cumulative number of linear convergence failures.

      @ida IDAGetNumLinConvFails *)
  val get_num_lin_conv_fails   : ('d, 'k) session -> int

  (** Returns the number of calls to the setup function.

      @ida IDAGetNumPrecEvals *)
  val get_num_prec_evals   : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the preconditioner solve
      function.

      @ida IDAGetNumPrecSolves *)
  val get_num_prec_solves  : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the Jacobian-vector
      setup function.

      @ida IDAGetNumJTSetupEvals
      @since 3.0.0 *)
  val get_num_jtsetup_evals : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the Jacobian-vector
      function.

      @ida IDAGetNumJtimesEvals *)
  val get_num_jtimes_evals : ('d, 'k) session -> int

  (** Returns the number of calls to the residual callback for
      finite difference Jacobian-vector product approximation. This counter is
      only updated if the default difference quotient function is used.

      @ida IDAGetNumLinResEvals *)
  val get_num_lin_res_evals    : ('d, 'k) session -> int

  (** {3:lowlevel Low-level solver manipulation}

      The {!init} and {!reinit} functions are the preferred way to set or
      change preconditioner functions. These low-level functions are provided
      for experts who want to avoid resetting internal counters and other
      associated side-effects. *)

  (** Change the preconditioner functions.

      @ida IDASetPreconditioner
      @ida IDALsPrecSolveFn
      @ida IDALsPrecSetupFn *)
   val set_preconditioner :
     ('d,'k) session
     -> ?setup:'d prec_setup_fn
     -> 'd prec_solve_fn
     -> unit

  (** Change the Jacobian-times-vector function.

      NB: the [jac_times_setup] argument is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

      @ida IDASetJacTimes
      @ida IDALsJacTimesVecFn *)
  val set_jac_times :
    ('d,'k) session
    -> ?jac_times_setup:'d jac_times_setup_fn
    -> 'd jac_times_vec_fn
    -> unit

  (** Remove a Jacobian-times-vector function and use the default
      implementation.

      @ida IDASetJacTimes
      @ida IDAJacTimesVecFn *)
  val clear_jac_times : ('d, 'k) session -> unit
end (* }}} *)

(** Create an IDA-specific linear solver from a generic matrix embedded
    solver.

    @ida IDASetLinearSolver
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

(** {2:init Initialization} *)

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

    @ida IDARootFn *)
type 'd rootsfn = float -> 'd -> 'd -> RealArray.t -> unit

(** Creates and initializes a session with the solver. The call
    {[init linsolv tol ~nlsolver ~nlsresfn ~lsolver f ~varid:varid ~roots:(nroots, g) t0 y0 y'0]}
    has as arguments:
    - [tol],     the integration tolerances,
    - [nlsolver], the solver to use to calculate integration steps
                  and initial conditions,
    - [nlsresfn], alternative residual function to use in
                  nonlinear system function evaluations,
    - [lsolver],  used by [nlsolver]s based on Newton interation,
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

    If an [nlsolver] is not specified, then the
    {{!Sundials_NonlinearSolver.Newton}Newton} module is used by default.
    The [nlsolver] must be of type
    {{!Sundials_NonlinearSolver.nonlinear_solver_type}RootFind}, otherwise an
    {!IllInput} exception is raised.

    The alternative residual function for nonlinear system function
    evaluations is only supported for Sundials >= 5.8.0.

    By default, the session is created using the context returned by
    {!Sundials.Context.default}, but this can be overridden by passing an
    optional [context] argument.

    @ida IDACreate
    @ida IDAInit
    @ida IDARootInit
    @ida IDASetLinearSolver
    @ida IDASetNonlinearSolver
    @ida IDASStolerances
    @ida IDASVtolerances
    @ida IDAWFtolerances
    @ida IDAEwtFn
    @ida IDASetId
    @ida IDASetNlsResFn *)
val init :
       ?context:Context.t
    -> ('d, 'kind) tolerance
    -> ?nlsolver: ('d, 'kind, ('d, 'kind) session, [`Nvec])
                    Sundials_NonlinearSolver.t
    -> ?nlsresfn:'d resfn
    -> lsolver:('d, 'kind) linear_solver
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

    @ida IDASetId *)
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

    @ida IDASetId *)
val set_id : ('d, 'k) session -> ('d,'k) Nvector.t -> unit

(** Indicates whether or not to ignore algebraic variables in the local error
    test. When ignoring algebraic variables ([true]), a [varid] vector must be
    specified either in the call or by a prior call to {!init} or {!set_id}.
    Suppressing local error tests for algebraic variables is {i discouraged}
    for DAE systems of index 1 and {i encouraged} for systems of index 2 or
    more.

    @ida IDASetId
    @ida IDASetSuppressAlg *)
val set_suppress_alg : ('d, 'k) session
                       -> ?varid:('d, 'k) Nvector.t -> bool -> unit

(** Computes the initial state vector for certain index-one problems.
    All components of $y$ are computed, using {% $\dot{y}$%}, to satisfy
    the constraint {% $F(t_0, y_0, \dot{y}_0) = 0$%}. If given, the
    [~y] vector is filled with the corrected values. The last argument is
    the first vale of $t$ at which a solution will be requested.
    A {!reinit} is required before calling this function after
    {!solve_normal} or {!solve_one_step}.

    @ida IDACalcIC (IDA_Y_INIT)
    @ida IDAGetConsistentIC
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

    @ida IDACalcIC (IDA_YA_YDP_INIT)
    @ida IDASetId
    @ida IDAGetConsistentIC
    @raise IdNotSet Variable ids have not been set (see {!set_id}). *)
val calc_ic_ya_yd' :
  ('d, 'k) session
  -> ?y:('d, 'k) Nvector.t
  -> ?y':('d, 'k) Nvector.t
  -> ?varid:('d, 'k) Nvector.t
  -> float
  -> unit

(** Specifies the positive constant in the nonlinear convergence test
    of the initial condition calculation.

    @ida IDASetNonlinConvCoefIC
    @ida IDACalcIC *)
val set_nonlin_conv_coef_ic : ('d, 'k) session -> float -> unit

(** Specifies the maximum number of steps taken in attempting to reach
    a given output time in the initial condition calculation.

    @ida IDASetMaxNumStepsIC
    @ida IDACalcIC (IDA_YA_YDP_INIT) *)
val set_max_num_steps_ic : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of approximate Jacobian or preconditioner
    evaluations allowed when the Newton iteration appears to be slowly
    converging.

    @ida IDASetMaxNumJacsIC
    @ida IDACalcIC *)
val set_max_num_jacs_ic : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of Newton iterations allowed in any one
    attempt to calculate initial conditions.

    @ida IDASetMaxNumItersIC
    @ida IDACalcIC *)
val set_max_num_iters_ic : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of linesearch backtracks allowed in any
    Newton iteration, when solving the initial conditions calculation
    problem.

    @ida IDASetMaxBacksIC
    @ida IDACalcIC
    @raise Config.NotImplementedBySundialsVersion Feature not available.
    @since 2.7.0 *)
val set_max_backs_ic : ('d, 'k) session -> int -> unit

(** Enables ([true]) or disables ([false]) the linesearch algorithm
    in the initial condition calculation.

    @ida IDASetLineSearchOffIC
    @ida IDACalcIC *)
val set_line_search_ic : ('d, 'k) session -> bool -> unit

(** Specifies a positive lower bound on the Newton step in the initial condition
    calculation.

    @ida IDASetStepToleranceIC
    @ida IDACalcIC *)
val set_step_tolerance_ic : ('d, 'k) session -> float -> unit

(** Returns the number of backtrack operations during {!calc_ic_ya_yd'} or
    {!calc_ic_y}.

    @ida IDAGetNumBacktrackOps *)
val get_num_backtrack_ops : ('d, 'k) session -> int

(** {2:solver Solution} *)

(** Values returned by the step functions. Failures are indicated by
    exceptions.

    @ida IDASolve *)
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

    @ida IDASolve (IDA_NORMAL)
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

    @ida IDASolve (IDA_ONE_STEP) *)
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

    @ida IDAGetDky
    @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
    @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
val get_dky : ('d, 'k) session -> ('d, 'k) Nvector.t -> float -> int -> unit

(** Reinitializes the solver with new parameters and state values. The
    values of the independent variable, i.e., the simulation time, the
    state variables, and the derivatives must be given.
    If the argument [~linsolv] is not given, the current linear solver
    remains unchanged. The argument [~roots] works similarly; pass
    {!no_roots} to disable root finding. The [~resfn] argument gives the
    possibility to change the residual function. Similarly for [~nlsresfn]
    and the alternative residual function for nonlinear system function
    iterations.

    @ida IDAReInit
    @ida IDASetLinearSolver
    @ida IDASetNonlinearSolver
    @ida IDASetNlsResFn *)
val reinit :
  ('d, 'k) session
  -> ?nlsolver:('d, 'k, ('d, 'k) session, [`Nvec]) Sundials_NonlinearSolver.t
  -> ?nlsresfn:'d resfn
  -> ?lsolver:('d, 'k) linear_solver
  -> ?roots:(int * 'd rootsfn)
  -> ?resfn:'d resfn
  -> float
  -> ('d, 'k) Nvector.t
  -> ('d, 'k) Nvector.t
  -> unit

(** {2:set Modifying the solver} *)

(** Set the integration tolerances.

    @ida IDASStolerances
    @ida IDASVtolerances
    @ida IDAWFtolerances
    @ida IDAEwtFn *)
val set_tolerances : ('d, 'k) session -> ('d, 'k) tolerance -> unit

(** Specifies a vector defining inequality constraints for each
    component of the solution vector [u].  See {!Sundials.Constraint}.

    @ida IDASetConstraints *)
val set_constraints : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Disables constraint checking.

    @ida IDASetConstraints *)
val clear_constraints : ('d, 'k) session -> unit

(** {3:set_main Main solver optional input functions} *)

(** Configure the default error handler to write messages to a file.
    By default it writes to Logfile.stderr.

    @ida IDASetErrFile *)
val set_error_file : ('d, 'k) session -> Logfile.t -> unit

(** Specifies a custom function for handling error messages.
    The handler must not fail: any exceptions are trapped and discarded.

    @ida IDASetErrHandlerFn
    @ida IDAErrHandlerFn *)
val set_err_handler_fn
  : ('d, 'k) session -> (Util.error_details -> unit) -> unit

(** Restores the default error handling function.

    @ida IDASetErrHandlerFn *)
val clear_err_handler_fn : ('d, 'k) session -> unit

(** Specifies the maximum order of the linear multistep method.

    @ida IDASetMaxOrd *)
val set_max_ord : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of steps taken in attempting to reach
    a given output time.

    @ida IDASetMaxNumSteps *)
val set_max_num_steps : ('d, 'k) session -> int -> unit

(** Specifies the initial step size.

    @ida IDASetInitStep *)
val set_init_step : ('d, 'k) session -> float -> unit

(** Specifies an upper bound on the magnitude of the step size.

    @ida IDASetMaxStep *)
val set_max_step : ('d, 'k) session -> float -> unit

(** Specifies the minimum absolute value of the step size.

    @ida IDASetMinStep
    @since 6.2.0 *)
val set_min_step: ('d, 'k) session -> float -> unit

(** Limits the value of the independent variable [t] when solving.
    By default no stop time is imposed.

    @ida IDASetStopTime *)
val set_stop_time : ('d, 'k) session -> float -> unit

(** Specifies the maximum number of error test failures permitted in attempting
    one step.

    @ida IDASetMaxErrTestFails *)
val set_max_err_test_fails : ('d, 'k) session -> int -> unit

(** {3:set_nls Linear solver interface optional input functions} *)

(** Specifies the parameter {% $\delta_\mathrm{cj}$ %} that determines the
    bounds on changes that require a linear solver setup call.
    If {% $\frac{\mathrm{cj}_\mathrm{current}}{\mathrm{cj}_\mathrm{previous}}
           < \frac{1 - \delta_\mathrm{cj}}{1 + \delta_\mathrm{cj}} $ %}
    or {% $\frac{\mathrm{cj}_\mathrm{current}}{\mathrm{cj}_\mathrm{previous}}
           > \frac{1 + \delta_\mathrm{cj}}{1 - \delta_\mathrm{cj}} $ %}, the
    linear setup function is called.
    If the argument is [< 0] or [>= 1], then the default value of 0.25 is
    used.

    @ida IDASetDeltaCjLSetup
    @since 6.2.0 *)
val set_delta_cj_lsetup : ('d, 'k) session -> float -> unit

(** Sets the factor by which the Krylov linear solver's convergence test
    constant is reduced from the Newton iteration test constant.
    This factor must be >= 0; passing 0 specifies the default (0.05).

    @ida IDASetEpsLin *)
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

    @ida IDASetLSNormFactor
    @since 5.4.0 *)
val set_ls_norm_factor : ('d, 'k) session -> float -> unit

(** Enables or disables scaling of the linear system solution to account
    for a change in {% $\gamma$ %} in the linear system.
    Linear solution scaling is enabled by default when a matrix-based
    linear solver is attached.

    @ida IDASetLinearSolutionScaling
    @since 5.2.0 *)
val set_linear_solution_scaling : ('d, 'k) session -> bool -> unit

(** Sets the increment factor ([dqincfac]) to use in the difference-quotient
    approximation.

    Specifically, the  product {% $Jv$ %} is approximated
    by {% $Jv = \frac{1}{\sigma}\left(
                  F(t, \tilde{y}, \tilde{y}') - F(t, y, y')
                \right)$ %}.
    where {% $\tilde{y} = y + \sigma v$ %},
          {% $\tilde{y}' = y' + c_j \sigma v$ %},
    {% $c_j$ %} is a BDF parameter proportional to the step size,
    {% $\sigma = \sqrt{N} \mathtt{dqincfac}$ %},
    and {% $N$ %} is the number of equations in the DAE system.

    @ida IDASetIncrementFactor *)
val set_increment_factor : ('d, 'k) session -> float -> unit

(** {3:set_nls Nonlinear solver interface optional input functions} *)

(** Specifies the maximum number of nonlinear solver iterations permitted per
    step.

    @ida IDASetMaxNonlinIters *)
val set_max_nonlin_iters : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver convergence failures
    permitted during one step.

    @ida IDASetMaxConvFails *)
val set_max_conv_fails : ('d, 'k) session -> int -> unit

(** Specifies the safety factor used in the nonlinear convergence test.

    @ida IDASetNonlinConvCoef *)
val set_nonlin_conv_coef : ('d, 'k) session -> float -> unit

(** {3:set_tsa Time step adaptivity optional input functions} *)

(** Specifies the lower {% $\eta_\mathrm{min\_fx}$ %} and
    upper {% $\eta_\mathrm{max\_fx}$ %} bounds in which the step size will
    remain unchanged. The default values are {% $\eta_\mathrm{min\_fx} = 1$ %}
    and {% $\eta_\mathrm{max\_fx} = 2$ %}.
    If the step size change factor satisfies,
    {% $\eta_\mathrm{min\_fx} < \eta < \eta_\mathrm{max\_fx}$ %},
    the current step size is retained.

    If [min < 0] or [min >= 1], the default minimum value is used.
    If [max < 1], the default maximum value is used.

    @cvode IDASetEtaFixedStepBounds
    @since 6.2.0 *)
val set_eta_fixed_step_bounds
      : ('d, 'k) session -> min:float -> max:float -> unit -> unit

(** Specifies the minimum step size reduction factor, {% $\eta_\mathrm{min}$ %}.
    If the argument is [<= 0] or [>= 1], the default value of 0.5 is used.

    @cvode IDASetEtaMin
    @since 6.2.0 *)
val set_eta_min                      : ('d, 'k) session -> float -> unit

(** Specifies the maximum step size growth factor, {% $\eta_\mathrm{max}$ %}.
    The default value is 2.

    @cvode IDASetEtaMax
    @since 6.2.0 *)
val set_eta_max                      : ('d, 'k) session -> float -> unit

(** Specifies the maximum step size reduction factor, {% $\eta_\mathrm{low}$ %}.
    If the argument is [< 0] or [> 1], the default value of 0.9 is used.

    @cvode IDASetEtaLow
    @since 6.2.0 *)
val set_eta_low                      : ('d, 'k) session -> float -> unit

(** Specifies the minimum step size reduction factor after an error test
    failure, {% $\eta_\mathrm{min\_ef}$ %}.
    If the argument is [<= 0] or [>= 1], the default value of 0.25 is used.

    @cvode IDASetEtaMinErrFail
    @since 6.2.0 *)
val set_eta_min_err_fail             : ('d, 'k) session -> float -> unit

(** Specifies the step size factor after a nonlinear solver convergence
    failure, {% $\eta_\mathrm{cf}$ %}.
    If the argument is [<= 0] or [>= 1], the default value of 0.25 is used.

    @cvode IDASetEtaConvFail
    @since 6.2.0 *)
val set_eta_conv_fail                : ('d, 'k) session -> float -> unit

(** {2:get Querying the solver (optional output functions)} *)

(** Returns the sizes of the real and integer workspaces.

    @ida IDAGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space          : ('d, 'k) session -> int * int

(** Returns the cumulative number of internal solver steps.

    @ida IDAGetNumSteps *)
val get_num_steps           : ('d, 'k) session -> int

(** Returns the number of calls to the residual function.

    @ida IDAGetNumResEvals *)
val get_num_res_evals       : ('d, 'k) session -> int

(** Returns the number of calls made to the linear solver's setup function.

    @ida IDAGetNumLinSolvSetups *)
val get_num_lin_solv_setups : ('d, 'k) session -> int

(** Returns the number of local error test failures that have occurred.

    @ida IDAGetNumErrTestFails *)
val get_num_err_test_fails  : ('d, 'k) session -> int

(** Returns the number of failed steps due to a nonlinear solver failure.

    @cvode IDAGetNumStepSolveFails
    @since 6.2.0 *)
val get_num_step_solve_fails : ('d, 'k) session -> int

(** Returns the integration method order used during the last internal step.

    @ida IDAGetLastOrder *)
val get_last_order          : ('d, 'k) session -> int

(** Returns the integration method order to be used on the next internal step.

    @ida IDAGetCurrentOrder *)
val get_current_order       : ('d, 'k) session -> int

(** Returns the integration step size taken on the last internal step.

    @ida IDAGetLastStep *)
val get_last_step           : ('d, 'k) session -> float

(** Returns the integration step size to be attempted on the next internal step.

    @ida IDAGetCurrentStep *)
val get_current_step        : ('d, 'k) session -> float

(** Returns the the value of the integration step size used on the first step.

    @ida IDAGetActualInitStep *)
val get_actual_init_step    : ('d, 'k) session -> float

(** Returns the the current internal time reached by the solver.

    @ida IDAGetCurrentTime *)
val get_current_time        : ('d, 'k) session -> float

(** Returns a suggested factor by which the user's tolerances should be scaled
    when too much accuracy has been requested for some internal step.

    @ida IDAGetTolScaleFactor *)
val get_tol_scale_factor : ('d, 'k) session -> float

(** Returns the solution error weights at the current time.

    @ida IDAGetErrWeights *)
val get_err_weights : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Returns the vector of estimated local errors.

    @ida IDAGetEstLocalErrors *)
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

    @ida IDAGetIntegratorStats *)
val get_integrator_stats    : ('d, 'k) session -> integrator_stats

(** Prints the integrator statistics on the given channel.

    @ida IDAGetIntegratorStats *)
val print_integrator_stats  : ('d, 'k) session -> out_channel -> unit


(** Returns the cumulative number of nonlinear (functional or Newton)
    iterations.

    @ida IDAGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('d, 'k) session -> int

(** Returns the cumulative number of nonlinear convergence failures.

    @ida IDAGetNumNonlinSolvConvFails *)
val get_num_nonlin_solv_conv_fails : ('d, 'k) session -> int

(** Returns both the numbers of nonlinear iterations performed [nniters] and
    nonlinear convergence failures [nncfails].

    @ida IDAGetNonlinSolvStats
    @return ([nniters], [nncfails]) *)
val get_nonlin_solv_stats : ('d, 'k) session -> int *int

(** Outputs all of the integrator, nonlinear solver, linear solver, and other
    statistics.

    @ida IDAPrintAllStats
    @since 6.2.0 *)
val print_all_stats
      : ('d, 'k) session -> Logfile.t -> Sundials.output_format -> unit

(** {2:roots Additional root-finding functions} *)

(** [set_root_direction s dir] specifies the direction of zero-crossings to be
    located and returned. [dir] may contain one entry for each root function.

    @ida IDASetRootDirection *)
val set_root_direction : ('d, 'k) session -> RootDirs.d array -> unit

(** Like {!set_root_direction} but specifies a single direction for all root
    functions.

    @ida IDASetRootDirection *)
val set_all_root_directions : ('d, 'k) session -> RootDirs.d -> unit

(** Disables issuing a warning if some root function appears to be identically
    zero at the beginning of the integration.

    @ida IDASetNoInactiveRootWarn *)
val set_no_inactive_root_warn : ('d, 'k) session -> unit

(** Returns the number of root functions. *)
val get_num_roots : ('d, 'k) session -> int

(** Fills an array showing which functions were found to have a root.

    @ida IDAGetRootInfo *)
val get_root_info : ('d, 'k) session -> Roots.t -> unit

(** Returns the cumulative number of calls made to the user-supplied root
    function g.

    @ida IDAGetNumGEvals *)
val get_num_g_evals : ('d, 'k) session -> int

(** {2:nls Advanced nonlinear solver functions}

    These advanced functions may be useful for writing customized nonlinear
    solver routines. *)

(** Returns the scalar {% $c_j$ %}, which is proportional to the inverse of
    the step size.

    @ida IDAGetCurrentCj
    @since 5.0.0 *)
val get_current_cj : ('d, 'k) session -> float

(** Returns the current {% $y$ %} vector. This vector provides direct access
    to the data within the integrator.

    @ida IDAGetCurrentY
    @since 5.0.0 *)
val get_current_y : ('d, 'k) session -> 'd

(** Returns the current {% $\dot{y}$ %} vector. This vector provides direct
    access to the data within the integrator.

    @ida IDAGetCurrentYp
    @since 5.0.0 *)
val get_current_yp : ('d, 'k) session -> 'd

(** Internal data required to construct the current nonlinear implicit
    system within a nonlinear solver. *)
type 'd nonlin_system_data = {
  tn     : float;
    (** Independent variable value {% $t_n$ %}. *)
  yypred : 'd;
    (** Predicted value of {% $y_{\mathit{pred}}$ %} at {% $t_n$ %}. This
        data must not be changed. *)
  yppred : 'd;
    (** Predicted value of {% $\dot{y}_{\mathit{pred}}$ %} at {% $t_n$ %}. This
        data must not be changed. *)
  yyn    : 'd;
    (** The vector {% $y_n$ %}. This data may not be current and may
        need to be filled. *)
  ypn    : 'd;
    (** The vector {% $\dot{y}_n$ %}. This data may not be current and may
        need to be filled. *)
  res   : 'd;
    (** The residual function evaluated at the current time and state,
        {% $F(t_n, y_n, \dot{y}_n)$ %}. * This data may not be current and
        may need to be filled. *)
  cj    : float;
    (** The scalar {% $c_j$ %} which is proportional to the inverse of
        the step size {% $\alpha$ %}. *)
}

(** Gives direct access to the internal data required to construct the
    current nonlinear system within a nonlinear solver. This
    function should be called inside the nonlinear system function.
    If the nonlinear solver uses the [lsetup] or [lsolve] functions, then
    the nonlinear solver system function must fill the [yyn], [ypn], and [res]
    vectors with, respectively:
    {% $\mathit{yyn} = y_{\mathit{pred}} + y_{\mathit{cor}}$ %},
    {% $\mathit{ypn} = \dot{y}_{\mathit{pred}} + \alpha\dot{y}_{\mathit{cor}}$ %},
    and
    {% $\mathit{res} = F(t_n, y_n, \dot{y}_n)$ %} where
    {% $y_{\mathit{cor}}$ %} is the
    first argument of the nonlinear solver system function. Within a custom
    linear solver, then the vectors [yyn], [ypn], and [res] are only current
    after an evaluation of the nonlinear system function.

    @ida IDAGetNonlinearSystemData
    @since 5.4.0 *)
val get_nonlin_system_data : ('d, 'k) session -> 'd nonlin_system_data

(** Computes the current {% $y$ %} vector from a correction vector.

    @ida IDAComputeY
    @since 5.0.0 *)
val compute_y
  : ('d, 'k) session -> ycor:('d, 'k) Nvector.t -> y:('d, 'k) Nvector.t -> unit

(** Computes the current {% $\dot{y}$ %} vector from a correction vector.

    @ida IDAComputeYp
    @since 5.0.0 *)
val compute_yp
  : ('d, 'k) session -> ycor:('d, 'k) Nvector.t -> yp:('d, 'k) Nvector.t -> unit

(** {2:exceptions Exceptions} *)

(** Raised on missing or illegal solver inputs. Also raised if an element
    of the error weight vector becomes zero during time stepping, or the
    linear solver initialization function failed, or a root was found both at
    [t] and very near [t].

    @ida <Constants_link.html> IDA_ILL_INPUT *)
exception IllInput

(** The requested time could not be reached in [mxstep] internal steps.
    See {!set_max_num_steps}

    @ida <Constants_link.html> IDA_TOO_MUCH_WORK *)
exception TooMuchWork

(** The requested accuracy could not be satisfied.

    @ida <Constants_link.html> IDA_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** Too many error test failures within a step. See {!set_max_err_test_fails}.

    @ida <Constants_link.html> IDA_ERR_FAIL *)
exception ErrFailure

(** Too many convergence test failures within a step,
    or Newton convergence failed. See {!set_max_conv_fails}.

    @ida <Constants_link.html> IDA_CONV_FAIL *)
exception ConvergenceFailure

(** Linear solver initialization failed.

    @ida <Constants_link.html> IDA_LINIT_FAIL *)
exception LinearInitFailure

(** Linear solver setup failed in an unrecoverable manner.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.PSetFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @ida IDAGetLastLinFlag
    @ida <Constants_link.html> IDA_LSETUP_FAIL *)
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

    @ida IDAGetLastLinFlag
    @ida <Constants_link.html> IDA_LSOLVE_FAIL *)
exception LinearSolveFailure of exn option

(** The nonlinear solver failed in a general way.

    @ida <Constants_link.html> CV_NLS_FAIL
    @since 5.0.0 *)
exception NonlinearSolverFailure

(** Nonlinear solver initialization failed.

    @ida <Constants_link.html> IDA_NLS_INIT_FAIL *)
exception NonlinearInitFailure

(** Nonlinear solver setup failed in an unrecoverable manner.

    @ida <Constants_link.html> IDA_NLS_SETUP_FAIL *)
exception NonlinearSetupFailure

(** Nonlinear solver setup failed in a recoverable manner.

    @ida <Constants_link.html> IDA_NLS_SETUP_RECVR *)
exception NonlinearSetupRecoverable

(** The residual function failed in an unrecoverable manner.

    @ida <Constants_link.html> IDA_RES_FAIL *)
exception ResFuncFailure

(** The residual function had a recoverable error when first called.

    @ida <Constants_link.html> IDA_FIRST_RES_FAIL *)
exception FirstResFuncFailure

(** Too many convergence test failures, or unable to estimate the initial step
    size, due to repeated recoverable errors in the residual function.

    @ida <Constants_link.html> IDA_REP_RES_ERR *)
exception RepeatedResFuncFailure

(** The rootfinding function failed.

    @ida <Constants_link.html> IDA_RTFUNC_FAIL *)
exception RootFuncFailure

(** No solution satisfying the inequality constraints could be found.

    @ida <Constants_link.html> IDA_CONSTR_FAIL *)
exception ConstraintFailure

(** Linesearch could not find a solution with a step larger than steptol in
    weighted RMS norm.

    @ida <Constants_link.html> IDA_LINESEARCH_FAIL *)
exception LinesearchFailure

(** A recoverable error occurred in a callback but no recovery was possible.

    @ida <Constants_link.html> IDA_NO_RECOVERY *)
exception NoRecovery

(** A component of the error weight vector, either for the input value or a
    corrected value, is zero.

    @ida <Constants_link.html> IDA_BAD_EWT *)
exception BadEwt

(** Raised by {!get_dky} for invalid order values.

    @ida IDAGetDky (IDA_BAD_K) *)
exception BadK

(** Raised by {!get_dky} for invalid time values.

    @ida IDAGetDky (IDA_BAD_T) *)
exception BadT

(** Variable ids are required but not set. See {!set_id}. *)
exception IdNotSet

(** A fused vector operation failed.

    @ida <Constants_link.html> IDA_VECTOROP_ERR *)
exception VectorOpErr

