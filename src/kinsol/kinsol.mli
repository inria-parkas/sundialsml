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

(** Solves nonlinear systems using Newton-Krylov techniques.

    This module solves numerically problems of the form
    {% $F(u) = 0$%} given an initial guess $u_0$.

    This documented interface is structured as follows.
    {ol
      {- {{:#linear}Linear solvers}}
      {- {{:#solver}Solver initialization and use}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#exceptions}Exceptions}}}

     @version VERSION()
     @author Timothy Bourke (Inria/ENS)
     @author Jun Inoue (Inria/ENS)
     @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** A session with the KINSOL solver.

    An example session with Kinsol ({openfile kinsol_skel.ml}): {[
#include "../../examples/ocaml/skeletons/kinsol_skel.ml"
    ]}

    @kinsol <Usage/index.html#a-skeleton-of-the-users-main-program> Skeleton of main program *)
type ('data, 'kind) session = ('data, 'kind) Kinsol_impl.session

(** Alias for sessions based on serial nvectors. *)
type 'kind serial_session = (Nvector_serial.data, 'kind) session
                            constraint 'kind = [>Nvector_serial.kind]

(** {2:linear Linear solvers} *)

(** Linear solvers used by Kinsol.

    @kinsol <Usage/index.html#linear-solver-specification-functions> Linear solver specification functions *)
type ('data, 'kind) linear_solver = ('data, 'kind) Kinsol_impl.linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type 'kind serial_linear_solver =
  (Nvector_serial.data, 'kind) linear_solver
  constraint 'kind = [>Nvector_serial.kind]

(** Workspaces with two temporary vectors. *)
type 'd double = 'd * 'd

(** Arguments common to Jacobian callback functions.

    @kinsol KINLsJacFn
    @kinsol KINLsPrecSolveFn
    @kinsol KINLsPrecSetupFn *)
type ('t, 'd) jacobian_arg = ('t, 'd) Kinsol_impl.jacobian_arg =
  {
    jac_u   : 'd;   (** The current unscaled iterate. *)
    jac_fu  : 'd;   (** The current value of the vector $F(u)$. *)
    jac_tmp : 't    (** Workspace data. *)
  }

(** System function that defines nonlinear problem. The call
    [sysfun u fval] must calculate $F(u)$ into [fval] using the current value
    vector [u].

     Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
     Any other exception is treated as an unrecoverable error.

    {warning [u] and [fval] should not be accessed after the function
             returns.}

    @kinsol KINSysFn *)
type 'data sysfn = 'data -> 'data -> unit

(** Direct Linear Solvers operating on dense, banded, and sparse matrices. *)
module Dls : sig (* {{{ *)
  include module type of Sundials_LinearSolver.Direct

  (** Callback functions that compute dense approximations to a Jacobian
      matrix. In the call [jac arg jm], [arg] is a {!jacobian_arg}
      with two work vectors and the computed Jacobian must be stored
      in [jm].

      The callback should load the [(i,j)]th entry of [jm] with
      {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of the
      [i]th equation with respect to the [j]th variable, evaluated at the
      values of [t] and [y] obtained from [arg]. Only nonzero elements need
      be loaded into [jm].

      {warning Neither the elements of [arg] nor the matrix [jm] should
               be accessed after the function has returned.}

      @kinsol KINLsJacFn *)
  type 'm jac_fn =
    (RealArray.t double, RealArray.t) jacobian_arg -> 'm -> unit

  (** Create a Kinsol-specific linear solver from a Jacobian approximation
      function and a generic direct linear solver.
      The Jacobian approximation function is optional for dense and banded
      solvers (if not given an internal difference quotient approximation is
      used), but must be provided for other solvers (or [Invalid_argument] is
      raised).

      @kinsol KINSetLinearSolver
      @kinsol KINSetJacFn *)
  val solver :
    ?jac:'m jac_fn ->
    ('m, RealArray.t, 'kind, [>`Dls]) LinearSolver.t ->
    'kind serial_linear_solver

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by a direct
      linear solver.

      @kinsol KINGetLinWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space : 'k serial_session -> int * int

  (** Returns the number of calls made by a direct linear solver to the
      Jacobian approximation function.

      @kinsol KINGetNumJacEvals *)
  val get_num_jac_evals : 'k serial_session -> int

  (** Returns the number of calls made by a direct linear solver to the
      user system function for computing the difference quotient approximation
      to the Jacobian.

      @kinsol KINGetNumLinFuncEvals *)
  val get_num_lin_func_evals : 'k serial_session -> int

end (* }}} *)

(** Scaled Preconditioned Iterative Linear Solvers. *)
module Spils : sig (* {{{ *)
  include module type of Sundials_LinearSolver.Iterative

  (** {3:precond Preconditioners} *)

  (** Arguments passed to the preconditioner solver function.

      @kinsol KINLsPrecSolveFn *)
  type 'data solve_arg =
    {
      uscale : 'data; (** Diagonal elements of the scaling matrix for [u]. *)
      fscale : 'data; (** Diagonal elements of the scaling matrix
                          for [fval]. *)
    }

  (** Callback functions that solve a linear system involving a
      preconditioner matrix. The call [prec_solve_fn jarg sarg v] must solve
      {% $Pz = r$%}, where [jarg] is a {!jacobian_arg} with no work vectors,
      [sarg] is a {!solve_arg} giving the scaling matrices, and
      [v] is initially the right-hand side vector $r$ and is later filled
      with the computed solution $z$.
      $P$ is a preconditioner matrix that approximates the system
      Jacobian {% $J = \frac{\partial F}{\partial u}$%}.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the elements of [jarg] or [sarg], nor [z] should be
               accessed after the function has returned.}

      @kinsol KINSetPreconditioner
      @kinsol KINLsPrecSolveFn *)
  type 'd prec_solve_fn =
    (unit, 'd) jacobian_arg
    -> 'd solve_arg
    -> 'd
    -> unit

  (** Callback functions that preprocess or evaluate Jacobian-related data
      need by {!prec_solve_fn}. In the call [prec_setup_fn jarg sarg],
      [jarg] is a {!jacobian_arg} with no work vectors and [sarg] is a
      {!solve_arg} giving the scaling matrices.

      The callback should raise an exception if unsuccessful.

      {warning The elements of [jarg] and [sarg] should not be accessed after
               the function has returned.}

      @kinsol KINLsPrecSetupFn
      @kinsol KINSetPreconditioner
   *)
  type 'd prec_setup_fn =
    (unit, 'd) jacobian_arg
    -> 'd solve_arg
    -> unit

  (** Specifies a preconditioner, including the type of preconditioning
      (none or right) and callback functions.
      The following functions and those in {!Kinsol_bbd} construct
      preconditioners.

      @kinsol KINSetPreconditioner
      @kinsol KINLsPrecSolveFn
      @kinsol KINLsPrecSetupFn *)
  type ('d, 'k) preconditioner = ('d, 'k) Kinsol_impl.SpilsTypes.preconditioner

  (** No preconditioning.  *)
  val prec_none : ('d, 'k) preconditioner

  (** Right preconditioning. The {!prec_setup_fn} should compute
      the right preconditioner matrix $P$ which is used to form
      the scaled preconditioned linear system
      {% $(D_F J(u) P^{-1} D_u^{-1} \cdot (D_u P x) = * -D_F F(u)$%}. *)
  val prec_right :
    ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> ('d, 'k) preconditioner

  (** {3:lsolvers Solvers} *)

  (** Callback functions that compute (an approximation to) the Jacobian
      times a vector. In the call [jac_times_vec_fn v jv u new_u], [v] is the
      vector multiplying the Jacobian, [jv] is the vector in which to store
      the result—{% $\mathtt{jv} = J\mathtt{v}$%}—, [u] is the current
      value of the dependent variable vector, and [new_u=true] indicates that
      the Jacobian data should be recomputed. Returning [false] requests an
      update of the Jacobian data at the next call.

      {warning [v], [jv], and [u] should not be accessed after the function
               has returned.}

      @kinsol KINLsJacTimesVecFn
   *)
  type 'data jac_times_vec_fn =
    'data      (* v *)
    -> 'data   (* jv *)
    -> 'data   (* u *)
    -> bool    (* new_u *)
    -> bool

  (** Create a Kinsol-specific linear solver from a generic iterative
      linear solver.

      The [jac_times_sys] argument specifies an alternative system
      function for use in the internal Jacobian-vector product difference
      quotient approximation. It is incorrect to specify both this argument
      and [jac_times_vec].

      NB: a [jac_times_setup_fn] is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

      NB: a [jac_times_sys] function is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 5.3.0.

      @kinsol KINSetLinearSolver
      @kinsol KINSetJacTimesVecFn
      @kinsol KINSetJacTimesVecSysFn *)
  val solver :
    ('m, 'd, 'k, [>`Iter]) LinearSolver.t
    -> ?jac_times_vec:'d jac_times_vec_fn
    -> ?jac_times_sys:'d sysfn
    -> ('d, 'k) preconditioner
    -> ('d, 'k) linear_solver

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by the spils
      linear solver.

      @kinsol KINGetLinWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space       : ('d, 'k) session -> int * int

  (** Returns the cumulative number of linear iterations.

      @kinsol KINGetNumLinIters *)
  val get_num_lin_iters    : ('d, 'k) session -> int

  (** Returns the cumulative number of linear convergence failures.

      @kinsol KINGetNumLinConvFails *)
  val get_num_lin_conv_fails   : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the setup function.

      @kinsol KINGetNumPrecEvals *)
  val get_num_prec_evals   : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the preconditioner solve
      function.

      @kinsol KINGetNumPrecSolves *)
  val get_num_prec_solves  : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the Jacobian-vector
      function.

      @kinsol KINGetNumJtimesEvals *)
  val get_num_jtimes_evals : ('d, 'k) session -> int

  (** Returns the number of calls to the system function for finite
      difference quotient Jacobian-vector product approximations. This
      counter is only updated if the default difference quotient function
      is used.

      @kinsol KINGetNumLinFuncEvals *)
  val get_num_lin_func_evals    : ('d, 'k) session -> int

  (** {3:lowlevel Low-level solver manipulation} *)

  (** Change the preconditioner functions.

      @kinsol KINSetPreconditioner
      @kinsol KINLsPrecSetupFn
      @kinsol KINLsPrecSolveFn *)
  val set_preconditioner :
    ('d, 'k) session
    -> ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> unit

  (** Change the Jacobian-times-vector function.

      @kinsol KINSetJacTimesVecFn
      @kinsol KINLsJacTimesVecFn *)
  val set_jac_times :
    ('d, 'k) session
    -> 'd jac_times_vec_fn
    -> unit

  (** Remove a Jacobian-times-vector function and use the default
      implementation.

      @kinsol KINSetJacTimesVecFn *)
  val clear_jac_times : ('d, 'k) session -> unit
end (* }}} *)

(** Create a Kinsol-specific linear solver from a generic matrix embedded
    solver.

    @kinsol KINSetLinearSolver
    @since 5.8.0 *)
val matrix_embedded_solver :
       (unit, 'data, 'kind, [>`MatE]) LinearSolver.t
    -> ('data, 'kind) linear_solver

(** {2:solver Solver initialization and use} *)

(** Orthogonalization routines of the QR factorization portion of Anderson
    acceleration.

    @kinsol KINSetOrthAA
    @since 6.0.0 *)
type orthaa =
  | MGS     (** Modified Gram Schmidt (default)
                {cconst KIN_ORTH_MGS} *)
  | ICWY    (** Inverse Compact WY Modified Gram Schmidt
                {cconst KIN_ORTH_ICWY} *)
  | CGS2    (** Classical Gram Schmidt with Reorthogonalization
                {cconst KIN_ORTH_CGS2} *)
  | DCGS2   (** Classical Gram Schmidt with Delayed Reorthogonlization
                {cconst KIN_ORTH_DCGS2} *)

(** Creates and initializes a session with the Kinsol solver. The call
    [init ~max_lin_iters:mli ~maa:maa ~linsolv:ls f tmpl] has as arguments:
     - [mli], the maximum number of nonlinear iterations allowed,
     - [maa], the size of the Anderson acceleration subspace for the
              {{!strategy}Picard} and {{!strategy}FixedPoint} strategies,
     - [orthaa], specifies the othogonalization routine to be used in the QR
                 factorization portion of Anderson acceleration,
     - [ls], the linear solver to use (required for the {{!strategy}Newton},
             {{!strategy}LineSearch}, and {{!strategy}Picard} strategies),
     - [f],       the system function of the nonlinear problem, and,
     - [tmpl]     a template to initialize the session (e.g., the
                  initial guess vector).

     By default, the session is created using the context returned by
     {!Sundials.Context.default}, but this can be overridden by passing
     an optional [context] argument.

     The [orthaa] argument is ignored in Sundials < 6.0.0.

     @kinsol KINCreate
     @kinsol KINInit
     @kinsol KINSetNumMaxIters
     @kinsol KINSetMAA
     @kinsol KINSetOrthAA *)
val init :
     ?context:Context.t
  -> ?max_iters:int
  -> ?maa:int
  -> ?orthaa:orthaa
  -> ?lsolver:('data, 'kind) linear_solver
  -> 'data sysfn
  -> ('data, 'kind) Nvector.t
  -> ('data, 'kind) session

(** Strategy used to solve the nonlinear system. *)
type strategy =
  | Newton            (** Basic Newton iteration. {cconst KIN_NONE} *)
  | LineSearch        (** Newton iteration with globalization.
                          {cconst KIN_LINESEARCH} *)
  | Picard            (** Picard iteration with Anderson Acceleration.
                          {cconst KIN_PICARD} *)
  | FixedPoint        (** Fixed-point iteration with Anderson Acceleration.
                          {cconst KIN_FP} *)

(** Results of nonlinear solution attempts. *)
type result =
  | Success           (** The scaled norm of $F(u)$ is less than [fnormtol].
                          See {!set_func_norm_tol}. {cconst KIN_SUCCESS} *)
  | InitialGuessOK    (** The initial guess already satisfies the system.
                          {cconst KIN_INITIAL_GUESS_OK} *)
  | StoppedOnStepTol  (** Stopped based on scaled step length. The
                          current iterate is an approximate solution, or the
                          algorithm stalled near an invalid solution, or
                          [scsteptol] is too large
                          (see {!set_scaled_step_tol }).
                          {cconst KIN_STEP_LT_STPTOL} *)

(** Computes an approximate solution to a nonlinear system. The call
    [solve s u strategy u_scale f_scale] has arguments:
    - [s], a solver session,
    - [u], an initial guess that is replaced with an approximate solution
           for $F(u) = 0$,
    - [strategy], strategy used to solve the nonlinear system,
    - [u_scale], the diagonal elements of the scaling matrix $D_u$ for
                 vector [u] chosen so that all $D_u u$ all have roughly the
                 same magnitude when [u] is close to a root of $F(u)$, and,
    - [f_scale], the diagonal elements of the scaling matrix $D_f$ for
                 $F(u)$ chosen so that all $D_f F(u)$ have roughtly the same
                 magnitude when [u] is not near a root of $F(u)$.

    The function either returns a {!result} or raises one of the exceptions
    listed below.

    @kinsol KINSol
    @raise MissingLinearSolver A linear solver is required but was not given.
    @raise IllInput Missing or illegal solver inputs.
    @raise LineSearchNonConvergence Line search could not find a suitable iterate.
    @raise MaxIterationsReached The maximum number of nonlinear iterations was reached.
    @raise MaxNewtonStepExceeded Five consecutive steps satisfied a scaled step length test.
    @raise LineSearchBetaConditionFailure  Line search could not satisfy the beta-condition.
    @raise LinearSolverNoRecovery The {!Spils.prec_solve_fn} callback raised {!Sundials.RecoverableFailure} but the preconditioner is already current.
    @raise LinearSolverInitFailure Linear solver initialization failed.
    @raise LinearSetupFailure Linear solver setup failed unrecoverably.
    @raise LinearSolveFailure Linear solver solution failed unrecoverably.
    @raise SystemFunctionFailure The {!sysfn} callback failed unrecoverably.
    @raise FirstSystemFunctionFailure The {!sysfn} callback raised {!Sundials.RecoverableFailure} when first called.
    @raise RepeatedSystemFunctionFailure  The {!sysfn} callback raised {!Sundials.RecoverableFailure} repeatedly. *)
val solve :
    ('d, 'k) session
    -> ('d, 'k) Nvector.t
    -> strategy
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) Nvector.t
    -> result

(** {2:set Modifying the solver (optional input functions)} *)

(** Specifies that an initial call to the preconditioner setup function
    should {i not} be made. This feature is useful when solving a sequence of
    problems where the final preconditioner values of one problem become the
    initial values for the next problem.

    @kinsol KINSetNoInitSetup *)
val set_no_init_setup : ('d, 'k) session -> unit

(** Specifies that an initial call to the preconditioner setup function
    should be made (the default).

    @kinsol KINSetNoInitSetup *)
val set_init_setup : ('d, 'k) session -> unit

(** Disables the nonlinear residual monitoring scheme that controls Jacobian
    updating. It only has an effect for the Dense and Band solvers.

    @kinsol KINSetNoResMon *)
val set_no_res_mon : 'k serial_session -> unit

(** Enables the nonlinear residual monitoring scheme that controls Jacobian
    updating. It only has an effect for the Dense and Band solvers.

    @kinsol KINSetNoResMon *)
val set_res_mon : 'k serial_session -> unit

(** Specifies the maximum number of nonlinear iterations between calls to the
    preconditioner setup function. Pass 0 to set the default (10).

    @kinsol KINSetMaxSetupCalls *)
val set_max_setup_calls : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear iterations between checks by the
    residual monitoring algorithm. Pass 0 to set the default (5). It only
    affects the Dense and Band solvers.

    @kinsol KINSetMaxSubSetupCalls *)
val set_max_sub_setup_calls : 'k serial_session -> int -> unit

(** The parameters {i gamma} and {i alpha} in the formula for the Eisenstat and
    Walker Choice 2 for {i eta}. Set either to [None] to specify its default
    value. The legal values are
    {% $0 < \mathtt{egamma} \leq 1.0 \wedge 1 < \mathtt{ealpha} \leq 2.0$%}.

    @kinsol <Mathematics_link.html#stopping-criteria-for-iterative-linear-solvers> Stopping criteria for iterative linear solvers *)
type eta_params = {
  egamma : float option; (** default = 0.9 *)
  ealpha : float option; (** default = 2.0 *)
}

(** The {i eta} parameter in the stopping criteria for the linear system solver.

    @kinsol <Mathematics_link.html#stopping-criteria-for-iterative-linear-solvers> Stopping criteria for iterative linear solvers *)
type eta_choice =
  | EtaChoice1                   (** Eisenstat and Walker Choice 1 *)
  | EtaChoice2 of eta_params     (** Eisenstat and Walker Choice 2 *)
  | EtaConstant of float option  (** Constant (default = 0.1) *)

(** Specifies the method for computing the value of the {i eta} coefficient used
    in the calculation of the linear solver convergence tolerance.

    @kinsol KINSetEtaForm
    @kinsol KINSetEtaConstValue
    @kinsol KINSetEtaParams *)
val set_eta_choice : ('d, 'k) session -> eta_choice -> unit

(** Specifies the constant value of {i omega} when using residual monitoring.
    Pass 0.0 to specify the default value (0.9). The legal values are
    {% $0 < \mathtt{omega} < 1.0 $%}.

    @kinsol KINSetResMonConstValue *)
val set_res_mon_const_value : ('d, 'k) session -> float -> unit

(** Specifies the minimum and maximum values in the formula for {i omega}.
    The legal values are
    {% $0 < \mathtt{omegamin} < \mathtt{omegamax} < 1.0$%}.

    @kinsol KINSetResMonParams
    @kinsol <Mathematics_link.html#residual-monitoring-for-modified-newton-method> Residual monitoring for Modified Newton method *)
val set_res_mon_params : ('d, 'k) session
                         -> ?omegamin:float
                         -> ?omegamax:float
                         -> unit
                         -> unit

(** Specifies that the scaled linear residual tolerance ({i epsilon})
    is not bounded from below.

    @kinsol KINSetNoMinEps
    @kinsol KINSetFuncNormTol *)
val set_no_min_eps : ('d, 'k) session -> unit

(** Specifies that the scaled linear residual tolerance ({i epsilon})
    is bounded from below. That is, the positive minimum value
    {% $0.01\mathtt{fnormtol}$%} is applied to {i epsilon}.

    @kinsol KINSetNoMinEps
    @kinsol KINSetFuncNormTol *)
val set_min_eps : ('d, 'k) session -> unit

(** Specifies the maximum allowable scaled length of the Newton step. Pass
    0.0 to specify the default value {% $1000\lVert u_0 \rVert_{D_u}$%},
    otherwise the given value must be greater than zero.

    @kinsol KINSetMaxNewtonStep *)
val set_max_newton_step : ('d, 'k) session -> float -> unit

(** Specifies the maximum number of beta-condition failures in the
    line search algorithm. Pass 0.0 to specify the default (10).

    @kinsol KINSetMaxBetaFails *)
val set_max_beta_fails : ('d, 'k) session -> float -> unit

(** Specifies the relative error in computing $F(u)$, which is used in the
    difference quotient approximation of the Jacobian-vector product. Pass
    0.0 to specify the default value
    ({% $\sqrt{\mathtt{unit\_roundoff}}$%}).

    @kinsol KINSetRelErrFunc *)
val set_rel_err_func : ('d, 'k) session -> float -> unit

(** Specifies the stopping tolerance on the scaled maximum norm.
    It must be greater than zero. Pass 0.0 to specify the default
    value ({% $\mathtt{unit\_roundoff}^\frac{1}{3}$%}).

    @kinsol KINSetFuncNormTol *)
val set_func_norm_tol : ('d, 'k) session -> float -> unit

(** Specifies the stopping tolerance on the minimum scaled step length, which
    must be greater than zero. Pass 0.0 to specify the default
    value ({% $\mathtt{unit\_roundoff}^\frac{1}{3}$%}).

    @kinsol KINSetScaledStepTol *)
val set_scaled_step_tol : ('d, 'k) session -> float -> unit

(** Specifies a vector defining inequality constraints for each
    component of the solution vector [u].  See {!Sundials.Constraint}.

    @kinsol KINSetConstraints *)
val set_constraints : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Changes the system function. Allows solutions of several problems of the
    same size but with different functions.

    @kinsol KINSetSysFunc
    @kinsol KINSysFn *)
val set_sys_func : ('d, 'k) session -> ('d -> 'd -> unit) -> unit

(** {3:info Logging and error handling} *)

(** Configure the default error handler to write messages to a file.
    By default it writes to Logfile.stderr.

    @kinsol KINSetErrFile *)
val set_error_file : ('d, 'k) session -> Logfile.t -> unit

(** Specifies a custom function for handling error messages.
    The handler must not fail: any exceptions are trapped and discarded.

    @kinsol KINSetErrHandlerFn
    @kinsol KINErrHandlerFn *)
val set_err_handler_fn
  : ('d, 'k) session -> (Util.error_details -> unit) -> unit

(** Restores the default error handling function.

    @kinsol KINSetErrHandlerFn *)
val clear_err_handler_fn : ('d, 'k) session -> unit

(** Increasing levels of verbosity for informational messages. *)
type print_level =
  | NoInformation     (** No information displayed. {cconst 0} *)
  | ShowScaledNorms   (** At each nonlinear iteration, display the scaled
                          Euclidean norm of the system function at the
                          current iterate, the scaled norm of the Newton step
                          (if no globalization strategy is used), and the number
                          of function evaluations performed so far.
                          {cconst 1} *)
  | ShowScaledDFNorm  (** Additionally display {% $\lVert F(u) \rVert_{DF}$%}
                          if no globalization strategy is used, and
                          {% $\lVert F(u)\rVert_{DF,\infty}$%}, otherwise.
                          {cconst 2} *)
  | ShowGlobalValues  (** Additionally display the values used by the global
                          strategy and statistical information for the linear
                          solver. {cconst 3} *)

(** Sets the level of verbosity of informational messages.

    @kinsol KINSetPrintLevel *)
val set_print_level : ('d, 'k) session -> print_level -> unit

(** Write informational (non-error) messages to the given file.
    By default they are written to Logfile.stdout.
    The optional argument is a convenience for invoking {!set_print_level}.

    @kinsol KINSetInfoFile *)
val set_info_file
      : ('d, 'k) session -> ?print_level:print_level -> Logfile.t -> unit

(** Specifies a custom function for handling informational (non-error) messages.
    The [error_code] field of {{!Sundials.Util.error_details}Util.error_details}
    is [0] for such messages.
    The handler must not fail: any exceptions are trapped and discarded.

    @kinsol KINSetInfoHandlerFn
    @kinsol KINInfoHandlerFn *)
val set_info_handler_fn
  : ('d, 'k) session -> (Util.error_details -> unit) -> unit

(** Restores the default information handling function.

    @kinsol KINSetErrHandlerFn *)
val clear_info_handler_fn : ('d, 'k) session -> unit

(** Specifies whether fixed-point iteration should return the newest
    iteration or the iteration consistent with the last function
    evaluation. The default values is false.

    @kinsol KINSetReturnNewest
    @since 5.8.0 *)
val set_return_newest : ('d, 'k) session -> bool -> unit

(** Sets the damping parameter for the fixed point or Picard iteration.

    To applying damping, the given [beta] value must be greater than 0
    and less than 1. Damping is disabled if [beta >= 1]. The default
    value is 1.

    @kinsol KINSetDamping
    @since 5.8.0 *)
val set_damping : ('d, 'k) session -> float -> unit

(** Set the Anderson acceleration damping parameter.

    To applying damping, the given [beta] value must be greater than 0
    and less than 1. Damping is disabled if [beta >= 1]. The default
    value is 1.

    @kinsol KINSetDampingAA
    @since 5.1.0 *)
val set_damping_aa : ('d, 'k) session -> float -> unit

(** Sets the number of iterations to delay the start of Anderson acceleration.
    The default value is 0.

    @kinsol KINSetDelayAA
    @since 5.8.0 *)
val set_delay_aa : ('d, 'k) session -> float -> unit

(** {2:get Querying the solver (optional output functions)} *)

(** Returns the sizes of the real and integer workspaces.

    @kinsol KINGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : ('d, 'k) session -> int * int

(** Returns the number of evaluations of the system function.

    @kinsol KINGetNumFuncEvals *)
val get_num_func_evals : ('d, 'k) session -> int

(** Returns the cumulative number of nonlinear iterations.

    @kinsol KINGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('d, 'k) session -> int

(** Returns the number of beta-condition failures.

    @kinsol KINGetNumBetaCondFails *)
val get_num_beta_cond_fails : ('d, 'k) session -> int

(** Returns the number of backtrack operations (step length adjustments)
    performed by the line search algorithm.

    @kinsol KINGetNumBacktrackOps *)
val get_num_backtrack_ops : ('d, 'k) session -> int

(** Returns the scaled Euclidiean {i l2} norm of the nonlinear system function
    $F(u)$ evaluated at the current iterate.

    @kinsol KINGetFuncNorm *)
val get_func_norm : ('d, 'k) session -> float

(** Returns the scaled Euclidiean {i l2} norm of the step used during the
    previous iteration.

    @kinsol KINGetStepLength *)
val get_step_length : ('d, 'k) session -> float

(** {2:exceptions Exceptions} *)

(** An input parameter was invalid.

    @kinsol <Constants_link.html> KIN_ILL_INPUT *)
exception IllInput

(** Line search could not find an iterate sufficiently distinct
    from the current one, or an iterate satisfying the sufficient decrease
    condition. The latter could mean that the current iterate is “close” to an
    approximate solution, but that the difference approximation of the
    matrix-vector product is inaccurate, or that [scsteptol]
    ({!set_scaled_step_tol}) is too large.

    @kinsol <Constants_link.html> KIN_LINESEARCH_NONCONV *)
exception LineSearchNonConvergence

(** The maximum number of nonlinear iterations has been reached.

    @kinsol <Constants_link.html> KIN_MAXITER_REACHED *)
exception MaxIterationsReached

(** Five consecutive steps exceed the maximum newton step.
    That is, the five steps satisfy the inequality
    {% $\\|D_u p\\|_{L2} > 0.99 \mathtt{mxnewtstep}$%},
    where $p$ denotes the current step and [mxnewtstep] is a scalar
    upper bound on the scaled step length (see {!set_max_newton_step}).
    It could be that {% $\\| D_F F(u)\\|_{L2}$%} is bounded from above by
    a positive value or that [mxnewtstep] is too small.

    @kinsol <Constants_link.html> KIN_MXNEWT_5X_EXCEEDED *)
exception MaxNewtonStepExceeded

(** The line search algorithm could not satisfy the “beta-condition” for
    [mxnbcf + 1] nonlinear iterations. The failures need not occur in
    consecutive iterations. They may indicate that the algorithm is making
    poor progress.

    @kinsol <Constants_link.html> KIN_LINESEARCH_BCFAIL *)
exception LineSearchBetaConditionFailure

(** The {!Spils.prec_solve_fn} callback raised {!Sundials.RecoverableFailure}
    but the preconditioner is already current.

    @kinsol <Constants_link.html> KIN_LINSOLV_NO_RECOVERY *)
exception LinearSolverNoRecovery

(** Linear solver initialization failed.

    @kinsol <Constants_link.html> KIN_LINIT_FAIL *)
exception LinearSolverInitFailure

(** The {!Spils.prec_setup_fn} callback failed unrecoverably.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.PSetFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @kinsol KINGetLastLinFlag
    @kinsol <Constants_link.html> KIN_LSETUP_FAIL *)
exception LinearSetupFailure of exn option

(** Either {!Spils.prec_solve_fn} failed unrecoverably or the linear solver
    encountered an error condition.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.ATimesFailure},
    {!Sundials_LinearSolver.PSolveFailure},
    {!Sundials_LinearSolver.GSFailure},
    {!Sundials_LinearSolver.QRSolFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @kinsol KINGetLastLinFlag
    @kinsol <Constants_link.html> KIN_LSOLVE_FAIL *)
exception LinearSolveFailure of exn option

(** The {!sysfn} callback failed unrecoverably.

    @kinsol <Constants_link.html> KIN_SYSFUNC_FAIL *)
exception SystemFunctionFailure

(** The {!sysfn} callback raised {!Sundials.RecoverableFailure} when first called.

    @kinsol <Constants_link.html> KIN_FIRST_SYSFUNC_FAIL *)
exception FirstSystemFunctionFailure

(** The {!sysfn} callback raised {!Sundials.RecoverableFailure} repeatedly.
    No recovery is possible.

    @kinsol <Constants_link.html> KIN_REPTD_SYSFUNC_ERR *)
exception RepeatedSystemFunctionFailure

(** A linear solver is required but was not specified. *)
exception MissingLinearSolver

(** A fused vector operation failed.

    @kinsol <Constants_link.html> KIN_VECTOROP_ERR *)
exception VectorOpErr

