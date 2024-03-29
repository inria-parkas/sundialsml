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
(*               User Documentation for CVODES v2.7.0                  *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Sensitivity analysis (forward and adjoint) and quadrature equations
    for CVODE.

    These submodules extend basic state functions, {% $\dot{y} = f(t, y)$%},
    with parameters $p$ and quadrature variables $y_Q$. The former exposes
    a new dependency, {% $\dot{y} = f(t, y, p)$%}, relative to which the
    sensitivity of a solution can be approximated. The latter effectively
    suppresses a dependency ($y_Q$ is disjoint from $y$),
    {% $\dot{y}_Q = f_Q(t, y, p)$%}, to reduce computation costs.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** Alias for Cvode sessions. The Cvodes submodules all work by ‘extending’
    an existing session created with {!Cvode.init}. *)
type ('data, 'kind) session = ('data, 'kind) Cvode.session

(** Integration of pure quadrature equations.

    Adds a vector $y_Q$ of $N_Q$ quadrature variables defined by
    {% $\frac{\mathrm{d} y_Q}{\mathrm{d}t} = f_Q(t, y, p)$%}. The values of
    these variables are calculated more efficiently since they are excluded
    from the nonlinear solution stage.

    @cvodes <Usage/SIM.html#integration-of-pure-quadrature-equations> Integration of pure quadrature equations
    @cvodes <Mathematics_link.html#pure-quadrature-integration> Pure Quadrature Integration *)
module Quadrature : sig (* {{{ *)
  (** {2:init Initialization and access} *)

  (** Functions defining quadrature variables. They are passed three
      arguments:
      - [t], the value of the independent variable, i.e., the simulation time,
      - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
      - [yQ'], a vector for storing the computed value of
               {% $\dot{y}_Q = f_Q(t, y)$%}.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning [y] and [yQ'] should not be accessed after the function
               returns.}

      @cvodes_quad CVQuadRhsFn *)
  type 'd quadrhsfn = float -> 'd -> 'd -> unit

  (** Activates the integration of quadrature equations. The vector
      gives the initial value of $y_Q$.

      @cvodes_quad CVodeQuadInit *)
  val init : ('d, 'k) Cvode.session -> 'd quadrhsfn
                -> ('d, 'k) Nvector.t -> unit

  (** Reinitializes the integration of quadrature equations. The vector
      gives a new value for $y_Q$.

      @cvodes_quad CVodeQuadReInit *)
  val reinit : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> unit

  (** Returns the quadrature solutions and time reached after a successful
      solver step. The given vector is filled with values calculated during
      either {!Cvode.solve_normal} or {!Cvode.solve_one_step} and the
      value of the independent variable is returned.

      @cvodes_quad CVodeGetQuad *)
  val get : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> float

  (** Returns the interpolated solution or derivatives of quadrature
      variables.

      [get_dky s dkyq t k] computes the [k]th derivative at time [t], i.e.,
      {% $\frac{d^\mathtt{k}y_Q(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
      and stores it in [dkyq]. The arguments must satisfy {% $t_n - h_u \leq
      \mathtt{t} \leq t_n$%}—where $t_n$ denotes {!Cvode.get_current_time}
      and $h_u$ denotes {!Cvode.get_last_step},— and
      {% $0 \leq \mathtt{k} \leq q_u$%}—where
      $q_u$ denotes {!Cvode.get_last_order}.

      This function may only be called after a successful return from either
      {!Cvode.solve_normal} or {!Cvode.solve_one_step}.

      @cvodes_quad CVodeGetQuadDky
      @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
      @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
  val get_dky : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> float -> int -> unit

  (** {2:tols Tolerances} *)

  (** Tolerances for calculating quadrature variables. *)
  type ('d, 'k) tolerance =
      NoStepSizeControl
      (** Do not use quadrature variables for step-size control (default). *)
    | SStolerances of float * float
      (** [(rel, abs)] : scalar relative and absolute tolerances. *)
    | SVtolerances of float * ('d, 'k) Nvector.t
      (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

  (** Specify how to use quadrature variables in step size control.

      @cvodes_quad CVodeSetQuadErrCon
      @cvodes_quad CVodeQuadSStolerances
      @cvodes_quad CVodeQuadSVtolerances *)
  val set_tolerances : ('d, 'k) Cvode.session -> ('d, 'k) tolerance -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the number of calls to the quadrature function.

      @cvodes_quad CVodeGetQuadNumRhsEvals *)
  val get_num_rhs_evals       : ('d, 'k) Cvode.session -> int

  (** Returns the number of local error test failures that have occurred
      due to quadrature variables.

      @cvodes_quad CVodeGetQuadNumErrTestFails *)
  val get_num_err_test_fails  : ('d, 'k) Cvode.session -> int

  (** Returns the quadrature error weights at the current time.

      @cvodes_quad CVodeGetQuadErrWeights *)
  val get_err_weights : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> unit

  (** Returns quadrature-related statistics. These are the
      number of calls to the quadrature function ([nfqevals]) and the number
      of error test failures due to quadrature variables ([nqetfails]).

      @cvodes_quad CVodeGetQuadStats
      @return ([nfqevals], [nqetfails]) *)
  val get_stats : ('d, 'k) Cvode.session -> int * int

  (** {2:exceptions Exceptions} *)

  (** Quadrature integration was not initialized.

      @cvodes <Constants_link.html> CV_NO_QUAD *)
  exception QuadNotInitialized

  (** The quadrature function failed in an unrecoverable manner.

      @cvodes <Constants_link.html> CV_QRHSFUNC_FAIL *)
  exception QuadRhsFuncFailure

  (** The quadrature function failed at the first call.

      @cvodes <Constants_link.html> CV_FIRST_QRHSFUNC_ERR *)
  exception FirstQuadRhsFuncFailure

  (** Convergence test failures occurred too many times due to repeated
      recoverable errors in the quadrature function. Also raised if the
      quadrature function had repeated recoverable errors during the
      estimation of an initial step size (if quadrature variables are
      included in error tests).

      @cvodes <Constants_link.html> CV_REPTD_QRHSFUNC_ERR *)
  exception RepeatedQuadRhsFuncFailure

  (** The quadrature function had a recoverable error, but no
      recovery was possible. This failure mode is rare, as it can occur only
      if the quadrature function fails recoverably after an error test failed
      while at order one.

      @cvodes <Constants_link.html> CV_UNREC_QRHSFUNC_ERR *)
  exception UnrecoverableQuadRhsFuncFailure
end (* }}} *)

(** (Forward) Sensitivity analysis of ODEs with respect to their parameters.

    Formalizes the dependency of a set of ODEs on $N_p$ parameters $p$ and
    calculates the sensitivities $s$ of the solution $y$ to a subset of
    {% $N_s \leq N_p$%} of those parameters;
    {% $s(t) = \frac{\partial y(t)}{\partial p}$%}.

    The {i solution sensitivity} with respect to a single parameter $p_i$
    satisfies the {i (forward) sensitivity equations}
    {% $\dot{s}_i = \frac{\partial f}{\partial y}s_i + \frac{\partial f}{\partial p_i}$%}
    and {% $s_i(t_0) = \frac{\partial y_0(p)}{\partial p_i}$%}, where
    $f$ and $y$ are from the $N$ equations of the original model.

    This documented interface is structured as follows.
    {ol
      {- {{:#init}Initialization} (including {{!Quadrature}Quadrature equations})}
      {- {{:#sensout}Output functions}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#exceptions}Exceptions}}}

    @cvodes <Usage/FSA.html#using-cvodes-for-forward-sensitivity-analysis> Using CVODES for Forward Sensitivity Analysis
    @cvodes <Usage/FSA.html#a-skeleton-of-the-users-main-program> Enhanced skeleton for sensitivity analysis
    @cvodes <Mathematics_link.html#forward-sensitivity-analysis> Forward sensitivity analysis *)
module Sensitivity : sig (* {{{ *)
  (** {2:init Initialization} *)

  (** Common arguments to {!sensrhsfn1} and {!sensrhsfn_all}.  *)
  type 'd sensrhsfn_args =
    {
      t : float;
      (** The value of the independent variable. *)

      y : 'd;
      (** The vector of dependent-variable values $y(t)$. *)

      y' : 'd;
      (** The value of the right-hand side of the state
          equations {% $\dot{y} = f(t, y, p)$%}. *)

      tmp : 'd Cvode.double;
      (** Temporary storage vectors. *)
    }

  (** Sensitivity functions that calculate the right-hand sides of
      all sensitivity equations.  They are passed the arguments:

      - [args], the current values of state variables (see {!sensrhsfn_args}),
      - [s], an array of vectors holding the current values of sensitivity
             variables, and,
      - [s'], an array of vectors to be filled with the derivatives
              of the sensitivity variables.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning The vectors in the function arguments should not
               be accessed after the function returns.}

      @cvodes_sens CVSensRhsFn *)
  type 'd sensrhsfn_all = 'd sensrhsfn_args
                        -> 'd array
                        -> 'd array
                        -> unit

  (** Sensitivity functions that calculate the right-hand side of a
      single sensitivity equation.  They are passed the arguments:

      - [i], the index of the sensitivity equation to compute,
      - [args], the current values of state variables
                (see {!sensrhsfn_args}),
      - [s], a vector holding the current value of the {i i}th
             sensitivity variable, and
      - [s'], a vector to be filled with the current value of the {i i}th
              sensitivity variable's derivative.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning The vectors in the function arguments should not
               be accessed after the function returns.}

      @cvodes_sens CVSensRhs1Fn *)
  and 'd sensrhsfn1 = int
                    -> 'd sensrhsfn_args
                    -> 'd
                    -> 'd
                    -> unit

  (** Specify a sensitivity function. *)
  type 'd sensrhsfn =
      AllAtOnce of 'd sensrhsfn_all option
      (** Calculate sensitivity functions all at once. The argument [None]
          specifies an internal difference quotient routine.

          @cvodes_sens CVodeSensInit
          @cvodes_sens CVSensRhsFn *)

    | OneByOne of 'd sensrhsfn1 option
      (** Calculate sensitivity functions one parameter at a time. The
          argument [None] specifies an internal difference quotient routine.

          @cvodes_sens CVodeSensInit1
          @cvodes_sens CVSensRhs1Fn *)

  (** Specifies a sensitivity solution method. The method specification may
      optionally include the nonlinear solver to be used for corrections.

      @cvodes_sens CVodeSensInit
      @cvodes_sens CVodeSensInit1
      @cvodes_sens CVodeSetNonlinearSolverSensSim
      @cvodes_sens CVodeSetNonlinearSolverSensStg
      @cvodes_sens CVodeSetNonlinearSolverSensStg1 *)
  type ('d, 'k) sens_method =
      Simultaneous of
        (('d, 'k, ('d, 'k) Cvode.session, [`Sens]) Sundials.NonlinearSolver.t) option
      (** Correct state and sensitivity variables at the same time.
          {cconst CV_SIMULTANEOUS} *)
    | Staggered of
        (('d, 'k, ('d, 'k) Cvode.session, [`Sens]) Sundials.NonlinearSolver.t) option
      (** The correction step for the sensitivity variables takes place at the
          same time for all sensitivity equations, but only after the
          correction of the state variables has converged and the state
          variables have passed the local error test.
          {cconst CV_STAGGERED} *)
    | Staggered1 of
        (('d, 'k, ('d, 'k) Cvode.session, [`Nvec]) Sundials.NonlinearSolver.t) option
      (** All corrections are done sequentially, first for the state variables
          and then for the sensitivity variables, one parameter at a time. If
          the sensitivity variables are not included in the error control,
          this approach is equivalent to [Staggered]. Note that this
          approach can only be used if the user-provided sensitivity
          right-hand side function is of type [OneByOne].
          {cconst CV_STAGGERED1} *)

  (** Specifies problem parameter information for sensitivity
      calculations.  Which fields are required varies as follows:
      - [pvals] is mandatory if {!Sensitivity.init} and/or
        {!Sensitivity.Quadrature.init} is not given a
        user-defined callback.  Otherwise, it is ignored.
      - [plist] is optional if [pvals] is given.  Otherwise, it is
        ignored.
      - [pbar] is always meaningful and optional.

      @cvodes_sens CVodeSetSensParams *)
  type sens_params = {
      pvals  : RealArray.t option;
      (** An array of $N_p$ parameters $p$ in $f(t, y, p)$.  If
          specified, this array is updated by the solver to pass
          perturbed parameter values to the original right-hand side
          and root functions.  Those functions must (re-)read
          parameter values from this array on every invocation. *)
      pbar   : RealArray.t option;
      (** An array of $N_s$ positive scaling factors. These are needed
          when estimating tolerances or using the internal difference
          quotient function. *)
      plist  : int array option;
      (** An array of $N_s$ ($< N_p$) indices into [pvals]
          specifying the parameters to analyze. If not specified,
          the first $N_s$ parameters are analyzed. *)
    }

  (** Empty [pvals], [plist], and [pbar] fields. *)
  val no_sens_params : sens_params

  (** Tolerances for calculating sensitivities. *)
  type ('d, 'k) tolerance =
      SStolerances of float * RealArray.t
      (** [(rel, abs)] : scalar relative and absolute tolerances. *)
    | SVtolerances of float * ('d, 'k) Nvector.t array
      (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
    | EEtolerances
      (** Calculate the integration tolerances for sensitivities
          from those for state variables and the scaling factors
          in {{!sens_params}pbar}. *)

  (** Activates the calculation of forward sensitivities. The call
      [init s tol sens_method ~sens_params fs ys0], has as arguments:
      - [s], a session created with {!Cvode.init},
      - [tol], the tolerances desired,
      - [sens_method], the solution method,
      - [sens_params], the parameter information (see {!sens_params}),
      - [fs], the sensitivity function, and,
      - [ys0], initial values of the sensitivities for each parameter.

      If [fs] is not given, an internal difference quotient routine
      is used.  In that case, or if the internal difference quotient
      routine will be specified in a subsequent call to
      {!Sensitivity.Quadrature.init}, then [sens_params] must be
      given with [pvals] set to non-[None].

      @cvodes_sens CVodeSensInit
      @cvodes_sens CVodeSensInit1
      @cvodes_sens CVodeSetSensParams
      @cvodes_sens CVodeSensSStolerances
      @cvodes_sens CVodeSensSVtolerances
      @cvodes_sens CVodeSensEEtolerances *)
  val init : ('d, 'k) Cvode.session
             -> ('d, 'k) tolerance
             -> ('d, 'k) sens_method
             -> ?sens_params:sens_params
             -> 'd sensrhsfn
             -> ('d, 'k) Nvector.t array
             -> unit

  (** Reinitializes the forward sensitivity computation.

      @cvodes_sens CVodeSensReInit *)
  val reinit : ('d, 'k) Cvode.session -> ('d, 'k) sens_method
                    -> ('d, 'k) Nvector.t array -> unit

  (** Deactivates forward sensitivity calculations without deallocating
      memory. Sensitivities can be reactivated with {!reinit}.

      @cvodes_sens CVodeSensToggleOff *)
  val toggle_off : ('d, 'k) Cvode.session -> unit

  (** Deactivates forward sensitivity calculations with memory
      deallocation.

      @cvodes_sens CVodeSensFree *)
  val turn_off : ('d, 'k) Cvode.session -> unit

  (** Support for quadrature sensitivity equations.

      Adds an additional vector {% $s_\mathit{Q}$%} of {% $N_\mathit{Q}$%}
      quadrature sensitivities
      {% $\frac{\mathrm{d} s_\mathit{Q}}{\mathrm{d}t}
              = f_{\mathit{QS}}(t, y, s, \dot{y}_Q, p)$%}.
      This mechanism allows, in particular, the calculation of the
      sensitivities of the ‘pure’ {!Cvodes.Quadrature} variables, $y_Q$.

      @cvodes <Usage/FSA.html#integration-of-quadrature-equations-depending-on-forward-sensitivities> Integration of quadratures equations depending on forward sensitivities *)
  module Quadrature : sig (* {{{ *)
    (** {2:init Initialization} *)

    (** Arguments to {!quadsensrhsfn}. *)
    type 'd quadsensrhsfn_args =
      {
        t : float;
        (** The value of the independent variable. *)

        y : 'd;
        (** The vector of dependent-variable values $y(t)$. *)

        s : 'd array;
        (** The array of sensitivity vectors. *)

        yq' : 'd;
        (** The value of the quadrature-right hand side {% $\dot{y}_Q$%}. *)

        tmp : 'd Cvode.double;
        (** Temporary storage vectors. *)
      }

    (** Functions defining quadrature sensitivities.
        They are passed the arguments:
        - [args], the current values of state and sensitivity variables,
                  and,
        - [sq'], an array of vectors for storing the computed values of
                {% $\dot{s}_\mathit{Q} = f_\mathit{QS}(t, y, s, \dot{y}_Q)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure}
        exception indicates a recoverable error. Any other exception is
        treated as an unrecoverable error.

        {warning The vectors in the function arguments should not
                 be accessed after the function returns.}

       @cvodes_sens CVodeQuadSensRhsFn *)
    type 'd quadsensrhsfn = 'd quadsensrhsfn_args -> 'd array -> unit

    (** Activate the integration of quadrature sensitivities.  The
        arguments are:
        - [~fqs], a function that computes the right-hand sides of the
          quadrature sensitivities, and
        - [q0], an array of vectors specifying initial values for the
          quadrature sensitivities.

        If [~fqs] is not given, an internal difference quotient
        routine is used.  In that case, {!Sensitivity.init} must
        have been invoked with a [sens_params] whose [pvals] is
        non-[None].

        @cvodes_sens CVodeQuadSensInit
        @raise QuadNotInitialized {!Quadrature.init} has not been called. *)
    val init : ('d, 'k) Cvode.session -> ?fqs:'d quadsensrhsfn
             -> ('d, 'k) Nvector.t array -> unit

    (** Reinitializes the quadrature sensitivity integration.

        @cvodes_sens CVodeQuadSensReInit *)
    val reinit : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array -> unit

    (** {2:tols Tolerance specification} *)

    (** Tolerances for calculating quadrature sensitivities. *)
    type ('d, 'k) tolerance =
        NoStepSizeControl
        (** Quadrature variables are not used for step-size control
            (the default). *)
      | SStolerances of float * RealArray.t
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('d, 'k) Nvector.t array
        (** [(rel, abs)] : scalar relative and vector absolute
            tolerances. *)
      | EEtolerances
        (** Calculate the integration tolerances for the
            quadrature sensitivities from those provided for
            the pure quadrature variables. *)

    (** Specify how to use quadrature sensitivities in step size control.

        @cvodes_sens CVodeSetQuadSensErrCon
        @cvodes_sens CVodeQuadSensSStolerances
        @cvodes_sens CVodeQuadSensSVtolerances
        @cvodes_sens CVodeQuadSensEEtolerances *)
    val set_tolerances : ('d, 'k) Cvode.session
                          -> ('d, 'k) tolerance -> unit

    (** {2:quadout Output functions} *)

    (** Returns the quadrature sensitivity solutions and time reached
        after a successful solver step. The given vectors are filled with
        values calculated during either {!Cvode.solve_normal} or
        {!Cvode.solve_one_step} and the value of the independent variable
        is returned.

        @cvodes_sens CVodeGetQuadSens *)
    val get : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array -> float

    (** Returns a single quadrature sensitivity vector after a successful
        solver step. Like {!get}, but the argument [i] specifies a specific
        vector to return.

        @cvodes_sens CVodeGetQuadSens1
        @raise BadIS The index is not in the allowed range. *)
    val get1 : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> int -> float

    (** Returns the interpolated solution or derivatives of the quadrature
        sensitivity solution.

        [get_dky s dksq t k] computes the [k]th derivative at time [t],
        i.e.,
        {% $\frac{d^\mathtt{k}s_\mathit{Q}(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
        and stores it in [dksq]. The arguments must satisfy {% $t_n - h_u
        \leq \mathtt{t} \leq t_n$%}—where $t_n$ denotes
        {!Cvode.get_current_time} and $h_u$ denotes
        {!Cvode.get_last_step},—and
        {% $0 \leq \mathtt{k} \leq q_u$%}—where $q_u$ denotes
        {!Cvode.get_last_order}.

        @cvodes_sens CVodeGetQuadSensDky
        @raise BadIS The index is not in the allowed range.
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array
                    -> float -> int -> unit

    (** Returns the interpolated solution or derivatives of a single
        quadrature sensitivity solution vector.
        [get_dky s dksq t k i] is like
        {!get_dky} but restricted to the [i]th sensitivity solution vector.

        @cvodes_sens CVodeGetQuadSensDky1
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky1 : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t
                     -> float -> int -> int -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the number of calls to the quadrature right-hand side
        function.

        @cvodes_sens CVodeGetQuadSensNumRhsEvals *)
    val get_num_rhs_evals       : ('d, 'k) Cvode.session -> int

    (** Returns the number of local error test failures due to quadrature
        variables.

        @cvodes_sens CVodeGetQuadSensNumErrTestFails *)
    val get_num_err_test_fails  : ('d, 'k) Cvode.session -> int

    (** Returns the quadrature error weights at the current time.

        @cvodes_sens CVodeGetQuadSensErrWeights *)
    val get_err_weights
          : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array -> unit

    (** Returns quadrature-related statistics. These are the
        number of calls to the quadrature function ([nfqevals]) and the
        number of error test failures due to quadrature variables
        ([nqetfails]).

        @cvodes_sens CVodeGetQuadSensStats
        @return ([nfqevals], [nqetfails]) *)
    val get_stats : ('d, 'k) Cvode.session -> int * int

    (** {2:exceptions Exceptions} *)

    (** Quadrature integration was not initialized.

        @cvodes <Constants_link.html> CV_NO_QUAD_SENS *)
    exception QuadSensNotInitialized

    (** The sensitivity quadrature function failed in an unrecoverable
        manner.

        @cvodes <Constants_link.html> CV_QSRHSFUNC_FAIL *)
    exception QuadSensRhsFuncFailure

    (** The sensitivity quadrature function failed at the first call.

        @cvodes <Constants_link.html> CV_FIRST_QSRHSFUNC_ERR *)
    exception FirstQuadSensRhsFuncFailure

    (** Convergence test failures occurred too many times due to repeated
        recoverable errors in the quadrature function. Also raised if the
        sensitivity quadrature function had repeated recoverable errors
        during the estimation of an initial step size (if quadrature
        variables are included in error tests).

        @cvodes <Constants_link.html> CV_REPTD_QSRHSFUNC_ERR *)
    exception RepeatedQuadSensRhsFuncFailure

    (** The sensitivity quadrature function had a recoverable error, but no
        recovery was possible. This failure mode is rare, as it can occur
        only if the quadrature function fails recoverably after an error
        test failed while at order one.

        @cvodes <Constants_link.html> CV_UNREC_QSRHSFUNC_ERR *)
    exception UnrecoverableQuadSensRhsFuncFailure
  end (* }}} *)

  (** {2:sensout Output functions} *)

  (** Returns the sensitivity solution vectors after a successful solver
      step. The given vectors are filled with values calculated during
      either {!Cvode.solve_normal} or {!Cvode.solve_one_step} and the
      value of the independent variable is returned.

      @cvodes_sens CVodeGetSens *)
  val get : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array -> float

  (** Returns the interpolated solution or derivatives of the sensitivity
      solution vectors.

      [get_dky s dks t k] computes the [k]th derivative at time [t], i.e.,
      {% $\frac{d^\mathtt{k}s(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
      and stores it in [dks]. The arguments must satisfy {% $t_n - h_u \leq
      \mathtt{t} \leq t_n$%}—where $t_n$ denotes {!Cvode.get_current_time}
      and $h_u$ denotes {!Cvode.get_last_step},— and
      {% $0 \leq \mathtt{k} \leq q_u$%}—where
      $q_u$ denotes {!Cvode.get_last_order}.

      This function may only be called after a successful return from either
      {!Cvode.solve_normal} or {!Cvode.solve_one_step}.

      @cvodes_sens CVodeGetSensDky
      @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
      @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
  val get_dky : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array
                  -> float -> int -> unit

  (** Returns a single sensitivity solution vector after a successful solver
      step. The given vector is filled with values calculated for the [i]th
      sensitivity vector during either {!Cvode.solve_normal} or
      {!Cvode.solve_one_step} and the value of the independent variable is
      returned.

      @cvodes_sens CVodeGetSens1
      @raise BadIS The index [i] is not in the allowed range. *)
  val get1 : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> int -> float

  (** Returns the interpolated solution or derivatives of a single
      sensitivity solution vector. [get_dky s dks t k i] is like {!get_dky}
      but restricted to the [i]th sensitivity solution vector.

      @cvodes_sens CVodeGetSensDky1
      @raise BadIS The index [i] is not in the allowed range.
      @raise BadK [k] is not in the range 0, 1, ..., [qlast].
      @raise BadT [t] is not in the allowed range. *)
  val get_dky1 : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t
                   -> float -> int -> int -> unit

  (** {2:set Modifying the solver (optional input functions)} *)

  (** Sets the integration tolerances for sensitivities.

      {b NB}: Unlike the other [set_tolerances] functions in [Cvodes], this
      one does {b not} call {!set_err_con} (which defaults to [false]).

      @cvodes_sens CVodeSensSStolerances
      @cvodes_sens CVodeSensSVtolerances
      @cvodes_sens CVodeSensEEtolerances *)
  val set_tolerances : ('d, 'k) Cvode.session -> ('d, 'k) tolerance -> unit

  (** Sets whether sensitivity variables are used in the error control
      mechanism. The default is [false].

      @cvodes_sens CVodeSetSensErrCon *)
  val set_err_con : ('d, 'k) Cvode.session -> bool -> unit

  (** A difference quotient strategy. See {!set_dq_method}. *)
  type dq_method = DQCentered (** {cconst CV_CENTERED} *)
                 | DQForward  (** {cconst CV_FORWARD} *)

  (** Sets the difference quotient strategy when sensitivity equations
      are computed internally by the solver rather than via callback. A
      method and [dqrhomax] value must be given. The latter determines when
      to switch between simultaneous or separate approximations of the two
      terms in the sensitivity right-hand side.

      @cvodes_sens CVodeSetSensDQMethod *)
  val set_dq_method : ('d, 'k) Cvode.session -> dq_method -> float -> unit

  (** Sets the maximum number of nonlinear solver iterations for
      sensitivity variables permitted per step.

      @cvodes_sens CVodeSetSensMaxNonlinIters *)
  val set_max_nonlin_iters : ('d, 'k) Cvode.session -> int -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the number of calls to the sensitivity function.

      @cvodes_sens CVodeGetSensNumRhsEvals *)
  val get_num_rhs_evals       : ('d, 'k) Cvode.session -> int

  (** Returns the number of calls to the right-hand side function due
      to the internal finite difference approximation of the sensitivity
      equations.

      @cvodes_sens CVodeGetNumRhsEvalsSens *)
  val get_num_rhs_evals_sens  : ('d, 'k) Cvode.session -> int

  (** Returns the number of local error test failures for the sensitivity
      variables that have occurred.

      @cvodes_sens CVodeGetSensNumErrTestFails *)
  val get_num_err_test_fails  : ('d, 'k) Cvode.session -> int

  (** Returns the number of calls made to the linear solver's setup function
      due to forward sensitivity calculations.

      @cvodes_sens CVodeGetSensNumLinSolvSetups *)
  val get_num_lin_solv_setups : ('d, 'k) Cvode.session -> int

  (** Summaries of sensitivity stats. *)
  type sensitivity_stats = {
      num_sens_evals :int;
        (** Number of calls to the sensitivity function. *)
      num_rhs_evals : int;
        (** Number of calls to the right-hand side function to
            calculate sensitivities. *)
      num_err_test_fails : int;
        (** Number of local error test failures for sensitivity variables. *)
      num_lin_solv_setups :int;
        (** Number of setups calls to the linear solver for sensitivity
            calculations. *)
    }

  (** Returns the sensitivity-related statistics as a group.

      @cvodes_sens CVodeGetSensStats *)
  val get_stats : ('d, 'k) Cvode.session -> sensitivity_stats

  (** Returns the sensitivity error weights at the current time.

      @cvodes_sens CVodeGetSensErrWeights
      @cvodes <Mathematics_link.html#equation-cvodes-errwt> IVP solution (W_i) *)
  val get_err_weights : ('d, 'k) Cvode.session
                          -> ('d, 'k) Nvector.t array -> unit

  (** Returns the cumulative number of nonlinear iterations performed for
      sensitivity calculations.

      @cvodes_sens CVodeGetSensNumNonlinSolvIters *)
  val get_num_nonlin_solv_iters : ('d, 'k) Cvode.session -> int

  (** Returns the cumulative number of nonlinear convergence failures
      during sensitivity calculations.

      @cvodes_sens CVodeGetSensNumNonlinSolvConvFails *)
  val get_num_nonlin_solv_conv_fails : ('d, 'k) Cvode.session -> int

  (** Returns both the numbers of nonlinear iterations performed [nniters]
      and nonlinear convergence failures [nncfails] during sensitivity
      calculations.

      @cvodes_sens CVodeGetSensNonlinSolvStats
      @return ([nniters], [nncfails]) *)
  val get_nonlin_solv_stats : ('d, 'k) Cvode.session -> int * int

  (** Returns the cumulative number of nonlinear (functional or Newton)
      iterations for each sensitivity equation separately in the [Staggered1]
      case.

      @cvodes_sens CVodeGetStgrSensNumNonlinSolvIters *)
  val get_num_stgr_nonlin_solv_iters : ('d, 'k) Cvode.session
                                       -> LintArray.t -> unit

  (** Returns the cumulative number of nonlinear convergence failures
      for each sensitivity equation separately in the [Staggered1] case.

      @cvodes_sens CVodeGetStgrSensNumNonlinSolvConvFails *)
  val get_num_stgr_nonlin_solv_conv_fails : ('d, 'k) Cvode.session
                                            -> LintArray.t -> unit

  (** Returns the current sensitivity state vector array. The vectors in the
      returned array provide direct access to the data within the integrator.

      @cvodes_sens CVodeGetCurrentStateSens
      @since 5.0.0 *)
  val get_current_state_sens : ('d, 'k) session -> 'd array

  (** Internal data required to construct the current nonlinear implicit
      system within a nonlinear solver. *)
  type 'd nonlin_system_data = {
    tn     : float;
      (** Independent variable value {% $t_n$ %}. *)
    yspred : 'd array;
      (** Predicted state vectors {% $\mathit{yS}_{i,\mathit{pred}}$ %}
          at {% $t_n$ %} for {% $i = 0,\ldots,N_s-1$ %}. This data must not
          be changed. *)
    ysn    : 'd array;
      (** State vectors {% $\mathit{yS}^n_i$ %}
          for {% $i = 0,\ldots,N_s-1$ %}.
          This data may not be current and may need to be filled. *)
    gamma  : float;
      (** Current value of {% $\gamma$ %}. *)
    rls1   : float;
        (** A scaling factor used to compute {% $\tilde{a}_n =
            \mathtt{rls1}\cdot\mathtt{zns1}$ %}. *)
    zns1   : 'd array;
        (** Vectors used to compute {% $\tilde{a}_n =
            \mathtt{rl1}\cdot\mathtt{zn1}$ %}. *)
  }

  (** Gives direct access to the internal data required to construct the
      current nonlinear system within a nonlinear solver. This
      function should be called inside the nonlinear system function.
      The vectors [ysn] are provided as additional workspace and do not need
      to be filled in. They are only current after an evaluation of the
      nonlinear system function.

      @cvodes_sens CVodeGetNonlinearSystemDataSens
      @since 5.4.0 *)
  val get_nonlin_system_data : ('d, 'k) session -> 'd nonlin_system_data

  (** Computes the current sensitivity vector for all sensitivities using
      the stored prediction and the supplied correction vectors from the
      nonlinear solver. The call
      [compute_state s i yscor ysn] computes
      {% $\mathit{yS}^n = \mathit{yS}_{\mathit{pred}}
                            + \mathit{yS}_{\mathit{cor}}$ %}.

      @cvodes_sens CVodeComputeStateSens
      @since 5.4.0 *)
  val compute_state : ('d, 'k) session
                      -> ('d, 'k) Nvector.t array
                      -> ('d, 'k) Nvector.t array
                      -> unit

  (** Computes the current sensitivity vector for the sensitivity at the given
      index using the stored prediction and the supplied correction vector
      from the nonlinear solver. The call
      [compute_state s i yscor1 ysn1] computes
      {% $\mathit{yS}^n_i = \mathit{yS}_{i,\mathit{pred}}
                            + \mathit{yS}_{i,\mathit{cor}}$ %}.

      @cvodes_sens CVodeComputeStateSens1
      @since 5.4.0 *)
  val compute_state1 : ('d, 'k) session
                       -> int
                       -> ('d, 'k) Nvector.t
                       -> ('d, 'k) Nvector.t
                       -> unit

  (** Returns the index of the current sensitivity solve when using
      the {{!sens_method}Staggered1} method.

      @cvodes_sens CVodeGetCurrentSensSolveIndex
      @since 5.0.0 *)
  val get_current_sens_solve_index : ('d, 'k) session -> int

  (** {2:exceptions Exceptions} *)

  (** Sensitivity analysis was not initialized.

      @cvodes <Constants_link.html> CV_NO_SENS *)
  exception SensNotInitialized

  (** The sensitivity function failed in an unrecoverable manner.

      @cvodes <Constants_link.html> CV_SRHSFUNC_FAIL *)
  exception SensRhsFuncFailure

  (** The sensitivity function had a recoverable error when first called.

      @cvodes <Constants_link.html> CV_FIRST_SRHSFUNC_ERR *)
  exception FirstSensRhsFuncFailure

  (** Too many convergence test failures, or unable to estimate the initial
      step size, due to repeated recoverable errors in the sensitivity
      function.

      @cvodes <Constants_link.html> CV_REPTD_SRHSFUNC_ERR *)
  exception RepeatedSensRhsFuncFailure

  (** The sensitivity function had a recoverable error, but no recovery was
      possible. This error can only occur after an error test failure at
      order one.

    @cvodes <Constants_link.html> CV_UNREC_SRHSFUNC_ERR *)
  exception UnrecoverableSensRhsFuncFailure

  (** The index passed to identify a particular sensitivity is invalid.

      @cvodes <Constants_link.html> CV_BAD_IS *)
  exception BadSensIdentifier
end (* }}} *)

(** (Adjoint) Sensitivity analysis of ODEs with respect to their parameters.

    Provides an alternative to forward sensitivity analysis, which can become
    prohibitively expensive. This technique does not calculate sensitivities,
    but rather gradients with respect to the parameters of a relatively few
    derived functionals of the solution, that is the gradient
    {% $\frac{\mathrm{d}G}{\mathrm{d}p}$%} of
    {% $G(p) = \int_{t_0}^T \! g(t, y, p)\,\mathrm{d}t$%}. The gradients
    are evaluated by first calculating forward and checkpointing certain
    intermediate state values, and then integrating backward to $t_0$.

    This documented interface is structured as follows.
    {ol
      {- {{:#fwd}Forward solution}}
      {- {{:#linear}Linear solvers}}
      {- {{:#bsolve}Backward solutions} (including {{!Quadrature}Quadrature equations})}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#exceptions}Exceptions}}}

    @cvodes <Usage/ADJ.html#using-cvodes-for-adjoint-sensitivity-analysis> Using CVODES for Adjoint Sensitivity Analysis
    @cvodes <Usage/ADJ.html#a-skeleton-of-the-users-main-program> Enhanced Skeleton for Adjoint Sensitivity Analysis
    @cvodes <Mathematics_link.html#adjoint-sensitivity-analysis> Adjoint sensitivity analysis *)
module Adjoint : sig (* {{{ *)
  (** A backward session with the CVODES solver. Multiple backward sessions
      may be associated with a single parent session.

      @cvodes <Usage/ADJ.html#adjoint-sensitivity-allocation-and-deallocation-functions> Adjoint sensitivity allocation and deallocation functions *)
  type ('data, 'kind) bsession = ('data, 'kind) Cvode_impl.AdjointTypes.bsession

  (** Alias for backward sessions based on serial nvectors. *)
  type 'kind serial_bsession = (Nvector_serial.data, 'kind) bsession
                               constraint 'kind = [>Nvector_serial.kind]

  (** {2:fwd Forward solution} *)

  (** Specifies the type of interpolation to use between checkpoints.

      @cvodes <Mathematics_link.html#checkpointing-scheme> Checkpointing scheme *)
  type interpolation = IPolynomial (** {cconst CV_POLYNOMIAL} *)
                     | IHermite    (** {cconst CV_HERMITE} *)

  (** Activates the forward-backward problem. The arguments specify the number
      of integration steps between consecutive checkpoints, and the type of
      variable-degree interpolation.

      @cvodes_adj CVodeAdjInit *)
  val init : ('d, 'k) Cvode.session -> int -> interpolation -> unit

  (** Integrates the forward problem over an interval and saves
      checkpointing data. The arguments are the next time at which a solution
      is desired ([tout]) and a vector to receive the computed result
      ([y]). The function returns a triple [tret, ncheck, sr]: the time
      reached by the solver, the cumulative number of checkpoints stored, and
      whether [tout] was reached. The solver takes internal steps until it
      has reached or just passed the [tout] parameter {cconst CV_NORMAL},
      it then interpolates to approximate [y(tout)].

      @cvodes_adj CVodeF
      @raise Cvode.IllInput           One of the inputs is invalid.
      @raise Cvode.TooMuchWork        Could not reach [tout] in [mxstep] steps
      @raise Cvode.TooMuchAccuracy    Could not satisfy the demanded accuracy
      @raise Cvode.ErrFailure         Too many error test failures.
      @raise Cvode.ConvergenceFailure Too many convergence test failures.
      @raise Cvode.LinearSetupFailure Unrecoverable failure in linear solver setup function.
      @raise Cvode.LinearSolveFailure Unrecoverable failure in linear solver solve function.
      @raise AdjointNotInitialized    The {!init} function has not been called. *)
  val forward_normal :
    ('d, 'k) Cvode.session
    -> float
    -> ('d, 'k) Nvector.t
    -> float * int * Cvode.solver_result

  (** Integrates the forward problem over an interval and saves
      checkpointing data. The arguments are the next time at which a solution
      is desired ([tout]) and a vector to receive the computed result
      ([yret]). The function returns a triple [tret, ncheck, sr]: the time
      reached by the solver, the cumulative number of checkpoints stored, and
      whether [tout] was reached. The solver takes one step
      {cconst CV_ONE_STEP} and returns the solution reached.

      @cvodes_adj CVodeF
      @raise Cvode.IllInput           One of the inputs is invalid.
      @raise Cvode.TooMuchWork        Could not reach [tout] in [mxstep] steps
      @raise Cvode.TooMuchAccuracy    Could not satisfy the demanded accuracy
      @raise Cvode.ErrFailure         Too many error test failures.
      @raise Cvode.ConvergenceFailure Too many convergence test failures.
      @raise Cvode.LinearSetupFailure Unrecoverable failure in linear solver setup function.
      @raise Cvode.LinearSolveFailure Unrecoverable failure in linear solver solve function.
      @raise AdjointNotInitialized    The {!init} function has not been called. *)
  val forward_one_step :
    ('d, 'k) Cvode.session
    -> float
    -> ('d, 'k) Nvector.t
    -> float * int * Cvode.solver_result

  (** {2:linear Linear solvers} *)

  (** Linear solvers used in backward problems.

      @cvodes_adj <Usage/SIM.html#cvodes-usage-sim-user-callable-lin-solv-init> Linear solver interface functions *)
  type ('data, 'kind) linear_solver =
          ('data, 'kind) Cvode_impl.AdjointTypes.linear_solver

  (** Alias for linear solvers that are restricted to serial nvectors. *)
  type 'kind serial_linear_solver =
    (Nvector_serial.data, 'kind) linear_solver
    constraint 'kind = [>Nvector_serial.kind]

  (** Workspaces with three temporary vectors. *)
  type 'd triple = 'd * 'd * 'd

  (** Arguments common to Jacobian callback functions.

      @cvodes_adj CVodeLsJacFnB
      @cvodes_adj CVodeJacTimesVecFnB
      @cvodes_adj CVodeLsPrecSolveFnB
      @cvodes_adj CVodeLsPrecSetupFnB *)
  type ('t, 'd) jacobian_arg = ('t, 'd) Cvode_impl.AdjointTypes.jacobian_arg =
    {
      jac_t   : float;        (** The independent variable. *)
      jac_y   : 'd;           (** The forward solution vector. *)
      jac_yb  : 'd;           (** The backward solution vector. *)
      jac_fyb : 'd;           (** The backward right-hand side function [fB]. *)
      jac_tmp : 't;           (** Temporary storage vectors. *)
    }

  (** Diagonal approximation of Jacobians by difference quotients. *)
  module Diag : sig (* {{{ *)
    (** A linear solver based on Jacobian approximation by difference
        quotients.

        @cvodes_adj CVDiagB *)
    val solver : ('data, 'kind) linear_solver

    (** Returns the sizes of the real and integer workspaces used by the
        Diagonal linear solver.

        @cvodes_adj CVDiagGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : ('d, 'k) bsession -> int * int

    (** Returns the number of calls made to the right-hand side
        function due to finite difference Jacobian approximation in the
        Diagonal linear solver.

        @cvodes_adj CVDiagGetNumRhsEvals *)
    val get_num_rhs_evals : ('d, 'k) bsession -> int
  end (* }}} *)

  (** Direct Linear Solvers operating on dense, banded, and sparse matrices. *)
  module Dls : sig (* {{{ *)
    include module type of Sundials_LinearSolver.Direct

    (** Callback functions that compute dense approximations to a Jacobian
        matrix without forward sensitivities. In the call [jac arg jm],
        [arg] is a {!jacobian_arg} with three work vectors and the
        computed Jacobian must be stored in [jm].

        The callback should load the [(i,j)]th entry of [jm] with
        {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
        the [i]th equation with respect to the [j]th variable, evaluated
        at the values of [t] and [y] obtained from [arg]. Only nonzero
        elements need be loaded into [jm].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jm] should
                 be accessed after the function has returned.}

        @cvodes_adj CVodeLsJacFnB *)
    type 'm jac_fn_no_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg -> 'm -> unit

    (** Callback functions that compute dense approximations to a Jacobian
        matrix with forward sensitivities. In the call [jac arg s jm],
        [arg] is a {!jacobian_arg} with three work vectors, [s] is an
        array of forward sensitivity vectors, and the computed Jacobian
        must be stored in [jm].

        The callback should load the [(i,j)]th entry of [jm] with
        {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
        the [i]th equation with respect to the [j]th variable, evaluated
        at the values of [t] and [y] obtained from [arg]. Only nonzero
        elements need be loaded into [jm].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg], [s] nor the matrix [jm]
                 should be accessed after the function has returned.}

        @cvodes_adj CVodeLsJacFnBS *)
    type 'm jac_fn_with_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg
      -> RealArray.t array
      -> 'm
      -> unit

    (** Callback functions that compute dense approximations to a Jacobian
        matrix.

        @cvodes_adj CVodeLsJacFnB
        @cvodes_adj CVodeLsJacFnBS *)
    type 'm jac_fn =
        NoSens of 'm jac_fn_no_sens
        (** Does not depend on forward sensitivities. *)
      | WithSens of 'm jac_fn_with_sens
        (** Depends on forward sensitivities. *)

    (** Function to compute the linear system matrix or an approximation to it
        without forward sensitivities. See {!linsys_fn} for details.

        @cvodes_adj CVLsLinSysFnB
        @since 5.0.0 *)
    type 'm linsys_fn_no_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg
      -> 'm
      -> bool
      -> float
      -> bool

    (** Function to compute the linear system matrix or an approximation to it
        with forward sensitivities. See {!linsys_fn} for details.

        @cvodes_adj CVLsLinSysFnBS
        @since 5.0.0 *)
    type 'm linsys_fn_with_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg
      -> RealArray.t array
      -> 'm
      -> bool
      -> float
      -> bool

    (** Function to compute the linear system matrix
        {% $M_B = I - \gamma_B J_B$ %}
        or an approximation of it for the backward problem. Offers an
        alternative to evaluating the Jacobian of the right-hand-side function.

        In addition to those shared with the Jacobian function, the arguments of
        this function are
        - [m], storage for the computed linear system matrix,
        - [jok], indicates whether the Jacobian-related data needs to be
                 updated, and
        - [gammab], the scalar in the formula above.

        The function should return true only if the Jacobian data was
        recomputed.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the Jacobian argument elements nor the matrix [m]
                 should be accessed after the function has returned.}

        @cvodes_adj CVLsLinSysFnB
        @cvodes_adj CVLsLinSysFnBS
        @since 5.0.0 *)
    type 'm linsys_fn =
        LNoSens of 'm linsys_fn_no_sens
        (** Does not depend on forward sensitivities. *)
      | LWithSens of 'm linsys_fn_with_sens
        (** Depends on forward sensitivities. *)

    (** Create a Cvodes-specific linear solver from a a Jacobian approximation
        function and a generic direct linear solver.

        The Jacobian approximation function is optional for dense and banded
        solvers (if not given an internal difference quotient approximation is
        used), but must be provided for other solvers (or [Invalid_argument]
        is raised).

        The [linsys] argument allows to override the standard linear system
        function that calls [jac] to compute {% $M_B$ %}. This feature is only
        available in Sundials >= 5.0.0.

        @cvodes_adj CVodeSetLinearSolverB
        @cvodes_adj CVodeSetJacFnB
        @cvodes_adj CVodeSetJacFnBS
        @cvodes_adj CVodeSetLinSysFnB
        @cvodes_adj CVodeSetLinSysFnBS *)
    val solver :
      ?jac:'m jac_fn ->
      ?linsys:'m linsys_fn ->
      ('m, RealArray.t, 'kind, 't) LinearSolver.t ->
      'kind serial_linear_solver

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by a direct
        linear solver.

        @cvodes CVodeGetLinWorkSpace
        @cvodes_adj CVodeGetAdjCVodeBmem
        @return ([real_size], [integer_size]) *)
    val get_work_space : 'k serial_bsession -> int * int

    (** Returns the number of calls made by a direct linear solver to the
        Jacobian approximation function.

        @cvodes CVodeGetNumJacEvals
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_jac_evals : 'k serial_bsession -> int

    (** Returns the number of calls to the right-hand side callback due to
        the finite difference Jacobian approximation.

        @cvodes CVodeGetNumLinRhsEvals
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_lin_rhs_evals : 'k serial_bsession -> int

  end (* }}} *)

  (** Scaled Preconditioned Iterative Linear Solvers.

      @cvodes_adj CVodePrecSolveFnB
      @cvodes_adj CVodePrecSetupFnB *)
  module Spils : sig (* {{{ *)
    include module type of Sundials_LinearSolver.Iterative

    (** {3:precond Preconditioners} *)

    (** Arguments passed to the preconditioner solver function.

        @cvodes_adj CVodePrecSolveFnB *)
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
        preconditioner matrix without forward sensitivities.
        In the call [prec_solve_fn jac arg z],
        - [jac] is a {!jacobian_arg} with no work vectors,
        - [arg] is {!prec_solve_arg} that specifies the linear system, and
        - [z] is computed to solve {% $P\mathtt{z} = \mathtt{arg.rhs}$%}.
        $P$ is a preconditioner matrix, which approximates, however crudely,
        the Newton matrix {% $M = I - \gamma J$%} where
        {% $J = \frac{\partial f}{\partial y}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac], [arg], and [z] should not be
        accessed after the function has returned.}

        @cvodes_adj CVodePrecSolveFnB *)
    type 'd prec_solve_fn =
      (unit, 'd) jacobian_arg
      -> 'd prec_solve_arg
      -> 'd
      -> unit

    (** Callback functions that solve a linear system involving a
        preconditioner matrix with forward sensitivities.
        In the call [prec_solve_fn jac arg s z],
        - [jac] is a {!jacobian_arg} with no work vectors,
        - [arg] is {!prec_solve_arg} that specifies the linear system,
        - [s] is an array of forward sensitivity vectors, and
        - [z] is computed to solve {% $P\mathtt{z} = \mathtt{arg.rhs}$%}.
        $P$ is a preconditioner matrix, which approximates, however crudely,
        the Newton matrix {% $M = I - \gamma J$%} where
        {% $J = \frac{\partial f}{\partial y}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac], [arg], [s], and [z] should not be
                 accessed after the function has returned.}

        @cvodes_adj CVodePrecSolveFnBS *)
    type 'd prec_solve_fn_with_sens =
      (unit, 'd) jacobian_arg
      -> 'd prec_solve_arg
      -> 'd array
      -> 'd
      -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related data
        needed by {!prec_solve_fn} without forward sensitivities.
        In the call [prec_setup_fn jac s jok gamma],
        - [jac] is a {!jacobian_arg} with no work vectors,
        - [jok] indicates whether any saved Jacobian-related data can be
                reused with the current value of [gamma], and
        - [gamma] is the scalar $\gamma$ in the Newton matrix
                  {% $M = I - \gamma J$%} where $J$ is the Jacobian
                  matrix.
        A function should return [true] if Jacobian-related data was
        updated and [false] if saved data was reused.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac] should not be accessed after the
                 function has returned.}

        @cvodes_adj CVodePrecSetupFnB *)
    type 'd prec_setup_fn =
      (unit, 'd) jacobian_arg
      -> bool
      -> float
      -> bool

    (** Callback functions that preprocess or evaluate Jacobian-related data
        needed by {!prec_solve_fn} with forward sensitivities.
        In the call [prec_setup_fn jac s jok gamma],
        - [jac] is a {!jacobian_arg} with no work vectors,
        - [s] is an array of forward sensitivity vectors,
        - [jok] indicates whether any saved Jacobian-related data can be
                reused with the current value of [gamma], and
        - [gamma] is the scalar $\gamma$ in the Newton matrix
                  {% $M = I - \gamma J$%} where $J$ is the Jacobian
                  matrix.
        A function should return [true] if Jacobian-related data was
        updated and [false] if saved data was reused.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac] should not be accessed after the
                 function has returned.}

        @cvodes_adj CVodePrecSetupFnBS *)
    type 'd prec_setup_fn_with_sens =
      (unit, 'd) jacobian_arg
      -> 'd array
      -> bool
      -> float
      -> bool

    (** Specifies a preconditioner, including the type of preconditioning
        (none, left, right, or both) and callback functions. The following
        functions and those in {!Banded} and {!Cvodes_bbd} construct
        preconditioners.

        The {!prec_solve_fn} is mandatory. The {!prec_setup_fn} can
        be omitted if not needed.

        @cvodes_adj CVodeSetPreconditionerB
        @cvodes_adj CVodePrecSolveFnB
        @cvodes_adj CVodePrecSetupFnB
        @cvodes_adj CVodePrecSolveFnBS
        @cvodes_adj CVodePrecSetupFnBS *)
    type ('d, 'k) preconditioner =
      ('d, 'k) Cvode_impl.AdjointTypes.SpilsTypes.preconditioner

    (** No preconditioning.  *)
    val prec_none : ('d, 'k) preconditioner

    (** Left preconditioning without forward sensitivities.
        {% $(P^{-1}A)x = P^{-1}b$ %}. *)
    val prec_left :
      ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

    (** Left preconditioning with forward sensitiviites.
        {% $(P^{-1}A)x = P^{-1}b$ %}. *)
    val prec_left_with_sens :
      ?setup:'d prec_setup_fn_with_sens
      -> 'd prec_solve_fn_with_sens
      -> ('d, 'k) preconditioner

    (** Right preconditioning with sensitivities.
        {% $(AP^{-1})Px = b$ %}. *)
    val prec_right :
      ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

    (** Right preconditioning without sensitivities.
        {% $(AP^{-1})Px = b$ %}. *)
    val prec_right_with_sens :
      ?setup:'d prec_setup_fn_with_sens
      -> 'd prec_solve_fn_with_sens
      -> ('d, 'k) preconditioner

    (** Left and right preconditioning with sensitivities.
        {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)
    val prec_both :
      ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

    (** Left and right preconditioning without sensitivities.
        {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)
    val prec_both_with_sens :
      ?setup:'d prec_setup_fn_with_sens
      -> 'd prec_solve_fn_with_sens
      -> ('d, 'k) preconditioner

    (** Banded preconditioners.  *)
    module Banded : sig (* {{{ *)

      (** The range of nonzero entries in a band matrix. *)
      type bandrange =
        { mupper : int; (** The upper half-bandwidth.  *)
          mlower : int; (** The lower half-bandwidth.  *) }

      (** A band matrix {!preconditioner} based on difference quotients.
          The call [prec_left br] instantiates a left preconditioner which
          generates a banded approximation to the Jacobian with [br.mlower]
          sub-diagonals and [br.mupper] super-diagonals.

          @cvodes_adj CVBandPrecInitB *)
      val prec_left : bandrange -> (Nvector_serial.data,
                                    [>Nvector_serial.kind]) preconditioner

      (** Like {!prec_left} but preconditions from the right.

          @cvodes_adj CVBandPrecInitB *)
      val prec_right : bandrange -> (Nvector_serial.data,
                                     [>Nvector_serial.kind]) preconditioner

      (** Like {!prec_left} but preconditions from both sides.

          @cvodes_adj CVBandPrecInitB *)
      val prec_both : bandrange -> (Nvector_serial.data,
                                    [>Nvector_serial.kind]) preconditioner

      (** {4:stats Banded statistics} *)

      (** Returns the sizes of the real and integer workspaces
          used by the banded preconditioner module.

          @cvodes_adj CVBandPrecGetWorkSpace
          @cvodes_adj CVodeGetAdjCVodeBmem
          @return ([real_size], [integer_size]) *)
      val get_work_space : 'k serial_bsession -> int * int

      (** Returns the number of calls to the right-hand side callback for the
          difference banded Jacobian approximation. This counter is only updated
          if the default difference quotient function is used.

          @cvodes_adj CVBandPrecGetNumRhsEvals
          @cvodes_adj CVodeGetAdjCVodeBmem *)
      val get_num_rhs_evals : 'k serial_bsession -> int
    end (* }}} *)

    (** {3:lsolvers Solvers} *)

    (** Callback functions that preprocess or evaluate Jacobian-related
        data needed by the jac_times_vec_fn. In the call
        [jac_times_setup_fn arg], [arg] is a {!jacobian_arg} with no
        work vectors.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [arg] should not be accessed after the
                 function has returned.}

        @cvodes_adj CVodeSetJacTimesB
        @cvodes_adj CVodeJacTimesSetupFnB *)
    type 'd jac_times_setup_fn_no_sens = (unit, 'd) jacobian_arg -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related
        data needed by the jac_times_vec_fn. In the call
        [jac_times_setup_fn arg s], [arg] is a {!jacobian_arg} with no
        work vectors and [s] is an array of forward sensitivity vectors.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [arg] should not be accessed after the
                 function has returned.}

        @cvodes_adj CVodeSetJacTimesBS
        @cvodes_adj CVodeJacTimesSetupFnBS *)
    type 'd jac_times_setup_fn_with_sens =
      (unit, 'd) jacobian_arg -> 'd array -> unit

    (** Callback functions that compute the Jacobian times a vector without
        forward sensitivities.
        In the call [jac_times_vec_fn arg v jv],
        - [arg] is a {!jacobian_arg} with one work vector,
        - [v] is the vector multiplying the Jacobian, and
        - [jv] is the vector in which to store the
               result—{% $\mathtt{jv} = J\mathtt{v}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor [v] or [jv] should be
                 accessed after the function has returned.}

        @cvodes_adj CVodeJacTimesVecFnB *)
    type 'd jac_times_vec_fn_no_sens =
      ('d, 'd) jacobian_arg
      -> 'd
      -> 'd
      -> unit

    (** Callback functions that compute the Jacobian times a vector with
        forward sensitivities.
        In the call [jac_times_vec_fn arg s v jv],
        - [arg] is a {!jacobian_arg} with one work vector,
        - [s] is an array of forward sensitivity vectors,
        - [v] is the vector multiplying the Jacobian, and
        - [jv] is the vector in which to store the
               result—{% $\mathtt{jv} = J\mathtt{v}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg], [s], [v], nor [jv] should be
                 accessed after the function has returned.}

        @cvodes_adj CVodeJacTimesVecFnBS *)
    type 'd jac_times_vec_fn_with_sens =
      ('d, 'd) jacobian_arg
      -> 'd array
      -> 'd
      -> 'd
      -> unit

    (** Callback functions that compute the Jacobian times a vector.

        @cvodes_adj CVodeJacTimesSetupFnB
        @cvodes_adj CVodeJacTimesSetupFnBS
        @cvodes_adj CVodeJacTimesVecFnB
        @cvodes_adj CVodeJacTimesVecFnBS *)
    type 'd jac_times_vec_fn =
      | NoSens of 'd jac_times_setup_fn_no_sens option
                  * 'd jac_times_vec_fn_no_sens
        (** Does not depend on forward sensitivities. *)
      | WithSens of 'd jac_times_setup_fn_with_sens option
                    * 'd jac_times_vec_fn_with_sens
        (** Depends on forward sensitivities. *)

    (** Create a Cvodes-specific linear solver from a generic iterative
        linear solver.

        The [jac_times_rhs] argument specifies an alternative right-hand-side
        function for use in the internal Jacobian-vector product difference
        quotient approximation. It is incorrect to specify both this argument
        and [jac_times_vec].

        NB: the [jac_times_setup] argument is not supported in
            {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

        NB: a [jac_times_rhs] function is not supported in
            {{!Sundials_Config.sundials_version}Config.sundials_version} < 5.3.0.

        @cvodes_adj CVodeSetLinearSolverB
        @cvodes_adj CVodeSetJacTimesB
        @cvodes_adj CVodeSetJacTimesBS
        @cvodes_adj CVodeSetJacTimesRhsFnB *)
    val solver :
      ('m, 'd, 'k, 'f) LinearSolver.t
      -> ?jac_times_vec:'d jac_times_vec_fn
      -> ?jac_times_rhs:'d Cvode.rhsfn
      -> ('d, 'k) preconditioner
      -> ('d, 'k) linear_solver

    (** {3:set Solver parameters} *)

    (** Sets the maximum number of time steps to wait before recomputation of
        the Jacobian or recommendation to update the preconditioner.
        If the integer argument is less than or equal to 0, a default value
        of 50 is used.

        @cvodes_adj CVodeSetJacEvalFrequency
        @since 4.0.0 *)
    val set_jac_eval_frequency : ('d, 'k) bsession -> int -> unit

    (** Specifies the frequency of calls to the linear solver setup routine.
        Positive values specify the number of time steps between setup calls,
        negative values force recomputation at each Newton step, and zero
        values reset to the default (20).

        @cvodes_adj CVodeSetLSetupFrequency
        @since 5.4.0 *)
    val set_lsetup_frequency : ('d, 'k) bsession -> int -> unit

    (** Sets the factor by which the Krylov linear solver's convergence test
        constant is reduced from the Newton iteration test constant.
        This factor must be >= 0; passing 0 specifies the default (0.05).

        @cvodes_adj CVodeSetEpsLinB *)
    val set_eps_lin : ('d, 'k) bsession -> float -> unit

    (** Sets the factor for converting from the integrator tolerance (WRMS
        norm) to the linear solver tolerance (L2 norm). That is,
        {% $\mathit{tol}_{\mathsf{L2}} =
            \mathit{fact}\cdot\mathit{tol}_{\mathsf{WRMS}}$ %}.
        The given value is used directly if it is greater than zero.
        If it is zero (the default), then the square root of the state
        vector length is used.
        If it is less than zero, then the square root of the dot product of a
        state vector full of ones with itself is used.

        @cvodes_adj CVodeSetLSNormFactorB
        @since 5.4.0 *)
    val set_ls_norm_factor : ('d, 'k) bsession -> float -> unit

    (** Enables or disables scaling of the linear system solution to account
        for a change in {% $\gamma$ %} in the linear system. Linear solution
        scaling is enabled by default when a matrix-based linear solver is
        attached.

        @cvodes_adj CVodeSetLinearSolutionScalingB
        @since 5.2.0 *)
    val set_linear_solution_scaling : ('d, 'k) bsession -> bool -> unit

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the spils
        linear solver.

        @cvodes_adj CVodeGetWorkSpace
        @cvodes_adj CVodeGetAdjCVodeBmem
        @return ([real_size], [integer_size]) *)
    val get_work_space       : ('d, 'k) bsession -> int * int

    (** Returns the cumulative number of linear iterations.

        @cvodes_adj CVodeGetNumLinIters
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_lin_iters    : ('d, 'k) bsession -> int

    (** Returns the cumulative number of linear convergence failures.

        @cvodes_adj CVodeGetNumLinConvFails
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_lin_conv_fails   : ('d, 'k) bsession -> int

    (** Returns the cumulative number of calls to the setup function with
        [jok=false].

        @cvodes_adj CVodeGetNumPrecEvals
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_prec_evals   : ('d, 'k) bsession -> int

    (** Returns the cumulative number of calls to the preconditioner solve
        function.

        @cvodes_adj CVodeGetNumPrecSolves
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_prec_solves  : ('d, 'k) bsession -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        setup function.

        @cvodes_adj CVodeGetNumJTSetupEvals
        @cvodes_adj CVodeGetAdjCVodeBmem
        @since 3.0.0 *)
    val get_num_jtsetup_evals : ('d, 'k) bsession -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        function.

        @cvodes_adj CVodeGetNumJtimesEvals
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_jtimes_evals : ('d, 'k) bsession -> int

    (** Returns the number of calls to the right-hand side callback for
        finite difference Jacobian-vector product approximation. This counter is
        only updated if the default difference quotient function is used.

        @cvodes_adj CVodeGetNumLinRhsEvals
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_lin_rhs_evals    : ('d, 'k) bsession -> int

    (** {3:lowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change preconditioner functions. These low-level functions are
        provided for experts who want to avoid resetting internal counters
        and other associated side-effects. *)

    (** Change the preconditioner functions without using forward
        sensitivities.

        @cvodes_adj CVodeSetPreconditionerB
        @cvodes_adj CVodePrecSolveFnB
        @cvodes_adj CVodePrecSetupFnB *)
    val set_preconditioner :
      ('d,'k) bsession
      -> ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> unit

    (** Change the preconditioner functions using forward sensitivities.

        @cvodes_adj CVodeSetPreconditionerBS
        @cvodes_adj CVodePrecSolveFnBS
        @cvodes_adj CVodePrecSetupFnBS *)
    val set_preconditioner_with_sens :
      ('d,'k) bsession
      -> ?setup:'d prec_setup_fn_with_sens
      -> 'd prec_solve_fn_with_sens
      -> unit

    (** Change the Jacobian-times-vector function.

        @cvodes_adj CVodeSetJacTimesVecFnB
        @cvodes_adj CVodeSetJacTimesVecFnBS
        @cvodes_adj CVodeJacTimesVecFnB
        @cvodes_adj CVodeJacTimesVecFnBS *)
    val set_jac_times :
      ('d,'k) bsession
      -> 'd jac_times_vec_fn
      -> unit

    (** Remove a Jacobian-times-vector function and use the default
        implementation.

        @cvodes_adj CVodeSetJacTimesVecFnB *)
    val clear_jac_times : ('d, 'k) bsession -> unit

  end (* }}} *)

  (** Create a CVode-specific linear solver from a generic matrix embedded
      solver.

      @cvodes_adj CVodeSetLinearSolver
      @since 5.8.0 *)
  val matrix_embedded_solver :
         (unit, 'data, 'kind, [>`MatE]) LinearSolver.t
      -> ('data, 'kind) linear_solver

  (** {2:bsolve Backward solutions} *)

  (** Arguments common to {!brhsfn_no_sens} and {!brhsfn_with_sens}.  *)
  type 'd brhsfn_args = 'd Cvode_impl.AdjointTypes.brhsfn_args =
    {
      t : float;
      (** The value of the independent variable. *)

      y : 'd;
      (** The vector of dependent-variable values $y(t)$. *)

      yb : 'd;
      (** The vector of backward dependent-variable values $y_B(t)$. *)
    }

  (** Backward functions without forward sensitivities. They are passed
      the arguments:
      - [args], the current values of forward and backward variables, and,
      - [yb'], a vector for storing the values
               {% $\dot{y}_B = f_B(t, y, y_B)$%}.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning The vectors in the function arguments should not
               be accessed after the function returns.}

      @cvodes_adj CVRhsFnB *)
  type 'd brhsfn_no_sens = 'd brhsfn_args -> 'd -> unit

  (** Backward functions with forward sensitivities. They are passed the
      arguments:
      - [args], the current values of state and backward sensitivity variables,
      - [s], an array holding the values of forward sensitivity vectors, and,
      - [yb'], a vector for storing the values
               {% $\dot{y}_B = f_B(t, y, y_S, y_B)$%}.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning The vectors in the function arguments should not
               be accessed after the function returns.}

      @cvodes_adj CVRhsFnBS *)
  type 'd brhsfn_with_sens = 'd brhsfn_args -> 'd array -> 'd -> unit

  (** Functions that evaluate the right-hand side of a backward ODE system
      with or without forward sensitivities. *)
  type 'd brhsfn =
      NoSens of 'd brhsfn_no_sens
        (** No dependency on forward sensitivities. *)
    | WithSens of 'd brhsfn_with_sens
        (** Dependency on forward sensitivities. *)

  (** Tolerance specifications. *)
  type ('d, 'k) tolerance =
    | SStolerances of float * float
      (** [(rel, abs)] : scalar relative and absolute tolerances. *)
    | SVtolerances of float * ('d, 'k) Nvector.t
      (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

  (** Creates and initializes a backward session attached to an existing
      (forward) session. The call
      {[init_backward s lmm iter tol fb tb0 yb0]} has as arguments:
      - [s], the parent (forward) session,
      - [lmm], the linear multistep method (see {!Cvode.lmm}),
      - [tol], the integration tolerances,
      - [nlsolver], the solver to use to calculate integration steps,
      - [lsolver],  used by [nlsolver]s based on Newton interation,
      - [fb], the backward right-hand side function,
      - [tb0], specifies the endpoint where final conditions are provided
               for the backward problem, which is normally the endpoint of
               forward integration, and,
      - [yb0], a vector of final values that also determines the number
               of equations.

      This function does everything necessary to initialize a backward
      session, i.e., it makes the calls referenced below. The
      {!backward_normal} and {!backward_one_step} functions may be called
      directly.

      If an [nlsolver] is not specified, then the
      {{!Sundials_NonlinearSolver.Newton}Newton} module is used by default.
      In this case only, [lsolver] defaults to {!Diag.solver} if not otherwise
      specified. Specifying an [nlsolver] that requires a linear solver without
      specifying an [lsolver] results in a {!Cvode.NonlinearInitFailure} (or
      {!Cvode.IllInput} for Sundials < 4.0.0) exception on the first call to
      {!backward_normal} or {!backward_one_step}.

      @cvodes_adj CVodeCreateB
      @cvodes_adj CVodeInitB
      @cvodes_adj CVodeInitBS
      @cvodes_adj CVodeSetLinearSolverB
      @cvodes_adj CVodeSetNonlinearSolverB
      @cvodes_adj CVodeSStolerancesB
      @cvodes_adj CVodeSVtolerancesB
      @cvodes_adj CVodeSVtolerancesB
      @raise AdjointNotInitialized The {!init} function has not been called.
      @raise BadFinalTime The final time is outside the interval over which the forward problem was solved. *)
  val init_backward :
       ('d, 'k) Cvode.session
    -> Cvode.lmm
    -> ('d, 'k) tolerance
    -> ?nlsolver : ('d, 'k, ('d, 'k) Cvode.session, [`Nvec]) Sundials.NonlinearSolver.t
    -> ?lsolver  : ('d, 'k) linear_solver
    -> 'd brhsfn
    -> float
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) bsession

  (** Support for backward quadrature equations that may or may
      not depend on forward sensitivities.

      @cvodes_adj <Usage/ADJ.html#backward-integration-of-quadrature-equations> Backward integration of quadrature equations *)
  module Quadrature : sig (* {{{ *)
    (** {2:init Initialization} *)

    (** Arguments common to {!bquadrhsfn_no_sens} and
       {!bquadrhsfn_with_sens}. *)
    type 'd bquadrhsfn_args =
      {
        t : float;
        (** The value of the independent variable. *)

        y : 'd;
        (** The vector of dependent-variable values $y(t)$. *)

        yb : 'd;
        (** The vector of backward dependent-variable values $y_B(t)$. *)
      }

    (** Functions defining backward quadrature variables without forward
        sensitivities.  These functions are passed the arguments:
        - [args], the current values of forward and backward variables, and,
        - [qb'], a vector for storing the computed value of
                 {% $\dot{y}_\mathit{QB} = f_\mathit{QB}(t, y, y_B)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure}
        exception indicates a recoverable error. Any other exception is
        treated as an unrecoverable error.

        {warning The vectors in the function arguments should not
                 be accessed after the function returns.}

        @cvodes_adj CVQuadRhsFnB *)
    type 'd bquadrhsfn_no_sens = 'd bquadrhsfn_args -> 'd -> unit

    (** Functions defining backward quadrature variables with forward
        sensitivities.  These functions are passed the arguments:
        - [args], the current values of forward and backward variables,
        - [s], an array holding the values of forward sensitivity vectors,
               and,
        - [qb'], a vector for storing the computed value of
               {% $\dot{y}_\mathit{QB} = f_\mathit{QB}(t, y, y_S, y_B)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure}
        exception indicates a recoverable error. Any other exception is
        treated as an unrecoverable error.

        {warning The vectors in the function arguments should not
                 be accessed after the function returns.}

        @cvodes_adj CVQuadRhsFnBS *)
    type 'd bquadrhsfn_with_sens =
      'd bquadrhsfn_args -> 'd array -> 'd -> unit

    (** These functions compute the quadrature equation right-hand side for
        the backward problem. *)
    type 'd bquadrhsfn =
        NoSens of 'd bquadrhsfn_no_sens
        (** Does not depend on forward sensitivities. *)
      | WithSens of 'd bquadrhsfn_with_sens
        (** Depends on forward sensitivities. *)

    (** This function activates the integration of quadrature equations.
        The arguments are the function that computes the right-hand side of
        the backward quadrature equations, and a vector giving the values
        of the quadrature variables at [tB0].

        @cvodes_adj CVodeQuadInitB
        @cvodes_adj CVodeQuadInitBS *)
    val init : ('d, 'k) bsession -> 'd bquadrhsfn
             -> ('d, 'k) Nvector.t -> unit

    (** This function reinitializes the integration of quadrature equations
        during the backward phase.

        @cvodes_adj CVodeQuadReInitB *)
    val reinit : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

    (** {2:tols Tolerance specification} *)

    (** Tolerances for calculating backward quadrature variables. *)
    type ('d, 'k) tolerance =
        NoStepSizeControl
        (** Quadrature variables are not used for step-size control
            (the default). *)
      | SStolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('d, 'k) Nvector.t
        (** [(rel, abs)] : scalar relative and vector absolute
            tolerances. *)

    (** Specify how to use quadrature variables in step size control.

        @cvodes_adj CVodeSetQuadErrCon
        @cvodes_adj CVodeQuadSStolerances
        @cvodes_adj CVodeQuadSVtolerances *)
    val set_tolerances : ('d, 'k) bsession -> ('d, 'k) tolerance -> unit

    (** {2:quadout Output functions} *)

    (** Returns the backward quadrature solutions and time reached
        after a successful solver step. The given vectors are filled with
        values calculated during either {!backward_normal} or
        {!backward_one_step} and the value of the independent variable
        is returned.

      @cvodes_adj CVodeGetQuadB *)
    val get : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> float

    (** {2:get Querying the solver (optional output functions)}

        @cvodes_adj <Usage/ADJ.html#optional-input-output-functions-for-backward-quadrature-integration> Optional input/output functions for backward quadrature integration *)

    (** Returns the number of calls to the backward quadrature right-hand
        side function.

        @cvodes_adj CVodeGetQuadNumRhsEvals
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_rhs_evals       : ('d, 'k) bsession -> int

    (** Returns the number of local error test failures due to quadrature
        variables.

        @cvodes_adj CVodeGetQuadNumErrTestFails
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_num_err_test_fails  : ('d, 'k) bsession -> int

    (** Returns the quadrature error weights at the current time.

        @cvodes_adj CVodeGetQuadErrWeights
        @cvodes_adj CVodeGetAdjCVodeBmem *)
    val get_err_weights : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

    (** Returns quadrature-related statistics. These are the
        number of calls to the quadrature function ([nfqevals]) and the
        number of error test failures due to quadrature variables
        ([nqetfails]).

        @cvodes_adj CVodeGetQuadStats
        @cvodes_adj CVodeGetAdjCVodeBmem
        @return ([nfqevals], [nqetfails]) *)
    val get_stats : ('d, 'k) bsession -> int * int
  end (* }}} *)

  (** Integrates a backward ODE system over an interval. The solver takes
      internal steps until it has reached or just passed the specified value.

      @cvodes_adj CVodeB (CV_NORMAL)
      @raise AdjointNotInitialized    The {!init} function has not been called.
      @raise NoBackwardProblem        The {!init_backward} function has not been called.
      @raise NoForwardCall            Neither {!forward_normal} nor {!forward_one_step} has been called.
      @raise Cvode.IllInput           One of the inputs is invalid.
      @raise Cvode.TooMuchWork        Could not reach [tout] in [mxstep] steps
      @raise Cvode.TooMuchAccuracy    Could not satisfy the demanded accuracy
      @raise Cvode.ErrFailure         Too many error test failures.
      @raise Cvode.ConvergenceFailure Too many convergence test failures.
      @raise Cvode.LinearSetupFailure Unrecoverable failure in linear solver setup function.
      @raise Cvode.LinearSolveFailure Unrecoverable failure in linear solver solve function.
      @raise BadOutputTime            The requested output time is outside the interval over which the forward problem was solved.
      @raise ForwardReinitializationFailed Reinitialization of the forward problem failed at the first checkpoint (corresponding to the initial time of the forward problem).
      @raise ForwardFail              An error occurred during the integration of the forward problem. *)
  val backward_normal : ('d, 'k) Cvode.session -> float -> unit

  (** Like {!backward_normal} but returns after one internal solver step.

      @cvodes_adj CVodeB (CV_ONE_STEP) *)
  val backward_one_step : ('d, 'k) Cvode.session -> float -> unit

  (** Fills the given vector with the solution of the backward ODE problem at
      the returned time, interpolating if necessary.

      @cvodes_adj CVodeGetB *)
  val get : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> float

  (** Returns the interpolated solution or derivatives.
      [get_dky s dkyb t k] computes the [k]th derivative of the backward
      function at time [t], i.e.,
      {% $\frac{d^\mathtt{k}y_B(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
      and stores it in [dkyb]. The arguments must satisfy
      {% $t_n - h_u \leq \mathtt{t} \leq t_n$%}—where $t_n$
      denotes {!get_current_time} and $h_u$ denotes {!get_last_step},—
      and {% $0 \leq \mathtt{k} \leq q_u$%}—where $q_u$ denotes
      {!get_last_order}.

      This function may only be called after a successful return from either
      {!backward_normal} or {!backward_one_step}.

      @cvodes_adj CVodeGetDky
      @cvodes_adj CVodeGetAdjIDABmem
      @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
      @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
  val get_dky
        : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> float -> int -> unit

  (** Fills the vector with the interpolated forward solution at the
      given time during a backward simulation.

      @cvodes_adj CVodeGetAdjY *)
  val get_y : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t  -> float -> unit

  (** Reinitializes the backward problem with new parameters and state
      values. The values of the independent variable, i.e., the simulation
      time, and the state variables must be given. It is also possible to
      change the solution method (and linear solver).

      @cvodes_adj CVodeReInitB
      @cvodes_adj CVodeSetLinearSolverB
      @cvodes_adj CVodeSetNonlinearSolverB
      @raise AdjointNotInitialized The {!init} function has not been called.
      @raise BadFinalTime The final time is not within the forward problem solution interval. *)
  val reinit :
    ('d, 'k) bsession
    -> ?nlsolver : ('d, 'k, ('d, 'k) Cvode.session, [`Nvec]) Sundials.NonlinearSolver.t
    -> ?lsolver  : ('d, 'k) linear_solver
    -> float
    -> ('d, 'k) Nvector.t
    -> unit

  (** {2:set Modifying the solver (optional input functions)} *)

  (** Cancels the storage of sensitivity checkpointing data during forward
      solution (with {!forward_normal} or {!forward_one_step}).

      @cvodes_adj CVodeAdjSetNoSensi *)
  val set_no_sensitivity : ('d, 'k) Cvode.session -> unit

  (** Sets the integration tolerances for the backward problem.

      @cvodes_adj CVodeSStolerancesB
      @cvodes_adj CVodeSVtolerancesB *)
  val set_tolerances : ('d, 'k) bsession -> ('d, 'k) tolerance -> unit

  (** Specifies the maximum order of the linear multistep method.

      @cvodes_adj CVodeSetMaxOrdB *)
  val set_max_ord : ('d, 'k) bsession -> int -> unit

  (** Specifies the maximum number of steps taken in attempting to reach
      a given output time.

      @cvodes_adj CVodeSetMaxNumStepsB *)
  val set_max_num_steps : ('d, 'k) bsession -> int -> unit

  (** Specifies the initial step size.

      @cvodes_adj CVodeSetInitStepB *)
  val set_init_step : ('d, 'k) bsession -> float -> unit

  (** Specifies a lower bound on the magnitude of the step size.

      @cvodes_adj CVodeSetMinStepB *)
  val set_min_step : ('d, 'k) bsession -> float -> unit

  (** Specifies an upper bound on the magnitude of the step size.

      @cvodes_adj CVodeSetMaxStepB *)
  val set_max_step : ('d, 'k) bsession -> float -> unit

  (** Indicates whether the BDF stability limit detection algorithm should be
      used.

      @cvodes CVodeSetStabLimDet
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val set_stab_lim_det : ('d, 'k) bsession -> bool -> unit

  (** Specifies a vector defining inequality constraints for each
      component of the solution vector [y].  See {!Sundials.Constraint}.

      @cvodes_adj CVodeSetConstraintsB *)
  val set_constraints : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

  (** Disables constraint checking.

      @cvodes_adj CVodeSetConstraints *)
  val clear_constraints : ('d, 'k) bsession -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the real and integer workspace sizes.

      @cvodes_adj CVodeGetWorkSpace
      @cvodes_adj CVodeGetAdjCVodeBmem
      @return ([real_size], [integer_size]) *)
  val get_work_space          : ('d, 'k) bsession -> int * int

  (** Returns the cumulative number of internal steps taken by the solver.

      @cvodes_adj CVodeGetNumSteps
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_num_steps           : ('d, 'k) bsession -> int

  (** Returns the number of calls to the backward right-hand side function.

      @cvodes_adj CVodeGetNumRhsEvals
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_num_rhs_evals       : ('d, 'k) bsession -> int

  (** Returns the number of calls made to the linear solver's setup function.

      @cvodes_adj CVodeGetNumLinSolvSetups
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_num_lin_solv_setups : ('d, 'k) bsession -> int

  (** Returns the number of local error test failures that have occurred.

      @cvodes_adj CVodeGetNumErrTestFails
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_num_err_test_fails  : ('d, 'k) bsession -> int

  (** Returns the integration method order used during the last internal step.

      @cvodes_adj CVodeGetLastOrder
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_last_order          : ('d, 'k) bsession -> int

  (** Returns the integration method order to be used on the next internal
      step.

      @cvodes_adj CVodeGetCurrentOrder
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_current_order       : ('d, 'k) bsession -> int

  (** Returns the integration step size taken on the last internal step.

      @cvodes_adj CVodeGetLastStep
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_last_step           : ('d, 'k) bsession -> float

  (** Returns the integration step size to be attempted on the next internal
      step.

      @cvodes_adj CVodeGetCurrentStep
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_current_step        : ('d, 'k) bsession -> float

  (** Returns the the value of the integration step size used on the first
      step.

      @cvodes_adj CVodeGetActualInitStep
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_actual_init_step    : ('d, 'k) bsession -> float

  (** Returns the the current internal time reached by the solver.

      @cvodes_adj CVodeGetCurrentTime
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_current_time        : ('d, 'k) bsession -> float

  (** Returns the number of order reductions dictated by the BDF stability
      limit detection algorithm.

      @cvodes_adj CVodeGetNumStabLimOrderReds
      @cvodes_adj CVodeGetAdjCVodeBmem
      @cvodes_adj <Mathematics_link.html#bdf-stability-limit-detection> BDF stability limit detection *)
  val get_num_stab_lim_order_reds : ('d, 'k) bsession -> int

  (** Returns a suggested factor by which the user's tolerances should be
      scaled when too much accuracy has been requested for some internal
      step.

      @cvodes_adj CVodeGetTolScaleFactor
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_tol_scale_factor : ('d, 'k) bsession -> float

  (** Returns the solution error weights at the current time.

      @cvodes_adj CVodeGetErrWeights
      @cvodes_adj CVodeGetAdjCVodeBmem
      @cvodes_adj <Mathematics_link.html#equation-cvodes-errwt> IVP solution (W_i) *)
  val get_err_weights : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

  (** Returns the vector of estimated local errors.

      @cvodes_adj CVodeGetEstLocalErrors
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_est_local_errors : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

  (** Returns the integrator statistics as a group.

      @cvodes_adj CVodeGetIntegratorStats
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_integrator_stats    : ('d, 'k) bsession -> Cvode.integrator_stats

  (** Prints the integrator statistics on the given channel.

      @cvodes_adj CVodeGetIntegratorStats
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val print_integrator_stats  : ('d, 'k) bsession -> out_channel -> unit

  (** Returns the cumulative number of nonlinear (functional or Newton)
      iterations.

      @cvodes_adj CVodeGetNumNonlinSolvIters
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_num_nonlin_solv_iters : ('d, 'k) bsession -> int

  (** Returns the cumulative number of nonlinear convergence failures.

      @cvodes_adj CVodeGetNumNonlinSolvConvFails
      @cvodes_adj CVodeGetAdjCVodeBmem *)
  val get_num_nonlin_solv_conv_fails : ('d, 'k) bsession -> int

  (** Returns both the numbers of nonlinear iterations performed [nniters] and
      nonlinear convergence failures [nncfails].

      @cvodes CVodeGetNonlinSolvStats
      @cvodes_adj CVodeGetAdjCVodeBmem
      @return ([nniters], [nncfails]) *)
  val get_nonlin_solv_stats : ('d, 'k) bsession -> int *int

  (** {2:exceptions Exceptions} *)

  (** Adjoint sensitivity analysis was not initialized.

      @cvodes <Constants_link.html> CV_NO_ADJ *)
  exception AdjointNotInitialized

  (** Neither {!forward_normal} nor {!forward_one_step} has been called.

      @cvodes <Constants_link.html> CV_NO_FWD *)
  exception NoForwardCall

  (** Reinitialization of the forward problem failed at the first checkpoint
      (corresponding to the initial time of the forward problem).

      @cvodes <Constants_link.html> CV_REIFWD_FAIL *)
  exception ForwardReinitFailure

  (** An error occured when integrating the forward problem from a
      checkpoint.

      @cvodes <Constants_link.html> CV_FWD_FAIL *)
  exception ForwardFailure

  (** No backward problem has been created.

      @cvodes <Constants_link.html> CV_NO_BCK *)
  exception NoBackwardProblem

  (** The final time was outside the interval over which the forward
      problem was solved.

      @cvodes <Constants_link.html> CV_BAD_TB0 *)
  exception BadFinalTime

end (* }}} *)

