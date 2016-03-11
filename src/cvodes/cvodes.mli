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

open Cvode_impl
open Sundials

(** Alias for Cvode sessions. The Cvodes submodules all work by ‘extending’
    an existing session created with {!Cvode.init}. *)
type ('data, 'kind) session = ('data, 'kind) Cvode.session

(** Integration of pure quadrature equations.
 
    Adds a vector $y_Q$ of $N_Q$ quadrature variables defined by
    {% $\frac{\mathrm{d} y_Q}{\mathrm{d}t} = f_Q(t, y, p)$%}. The values of
    these variables are calculated more efficiently since they are excluded
    from the nonlinear solution stage.

    @cvodes <node5#SECTION00570000000000000000> Integration of pure quadrature equations
    @cvodes <node3#s:quad> Pure quadrature integration *)
module Quadrature :
  sig (* {{{ *)
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

        @cvodes <node5#ss:user_fct_quad> CVQuadRhsFn *)
    type 'd quadrhsfn = float -> 'd -> 'd -> unit

    (** Activates the integration of quadrature equations. The vector
        gives the initial value of $y_Q$.

        @cvodes <node5#ss:quad_malloc> CVodeQuadInit *)
    val init : ('d, 'k) Cvode.session -> 'd quadrhsfn
                  -> ('d, 'k) Nvector.t -> unit

    (** Reinitializes the integration of quadrature equations. The vector
        gives a new value for $y_Q$.

        @cvodes <node5#ss:quad_malloc> CVodeQuadReInit *)
    val reinit : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> unit

    (** Returns the quadrature solutions and time reached after a successful
        solver step. The given vector is filled with values calculated during
        either {!Cvode.solve_normal} or {!Cvode.solve_one_step} and the
        value of the independent variable is returned.

        @cvodes <node5#ss:quad_get> CVodeGetQuad *)
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

        @cvodes <node5#ss:quad_get> CVodeGetQuadDky
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

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
    val set_tolerances : ('d, 'k) Cvode.session -> ('d, 'k) tolerance -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the number of calls to the quadrature function.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals *)
    val get_num_rhs_evals       : ('d, 'k) Cvode.session -> int

    (** Returns the number of local error test failures that have occurred
        due to quadrature variables.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails *)
    val get_num_err_test_fails  : ('d, 'k) Cvode.session -> int

    (** Returns the quadrature error weights at the current time.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights *)
    val get_err_weights : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> unit

    (** Returns quadrature-related statistics. These are the
        number of calls to the quadrature function ([nfqevals]) and the number
        of error test failures due to quadrature variables ([nqetfails]).

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats
        @return ([nfqevals], [nqetfails]) *)
    val get_stats : ('d, 'k) Cvode.session -> int * int

    (** {2:exceptions Exceptions} *)

    (** Quadrature integration was not initialized.

        @cvodes <node5#ss:quad_get> CV_NO_QUAD *)
    exception QuadNotInitialized

    (** The quadrature function failed in an unrecoverable manner.

        @cvodes <node5#SECTION00572000000000000000> CV_QRHSFUNC_FAIL *)
    exception QuadRhsFuncFailure

    (** The quadrature function failed at the first call.

        @cvodes <node5#SECTION00572000000000000000> CV_FIRST_QRHSFUNC_ERR *)
    exception FirstQuadRhsFuncFailure

    (** Convergence test failures occurred too many times due to repeated
        recoverable errors in the quadrature function. Also raised if the
        quadrature function had repeated recoverable errors during the
        estimation of an initial step size (if quadrature variables are
        included in error tests).

        @cvodes <node5#SECTION00572000000000000000> CV_REPTD_QRHSFUNC_ERR *)
    exception RepeatedQuadRhsFuncFailure

    (** The quadrature function had a recoverable error, but no
        recovery was possible. This failure mode is rare, as it can occur only
        if the quadrature function fails recoverably after an error test failed
        while at order one.

        @cvodes <node5#SECTION00572000000000000000> CV_UNREC_QRHSFUNC_ERR *)
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

    @cvodes <node6#ss:forward_usage> Enhanced skeleton for sensitivity analysis
    @cvodes <node6#s:forward> Using CVODES for Forward Sensitivity Analysis
    @cvodes <node3#ss:fwd_sensi> Forward sensitivity analysis *)
module Sensitivity :
  sig (* {{{ *)
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

        tmp : 'd double;
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

        @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn *)
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

        @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn *)
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

            @cvodes <node6#ss:sensi_malloc> CVodeSensInit
            @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn *)

      | OneByOne of 'd sensrhsfn1 option
        (** Calculate sensitivity functions one parameter at a time. The
            argument [None] specifies an internal difference quotient routine.

            @cvodes <node6#ss:sensi_malloc> CVodeSensInit1
            @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn *)

    (** Specifies a sensitivity solution method.

        @cvodes <node6#ss:sensi_malloc> CVodeSensInit
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit1 *)
    type sens_method =
        Simultaneous
        (** Correct state and sensitivity variables at the same time.
            If [Newton] was selected as the nonlinear system
            solution method, this amounts to performing a modified Newton
            iteration on the combined nonlinear system.
            {cconst CV_SIMULTANEOUS} *)
      | Staggered
        (** The correction step for the sensitivity variables takes place at the
            same time for all sensitivity equations, but only after the
            correction of the state variables has converged and the state
            variables have passed the local error test.
            {cconst CV_STAGGERED} *)
      | Staggered1
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

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensParams *)
    type sens_params = {
        pvals  : Sundials.RealArray.t option;
        (** An array of $N_p$ parameters $p$ in $f(t, y, p)$.  If
            specified, this array is updated by the solver to pass
            perturbed parameter values to the original right-hand side
            and root functions.  Those functions must (re-)read
            parameter values from this array on every invocation. *)
        pbar   : Sundials.RealArray.t option;
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
        SStolerances of float * Sundials.RealArray.t
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('d, 'k) Nvector.t array
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
      | EEtolerances
        (** Calculate the integration tolerances for sensitivities
            from those for state variables and the scaling factors
            in {{!sens_params}pbar}. *)

    (** Activates the calculation of forward sensitivities. The call
        [init s tol sm sp fs ys0], has as arguments:
        - [s], a session created with {!Cvode.init},
        - [tol], the tolerances desired,
        - [sens_method], the solution method,
        - [~sens_params], the parameter information (see {!sens_params}),
        - [~fs], the sensitivity function, and,
        - [ys0], initial values of the sensitivities for each parameter.

        If [~fs] is not given, an internal difference quotient routine
        is used.  In that case, or if the internal difference quotient
        routine will be specified in a subsequent call to
        {!Sensitivity.Quadrature.init}, then [sens_params] must be
        given with [pvals] set to non-[None].

        @cvodes <node6#ss:sensi_malloc> CVodeSensInit
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit1
        @cvodes <node6#ss:sens_optional_input> CVodeSetSensParams
        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSStolerances
        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSVtolerances
        @cvodes <node6#sss:cvfwdtolerances> CVodeSensEEtolerances *)
    val init : ('d, 'k) Cvode.session
               -> ('d, 'k) tolerance
               -> sens_method
               -> ?sens_params:sens_params
               -> 'd sensrhsfn
               -> ('d, 'k) Nvector.t array
               -> unit

    (** Reinitializes the forward sensitivity computation.

        @cvodes <node6#ss:sensi_malloc> CVodeSensReInit *)
    val reinit : ('d, 'k) Cvode.session -> sens_method
                      -> ('d, 'k) Nvector.t array -> unit

    (** Deactivates forward sensitivity calculations without deallocating
        memory. Sensitivities can be reactivated with {!reinit}.

        @cvodes <node6#ss:sensi_malloc> CVodeSensToggleOff *)
    val toggle_off : ('d, 'k) Cvode.session -> unit

    (** Support for quadrature sensitivity equations.

        Adds an additional vector {% $s_\mathit{Q}$%} of {% $N_\mathit{Q}$%}
        quadrature sensitivities
        {% $\frac{\mathrm{d} s_\mathit{Q}}{\mathrm{d}t}
                = f_{\mathit{QS}}(t, y, s, \dot{y}_Q, p)$%}.
        This mechanism allows, in particular, the calculation of the
        sensitivities of the ‘pure’ {!Cvodes.Quadrature} variables, $y_Q$.

        @cvodes <node3#SECTION00364000000000000000> Quadratures depending on forward sensitivities
        @cvodes <node6#SECTION00640000000000000000> Integration of quadrature equations depending on forward sensitivities *)
    module Quadrature :
      sig (* {{{ *)
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

            tmp : 'd double;
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

           @cvodes <node6#ss:user_fct_quad_sens> CVodeQuadSensRhsFn *)
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

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensInit
            @raise QuadNotInitialized {!Quadrature.init} has not been called. *)
        val init : ('d, 'k) Cvode.session -> ?fqs:'d quadsensrhsfn
                 -> ('d, 'k) Nvector.t array -> unit

        (** Reinitializes the quadrature sensitivity integration.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensReInit *)
        val reinit : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array -> unit

        (** {2:tols Tolerance specification} *)

        (** Tolerances for calculating quadrature sensitivities. *)
        type ('d, 'k) tolerance =
            NoStepSizeControl
            (** Quadrature variables are not used for step-size control
                (the default). *)
          | SStolerances of float * Sundials.RealArray.t
            (** [(rel, abs)] : scalar relative and absolute tolerances. *)
          | SVtolerances of float * ('d, 'k) Nvector.t array
            (** [(rel, abs)] : scalar relative and vector absolute
                tolerances. *)
          | EEtolerances
            (** Calculate the integration tolerances for the
                quadrature sensitivities from those provided for
                the pure quadrature variables. *)

        (** Specify how to use quadrature sensitivities in step size control.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeSetQuadSensErrCon
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSStolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSVtolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensEEtolerances *)
        val set_tolerances : ('d, 'k) Cvode.session
                              -> ('d, 'k) tolerance -> unit

        (** {2:quadout Output functions} *)

        (** Returns the quadrature sensitivity solutions and time reached
            after a successful solver step. The given vectors are filled with
            values calculated during either {!Cvode.solve_normal} or
            {!Cvode.solve_one_step} and the value of the independent variable
            is returned.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens *)
        val get : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array -> float

        (** Returns a single quadrature sensitivity vector after a successful
            solver step. Like {!get}, but the argument [i] specifies a specific
            vector to return.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens1
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

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky
            @raise BadIS The index is not in the allowed range.
            @raise BadK [k] is not in the range 0, 1, ..., [qlast].
            @raise BadT [t] is not in the allowed range. *)
        val get_dky : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array
                        -> float -> int -> unit

        (** Returns the interpolated solution or derivatives of a single
            quadrature sensitivity solution vector.
            [get_dky s dksq t k i] is like
            {!get_dky} but restricted to the [i]th sensitivity solution vector.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky1
            @raise BadK [k] is not in the range 0, 1, ..., [qlast].
            @raise BadT [t] is not in the allowed range. *)
        val get_dky1 : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t
                         -> float -> int -> int -> unit

        (** {2:get Querying the solver (optional output functions)} *)

        (** Returns the number of calls to the quadrature right-hand side
            function.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumRhsEvals *)
        val get_num_rhs_evals       : ('d, 'k) Cvode.session -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumErrTestFails *)
        val get_num_err_test_fails  : ('d, 'k) Cvode.session -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensErrWeights *)
        val get_err_weights
              : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array -> unit

        (** Returns quadrature-related statistics. These are the
            number of calls to the quadrature function ([nfqevals]) and the
            number of error test failures due to quadrature variables
            ([nqetfails]).
         
            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensStats
            @return ([nfqevals], [nqetfails]) *)
        val get_stats : ('d, 'k) Cvode.session -> int * int

        (** {2:exceptions Exceptions} *)

        (** Quadrature integration was not initialized.

            @cvodes <node5#SECTION00642000000000000000> CV_NO_QUAD_SENS *)
        exception QuadSensNotInitialized

        (** The sensitivity quadrature function failed in an unrecoverable
            manner.

            @cvodes <node6#SECTION00642000000000000000> CV_QSRHSFUNC_FAIL *)
        exception QuadSensRhsFuncFailure

        (** The sensitivity quadrature function failed at the first call.

            @cvodes <node6#SECTION00642000000000000000> CV_FIRST_QSRHSFUNC_ERR *)
        exception FirstQuadSensRhsFuncFailure

        (** Convergence test failures occurred too many times due to repeated
            recoverable errors in the quadrature function. Also raised if the
            sensitivity quadrature function had repeated recoverable errors
            during the estimation of an initial step size (if quadrature
            variables are included in error tests).
          
            @cvodes <node6#SECTION00642000000000000000> CV_REPTD_QSRHSFUNC_ERR *)
        exception RepeatedQuadSensRhsFuncFailure

        (** The sensitivity quadrature function had a recoverable error, but no
            recovery was possible. This failure mode is rare, as it can occur
            only if the quadrature function fails recoverably after an error
            test failed while at order one.
          
            @cvodes <node6#SECTION00642000000000000000> CV_UNREC_QSRHSFUNC_ERR *)
        exception UnrecoverableQuadSensRhsFuncFailure
      end (* }}} *)

    (** {2:sensout Output functions} *)

    (** Returns the sensitivity solution vectors after a successful solver
        step. The given vectors are filled with values calculated during
        either {!Cvode.solve_normal} or {!Cvode.solve_one_step} and the
        value of the independent variable is returned.

        @cvodes <node6#ss:sensi_get> CVodeGetSens *)
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

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky
        @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
        @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
    val get_dky : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t array
                    -> float -> int -> unit

    (** Returns a single sensitivity solution vector after a successful solver
        step. The given vector is filled with values calculated for the [i]th
        sensitivity vector during either {!Cvode.solve_normal} or
        {!Cvode.solve_one_step} and the value of the independent variable is
        returned.

        @cvodes <node6#ss:sensi_get> CVodeGetSens1
        @raise BadIS The index [i] is not in the allowed range. *)
    val get1 : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t -> int -> float

    (** Returns the interpolated solution or derivatives of a single
        sensitivity solution vector. [get_dky s dks t k i] is like {!get_dky}
        but restricted to the [i]th sensitivity solution vector.

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky1
        @raise BadIS The index [i] is not in the allowed range.
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky1 : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t
                     -> float -> int -> int -> unit

    (** {2:set Modifying the solver (optional input functions)} *)

    (** Sets the integration tolerances for sensitivities.

        {b NB}: Unlike the other [set_tolerances] functions in [Cvodes], this
        one does {b not} call {!set_err_con} (which defaults to [false]).

        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSStolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensSVtolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensEEtolerances *)
    val set_tolerances : ('d, 'k) Cvode.session -> ('d, 'k) tolerance -> unit

    (** Sets whether sensitivity variables are used in the error control
        mechanism. The default is [false].

        @cvodes <node5#ss:sens_optional_input> CVodeSetSensErrCon *)
    val set_err_con : ('d, 'k) Cvode.session -> bool -> unit

    (** A difference quotient strategy. See {!set_dq_method}. *)
    type dq_method = DQCentered (** {cconst CV_CENTERED} *)
                   | DQForward  (** {cconst CV_FORWARD} *)

    (** Sets the difference quotient strategy when sensitivity equations
        are computed internally by the solver rather than via callback. A
        method and [dqrhomax] value must be given. The latter determines when
        to switch between simultaneous or separate approximations of the two
        terms in the sensitivity right-hand side.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensDQMethod
        @cvodes <node3#ss:fwd_sensi> Forward Sensitivity Analysis *)
    val set_dq_method : ('d, 'k) Cvode.session -> dq_method -> float -> unit

    (** Sets the maximum number of nonlinear solver iterations for
        sensitivity variables permitted per step.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensMaxNonlinIters *)
    val set_max_nonlin_iters : ('d, 'k) Cvode.session -> int -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the number of calls to the sensitivity function.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumRhsEvals *)
    val get_num_rhs_evals       : ('d, 'k) Cvode.session -> int

    (** Returns the number of calls to the right-hand side function due
        to the internal finite difference approximation of the sensitivity
        equations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetNumRhsEvalsSens *)
    val get_num_rhs_evals_sens  : ('d, 'k) Cvode.session -> int

    (** Returns the number of local error test failures for the sensitivity
        variables that have occurred.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumErrTestFails *)
    val get_num_err_test_fails  : ('d, 'k) Cvode.session -> int

    (** Returns the number of calls made to the linear solver's setup function
        due to forward sensitivity calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumLinSolvSetups *)
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

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensStats *)
    val get_stats : ('d, 'k) Cvode.session -> sensitivity_stats

    (** Returns the sensitivity error weights at the current time.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensErrWeights
        @cvodes <node3#e:errwt> IVP solution (W_i, Eq. (2.7)) *)
    val get_err_weights : ('d, 'k) Cvode.session
                            -> ('d, 'k) Nvector.t array -> unit

    (** Returns the number of nonlinear iterations performed for sensitivity
        calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvIters *)
    val get_num_nonlin_solv_iters : ('d, 'k) Cvode.session -> int

    (** Returns the number of nonlinear convergence failures that have occurred
        during sensitivity calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvConvFails *)
    val get_num_nonlin_solv_conv_fails : ('d, 'k) Cvode.session -> int

    (** Returns both the numbers of nonlinear iterations performed [nniters]
        and nonlinear convergence failures [nncfails] during sensitivity
        calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNonlinSolvStats
        @return ([nniters], [nncfails]) *)
    val get_nonlin_solv_stats : ('d, 'k) Cvode.session -> int * int

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed for each sensitivity equation separately in the [Staggered1]
        case.

        @cvodes <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvIters *)
    val get_num_stgr_nonlin_solv_iters : ('d, 'k) Cvode.session
                                         -> Sundials.LintArray.t -> unit

    (** Returns the number of nonlinear convergence failures that have occurred
        for each sensitivity equation separately in the [Staggered1] case.

        @cvodes <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvConvFails *)
    val get_num_stgr_nonlin_solv_conv_fails : ('d, 'k) Cvode.session
                                              -> Sundials.LintArray.t -> unit

    (** {2:exceptions Exceptions} *)

    (** Sensitivity analysis was not initialized.

        @cvodes <node5#ss:sensi_get> CV_NO_SENS *)
    exception SensNotInitialized

    (** The sensitivity function failed in an unrecoverable manner.

        @cvodes <node6#SECTION00623000000000000000> CV_SRHSFUNC_FAIL *)
    exception SensRhsFuncFailure

    (** The sensitivity function had a recoverable error when first called.

        @cvodes <node6#SECTION00623000000000000000> CV_FIRST_SRHSFUNC_ERR *)
    exception FirstSensRhsFuncFailure

    (** Too many convergence test failures, or unable to estimate the initial
        step size, due to repeated recoverable errors in the sensitivity
        function.

        @cvodes <node6#SECTION00623000000000000000> CV_REPTD_SRHSFUNC_ERR *)
    exception RepeatedSensRhsFuncFailure

    (** The sensitivity function had a recoverable error, but no recovery was
        possible. This error can only occur after an error test failure at
        order one.
     
      @cvodes <node6#SECTION00623000000000000000> CV_UNREC_SRHSFUNC_ERR *)
    exception UnrecoverableSensRhsFuncFailure

    (** The index passed to identify a particular sensitivity is invalid.

        @cvodes <node6> CV_BAD_IS *)
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
  
    @cvodes <node7#ss:skeleton_adj> Enhanced Skeleton for Adjoint Sensitivity Analysis
    @cvodes <node7#s:adjoint> Using CVODES for Adjoint Sensitivity Analysis
    @cvodes <node3#ss:adj_sensi> Adjoint sensitivity analysis *)
module Adjoint :
  sig (* {{{ *)
    (** A backward session with the CVODES solver. Multiple backward sessions
        may be associated with a single parent session.

        @cvodes <node7#sss:cvinitb> Backward problem initialization functions *)
    type ('data, 'kind) bsession = ('data, 'kind) AdjointTypes.bsession

    (** Alias for backward sessions based on serial nvectors. *)
    type 'kind serial_bsession = (Nvector_serial.data, 'kind) bsession
                                 constraint 'kind = [>Nvector_serial.kind]

    (** {2:fwd Forward solution} *)

    (** Specifies the type of interpolation to use between checkpoints.

        @cvodes <node3#ss:checkpointing> Checkpointing scheme *)
    type interpolation = IPolynomial (** {cconst CV_POLYNOMIAL} *)
                       | IHermite    (** {cconst CV_HERMITE} *)

    (** Activates the forward-backward problem. The arguments specify the number
        of integration steps between consecutive checkpoints, and the type of
        variable-degree interpolation.

        @cvodes <node7#sss:cvadjinit> CVodeAdjInit *)
    val init : ('d, 'k) Cvode.session -> int -> interpolation -> unit

    (** Integrates the forward problem over an interval and saves
        checkpointing data. The arguments are the next time at which a solution
        is desired ([tout]) and a vector to receive the computed result
        ([y]). The function returns a triple [tret, ncheck, sr]: the time
        reached by the solver, the cumulative number of checkpoints stored, and
        whether [tout] was reached. The solver takes internal steps until it
        has reached or just passed the [tout] parameter {cconst CV_NORMAL},
        it then interpolates to approximate [y(tout)].

        @cvodes <node7#sss:cvsolvef> CVodeF
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

        @cvodes <node7#sss:cvsolvef> CVodeF
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

        @cvodes <node7#sss:lin_solv_b> Linear Solver Initialization Functions *)
    type ('data, 'kind) linear_solver =
            ('data, 'kind) AdjointTypes.linear_solver

    (** Alias for linear solvers that are restricted to serial nvectors. *)
    type 'kind serial_linear_solver = (Nvector_serial.data, 'kind) linear_solver
                                      constraint 'kind = [>Nvector_serial.kind]

    (** Workspaces with three temporary vectors. *)
    type 'd triple = 'd * 'd * 'd

    (** Arguments common to Jacobian callback functions.

        @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB
        @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB 
        @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
        @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
        @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
    type ('t, 'd) jacobian_arg = ('t, 'd) Cvode_impl.AdjointTypes.jacobian_arg =
      {
        jac_t   : float;        (** The independent variable. *)
        jac_y   : 'd;           (** The forward solution vector. *)
        jac_yb  : 'd;           (** The backward solution vector. *)
        jac_fyb : 'd;           (** The backward right-hand side function [fB]. *)
        jac_tmp : 't;           (** Temporary storage vectors. *)
      }

    (** The range of nonzero entries in a band matrix. *)
    type bandrange = Cvode_impl.bandrange =
      { mupper : int; (** The upper half-bandwidth.  *)
        mlower : int; (** The lower half-bandwidth.  *) }

    (** Diagonal approximation of Jacobians by difference quotients. *)
    module Diag :
      sig (* {{{ *)
        (** A linear solver based on Jacobian approximation by difference
            quotients.

            @cvodes <node7#sss:lin_solv_b> CVDiagB *)
        val solver : ('data, 'kind) linear_solver

        (** Returns the sizes of the real and integer workspaces used by the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetWorkSpace
            @return ([real_size], [integer_size]) *)
        val get_work_space : ('d, 'k) bsession -> int * int

        (** Returns the number of calls made to the right-hand side
            function due to finite difference Jacobian approximation in the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetNumRhsEvals *)
        val get_num_rhs_evals : ('d, 'k) bsession -> int
      end (* }}} *)

    (** Direct Linear Solvers operating on dense and banded matrices. *)
    module Dls :
      sig (* {{{ *)

        (** Callback functions that compute dense approximations to a Jacobian
            matrix without forward sensitivities. In the call
            [dense_jac_fn arg jac], [arg] is a {!jacobian_arg} with three work
            vectors and the computed Jacobian must be stored in [jac].

            The callback should load the [(i,j)]th entry of [jac] with
            {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
            the [i]th equation with respect to the [j]th variable, evaluated
            at the values of [t] and [y] obtained from [arg]. Only nonzero
            elements need be loaded into [jac].

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning Neither the elements of [arg] nor the matrix [jac] should
                     be accessed after the function has returned.}

            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
        type dense_jac_fn_no_sens
          = (RealArray.t triple, RealArray.t) jacobian_arg
            -> Dls.DenseMatrix.t
            -> unit

        (** Callback functions that compute dense approximations to a Jacobian
            matrix with forward sensitivities. In the call
            [dense_jac_fn arg s jac], [arg] is a {!jacobian_arg} with three work
            vectors, [s] is an array of forward sensitivity vectors, and the
            computed Jacobian must be stored in [jac].

            The callback should load the [(i,j)]th entry of [jac] with
            {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
            the [i]th equation with respect to the [j]th variable, evaluated
            at the values of [t] and [y] obtained from [arg]. Only nonzero
            elements need be loaded into [jac].

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning Neither the elements of [arg], [s] nor the matrix [jac]
                     should be accessed after the function has returned.}

            @nocvodes <node7#ss:densejac_bs> CVDlsDenseJacFnBS *)
        type dense_jac_fn_with_sens
          = (RealArray.t triple, RealArray.t) jacobian_arg
            -> RealArray.t array
            -> Dls.DenseMatrix.t
            -> unit

        (** Callback functions that compute dense approximations to a Jacobian
            matrix.

            @cvodes <node5#ss:sjacFnB> CVDlsDenseJacFnB
            @nocvodes <node5#ss:sjacFnBS> CVDlsDenseJacFnBS *)
        type dense_jac_fn =
            DenseNoSens of dense_jac_fn_no_sens
            (** Does not depend on forward sensitivities. *)
          | DenseWithSens of dense_jac_fn_with_sens
            (** Depends on forward sensitivities. *)

        (** A direct linear solver on dense matrices. The optional argument
            specifies a callback function for computing an approximation to the
            Jacobian matrix. If this argument is omitted, then a default
            implementation based on difference quotients is used.

            @cvodes <node7#sss:lin_solv_b> CVDenseB
            @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @nocvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnBS
            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB
            @nocvodes <node7#ss:densejac_bs> CVDlsDenseJacFnBS *)
        val dense : ?jac:dense_jac_fn -> unit -> 'k serial_linear_solver

        (** A direct linear solver on dense matrices using LAPACK. See {!dense}.
            Only available if {!Sundials.lapack_enabled}.

            @raise Sundials.NotImplementedBySundialsVersion Solver not available.
            @cvodes <node7#sss:lin_solv_b> CVLapackDenseB
            @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
        val lapack_dense : ?jac:dense_jac_fn -> unit -> 'k serial_linear_solver

        (** Callback functions that compute banded approximations to
            a Jacobian matrix without forward sensitivities. In the call
            [band_jac_fn {mupper; mlower} arg jac],
            - [mupper] is the upper half-bandwidth of the Jacobian,
            - [mlower] is the lower half-bandwidth of the Jacobian,
            - [arg] is a {!jacobian_arg} with three work vectors, and,
            - [jac] is storage for the computed Jacobian.

            The callback should load the [(i,j)]th entry of [jac] with
            {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
            the [i]th equation with respect to the [j]th variable, evaluated
            at the values of [t] and [y] obtained from [arg]. Only nonzero
            elements need be loaded into [jac].

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning Neither the elements of [arg] nor the matrix [jac] should
                     be accessed after the function has returned.}

            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        type band_jac_fn_no_sens
          = bandrange
            -> (RealArray.t triple, RealArray.t) jacobian_arg
            -> Dls.BandMatrix.t
            -> unit

        (** Callback functions that compute banded approximations to
            a Jacobian matrix with forward sensitivities. In the call
            [band_jac_fn {mupper; mlower} arg s jac],
            - [mupper] is the upper half-bandwidth of the Jacobian,
            - [mlower] is the lower half-bandwidth of the Jacobian,
            - [arg] is a {!jacobian_arg} with three work vectors,
            - [s] is an array of forward sensitivity vectors, and
            - [jac] is storage for the computed Jacobian.

            The callback should load the [(i,j)]th entry of [jac] with
            {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
            the [i]th equation with respect to the [j]th variable, evaluated
            at the values of [t] and [y] obtained from [arg]. Only nonzero
            elements need be loaded into [jac].

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning Neither the elements of [arg], [s], nor the matrix [jac]
                     should be accessed after the function has returned.}

            @nocvodes <node7#ss:bandjac_bs> CVDlsBandJacFnBS *)
        type band_jac_fn_with_sens
          = bandrange
            -> (RealArray.t triple, RealArray.t) jacobian_arg
            -> RealArray.t array
            -> Dls.BandMatrix.t
            -> unit

        (** Callback functions that compute banded approximations to a Jacobian
            matrix.

            @cvodes <node5#ss:bandjac_b> CVDlsBandJacFnB
            @nocvodes <node5#ss:bandjac_bs> CVDlsBandJacFnBS *)
        type band_jac_fn =
            BandNoSens of band_jac_fn_no_sens
            (** Does not depend on forward sensitivities. *)
          | BandWithSens of band_jac_fn_with_sens
            (** Depends on forward sensitivities. *)

        (** A direct linear solver on banded matrices. The optional argument
            specifies a callback function for computing an approximation to the
            Jacobian matrix. If this argument is omitted, then a default
            implementation based on difference quotients is used. The other
            argument gives the width of the bandrange.

            @cvodes <node7#sss:lin_solv_b> CVBandB
            @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        val band : ?jac:band_jac_fn -> bandrange -> 'k serial_linear_solver

        (** A direct linear solver on banded matrices using LAPACK. See {!band}.
            Only available if {!Sundials.lapack_enabled}.

            @raise Sundials.NotImplementedBySundialsVersion Solver not available.
            @cvodes <node7#sss:lin_solv_b> CVLapackBandB
            @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        val lapack_band
              : ?jac:band_jac_fn -> bandrange -> 'k serial_linear_solver

        (** {3:stats Solver statistics} *)

        (** Returns the sizes of the real and integer workspaces used by a direct
            linear solver.

            @cvode <node5#sss:optout_dls> CVDlsGetWorkSpace
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
            @return ([real_size], [integer_size]) *)
        val get_work_space : 'k serial_bsession -> int * int

        (** Returns the number of calls made by a direct linear solver to the
            Jacobian approximation function.

            @cvode <node5#sss:optout_dls> CVDlsGetNumJacEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_jac_evals : 'k serial_bsession -> int

        (** Returns the number of calls to the right-hand side callback due to
            the finite difference Jacobian approximation.

            @cvode <node5#sss:optout_dls> CVDlsGetNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_rhs_evals : 'k serial_bsession -> int

        (** {3:lowlevel Low-level solver manipulation}

            The {!init} and {!reinit} functions are the preferred way to set or
            change a Jacobian function. These low-level functions are provided for
            experts who want to avoid resetting internal counters and other
            associated side-effects. *)

        (** Change the dense Jacobian function.
       
            @cvode <node5#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @nocvode <node5#SECTION00728200000000000000> CVDlsSetDenseJacFnBS *)
        val set_dense_jac_fn : 'k serial_bsession -> dense_jac_fn -> unit

        (** Remove a dense Jacobian function and use the default
            implementation.

            @cvode <node5#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @nocvode <node5#SECTION00728200000000000000> CVDlsSetDenseJacFnBS *)
        val clear_dense_jac_fn : 'k serial_bsession -> unit

        (** Change the band Jacobian function.

            @cvode <node5#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @nocvode <node5#SECTION00728300000000000000> CVDlsSetBandJacFnBS *)
        val set_band_jac_fn : 'k serial_bsession -> band_jac_fn -> unit

        (** Remove a banded Jacobian function and use the default
            implementation.

            @cvode <node5#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @nocvode <node5#SECTION00728300000000000000> CVDlsSetBandJacFnBS *)
        val clear_band_jac_fn : 'k serial_bsession -> unit
      end (* }}} *)

    (** Sparse Linear Solvers.

        @nocvodes <node> The SLS modules *)
    module Sls :
      sig (* {{{ *)
        (** Callback functions that compute sparse approximations to a Jacobian
            matrix without forward sensitivites. In the call [sparse_jac_fn arg
            jac], [arg] is a {!Cvodes.Adjoint.jacobian_arg} with three work
            vectors and the computed Jacobian must be stored in [jac].

            The callback should load the [(i,j)]th entry of [jac] with
            {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
            the [i]th equation with respect to the [j]th variable, evaluated at
            the values of [t] and [y] obtained from [arg]. Only nonzero elements
            need be loaded into [jac].

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning Neither the elements of [arg] nor the matrix [jac] should
                     be accessed after the function has returned.}

            @nocvodes <node5#ss:sjacFnB> CVSlsSparseJacFnB *)
        type sparse_jac_fn_no_sens =
          (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg
          -> Sls.SparseMatrix.t -> unit

        (** Callback functions that compute sparse approximations to a Jacobian
            matrix with forward sensitivities. In the call [sparse_jac_fn arg s
            jac], [arg] is a {!Cvodes.Adjoint.jacobian_arg} with three work
            vectors, [s] is an array of forward sensitivity vectors, and the
            computed Jacobian must be stored in [jac].

            The callback should load the [(i,j)]th entry of [jac] with
            {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
            the [i]th equation with respect to the [j]th variable, evaluated at
            the values of [t] and [y] obtained from [arg]. Only nonzero elements
            need be loaded into [jac].

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning Neither the elements of [arg] nor the matrix [jac] should
                     be accessed after the function has returned.}

            @nocvodes <node5#ss:sjacFnBS> CVSlsSparseJacFnBS *)
        type sparse_jac_fn_with_sens =
          (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg
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

        (** KLU sparse-direct linear solver module for CVODES (requires KLU).

            @nocvodes <node5#sss:cvklu> The KLU Solver *)
        module Klu : sig (* {{{ *)

          (** A direct linear solver on sparse matrices. In the call,
              [klu jfn nnz], [jfn] is a callback function that computes an
              approximation to the Jacobian matrix and [nnz] is the maximum
              number of nonzero entries in that matrix.

              @raise Sundials.NotImplementedBySundialsVersion Solver not available.
              @nocvodes <node5#sss:lin_solv_init> CVKLUB
              @nocvodes <node5#sss:optin_sls> CVSlsSetSparseJacFnB
              @nocvodes <node5#sss:optin_sls> CVSlsSetSparseJacFnBS
              @nocvodes <node5#ss:sjacFnB> CVSlsSparseJacFnB
              @nocvodes <node5#ss:sjacFnBS> CVSlsSparseJacFnBS *)
          val solver : sparse_jac_fn -> int -> 'k serial_linear_solver

          (** The ordering algorithm used for reducing fill. *)
          type ordering = Cvode.Sls.Klu.ordering =
               Amd      (** Approximate minimum degree permutation. *)
             | ColAmd   (** Column approximate minimum degree permutation. *)
             | Natural  (** Natural ordering. *)

          (** Sets the ordering algorithm used to minimize fill-in.

              @nocvodes <node5#ss:sls_optin> CVKLUSetOrdering
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val set_ordering : 'k serial_bsession -> ordering -> unit

          (** Reinitializes the Jacobian matrix memory and flags.
              In the call, [reinit s n nnz realloc], [n] is the number of system
              state variables, and [nnz] is the number of non-zeroes in the
              Jacobian matrix. New symbolic and numeric factorizations will be
              completed at the next solver step. If [realloc] is true, the
              Jacobian matrix will be reallocated based on [nnz].

              @nocvodes <node5#ss:sls_optin> CVKLUReInit
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val reinit : 'k serial_bsession -> int -> int -> bool -> unit

          (** Returns the number of calls made by a sparse linear solver to the
              Jacobian approximation function.

              @nocvodes <node5#sss:optout_sls> CVSlsGetNumJacEvals
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val get_num_jac_evals : 'k serial_bsession -> int

        end (* }}} *)

        (** SuperLU_MT sparse-direct linear solver module for CVODES
            (requires SuperLU_MT).

            @nocvodes <node5#sss:cvsuperlumt> The SuperLUMT Solver *)
        module Superlumt : sig (* {{{ *)

          (** A direct linear solver on sparse matrices. In the call,
              [superlumt jfn nnz], [jfn] specifies a callback function that
              computes an approximation to the Jacobian matrix and [nnz] is the
              maximum number of nonzero entries in that matrix.

              @raise Sundials.NotImplementedBySundialsVersion Solver not available.
              @nocvodes <node5#sss:lin_solv_init> CVSuperLUMTB
              @nocvodes <node5#sss:optin_sls> CVSlsSetSparseJacFnB
              @nocvodes <node5#sss:optin_sls> CVSlsSetSparseJacFnBS
              @nocvodes <node5#ss:sjacFnB> CVSlsSparseJacFnB
              @nocvodes <node5#ss:sjacFnBS> CVSlsSparseJacFnBS *)
          val solver : sparse_jac_fn -> nnz:int -> nthreads:int
                            -> 'k serial_linear_solver

          (** The ordering algorithm used for reducing fill. *)
          type ordering = Cvode.Sls.Superlumt.ordering =
               Natural       (** Natural ordering. *)
             | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
             | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
             | ColAmd        (** Column approximate minimum degree permutation. *)

          (** Sets the ordering algorithm used to minimize fill-in.

              @nocvodes <node5#ss:sls_optin> CVSuperLUMTSetOrdering
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val set_ordering : 'k serial_bsession -> ordering -> unit

          (** Returns the number of calls made by a sparse linear solver to the
              Jacobian approximation function.

              @nocvodes <node5#sss:optout_sls> CVSlsGetNumJacEvals
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val get_num_jac_evals : 'k serial_bsession -> int

        end (* }}} *)
      end (* }}} *)

    (** Scaled Preconditioned Iterative Linear Solvers.

        @cvodes <node7#ss:optional_output_b> Optional output functions for the backward problem.
        @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
        @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
    module Spils :
      sig (* {{{ *)
        (** Arguments passed to the preconditioner solver function.

            @cvode <node7#ss:psolve_b> CVSpilsPrecSolveFnB *)
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
            - [jac] is a {!jacobian_arg} with one work vector,
            - [arg] is {!prec_solve_arg} that specifies the linear system, and
            - [z] is computed to solve {% $P\mathtt{z} = \mathtt{arg.rhs}$%}.
            $P$ is a preconditioner matrix, which approximates, however crudely,
            the Newton matrix {% $M = I - \gamma J$%} where
            {% $J = \frac{\partial f}{\partial y}$%}.

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning The elements of [jac], [arg], and [z] should not be
            accessed after the function has returned.}

            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB *)
        type 'd prec_solve_fn =
          ('d, 'd) jacobian_arg
          -> 'd prec_solve_arg
          -> 'd
          -> unit

        (** Callback functions that solve a linear system involving a
            preconditioner matrix with forward sensitivities.
            In the call [prec_solve_fn jac arg s z],
            - [jac] is a {!jacobian_arg} with one work vector,
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

            @nocvodes <node7#ss:psolve_bs> CVSpilsPrecSolveFnBS *)
        type 'd prec_solve_fn_with_sens =
          ('d, 'd) jacobian_arg
          -> 'd prec_solve_arg
          -> 'd array
          -> 'd
          -> unit

        (** Callback functions that preprocess or evaluate Jacobian-related data
            needed by {!prec_solve_fn} without forward sensitivities.
            In the call [prec_setup_fn jac s jok gamma],
            - [jac] is a {!jacobian_arg} with three work vectors,
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

            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
        type 'd prec_setup_fn =
          ('d triple, 'd) jacobian_arg
          -> bool
          -> float
          -> bool

        (** Callback functions that preprocess or evaluate Jacobian-related data
            needed by {!prec_solve_fn} with forward sensitivities.
            In the call [prec_setup_fn jac s jok gamma],
            - [jac] is a {!jacobian_arg} with three work vectors,
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

            @nocvodes <node7#ss:psetup_bs> CVSpilsPrecSetupFnBS *)
        type 'd prec_setup_fn_with_sens =
          ('d triple, 'd) jacobian_arg
          -> 'd array
          -> bool
          -> float
          -> bool

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

            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB *)
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

            @nocvodes <node7#ss:jtimesv_bs> CVSpilsJacTimesVecFnBS *)
        type 'd jac_times_vec_fn_with_sens =
          ('d, 'd) jacobian_arg
          -> 'd array
          -> 'd
          -> 'd
          -> unit

        (** Callback functions that compute the Jacobian times a vector.

            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
            @nocvodes <node7#ss:jtimesv_bs> CVSpilsJacTimesVecFnBS *)
        type 'd jac_times_vec_fn =
          | NoSens of 'd jac_times_vec_fn_no_sens
            (** Does not depend on forward sensitivities. *)
          | WithSens of 'd jac_times_vec_fn_with_sens
            (** Depends on forward sensitivities. *)

        (** Specifies a preconditioner, including the type of preconditioning
            (none, left, right, or both) and callback functions. The following
            functions and those in {!Banded} and {!Cvodes_bbd} construct
            preconditioners.

            The {!prec_solve_fn} is mandatory. The {!prec_setup_fn} can
            be omitted if not needed.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB
            @nocvodes <node7#ss:psolve_bs> CVSpilsPrecSolveFnBS
            @nocvodes <node7#ss:psetup_bs> CVSpilsPrecSetupFnBS *)
        type ('d, 'k) preconditioner =
          ('d, 'k) AdjointTypes.SpilsTypes.preconditioner

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
        module Banded : sig

          (** A band matrix {!preconditioner} based on difference quotients.
              The call [prec_left br] instantiates a left preconditioner which
              generates a banded approximation to the Jacobian with [br.mlower]
              sub-diagonals and [br.mupper] super-diagonals.

              @cvode <node7#SECTION00741000000000000000> CVBandPrecInitB *)
          val prec_left : bandrange -> (Nvector_serial.data,
                                        [>Nvector_serial.kind]) preconditioner

          (** Like {!prec_left} but preconditions from the right.

              @cvode <node7#SECTION00741000000000000000> CVBandPrecInitB *)
          val prec_right : bandrange -> (Nvector_serial.data,
                                         [>Nvector_serial.kind]) preconditioner

          (** Like {!prec_left} but preconditions from both sides.

              @cvode <node7#SECTION00741000000000000000> CVBandPrecInitB *)
          val prec_both : bandrange -> (Nvector_serial.data,
                                        [>Nvector_serial.kind]) preconditioner

          (** {4:stats Banded statistics} *)

          (** Returns the sizes of the real and integer workspaces
              used by the banded preconditioner module.

              @cvodes <node5#sss:cvbandpre> CVBandPrecGetWorkSpace
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
              @return ([real_size], [integer_size]) *)
          val get_work_space : 'k serial_bsession -> int * int

          (** Returns the number of calls to the right-hand side callback for the
              difference banded Jacobian approximation. This counter is only updated
              if the default difference quotient function is used.

              @cvodes <node5#sss:cvbandpre> CVBandPrecGetNumRhsEvals
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val get_num_rhs_evals : 'k serial_bsession -> int
        end

        (** {3:lsolvers Solvers} *)

        (** Krylov iterative solver using the scaled preconditioned generalized
            minimum residual (GMRES) method.
            In the call [spgmr ~maxl:maxl ~jac_times_vec:jtv prec],
            - [maxl] is the maximum dimension of the Krylov subspace
                     (defaults to 5),
            - [jtv] computes an approximation to the product between the
                    Jacobian matrix and a vector, and
            - [prec] is a {!preconditioner}.
            
            If the {!jac_times_vec_fn} is omitted, a default implementation
            based on difference quotients is used.

            @cvodes <node7#sss:lin_solv_b> CVSpgmrB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @nocvodes <node7> CVSpilsSetPreconditionerBS
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @nocvodes <node7> CVSpilsSetJacTimesVecFnBS *)
        val spgmr :
          ?maxl:int
          -> ?jac_times_vec:'d jac_times_vec_fn
          -> ('d, 'k) preconditioner
          -> ('d, 'k) linear_solver

        (** Krylov iterative solver using the scaled preconditioned biconjugate
            stabilized (Bi-CGStab) method.
            In the call [spbcg ~maxl:maxl ~jac_times_vec:jtv prec],
            - [maxl] is the maximum dimension of the Krylov subspace
                     (defaults to 5),
            - [jtv] computes an approximation to the product between the
                    Jacobian matrix and a vector, and
            - [prec] is a {!preconditioner}.
            
            If the {!jac_times_vec_fn} is omitted, a default implementation
            based on difference quotients is used.

            @cvodes <node7#sss:lin_solv_b> CVSpbcgB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @nocvodes <node7> CVSpilsSetPreconditionerBS
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB
            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @nocvodes <node7> CVSpilsSetJacTimesVecFnBS *)
        val spbcg :
          ?maxl:int
          -> ?jac_times_vec:'d jac_times_vec_fn
          -> ('d, 'k) preconditioner
          -> ('d, 'k) linear_solver

        (** Krylov iterative with the scaled preconditioned transpose-free
            quasi-minimal residual (SPTFQMR) method.
            In the call [sptfqmr ~maxl:maxl ~jac_times_vec:jtv prec],
            - [maxl] is the maximum dimension of the Krylov subspace
                     (defaults to 5),
            - [jtv] computes an approximation to the product between the
                    Jacobian matrix and a vector, and
            - [prec] is a {!preconditioner}.
            
            If the {!jac_times_vec_fn} is omitted, a default implementation
            based on difference quotients is used.

            @cvodes <node7#sss:lin_solv_b> CVSptfqmrB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @nocvodes <node7> CVSpilsSetPreconditionerBS
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB
            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @nocvodes <node7> CVSpilsSetJacTimesVecFnBS *)
        val sptfqmr :
          ?maxl:int
          -> ?jac_times_vec:'d jac_times_vec_fn
          -> ('d, 'k) preconditioner
          -> ('d, 'k) linear_solver

        (** {3:set Solver parameters} *)

        (** Sets the Gram-Schmidt orthogonalization to be used with the
            Spgmr {!linear_solver}.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetGSTypeB *)
        val set_gs_type : ('d, 'k) bsession -> Spils.gramschmidt_type -> unit

        (** Sets the factor by which the Krylov linear solver's convergence test
            constant is reduced from the Newton iteration test constant.
            This factor must be >= 0; passing 0 specifies the default (0.05).

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetEpsLinB *)
        val set_eps_lin : ('d, 'k) bsession -> float -> unit

        (** Resets the maximum Krylov subspace dimension for the Bi-CGStab and
            TFQMR methods. A value <= 0 specifies the default (5.0).

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB *)
        val set_maxl : ('d, 'k) bsession -> int option -> unit

        (** {3:stats Solver statistics} *)

        (** Returns the sizes of the real and integer workspaces used by the spils
            linear solver.

            @cvodes <node5#sss:optout_spils> CVSpilsGetWorkSpace
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
            @return ([real_size], [integer_size]) *)
        val get_work_space       : ('d, 'k) bsession -> int * int

        (** Returns the cumulative number of linear iterations.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumLinIters
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_lin_iters    : ('d, 'k) bsession -> int

        (** Returns the cumulative number of linear convergence failures.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumConvFails
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_conv_fails   : ('d, 'k) bsession -> int

        (** Returns the cumulative number of calls to the setup function with
            [jok=false].

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_prec_evals   : ('d, 'k) bsession -> int

        (** Returns the cumulative number of calls to the preconditioner solve
            function.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecSolves
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_prec_solves  : ('d, 'k) bsession -> int

        (** Returns the cumulative number of calls to the Jacobian-vector
            function.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumJtimesEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_jtimes_evals : ('d, 'k) bsession -> int

        (** Returns the number of calls to the right-hand side callback for
            finite difference Jacobian-vector product approximation. This counter is
            only updated if the default difference quotient function is used.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_rhs_evals    : ('d, 'k) bsession -> int

        (** {3:lowlevel Low-level solver manipulation}

            The {!init} and {!reinit} functions are the preferred way to set or
            change preconditioner functions. These low-level functions are
            provided for experts who want to avoid resetting internal counters
            and other associated side-effects. *)

        (** Change the preconditioner functions without using forward
            sensitivities.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
        val set_preconditioner :
          ('d,'k) bsession
          -> ?setup:'d prec_setup_fn
          -> 'd prec_solve_fn
          -> unit

        (** Change the preconditioner functions using forward sensitivities.

            @nocvodes <node7> CVSpilsSetPreconditionerBS
            @nocvodes <node7#ss:psolve_bs> CVSpilsPrecSolveFnBS
            @nocvodes <node7#ss:psetup_bs> CVSpilsPrecSetupFnBS *)
        val set_preconditioner_with_sens :
          ('d,'k) bsession
          -> ?setup:'d prec_setup_fn_with_sens
          -> 'd prec_solve_fn_with_sens
          -> unit

        (** Change the Jacobian-times-vector function.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @nocvodes <node7> CVSpilsSetJacTimesVecFnBS
            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
            @nocvodes <node7#ss:jtimesv_bs> CVSpilsJacTimesVecFnBS *)
        val set_jac_times_vec_fn :
          ('d,'k) bsession
          -> 'd jac_times_vec_fn
          -> unit

        (** Remove a Jacobian-times-vector function and use the default
            implementation.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB *)
        val clear_jac_times_vec_fn : ('d, 'k) bsession -> unit

        (** Change the preconditioning direction without modifying
            callback functions. If the preconditioning type is changed from
            {{!Spils.preconditioning_type}Spils.PrecNone}
            then {!set_preconditioner} must be called to install the necessary
            callbacks.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPrecTypeB
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val set_prec_type : ('d, 'k) bsession
                            -> Spils.preconditioning_type -> unit
      end (* }}} *)

    (* TODO: Add alternate linear solvers? *)

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

        @cvodes <node7#ss:ODErhs_b> CVRhsFnB
        @cvodes <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis *)
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

        @cvodes <node7#ss:ODErhs_bs> CVRhsFnBS
        @cvodes <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis *)
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

    (** Choice of method for solving non-linear systems that arise in solver
        formulas.

        @cvodes <node7#sss:cvinitb> CVodeCreateB *)
    type ('data, 'kind) iter =
      | Newton of ('data, 'kind) linear_solver
        (** Newton iteration with a given linear solver *)
      | Functional
        (** Functional iteration (non-stiff systems only) *)

    (** Creates and initializes a backward session attached to an existing
        (forward) session. The call
        {[init_backward s lmm iter tol fb tb0 yb0]} has as arguments:
        - [s], the parent (forward) session,
        - [lmm], the linear multistep method (see {!Cvode.lmm}),
        - [iter], either functional or Newton iteration (see {!iter}),
        - [tol], the integration tolerances,
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

        @cvodes <node7#sss:cvinitb> CVodeCreateB
        @cvodes <node7#sss:cvinitb> CVodeInitB
        @cvodes <node7#sss:cvinitb> CVodeInitBS
        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB 
        @raise AdjointNotInitialized The {!init} function has not been called.
        @raise BadFinalTime The final time is outside the interval over which the forward problem was solved. *)
    val init_backward :
         ('d, 'k) Cvode.session
      -> Cvode.lmm
      -> ('d, 'k) iter
      -> ('d, 'k) tolerance
      -> 'd brhsfn
      -> float
      -> ('d, 'k) Nvector.t
      -> ('d, 'k) bsession

    (** Support for backward quadrature equations that may or may
        not depend on forward sensitivities.

        @cvodes <node7#SECTION007210000000000000000> Backward integration of quadrature equations *)
    module Quadrature :
      sig (* {{{ *)
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

            @cvodes <node7#ss:ODErhs_quad_b> CVQuadRhsFnB *)
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

            @cvodes <node7#ss:ODErhs_quad_sens_B> CVQuadRhsFnBS *)
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

            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitB
            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitBS *)
        val init : ('d, 'k) bsession -> 'd bquadrhsfn
                 -> ('d, 'k) Nvector.t -> unit

        (** This function reinitializes the integration of quadrature equations
            during the backward phase.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadReInitB *)
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

            @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
        val set_tolerances : ('d, 'k) bsession -> ('d, 'k) tolerance -> unit

        (** {2:quadout Output functions} *)

        (** Returns the backward quadrature solutions and time reached
            after a successful solver step. The given vectors are filled with
            values calculated during either {!backward_normal} or
            {!backward_one_step} and the value of the independent variable
            is returned.

          @cvodes <node7#sss:quad_get_b> CVodeGetQuadB *)
        val get : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> float

        (** {2:get Querying the solver (optional output functions)}

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration *)

        (** Returns the number of calls to the backward quadrature right-hand
            side function.

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_rhs_evals       : ('d, 'k) bsession -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_err_test_fails  : ('d, 'k) bsession -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_err_weights : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

        (** Returns quadrature-related statistics. These are the
            number of calls to the quadrature function ([nfqevals]) and the
            number of error test failures due to quadrature variables
            ([nqetfails]).

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
            @return ([nfqevals], [nqetfails]) *)
        val get_stats : ('d, 'k) bsession -> int * int
      end (* }}} *)

    (** Integrates a backward ODE system over an interval. The solver takes
        internal steps until it has reached or just passed the specified value.

        @cvodes <node7#sss:cvsolveb> CVodeB (CV_NORMAL)
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

        @cvodes <node7#sss:cvsolveb> CVodeB (CV_ONE_STEP) *)
    val backward_one_step : ('d, 'k) Cvode.session -> float -> unit

    (** Fills the given vector with the solution of the backward ODE problem at
        the returned time, interpolating if necessary.

        @cvodes <node7#sss:cvsolveb> CVodeGetB *)
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

        @cvodes <node5#ss:optional_dky> CVodeGetDky
        @cvodes <node5#ss:optional_dky> CVodeGetAdjIDABmem
        @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
        @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
    val get_dky
          : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> float -> int -> unit

    (** Fills the vector with the interpolated forward solution at the
        given time during a backward simulation.

        @nocvodes <node5#ss:get_adjy> CVodeGetAdjY *)
    val get_y : ('d, 'k) Cvode.session -> ('d, 'k) Nvector.t  -> float -> unit

    (** Reinitializes the backward problem with new parameters and state
        values. The values of the independent variable, i.e., the simulation
        time, and the state variables must be given. It is also possible to
        change the solution method (and linear solver).

        @cvodes <node7#sss:cvinitb> CVodeReInitB
        @raise AdjointNotInitialized The {!init} function has not been called.
        @raise BadFinalTime The final time is not within the forward problem solution interval. *)
    val reinit :
      ('d, 'k) bsession
      -> ?iter_type:('d, 'k) iter
      -> float
      -> ('d, 'k) Nvector.t
      -> unit

    (** {2:set Modifying the solver (optional input functions)} *)

    (** Cancels the storage of sensitivity checkpointing data during forward
        solution (with {!forward_normal} or {!forward_one_step}).

        @cvodes <node7#SECTION00727000000000000000> CVodeAdjSetNoSensi *)
    val set_no_sensitivity : ('d, 'k) Cvode.session -> unit

    (** Sets the integration tolerances for the backward problem.

        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB *)
    val set_tolerances : ('d, 'k) bsession -> ('d, 'k) tolerance -> unit

    (** Specifies the maximum order of the linear multistep method.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxOrdB *)
    val set_max_ord : ('d, 'k) bsession -> int -> unit

    (** Specifies the maximum number of steps taken in attempting to reach
        a given output time.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxNumStepsB *)
    val set_max_num_steps : ('d, 'k) bsession -> int -> unit

    (** Specifies the initial step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetInitStepB *)
    val set_init_step : ('d, 'k) bsession -> float -> unit

    (** Specifies a lower bound on the magnitude of the step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetMinStepB *)
    val set_min_step : ('d, 'k) bsession -> float -> unit

    (** Specifies an upper bound on the magnitude of the step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxStepB *)
    val set_max_step : ('d, 'k) bsession -> float -> unit

    (** Indicates whether the BDF stability limit detection algorithm should be
        used.

        @cvode <node7#ss:optional_input_b> CVodeSetStabLimDet *)
    val set_stab_lim_det : ('d, 'k) bsession -> bool -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the real and integer workspace sizes.

        @cvodes <node5#sss:optout_main> CVodeGetWorkSpace
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @return ([real_size], [integer_size]) *)
    val get_work_space          : ('d, 'k) bsession -> int * int

    (** Returns the cumulative number of internal steps taken by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetNumSteps
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_steps           : ('d, 'k) bsession -> int

    (** Returns the number of calls to the backward right-hand side function.

        @cvodes <node5#sss:optout_main> CVodeGetNumRhsEvals
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_rhs_evals       : ('d, 'k) bsession -> int

    (** Returns the number of calls made to the linear solver's setup function.

        @cvodes <node5#sss:optout_main> CVodeGetNumLinSolvSetups
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_lin_solv_setups : ('d, 'k) bsession -> int

    (** Returns the number of local error test failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumErrTestFails
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_err_test_fails  : ('d, 'k) bsession -> int

    (** Returns the integration method order used during the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastOrder
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_last_order          : ('d, 'k) bsession -> int

    (** Returns the integration method order to be used on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentOrder
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_current_order       : ('d, 'k) bsession -> int

    (** Returns the integration step size taken on the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_last_step           : ('d, 'k) bsession -> float

    (** Returns the integration step size to be attempted on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_current_step        : ('d, 'k) bsession -> float

    (** Returns the the value of the integration step size used on the first
        step.

        @cvodes <node5#sss:optout_main> CVodeGetActualInitStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_actual_init_step    : ('d, 'k) bsession -> float

    (** Returns the the current internal time reached by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentTime
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_current_time        : ('d, 'k) bsession -> float

    (** Returns the number of order reductions dictated by the BDF stability
        limit detection algorithm.

        @cvodes <node5#sss:optout_main> CVodeGetNumStabLimOrderReds
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @cvodes <node3#s:bdf_stab> BDF stability limit detection *)
    val get_num_stab_lim_order_reds : ('d, 'k) bsession -> int

    (** Returns a suggested factor by which the user's tolerances should be
        scaled when too much accuracy has been requested for some internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetTolScaleFactor
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_tol_scale_factor : ('d, 'k) bsession -> float

    (** Returns the solution error weights at the current time.

        @cvodes <node5#sss:optout_main> CVodeGetErrWeights
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @cvodes <node3#ss:ivp_sol> IVP solution (W_i) *)
    val get_err_weights : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

    (** Returns the vector of estimated local errors.

        @cvodes <node5#sss:optout_main> CVodeGetEstLocalErrors
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_est_local_errors : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

    (** Returns the integrator statistics as a group.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_integrator_stats    : ('d, 'k) bsession -> Cvode.integrator_stats

    (** Prints the integrator statistics on the given channel.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val print_integrator_stats  : ('d, 'k) bsession -> out_channel -> unit

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvIters
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_nonlin_solv_iters : ('d, 'k) bsession -> int

    (** Returns the number of nonlinear convergence failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_nonlin_solv_conv_fails : ('d, 'k) bsession -> int

    (** Returns both the numbers of nonlinear iterations performed [nniters] and
        nonlinear convergence failures [nncfails].

        @cvode <node5#sss:optout_main> CVodeGetNonlinSolvStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @return ([nniters], [nncfails]) *)
    val get_nonlin_solv_stats : ('d, 'k) bsession -> int *int

    (** {2:exceptions Exceptions} *)

    (** Adjoint sensitivity analysis was not initialized.

        @cvodes <node7#sss:cvsolvef> CV_NO_ADJ *)
    exception AdjointNotInitialized

    (** Neither {!forward_normal} nor {!forward_one_step} has been called.

        @cvodes <node7#sss:cvsolveb> CV_NO_FWD *)
    exception NoForwardCall

    (** Reinitialization of the forward problem failed at the first checkpoint
        (corresponding to the initial time of the forward problem).

        @cvodes <node7#sss:cvsolveb> CV_REIFWD_FAIL *)
    exception ForwardReinitFailure

    (** An error occured when integrating the forward problem from a
        checkpoint.

        @cvodes <node7#sss:cvsolveb> CV_FWD_FAIL *)
    exception ForwardFailure

    (** No backward problem has been created.

        @cvodes <node7#sss:cvsolveb> CV_NO_BCK *)
    exception NoBackwardProblem

    (** The final time was outside the interval over which the forward
        problem was solved.

        @cvodes <node7#sss:cvinitb> CV_BAD_TB0 *)
    exception BadFinalTime

  end (* }}} *)

