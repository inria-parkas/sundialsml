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
    {% $\dot{y}_Q = f_Q(y, t, p)$%}, to reduce computation costs.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS) *)

open Cvode_impl
open Sundials

(** Alias for Cvode sessions. The Cvodes submodules all work by ‘extending’
    an existing session created with {!Cvode.init}. *)
type ('data, 'kind) session = ('data, 'kind) Cvode.session

(** Integration of pure quadrature equations.
 
    Adds an additional vector $y_Q$ of $N_Q$ quadrature variables defined by
    {% $\frac{\mathrm{d} y_Q}{\mathrm{d}t} = f_Q(t, y, p)$%}. The values of
    these variables are calculated more efficiently since they are excluded
    from the nonlinear solution stage.

    An example session with Cvodes using quadrature variables
    ({openfile cvodes_quad_skel.ml}): {[
#include "examples/ocaml/skeletons/cvodes_quad_skel.ml"
    ]}

    @cvodes <node5#SECTION00570000000000000000> Integration of pure quadrature equations *)
module Quadrature :
  sig
    (** {2:init Initialization and access} *)

    (** Functions defining quadrature variables. They are passed three
        arguments:
        - [t], the value of the independent variable, i.e., the simulation time,
        - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
        - [dyq], a vector for storing the computed value of
                 {% $\dot{y}_Q = f_Q(t, y)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure} exception
        indicates a recoverable error. Any other exception is treated as an
        unrecoverable error.

        {warning [y] and [dyq] should not be accessed after the function
                 returns.}

        @cvodes <node5#ss:user_fct_quad> CVQuadRhsFn *)
    type 'a quadrhsfn = float -> 'a -> 'a -> unit

    (** Activates the integration of quadrature equations. The vector
        gives the initial value of $y_Q$.

        @cvodes <node5#ss:quad_malloc> CVodeQuadInit *)
    val init : ('a, 'k) Cvode.session -> 'a quadrhsfn
                  -> ('a, 'k) Nvector.t -> unit

    (** Reinitialize the integration of quadrature equations. The vector
        gives a new value for $y_Q$.

        @cvodes <node5#ss:quad_malloc> CVodeQuadReInit *)
    val reinit : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t -> unit

    (** Returns the quadrature solutions and time reached after a successful
        solver step. The given vector is filled with values calculated during
        either {!Cvode.solve_normal} or {!Cvode.solve_one_step} and the
        value of the independent variable is returned.

        @cvodes <node5#ss:quad_get> CVodeGetQuad *)
    val get : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t -> float

    (** Returns the interpolated solution or derivatives of quadrature
        variables.

        [get_dky s dky t k] computes the [k]th derivative at time [t], i.e.,
        {% $\frac{d^\mathtt{k}y_Q(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
        and stores it in [dky]. The arguments must satisfy {% $t_n - h_u \leq
        \mathtt{t} \leq t_n$%}—where $t_n$ denotes {!Cvode.get_current_time}
        and $h_u$ denotes {!Cvode.get_last_step},— and
        {% $0 \leq \mathtt{k} \leq q_u$%}—where
        $q_u$ denotes {!Cvode.get_last_order}.

        This function may only be called after a successful return from either
        {!Cvode.solve_normal} or {!Cvode.solve_one_step}.

        @cvodes <node5#ss:quad_get> CVodeGetQuadDky
        @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
        @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
    val get_dky : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t -> float -> int -> unit

    (** {2:tols Tolerances} *)

    (** Tolerances for calculating quadrature variables. *)
    type ('a, 'k) tolerance =
        NoStepSizeControl
        (** Do not use quadrature variables for step-size control (default). *)
      | SStolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('a, 'k) Nvector.t
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

    (** Specify how to use quadrature variables in step size control.

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
    val set_tolerances : ('a, 'k) Cvode.session -> ('a, 'k) tolerance -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the number of calls to the quadrature function.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals *)
    val get_num_rhs_evals       : ('a, 'k) Cvode.session -> int

    (** Returns the number of local error test failures that have occurred
        due to quadrature variables.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails *)
    val get_num_err_test_fails  : ('a, 'k) Cvode.session -> int

    (** Returns the quadrature error weights at the current time.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights *)
    val get_err_weights : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t -> unit

    (** Returns quadrature-related statistics. These are the
        number of calls to the quadrature function ([nfqevals]) and the number
        of error test failures due to quadrature variables ([nqetfails]).

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats
        @return ([nfqevals], [nqetfails]) *)
    val get_stats : ('a, 'k) Cvode.session -> int * int

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
  end

(** (Forward) Sensitivity analysis of ODEs with respect to their parameters.
 
    Formalizes the dependence of a set of ODEs on $N_p$ parameters $p$ and
    calculates the sensitivities $s$ of the solution $y$ to a subset of
    {% $N_s \leq N_p$%} of those parameters;
    {% $s(t) = \frac{\partial y(t)}{\partial p}$%}.

    The {i solution sensitivity} with respect to a single parameter $p_i$
    satisfies the {i (forward) sensitivity equations}
    {% $\dot{s}_i = \frac{\partial f}{\partial y}s_i
                                + \frac{\partial f}{\partial p_i}$%} and
    {% $s_i(t_0) = \frac{\partial y_0(p)}{\partial p_i}$%}, where
    $f$ and $y$ are from the $N$ equations of the original model.

    This documented interface is structured as follows.
    {ol
      {- {{:#init}Initialization} (including {{!Quadrature}Quadrature equations})}
      {- {{:#sensout}Output functions}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#exceptions}Exceptions}}}

    An example session with Cvodes using sensitivity analysis
    ({openfile cvodes_sens_skel.ml}): {[
#include "examples/ocaml/skeletons/cvodes_sens_skel.ml"
    ]}

    @cvodes <node6#ss:forward_usage> Enhanced skeleton for sensitivity analysis *)
module Sensitivity :
  sig
    (** {2:init Initialization} *)

    (** Sensitivity functions that calculate the right-hand sides of all
        sensitivity equations. They are passed the arguments:
        - [t], the value of the independent variable, i.e., the simulation time,
        - [y], the vector of dependent-variable values, i.e., $y(t)$,
        - [dy], the value of the right-hand side of the state
                  equations {% $\dot{y} = f(t, y, p)$%},
        - [s], the array of sensitivity vectors,
        - [ds], an array of vectors to be filled with the values of
                   {% $\dot{s}_i$%} for all $i$, and,
        - [tmp1] and [tmp2], temporary storage vectors.

        Within the function, raising a {!Sundials.RecoverableFailure} exception
        indicates a recoverable error. Any other exception is treated as an
        unrecoverable error.

        {warning Neither [y], [ds], [tmp1], [tmp2], nor the elements of [s]
                 and [dys] should be accessed after the function returns.}

        @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn *)
    type 'a sensrhsfn_all =
      float           (* t *)
      -> 'a           (* y *)
      -> 'a           (* dy *)
      -> 'a array     (* ys *)
      -> 'a array     (* dys *)
      -> 'a           (* tmp1 *)
      -> 'a           (* tmp2 *)
      -> unit

    (** Sensitivity functions that calculate the right-hand side of a single
        sensitivity equation. They are passed the arguments:
        - [t], the value of the independent variable, i.e., the simulation time,
        - [y], the vector of dependent-variable values, i.e., $y(t)$,
        - [dy], the value of the right-hand side of the state
                  equations {% $\dot{y} = f(t, y, p)$%},
        - [i], the index of the sensitivity equation to compute,
        - [si], the {i is}th sensitivity vector,
        - [ds], a vector to be filled with the values of {% $\dot{s}_i$%}, and,
        - [tmp1] and [tmp2], temporary storage vectors.

        Within the function, raising a {!Sundials.RecoverableFailure} exception
        indicates a recoverable error. Any other exception is treated as an
        unrecoverable error.

        {warning [y], [dy], [si], [ds], [tmp1], and [tmp2] should not be
                 accessed after the function returns.}

        @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn *)
    type 'a sensrhsfn1 =
      float           (* t *)
      -> 'a           (* y *)
      -> 'a           (* dy *)
      -> int          (* i *)
      -> 'a           (* ys *)
      -> 'a           (* dys *)
      -> 'a           (* tmp1 *)
      -> 'a           (* tmp2 *)
      -> unit

    (** Specify a sensitivity function. *)
    type 'a sensrhsfn =
        AllAtOnce of 'a sensrhsfn_all option
        (** Calculate sensitivity functions all at once. The argument [None]
            specifies an internal difference quotient routine.

            @cvodes <node6#ss:sensi_malloc> CVodeSensInit
            @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn *)

      | OneByOne of 'a sensrhsfn1 option
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

    (** Specifies problem parameter information for sensitivity calculations.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensParams *)
    type sens_params = {
        pvals  : Sundials.RealArray.t option;
        (** An array of $N_p$ parameters $p$ in $f(t, y, p)$. If specified,
            this array is updated by the solver to pass parameter values to the
            original right-hand side and root functions. *)
        pbar   : Sundials.RealArray.t option;
        (** An array of $N_s$ positive scaling factors. These are only needed
            when estimating tolerances or using the internal difference
            quotient function. *)
        plist  : int array option;
        (** An array of $N_s$ ($< N_p$) indices specifying the parameters
            to analyze. If not specified, the first $N_s$ parameters are
            analyzed. *)
      }

    (** Empty pvals, plist, and pbar. *)
    val no_sens_params : sens_params

    (** Tolerances for calculating sensitivities. *)
    type ('a, 'k) tolerance =
        SStolerances of float * Sundials.RealArray.t
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('a, 'k) Nvector.t array
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
      | EEtolerances
        (** Calculate the integration tolerances for sensitivities
            from those for state variables and the scaling factors
            in {{!sens_params}pbar}. *)

    (** Activates the calculation of forward sensitivities. The call
        [init s tol sm sp fs ys0], has as arguments:
        - [s], a session created with {!Cvode.init},
        - [tol], the tolerances desired,
        - [sm], the solution method,
        - [sp], the parameter information,
        - [fs], the sensitivity function, and,
        - [ys0], initial values of the sensitivities for each parameter.

        @cvodes <node6#ss:sensi_malloc> CVodeSensInit
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit1
        @cvodes <node6#ss:sens_optional_input> CVodeSetSensParams
        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSStolerances
        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSVtolerances
        @cvodes <node6#sss:cvfwdtolerances> CVodeSensEEtolerances *)
    val init : ('a, 'k) Cvode.session
               -> ('a, 'k) tolerance
               -> sens_method
               -> sens_params
               -> 'a sensrhsfn
               -> ('a, 'k) Nvector.t array
               -> unit

    (** Reinitializes the forward sensitivity computation.

        @cvodes <node6#ss:sensi_malloc> CVodeSensReInit *)
    val reinit : ('a, 'k) Cvode.session -> sens_method
                      -> ('a, 'k) Nvector.t array -> unit

    (** Deactivates forward sensitivity calculations without deallocating
        memory. Sensitivities can be reactivated with {!reinit}.

        @cvodes <node6#ss:sensi_malloc> CVodeSensToggleOff *)
    val toggle_off : ('a, 'k) Cvode.session -> unit

    (** Support for quadrature equations that depend not only on
        state variables but also on forward sensitivities.

        Adds an additional vector {% $y_\mathit{QS}$%} of {% $N_\mathit{QS}$%}
        quadrature variables defined by
        {% $\frac{\mathrm{d} y_\mathit{QS}}{\mathrm{d}t}
                = f_{\mathit{QS}}(t, y, s, \dot{y}_Q, p)$%}.
        The values of these variables are calculated more efficiently since they
        are excluded from the nonlinear solution stage. While the sensitivities
        of the ‘pure’ {!Cvodes.Quadrature} variables, $y_Q$, are not provided
        directly, they can be calculated using this more general mechanism.

        @cvodes <node3#SECTION00364000000000000000> Quadratures depending on forward sensitivities
        @cvodes <node6#SECTION00640000000000000000> Integration of quadrature equations depending on forward sensitivities *)
    module Quadrature :
      sig
        (** {2:init Initialization} *)

        (** Functions defining sensitivity-dependent quadrature variables.
            The call [fQS t y ys dyq dyqs tmp1 tmp2] has arguments:
            - [t], the value of the independent variable, i.e., the simulation time,
            - [y], the vector of dependent-variable values, i.e., $y(t)$,
            - [s], the array of sensitivity vectors,
            - [dyq], the value of the quadrature right-hand side, i.e.,
                     {% $\dot{y}_Q$%},
            - [dyqs], an array of vectors for storing the computed values of
              {% $\dot{y}_\mathit{QS} = f_\mathit{QS}(t, y, s, \dot{y}_q)$%}, and,
            - [tmp1] and [tmp2], temporary storage vectors.

            Within the function, raising a {!Sundials.RecoverableFailure}
            exception indicates a recoverable error. Any other exception is
            treated as an unrecoverable error.

            {warning Neither of [y], [dyq], [tmp1], [tmp1], nor the elements
                     of [s] or [dyqs] should be accessed after the function
                     returns.}

           @cvodes <node6#ss:user_fct_quad_sens> CVodeQuadSensRhsFn *)
        type 'a quadsensrhsfn =
           float          (* t *)
           -> 'a          (* y *)
           -> 'a array    (* s *)
           -> 'a          (* dyq *)
           -> 'a array    (* dyqs *)
           -> 'a          (* tmp1 *)
           -> 'a          (* tmp2 *)
           -> unit

        (** Activate the integration of quadrature equations that depend on
            sensitivities. The right-hand sides of the sensitivity-dependent
            quadrature equations are computed with [~fQS] if given, and
            otherwise using an internal implementation based on difference
            quotients. An array of vectors specifies initial values for
            the quadrature equations.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensInit *)
        val init : ('a, 'k) Cvode.session -> ?fQS:'a quadsensrhsfn
                 -> ('a, 'k) Nvector.t array -> unit

        (** Reinitializes the sensitivity-dependent quadrature integration.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensReInit *)
        val reinit : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t array -> unit

        (** {2:tols Tolerance specification} *)

        (** Tolerances for calculating sensitivity-dependent quadrature
            variables. *)
        type ('a, 'k) tolerance =
            NoStepSizeControl
            (** Quadrature variables are not used for step-size control
                (the default). *)
          | SStolerances of float * Sundials.RealArray.t
            (** [(rel, abs)] : scalar relative and absolute tolerances. *)
          | SVtolerances of float * ('a, 'k) Nvector.t array
            (** [(rel, abs)] : scalar relative and vector absolute
                tolerances. *)
          | EEtolerances
            (** Calculate the integration tolerances for the
                sensitivity-dependent quadratures from those provided for
                the pure quadrature variables. *)

        (** Specify how to use quadrature variables in step size control.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeSetQuadSensErrCon
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSStolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSVtolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensEEtolerances *)
        val set_tolerances : ('a, 'k) Cvode.session
                              -> ('a, 'k) tolerance -> unit

        (** {2:quadout Output functions} *)

        (** Returns the quadrature sensitivity solutions and time reached
            after a successful solver step. The given vectors are filled with
            values calculated during either {!Cvode.solve_normal} or
            {!Cvode.solve_one_step} and the value of the independent variable
            is returned.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens *)
        val get : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t array -> float

        (** Returns a single quadrature sensitivity vector after a successful
            solver step. Like {!get}, but the argument [i] specifies a specific
            vector to return.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens1
            @raise BadIS The index is not in the allowed range. *)
        val get1 : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t -> int -> float

        (** Returns the interpolated solution or derivatives of the quadrature
            sensitivity solution.

            [get_dky s dkyqs t k] computes the [k]th derivative at time [t],
            i.e.,
            {% $\frac{d^\mathtt{k}y_\mathit{QS}(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
            and stores it in [dkyqs]. The arguments must satisfy {% $t_n - h_u
            \leq \mathtt{t} \leq t_n$%}—where $t_n$ denotes
            {!Cvode.get_current_time} and $h_u$ denotes
            {!Cvode.get_last_step},—and
            {% $0 \leq \mathtt{k} \leq q_u$%}—where $q_u$ denotes
            {!Cvode.get_last_order}.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky
            @raise BadIS The index is not in the allowed range.
            @raise BadK [k] is not in the range 0, 1, ..., [qlast].
            @raise BadT [t] is not in the allowed range. *)
        val get_dky : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t array
                        -> float -> int -> unit

        (** Returns the interpolated solution or derivatives of a single
            quadrature sensitivity solution vector.
            [get_dky s dkys t k i] is like
            {!get_dky} but restricted to the [i]th sensitivity solution vector.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky1
            @raise BadK [k] is not in the range 0, 1, ..., [qlast].
            @raise BadT [t] is not in the allowed range. *)
        val get_dky1 : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t
                         -> float -> int -> int -> unit

        (** {2:get Querying the solver (optional output functions)} *)

        (** Returns the number of calls to the quadrature right-hand side
            function.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumRhsEvals *)
        val get_num_rhs_evals       : ('a, 'k) Cvode.session -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumErrTestFails *)
        val get_num_err_test_fails  : ('a, 'k) Cvode.session -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensErrWeights *)
        val get_err_weights
              : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t array -> unit

        (** Returns quadrature-related statistics. These are the
            number of calls to the quadrature function ([nfqevals]) and the
            number of error test failures due to quadrature variables
            ([nqetfails]).
         
            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensStats
            @return ([nfqevals], [nqetfails]) *)
        val get_stats : ('a, 'k) Cvode.session -> int * int

        (** {2:exceptions Exceptions} *)

        (** Quadrature integration was not initialized.

            @cvodes <node5#SECTION00642000000000000000> CV_NO_QUAD_SENS *)
        exception QuadSensNotInitialized

        (** The quadrature function failed in an unrecoverable manner.

            @cvodes <node6#SECTION00642000000000000000> CV_QSRHSFUNC_FAIL *)
        exception QuadSensRhsFuncFailure

        (** The quadrature function failed at the first call.

            @cvodes <node6#SECTION00642000000000000000> CV_FIRST_QSRHSFUNC_ERR *)
        exception FirstQuadSensRhsFuncFailure

        (** Convergence test failures occurred too many times due to repeated
            recoverable errors in the quadrature function. Also raised if the
            quadrature function had repeated recoverable errors during the
            estimation of an initial step size (if quadrature variables are
            included in error tests).
          
            @cvodes <node6#SECTION00642000000000000000> CV_REPTD_QSRHSFUNC_ERR *)
        exception RepeatedQuadSensRhsFuncFailure

        (** The quadrature function had a recoverable error, but no
            recovery was possible. This failure mode is rare, as it can occur
            only if the quadrature function fails recoverably after an error
            test failed while at order one.
          
            @cvodes <node6#SECTION00642000000000000000> CV_UNREC_QSRHSFUNC_ERR *)
        exception UnrecoverableQuadSensRhsFuncFailure
      end

    (** {2:sensout Output functions} *)

    (** Returns the sensitivity solution vectors after a successful solver
        step. The given vectors are filled with values calculated during
        either {!Cvode.solve_normal} or {!Cvode.solve_one_step} and the
        value of the independent variable is returned.

        @cvodes <node6#ss:sensi_get> CVodeGetSens *)
    val get : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t array -> float

    (** Returns the interpolated solution or derivatives of the sensitivity
        solution vectors.

        [get_dky s dkys t k] computes the [k]th derivative at time [t], i.e.,
        {% $\frac{d^\mathtt{k}s(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
        and stores it in [dkys]. The arguments must satisfy {% $t_n - h_u \leq
        \mathtt{t} \leq t_n$%}—where $t_n$ denotes {!Cvode.get_current_time}
        and $h_u$ denotes {!Cvode.get_last_step},— and
        {% $0 \leq \mathtt{k} \leq q_u$%}—where
        $q_u$ denotes {!Cvode.get_last_order}.

        This function may only be called after a successful return from either
        {!Cvode.solve_normal} or {!Cvode.solve_one_step}.

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky
        @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
        @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
    val get_dky : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t array
                    -> float -> int -> unit

    (** Returns a single sensitivity solution vector after a successful solver
        step. The given vector is filled with values calculated for the [i]th
        sensitivity vector during either {!Cvode.solve_normal} or
        {!Cvode.solve_one_step} and the value of the independent variable is
        returned.

        @cvodes <node6#ss:sensi_get> CVodeGetSens1
        @raise BadIS The index [i] is not in the allowed range. *)
    val get1 : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t -> int -> float

    (** Returns the interpolated solution or derivatives of a single
        sensitivity solution vector. [get_dky s dkys t k i] is like {!get_dky}
        but restricted to the [i]th sensitivity solution vector.

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky1
        @raise BadIS The index [i] is not in the allowed range.
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky1 : ('a, 'k) Cvode.session -> ('a, 'k) Nvector.t
                     -> float -> int -> int -> unit

    (** {2:set Modifying the solver (optional input functions)} *)

    (** Sets the integration tolerances for sensitivities.

        {b NB}: Unlike the other [set_tolerances] functions in [Cvodes], this
        one does {b not} call {!set_err_con} (which defaults to [false]).

        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSStolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensSVtolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensEEtolerances *)
    val set_tolerances : ('a, 'k) Cvode.session -> ('a, 'k) tolerance -> unit

    (** Sets whether sensitivity variables are used in the error control
        mechanism. The default is [false].

        @cvodes <node5#ss:sens_optional_input> CVodeSetSensErrCon *)
    val set_err_con : ('a, 'k) Cvode.session -> bool -> unit

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
    val set_dq_method : ('a, 'k) Cvode.session -> dq_method -> float -> unit

    (** Sets the maximum number of nonlinear solver iterations for
        sensitivity variables permitted per step.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensMaxNonlinIters *)
    val set_max_nonlin_iters : ('a, 'k) Cvode.session -> int -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the number of calls to the sensitivity function.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumRhsEvals *)
    val get_num_rhs_evals       : ('a, 'k) Cvode.session -> int

    (** Returns the number of calls to the right-hand side function due
        to the internal finite difference approximation of the sensitivity
        equations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetNumRhsEvalsSens *)
    val get_num_rhs_evals_sens  : ('a, 'k) Cvode.session -> int

    (** Returns the number of local error test failures for the sensitivity
        variables that have occurred.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumErrTestFails *)
    val get_num_err_test_fails  : ('a, 'k) Cvode.session -> int

    (** Returns the number of calls made to the linear solver's setup function
        due to forward sensitivity calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumLinSolvSetups *)
    val get_num_lin_solv_setups : ('a, 'k) Cvode.session -> int

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
    val get_stats : ('a, 'k) Cvode.session -> sensitivity_stats

    (** Returns the sensitivity error weights at the current time.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensErrWeights
        @cvodes <node3#e:errwt> IVP solution (W_i, Eq. (2.7)) *)
    val get_err_weights : ('a, 'k) Cvode.session
                            -> ('a, 'k) Nvector.t array -> unit

    (** Returns the number of nonlinear iterations performed for sensitivity
        calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvIters *)
    val get_num_nonlin_solv_iters : ('a, 'k) Cvode.session -> int

    (** Returns the number of nonlinear convergence failures that have occurred
        during sensitivity calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvConvFails *)
    val get_num_nonlin_solv_conv_fails : ('a, 'k) Cvode.session -> int

    (** Returns both the numbers of nonlinear iterations performed [nniters]
        and nonlinear convergence failures [nncfails] during sensitivity
        calculations.

        @cvodes <node6#ss:sens_optional_output> CVodeGetSensNonlinSolvStats
        @return ([nniters], [nncfails]) *)
    val get_nonlin_solv_stats : ('a, 'k) Cvode.session -> int * int

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed for each sensitivity equation separately in the [Staggered1]
        case.

        @cvodes <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvIters *)
    val get_num_stgr_nonlin_solv_iters : ('a, 'k) Cvode.session
                                         -> Sundials.LintArray.t -> unit

    (** Returns the number of nonlinear convergence failures that have occurred
        for each sensitivity equation separately in the [Staggered1] case.

        @cvodes <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvConvFails *)
    val get_num_stgr_nonlin_solv_conv_fails : ('a, 'k) Cvode.session
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
  end

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
  
    An example session with Cvodes using sensitivity analysis
    ({openfile cvodes_adj_skel.ml}): {[
#include "examples/ocaml/skeletons/cvodes_adj_skel.ml"
    ]}

    @cvodes <node7#ss:skeleton_adj> Enhanced Skeleton for Adjoint Sensitivity Analysis *)
module Adjoint :
  sig
    (** A backward session with the CVODES solver. Multiple backward sessions
        may be associated with a single parent session.

        @cvodes <node7#sss:cvinitb> Backward problem initialization functions *)
    type ('data, 'kind) bsession = ('data, 'kind) AdjointTypes.bsession

    (** Alias for backward sessions based on serial nvectors. *)
    type serial_bsession = (Nvector_serial.data, Nvector_serial.kind) bsession

    (** {2:fwd Forward solution} *)

    (** Specifies the type of interpolation to use between checkpoints.

        @cvodes <node3#ss:checkpointing> Checkpointing scheme *)
    type interpolation = IPolynomial (** {cconst CV_POLYNOMIAL} *)
                       | IHermite    (** {cconst CV_HERMITE} *)

    (** Activates the forward-backward problem. The arguments specify the number
        of integration steps between consecutive checkpoints, and the type of
        variable-degree interpolation.

        @cvodes <node7#sss:cvadjinit> CVodeAdjInit *)
    val init : ('a, 'k) Cvode.session -> int -> interpolation -> unit

    (** Integrates the forward problem over an interval and saves
        checkpointing data. The arguments are the next time at which a solution
        is desired ([tout]) and a vector for storing the computed result
        ([yret]). The function returns a triple [tret, ncheck, sr]: the time
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
      ('a, 'k) Cvode.session
      -> float
      -> ('a, 'k) Nvector.t
      -> float * int * Cvode.solver_result

    (** Integrates the forward problem over an interval and saves
        checkpointing data. The arguments are the next time at which a solution
        is desired ([tout]) and a vector for storing the computed result
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
      ('a, 'k) Cvode.session
      -> float
      -> ('a, 'k) Nvector.t
      -> float * int * Cvode.solver_result

    (** {2:linear Linear solvers} *)

    (** Linear solvers used in backward problems.

        @cvodes <node7#sss:lin_solv_b> Linear Solver Initialization Functions *)
    type ('data, 'kind) linear_solver =
            ('data, 'kind) AdjointTypes.linear_solver

    (** Alias for linear solvers that are restricted to serial nvectors. *)
    type serial_linear_solver =
            (Nvector_serial.data, Nvector_serial.kind) linear_solver

    (** Workspaces with three temporary vectors. *)
    type 'a triple = 'a * 'a * 'a

    (** Arguments common to Jacobian callback functions.

        @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB
        @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB 
        @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
        @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
        @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
    type ('t, 'a) jacobian_arg =
      {
        jac_t   : float;        (** The independent variable. *)
        jac_y   : 'a;           (** The forward solution vector. *)
        jac_yb  : 'a;           (** The backward dependent variable vector. *)
        jac_fyb : 'a;           (** The backward right-hand side function [fB]. *)
        jac_tmp : 't            (** Workspace data. *)
      }

    (** The range of nonzero entries in a band matrix. *)
    type bandrange = Cvode_impl.bandrange =
      { mupper : int; (** The upper half-bandwidth.  *)
        mlower : int; (** The lower half-bandwidth.  *) }

    (** Diagonal approximation of Jacobians by difference quotients. *)
    module Diag :
      sig
        (** A linear solver based on Jacobian approximation by difference
            quotients.

            @cvodes <node7#sss:lin_solv_b> CVDiagB *)
        val solver : ('data, 'kind) linear_solver

        (** Returns the sizes of the real and integer workspaces used by the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetWorkSpace
            @return ([real_size], [integer_size]) *)
        val get_work_space : ('a, 'k) bsession -> int * int

        (** Returns the number of calls made to the right-hand side
            function due to finite difference Jacobian approximation in the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetNumRhsEvals *)
        val get_num_rhs_evals : ('a, 'k) bsession -> int
      end


    (** Direct Linear Solvers operating on dense and banded matrices. *)
    module Dls :
      sig

        (** Callback functions that compute dense approximations to a Jacobian
            matrix. In the call [dense_jac_fn arg jac], [arg] is a
            {!jacobian_arg} with three work vectors and the computed Jacobian
            must be stored in [jac].

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
        type dense_jac_fn = (RealArray.t triple, RealArray.t) jacobian_arg
                                -> Dls.DenseMatrix.t -> unit

        (** A direct linear solver on dense matrices. The optional argument
            specifies a callback function for computing an approximation to the
            Jacobian matrix. If this argument is omitted, then a default
            implementation based on difference quotients is used.

            @cvodes <node7#sss:lin_solv_b> CVDenseB
            @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
        val dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

        (** A direct linear solver on dense matrices using LAPACK. See {!dense}.
            Only available if {!Sundials.lapack_enabled}.

            @cvodes <node7#sss:lin_solv_b> CVLapackDenseB
            @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
        val lapack_dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

        (** Callback functions that compute banded approximations to
            a Jacobian matrix. In the call
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
        type band_jac_fn = bandrange
                            -> (RealArray.t triple, RealArray.t) jacobian_arg
                            -> Dls.BandMatrix.t -> unit

        (** A direct linear solver on banded matrices. The optional argument
            specifies a callback function for computing an approximation to the
            Jacobian matrix. If this argument is omitted, then a default
            implementation based on difference quotients is used. The other
            argument gives the width of the bandrange.

            @cvodes <node7#sss:lin_solv_b> CVBandB
            @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        val band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

        (** A direct linear solver on banded matrices using LAPACK. See {!band}.
            Only available if {!Sundials.lapack_enabled}.

            @cvodes <node7#sss:lin_solv_b> CVLapackBandB
            @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        val lapack_band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

        (** {3:stats Solver statistics} *)

        (** Returns the sizes of the real and integer workspaces used by a direct
            linear solver.

            @cvode <node5#sss:optout_dls> CVDlsGetWorkSpace
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
            @return ([real_size], [integer_size]) *)
        val get_work_space : serial_bsession -> int * int

        (** Returns the number of calls made by a direct linear solver to the
            Jacobian approximation function.

            @cvode <node5#sss:optout_dls> CVDlsGetNumJacEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_jac_evals : serial_bsession -> int

        (** Returns the number of calls to the right-hand side callback due to
            the finite difference Jacobian approximation.

            @cvode <node5#sss:optout_dls> CVDlsGetNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_rhs_evals : serial_bsession -> int

        (** {3:lowlevel Low-level solver manipulation}

            The {!init} and {!reinit} functions are the preferred way to set or
            change a Jacobian function. These low-level functions are provided for
            experts who want to avoid resetting internal counters and other
            associated side-effects. *)

        (** Change the dense Jacobian function.
       
            @cvode <node5#SECTION00728200000000000000> CVDlsSetDenseJacFnB *)
        val set_dense_jac_fn : serial_bsession -> dense_jac_fn -> unit

        (** Remove a dense Jacobian function and use the default
            implementation.

            @cvode <node5#SECTION00728200000000000000> CVDlsSetDenseJacFnB *)
        val clear_dense_jac_fn : serial_bsession -> unit

        (** Change the band Jacobian function.

            @cvode <node5#SECTION00728300000000000000> CVDlsSetBandJacFnB *)
        val set_band_jac_fn : serial_bsession -> band_jac_fn -> unit

        (** Remove a banded Jacobian function and use the default
            implementation.

            @cvode <node5#SECTION00728300000000000000> CVDlsSetBandJacFnB *)
        val clear_band_jac_fn : serial_bsession -> unit
      end

    (** Scaled Preconditioned Iterative Linear Solvers.

        @cvodes <node7#ss:optional_output_b> Optional output functions for the backward problem.
        @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
        @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
    module Spils :
      sig
        (** Arguments passed to the preconditioner solver function.

            @cvode <node7#ss:psolve_b> CVSpilsPrecSolveFnB *)
        type 'a prec_solve_arg =
          {
            rhs   : 'a;         (** Right-hand side vector of the linear system. *)
            gamma : float;      (** Scalar $\gamma$ in the Newton
                                    matrix given by $M = I - \gamma J$. *)
            delta : float;      (** Input tolerance for iterative methods. *)
            left  : bool;       (** [true] for left preconditioning and
                                    [false] for right preconditioning. *)
          }

        (** Callback functions that solve a linear system involving a
            preconditioner matrix. In the call [prec_solve_fn jac arg z],
            [jac] is a {!jacobian_arg} with one work vector, [arg] is
            a {!prec_solve_arg} that specifies the linear system, and [z] is
            computed to solve {% $P\mathtt{z} = \mathtt{arg.rhs}$%}.
            $P$ is a preconditioner matrix, which approximates, however crudely,
            the Newton matrix {% $M = I - \gamma J$%} where
            {% $J = \frac{\partial f}{\partial y}$%}.

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning The elements of [jac], [arg], and [z] should not be
            accessed after the function has returned.}

            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB *)
        type 'a prec_solve_fn =
          ('a, 'a) jacobian_arg
          -> 'a prec_solve_arg
          -> 'a
          -> unit

        (** Callback functions that preprocess or evaluate Jacobian-related data
            need by {!prec_solve_fn}. In the call [prec_setup_fn jac jok gamma],
            [jac] is a {!jacobian_arg} with one work vector, [jok] indicates
            whether any saved Jacobian-related data can be reused with the
            current value of [gamma], and [gamma] is the scalar $\gamma$ in the
            Newton matrix {% $M = I - \gamma J$%} where $J$ is the Jacobian
            matrix. A function should return [true] if Jacobian-related data was
            updated and [false] if saved data was reused.

            Raising {!Sundials.RecoverableFailure} indicates a recoverable
            error. Any other exception is treated as an unrecoverable error.

            {warning The elements of [jac] should not be accessed after the
                     function has returned.}

            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
        type 'a prec_setup_fn =
          ('a triple, 'a) jacobian_arg
          -> bool
          -> float
          -> bool

        (** Callback functions that compute the Jacobian times a vector. In the
            call [jac_times_vec_fn arg v jv], [arg] is a {!jacobian_arg} with one
            work vector, [v] is the vector multiplying the Jacobian, and [jv] is
            the vector in which to store the
            result—{% $\mathtt{jv} = J\mathtt{v}$%}.
          
            Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
            Any other exception is treated as an unrecoverable error.

            {warning Neither the elements of [arg] nor [v] or [jv] should be
                     accessed after the function has returned.}

            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB *)
        type 'a jac_times_vec_fn =
          ('a, 'a) jacobian_arg
          -> 'a
          -> 'a
          -> unit

        (** Specifies a preconditioner, including the type of preconditioning
            (none, left, right, or both) and callback functions. The following
            functions and those in {!Banded} and {!Cvode_bbd} construct
            preconditioners.

            The {!prec_solve_fn} is usually mandatory. The {!prec_setup_fn} can
            be omitted if not needed. If the {!jac_times_vec_fn} is omitted, a
            default implementation based on difference quotients is used.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB
            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB *)
        type ('a, 'k) preconditioner =
          ('a, 'k) AdjointTypes.SpilsTypes.preconditioner

        (** No preconditioning.  *)
        val prec_none : ('a, 'k) preconditioner

        (** Left preconditioning. {% $(P^{-1}A)x = P^{-1}b$ %}. *)
        val prec_left :
          ?setup:'a prec_setup_fn
          -> ?jac_times_vec:'a jac_times_vec_fn
          -> 'a prec_solve_fn
          -> ('a, 'k) preconditioner

        (** Right preconditioning. {% $(AP^{-1})Px = b$ %}. *)
        val prec_right :
          ?setup:'a prec_setup_fn
          -> ?jac_times_vec:'a jac_times_vec_fn
          -> 'a prec_solve_fn
          -> ('a, 'k) preconditioner

        (** Left and right preconditioning.
            {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)
        val prec_both :
          ?setup:'a prec_setup_fn
          -> ?jac_times_vec:'a jac_times_vec_fn
          -> 'a prec_solve_fn
          -> ('a, 'k) preconditioner

        (** Banded preconditioners.  *)
        module Banded : sig

          (** A band matrix {!preconditioner} based on difference quotients.
              The call [prec_left br] instantiates a left preconditioner which
              generates a banded approximation to the Jacobian with [br.mlower]
              sub-diagonals and [br.mupper] super-diagonals.

              @cvode <node7#SECTION00741000000000000000> CVBandPrecInitB *)
          val prec_left :
            ?jac_times_vec:(RealArray.t jac_times_vec_fn)
            -> bandrange
            -> (Nvector_serial.data, Nvector_serial.kind) preconditioner

          (** Like {!prec_left} but preconditions from the right.

              @cvode <node7#SECTION00741000000000000000> CVBandPrecInitB *)
          val prec_right :
            ?jac_times_vec:(RealArray.t jac_times_vec_fn)
            -> bandrange
            -> (Nvector_serial.data, Nvector_serial.kind) preconditioner

          (** Like {!prec_left} but preconditions from both sides.

              @cvode <node7#SECTION00741000000000000000> CVBandPrecInitB *)
          val prec_both :
            ?jac_times_vec:(RealArray.t jac_times_vec_fn)
            -> bandrange
            -> (Nvector_serial.data, Nvector_serial.kind) preconditioner

          (** {4:stats Banded statistics} *)

          (** Returns the sizes of the real and integer workspaces
              used by the banded preconditioner module.

              @cvodes <node5#sss:cvbandpre> CVBandPrecGetWorkSpace
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
              @return ([real_size], [integer_size]) *)
          val get_work_space : serial_bsession -> int * int

          (** Returns the number of calls to the right-hand side callback for the
              difference banded Jacobian approximation. This counter is only updated
              if the default difference quotient function is used.

              @cvodes <node5#sss:cvbandpre> CVBandPrecGetNumRhsEvals
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val get_num_rhs_evals : serial_bsession -> int
        end

        (** {3:lsolvers Solvers} *)

        (** Krylov iterative solver using the scaled preconditioned generalized
            minimum residual (GMRES) method. In the call
            [spgmr ~maxl:maxl prec], [maxl] is the maximum dimension of the
            Krylov subspace (defaults to 5), and [prec] is a {!preconditioner}.

            @cvodes <node7#sss:lin_solv_b> CVSpgmrB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB *)
        val spgmr : ?maxl:int -> ('a, 'k) preconditioner
                      -> ('a, 'k) linear_solver

        (** Krylov iterative solver using the scaled preconditioned biconjugate
            stabilized (Bi-CGStab) method. In the call [spbcg ~maxl:maxl prec],
            [maxl] is the maximum dimension of the Krylov subspace (defaults to
            5), and [prec] is a {!preconditioner}.

            @cvodes <node7#sss:lin_solv_b> CVSpbcgB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB *)
        val spbcg : ?maxl:int -> ('a, 'k) preconditioner
                      -> ('a, 'k) linear_solver

        (** Krylov iterative with the scaled preconditioned transpose-free
            quasi-minimal residual (SPTFQMR) method.
            In the call [sptfqmr ~maxl:maxl prec], [maxl] is the maximum dimension
            of the Krylov subspace (defaults to 5), and [prec] is a
            {!preconditioner}.

            @cvodes <node7#sss:lin_solv_b> CVSptfqmrB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB *)
        val sptfqmr : ?maxl:int -> ('a, 'k) preconditioner
                      -> ('a, 'k) linear_solver

        (** {3:set Solver parameters} *)

        (** The type of Gram-Schmidt orthogonalization.

            @cvode <node9#ss:spgmr> ModifiedGS/ClassicalGS *)
        type gramschmidt_type = Spils.gramschmidt_type =
          | ModifiedGS   (** Modified Gram-Schmidt orthogonalization
                             {cconst MODIFIED_GS} *)
          | ClassicalGS  (** Classical Gram Schmidt orthogonalization
                             {cconst CLASSICAL_GS} *)

        (** Sets the Gram-Schmidt orthogonalization to be used with the
            Spgmr {!linear_solver}.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetGSTypeB *)
        val set_gs_type : ('a, 'k) bsession -> gramschmidt_type -> unit

        (** Sets the factor by which the Krylov linear solver's convergence test
            constant is reduced from the Newton iteration test constant.
            This factor must be >= 0; passing 0 specifies the default (0.05).

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetEpsLinB *)
        val set_eps_lin : ('a, 'k) bsession -> float -> unit

        (** Resets the maximum Krylov subspace dimension for the Bi-CGStab and
            TFQMR methods. A value <= 0 specifies the default (5.0).

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB *)
        val set_maxl : ('a, 'k) bsession -> int option -> unit

        (** {3:stats Solver statistics} *)

        (** Returns the sizes of the real and integer workspaces used by the spils
            linear solver.

            @cvodes <node5#sss:optout_spils> CVSpilsGetWorkSpace
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
            @return ([real_size], [integer_size]) *)
        val get_work_space       : ('a, 'k) bsession -> int * int

        (** Returns the cumulative number of linear iterations.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumLinIters
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_lin_iters    : ('a, 'k) bsession -> int

        (** Returns the cumulative number of linear convergence failures.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumConvFails
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_conv_fails   : ('a, 'k) bsession -> int

        (** Returns the cumulative number of calls to the setup function with
            [jok=false].

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_prec_evals   : ('a, 'k) bsession -> int

        (** Returns the cumulative number of calls to the preconditioner solve
            function.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecSolves
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_prec_solves  : ('a, 'k) bsession -> int

        (** Returns the cumulative number of calls to the Jacobian-vector
            function.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumJtimesEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_jtimes_evals : ('a, 'k) bsession -> int

        (** Returns the number of calls to the right-hand side callback for
            finite difference Jacobian-vector product approximation. This counter is
            only updated if the default difference quotient function is used.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_rhs_evals    : ('a, 'k) bsession -> int

        (** {3:lowlevel Low-level solver manipulation}

            The {!init} and {!reinit} functions are the preferred way to set or
            change preconditioner functions. These low-level functions are
            provided for experts who want to avoid resetting internal counters
            and other associated side-effects. *)

        (** Change the preconditioner functions.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
        val set_preconditioner :
          ('a,'k) bsession
          -> ?setup:'a prec_setup_fn
          -> 'a prec_solve_fn
          -> unit

        (** Change the Jacobian-times-vector function.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB *)
        val set_jac_times_vec_fn :
          ('a,'k) bsession
          -> 'a jac_times_vec_fn
          -> unit

        (** Remove a Jacobian-times-vector function and use the default
            implementation.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB *)
        val clear_jac_times_vec_fn : ('a, 'k) bsession -> unit

        (** Change the preconditioning direction without modifying
            callback functions. If the preconditioning type is changed from
            {{!Spils.preconditioning_type}Spils.PrecNone}
            then {!set_preconditioner} must be called to install the necessary
            callbacks.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPrecTypeB
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val set_prec_type : ('a, 'k) bsession
                            -> Spils.preconditioning_type -> unit
      end

    (* TODO: Add alternate linear solvers? *)

    (** {2:bsolve Backward solutions} *)

    (** Backward functions without forward sensitivities. They are passed
        the arguments:
        - [t], the value of the independent variable, i.e., the simulation time,
        - [y], the vector of dependent-variable values, i.e., $y(t)$,
        - [yb], the vector of backward dependent-variable values, i.e., $y_B(t)$,
        - [dyb], a vector for storing the values
                 {% $\dot{y}_B = f_B(t, y, y_B)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure} exception
        indicates a recoverable error. Any other exception is treated as an
        unrecoverable error.

        {warning [y], [yb], and [dyb] should not be accessed after the function
                 returns.}

        @cvodes <node7#ss:ODErhs_b> CVRhsFnB
        @cvodes <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis *)
    type 'a brhsfn_no_sens =
      float    (* t *)
      -> 'a    (* y *)
      -> 'a    (* yb *)
      -> 'a    (* dyb *)
      -> unit

    (** Backward functions with forward sensitivities. They are passed the
        arguments:
        - [t], the value of the independent variable, i.e., the simulation time,
        - [y], the vector of dependent-variable values, i.e., $y(t)$,
        - [s], the array of forward sensitivity vectors,
        - [yb], the vector of backward dependent-variable values,
                i.e., $y_B(t)$,
        - [dyb], a vector for storing the values
                 {% $\dot{y}_B = f_B(t, y, y_S, y_B)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure} exception
        indicates a recoverable error. Any other exception is treated as an
        unrecoverable error.

        {warning Neither [y], [yb], [dyb], nor the elements of [s] should be
                 accessed after the function returns.}

        @cvodes <node7#ss:ODErhs_bs> CVRhsFnBS
        @cvodes <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis *)
    type 'a brhsfn_with_sens =
      float        (* t *)
      -> 'a        (* y *)
      -> 'a array  (* ys *)
      -> 'a        (* yb *)
      -> 'a        (* dyb *)
      -> unit

    (** Functions that evaluate the right-hand side of a backward ODE system
        with or without forward sensitivities. *)
    type 'a brhsfn =
        NoSens of 'a brhsfn_no_sens
          (** No dependency on forward sensitivities. *)
      | WithSens of 'a brhsfn_with_sens
          (** Dependency on forward sensitivities. *)

    (** Tolerance specifications. *)
    type ('a, 'k) tolerance =
      | SStolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('a, 'k) Nvector.t
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
         ('a, 'k) Cvode.session
      -> Cvode.lmm
      -> ('a, 'k) iter
      -> ('a, 'k) tolerance
      -> 'a brhsfn
      -> float
      -> ('a, 'k) Nvector.t
      -> ('a, 'k) bsession

    (** Support for backward quadrature equations that may or may
        not depend on forward sensitivities.

        @cvodes <node7#SECTION007210000000000000000> Backward integration of quadrature equations *)
    module Quadrature :
      sig
        (** {2:init Initialization} *)

        (** Functions defining backward quadrature variables without forward
            sensitivities.
            The call [fBQS t y yb dqb] has arguments:
            - [t], the value of the independent variable, i.e.,
                   the simulation time,
            - [y], the vector of dependent-variable values, i.e., $y(t)$,
            - [yb], the vector of backward dependent-variable values,
                    i.e., $y_B(t)$,
            - [dqb], a vector for storing the computed value of
                     {% $\dot{y}_\mathit{QB} = f_\mathit{QB}(t, y, y_B)$%}.

            Within the function, raising a {!Sundials.RecoverableFailure}
            exception indicates a recoverable error. Any other exception is
            treated as an unrecoverable error.

            {warning [y], [yb], and [dqb] should not be accessed after the
                     function returns.}

            @cvodes <node7#ss:ODErhs_quad_b> CVQuadRhsFnB *)
        type 'a bquadrhsfn_no_sens =
          float       (* t *)
          -> 'a       (* y *)
          -> 'a       (* yb *)
          -> 'a       (* dqb *)
          -> unit

        (** Functions defining backward quadrature variables with forward
            sensitivities.
            The call [fBQS t y ys yb dqb] has arguments:
            - [t], the value of the independent variable, i.e.,
                   the simulation time,
            - [y], the vector of dependent-variable values, i.e., $y(t)$,
            - [s], the array of forward sensitivity vectors,
            - [yb], the vector of backward dependent-variable values,
                    i.e., $y_B(t)$,
            - [dqb], a vector for storing the computed value of
                   {% $\dot{y}_\mathit{QB} = f_\mathit{QB}(t, y, y_S, y_B)$%}.

            Within the function, raising a {!Sundials.RecoverableFailure}
            exception indicates a recoverable error. Any other exception is
            treated as an unrecoverable error.

            {warning Neither [y], [yb], [dqb], nor the elements of [s] should
                     be accessed after the function returns.}

            @cvodes <node7#ss:ODErhs_quad_sens_B> CVQuadRhsFnBS *)
        type 'a bquadrhsfn_with_sens =
          float        (* t *)
          -> 'a        (* y *)
          -> 'a array  (* ys *)
          -> 'a        (* yb *)
          -> 'a        (* qbdot *)
          -> unit

        (** These functions compute the quadrature equation right-hand side for
            the backward problem. *)
        type 'a bquadrhsfn =
            NoSens of 'a bquadrhsfn_no_sens
            (** Doesn't depend on forward sensitivities. *)
          | WithSens of 'a bquadrhsfn_with_sens
            (** Depends on forward sensitivities. *)

        (** This function, [init s fQB yQB0], activates integration of
            quadrature equations, with or without sensitivities, where [fQB]
            computes the right-hand side of the backward quadrature equations,
            and [yQB0] contains the values of the quadrature variables at [tB0].

            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitB
            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitBS *)
        val init : ('a, 'k) bsession -> 'a bquadrhsfn
                 -> ('a, 'k) Nvector.t -> unit

        (** This function reinitializes the integration of quadrature equations
            during the backward phase.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadReInitB *)
        val reinit : ('a, 'k) bsession -> ('a, 'k) Nvector.t -> unit

        (** {2:tols Tolerance specification} *)

        (** Tolerances for calculating backward quadrature variables. *)
        type ('a, 'k) tolerance =
            NoStepSizeControl
            (** Quadrature variables are not used for step-size control
                (the default). *)
          | SStolerances of float * float
            (** [(rel, abs)] : scalar relative and absolute tolerances. *)
          | SVtolerances of float * ('a, 'k) Nvector.t
            (** [(rel, abs)] : scalar relative and vector absolute
                tolerances. *)

        (** Specify how to use quadrature variables in step size control.

            @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
        val set_tolerances : ('a, 'k) bsession -> ('a, 'k) tolerance -> unit

        (** {2:quadout Output functions} *)

        (** Returns the backward quadrature solutions and time reached
            after a successful solver step. The given vectors are filled with
            values calculated during either {!backward_normal} or
            {!backward_one_step} and the value of the independent variable
            is returned.

          @cvodes <node7#sss:quad_get_b> CVodeGetQuadB *)
        val get : ('a, 'k) bsession -> ('a, 'k) Nvector.t -> float

        (** {2:get Querying the solver (optional output functions)}

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration *)

        (** Returns the number of calls to the backward quadrature right-hand
            side function.

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_rhs_evals       : ('a, 'k) bsession -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_err_test_fails  : ('a, 'k) bsession -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_err_weights : ('a, 'k) bsession -> ('a, 'k) Nvector.t -> unit

        (** Returns quadrature-related statistics. These are the
            number of calls to the quadrature function ([nfqevals]) and the
            number of error test failures due to quadrature variables
            ([nqetfails]).

            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_stats : ('a, 'k) bsession -> int * int
      end

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
    val backward_normal : ('a, 'k) Cvode.session -> float -> unit

    (** Like {!backward_normal} but returns after one internal solver step.

        @cvodes <node7#sss:cvsolveb> CVodeB (CV_ONE_STEP) *)
    val backward_one_step : ('a, 'k) Cvode.session -> float -> unit

    (** Fills the given vector with the solution of the backward ODE problem at
        the returned time, interpolating if necessary.

        @cvodes <node7#sss:cvsolveb> CVodeGetB *)
    val get : ('a, 'k) bsession -> ('a, 'k) Nvector.t -> float

    (** Returns the interpolated solution or derivatives.
        [get_dky s dky t k] computes the [k]th derivative of the backward
        function at time [t], i.e.,
        {% $\frac{d^\mathtt{k}y_B(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
        and stores it in [dky]. The arguments must satisfy
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
          : ('a, 'k) bsession -> ('a, 'k) Nvector.t -> float -> int -> unit

    (** Reinitializes the backward problem with new parameters and state
        values. The values of the independent variable, i.e., the simulation
        time, and the state variables must be given.

        @cvodes <node7#sss:cvinitb> CVodeReInitB
        @raise AdjointNotInitialized The {!init} function has not been called.
        @raise BadFinalTime The final time is not within the forward problem solution interval. *)
    val reinit :
      ('a, 'k) bsession
      -> ?iter_type:('a, 'k) iter
      -> float
      -> ('a, 'k) Nvector.t
      -> unit

    (** {2:set Modifying the solver (optional input functions)} *)

    (** Cancels the storage of sensitivity checkpointing data during forward
        solution (with {!forward_normal} or {!forward_one_step}).

        @cvodes <node7#SECTION00727000000000000000> CVodeAdjSetNoSensi *)
    val set_no_sensitivity : ('a, 'k) Cvode.session -> unit

    (** Sets the integration tolerances for the backward problem.

        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB *)
    val set_tolerances : ('a, 'k) bsession -> ('a, 'k) tolerance -> unit

    (** Specifies the maximum order of the linear multistep method.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxOrdB *)
    val set_max_ord : ('a, 'k) bsession -> int -> unit

    (** Specifies the maximum number of steps taken in attempting to reach
        a given output time.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxNumStepsB *)
    val set_max_num_steps : ('a, 'k) bsession -> int -> unit

    (** Specifies the initial step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetInitStepB *)
    val set_init_step : ('a, 'k) bsession -> float -> unit

    (** Specifies a lower bound on the magnitude of the step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetMinStepB *)
    val set_min_step : ('a, 'k) bsession -> float -> unit

    (** Specifies an upper bound on the magnitude of the step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxStepB *)
    val set_max_step : ('a, 'k) bsession -> float -> unit

    (** Indicates whether the BDF stability limit detection algorithm should be
        used.

        @cvode <node7#ss:optional_input_b> CVodeSetStabLimDet *)
    val set_stab_lim_det : ('a, 'k) bsession -> bool -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the real and integer workspace sizes.

        @cvodes <node5#sss:optout_main> CVodeGetWorkSpace
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @return ([real_size], [integer_size]) *)
    val get_work_space          : ('a, 'k) bsession -> int * int

    (** Returns the cumulative number of internal steps taken by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetNumSteps
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_steps           : ('a, 'k) bsession -> int

    (** Returns the number of calls to the backward right-hand side function.

        @cvodes <node5#sss:optout_main> CVodeGetNumRhsEvals
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_rhs_evals       : ('a, 'k) bsession -> int

    (** Returns the number of calls made to the linear solver's setup function.

        @cvodes <node5#sss:optout_main> CVodeGetNumLinSolvSetups
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_lin_solv_setups : ('a, 'k) bsession -> int

    (** Returns the number of local error test failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumErrTestFails
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_err_test_fails  : ('a, 'k) bsession -> int

    (** Returns the integration method order used during the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastOrder
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_last_order          : ('a, 'k) bsession -> int

    (** Returns the integration method order to be used on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentOrder
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_current_order       : ('a, 'k) bsession -> int

    (** Returns the integration step size taken on the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_last_step           : ('a, 'k) bsession -> float

    (** Returns the integration step size to be attempted on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_current_step        : ('a, 'k) bsession -> float

    (** Returns the the value of the integration step size used on the first
        step.

        @cvodes <node5#sss:optout_main> CVodeGetActualInitStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_actual_init_step    : ('a, 'k) bsession -> float

    (** Returns the the current internal time reached by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentTime
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_current_time        : ('a, 'k) bsession -> float

    (** Returns the number of order reductions dictated by the BDF stability
        limit detection algorithm.

        @cvodes <node5#sss:optout_main> CVodeGetNumStabLimOrderReds
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @cvodes <node3#s:bdf_stab> BDF stability limit detection *)
    val get_num_stab_lim_order_reds : ('a, 'k) bsession -> int

    (** Returns a suggested factor by which the user's tolerances should be
        scaled when too much accuracy has been requested for some internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetTolScaleFactor
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_tol_scale_factor : ('a, 'k) bsession -> float

    (** Returns the solution error weights at the current time.

        @cvodes <node5#sss:optout_main> CVodeGetErrWeights
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @cvodes <node3#ss:ivp_sol> IVP solution (W_i) *)
    val get_err_weights : ('a, 'k) bsession -> ('a, 'k) Nvector.t -> unit

    (** Returns the vector of estimated local errors.

        @cvodes <node5#sss:optout_main> CVodeGetEstLocalErrors
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_est_local_errors : ('a, 'k) bsession -> ('a, 'k) Nvector.t -> unit

    (** Returns the integrator statistics as a group.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_integrator_stats    : ('a, 'k) bsession -> Cvode.integrator_stats

    (** Prints the integrator statistics on the given channel.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val print_integrator_stats  : ('a, 'k) bsession -> out_channel -> unit

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvIters
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_nonlin_solv_iters : ('a, 'k) bsession -> int

    (** Returns the number of nonlinear convergence failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_num_nonlin_solv_conv_fails : ('a, 'k) bsession -> int

    (** Returns both the numbers of nonlinear iterations performed [nniters] and
        nonlinear convergence failures [nncfails].

        @cvode <node5#sss:optout_main> CVodeGetNonlinSolvStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
    val get_nonlin_solv_stats : ('a, 'k) bsession -> int *int

    (** {2:exceptions Exceptions} *)

    (** Adjoint sensitivity analysis was not initialized.

        @cvodes <node7#sss:cvsolvef> CV_NO_ADJ *)
    exception AdjointNotInitialized

    (** Neither {!forward_normal} nor {!forward_one_step} has previously been
        called.

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

  end

