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
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for CVODES v2.7.0                  *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Sensitivity analysis (forward and adjoint) and quadrature equations
    for CVODE.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

open Cvode_impl
open Sundials

type ('data, 'kind) session = ('data, 'kind) Cvode.session
type ('data, 'kind) nvector = ('data, 'kind) Nvector.t

(** {2:quad Quadrature Equations} *)

module Quadrature :
  sig
    (** A skeleton of an enhanced main program:
        + {b Initialize a session [s] per the skeleton at
           {!Cvode.session}}
        {[...]}
        The vector of initial values should not include values for the
        quadrature variables.
        + {b Set vector of quadrature variables }
        {[let yQ = Cvode.RealArray.of_array [| 0.0; 0.0 |] ]}
        The length of this vector determines the number of quadrature variables.
        + {b Initialize quadrature integration}
        {[init s fQ yQ]}
        + {b Specify integration tolerances (optional)}, e.g.
        {[set_tolerances s SStolerances (reltol, abstol)]}
        + {b Advance the solution in time as per normal}
        + {b Extract quadrature variables}
        {[get s yQ]}
        + {b Get quadrature optional outputs}
        {[let nre = get_num_rhs_evals s in ...]}
        Call any of the [get_*] functions to examine solver statistics.

        @cvodes <node5#SECTION00570000000000000000> Integration of pure quadrature equations *)

    (** {3:quadexcept Exceptions} *)

    (** Quadrature integration was not initialized.

        @cvodes <node5#ss:quad_get> CV_NO_QUAD *)
    exception QuadNotInitialized

    (** The quadrature right-hand side function failed in an unrecoverable
        manner.

        @cvodes <node5#SECTION00572000000000000000> CV_QRHSFUNC_FAIL *)
    exception QuadRhsFuncFailure

    (** The quadrature right-hand side function failed at the first call.

        @cvodes <node5#SECTION00572000000000000000> CV_FIRST_QRHSFUNC_ERR *)
    exception FirstQuadRhsFuncFailure

    (** Convergence test failures occurred too many times due to repeated
        recoverable errors in the quadrature right-hand side function. This
        value will also be returned if the quadrature right-hand side function
        had repeated recoverable errors during the estimation of an initial
        step size (assuming the quadrature variables are included in the error
        tests).

        @cvodes <node5#SECTION00572000000000000000> CV_REPTD_QRHSFUNC_ERR *)
    exception RepeatedQuadRhsFuncFailure

    (** The quadrature right-hand function had a recoverable error, but no
        recovery was possible. This failure mode is rare, as it can occur only
        if the quadrature right-hand side function fails recoverably after an
        error test failed while at order one.

        @cvodes <node5#SECTION00572000000000000000> CV_UNREC_QRHSFUNC_ERR *)
    exception UnrecoverableQuadRhsFuncFailure

    (** {3:quadinit Initialization} *)

    (** This function, [f_Q t y yQdot], computes the quadrature equation
        right-hand side [yQdot] for a given value of the independent variable
        [t] and state vector [y].

        @cvodes <node5#ss:user_fct_quad> CVQuadRhsFn *)
    type 'a quadrhsfn = float -> 'a -> 'a -> unit

    (** Activates the integration of quadrature equations.

        @cvodes <node5#ss:quad_malloc> CVodeQuadInit *)
    val init : ('a, 'k) session -> 'a quadrhsfn -> ('a, 'k) nvector -> unit

    (** Reinitialize the integration of quadrature equations.

        @cvodes <node5#ss:quad_malloc> CVodeQuadReInit *)
    val reinit : ('a, 'k) session -> ('a, 'k) nvector -> unit

    (** {3:quadtol Tolerance specification} *)

    type ('a, 'k) tolerance =
        NoStepSizeControl
        (** Do not use quadrature variables for step-size control (default). *)
      | SStolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('a, 'k) nvector
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

    (** Specify whether and how quadrature variables should be used in the step
        size control mechanism.

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
    val set_tolerances : ('a, 'k) session -> ('a, 'k) tolerance -> unit

    (** {3:quadout Output Functions} *)

    (** [tret = get s yq] fills [yq] with the quadrature solution vector after a
        successful return from {!Cvode.solve_normal} or
        {!Cvode.solve_one_step}, and returns the time reached by the
        solver.

        @cvodes <node5#ss:quad_get> CVodeGetQuad *)
    val get : ('a, 'k) session -> ('a, 'k) nvector -> float

    (** [tret = get_dky s t k dkyq] fills [dkyq] with the derivatives of the
        quadrature solution vector after a successful return from
        {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The
        time requested, [t], must fall within the interval defined by the last
        successful step ({!Cvode.get_last_step}). The requested order,
        [k], must be less than or equal to the value returned by
        {!Cvode.get_last_order}.

        @cvodes <node5#ss:quad_get> CVodeGetQuadDky
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit

    (** {3:quadoptout Optional Output Functions} *)

    (** Returns the number of calls to the user's quadrature right-hand side
        function.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals *)
    val get_num_rhs_evals       : ('a, 'k) session -> int

    (** Returns the number of local error test failures due to quadrature
        variables.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails *)
    val get_num_err_test_fails  : ('a, 'k) session -> int

    (** Returns the quadrature error weights at the current time.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights *)
    val get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit

    (** [nfqevals, nqetfails = get_stats s] returns
        - [fqevals], the number of calls to the user's quadrature function, and,
        - [nqetfails], the number of error test failures due to quadrature
          variables.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats *)
    val get_stats : ('a, 'k) session -> int * int

  end

(** {2:sens (Forward) Sensitivity Analysis} *)

module Sensitivity :
  sig
    (** A skeleton of an enhanced main program:
        + {b Initialize a session [s] per the skeleton at
           {!Cvode.session} or {!Quadrature.init}}
        {[...]}
        + {b Define the sensitivity problem}
        {[let p = Cvode.RealArray.make np in
let sp = { pvals = Some p; pbar = ...; plist = ... }]}
        + {b Set sensitivity initial conditions }
        {[let yS0 = Array.init ns (fun _ -> RealArray.init neq 0.0) in]}
        + {b Activate sensitivity calculations}
        {[init s (SStolerances ...) Simultaneous sp fS yS0]}
        + {b Set optional inputs}
        {[set_dq_method s ...]}
        + {b Advance the solution in time as per normal}
        + {b Extract sensitivity solution}
        {[let t = get s yS]}

        @cvodes <node6#ss:forward_usage> Enhanced skeleton for sensitivity analysis *)

    (** {3:sensexcept Exceptions} *)

    (** Sensitivity analysis was not initialized.

        @cvodes <node5#ss:sensi_get> CV_NO_SENS *)
    exception SensNotInitialized

    (** @cvodes <node6#SECTION00623000000000000000> CV_SRHSFUNC_FAIL *)
    exception SensRhsFuncFailure

    (** @cvodes <node6#SECTION00623000000000000000> CV_FIRST_SRHSFUNC_ERR *)
    exception FirstSensRhsFuncFailure

    (** @cvodes <node6#SECTION00623000000000000000> CV_REPTD_SRHSFUNC_ERR *)
    exception RepeatedSensRhsFuncFailure

    (** @cvodes <node6#SECTION00623000000000000000> CV_UNREC_SRHSFUNC_ERR *)
    exception UnrecoverableSensRhsFuncFailure

    (** @cvodes <node6> CV_BAD_IS *)
    exception BadSensIdentifier

    (** {3:sensinit Initialization} *)

    type 'a sensrhsfn =
        AllAtOnce of 'a sensrhsfn_all option
        (** Computes the sensitivity right-hand side for all
            sensitivity equations at once.  See {!sensrhsfn_all} for
            details.

            @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn
            @cvodes <node6#ss:sensi_malloc> CVodeSensInit
          *)
      | OneByOne of 'a sensrhsfn1 option
        (** Computes the sensitivity right-hand side one sensitivity
            parameter at a time.  See {!sensrhsfn1} for details.

            @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn
            @cvodes <node6#ss:sensi_malloc> CVodeSensInit1
          *)

    (** This function, [fS t y ydot yS ySdot tmp1 tmp2], computes the
        sensitivity right-hand side for all sensitivity equations at
        once, given
        - [t], the current value of the independent variable,
        - [y], the current value of the state vector,
        - [ydot], the current value of the right-hand side of the
        state equations,
        - [yS], the current values of the sensitivity vectors,
        - [ySdot], the sensitivity right-hand side vectors must be stored
        here,
        - [tmp1], and [tmp2] can be used as temporary storage.

        If a function is not given ([None]) then the default internal
        difference quotient sensitivity right-hand side routine is used.

        See also {!sensrhsfn}.

        @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit
      *)
    and 'a sensrhsfn_all =
      float           (* t *)
      -> 'a          (* y *)
      -> 'a          (* ydot *)
      -> 'a array    (* yS *)
      -> 'a array    (* ySdot *)
      -> 'a          (* tmp1 *)
      -> 'a          (* tmp2 *)
      -> unit

    (** This function, [fS t y ydot iS yS ySdot tmp1 tmp2], computes the
        sensitivity right-hand side one sensitivity parameter at a time,
        given
        - [t], the current value of the independent variable,
        - [y], the current value of the state vector,
        - [ydot], the current value of the right-hand side of the
        state equations,
        - [iS], the index of the parameter for which the sensitivity
        right-hand side must be computed,
        - [yS], the current value of the [iS]th sensitivity vector,
        - [ySdot], the [iS]th sensitivity right-hand side vector must be
        stored here,
        - [tmp1], and [tmp2] can be used as temporary storage.

        If a function is not given ([None]) then the default internal
        difference quotient sensitivity right-hand side routine is used.

        See also {!sensrhsfn}.

        @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit1
      *)
    and 'a sensrhsfn1 =
      float           (* t *)
      -> 'a          (* y *)
      -> 'a          (* ydot *)
      -> int         (* iS *)
      -> 'a          (* yS *)
      -> 'a          (* ySdot *)
      -> 'a          (* tmp1 *)
      -> 'a          (* tmp2 *)
      -> unit

    (** Specifies a sensitivity solution method.

        @cvodes <node6#ss:sensi_malloc> CVodeSensInit
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit1
      *)
    type sens_method =
        Simultaneous
        (** Correct state and sensitivity variables at the same time.
            If [Newton] was selected as the nonlinear system
            solution method, this amounts to performing a modified Newton
            iteration on the combined nonlinear system (CV_SIMULTANEOUS). *)
      | Staggered
        (** The correction step for the sensitivity variables takes place at the
            same time for all sensitivity equations, but only after the
            correction of the state variables has converged and the state
            variables have passed the local error test (CV_STAGGERED). *)
      | Staggered1
        (** All corrections are done sequentially, first for the state variables
            and then for the sensitivity variables, one parameter at a time. If
            the sensitivity variables are not included in the error control,
            this approach is equivalent to [Staggered]. Note that this
            approach can only be used if the user-provided sensitivity
            right-hand side function is of type [OneByOne] (CV_STAGGERED1). *)

    (** Used for specifying problem parameter information for
        sensitivity calculations.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensParams *)
    type sens_params = {
        pvals  : Sundials.RealArray.t option;
        (** The parameters used to evaluate {i f(t, y, p)}. *)
        pbar   : Sundials.RealArray.t option;
        (** An array of {i ns} positive scaling factors. *)
        plist  : int array option;
        (** An array of non-negative indices to specify which components
            to use in estimating the sensitivity equations. *)
      }

    val no_sens_params : sens_params

    type ('a, 'k) tolerance =
        SStolerances of float * Sundials.RealArray.t
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('a, 'k) nvector array
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
      | EEtolerances
        (** Calculate the integration tolerances for sensitivities
            based on those for state variables and the scaling factors
            (i.e. [pbar] in {!sens_params}). *)

    (** This function, [init s tol ism ps fS yS0], activates the
        forward sensitivity computation, where [ism] selects the
        sensitivity solution method, [ps] gives problem parameter
        information, [fS] computes the sensitivity right-hand sides,
        and [yS0] gives the initial values of the sensitivities.

        Note that any array specified by [ps.pvals] is used to pass
        parameter information during problem solution, that is, the
        library will normally write to it from time to time.

        @cvodes <node6#ss:sensi_malloc> CVodeSensInit
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit1 *)
    val init : ('a, 'k) session
               -> ('a, 'k) tolerance
               -> sens_method
               -> sens_params
               -> 'a sensrhsfn
               -> ('a, 'k) nvector array
               -> unit

    (** This function reinitializes the forward sensitivity computation.

        @cvodes <node6#ss:sensi_malloc> CVodeSensReInit *)
    val reinit : ('a, 'k) session -> sens_method
                      -> ('a, 'k) nvector array -> unit

    (** Deactivates forward sensitivity calculations without deallocating
        memory. Sensitivities can be reactivated with {!reinit}.

        @cvodes <node6#ss:sensi_malloc> CVodeSensToggleOff *)
    val toggle_off : ('a, 'k) session -> unit

    (** {3:sensout Output Functions} *)

    (** [tret = get s ys] fills [ys] with the sensitivity solution vectors after
        a successful return from {!Cvode.solve_normal} or
        {!Cvode.solve_one_step}, and returns the time reached by the solver.

        @cvodes <node6#ss:sensi_get> CVodeGetSens *)
    val get : ('a, 'k) session -> ('a, 'k) nvector array -> float

    (** [tret = get_dky s t k dkys] fills [dkys] with the
        derivatives of the sensitivity solution vectors after a
        successful return from {!Cvode.solve_normal} or
        {!Cvode.solve_one_step}. The time requested, [t], must fall
        within the interval defined by the last successful step
        ({!Cvode.get_last_step}). The requested order, [k], must be
        less than or equal to the value returned by
        {!Cvode.get_last_order}.

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range.
    *)
    val get_dky : ('a, 'k) session -> float -> int
                          -> ('a, 'k) nvector array -> unit

    (** [tret = get1 s i ys] fills [ys] with the [i]th sensitivity
        solution vector after a successful return from
        {!Cvode.solve_normal} or {!Cvode.solve_one_step}, and returns
        the time reached by the solver.

        @cvodes <node6#ss:sensi_get> CVodeGetSens1
        @raise BadIS The index [i] is not in the allowed range.
    *)
    val get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float

    (** [tret = get_dky1 s t k i dkys] fills [dkys] with the
        derivatives of the [i]th sensitivity solution vector after a
        successful return from {!Cvode.solve_normal} or
        {!Cvode.solve_one_step}. The time requested, [t], must fall
        within the interval defined by the last successful step
        ({!Cvode.get_last_step}). The requested order, [k], must be
        less than or equal to the value returned by
        {!Cvode.get_last_order}.

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky1
        @raise BadIS The index [i] is not in the allowed range.
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky1 : ('a, 'k) session -> float -> int -> int
                          -> ('a, 'k) nvector -> unit

    (** {3:sensoptin Optional Input Functions} *)

    (** Specify the integration tolerances for sensitivities.

        {b NB}: Unlike the other [set_tolerances] functions in
        [Cvodes], this one does {b not} call {!set_err_con} (which
        defaults to [false]).

        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSStolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensSVtolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensEEtolerances *)
    val set_tolerances : ('a, 'k) session -> ('a, 'k) tolerance -> unit

    (** Set whether sensitivity variables should be used in the error control
        mechanism (the default is [false]). 

        @cvodes <node5#ss:sens_optional_input> CVodeSetSensErrCon *)
    val set_err_con : ('a, 'k) session -> bool -> unit

    type dq_method = DQCentered (** CV_CENTERED *)
                   | DQForward  (** CV_FORWARD *)

    (** [set_dq_method s dqtype dqrhomax] specifies the difference quotient
        strategy in the case in which the right-hand side of the sensitivity
        equations is to be computed by CVODES; [dqrhomax] is used in deciding
        the switching between simultaneous or separate approximations of the two
        terms in the sensitivity right-hand side.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensDQMethod
        @cvodes <node3#ss:fwd_sensi> Forward Sensitivity Analysis *)
    val set_dq_method : ('a, 'k) session -> dq_method -> float -> unit

    (** Specifies the maximum number of nonlinear solver iterations for
        sensitivity variables permitted per step.

        @cvode <node5#ss:sens_optional_input> CVodeSetSensMaxNonlinIters *)
    val set_max_nonlin_iters : ('a, 'k) session -> int -> unit

    (** {3:sensoptout Optional Output Functions} *)

    (** Returns the number of calls to the sensitivity right-hand side function.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumRhsEvals *)
    val get_num_rhs_evals       : ('a, 'k) session -> int

    (** Returns the number of calls to the user's right-hand side function due
        to the internal finite difference approximation of the sensitivity
        right-hand sides.

        @cvode <node6#ss:sens_optional_output> CVodeGetNumRhsEvalsSens *)
    val get_num_rhs_evals_sens  : ('a, 'k) session -> int

    (** Returns the number of local error test failures for the sensitivity
        variables that have occurred.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumErrTestFails *)
    val get_num_err_test_fails  : ('a, 'k) session -> int

    (** Returns the number of calls made to the linear solver's setup function
        due to forward sensitivity calculations.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumLinSolvSetups *)
    val get_num_lin_solv_setups : ('a, 'k) session -> int

    type sensitivity_stats = {
        num_rhs_evals : int;
        num_sens_evals :int;
        num_err_test_fails : int;
        num_lin_solv_setups :int;
      }

    (** Returns all of the sensitivity-related solver statistics as a group.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensStats *)
    val get_stats : ('a, 'k) session -> sensitivity_stats

    (** Returns the sensitivity error weight vectors at the current time.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensErrWeights
        @cvode <node3#e:errwt> Eq. (2.7) IVP solution (W_i) *)
    val get_err_weights : ('a, 'k) session -> ('a, 'k) nvector array -> unit

    (** Returns the number of nonlinear iterations performed for sensitivity
        calculations.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvIters *)
    val get_num_nonlin_solv_iters : ('a, 'k) session -> int

    (** Returns the number of nonlinear convergence failures that have occurred
        for sensitivity calculations.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvConvFails *)
    val get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int

    (** [nni, ncfn = get_nonlin_solv_stats s] returns the sensitivity-related
        nonlinear solver statistics as a group, where [nni] is the number of
        nonlinear iterations performed for sensitivity calculations, and [ncfn]
        is the number of nonlinear convergence failures that have occurred for
        sensitivity calculations.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNonlinSolvStats
     *)
    val get_nonlin_solv_stats : ('a, 'k) session -> int * int

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed for each sensitivity equation separately, in the [Staggered1]
        case.

      @cvode <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvIters *)
    val get_num_stgr_nonlin_solv_iters : ('a, 'k) session
                                         -> Sundials.LintArray.t -> unit

    (** Returns the number of nonlinear convergence failures that have occurred
        for each sensitivity equation separately, in the [Staggered1] case.

      @cvode <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvConvFails *)
    val get_num_stgr_nonlin_solv_conv_fails : ('a, 'k) session
                                              -> Sundials.LintArray.t -> unit

    (** {2:quadsens Quadrature Equations} *)

    (**
       Support for integration of quadrature equations that depends not only on
       the state variables but also on forward sensitivities.
     *)

    module Quadrature :
      sig
        (** A skeleton of an enhanced main program:
            + {b Initialize a session [s] per the skeleton at {!Sensitivity.init}}
            {[...]}
            + {b Set initial values of quadrature variables}
            {[let yQS = Array.init ns (fun _ -> RealArray.of_array [0.0; 0.0]) in]}
            + {b Initialize sensitivity-dependent quadrature integration}
            {[init s fQS yQS]}
            + {b Set optional inputs}
            {[set_tolerances s ...]}
            + {b Advance the solution in time as per normal}
            + {b Extract sensitivity-dependent quadrature variables}
            {[let t = get s yQS]}
            + {b Get sensitivity-dependent optional outputs}
            {[let e, f = get_stats s]}

            @cvodes <node6#SECTION00640000000000000000> Integration of quadrature equations depending on forward sensitivities *)

        (** {3:quadsensexcept Exceptions} *)

        (** Quadrature integration was not initialized.

            @cvodes <node5#SECTION00642000000000000000> CV_NO_QUAD_SENS *)
        exception QuadSensNotInitialized

        (** @cvodes <node6#SECTION00642000000000000000> CV_QSRHSFUNC_FAIL *)
        exception QuadSensRhsFuncFailure

        (** @cvodes <node6#SECTION00642000000000000000> CV_FIRST_QSRHSFUNC_ERR *)
        exception FirstQuadSensRhsFuncFailure

        (** @cvodes <node6#SECTION00642000000000000000> CV_REPTD_QSRHSFUNC_ERR *)
        exception RepeatedQuadSensRhsFuncFailure

        (** @cvodes <node6#SECTION00642000000000000000> CV_UNREC_QSRHSFUNC_ERR *)
        exception UnrecoverableQuadSensRhsFuncFailure

        (** {3:quadsensinit Initialization} *)

        (** This function, [fQS t y yS yQdot rhsvalQs tmp1 tmp2], computes the
            sensitivity quadrature equation right-hand side given
            - [t], the current value of the independent variable,
            - [y], the current value of the state vector,
            - [yS], an array of dependent sensitivity vectors,
            - [yQdot], the current value of the quadrature right-hand side,
            - [rhsvalQs], the right-hand side vectors must be stored here,
            - [tmp1], and [tmp2] can be used as temporary storage.

           @cvodes <node6#ss:user_fct_quad_sens> CVodeQuadSensRhsFn *)
        type 'a quadsensrhsfn =
           float          (* t *)
           -> 'a          (* y *)
           -> 'a array    (* yS *)
           -> 'a          (* yQdot *)
           -> 'a array    (* rhsvalQs *)
           -> 'a          (* tmp1 *)
           -> 'a          (* tmp2 *)
           -> unit

        (** This function, [init s fQS yQS0], activates integration of
            quadrature equations depending on sensitivities, where
            [fQS] computes the right-hand side of the
            sensitivity-dependent quadrature equations, and [yQS0]
            contains the initial values of sensitivity-dependent
            quadratures.  When no [fQB] is supplied, the solver uses
            an internal implementation based on difference quotients.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensInit *)
        val init : ('a, 'k) session -> 'a quadsensrhsfn option
                 -> ('a, 'k) nvector array -> unit

        (** This function reinitializes the forward sensitivity computation.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensReInit *)
        val reinit : ('a, 'k) session -> ('a, 'k) nvector array -> unit

        (** {3:quadsenstol Tolerance specification} *)

        type ('a, 'k) tolerance =
            NoStepSizeControl
            (** Do not use quadrature variables for step-size control
                (default). *)
          | SStolerances of float * Sundials.RealArray.t
            (** [(rel, abs)] : scalar relative and absolute tolerances. *)
          | SVtolerances of float * ('a, 'k) nvector array
            (** [(rel, abs)] : scalar relative and vector absolute
                tolerances. *)
          | EEtolerances
            (** Estimate the tolerances for the sensitivity-dependent
                quadratures from those provided for the pure quadrature
                variables. *)

        (** Specify whether and how quadrature variables should be used in the
            step size control mechanism.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeSetQuadSensErrCon
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSStolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSVtolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensEEtolerances *)
        val set_tolerances : ('a, 'k) session -> ('a, 'k) tolerance -> unit

        (** {3:quadsensout Output Functions} *)

        (** [tret = get s yqs] fills [yqs] with quadrature solution vectors
            after a successful return from {!Cvode.solve_normal} or
            {!Cvode.solve_one_step}, and returns the time reached by
            the solver.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens *)
        val get : ('a, 'k) session -> ('a, 'k) nvector array -> float

        (**
          [tret = get s i yqs] fills [yqs] with the [i]th quadrature solution
          vector after a successful return from {!Cvode.solve_normal}
          or {!Cvode.solve_one_step}, and returns the time reached by
          the solver.

          @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens1
          @raise BadIS The index [i] is not in the allowed range. *)
        val get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float

        (**
          [tret = get_dky s t k dkyqs] fills [dkyqs] with the derivatives of the
          quadrature solution vectors after a successful return from
          {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The
          time requested, [t], must fall within the interval defined by the last
          successful step ({!Cvode.get_last_step}). The requested order,
          [k], must be less than or equal to the value returned by
          {!Cvode.get_last_order}.

          @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky
          @raise BadIS The index is not in the allowed range.
          @raise BadK [k] is not in the range 0, 1, ..., [qlast].
          @raise BadT [t] is not in the allowed range. *)
        val get_dky : ('a, 'k) session -> float -> int
                          -> ('a, 'k) nvector array -> unit

        (** [tret = get_dky s t k i dkyqs] fills [dkyqs] with the derivatives of
            the [i]th quadrature solution vector after a successful return from
            {!Cvode.solve_normal} or {!Cvode.solve_one_step}.
            The time requested, [t], must fall within the interval defined by
            the last successful step ({!Cvode.get_last_step}). The
            requested order, [k], must be less than or equal to the value
            returned by {!Cvode.get_last_order}.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky1
            @raise BadK [k] is not in the range 0, 1, ..., [qlast].
            @raise BadT [t] is not in the allowed range. *)
        val get_dky1 : ('a, 'k) session -> float -> int -> int
                            -> ('a, 'k) nvector -> unit

        (** {3:quadsensoptout Optional Output Functions} *)

        (** Returns the number of calls to the user's quadrature right-hand side
            function.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumRhsEvals *)
        val get_num_rhs_evals       : ('a, 'k) session -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumErrTestFails *)
        val get_num_err_test_fails  : ('a, 'k) session -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensErrWeights *)
        val get_err_weights : ('a, 'k) session -> ('a, 'k) nvector array -> unit

        (** [nfqevals, nqetfails = get_stats s] returns
            - [fqevals], the number of calls to the user's quadrature function,
            and,
            - [nqetfails], the number of error test failures due to quadrature
            variables.

          @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensStats *)
        val get_stats : ('a, 'k) session -> int * int

      end
  end

(** {2:adj Adjoint Sensitivity Analysis} *)

module Adjoint :
  sig
    (**
        A skeleton of an enhanced main program:
        + {b Initialize a session [s] per the skeleton at {!Cvode.init}}
        {[...]}
        Adding quadrature variables using {!Quadrature.init} if desired.
        + {b Initialize the adjoint computation}
        {[init s nsteps IHermite]}
        + {b Integrate forward problem}
        {[let t, ncheck, r = forward_normal s tout y0]}
        + {b Setup the backward problem and attach a linear solver}
        {[let yB0 = RealArray.of_list [0.0; 0.0; ...]
let bs = init_backward s lmm (Newton ...) (SStolerances ...) fB tB0 yB0]}
        + {b Set optional inputs}
        {[set_max_ord bs ...]}
        + {b Initialize quadrature calculation}
        {[Quadrature.init bs fQb yQB0]}
        + {b Integrate backward problem}
        {[backward_normal s tB]}
        + {b Extract quadrature variables}
        {[let t = Quadrature.get s yQS]}

        @cvodes <node7#ss:skeleton_adj> Enhanced Skeleton for Adjoint Sensitivity Analysis
    *)

    (** {3:adjexcept Exceptions} *)

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

    (** The final time [tB0] was outside the interval over which the forward
        problem was solved.

        @cvodes <node7#sss:cvinitb> CV_BAD_TB0 *)
    exception BadFinalTime

    (** {3:adjfwd Forward Solutions} *)

    (** {4:adjfwdinit Initialization} *)

    (** Specifies the type of interpolation.

        @cvodes <node3#ss:checkpointing> Checkpointing scheme *)
    type interpolation = IPolynomial (** CV_POLYNOMIAL *)
                       | IHermite    (** CV_HERMITE *)

    (** [init s nd interp] initializes the forward-backward problem with [nd]
        integration steps between consecutive checkpoints and variable-degree
        interpolation according to [interp]. This function must be called before
        either {!forward_normal} or {!forward_one_step}.

        @cvodes <node7#ss:cvadjinit> CVodeAdjInit *)
    val init : ('a, 'k) session -> int -> interpolation -> unit

    (** {4:adjfwdintegration Forward Integration} *)

    (** [tret, ncheck, sr = forward_normal s tout yret] integrates the forward
        problem over an interval and saves checkpointing data. The function
        takes as arguments the next time at which a solution is desired
        ([tout]), a vector for storing the computed result ([yret]), and returns
        the time reached by the solver ([tret]), the number of checkpoints
        stored so far ([ncheck]), and whether the solver reached [tout] or not
        ([sr]).

        This call asks the solver to take internal steps until it has reached or
        just passed the [tout] parameter ([CV_NORMAL]). The solver then
        interpolates in order to return an approximate value of [y(tout)].

        @cvodes <node7#sss:cvsolvef> CVodeF
        @raise Cvode.IllInput           One of the inputs is invalid.
        @raise Cvode.TooMuchWork        Could not reach [tout] in [mxstep] steps
        @raise Cvode.TooMuchAccuracy    Could not satisfy the demanded accuracy
        @raise Cvode.ErrFailure         Too many error test failures.
        @raise Cvode.ConvergenceFailure Too many convergence test failures.
        @raise Cvode.LinearSetupFailure Unrecoverable failure in linear solver setup function.
        @raise Cvode.LinearSolveFailure Unrecoverable failure in linear solver solve function.
        @raise AdjointNotInitialized    The [init] function has not previously been called. *)
    val forward_normal :
      ('a, 'k) session
      -> float
      -> ('a, 'k) nvector
      -> float * int * Cvode.solver_result

    (** [tret, ncheck, sr = forward_normal s tout yret] integrates the forward
        problem over an interval and saves checkpointing data. The function
        takes as arguments the next time at which a solution is desired
        ([tout]), a vector for storing the computed result ([yret]), and returns
        the time reached by the solver ([tret]) and the number of checkpoints
        stored so far ([ncheck]), and whether the solver reached [tout] or not
        ([sr]).

        This call asks the solver to take one internal step and to return the
        solution at the point reached by that step ([CV_ONE_STEP]).

        @cvodes <node7#sss:cvsolvef> CVodeF
        @raise Cvode.IllInput           One of the inputs is invalid.
        @raise Cvode.TooMuchWork        Could not reach [tout] in [mxstep] steps
        @raise Cvode.TooMuchAccuracy    Could not satisfy the demanded accuracy
        @raise Cvode.ErrFailure         Too many error test failures.
        @raise Cvode.ConvergenceFailure Too many convergence test failures.
        @raise Cvode.LinearSetupFailure Unrecoverable failure in linear solver setup function.
        @raise Cvode.LinearSolveFailure Unrecoverable failure in linear solver solve function.
        @raise AdjointNotInitialized    The [init] function has not previously been called. *)
    val forward_one_step :
      ('a, 'k) session
      -> float
      -> ('a, 'k) nvector
      -> float * int * Cvode.solver_result

    (** {3:adjbwd Backward Problems} *)

    (** Identifies a backward problem. *)
    type ('data, 'kind) bsession = ('data, 'kind) AdjointTypes.bsession
    type serial_bsession = (RealArray.t, Nvector_serial.kind) bsession

    (** {4:adjbwdinit Initialization} *)


    (** These functions evaluate the right-hand side of the backward ODE system
        with or without a dependence on forward sensitivities. *)
    type 'a brhsfn =
        NoSens of 'a brhsfn_no_sens
        (** Doesn't depend on forward sensitivities.  See
            {!brhsfn_no_sens} for details.  *)
      | WithSens of 'a brhsfn_with_sens
        (** Depends on forward sensitivities.  See {!brhsfn_with_sens}
            for details.  *)

    (** Backward rhs function that doesn't depend on forward sensitivities.

        See also {!brhsfn}.

        @cvodes <node7#ss:ODErhs_b> CVRhsFnB
        @cvodes <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis *)
    and 'a brhsfn_no_sens =
      float    (* t *)
      -> 'a    (* y *)
      -> 'a    (* yb *)
      -> 'a    (* ybdot *)
      -> unit

    (** Backward rhs function that depends on forward sensitivities.

        See also {!brhsfn}.

        @cvodes <node7#ss:ODErhs_bs> CVRhsFnBS
        @cvodes <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis *)
    and 'a brhsfn_with_sens =
      float        (* t *)
      -> 'a        (* y *)
      -> 'a array  (* ys *)
      -> 'a        (* yb *)
      -> 'a        (* ybdot *)
      -> unit

    type 'a single_tmp = 'a
    type 'a triple_tmp = 'a * 'a * 'a

    (** Arguments common to all Jacobian callback functions.

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
        jac_fyb : 'a;           (** The backward right-hand side function fB. *)
        jac_tmp : 't            (** Workspace data,
                                    either {!single_tmp} or {!triple_tmp}. *)
      }

    type bandrange = Cvode_impl.bandrange =
      { mupper : int; (** The upper half-bandwidth.  *)
        mlower : int; (** The lower half-bandwidth.  *) }

    (** Specify a linear solver.

        The Lapack solvers require that both Sundials and the OCaml interface
        were built to link with a LAPACK library.

        @cvodes <node7#sss:lin_solv_b> Linear Solver Initialization Functions *)
    type ('data, 'kind) linear_solver = ('data, 'kind) AdjointTypes.linear_solver
    type serial_linear_solver = (RealArray.t, Nvector_serial.kind) linear_solver

    type ('a, 'k) tolerance =
      | SStolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('a, 'k) nvector
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

    (** Specify a solution method.

        @cvodes <node7#sss:cvinitb> CVodeCreateB
     *)
    type ('data, 'kind) iter =
      | Newton of ('data, 'kind) linear_solver
        (** Newton iteration with a given linear solver *)
      | Functional
        (** Functional iteration (non-stiff systems only) *)

    (** [init_backward s lmm iter tol fB tB0 yB0] adds and initializes a
        backward problem that may or may not depend on forward sensitivities,
        where
        - [s] is the parent session (going forward),
        - [lmm]     specifies the linear multistep method, see {!Cvode.lmm},
        - [iter]    specifies either functional iteration or Newton iteration
                    with a specific linear solver, see {!iter},
        - [tol]     specifies the tolerances, see {!tolerance},
        - [fB]      computes the right-hand side of the backward ODE problem,
        - [tB0]     specifies the endpoint where final conditions are provided
                    for the backward problem, normally equal to the endpoint of
                    the forward integration, and,
        - [yB0]     is the final value of the backward problem.

        @cvodes <node7#sss:cvinitb> CVodeCreateB
        @cvodes <node7#sss:cvinitb> CVodeInitB
        @cvodes <node7#sss:cvinitb> CVodeInitBS
        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB 
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise BadFinalTime      The final time is outside the interval over which the forward problem was solved. *)
    val init_backward :
         ('a, 'k) session
      -> Cvode.lmm
      -> ('a, 'k) iter
      -> ('a, 'k) tolerance
      -> 'a brhsfn
      -> float
      -> ('a, 'k) nvector
      -> ('a, 'k) bsession

    (** Reinitialize the backward problem.

        @cvodes <node7#sss:cvinitb> CVodeReInitB
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise BadFinalTime      The final time is outside the interval over which the forward problem was solved. *)
    val reinit : ('a, 'k) bsession -> float -> ('a, 'k) nvector -> unit

    (** {4:adjbwdintegration Backward Integration} *)

    (** [backward_normal s tbout] integrates the backward ODE
        problem. The function takes internal steps until it has
        reached or just passed the user-specified value [tbout]
        ([CV_NORMAL]). The solver then interpolates in order to return
        an approximate value of [y(tbout)] when {!get} is called.

        @cvodes <node7#sss:cvsolveb> CVodeB
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise NoBackwardProblem        The [init_backward] function has not previously been called.
        @raise NoForwardCall            Neither [forward_normal] nor [forward_one_step] has previously been called.
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
    val backward_normal : ('a, 'k) session -> float -> unit

    (** [backward_one_step s tbout] integrates the backward ODE problem. The
        function takes one internal step ([CV_ONE_STEP]).

        @cvodes <node7#sss:cvsolveb> CVodeB
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise NoBackwardProblem        The [init_backward] function has not previously been called.
        @raise NoForwardCall            Neither [forward_normal] nor [forward_one_step] has previously been called.
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
    val backward_one_step : ('a, 'k) session -> float -> unit

    (** [tret = get bs yb] returns the solution of the backward ODE problem
        in [yb] at time [tret].

        @cvodes <node7#sss:cvsolveb> CVodeGetB *)
    val get : ('a, 'k) bsession -> ('a, 'k) nvector -> float

    (** [tret = get_dky s t k dkys] fills [dkys] with the derivatives of the
        sensitivity solution vectors after a successful return from
        {!backward_normal} or {!backward_one_step}. The time requested, [t],
        must fall within the interval defined by the last successful step
        ({!get_last_step}). The requested order, [k], must be less than or equal
        to the value returned by {!get_last_order}.

        @cvodes <node5#ss:optional_dky> CVodeGetDky
        @cvodes <node5#ss:optional_dky> CVodeGetAdjIDABmem
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky : ('a, 'k) bsession -> float -> int -> ('a, 'k) nvector -> unit

    (** {4:adjbwdoptout Optional Output Functions} *)

    (** Instructs {!forward_normal} and {!forward_one_step} not to save
        checkpointing data for forward sensitivities anymore.

        @cvodes <node7#SECTION00727000000000000000> CVodeAdjSetNoSensi *)
    val set_no_sensitivity : ('a, 'k) session -> unit

    (** Specify the integration tolerances for the backward problem.

        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB *)
    val set_tolerances : ('a, 'k) bsession -> ('a, 'k) tolerance -> unit

    (** Specifies the maximum order of the linear multistep method.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxOrdB *)
    val set_max_ord : ('a, 'k) bsession -> int -> unit

    (** Specifies the maximum number of steps to be taken by the solver in its
        attempt to reach the next output time.

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

    (** {4:adjbwddiag Diagonal Approximation} *)

    module Diag :
      sig
        (** Diagonal approximation of the Jacobian by difference quotients.

            @cvodes <node7#sss:lin_solv_b> CVDiagB *)
        val solver : ('data, 'kind) linear_solver

        (** Returns the number of calls made to the user-supplied right-hand
            side function due to finite difference Jacobian approximation in the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetWorkSpace
            @return ([real_size], [integer_size]) *)
        val get_work_space : ('a, 'k) bsession -> int * int

        (** Returns the number of calls made to the user-supplied right-hand
            side function due to finite difference Jacobian approximation in the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetNumRhsEvals *)
        val get_num_rhs_evals : ('a, 'k) bsession -> int
      end

    (** {4:adjbwddirect Direct Linear Solver} *)

    module Dls :
      sig
        (** This function computes the dense Jacobian of the backward problem
            (or an approximation to it).

            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
        type dense_jac_fn = (RealArray.t triple_tmp, RealArray.t) jacobian_arg
                                -> Dls.DenseMatrix.t -> unit

        (** This function computes the banded Jacobian of the backward problem
            (or an approximation to it).

            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        type band_jac_fn = bandrange
                            -> (RealArray.t triple_tmp, RealArray.t) jacobian_arg
                            -> Dls.BandMatrix.t -> unit

        (** Direct linear solver with dense matrix.  The optional argument
            specifies a callback function that computes an approximation to the
            Jacobian matrix (see {!dense_jac_fn} for details).  If this
            argument is [None], then CVODE uses a default implementation based
            on difference quotients.  See also {!Dls}.

            @cvodes <node7#sss:lin_solv_b> CVDenseB
            @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
        val dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

        (** Direct linear solver with dense matrix, using LAPACK.  The argument
            is the same as [Dense].  See also {!Dls}.

            @cvodes <node7#sss:lin_solv_b> CVLapackDenseB
            @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
            @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
        val lapack_dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

        (** Direct linear solver with banded matrix.  The arguments specify the
            width of the band ({!bandrange}) and an optional Jacobian
            function ({!band_jac_fn}).  If the Jacobian function is [None],
            CVODES uses an internal implementation based on difference
            quotients. See also {!Dls}.

            @cvodes <node7#sss:lin_solv_b> CVBandB
            @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        val band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

        (** Direct linear solver with banded matrix using LAPACK.  The arguments
            are the same as [Band].

            @cvodes <node7#sss:lin_solv_b> CVLapackBandB
            @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
            @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
        val lapack_band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

      end

    (** {4:adjbwdspils Scaled Preconditioned Iterative Linear Solvers (SPILS)} *)

    module Spils :
      sig
        (** Scaled Preconditioned Iterative Linear Solvers (SPILS)

        @cvodes <node7#ss:optional_output_b> Optional output functions for the backward problem. *)

        type gramschmidt_type = Spils.gramschmidt_type =
          | ModifiedGS
          | ClassicalGS

        (** Arguments passed to the preconditioner solve callback
            function.  See {!prec_solve_fn}.

            @cvode <node7#ss:psolve_b> CVSpilsPrecSolveFnB *)
        type 'a prec_solve_arg = 'a AdjointTypes.SpilsTypes.prec_solve_arg =
          {
            rhs   : 'a;        (** The right-hand side vector, {i r}, of the
                                   linear system. *)
            gamma : float;     (** The scalar {i g} appearing in the Newton
                                   matrix given by M = I - {i g}J. *)
            delta : float;     (** Input tolerance to be used if an iterative method
                                   is employed in the solution. *)

            left  : bool;      (** Indicates whether to use the left preconditioner
                                   ([true]) or the right one ([false]). *)
          }

        (** This function solves the preconditioning system {i Pz = r} for
            the backward problem.

            See also {!callbacks}.

            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB *)
        type 'a prec_solve_fn =
          ('a single_tmp, 'a) jacobian_arg
          -> 'a prec_solve_arg
          -> 'a
          -> unit

        (** This function preprocesses and/or evaluates Jacobian-related
            data needed by the preconditioner for the backward problem.

            See also {!callbacks}.

            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
        type 'a prec_setup_fn =
          ('a triple_tmp, 'a) jacobian_arg
          -> bool
          -> float
          -> bool

        (** This function computes the action of the Jacobian for the
            backward problem on a given vector.

            See also {!callbacks}.

            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB *)
        type 'a jac_times_vec_fn =
          ('a single_tmp, 'a) jacobian_arg
          -> 'a
          -> 'a
          -> unit

        (** A type that specifies a preconditioner, including the side
            from which to precondition (none, left, right, or both)
            and callbacks or other data needed to solve the
            preconditioning system $Pz = r$.

            Generally, three callbacks are needed in preconditioning:

            - [solve], the main function that solves the
              preconditioning system $Pz = r$, where $P$ is a
              preconditioning matrix chosen by the user.  See
              {!prec_solve_fn} for details.
            - [setup], which preprocesses and/or evaluates
              Jacobian-related data needed by [solve].  It can be
              omitted if there are no such data.  See {!prec_setup_fn}
              for details.
            - [jac_times_vec], which multiplies the system Jacobian to
              a given vector.  See {!jac_times_vec_fn} for details.
              If the user doesn't give such a function, CVODE uses a
              default implementation based on difference quotients.

            Like the {!linear_solver}, there are several functions
            which construct preconditioners.  The simples is
            {!prec_none}, which does no preconditioning.  This module
            also provides {!prec_left}, {!prec_right}, and
            {!prec_both}, which take a [solve] function and optionally
            [setup] and [jac_times_vec] which together define an
            arbitrary, user-defined preconditioner.

            The {!Banded} module gives access to CVODE's banded
            preconditioner, while {!Cvode_bbd} contains the parallel
            band-block diagonal preconditioners.


            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetJacTimesVecFnB
            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB
            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
          *)
        type ('a, 'k) preconditioner =
          ('a, 'k) AdjointTypes.SpilsTypes.preconditioner

        (** {!preconditioner} restricted to serial nvectors. *)
        type serial_preconditioner =
          (Nvector_serial.data, Nvector_serial.kind) preconditioner

        (** See {!preconditioner}.  *)
        val prec_none : ('a, 'k) preconditioner

        (** See {!preconditioner}. *)
        val prec_left :
          ?setup:'a prec_setup_fn
          -> ?jac_times_vec:'a jac_times_vec_fn
          -> 'a prec_solve_fn
          -> ('a, 'k) preconditioner

        (** See {!preconditioner}. *)
        val prec_right :
          ?setup:'a prec_setup_fn
          -> ?jac_times_vec:'a jac_times_vec_fn
          -> 'a prec_solve_fn
          -> ('a, 'k) preconditioner

        (** See {!preconditioner}. *)
        val prec_both :
          ?setup:'a prec_setup_fn
          -> ?jac_times_vec:'a jac_times_vec_fn
          -> 'a prec_solve_fn
          -> ('a, 'k) preconditioner

        (** Krylov iterative solver with the scaled preconditioned
            GMRES method.  The arguments specify the maximum dimension
            of the Krylov subspace and preconditioning type
            ({!Spils.preconditioning_type}) and the preconditioner
            callback functions ({!callbacks}).

            @cvodes <node7#sss:lin_solv_b> CVSpgmrB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
        val spgmr : ?maxl:int -> ('a, 'k) preconditioner
                      -> ('a, 'k) linear_solver

        (** Krylov iterative solver with the scaled preconditioned Bi-CGStab
            method. The arguments are the same as [Spgmr].

            @cvodes <node7#sss:lin_solv_b> CVSpbcgB
            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
        val spbcg : ?maxl:int -> ('a, 'k) preconditioner
                      -> ('a, 'k) linear_solver

        (** Krylov iterative with the scaled preconditioned TFQMR method.  The
            arguments are the same as [Spgmr].  See also {!Spils}. *)
        val sptfqmr : ?maxl:int -> ('a, 'k) preconditioner
                      -> ('a, 'k) linear_solver

        (** {5:llsolvermanip Low-level solver manipulation} *)

        (** This type is used only for low-level solver manipulation.
            Information about the type of preconditioning to be done,
            without any of the necessary callbacks to make it happen:
            [PrecNone], [PrecLeft ()], [PrecRight ()], or [PrecBoth ()].
        *)
        type preconditioning_type =
          | PrecNone
          | PrecLeft
          | PrecRight
          | PrecBoth

        (** This function changes the type of preconditioning without
            affecting the preconditioning callbacks.  If the
            preconditioning type is changed from [PrecNone] to
            something else, then {!set_prec_callbacks} must be called
            to install the necessary callbacks.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPrecTypeB
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val set_prec_type : ('a, 'k) bsession -> preconditioning_type -> unit

        (** {5:adjbwdspilsoptin Optional Input Functions} *)

        (** Sets the Gram-Schmidt orthogonalization to be used with the
            Spgmr {!linear_solver}.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetGSTypeB
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val set_gs_type : ('a, 'k) bsession -> gramschmidt_type -> unit

        (** [set_eps_lin eplifac] sets the factor by which the Krylov linear
            solver's convergence test constant is reduced from the Newton
            iteration test constant. [eplifac]  must be >= 0. Passing a value of
            0 specifies the default (which is 0.05).

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetEpsLinB
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val set_eps_lin : ('a, 'k) bsession -> float -> unit

        (** [set_maxl maxl] resets the maximum Krylov subspace dimension for the
            Bi-CGStab or TFQMR methods. [maxl] is the maximum dimension of the
            Krylov subspace.  A value of [None] (or [maxl] <= 0) specifies the
            default of 5.0.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val set_maxl : ('a, 'k) bsession -> int option -> unit

        (** {5:adjbwdspilsoptout Optional Output Functions} *)

        (** Returns the sizes of the real and integer workspaces used by the
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

        (** Returns the number of preconditioner evaluations, i.e., the number
            of calls made to psetup with jok = [false].

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_prec_evals   : ('a, 'k) bsession -> int

        (** Returns the cumulative number of calls made to the preconditioner
            solve function, psolve.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecSolves
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_prec_solves  : ('a, 'k) bsession -> int

        (** Returns the cumulative number of calls made to the Jacobian-vector
            function, jtimes.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumJtimesEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_jtimes_evals : ('a, 'k) bsession -> int

        (** Returns the number of calls to the user right-hand side function for
            finite difference Jacobian-vector product approximation. This
            counter is only updated if the default difference quotient function
            is used.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
        val get_num_rhs_evals    : ('a, 'k) bsession -> int

        module Banded : sig

          (** Denotes CVODE's internal band matrix {!preconditioner} based
              on difference quotients.  [prec_left br] creates a left
              preconditioner which generates a banded approximation to the
              Jacobian with [br.mlower] sub-diagonals and [br.mupper]
              super-diagonals.

              @cvode <node5#sss:cvbandpre> CVBandPrecInit
          *)
          val prec_left :
            ?jac_times_vec:(RealArray.t jac_times_vec_fn)
            -> bandrange
            -> serial_preconditioner

          (** Like {!prec_left} but preconditions from the right.

              @cvode <node5#sss:cvbandpre> CVBandPrecInit
          *)
          val prec_right :
            ?jac_times_vec:(RealArray.t jac_times_vec_fn)
            -> bandrange
            -> serial_preconditioner

          (** Like {!prec_left} but preconditions from both sides.

              @cvode <node5#sss:cvbandpre> CVBandPrecInit
          *)
          val prec_both :
            ?jac_times_vec:(RealArray.t jac_times_vec_fn)
            -> bandrange
            -> serial_preconditioner

          (** Returns the sizes of the real and integer workspaces
              used by the serial banded preconditioner module.

              @cvodes <node5#sss:cvbandpre> CVBandPrecGetWorkSpace
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
              @return ([real_size], [integer_size]) *)
          val get_work_space : serial_bsession -> int * int

          (** Returns the number of calls made to the user-supplied
              right-hand side function due to finite difference banded
              Jacobian approximation in the banded preconditioner
              setup function.

              @cvodes <node5#sss:cvbandpre> CVBandPrecGetNumRhsEvals
              @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem *)
          val get_num_rhs_evals : serial_bsession -> int
        end
      end

    (** {4:adjbwdout Output} *)

    (** Returns the real and integer workspace sizes.

        @cvodes <node5#sss:optout_main> CVodeGetWorkSpace
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @return ([real_size], [integer_size]) *)
    val get_work_space          : ('a, 'k) bsession -> int * int

    (** Returns the cumulative number of internal steps taken by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetNumSteps
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_num_steps           : ('a, 'k) bsession -> int

    (** Returns the number of calls to the user's right-hand side function.

        @cvodes <node5#sss:optout_main> CVodeGetNumRhsEvals
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_num_rhs_evals       : ('a, 'k) bsession -> int

    (** Returns the number of calls made to the linear solver's setup function.

        @cvodes <node5#sss:optout_main> CVodeGetNumLinSolvSetups
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_num_lin_solv_setups : ('a, 'k) bsession -> int

    (** Returns the number of local error test failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumErrTestFails
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_num_err_test_fails  : ('a, 'k) bsession -> int

    (** Returns the integration method order used during the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastOrder
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_last_order          : ('a, 'k) bsession -> int

    (** Returns the integration method order to be used on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentOrder
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_current_order       : ('a, 'k) bsession -> int

    (** Returns the integration step size taken on the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_last_step           : ('a, 'k) bsession -> float

    (** Returns the integration step size to be attempted on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_current_step        : ('a, 'k) bsession -> float

    (** Returns the the value of the integration step size used on the first
        step.

        @cvodes <node5#sss:optout_main> CVodeGetActualInitStep
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_actual_init_step    : ('a, 'k) bsession -> float

    (** Returns the the current internal time reached by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentTime
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_current_time        : ('a, 'k) bsession -> float

    (** Returns the number of order reductions dictated by the BDF stability
        limit detection algorithm.

        @cvodes <node5#sss:optout_main> CVodeGetNumStabLimOrderReds
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @cvodes <node3#s:bdf_stab> BDF stability limit detection
     *)
    val get_num_stab_lim_order_reds : ('a, 'k) bsession -> int

    (** Returns a suggested factor by which the user's tolerances should be
        scaled when too much accuracy has been requested for some internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetTolScaleFactor
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_tol_scale_factor : ('a, 'k) bsession -> float

    (** Returns the solution error weights at the current time.

        @cvodes <node5#sss:optout_main> CVodeGetErrWeights
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
        @cvodes <node3#ss:ivp_sol> IVP solution (W_i) *)
    val get_err_weights : ('a, 'k) bsession -> ('a, 'k) nvector -> unit

    (** Returns the vector of estimated local errors.

        @cvodes <node5#sss:optout_main> CVodeGetEstLocalErrors
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_est_local_errors : ('a, 'k) bsession -> ('a, 'k) nvector -> unit

    (** Returns the integrator statistics as a group.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_integrator_stats    : ('a, 'k) bsession -> Cvode.integrator_stats

    (** Convenience function that calls get_integrator_stats and prints the
        results to stdout.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val print_integrator_stats  : ('a, 'k) bsession -> unit

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvIters
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_num_nonlin_solv_iters : ('a, 'k) bsession -> int

    (** Returns the number of nonlinear convergence failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_num_nonlin_solv_conv_fails : ('a, 'k) bsession -> int

    (** [nniters, nncfails = get_nonlin_solv_stats s] returns both the
        numbers of nonlinear iterations performed [nniters] and of
        nonlinear convergence failures that have occurred [nncfails].

        @cvode <node5#sss:optout_main> CVodeGetNonlinSolvStats
        @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
     *)
    val get_nonlin_solv_stats : ('a, 'k) bsession -> int *int

    (** {2:adjquad Quadrature Equations} *)

    (** Support for integration of backward quadrature equations that may or may
        not depend on forward sensitivities. *)

    module Quadrature :
      sig
        (** {3:adjquadinit Initialization} *)

        (** These functions compute the quadrature equation right-hand side for
            the backward problem. *)
        type 'a bquadrhsfn =
            NoSens of 'a bquadrhsfn_no_sens
            (** Doesn't depend on forward sensitivities.  See
                {!bquadrhsfn_no_sens} for details. *)
          | WithSens of 'a bquadrhsfn_with_sens
            (** Depends on forward sensitivities.  See
                {!bquadrhsfn_with_sens} for details. *)

        (** Quadrature rhs that doesn't depend on forward
            sensitivities.

            See also {!bquadrhsfn}.

            @cvodes <node7#ss:ODErhs_quad_b> CVQuadRhsFnB *)
        and 'a bquadrhsfn_no_sens = float -> 'a -> 'a -> 'a -> unit

        (** Quadrature rhs that depends on forward sensitivities.

            See also {!bquadrhsfn}.

            @cvodes <node7#ss:ODErhs_quad_sens_B> CVQuadRhsFnBS
          *)
        and 'a bquadrhsfn_with_sens =
          float        (* t *)
          -> 'a        (* y *)
          -> 'a array  (* ys *)
          -> 'a        (* yb *)
          -> 'a        (* qbdot *)
          -> unit

        (** This function, [init s fQB yQB0], activates integration of
            quadrature equations, with or without sensitivities, where [fQB]
            computes the right-hand side of the backward quadrature equations,
            and [yQB0] contains the values of the quadrature variables at [tB0].

            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitB
            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitBS *)
        val init : ('a, 'k) bsession -> 'a bquadrhsfn
                 -> ('a, 'k) nvector -> unit

        (** This function reinitializes the integration of quadrature equations
            during the backward phase.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadReInitB *)
        val reinit : ('a, 'k) bsession -> ('a, 'k) nvector -> unit

        (** {3:adjextraction Extraction function} *)

        (**
          [tret = get s w yqs] fills [yqs] with the quadrature solution vector
          after a successful return from {!backward_normal} or
          {!backward_one_step}, and returns the time reached by the solver.

          @cvodes <node7#sss:quad_get_b> CVodeGetQuadB
         *)
        val get : ('a, 'k) bsession -> ('a, 'k) nvector -> float

        (** {3:adjquadoptin Optional Input Functions} *)

        type ('a, 'k) tolerance =
            NoStepSizeControl
            (** Do not use quadrature variables for step-size control
                (default). *)
          | SStolerances of float * float
            (** [(rel, abs)] : scalar relative and absolute tolerances. *)
          | SVtolerances of float * ('a, 'k) nvector
            (** [(rel, abs)] : scalar relative and vector absolute
                tolerances. *)

        (** Specify whether and how quadrature variables should be used in the
            step size control mechanism.

            @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
        val set_tolerances : ('a, 'k) bsession -> ('a, 'k) tolerance -> unit

        (** {3:adjquadoptout Optional Output Functions} *)

        (** Returns the number of calls to the user's quadrature right-hand side
            function.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
         *)
        val get_num_rhs_evals       : ('a, 'k) bsession -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
         *)
        val get_num_err_test_fails  : ('a, 'k) bsession -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
         *)
        val get_err_weights : ('a, 'k) bsession -> ('a, 'k) nvector -> unit

        (** [nfqevals, nqetfails = get_stats s] returns
            - [fqevals], the number of calls to the user's quadrature function, and,
            - [nqetfails], the number of error test failures due to quadrature variables.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats
            @cvodes <node7#ss:optional_output_b> CVodeGetAdjCVodeBmem
         *)
        val get_stats : ('a, 'k) bsession -> int * int
      end
  end

