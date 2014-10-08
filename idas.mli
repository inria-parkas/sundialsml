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
(*               User Documentation for IDAS v1.1.0                    *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Sensitivity analysis (forward and adjoint) and quadrature equations for
    IDA.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

open Ida_impl
open Sundials

type ('data, 'kind) session = ('data, 'kind) Ida.session
type ('data, 'kind) nvector = ('data, 'kind) Nvector.t

(** {2:quad Quadrature Equations} *)

module Quadrature :
  sig
    (** A skeleton of an enhanced main program:
        + {b Initialize a session [s] per the skeleton at
           {!Ida.session}}
        {[...]}
        The vector of initial values should not include values for the
        quadrature variables.
        + {b Set vector of quadrature variables }
        {[let yQ = Ida.RealArray.of_array [| 0.0; 0.0 |] ]}
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

        @idas <node5#SECTION00570000000000000000> Integration of pure quadrature equations
     *)

    (** {3:quadexcept Exceptions} *)

    (** Quadrature integration was not initialized.

        @idas <node5#ss:quad_get> IDA_NO_QUAD *)
    exception QuadNotInitialized

    (** The quadrature right-hand side function failed in an unrecoverable
        manner.

        @idas <node5#SECTION00572000000000000000> IDA_QRHS_FAIL *)
    exception QuadRhsFuncFailure

    (** The quadrature right-hand side function failed in an
        unrecoverable manner on the first call.

        @idas <node5#SECTION00572000000000000000> IDA_FIRST_QRHS_ERR *)
    exception FirstQuadRhsFuncFailure

    (** Convergence test failures occurred too many times due to
        repeated recoverable errors in the quadrature right-hand side
        function. This value will also be returned if the quadrature
        right-hand side function had repeated recoverable errors
        during the estimation of an initial step size (assuming the
        quadrature variables are included in the error tests).

        @idas <node5#SECTION00572000000000000000> IDA_REP_QRHS_ERR *)
    exception RepeatedQuadRhsFuncFailure

    (** {3:quadinit Initialization} *)

    (** This function, [fQ t y y' rhsQ], computes the quadrature
        equation right-hand side [rhsQ] for a given value of the
        independent variable [t] and state vectors [y] and [y'].

        @idas <node5#ss:user_fct_quad> IDAQuadRhsFn *)
    type 'a quadrhsfn = float -> 'a -> 'a -> 'a -> unit

    (** Activates the integration of quadrature equations.

        @idas <node5#ss:quad_init> IDAQuadInit *)
    val init : ('a, 'b) session -> 'a quadrhsfn -> ('a, 'b) nvector -> unit

    (** Reinitialize the integration of quadrature equations.

        @idas <node5#ss:quad_init> IDAQuadReInit *)
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

        @idas <node5#ss:quad_optional_input> IDASetQuadErrCon
        @idas <node5#ss:quad_optional_input> IDAQuadSStolerances
        @idas <node5#ss:quad_optional_input> IDAQuadSVtolerances *)
    val set_tolerances : ('a, 'b) session -> ('a, 'b) tolerance -> unit

    (** {3:quadout Output Functions} *)

    (** [tret = get s yq] fills [yq] with the quadrature solution
        vector after a successful return from {!Ida.solve_normal} or
        {!Ida.solve_one_step}, and returns the time reached by the
        solver.

        @idas <node5#ss:quad_get> IdaGetQuad *)
    val get : ('a, 'k) session -> ('a, 'k) nvector -> float

    (** [tret = get_dky s t k dkyq] fills [dkyq] with the derivatives
        of the quadrature solution vector after a successful return
        from {!Ida.solve_normal} or {!Ida.solve_one_step}. The time
        requested, [t], must fall within the interval defined by the
        last successful step ({!Ida.get_last_step}). The requested
        order, [k], must be less than or equal to the value returned
        by {!Ida.get_last_order}.

        @idas <node5#ss:quad_get> IdaGetQuadDky
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit

    (** {3:quadoptout Optional Output Functions} *)

    (** Returns the number of calls to the user's quadrature
        right-hand side function.

        @idas <node5#ss:quad_optional_output> IDAGetQuadNumRhsEvals *)
    val get_num_rhs_evals : ('a, 'k) session -> int

    (** Returns the number of local error test failures due to quadrature
        variables.

        @idas <node5#ss:quad_optional_output> IDAGetQuadNumErrTestFails *)
    val get_num_err_test_fails : ('a, 'k) session -> int

    (** Returns the quadrature error weights at the current time.

        @idas <node5#ss:quad_optional_output> IDAGetQuadErrWeights *)
    val get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit

    (** [nfqevals, nqetfails = get_stats s] returns
        - [fqevals], the number of calls to the user's quadrature function, and,
        - [nqetfails], the number of error test failures due to quadrature
          variables.

        @idas <node5#ss:quad_optional_output> IDAGetQuadStats *)
    val get_stats : ('a, 'k) session -> int * int

  end

(** {2:sens (Forward) Sensitivity Analysis} *)

module Sensitivity :
  sig
    (** A skeleton of an enhanced main program:
        + {b Initialize a session [s] per the skeleton at
           {!Ida.session} or {!Idas.Quadrature.init}}
        {[...]}
          Initial value correction, if needed, should be held off
          until sensitivity calculations are activated below.
        + {b Define the sensitivity problem}
        {[let p = Ida.RealArray.make np in
let sp = { pvals = Some p; pbar = ...; plist = ... }]}
        + {b Set sensitivity initial conditions }
        {[let yS0 = Array.init ns (fun _ -> RealArray.init neq 0.0)
let yS'0 = Array.init ns (fun _ -> RealArray.init neq 0.0)]}
        + {b Activate sensitivity calculations}
        {[init s (SStolerances ...) Simultaneous sp fS yS0 y'S0]}
        + {b Correct initial values (optional)}
        {[calc_ic_y s tout1]}
        + {b Set optional inputs}
        {[set_dq_method s ...]}
        + {b Advance the solution in time as per normal}
        + {b Extract sensitivity solution}
        {[let t = get s yS]}

        @idas <node6#SECTION00610000000000000000> Enhanced skeleton for sensitivity analysis *)

    (** {3:sensexcept Exceptions} *)

    (** Forward sensitivity analysis was not initialized.

        @idas <node6#ss:sensi_get> IDA_NO_SENS *)
    exception SensNotInitialized

    (** The sensitivity residual function failed in an unrecoverable manner.

        @idas <node6#SECTION00624000000000000000> IDA_SRES_FAIL *)
    exception SensResFuncFailure

    (** The user's sensitivity residual function repeatedly returned a
        recoverable error flag, but the solver was unable to recover.

        @idas <node6#SECTION00624000000000000000> IDA_REP_SRES_ERR *)
    exception RepeatedSensResFuncFailure

    (** The sensitivity identifier is not valid.  This happens, for
        example, if you have [3] sensitivity variables and request the
        value of sensitivity variable number [5] in {!get_dky1}.

        @idas <node6#SECTION00625000000000000000> IDA_BAD_IS *)
    exception BadSensIdentifier

    (** {3:sensinit Initialization} *)

    (** This function, [fS t y y' res yS y'S resS tmp1 tmp2 tmp3],
        computes the sensitivity residual for all sensitivity equations,
        given
        - [t], the current value of the independent variable,
        - [y], the current value of the state vector,
        - [y'], the current value of the right-hand side of the
                state equations,
        - [yS], the current values of the sensitivity vectors,
        - [y'S], the current values of the derivatives of the
                 sensitivity vectors,
        - [resS], the sensitivity residual vectors.
        - [tmp1], [tmp2], [tmp3] can be used as temporary storage.

        For each [i], this function must compute the vector
        {i (dF/dy)s_i(t)+(dF/dy')s'_i(t)+(dF/dp_i)} and store it in
        {[resS[i]]}.

        If a function is not given during {!init}, then the default
        internal difference-quotient implementation is used.

        @idas <node6#s:user_fct_fwd> IDASensResFn  *)
    type 'a sensresfn =
      float                            (* t *)
      -> 'a                            (* y *)
      -> 'a                            (* y' *)
      -> 'a                            (* resval *)
      -> 'a array                      (* yS *)
      -> 'a array                      (* y'S *)
      -> 'a array                      (* resvalS *)
      -> 'a                            (* tmp1 *)
      -> 'a                            (* tmp2 *)
      -> 'a                            (* tmp3 *)
      -> unit

    (** Specifies a sensitivity solution method.

        @idas <node6#ss:sensi_init> IDASensInit
      *)
    type sens_method =
        Simultaneous
        (** The state and sensitivity variables are corrected at the
            same time.

            @idas <node6#ss:sensi_init> IDA_SIMULTANEOUS
        *)
      | Staggered
        (** The correction step for the sensitivity variables takes
            place at the same time for all sensitivity equations, but
            only after the correction of the state variables has
            converged and the state variables have passed the local
            error test.

            @idas <node6#ss:sensi_init> IDA_STAGGERED
        *)

    (** Used for specifying problem parameter information for
        sensitivity calculations.

        @idas <node6#ss:sens_optional_input> IDASetSensParams *)
    type sens_params = {
      pvals : Sundials.RealArray.t option;
        (** The parameters used to evaluate {i F(t, y, y', p)}. *)
      pbar : Sundials.RealArray.t option;
        (** An array of {i ns} positive scaling factors. *)
      plist : int array option;
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

    (** This function, [init s tol ism ps fS yS0 y'S0], activates the
        forward sensitivity computation, where [ism] selects the
        sensitivity solution method, [ps] gives problem parameter
        information, [fS] computes the sensitivity residuals, [yS0]
        gives the initial values of the sensitivities, and [y'S0]
        gives the initial values of the derivatives of sensitivities.

        Note that any array specified by [ps.pvals] is used to pass
        parameter information during problem solution, that is, the
        library will normally write to it from time to time.

        @idas <node6#ss:sensi_init> IDASensInit *)
    val init :
      ('a, 'b) session ->
      ('a, 'b) tolerance ->
      sens_method ->
      sens_params ->
      'a sensresfn option ->
      ('a, 'b) nvector array -> ('a, 'b) nvector array -> unit

    (** This function reinitializes the forward sensitivity computation.

        @idas <node6#ss:sensi_init> IDASensReInit *)
    val reinit :
      ('a, 'b) session ->
      sens_method -> ('a, 'b) nvector array -> ('a, 'b) nvector array -> unit

    (** Like {!Ida.calc_ic_ya_yd'}, but has extra output parameters
        [ys] and [y's] for receiving corrected sensitivities.

        @ida <node5#ss:idacalcic> IDACalcIC
        @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
        @idas <node6#sss:sens_optout_iccalc> IDAGetSensConsistentIC
      *)
    val calc_ic_ya_yd' :
      ('a, 'b) session ->
      ?y:('a, 'b) nvector ->
      ?y':('a, 'b) nvector ->
      ?ys:('a, 'b) nvector array ->
      ?y's:('a, 'b) nvector array -> ('a, 'b) nvector -> float -> unit

    (** Like {!Ida.calc_ic_y}, but has an extra output parameter [ys]
        for receiving corrected sensitivities.
        @ida <node5#ss:idacalcic> IDACalcIC
        @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
        @idas <node6#sss:sens_optout_iccalc> IDAGetSensConsistentIC
     *)
    val calc_ic_y :
      ('a, 'b) session ->
      ?y:('a, 'b) nvector -> ?ys:('a, 'b) nvector array -> float -> unit

    (** Deactivates forward sensitivity calculations without deallocating
        memory. Sensitivities can be reactivated with {!reinit}.

        @idas <node6#ss:sensi_init> IDASensToggleOff *)
    val toggle_off : ('a, 'k) session -> unit

    (** {3:sensout Output Functions} *)

    (** [tret = get s ys] fills [ys] with the sensitivity solution
        vectors after a successful return from {!Ida.solve_normal} or
        {!Ida.solve_one_step}, and returns the time reached by the
        solver.

        @idas <node6#ss:sensi_get> IDAGetSens
      *)
    val get : ('a, 'b) session -> ('a, 'b) nvector array -> float

    (** [tret = get_dky s t k dkys] fills [dkys] with the
        derivatives of the sensitivity solution vectors after a
        successful return from {!Ida.solve_normal} or
        {!Ida.solve_one_step}. The time requested, [t], must fall
        within the interval defined by the last successful step
        ({!Ida.get_last_step}). The requested order, [k], must be less
        than or equal to the value returned by {!Ida.get_last_order}.

        @idas <node6#ss:sensi_get> IDAGetSensDky
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range.
    *)
    val get_dky :
      ('a, 'b) session -> float -> int -> ('a, 'b) nvector array -> unit

    (** [tret = get1 s i ys] fills [ys] with the [i]th sensitivity
        solution vector after a successful return from
        {!Ida.solve_normal} or {!Ida.solve_one_step}, and returns
        the time reached by the solver.

        @idas <node6#ss:sensi_get> IDAGetSens1
        @raise BadIS The index [i] is not in the allowed range. *)
    val get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float

    (** [tret = get_dky1 s t k i dkys] fills [dkys] with the
        derivatives of the [i]th sensitivity solution vector after a
        successful return from {!Ida.solve_normal} or
        {!Ida.solve_one_step}. The time requested, [t], must fall
        within the interval defined by the last successful step
        ({!Ida.get_last_step}). The requested order, [k], must be
        less than or equal to the value returned by
        {!Ida.get_last_order}.

        @idas <node6#ss:sensi_get> IDAGetSensDky1
        @raise BadIS The index [i] is not in the allowed range.
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range.  *)
    val get_dky1 :
      ('a, 'k) session -> float -> int -> int -> ('a, 'k) nvector -> unit

    (** {3:sensoptin Optional Input Functions} *)

    (** Specify the integration tolerances for sensitivities.

        {b NB}: Unlike the other [set_tolerances] functions in [Idas],
        this one does {b not} call {!set_err_con} (which defaults to
        [false]).

        @idas <node6#ss:quad_sens_optional_input> IDASensSStolerances
        @idas <node6#ss:quad_sens_optional_input> IDASensSVtolerances
        @idas <node6#ss:quad_sens_optional_input> IDASensEEtolerances *)
    val set_tolerances : ('a, 'b) session -> ('a, 'b) tolerance -> unit

    (** Set whether sensitivity variables should be used in the error
        control mechanism (the default is [false]).

        @idas <node6#ss:sens_optional_input> IDASetSensErrCon *)
    val set_err_con : ('a, 'k) session -> bool -> unit

    (** Difference quotient strategy accepted by {!set_dq_method} below. *)
    type dq_method = DQCentered (** IDA_CENTERED *)
                   | DQForward  (** IDA_FORWARD *)

    (** [set_dq_method s dqtype dqrhomax] specifies the difference quotient
        strategy in the case in which the right-hand side of the sensitivity
        equations is to be computed by IDAS; [dqrhomax] is used in deciding
        the switching between simultaneous or separate approximations of the two
        terms in the sensitivity residual.

        @idas <node6#ss:sens_optional_input> IDASetSensDQMethod
        @idas <node3#ss:fwd_sensi> Forward Sensitivity Analysis *)
    val set_dq_method : ('a, 'k) session -> dq_method -> float -> unit

    (** Specifies the maximum number of nonlinear solver iterations for
        sensitivity variables permitted per step.

        @idas <node6#ss:sens_optional_input> IDASetSensMaxNonlinIters *)
    val set_max_nonlin_iters : ('a, 'k) session -> int -> unit

    (** {3:sensoptout Optional Output Functions} *)

    (** Returns the number of calls to the sensitivity residual function.

        @idas <node6#ss:sens_optional_output> IDAGetSensNumResEvals *)
    val get_num_res_evals : ('a, 'k) session -> int

    (** Returns the number of calls to the user's residual function due
        to the internal finite difference approximation of the sensitivity
        residual.

        @idas <node6#ss:sens_optional_output> IDAGetNumResEvalsSens *)
    val get_num_res_evals_sens : ('a, 'k) session -> int

    (** Returns the number of local error test failures for the sensitivity
        variables that have occurred.

        @idas <node6#ss:sens_optional_output> IDAGetSensNumErrTestFails *)
    val get_num_err_test_fails : ('a, 'k) session -> int

    (** Returns the number of calls made to the linear solver's setup function
        due to forward sensitivity calculations.

        @idas <node6#ss:sens_optional_output> IDAGetSensNumLinSolvSetups *)
    val get_num_lin_solv_setups : ('a, 'k) session -> int

    (** Return type of {!get_stats}. *)
    type sensitivity_stats = {
      num_res_evals : int;
      num_sens_evals : int;
      num_err_test_fails : int;
      num_lin_solv_setups : int;
    }

    (** Returns all of the sensitivity-related solver statistics as a group.

        @idas <node6#ss:sens_optional_output> IDAGetSensStats *)
    val get_stats : ('a, 'k) session -> sensitivity_stats

    (** Returns the sensitivity error weight vectors at the current time.

        @idas <node6#ss:sens_optional_output> IDAGetSensErrWeights
        @idas <node3#e:errwt> Eq. (2.7) IVP solution (W_i) *)
    val get_err_weights : ('a, 'k) session -> ('a, 'k) nvector array -> unit

    (** Returns the number of nonlinear iterations performed for sensitivity
        calculations.

        @idas <node6#ss:sens_optional_output> IDAGetSensNumNonlinSolvIters *)
    val get_num_nonlin_solv_iters : ('a, 'k) session -> int


    (** Returns the number of nonlinear convergence failures that have occurred
        for sensitivity calculations.

        @idas <node6#ss:sens_optional_output> IDAGetSensNumNonlinSolvConvFails *)
    val get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int

    (** [nni, ncfn = get_nonlin_solv_stats s] returns the sensitivity-related
        nonlinear solver statistics as a group, where [nni] is the number of
        nonlinear iterations performed for sensitivity calculations, and [ncfn]
        is the number of nonlinear convergence failures that have occurred for
        sensitivity calculations.

        @idas <node6#ss:sens_optional_output> IDAGetSensNonlinSolvStats
     *)
    val get_nonlin_solv_stats : ('a, 'k) session -> int * int

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

            @idas <node6#SECTION00640000000000000000> Integration of quadrature equations depending on forward sensitivities *)

        (** {3:quadsensexcept Exceptions} *)

        (** Quadrature integration was not initialized.

            @idas <node6#ss:quad_sens_get> IDA_NO_QUADSENS *)
        exception QuadSensNotInitialized

        (** The sensitivity quadrature right-hand side function failed
            in an unrecoverable manner.

            @idas <node5#SECTION00642000000000000000> IDA_QSRHS_FAIL
          *)
        exception QuadSensRhsFuncFailure

        (** The user-provided sensitivity-dependent quadrature right-
            hand side function failed in an unrecoverable manner on
            the first call.

            @idas <node5#SECTION00642000000000000000> IDA_FIRST_QSRHS_ERR *)
        exception FirstQuadSensRhsFuncFailure

        (** The user-provided sensitivity-dependent quadrature right-
            hand side repeatedly returned a recoverable error flag,
            but the solver was unable to recover.

            @idas <node6#SECTION00642000000000000000> IDA_REP_QSRHS_ERR
          *)
        exception RepeatedQuadSensRhsFuncFailure

        (** {3:quadsensinit Initialization} *)

        (** This function, [fQS t y yS yQdot rhsvalQs tmp1 tmp2], computes the
            sensitivity quadrature equation right-hand side given
            - [t], the current value of the independent variable,
            - [y], the current value of the state vector,
            - [y'], the current value of the derivative of the state vector,
            - [yS], an array of dependent sensitivity vectors,
            - [y'S], an array of derivatives of dependent sensitivity vectors,
            - [rrQ], the current value of the quadrature residual,
            - [rhsvalQs], the right-hand side vectors must be stored here,
            - [tmp1], [tmp2], and [tmp3] can be used as temporary storage.

           @idas <node6#ss:user_fct_quad_sens> IDAQuadSensRhsFn *)
        type 'a quadsensrhsfn =
          float          (* t *)
          -> 'a          (* y *)
          -> 'a          (* y' *)
          -> 'a array    (* yS *)
          -> 'a array    (* y'S *)
          -> 'a          (* rrQ *)
          -> 'a array    (* rhsvalQs *)
          -> 'a          (* tmp1 *)
          -> 'a          (* tmp2 *)
          -> 'a          (* tmp3 *)
          -> unit

        (** This function, [init s fQS yQS0], activates integration of
            quadrature equations depending on sensitivities, where
            [fQS] computes the right-hand side of the
            sensitivity-dependent quadrature equations, and [yQS0]
            contains the initial values of sensitivity-dependent
            quadratures.  When no [fQB] is supplied, the solver uses
            an internal implementation based on difference quotients.

            @idas <node6#ss:quad_sens_init> IDAQuadSensInit *)
        val init :
          ('a, 'b) session -> ?fQS:'a quadsensrhsfn
          -> ('a, 'b) nvector array -> unit

        (** This function reinitializes the forward sensitivity computation.

            @idas <node6#ss:quad_sens_init> IDAQuadSensReInit *)
        val reinit : ('a, 'b) session -> ('a, 'b) nvector array -> unit


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

            @idas <node6#ss:quad_sens_optional_input> IDASetQuadSensErrCon
            @idas <node6#ss:quad_sens_optional_input> IDAQuadSensSStolerances
            @idas <node6#ss:quad_sens_optional_input> IDAQuadSensSVtolerances
            @idas <node6#ss:quad_sens_optional_input> IDAQuadSensEEtolerances *)
        val set_tolerances : ('a, 'k) session -> ('a, 'k) tolerance -> unit

        (** {3:quadsensout Output Functions} *)

        (** [tret = get s yqs] fills [yqs] with quadrature solution vectors
            after a successful return from {!Ida.solve_normal} or
            {!Ida.solve_one_step}, and returns the time reached by
            the solver.

            @idas <node6#ss:quad_sens_get> IDAGetQuadSens *)
        val get : ('a, 'b) session -> ('a, 'b) nvector array -> float

        (** [tret = get s i yqs] fills [yqs] with the [i]th quadrature
            solution vector after a successful return from
            {!Ida.solve_normal} or {!Ida.solve_one_step}, and
            returns the time reached by the solver.

            @idas <node6#ss:quad_sens_get> IDAGetQuadSens1
            @raise BadIS The index [i] is not in the allowed range.
          *)
        val get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float

        (** [tret = get_dky s t k dkyqs] fills [dkyqs] with the
            derivatives of the quadrature solution vectors after a
            successful return from {!Ida.solve_normal} or
            {!Ida.solve_one_step}. The time requested, [t], must
            fall within the interval defined by the last successful
            step ({!Ida.get_last_step}). The requested order, [k],
            must be less than or equal to the value returned by
            {!Ida.get_last_order}.

            @idas <node6#ss:quad_sens_get> IDAGetQuadSensDky
            @raise BadIS The index is not in the allowed range.
            @raise BadK [k] is not in the range 0, 1, ..., [qlast].
            @raise BadT [t] is not in the allowed range. *)
        val get_dky :
          ('a, 'b) session -> float -> int -> ('a, 'b) nvector array -> unit

        (** [tret = get_dky s t k i dkyqs] fills [dkyqs] with the derivatives of
            the [i]th quadrature solution vector after a successful return from
            {!Ida.solve_normal} or {!Ida.solve_one_step}.
            The time requested, [t], must fall within the interval defined by
            the last successful step ({!Ida.get_last_step}). The
            requested order, [k], must be less than or equal to the value
            returned by {!Ida.get_last_order}.

            @idas <node6#ss:quad_sens_get> IDAGetQuadSensDky1
            @raise BadK [k] is not in the range 0, 1, ..., [qlast].
            @raise BadT [t] is not in the allowed range. *)
        val get_dky1 :
          ('a, 'k) session -> float -> int -> int -> ('a, 'k) nvector -> unit

        (** Returns the number of calls to the user's quadrature right-hand side
            function.

            @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensNumRhsEvals *)
        val get_num_rhs_evals : ('a, 'k) session -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensNumErrTestFails *)
        val get_num_err_test_fails : ('a, 'k) session -> int

        (** Returns the quadrature error weights at the current time.

            @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensErrWeights *)
        val get_err_weights :
          ('a, 'b) session -> ('a, 'b) nvector array -> unit

        (** [nfqevals, nqetfails = get_stats s] returns
            - [fqevals], the number of calls to the user's quadrature function,
            and,
            - [nqetfails], the number of error test failures due to quadrature
            variables.

          @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensStats *)
        val get_stats : ('a, 'k) session -> int * int
      end
  end


(** {2:adj Adjoint Sensitivity Analysis} *)

module Adjoint :
  sig
    (** A skeleton of an enhanced main program:
        + {b Initialize a session [s] per the skeleton at {!Ida.init}}
          {[...]}
          Adding quadrature variables using {!Quadrature.init} if desired.
          Initial value correction, if needed, should be held off
          until sensitivity calculations are activated below.
        + {b Initialize the adjoint computation}
          {[init s nsteps IHermite]}
        + {b Integrate forward problem}
          {[let t, ncheck, r = forward_normal s tout y0]}
        + {b Setup the backward problem and attach a linear solver}
          {[let yB0  = RealArray.of_list [0.0; 0.0; ...]
          let yB'0 = RealArray.of_list [0.0; 0.0; ...]
let bs = init_backward s (Spils.spgmr ...) (SStolerances ...) (NoSens fB) tB0 yB0 yB'0]}
        + {b Set optional inputs}
          {[set_max_ord bs ...]}
        + {b Initialize quadrature calculation}
          {[Quadrature.init bs fQb yQB0]}
        + {b Integrate backward problem}
          {[backward_normal s tB]}
        + {b Extract quadrature variables}
          {[let t = Quadrature.get s yQS]}

        @idas <node7#ss:skeleton_adj> Enhanced Skeleton for Adjoint Sensitivity Analysis
    *)

    (** Adjoint sensitivity analysis was not initialized.

        @idas <node7#sss:idasolvef> IDA_NO_ADJ *)
    exception AdjointNotInitialized

    (** Neither {!forward_normal} nor {!forward_one_step} has previously been
        called.

        @idas <node7#sss:idasolveb> IDA_NO_FWD *)
    exception NoForwardCall

    (** Reinitialization of the forward problem failed at the first checkpoint
        (corresponding to the initial time of the forward problem).

        @idas <node7#sss:idasolveb> IDA_REIFWD_FAIL *)
    exception ForwardReinitFailure

    (** An error occured during the integration of the forward problem.

        @idas <node7#sss:idasolveb> IDA_FWD_FAIL *)
    exception ForwardFailure

    (** No backward problem has been created.

        @idas <node7#sss:idasolveb> IDA_NO_BCK *)
    exception NoBackwardProblem

    (** The final time [tB0] was outside the interval over which the forward
        problem was solved.

        @idas <node7#sss:idainitb> IDA_BAD_TB0 *)
    exception BadFinalTime

    (** {3:adjfwd Forward Solutions} *)

    (** {4:adjfwdinit Initialization} *)

    (** Specifies the type of interpolation.

        @idas <node3#ss:checkpointing> Checkpointing scheme *)
    type interpolation = IPolynomial    (** IDA_POLYNOMIAL *)
                       | IHermite       (** IDA_HERMITE *)

    (** [init s nd interp] initializes the forward-backward problem with [nd]
        integration steps between consecutive checkpoints and variable-degree
        interpolation according to [interp]. This function must be called before
        either {!forward_normal} or {!forward_one_step}.

        @idas <node7#sss:idaadjinit> IDAAdjInit *)
    val init : ('a, 'k) session -> int -> interpolation -> unit

    (** {4:adjfwdintegration Forward Integration} *)
    (** [tret, ncheck, sr = forward_normal s tout y y'] integrates the
        forward problem over an interval and saves checkpointing
        data. The function takes as arguments the next time at which a
        solution is desired ([tout]), two vectors for storing the
        computed result ([y] and [y']), and returns the time reached
        by the solver ([tret]), the number of checkpoints stored so
        far ([ncheck]), and whether the solver reached [tout] or not
        ([sr]).

        This call asks the solver to take internal steps until it has
        reached or just passed the [tout] parameter
        ([IDA_NORMAL]). The solver then interpolates in order to
        return an approximate value of [y] and [y'] at time [tout].

        @idas <node7#sss:idasolvef> IDASolveF
        @raise Ida.IllInput           One of the inputs is invalid.
        @raise Ida.TooMuchWork        Could not reach [tout] in [mxstep] steps
        @raise Ida.TooMuchAccuracy    Could not satisfy the demanded accuracy
        @raise Ida.ErrFailure         Too many error test failures.
        @raise Ida.ConvergenceFailure Too many convergence test failures.
        @raise Ida.LinearSetupFailure Unrecoverable failure in linear solver setup function.
        @raise Ida.LinearSolveFailure Unrecoverable failure in linear solver solve function.
        @raise AdjointNotInitialized    The [init] function has not previously been called. *)
    val forward_normal :
      ('a, 'k) session ->
      float ->
      ('a, 'k) nvector ->
      ('a, 'k) nvector -> float * int * Ida.solver_result

    (** [tret, ncheck, sr = forward_normal s tout y y'] integrates the
        forward problem over an interval and saves checkpointing
        data. The function takes as arguments the next time at which a
        solution is desired ([tout]), two vectors for storing the
        computed result ([y] and [y']), and returns the time reached
        by the solver ([tret]), the number of checkpoints stored so
        far ([ncheck]), and whether the solver reached [tout] or not
        ([sr]).

        This call asks the solver to take one internal step and to
        return the solution at the point reached by that step
        ([IDA_ONE_STEP]).

        @idas <node7#sss:idasolvef> IDASolveF
        @raise Ida.IllInput           One of the inputs is invalid.
        @raise Ida.TooMuchWork        Could not reach [tout] in [mxstep] steps
        @raise Ida.TooMuchAccuracy    Could not satisfy the demanded accuracy
        @raise Ida.ErrFailure         Too many error test failures.
        @raise Ida.ConvergenceFailure Too many convergence test failures.
        @raise Ida.LinearSetupFailure Unrecoverable failure in linear solver setup function.
        @raise Ida.LinearSolveFailure Unrecoverable failure in linear solver solve function.
        @raise AdjointNotInitialized    The [init] function has not previously been called. *)
    val forward_one_step :
      ('a, 'k) session ->
      float ->
      ('a, 'k) nvector ->
      ('a, 'k) nvector -> float * int * Ida.solver_result

    (** {3:adjbwd Backward Problems} *)

    (** Identifies a backward problem. *)
    type ('a, 'k) bsession = ('a, 'k) AdjointTypes.bsession
    type serial_bsession = (RealArray.t, Nvector_serial.kind) bsession

    (** {4:adjbwdinit Initialization} *)

    (** These functions evaluate the residual function of the backward
        DAE system with or without a dependence on forward
        sensitivities.

        @idas <node7#ss:ODErhs_b> IDAResFnB
        @idas <node7#ss:ODErhs_bs> IDAResFnBS
        @idas <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis
      *)
    type 'a bresfn =
      | NoSens of 'a bresfn_no_sens
        (** Doesn't depend on forward sensitivities.  See
            {!bresfn_no_sens} for details. *)
      | WithSens of 'a bresfn_with_sens
        (** Depends on forward senstivites.  See {!bresfn_with_sens} for
            details. *)

    (** Backward DAE residual function that doesn't depend on forward
        sensitivities.

        @idas <node7#ss:ODErhs_b> IDAResFnB
        @idas <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis *)
    and 'a bresfn_no_sens =
      float             (* t *)
      -> 'a             (* y *)
      -> 'a             (* y' *)
      -> 'a             (* yB *)
      -> 'a             (* y'B *)
      -> 'a             (* resvalB *)
      -> unit

    (** Backward DAE residual function that depends on forward
        sensitivities.

        @idas <node7#ss:ODErhs_bs> IDAResFnBS
        @idas <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis *)
    and 'a bresfn_with_sens =
      float          (* t *)
      -> 'a          (* y *)
      -> 'a          (* y' *)
      -> 'a array    (* yS *)
      -> 'a array    (* y'S *)
      -> 'a          (* yB *)
      -> 'a          (* y'B *)
      -> 'a          (* resvalB *)
      -> unit

    type 'a single_tmp = 'a
    type 'a triple_tmp = 'a * 'a * 'a

    (** Arguments common to all Jacobian callback functions.

        @idas <node7#ss:densejac_b> IDADlsDenseJacFnB
        @idas <node7#ss:bandjac_b> IDADlsBandJacFnB
        @idas <node7#ss:jactimesvec_b> IDASpilsJacTimesVecFnB
        @idas <node7#ss:psolve_b> IDASpilsPrecSolveFnB
        @idas <node7#ss:psetup_b> IDASpilsPrecSetupFnB *)
    type ('t, 'a) jacobian_arg =
      {
        jac_t : float;
        jac_y : 'a;
        jac_y' : 'a;
        jac_yb : 'a;
        jac_y'b : 'a;
        jac_resb : 'a;
        jac_coef : float;
        jac_tmp : 't;
      }

    type bandrange = Ida.bandrange =
      { mupper : int; (** The upper half-bandwidth.  *)
        mlower : int; (** The lower half-bandwidth.  *) }

    (** Specify which variables are algebraic and which variables are
        differential, needed for {!set_suppress_alg}.  This function must
        not be called if you already called {!calc_ic_ya_yd'}.

        The SUNDIALS manual is not clear about whether it's safe to change the
        variable types after you've already set it.

        [set_var_types] corresponds to [IDASetIdB] in the C interface,
        and an alias {!set_id} is also available in this binding.  We
        prefer the more descriptive name {!set_var_types}, however.

        @idas <node7#ss:optional_input_b> IDASetIdB
    *)
    val set_var_types : ('a, 'b) bsession -> ('a, 'b) nvector -> unit

    (** An unpreferred alias for {!set_var_types}.  SUNDIALS calls
        variable types by the cryptic name "Id", and this OCaml
        binding preserves this alternative naming to help users
        transition from the C interface.

        @idas <node7#ss:optional_input_b> IDASetIdB
    *)
    val set_id : ('a, 'b) bsession -> ('a, 'b) nvector -> unit

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

        @idas <node7#ss:optional_input_b> IDASetSuppressAlgB
    *)
    val set_suppress_alg : ('a, 'b) bsession -> bool -> unit

    (** This function provides the same functionality for backward
        problems as {!Ida.calc_ic_ya_yd'} provides for forward
        problems: compute the algebraic components of [yB] and
        differential components of [y'B], given the differential
        components of [yB].  If [bs] is a backward session initialized
        with a call to {!init_backward} at time value [tB0] and a
        residual function independent of forward sensitivities, then
        [calc_ic bs tBout1 yB0 y'B0 yS0 y'S0] corrects the initial
        values [yB0] and [y'B0] at time [tB0] for the backward
        problem, where:
        - [tBout1] is the first value of [t] at which a solution will
          be requested (with {!backward_normal} or
          {!backward_one_step}). This value is needed here only to
          determine the direction of integration and rough scale in
          the independent variable [t].
        - [yB0] is the forward solution at final time [tB0].
        - [y'B0] is the forward derivative solution at final time [tB0].

        The optional arguments [~yb] and [~y'b], if present, will be
        filled with the corrected values of [yB] and [y'B],
        respectively.  The storage for [~yb] and [~y'b] should not
        overlap with each other, but they can overlap with [yB] and/or
        [y'B].

        Calling this function is optional. It is only necessary when
        the initial conditions do not satisfy the adjoint system.

        This function works if the residual function is independent of
        forward sensitivities.  If the residual function depends on
        sensitivities, {!calc_ic_sens} should be used instead.

        @idas <node7#sss:idacalcicB> IDACalcICB
        @idas <node7#ss:optional_output_b> IDAGetConsistentICB
      *)
    val calc_ic :
      ('a, 'b) bsession ->
      ?yb:('a, 'b) nvector ->
      ?y'b:('a, 'b) nvector ->
      float -> ('a, 'b) nvector -> ('a, 'b) nvector -> unit

    (** This function provides the same functionality for backward
        problems as {!Ida.calc_ic_ya_yd'} provides for forward
        problems: compute the algebraic components of [yB] and
        differential components of [y'B], given the differential
        components of [yB].  If [bs] is a backward session initialized
        with a call to {!init_backward} at time value [tB0] and a
        sensitivity-dependent residual function, then [calc_ic bs
        tBout1 yB0 y'B0 yS0 y'S0] corrects the initial values [yB0]
        and [y'B0] at time [tB0] for the backward problem, where:
        - [tBout1] is the first value of [t] at which a solution will
          be requested (with {!backward_normal} or
          {!backward_one_step}). This value is needed here only to
          determine the direction of integration and rough scale in
          the independent variable [t].
        - [yB0] is the forward solution at final time [tB0].
        - [y'B0] is the forward derivative solution at final time [tB0].
        - [yS0] is an array of vectors containing the sensitivities
          of the forward solution at final time [tB0].
        - [y'S0] is an array of vectors containing the sensitivities
          of the forward derivative solution at final time [tB0].

        The optional arguments [~yb] and [~y'b], if present, will be
        filled with the corrected values of [yB] and [y'B],
        respectively.  The storage for [~yb] and [~y'b] should not
        overlap with each other, but they can overlap with [yB] and/or
        [y'B].

        Calling this function is optional. It is only necessary when
        the initial conditions do not satisfy the adjoint system.

        This function works if the residual function is independent of
        forward sensitivities.  If the residual function depends on
        sensitivities, {!calc_ic_sens} should be used instead.

        @idas <node7#sss:idacalcicB> IDACalcICBS
        @idas <node7#ss:optional_output_b> IDAGetConsistentICB
      *)
    val calc_ic_sens :
      ('a, 'b) bsession ->
      ?yb:('a, 'b) nvector ->
      ?y'b:('a, 'b) nvector ->
      float ->
      ('a, 'b) nvector ->
      ('a, 'b) nvector ->
      ('a, 'b) nvector array -> ('a, 'b) nvector array -> unit

    (** Specify a linear solver.

        @idas <node7#sss:lin_solv_b> Linear Solver Initialization Functions *)
    type ('data, 'kind) linear_solver =
      ('data, 'kind) AdjointTypes.linear_solver
    type serial_linear_solver =
      (RealArray.t, Nvector_serial.kind) linear_solver

    type ('a, 'k) tolerance =
        SStolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVtolerances of float * ('a, 'k) nvector
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

    (** [init_backward s linsolv tol fB tB0 yB0 yB'0] adds and initializes a
        backward problem that may or may not depend on forward sensitivities,
        where
        - [s] is the parent session (going forward),
        - [linsolv] specifies a linear solver,
        - [tol]     specifies the tolerances, see {!tolerance},
        - [fB]      computes the right-hand side of the backward ODE problem,
        - [tB0]     specifies the endpoint where final conditions are provided
                    for the backward problem, normally equal to the endpoint of
                    the forward integration, and,
        - [yB0]     is the final value of the backward problem.

        @idas <node7#sss:idainitb> IDACreateB
        @idas <node7#sss:idainitb> IDAInitB
        @idas <node7#sss:idainitb> IDAInitBS
        @idas <node7#sss:idatolerances_b> IDASStolerancesB
        @idas <node7#sss:idatolerances_b> IDASVtolerancesB 
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise BadFinalTime      The final time is outside the interval over which the forward problem was solved. *)
    val init_backward :
         ('a, 'k) session
      -> ('a, 'k) linear_solver
      -> ('a, 'k) tolerance
      -> 'a bresfn
      -> float
      -> ('a, 'k) nvector
      -> ('a, 'k) nvector
      -> ('a, 'k) bsession

    (** Reinitialize the backward problem.

        @idas <node7#sss:idainitb> IDAReInitB
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise BadFinalTime      The final time is outside the interval over which the forward problem was solved. *)
    val reinit :
      ('a, 'b) bsession ->
      ?linsolv:('a, 'b) linear_solver ->
      float -> ('a, 'b) nvector -> ('a, 'b) nvector -> unit

    (** {4:adjbwdintegration Backward Integration} *)

    (** [backward_normal s tbout] integrates the backward DAE
        problem. The function takes internal steps until it has
        reached or just passed the user-specified value [tbout]
        ([IDA_NORMAL]). The solver then interpolates in order to
        return an approximate value of [y(tbout)] when {!get} is
        called.

        @idas <node7#sss:idasolveb> IDASolveB
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise NoBackwardProblem        The [init_backward] function has not previously been called.
        @raise NoForwardCall            Neither [forward_normal] nor [forward_one_step] has previously been called.
        @raise Ida.IllInput           One of the inputs is invalid.
        @raise Ida.TooMuchWork        Could not reach [tout] in [mxstep] steps
        @raise Ida.TooMuchAccuracy    Could not satisfy the demanded accuracy
        @raise Ida.ErrFailure         Too many error test failures.
        @raise Ida.ConvergenceFailure Too many convergence test failures.
        @raise Ida.LinearSetupFailure Unrecoverable failure in linear solver setup function.
        @raise Ida.LinearSolveFailure Unrecoverable failure in linear solver solve function.
        @raise BadOutputTime            The requested output time is outside the interval over which the forward problem was solved.
        @raise ForwardReinitializationFailed Reinitialization of the forward problem failed at the first checkpoint (corresponding to the initial time of the forward problem).
        @raise ForwardFail              An error occurred during the integration of the forward problem. *)
    val backward_normal : ('a, 'k) session -> float -> unit

    (** [backward_one_step s tbout] integrates the backward DAE problem. The
        function takes one internal step ([IDA_ONE_STEP]).

        @idas <node7#sss:idasolveb> IDASolveB
        @raise AdjointNotInitialized    The [init] function has not previously been called.
        @raise NoBackwardProblem        The [init_backward] function has not previously been called.
        @raise NoForwardCall            Neither [forward_normal] nor [forward_one_step] has previously been called.
        @raise Ida.IllInput           One of the inputs is invalid.
        @raise Ida.TooMuchWork        Could not reach [tout] in [mxstep] steps
        @raise Ida.TooMuchAccuracy    Could not satisfy the demanded accuracy
        @raise Ida.ErrFailure         Too many error test failures.
        @raise Ida.ConvergenceFailure Too many convergence test failures.
        @raise Ida.LinearSetupFailure Unrecoverable failure in linear solver setup function.
        @raise Ida.LinearSolveFailure Unrecoverable failure in linear solver solve function.
        @raise BadOutputTime            The requested output time is outside the interval over which the forward problem was solved.
        @raise ForwardReinitializationFailed Reinitialization of the forward problem failed at the first checkpoint (corresponding to the initial time of the forward problem).
        @raise ForwardFail              An error occurred during the integration of the forward problem. *)
    val backward_one_step : ('a, 'k) session -> float -> unit

    (** [tret = get bs yB y'B] returns the solution of the backward DAE problem
        in [yB] and [y'B] at time [tret].

        @idas <node7#sss:idasolveb> IDAGetB *)
    val get :
      ('a, 'b) bsession -> ('a, 'b) nvector -> ('a, 'b) nvector -> float

    (** [tret = get_dky s t k dkys] fills [dkys] with the derivatives of the
        sensitivity solution vectors after a successful return from
        {!backward_normal} or {!backward_one_step}. The time requested, [t],
        must fall within the interval defined by the last successful step
        ({!get_last_step}). The requested order, [k], must be less than or equal
        to the value returned by {!get_last_order}.

        @idas <node5#ss:optional_dky> IDAGetDky
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky :
      ('a, 'b) bsession -> float -> int -> ('a, 'b) Ida.nvector -> unit

    (** {4:adjbwdoptout Optional Output Functions} *)

    (** Instructs {!forward_normal} and {!forward_one_step} not to save
        checkpointing data for forward sensitivities anymore.

        @idas <node7#SECTION00728000000000000000> IDAAdjSetNoSensi *)
    val set_no_sensitivity : ('a, 'k) session -> unit

    (** Specifies the maximum order of the linear multistep method.

        @idas <node7#ss:optional_input_b> IDASetMaxOrdB *)
    val set_max_ord : ('a, 'b) bsession -> int -> unit

    (** Specifies the maximum number of steps to be taken by the solver in its
        attempt to reach the next output time.

        @idas <node7#ss:optional_input_b> IDASetMaxNumStepsB *)
    val set_max_num_steps : ('a, 'b) bsession -> int -> unit

    (** Specifies the initial step size.

        @idas <node7#ss:optional_input_b> IDASetInitStepB *)
    val set_init_step : ('a, 'b) bsession -> float -> unit

    (** Specifies an upper bound on the magnitude of the step size.

        @idas <node7#ss:optional_input_b> IDASetMaxStepB *)
    val set_max_step : ('a, 'b) bsession -> float -> unit

    (** {4:adjbwddirect Direct Linear Solver} *)

    module Dls :
      sig
        (** This module provides the same functionality as {!Ida.Dls}, but
            for backward problems specifically.  Note sundials (as of
            2.5.0) doesn't support LAPACK for backward problems.  *)

        (** This function computes the dense Jacobian of the backward problem
            (or an approximation to it).

            @idas <node7#ss:densejac_b> IDADlsDenseJacFnB *)
        type dense_jac_fn =
            (RealArray.t triple_tmp, RealArray.t) jacobian_arg ->
            Dls.DenseMatrix.t -> unit

        (** This function computes the banded Jacobian of the backward problem
            (or an approximation to it).

            @idas <node7#ss:bandjac_b> IDADlsBandJacFnB *)
        type band_jac_fn =
            bandrange ->
            (RealArray.t triple_tmp, RealArray.t) jacobian_arg ->
            Dls.BandMatrix.t -> unit

        (** Direct linear solver with dense matrix.  The optional argument
            specifies a callback function that computes an approximation to the
            Jacobian matrix (see {!dense_jac_fn} for details).  If this
            argument is [None], then IDA uses a default implementation based
            on difference quotients.  See also {!Dls}.

            @idas <node7#sss:lin_solv_b> IDADenseB
            @idas <node7#SECTION00729200000000000000> IDADlsSetDenseJacFnB
            @idas <node7#ss:densejac_b> IDADlsDenseJacFnB *)
        val dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

        (** Direct linear solver with banded matrix.  The arguments specify the
            width of the band ({!bandrange}) and an optional Jacobian
            function ({!bband_jac_fn}).  If the Jacobian function is [None],
            IDAS uses an internal implementation based on difference
            quotients. See also {!Dls}.

            @idas <node7#sss:lin_solv_b> IDABandB
            @idas <node7#SECTION00729300000000000000> IDADlsSetBandJacFnB
            @idas <node7#ss:bandjac_b> IDADlsBandJacFnB *)
        val band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver
      end

    (** {4:adjbwdspils Scaled Preconditioned Iterative Linear Solvers (SPILS)} *)

    module Spils :
      sig
        (** Scaled Preconditioned Iterative Linear Solvers (SPILS)

            @idas <node7#ss:optional_output_b> Optional output functions for the backward problem. *)

        type gramschmidt_type =
          Spils.gramschmidt_type =
          | ModifiedGS
          | ClassicalGS

        (** Called like [prec_solve_fn arg r z delta] to solve the
            linear system {i P}[z] = [r], where {i P} is the (left)
            preconditioner matrix.
            - [arg] supplies the basic problem data as a {!jacobian_arg}.
            - [r] is the right-hand side vector.
            - [z] is the vector in which the result must be stored.
            - [delta] is an input tolerance.

            If set to [None] then no preconditioning is performed, and
            [prec_setup_fn] and [jac_times_vec_fn] are ignored.

            {i P} should approximate, at least crudely, the system
            Jacobian matrix {i J = dF/dy + c * dF/dy'} where {i F} is
            the residual function and {i c} is [arg.jac_coef].

            [delta] is an input tolerance to be used if an iterative
            method is employed in the solution.  In that case, the
            residual vector res = [r]
            - {i P} [z] of the system should be made less than [delta] in weighted
            l2 norm, i.e. [sqrt (sum over i ((res.{i} * ewt.{i})^2)) <
            delta], where the vector ewt can be obtained through
            {!get_err_weights}.

            This function can raise {!Sundials.RecoverableFailure} to
            instruct the integrator to retry with a different step
            size.  Raising any other kind of exception aborts the
            integrator.

            {b NB:} [r], [z], and the elements of [arg] must no longer
            be accessed after [prec_solve_fn] has returned, i.e. if
            their values are needed outside of the function call, then
            they must be copied to separate physical structures.

            @idas <node7#ss:psolve_b> IDASpilsPrecSolveFnB
          *)
        and 'a prec_solve_fn =
          ('a single_tmp, 'a) jacobian_arg
          -> 'a
          -> 'a
          -> float
          -> unit

        (** A function that preprocesses and/or evaluates any
            Jacobian-related data needed by [prec_solve_fn] above.

            The sole argument to this function specifies the basic
            problem data as a {!jacobian_arg}.

            Note that unlike in CVODES, whatever data this function
            computes has to be recomputed every time it is called.

            This function can raise {!Sundials.RecoverableFailure} to
            instruct the integrator to retry with a different step
            size.  Raising any other kind of exception aborts the
            integrator.

            {b NB:} The elements of [jac] must no longer be accessed
            after [psetup] has returned a result, i.e. if their values
            are needed outside of the function call, then they must be
            copied to a separate physical structure.

            @idas <node7#ss:psetup_b> IDASpilsPrecSetupFnB
          *)
        and 'a prec_setup_fn = ('a triple_tmp, 'a) jacobian_arg -> unit

        (** Specifies a Jacobian-times-vector function.
            [jac_times_vec_fn arg v jv] should compute the
            matrix-vector product {i J}[v], where {i J} is the system
            Jacobian.

            - [arg] provides the data necessary to compute the Jacobian.
            - [v] is the vector by which the Jacobian must be multiplied.
            - [jv] is the vector in which the result must be stored.

            The Jacobian {i J} (which is not explicitly constructed)
            has ({i i,j}) entry {i dFi/dyj + c*dFi/dy'j} where {i F}
            is the residual function, i.e. the partial derivative of
            the [i]-th equation with respect to the [j]-th component
            of the non-derivative vector.  [c] is the [jac_coef]
            field of [arg] (see {!jacobian_arg}).  See the [Dense]
            {!linear_solver} for a more detailed explanation.

            {b NB:} The elements of [jac], [v], and [Jv] must no
            longer be accessed after [psolve] has returned a result,
            i.e. if their values are needed outside of the function
            call, then they must be copied to separate physical
            structures.

            Raising any kind of exception (including
            {!Sundials.RecoverableFailure}) from this function results
            in the integrator being aborted.

            @idas <node7#ss:jactimesvec_b> IDASpilsJacTimesVecFnB
          *)
        and 'a jac_times_vec_fn =
          ('a single_tmp, 'a) jacobian_arg
          -> 'a
          -> 'a
          -> unit

        (** Specifies a preconditioner, including the type of
            preconditioning to be done (none or left), and a set of
            three callbacks if applicable:

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
              If the user doesn't give such a function, IDA uses a
              default implementation based on difference quotients.

            Like the {!linear_solver}, there are several functions
            which construct preconditioners.  The simplest is
            {!prec_none}, which does no preconditioning.  Arbitrary
            user-defined preconditioners can be constructed through
            {!prec_left}, which takes user-defined [solve], [setup],
            and [jac_times_vec], with the last two optional.
            {!Idas_bbd} gives access to the parallel band-block
            diagonal preconditioners that come with IDAS.

            @idas <node7#SECTION00729400000000000000> IDASpilsSetPreconditionerB
            @idas <node7#SECTION00729400000000000000> IDASpilsSetJacTimesVecFnB
            @idas <node7#ss:psolve_b> IDASpilsPrecSolveFnB
            @idas <node7#ss:psetup_b> IDASpilsPrecSetupFnB
            @idas <node7#ss:jactimesvec_b> IDASpilsJacTimesVecFnB
        *)
        type ('a, 'k) preconditioner =
          ('a, 'k) AdjointTypes.SpilsTypes.preconditioner

        (** See {!preconditioner}.  *)
        val prec_none : ('a, 'k) preconditioner

        (** See {!preconditioner}. *)
        val prec_left :
          ?setup:'a prec_setup_fn
          -> ?jac_times_vec:'a jac_times_vec_fn
          -> 'a prec_solve_fn
          -> ('a, 'k) preconditioner

        (** Krylov iterative linear solver with the scaled
            preconditioned GMRES method.  Called like [spgmr
            ~maxl:maxl ~max_restarts:maxr prec], where:

            - [~maxl] is the maximum dimension of the Krylov subspace.
              Defaults to [5].
            - [~max_restarts] is the maximum number of restarts.
              Defaults to [5].  Passing [0] disables restarts.
            - [prec] is a preconditioner.  See {!preconditioner}.

            @idas <node7#sss:lin_solv_b> IDASpgmrB
            @idas <node7#SECTION00729400000000000000> IDASpilsSetPreconditionerB
            @idas <node7#ss:psolve_b> IDASpilsPrecSolveFnB
            @idas <node7#ss:psetup_b> IDASpilsPrecSetupFnB
            @idas <node7#ss:jactimesvec_b> IDASpilsJacTimesVecFnB
         *)
        val spgmr : ?maxl:int -> ?max_restarts:int
                  -> ('a, 'k) preconditioner -> ('a, 'k) linear_solver

        (** Krylov iterative solver with the scaled preconditioned
            Bi-CGStab method.  The arguments are the same as [Spgmr],
            except the maximum number of restarts ([~max_restarts])
            cannot be specified.

            @idas <node7#sss:lin_solv_b> IDASpbcgB
            @idas <node7#SECTION00729400000000000000> IDASpilsSetPreconditionerB
            @idas <node7#ss:psolve_b> IDASpilsPrecSolveFnB
            @idas <node7#ss:psetup_b> IDASpilsPrecSetupFnB
            @idas <node7#ss:jactimesvec_b> IDASpilsJacTimesVecFnB
         *)
        val spbcg : ?maxl:int -> ('a, 'k) preconditioner
                  -> ('a, 'k) linear_solver

        (** Krylov iterative with the scaled preconditioned TFQMR
            method.  The arguments are the same as [Spgmr], except the
            maximum number of restarts ([~max_restarts]) cannot be
            specified.

            @idas <node7#sss:lin_solv_b> IDASptfqmrB
            @idas <node7#SECTION00729400000000000000> IDASpilsSetPreconditionerB
            @idas <node7#ss:psolve_b> IDASpilsPrecSolveFnB
            @idas <node7#ss:psetup_b> IDASpilsPrecSetupFnB
            @idas <node7#ss:jactimesvec_b> IDASpilsJacTimesVecFnB
         *)
        val sptfqmr : ?maxl:int -> ('a, 'k) preconditioner
                    -> ('a, 'k) linear_solver

        (** {5:adjbwdspilsoptin Optional Input Functions} *)

        (** Sets the Gram-Schmidt orthogonalization to be used with the
            Spgmr {!linear_solver}.

            @idas <node7#SECTION00729400000000000000> IDASpilsSetGSTypeB *)
        val set_gs_type :
          ('a, 'k) bsession -> Spils.gramschmidt_type -> unit

        (** [set_eps_lin eplifac] sets the factor by which the Krylov linear
            solver's convergence test constant is reduced from the Newton
            iteration test constant. [eplifac]  must be >= 0. Passing a value of
            0 specifies the default (which is 0.05).

            @idas <node7#SECTION00729400000000000000> IDASpilsSetEpsLinB *)
        val set_eps_lin : ('a, 'k) bsession -> float -> unit

        (** [set_maxl maxl] resets the maximum Krylov subspace dimension for the
            Bi-CGStab or TFQMR methods. [maxl] is the maximum dimension of the
            Krylov subspace.  A value of [None] (or [maxl] <= 0) specifies the
            default of 5.0.

            @idas <node7#SECTION00729400000000000000> IDASpilsSetMaxlB *)
        val set_maxl : ('a, 'b) bsession -> int option -> unit

        (** {5:adjbwdspilsoptout Optional Output Functions} *)

        (** Returns the sizes of the real and integer workspaces used by the
            linear solver.

            @idas <node5#sss:optout_spils> IDASpilsGetWorkSpace
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
            @return ([real_size], [integer_size])
         *)
        val get_work_space : ('a, 'b) bsession -> int * int

        (** Returns the cumulative number of linear iterations.

            @idas <node5#sss:optout_spils> IDASpilsGetNumLinIters
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_lin_iters : ('a, 'b) bsession -> int

        (** Returns the cumulative number of linear convergence failures.

            @idas <node5#sss:optout_spils> IDASpilsGetNumConvFails
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_conv_fails : ('a, 'b) bsession -> int

        (** Returns the number of preconditioner evaluations, i.e., the number
            of calls made to psetup with jok = [false].

            @idas <node5#sss:optout_spils> IDASpilsGetNumPrecEvals
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_prec_evals : ('a, 'b) bsession -> int

        (** Returns the cumulative number of calls made to the preconditioner
            solve function, psolve.

            @idas <node5#sss:optout_spils> IDASpilsGetNumPrecSolves
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_prec_solves : ('a, 'b) bsession -> int

        (** Returns the cumulative number of calls made to the Jacobian-vector
            function, jtimes.

            @idas <node5#sss:optout_spils> IDASpilsGetNumJtimesEvals
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_jtimes_evals : ('a, 'b) bsession -> int

        (** Returns the number of calls to the user residual function
            for finite difference Jacobian-vector product
            approximation. This counter is only updated if the default
            difference quotient function is used.

            @idas <node5#sss:optout_spils> IDASpilsGetNumRhsEvals
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_res_evals : ('a, 'b) bsession -> int
      end

    (** {4:adjbwdout Output} *)

    (** Returns the real and integer workspace sizes.

        @idas <node5#sss:optout_main> IDAGetWorkSpace
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
        @return ([real_size], [integer_size]) *)
    val get_work_space : ('a, 'b) bsession -> int * int

    (** Returns the cumulative number of internal steps taken by the solver.

        @idas <node5#sss:optout_main> IDAGetNumSteps
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
      *)
    val get_num_steps : ('a, 'b) bsession -> int

    (** Returns the number of calls to the user's residual function.

        @idas <node5#sss:optout_main> IDAGetNumRhsEvals
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
      *)
    val get_num_res_evals : ('a, 'b) bsession -> int

    (** Returns the number of calls made to the linear solver's setup function.

        @idas <node5#sss:optout_main> IDAGetNumLinSolvSetups
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_num_lin_solv_setups : ('a, 'b) bsession -> int

    (** Returns the number of local error test failures that have occurred.

        @idas <node5#sss:optout_main> IDAGetNumErrTestFails
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_num_err_test_fails : ('a, 'b) bsession -> int

    (** Returns the integration method order used during the last internal step.

        @idas <node5#sss:optout_main> IDAGetLastOrder
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_last_order : ('a, 'b) bsession -> int

    (** Returns the integration method order to be used on the next internal
        step.

        @idas <node5#sss:optout_main> IDAGetCurrentOrder
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_current_order : ('a, 'b) bsession -> int

    (** Returns the integration step size taken on the last internal step.

        @idas <node5#sss:optout_main> IDAGetLastStep
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_last_step : ('a, 'b) bsession -> float

    (** Returns the integration step size to be attempted on the next internal
        step.

        @idas <node5#sss:optout_main> IDAGetCurrentStep
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_current_step : ('a, 'b) bsession -> float

    (** Returns the the value of the integration step size used on the first
        step.

        @idas <node5#sss:optout_main> IDAGetActualInitStep
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_actual_init_step : ('a, 'b) bsession -> float

    (** Returns the the current internal time reached by the solver.

        @idas <node5#sss:optout_main> IDAGetCurrentTime
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_current_time : ('a, 'b) bsession -> float

    (** Returns a suggested factor by which the user's tolerances should be
        scaled when too much accuracy has been requested for some internal
        step.

        @idas <node5#sss:optout_main> IDAGetTolScaleFactor
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_tol_scale_factor : ('a, 'b) bsession -> float

    (** Returns the solution error weights at the current time.

        @idas <node5#sss:optout_main> IDAGetErrWeights
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
        @idas <node3#ss:ivp_sol> IVP solution (W_i) *)
    val get_err_weights : ('a, 'b) bsession -> ('a, 'b) Ida.nvector -> unit

    (** Returns the vector of estimated local errors.

        @idas <node5#sss:optout_main> IDAGetEstLocalErrors
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_est_local_errors :
      ('a, 'b) bsession -> ('a, 'b) Ida.nvector -> unit

    (** Returns the integrator statistics as a group.

        @idas <node5#sss:optout_main> IDAGetIntegratorStats
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_integrator_stats : ('a, 'b) bsession -> Ida.integrator_stats

    (** Convenience function that calls get_integrator_stats and prints the
        results to stdout.

        @idas <node5#sss:optout_main> IDAGetIntegratorStats
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val print_integrator_stats : ('a, 'b) bsession -> unit

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed.

        @idas <node5#sss:optout_main> IDAGetNumNonlinSolvIters
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_num_nonlin_solv_iters : ('a, 'b) bsession -> int

    (** Returns the number of nonlinear convergence failures that have occurred.

        @idas <node5#sss:optout_main> IDAGetNumNonlinSolvConvFails
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_num_nonlin_solv_conv_fails : ('a, 'b) bsession -> int

    (** [nniters, nncfails = get_nonlin_solv_stats s] returns both the
        numbers of nonlinear iterations performed [nniters] and of
        nonlinear convergence failures that have occurred [nncfails].

        @idas <node5#sss:optout_main> IDAGetNonlinSolvStats
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
     *)
    val get_nonlin_solv_stats : ('a, 'b) bsession -> int * int

    (** {2:adjquad Quadrature Equations} *)

    module Quadrature :
      sig

        (** {3:adjquadinit Initialization} *)

        (** These functions compute the quadrature equation right-hand side for
            the backward problem. *)
        type 'a bquadrhsfn =
          | NoSens of 'a bquadrhsfn_no_sens
          (** Doesn't depend on forward sensitivities.  See
              {!bquadrhsfn_no_sens} for details.  *)
          | WithSens of 'a bquadrhsfn_with_sens
          (** Depends on forward sensitivities.  See
              {!bquadrhsfn_with_sens} for details.  *)

        (** A quadrature rhs function that does not depend on forward
            sensitivities.  The function is called as [f t y y' yB y'B
            rhsvalBQ], where
            - [t] is the current value of the independent variable.
            - [y] is the current value of the forward solution vector.
            - [y'] is the current value of the forward derivative solution vector.
            - [yB] is the current value of the backward dependent variable vector.
            - [y'B] is the current value of the backward dependent derivative vector.
            - [rhsvalBQ] is the output vector containing the residual for the backward quadrature equations.

            This function can raise {!Sundials.RecoverableFailure} to
            instruct the integrator to retry with a different step
            size.  Raising any other kind of exception aborts the
            integrator.

            See also {!bquadrhsfn}.

            @idas <node7#sss:rhs_quad_B> IDAQuadRhsFnB *)
        and 'a bquadrhsfn_no_sens =
          float
          -> 'a
          -> 'a
          -> 'a
          -> 'a
          -> 'a
          -> unit

        (** A quadrature rhs function that depends on forward
            sensitivities.  The function is called as [f t y y' yS y'S
            yB y'B rhsvalQBS], where
            - [t] is the current value of the independent variable.
            - [y] is the current value of the forward solution vector.
            - [yS] is an array of vectors containing the sensitivities of the forward solution.
            - [y'S] is an array of vectors containing the sensitivities of the forward derivative solution.
            - [y'] is the current value of the forward derivative solution vector.
            - [yB] is the current value of the backward dependent variable vector.
            - [y'B] is the current value of the backward dependent derivative vector.
            - [rhsvalBQ] is the output vector containing the residual for the backward quadrature equations.

            This function can raise {!Sundials.RecoverableFailure} to
            instruct the integrator to retry with a different step
            size.  Raising any other kind of exception aborts the
            integrator.

            See also {!bquadrhsfn}.

            @idas <node7#sss:rhs_quad_sens_B> IDAQuadRhsFnBS *)
        and 'a bquadrhsfn_with_sens =
          float
          -> 'a
          -> 'a
          -> 'a array
          -> 'a array
          -> 'a
          -> 'a
          -> 'a
          -> unit

        (** This function, [init s fQB yQB0], activates integration of
            quadrature equations, with or without sensitivities, where [fQB]
            computes the right-hand side of the backward quadrature equations,
            and [yQB0] contains the values of the quadrature variables at [tB0].

            @idas <node7#SECTION007211000000000000000> IDAQuadInitB
            @idas <node7#SECTION007211000000000000000> IDAQuadInitBS *)
        val init :
          ('a, 'b) bsession -> 'a bquadrhsfn -> ('a, 'b) nvector -> unit

        (** This function reinitializes the integration of quadrature equations
            during the backward phase.

            @idas <node7#SECTION007211000000000000000> IDAQuadReInitB *)
        val reinit : ('a, 'b) bsession -> ('a, 'b) nvector -> unit

        (** {3:adjextraction Extraction function} *)

        (** [tret = get s w yqs] fills [yqs] with the quadrature
            solution vector after a successful return from
            {!backward_normal} or {!backward_one_step}, and returns
            the time reached by the solver.

            @idas <node7#sss:quad_get_b> IDAGetQuadB
         *)
        val get : ('a, 'b) bsession -> ('a, 'b) nvector -> float

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

            @idas <node7#sss:quad_optional_input_B> IDASetQuadErrConB
            @idas <node7#sss:quad_optional_input_B> IDAQuadSStolerancesB
            @idas <node7#sss:quad_optional_input_B> IDAQuadSVtolerancesB *)
        val set_tolerances : ('a, 'b) bsession -> ('a, 'b) tolerance -> unit

        (** {3:adjquadoptout Optional Output Functions} *)

        (** Returns the number of calls to the user's quadrature right-hand side
            function.

            @idas <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @idas <node5#ss:quad_optional_output> IDAGetQuadNumRhsEvals
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_rhs_evals : ('a, 'b) bsession -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @idas <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @idas <node5#ss:quad_optional_output> IDAGetQuadNumErrTestFails
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_num_err_test_fails : ('a, 'b) bsession -> int

        (** Returns the quadrature error weights at the current time.

            @idas <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @idas <node5#ss:quad_optional_output> IDAGetQuadErrWeights
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_err_weights : ('a, 'b) bsession -> ('a, 'b) nvector -> unit

        (** [nfqevals, nqetfails = get_stats s] returns
            - [fqevals], the number of calls to the user's quadrature function, and,
            - [nqetfails], the number of error test failures due to quadrature variables.

            @idas <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration
            @idas <node5#ss:quad_optional_output> IDAGetQuadStats
            @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
         *)
        val get_stats : ('a, 'b) bsession -> int * int
      end
  end
