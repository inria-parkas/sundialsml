(***********************************************************************)
(*                                                                     *)
(*               OCaml interface to (serial) Sundials                  *)
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

(* TODO: add "skeletons" for each solver extension *)

(* TODO: Verify that all functionalities are present and consistent...
         and likewise for CVODE and IDA *)
(* TODO: redo heading links {3:... } *)

(** Abstract nvector interface to the CVODES solver.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

type 'a session = 'a Cvode_nvector.session
type 'a nvector = 'a Cvode_nvector.nvector

(** {2:quadrature Integration of pure quadrature equations} *)

module Quadrature :
  sig

    (** {3:quadexcept Exceptions} *)

    (* TODO: add @raise everywhere? *)
    (** Quadrature integration was not initialized.

        @cvodes <node5#ss:quad_get> CV_NO_QUAD *)
    exception QuadNotInitialized

    (** The quadrature right-hand side function failed in an unrecoverable
        manner.
     
        @cvodes <node5#SECTION00572000000000000000> CV_QRHSFUNC_FAIL *)
    exception QuadRhsFuncFailure

    (** The quadrature right-hand side function failed at the first call.
      
        @cvodes <node5#SECTION00572000000000000000> CV_FIRST_QRHSFUNC_ERR *)
    exception FirstQuadRhsFuncErr

    (** Convergence test failures occurred too many times due to repeated
        recoverable errors in the quadrature right-hand side function. This
        value will also be returned if the quadrature right-hand side function
        had repeated recoverable errors during the estimation of an initial step
        size (assuming the quadrature variables are included in the error
        tests).

        @cvodes <node5#SECTION00572000000000000000> CV_REPTD_QRHSFUNC_ERR *)
    exception RepeatedQuadRhsFuncErr

    (** The quadrature right-hand function had a recoverable error, but no
        recovery was possible. This failure mode is rare, as it can occur only
        if the quadrature right-hand side function fails recoverably after an
        error test failed while at order one.
      
        @cvodes <node5#SECTION00572000000000000000> CV_UNREC_QRHSFUNC_ERR *)
    exception UnrecoverableQuadRhsFuncErr

    (** {3:quadinit Initialization} *)

    (** This function, [f_Q t y yQdot], computes the quadrature equation
        right-hand side [yQdot] for a given value of the independent variable
        [t] and state vector [y].

        @cvodes <node5#ss:user_fct_quad> CVQuadRhsFn *)
    type 'a quadrhsfn = float -> 'a -> 'a -> unit

    (** Activates the integration of quadrature equations.

        @cvodes <node5#ss:quad_malloc> CVodeQuadInit *)
    val init : 'a session -> 'a quadrhsfn -> 'a nvector -> unit

    (** Reinitialize the integration of quadrature equations.

        @cvodes <node5#ss:quad_malloc> CVodeQuadReInit *)
    val reinit : 'a session -> 'a nvector -> unit

    (** {3:quadtol Tolerance specification} *)

    type 'a tolerance =
        NoStepSizeControl
        (** Do not use quadrature variables for step-size control (default). *)
      | SSTolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVTolerances of float * 'a nvector
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

    (** Specify whether and how quadrature variables should be used in the step
        size control mechanism.

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
    val set_tolerances : 'a session -> 'a tolerance -> unit

    (** {3:quadout Output functions} *)

    (** [tret = get s yq] fills [yq] with the quadrature solution vector after a
        successful return from {!Cvode.solve_normal} or {!Cvode.solve_one_step},
        and returns the time reached by the solver.

        @cvodes <node5#ss:quad_get> CVodeGetQuad *)
    val get : 'a session -> 'a nvector -> float

    (** [tret = get_dky s t k dkyq] fills [dkyq] with the derivatives of the
        quadrature solution vector after a successful return from
        {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested,
        [t], must fall within the interval defined by the last successful step
        ({!Cvode.get_last_step}). The requested order, [k], must be less than or
        equal to the value returned by {!Cvode.get_last_order}.

        @cvodes <node5#ss:quad_get> CVodeGetQuadDky *)
    val get_dky : 'a session -> float -> int -> 'a nvector -> unit

    (** {3 Optional output functions} *)

    (** Returns the number of calls to the user's quadrature right-hand side
        function.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals *)
    val get_num_rhs_evals       : 'a session -> int

    (** Returns the number of local error test failures due to quadrature
        variables.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails *)
    val get_num_err_test_fails  : 'a session -> int

    (** Returns the quadrature error weights at the current time.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights *)
    val get_err_weights : 'a session -> 'a nvector -> unit

    (** [nfqevals, nqetfails = get_stats s] returns
        - [fqevals], the number of calls to the user's quadrature function, and,
        - [nqetfails], the number of error test failures due to quadrature variables.

        @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats *)
    val get_stats : 'a session -> int * int

  end

(** {2:sens (Forward) Sensitivity Analysis} *)

module Sensitivity :
  sig
    (** {3:senstol Tolerance specification} *)

    type 'a tolerance =
        SSTolerances of float * Sundials.real_array
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVTolerances of float * 'a nvector array
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
      | EETolerances
        (** Calculate the integration tolerances for sensitivities based
            on those for state variables and the scaling factors. *)

    (** {3:sensexcept Exceptions} *)

    (* TODO: add @raise everywhere? *)
    (** Sensitivity analysis was not initialized.

        @cvodes <node5#ss:sensi_get> CV_NO_SENS *)
    exception SensNotInitialized

    (** @cvodes <node6#SECTION00623000000000000000> CV_SRHSFUNC_FAIL *)
    exception SensRhsFuncFailure

    (** @cvodes <node6#SECTION00623000000000000000> CV_FIRST_SRHSFUNC_ERR *)
    exception FirstSensRhsFuncErr

    (** @cvodes <node6#SECTION00623000000000000000> CV_REPTD_SRHSFUNC_ERR *)
    exception RepeatedSensRhsFuncErr

    (** @cvodes <node6#SECTION00623000000000000000> CV_UNREC_SRHSFUNC_ERR *)
    exception UnrecoverableSensRhsFuncErr

    (** @cvodes <node6> CV_BAD_IS *)
    exception BadIS

    (** {3:sensinit Initialization} *)

    type sens_method =
        Simultaneous
        (** Correct state and sensitivity variables at the same time.
            If {!Cvode_nvector.Newton} was selected as the nonlinear system
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
            this approach is equivalent to {!Staggererd}. Note that this
            approach can only be used if the user-provided sensitivity
            right-hand side function is of type {!OneByOne} (CV_STAGGERED1). *)

    type 'a sensrhsfn =
        AllAtOnce of
         (float           (* t *)
           -> 'a          (* y *)
           -> 'a          (* ydot *)
           -> 'a array    (* yS *)
           -> 'a array    (* ySdot *)
           -> 'a          (* tmp1 *)
           -> 'a          (* tmp2 *)
           -> unit)
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
 
            @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn *)
      | OneByOne of
         (float           (* t *)
           -> 'a          (* y *)
           -> 'a          (* ydot *)
           -> int         (* iS *)
           -> 'a          (* yS *)
           -> 'a          (* ySdot *)
           -> 'a          (* tmp1 *)
           -> 'a          (* tmp2 *)
           -> unit)
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

           @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn *)

    (** Used for specifying problem parameter information for
        sensitivity calculations.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensParams *)
    type sens_params = {
        pvals  : Sundials.real_array option;
        (** The parameters used to evaluate {i f(t, y, p)}. *)
        pbar   : Sundials.real_array option;
        (** An array of {i ns} positive scaling factors. *)
        plist  : int array option;
        (** An array of non-negative indices to specify which components
            to use in estimating the sensitivity equations. *)
      }

    val no_sens_params : sens_params

    (** This function, [init s ism ps fS yS0], activates the forward sensitivity
        computation, where [ism] selects the sensitivity solution method, [ps]
        gives problem parameter information, [fS] computes the sensitivity
        right-hand sides, and [yS0] gives the initial values of the
        sensitivities.
        
        Note that any array specified by [ps.pvals] is used to pass parameter
        information during problem solution, that is, the library will normally
        write to it from time to time.

        @cvodes <node6#ss:sensi_malloc> CVodeSensInit
        @cvodes <node6#ss:sensi_malloc> CVodeSensInit1 *)
    val init : 'a session
               -> 'a tolerance
               -> sens_method
               -> sens_params
               -> 'a sensrhsfn
               -> 'a nvector array
               -> unit

    (** This function reinitializes the forward sensitivity computation.

        @cvodes <node6#ss:sensi_malloc> CVodeSensReInit *)
    val reinit : 'a session -> sens_method -> 'a nvector array -> unit

    (** Deactivates forward sensitivity calculations without deallocating
        memory. Sensitivities can be reactivated with {!reinit}.

        @cvodes <node6#ss:sensi_malloc> CVodeSensToggleOff *)
    val toggle_off : 'a session -> unit

    (** {3:output Output functions} *)

    (** [tret = get s ys] fills [ys] with the sensitivity solution vectors after
        a successful return from {!Cvode.solve_normal} or
        {!Cvode.solve_one_step}, and returns the time reached by the solver.

        @cvodes <node6#ss:sensi_get> CVodeGetSens *)
    val get : 'a session -> 'a nvector array -> float

    (** [tret = get_dky s t k dkys] fills [dkys] with the derivatives of the
        sensitivity solution vectors after a successful return from
        {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested,
        [t], must fall within the interval defined by the last successful step
        ({!Cvode.get_last_step}). The requested order, [k], must be less than or
        equal to the value returned by {!Cvode.get_last_order}.

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky *)
    val get_dky : 'a session -> float -> int -> 'a nvector array -> unit

    (** [tret = get s i ys] fills [ys] with the [i]th sensitivity solution
        vector after a successful return from {!Cvode.solve_normal} or
        {!Cvode.solve_one_step}, and returns the time reached by the solver.

        @cvodes <node6#ss:sensi_get> CVodeGetSens1 *)
    val get' : 'a session -> int -> 'a nvector -> float

    (** [tret = get_dky s t k i dkys] fills [dkys] with the derivatives of the
        [i]th sensitivity solution vector after a successful return from
        {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested,
        [t], must fall within the interval defined by the last successful step
        ({!Cvode.get_last_step}). The requested order, [k], must be less than or
        equal to the value returned by {!Cvode.get_last_order}.

        @cvodes <node6#ss:sensi_get> CVodeGetSensDky1 *)
    val get_dky' : 'a session -> float -> int -> int -> 'a nvector -> unit

    (** {3 Optional input functions} *)

    (** Specify the integration tolerances for sensitivities.

        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSStolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensSVtolerances
        @cvodes <node6#ss:cvfwdtolerances> CVodeSensEEtolerances *)
    val set_tolerances : 'a session -> 'a tolerance -> unit

    (** Set whether sensitivity variables should be used in the error control
        mechanism (the default is [true]). 

        @cvodes <node5#ss:sens_optional_input> CVodeSetSensErrCon *)
    val set_err_con : 'a session -> bool -> unit

    type dq_method = DQCentered (* CV_CENTERED *)
                   | DQForward  (* CV_FORWARD *)

    (** [set_dq_method s dqtype dqrhomax] specifies the difference quotient
        strategy in the case in which the right-hand side of the sensitivity
        equations is to be computed by CVODES; [dqrhomax] is used in deciding
        the switching between simultaneous or separate approximations of the two
        tersm in the sensitivity right-hand side.

        @cvodes <node6#ss:sens_optional_input> CVodeSetSensDQMethod
        @cvodes <node3#ss:fwd_sensi> Forward Sensitivity Analysis *)
    val set_dq_method : 'a session -> dq_method -> float -> unit

    (** Specifies the maximum number of nonlinear solver iterations for
        sensitivity variables permitted per step.

        @cvode <node5#ss:sens_optional_input> CVodeSetSensMaxNonlinIters *)
    val set_max_nonlin_iters : 'a session -> int -> unit

    (** {3 Optional output functions} *)

    (** Returns the number of calls to the sensitivity right-hand side function.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumRhsEvals *)
    val get_num_sens_evals       : 'a session -> int

    (** Returns the number of calls to the user's right-hand side function due
        to the internal finite difference approximation of the sensitivity
        right-hand sides.

        @cvode <node6#ss:sens_optional_output> CVodeGetNumRhsEvalsSens *)
    val get_num_rhs_evals      : 'a session -> int

    (** Returns the number of local error test failures for the sensitivity
        variables that have occurred.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumErrTestFails *)
    val get_num_err_test_fails  : 'a session -> int

    (** Returns the number of calls made to the linear solver's setup function
        due to forward sensitivity calculations.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumLinSolvSetups *)
    val get_num_lin_solv_setups : 'a session -> int

    type sensitivity_stats = {
        num_rhs_evals : int;
        num_sens_evals :int;
        num_err_test_fails : int;
        num_lin_solv_setups :int;
      }

    (** Returns all of the sensitivity-related solver statistics as a group.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensStats *)
    val get_stats : 'a session -> sensitivity_stats

    (** Returns the sensitivity error weight vectors at the current time.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensErrWeights
        @cvode <node3#e:errwt> Eq. (2.7) IVP solution (W_i) *)
    val get_err_weights : 'a session -> 'a nvector array -> unit

    (** Returns the number of nonlinear iterations performed for sensitivity
        calculations.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvIters *)
    val get_num_nonlin_solv_iters : 'a session -> int

    (** Returns the number of nonlinear convergence failures that have occurred
        for sensitivity calculations.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvConvFails *)
    val get_num_nonlin_solv_conv_fails : 'a session -> int

    type nonlin_stats = {
        num_nonlin_solv_iters : int;
        num_nonlin_solv_conv_fails : int;
      }

    (** Returns the sensitivity-related nonlinear solver statistics as a group.

        @cvode <node6#ss:sens_optional_output> CVodeGetSensNonlinSolvStats
     *)
    val get_nonlin_solv_stats : 'a session -> nonlin_stats

    (** Returns the number of nonlinear (functional or Newton) iterations
        performed for each sensitivity equation separately, in the {!Staggered1}
        case.

      @cvode <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvIters *)
    val get_num_stgr_nonlin_solv_iters : 'a session

                                         -> Sundials.lint_array -> unit

    (** Returns the number of nonlinear convergence failures that have occurred
        for each sensitivity equation separately, in the {!Staggered1} case.

      @cvode <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvConvFails *)
    val get_num_stgr_nonlin_solv_conv_fails : 'a session
                                              -> Sundials.lint_array -> unit

    (** {3 Integration of quadrature equations depending on forward sensitivities} *)
    module Quadrature :
      sig

        (** {4:quadsensexcept Exceptions} *)

        (* TODO: add @raise everywhere? *)
        (** Quadrature integration was not initialized.

            @cvodes <node5#SECTION00642000000000000000> CV_NO_QUAD_SENS *)
        exception QuadSensNotInitialized

        (** @cvodes <node6#SECTION00642000000000000000> CV_QSRHSFUNC_FAIL *)
        exception QuadSensRhsFuncFailure

        (** @cvodes <node6#SECTION00642000000000000000> CV_FIRST_QSRHSFUNC_ERR *)
        exception FirstQuadSensRhsFuncErr

        (** @cvodes <node6#SECTION00642000000000000000> CV_REPTD_QSRHSFUNC_ERR *)
        exception RepeatedQuadSensRhsFuncErr

        (** @cvodes <node6#SECTION00642000000000000000> CV_UNREC_QSRHSFUNC_ERR *)
        exception UnrecoverableQuadSensRhsFuncErr

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
            quadrature equations depending on sensitivities, where [fQS]
            computes the right-hand side of the sensitivity-dependent quadrature
            equations, and [yQS0] contains the initial values of
            sensitivity-dependent quadratures.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensInit *)
        val init : 'a session -> 'a quadsensrhsfn -> 'a nvector array -> unit

        (** This function reinitializes the forward sensitivity computation.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadSensReInit *)
        val reinit : 'a session -> 'a nvector array -> unit

        (** {4:sensquad Tolerance specification} *)

        type 'a tolerance =
            NoStepSizeControl
            (** Do not use quadrature variables for step-size control (default). *)
          | SSTolerances of float * Sundials.real_array
            (** [(rel, abs)] : scalar relative and absolute tolerances. *)
          | SVTolerances of float * 'a nvector array
            (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
          | EETolerances
            (** Estimate the tolerances for the sensitivity-dependent
                quadratures from those provided for the pure quadrature
                variables. *)

        (** Specify whether and how quadrature variables should be used in the
            step size control mechanism.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeSetQuadSensErrCon
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSStolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSVtolerances
            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensEEtolerances *)
        val set_tolerances : 'a session -> 'a tolerance -> unit

        (** {4:extraction Extraction functions} *)

        (** [tret = get s yqs] fills [yqs] with quadrature solution vectors
            after a successful return from {!Cvode.solve_normal} or
            {!Cvode.solve_one_step}, and returns the time reached by the solver.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens *)
        val get : 'a session -> 'a nvector array -> float

        (**
          [tret = get s i yqs] fills [yqs] with the [i]th quadrature solution
          vector after a successful return from {!Cvode.solve_normal} or
          {!Cvode.solve_one_step}, and returns the time reached by the solver.

          @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens1
         *)
        val get1 : 'a session -> int -> 'a nvector -> float

        (**
          [tret = get_dky s t k dkyqs] fills [dkyqs] with the derivatives of the
          quadrature solution vectors after a successful return from
          {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested, [t],
          must fall within the interval defined by the last successful step
          ({!Cvode.get_last_step}). The requested order, [k], must be less than or
          equal to the value returned by {!Cvode.get_last_order}.

          @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky
         *)
        val get_dky : 'a session -> float -> int -> 'a nvector array -> unit

        (** [tret = get_dky s t k i dkyqs] fills [dkyqs] with the derivatives of
            the [i]th quadrature solution vector after a successful return from
            {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time
            requested, [t], must fall within the interval defined by the last
            successful step ({!Cvode.get_last_step}). The requested order, [k],
            must be less than or equal to the value returned by
            {!Cvode.get_last_order}.

            @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky1 *)
        val get_dky1 : 'a session -> float -> int -> int -> 'a nvector array -> unit

        (** {4 Optional output functions} *)

        (** Returns the number of calls to the user's quadrature right-hand side
            function.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumRhsEvals *)
        val get_num_rhs_evals       : 'a session -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumErrTestFails *)
        val get_num_err_test_fails  : 'a session -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensErrWeights *)
        val get_err_weights : 'a session -> 'a nvector array -> unit

        (** [nfqevals, nqetfails = get_stats s] returns
            - [fqevals], the number of calls to the user's quadrature function,
            and,
            - [nqetfails], the number of error test failures due to quadrature
            variables.

          @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensStats *)
        val get_stats : 'a session -> int * int

      end
    end

(** {2:adjoint Adjoint Sensitivity Analysis} *)

module Adjoint :
  sig

    (** {3:adjexcept Exceptions} *)

    (* TODO: add @raise everywhere? *)
    (* TODO: check that all return flags have been considered; here and
       elsewhere *)
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
    exception ForwardReinitializationFailed

    (** An error occured during the integration of the forward problem.

        @cvodes <node7#sss:cvsolveb> CV_FWD_FAIL *)
    exception ForwardFailed

    (** No backward problem has been created.

        @cvodes <node7#sss:cvsolveb> CV_NO_BCK *)
    exception NoBackwardProblem

    (** The final time [tB0] was outside the interval over which the forward
        problem was solved.

        @cvodes <node7#sss:cvinitb> CV_BAD_TB0 *)
    exception BadTB0

    (** {3:adjinit Initialization and Forward Solution} *)

    (** Specifies the type of interpolation.

        @cvodes <node3#ss:checkpointing> Checkpointing scheme *)
    type interpolation = IPolynomial (* CV_POLYNOMIAL *)
                       | IHermite    (* CV_HERMITE *)

    (** [init s nd interp] initializes the forward-backward problem with [nd]
        integration steps between consecutive checkpoints and variable-degree
        interpolation according to [interp]. This function must be called before
        either {!forward_normal} or {!forward_one_step}.

        @cvodes <node7#ss:cvadjinit> CVodeAdjInit *)
    val init : 'a session -> int -> interpolation -> unit

    (** {3:adjforward Forward integration functions} *)

    (** [tret, ncheck = forward_normal s tout yret] integrates the forward
        problem over an interval and saves checkpointing data. The function
        takes as arguments the next time at which a solution is desired
        ([tout]), a vector for storing the computed result ([yret]), and returns
        the time reached by the solver ([tret]) and the number of checkpoints
        stored so far ([ncheck]).

        This call asks the solver to take internal steps until it has reached or
        just passed the [tout] parameter ([CV_NORMAL]). The solver then
        interpolates in order to return an approximate value of [y(tout)].

        @cvodes <node7#sss:cvsolvef> CVodeF
        TODO: list of exceptions raised. *)
    val forward_normal :
      'a session
      -> float
      -> 'a nvector
      -> float * int

    (** [tret, ncheck = forward_normal s tout yret] integrates the forward
        problem over an interval and saves checkpointing data. The function
        takes as arguments the next time at which a solution is desired
        ([tout]), a vector for storing the computed result ([yret]), and returns
        the time reached by the solver ([tret]) and the number of checkpoints
        stored so far ([ncheck]).

        This call asks the solver to take one internal step and to return the
        solution at the point reached by that step ([CV_ONE_STEP]).

        @cvodes <node7#sss:cvsolvef> CVodeF
        TODO: list of exceptions raised. *)
    val forward_one_step :
      'a session
      -> float
      -> 'a nvector
      -> float * int

    (** {3:adjbackinit Initialization of Backward Problems} *)

    (** These functions evaluate the right-hand side of the backward ODE system
        with or without a dependence on forward sensitivies.

        @cvodes <node7#ss:ODErhs_b> CVRhsFnB
        @cvodes <node7#ss:ODErhs_bs> CVRhsFnBS
        @cvodes <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis
        @cvodes <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis *)
    type 'a brhsfn =
        BackBasic of (float       (* t *)
                       -> 'a      (* y *)
                       -> 'a      (* yb *)
                       -> 'a      (* ybdot *)
                       -> unit)
      | BackWithSens of (float          (* t *)
                          -> 'a         (* y *)
                          -> 'a array   (* ys *)
                          -> 'a         (* yb *)
                          -> 'a         (* ybdot *)
                          -> unit)

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

    (** Specify a solution method.

        @cvodes <node7#sss:cvinitb> CVodeCreateB
     *)
    type 'a iter =
      | Newton of 'a linear_solver (** Newton iteration with a given linear solver *)
      | Functional              (** Functional iteration (non-stiff systems only) *)

    (** Specify a linear solver.

        The Lapack solvers require that both Sundials and the OCaml interface
        were built to link with a LAPACK library.

        The Banded Krylov solvers imply an additional call to
        {{:CVODE_DOC_ROOT(node7#SECTION00741000000000000000)} CVBandPrecInitB}.

        @cvodes <node7#sss:lin_solv_b> Linear Solver Initialization Functions *)
    and 'a linear_solver =
      | Diag
      (** Diagonal approximation of the Jacobian by difference quotients.

          @cvodes <node7#sss:lin_solv_b> CVDiagB *)
(*
      (* TODO: not for nvectors *)
      | Dense of bdense_jac_fn option
      (** Direct linear solver with dense matrix.  The optional argument specifies
          a callback function that computes an approximation to the Jacobian matrix
          (see {!bdense_jac_fn} for details).  If this argument is [None], then
          CVODE uses a default implementation based on difference quotients.  See
          also {!Dls}.

          @cvodes <node7#sss:lin_solv_b> CVDenseB
          @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
          @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
      (* TODO: not for nvectors *)
      | LapackDense of bdense_jac_fn option
      (** Direct linear solver with dense matrix, using LAPACK.  The argument is
          the same as [Dense].  See also {!Dls}.

          @cvodes <node7#sss:lin_solv_b> CVLapackDenseB
          @cvodes <node7#SECTION00728200000000000000> CVDlsSetDenseJacFnB
          @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
      (* TODO: not for nvectors *)
      | Band of bandrange * bband_jac_fn option
      (** Direct linear solver with banded matrix.  The arguments specify the width
          of the band ({!bandrange}) and an optional Jacobian function
          ({!bband_jac_fn}).  If the Jacobian function is [None], CVODES uses an
          internal implementation based on difference quotients.  See also {!Dls}.

          @cvodes <node7#sss:lin_solv_b> CVBandB
          @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
          @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
      (* TODO: not for nvectors *)
      | LapackBand of bandrange * bband_jac_fn option
      (** Direct linear solver with banded matrix using LAPACK.  The arguments
          are the same as [Band].

          @cvodes <node7#sss:lin_solv_b> CVLapackBandB
          @cvodes <node7#SECTION00728300000000000000> CVDlsSetBandJacFnB
          @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
*)
      | Spgmr of spils_params * 'a spils_callbacks
      (** Krylov iterative solver with the scaled preconditioned GMRES method.  The
          arguments specify the maximum dimension of the Krylov subspace and
          preconditioning type ({!spils_params}) and the preconditioner callback
          functions ({!spils_callbacks}).  See also {!Spils}.

          @cvodes <node7#sss:lin_solv_b> CVSpgmrB
          @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
          @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
          @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
      | Spbcg of spils_params * 'a spils_callbacks
      (** Krylov iterative solver with the scaled preconditioned Bi-CGStab method.
          The arguments are the same as [Spgmr].  See also {!Spils}.

          @cvodes <node7#sss:lin_solv_b> CVSpbcgB
          @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
          @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
          @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
      | Sptfqmr of spils_params * 'a spils_callbacks
(*
      (** Krylov iterative with the scaled preconditioned TFQMR method.  The
          arguments are the same as [Spgmr].  See also {!Spils}.

          @cvodes <node7#sss:lin_solv_b> CVSptfqmrB
          @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPreconditionerB
          @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
          @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)
      (* TODO: not for nvectors *)
      | BandedSpgmr of spils_params * bandrange
      (** Same as Spgmr (the Krylov iterative solver with scaled preconditioned
          GMRES), but the preconditioner is set to CVODE's internal implementation
          using a banded matrix of difference quotients.  The arguments specify the
          maximum dimension of the Krylov subspace and preconditioning type
          ({!spils_params}), along with the width of the band matrix
          ({!bandrange}).

          @cvodes <node7#sss:lin_solv_b> CVSpgmrB
          @cvodes <node7#SECTION00741000000000000000> CVBandPrecInitB *)
      (* TODO: not for nvectors *)
      | BandedSpbcg of spils_params * bandrange
      (** Same as Spbcg (the Krylov iterative solver with scaled preconditioned
          Bi-CGStab), but the preconditioner is set to CVODE's internal
          implementation using a banded matrix of difference quotients.  The
          arguments are the same as [BandedSpgmr].

          @cvodes <node7#sss:lin_solv_b> CVSpbcgB
          @cvodes <node7#SECTION00741000000000000000> CVBandPrecInitB *)
      (* TODO: not for nvectors *)
      | BandedSptfqmr of spils_params * bandrange
      (** Same as Spbcg (the Krylov iterative solver with scaled preconditioned
          Bi-CGStab), but the preconditioner is set to CVODE's internal
          implementation using a banded matrix of difference quotients.  The
          arguments are the same as [BandedSpgmr].

          @cvodes <node7#sss:lin_solv_b> CVSptfqmrB
          @cvodes <node7#SECTION00741000000000000000> CVBandPrecInitB *)
*)

(*
(* TODO: not for nvectors *)
    (** This function computes the dense Jacobian of the backward problem (or an
        approximation to it).

        @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB *)
    and 'a bdense_jac_fn =
      ('a triple_tmp, 'a) jacobian_arg -> Dls.DenseMatrix.t -> unit

    (** This function computes the banded Jacobian of the backward problem (or
        an approximation to it).

        @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB *)
    and 'a bband_jac_fn =
      bandrange -> ('a triple_tmp, 'a) jacobian_arg -> Dls.DenseMatrix.t -> unit

    (** The range of nonzero entries in a band matrix.  *)
    and bandrange = { mupper : int; (** The upper half-bandwidth.  *)
                      mlower : int; (** The lower half-bandwidth.  *) }
*)

    (** Common parameters for Krylov subspace linear solvers.  *)
    and spils_params = { maxl : int option;
                         (** Maximum dimension of the Krylov subspace
                             to be used.  Pass [None] to use the default
                             value [5]. *)
                         prec_type : Spils.preconditioning_type;
                         (** The type of preconditioning to be done.  *) }

    (** Callbacks for Krylov subspace linear solvers.  Ignored if the
        {!Spils.preconditioning_type} is set to [PrecNone].  In that case, you
        should use {!spils_no_precond} as [spils_callbacks].  *)
    and 'a spils_callbacks =
      {
        prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                         -> 'a -> unit) option;
        (** This function solves the preconditioning system {i Pz = r} for the
            backward problem.

            @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB *)

        prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> bool -> float
                         -> bool) option;
        (** This function preprocesses and/or evaluates Jacobian-related data needed
            by the preconditioner for the backward problem.

            @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB *)

        jac_times_vec_fn : (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a
                            -> unit) option;
        (** This function computes the action of the Jacobian for the backward
            problem on a given vector.

            @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB *)
      }

    (** Arguments passed to the preconditioner solve callback function.  See
        [prec_solve_fn] in {!spils_callbacks}.

        @cvode <node7#ss:psolveFn> CVSpilsPrecSolveFnB *)
    and 'a prec_solve_arg =
      {
        rvec   : 'a;       (** The right-hand side vector, {i r}, of the
                               linear system. *)
        gamma : float;     (** The scalar {i g} appearing in the Newton
                               matrix given by M = I - {i g}J. *)
        delta : float;     (** Input tolerance to be used if an iterative method
                               is employed in the solution. *)

        left  : bool;      (** Indicates whether to use the left preconditioner
                               ([true]) or the right one ([false]). *)
      }

    type 'a tolerance =
      | SSTolerances of float * float
        (** [(rel, abs)] : scalar relative and absolute tolerances. *)
      | SVTolerances of float * 'a nvector
        (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

    (** Identifies a backward problem. *)
    type 'a bsession

    (** [init_backward s lmm iter tol fB tB0 yB0] adds and initializes a
        backward problem that may or may not depend on forward sensitivities,
        where
        - [s] is the parent session (going forward),
        - [lmm]     specifies the linear multistep method, see {!Cvode.lmm},
        - [iter]    specifies either functional iteration or Newton iteration
                    with a specific linear solver, see {!iter},
        - [tol]     specifies the tolerances, see {!tolerance},
        - [fB]      computes the right-hand side of the backward ODE problem,
        - [nroots]  specifies the number of root functions (zero-crossings).
        - [tB0]     specifies the endpoint where final conditions are provided
                    for the backward problem, normally equal to the endpoint of
                    the forward integration, and,
        - [yB0]     is the final value of the backward problem.

        @cvodes <node7#sss:cvinitb> CVodeCreateB
        @cvodes <node7#sss:cvinitb> CVodeInitB
        @cvodes <node7#sss:cvinitb> CVodeInitBS
        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB *)
    val init_backward :
         'a session
      -> Cvode.lmm
      -> 'a iter
      -> 'a tolerance
      -> 'a brhsfn
      -> float
      -> 'a nvector
      -> 'a bsession

    (** Reinitialize the backward problem.

        @cvodes <node7#sss:cvinitb> CVodeReInitB *)
    val reinit : 'a bsession -> float -> 'a nvector -> unit

    (** {3:adjbackward Backward integration functions} *)

    (** [backward_normal s tbout] integrates the backward ODE problem. The
        function takes internal steps until it has reached or just passed the
        user-specified value [tbout] ([CV_NORMAL]). The solver then interpolates
        in order to return an approximate value of [y(tbout)].

        @cvodes <node7#sss:cvsolveb> CVodeB
        TODO: list of exceptions raised. *)
    val backward_normal : 'a session -> float -> unit

    (** [backward_one_step s tbout] integrates the backward ODE problem. The
        function takes one internal step ([CV_ONE_STEP]).

        @cvodes <node7#sss:cvsolveb> CVodeB
        TODO: list of exceptions raised. *)
    val backward_one_step : 'a session -> float -> unit

    (** [tret = get bs yb] returns the solution of the backward ODE problem
        in [yb] at time [tret].

        @cvodes <node7#sss:cvsolveb> CVodeGetB *)
    val get : 'a bsession -> 'a nvector -> float

    (** [tret = get_dky s t k dkys] fills [dkys] with the derivatives of the
        sensitivity solution vectors after a successful return from
        {!backward_normal} or {!backward_one_step}. The time requested, [t],
        must fall within the interval defined by the last successful step
        ({!get_last_step}). The requested order, [k], must be less than or equal
        to the value returned by {!get_last_order}.

        @cvodes <node5#ss:optional_dky> CVodeGetDky
     *)
    val get_dky : 'a bsession -> float -> int -> 'a nvector -> unit

    (** {3 Optional output functions} *)

    (** Instructs {!forward_normal} and {!forward_one_step} not to save
        checkpointing data for forward sensitivities anymore.

        @cvodes <node7#SECTION00727000000000000000> CVodeAdjSetNoSensi *)
    val set_no_sensitivity : 'a session -> unit

    (** Specify the integration tolerances for the backward problem.

        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB *)
    val set_tolerances : 'a bsession -> 'a tolerance -> unit

    (** Specifies the maximum order of the linear multistep method.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxOrdB *)
    val set_max_ord : 'a bsession -> int -> unit

    (** Specifies the maximum number of steps to be taken by the solver in its
        attempt to reach the next output time.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxNumStepsB *)
    val set_max_num_steps : 'a bsession -> int -> unit

    (** Specifies the initial step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetInitStepB *)
    val set_init_step : 'a bsession -> float -> unit

    (** Specifies a lower bound on the magnitude of the step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetMinStepB *)
    val set_min_step : 'a bsession -> float -> unit

    (** Specifies an upper bound on the magnitude of the step size.

        @cvodes <node7#ss:optional_input_b> CVodeSetMaxStepB *)
    val set_max_step : 'a bsession -> float -> unit

    (** Indicates whether the BDF stability limit detection algorithm should be
        used.

        @cvode <node7#ss:optional_input_b> CVodeSetStabLimDet *)
    val set_stab_lim_det : 'a bsession -> bool -> unit

    (** {3 Optional input functions} *)

    (** {3 Diagonal approximation} *)

    (** Get optional inputs for the linear solver that gives diagonal approximations
        of the Jacobian matrix.

        @cvodes <node7#ss:optional_output_b> Optional output functions for
                                             the backward problem. *)
    module Diag :
      sig
        (** {4 Optional input functions} *)

        (** Returns the number of calls made to the user-supplied right-hand
            side function due to finite difference Jacobian approximation in the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetWorkSpace
            @return ([real_size], [integer_size]) *)
        val get_work_space : 'a bsession -> int * int

        (** Returns the number of calls made to the user-supplied right-hand
            side function due to finite difference Jacobian approximation in the
            Diagonal linear solver.

            @cvodes <node5#sss:optout_diag> CVDiagGetNumRhsEvals *)
        val get_num_rhs_evals : 'a bsession -> int
      end

    (** {3 Scaled Preconditioned Iterative Linear Solvers (SPILS)} *)

    (** Set optional outputs, and get optional inputs for the Scaled
        Preconditioned Iterative Linear Solvers: SPGMR, SPBCG, SPTFQMR.

        @cvodes <node7#ss:optional_output_b> Optional output functions for
                                             the backward problem. *)
    module Spils :
      sig
        (** {4 Optional input functions} *)

        (** This function resets the type of preconditioning to be used using a
            value of type {!Spils.preconditioning_type}.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetPrecTypeB *)
        val set_prec_type : 'a bsession -> Spils.preconditioning_type -> unit

        (** Sets the Gram-Schmidt orthogonalization to be used with the
            Spgmr {!linear_solver}.

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetGSTypeB *)
        val set_gs_type : 'a bsession -> Spils.gramschmidt_type -> unit

        (** [set_eps_lin eplifac] sets the factor by which the Krylov linear
            solver's convergence test constant is reduced from the Newton
            iteration test constant. [eplifac]  must be >= 0. Passing a value of
            0 specifies the default (which is 0.05).

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetEpsLinB *)
        val set_eps_lin : 'a bsession -> float -> unit

        (** [set_maxl maxl] resets the maximum Krylov subspace dimension for the
            Bi-CGStab or TFQMR methods. [maxl] is the maximum dimension of the
            Krylov subspace. a value of [None] (or [maxl] <= 0) specifies the
            default of 5.0).

            @cvodes <node7#SECTION00728400000000000000> CVSpilsSetMaxlB *)
        val set_maxl : 'a bsession -> int option -> unit

        (** {4 Optional output functions} *)

        (** Returns the sizes of the real and integer workspaces used by the
            linear solver.

            @cvodes <node5#sss:optout_spils> CVSpilsGetWorkSpace
            @return ([real_size], [integer_size]) *)
        val get_work_space       : 'a bsession -> int * int

        (** Returns the cumulative number of linear iterations.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumLinIters *)
        val get_num_lin_iters    : 'a bsession -> int

        (** Returns the cumulative number of linear convergence failures.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumConvFails *)
        val get_num_conv_fails   : 'a bsession -> int

        (** Returns the number of preconditioner evaluations, i.e., the number
            of calls made to psetup with jok = [false] (see
            {!set_preconditioner}).

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecEvals *)
        val get_num_prec_evals   : 'a bsession -> int

        (** Returns the cumulative number of calls made to the preconditioner
            solve function, psolve (see {!set_preconditioner}).

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumPrecSolves *)
        val get_num_prec_solves  : 'a bsession -> int

        (** Returns the cumulative number of calls made to the Jacobian-vector
            function, jtimes (see {! set_jac_times_vec_fn}).

          @cvodes <node5#sss:optout_spils> CVSpilsGetNumJtimesEvals *)
        val get_num_jtimes_evals : 'a bsession -> int

        (** Returns the number of calls to the user right-hand side function for
            finite difference Jacobian-vector product approximation. This
            counter is only updated if the default difference quotient function
            is used.

            @cvodes <node5#sss:optout_spils> CVSpilsGetNumRhsEvals *)
        val get_num_rhs_evals    : 'a bsession -> int
      end

    (** {3 Banded preconditioner} *)

    (** Get optional outputs for the banded preconditioner module of the
        Scaled Preconditioned Iterative Linear Solvers: SPGMR, SPBCG, SPTFQMR.

        @cvodes <node7#ss:optional_output_b> Optional output functions for
                                             the backward problem. *)
    module BandPrec :
      sig
        (** {4 Optional output functions} *)

        (** Returns the sizes of the real and integer workspaces used by the
            serial banded preconditioner module.

            @cvodes <node5#sss:cvbandpre> CVBandPrecGetWorkSpace
            @return ([real_size], [integer_size]) *)
        val get_work_space : 'a bsession -> int * int

        (** Returns the number of calls made to the user-supplied right-hand
            side function due to finite difference banded Jacobian approximation
            in the banded preconditioner setup function.

            @cvodes <node5#sss:cvbandpre> CVBandPrecGetNumRhsEvals *)
        val get_num_rhs_evals : 'a bsession -> int
      end

    (** {3 Output} *)

    (** Returns the real and integer workspace sizes.

        @cvodes <node5#sss:optout_main> CVodeGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space          : 'a bsession -> int * int

    (** Returns the cumulative number of internal steps taken by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetNumSteps *)
    val get_num_steps           : 'a bsession -> int

    (** Returns the number of calls to the user's right-hand side function.

        @cvodes <node5#sss:optout_main> CVodeGetNumRhsEvals *)
    val get_num_rhs_evals       : 'a bsession -> int

    (** Returns the number of calls made to the linear solver's setup function.

        @cvodes <node5#sss:optout_main> CVodeGetNumLinSolvSetups *)
    val get_num_lin_solv_setups : 'a bsession -> int

    (** Returns the number of local error test failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumErrTestFails *)
    val get_num_err_test_fails  : 'a bsession -> int

    (** Returns the integration method order used during the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastOrder *)
    val get_last_order          : 'a bsession -> int

    (** Returns the integration method order to be used on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentOrder *)
    val get_current_order       : 'a bsession -> int

    (** Returns the integration step size taken on the last internal step.

        @cvodes <node5#sss:optout_main> CVodeGetLastStep *)
    val get_last_step           : 'a bsession -> float

    (** Returns the integration step size to be attempted on the next internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentStep *)
    val get_current_step        : 'a bsession -> float

    (** Returns the the value of the integration step size used on the first
        step.

        @cvodes <node5#sss:optout_main> CVodeGetActualInitStep *)
    val get_actual_init_step    : 'a bsession -> float

    (** Returns the the current internal time reached by the solver.

        @cvodes <node5#sss:optout_main> CVodeGetCurrentTime *)
    val get_current_time        : 'a bsession -> float

    (** Returns the number of order reductions dictated by the BDF stability
        limit detection algorithm.

        @cvodes <node5#sss:optout_main> CVodeGetNumStabLimOrderReds
        @cvodes <node3#s:bdf_stab> BDF stability limit detection
     *)
    val get_num_stab_lim_order_reds : 'a bsession -> int

    (** Returns a suggested factor by which the user's tolerances should be
        scaled when too much accuracy has been requested for some internal
        step.

        @cvodes <node5#sss:optout_main> CVodeGetTolScaleFactor *)
    val get_tol_scale_factor : 'a bsession -> float

    (** Returns the solution error weights at the current time.

        @cvodes <node5#sss:optout_main> CVodeGetErrWeights
        @cvodes <node3#ss:ivp_sol> IVP solution (W_i) *)
    val get_err_weights : 'a bsession -> 'a nvector -> unit

    (** Returns the vector of estimated local errors.

        @cvodes <node5#sss:optout_main> CVodeGetEstLocalErrors *)
    val get_est_local_errors : 'a bsession -> 'a nvector -> unit

    (** Returns the integrator statistics as a group.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats *)
    val get_integrator_stats    : 'a bsession -> Cvode.integrator_stats

    (** Convenience function that calls get_integrator_stats and prints the
        results to stdout.

        @cvodes <node5#sss:optout_main> CVodeGetIntegratorStats *)
    val print_integrator_stats  : 'a bsession -> unit


    (** Returns the number of nonlinear (functional or Newton) iterations
        performed.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvIters *)
    val get_num_nonlin_solv_iters : 'a bsession -> int

    (** Returns the number of nonlinear convergence failures that have occurred.

        @cvodes <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails *)
    val get_num_nonlin_solv_conv_fails : 'a bsession -> int

    (** {3 Backward integration of quadrature equations} *)
    module Quadrature :
      sig
        (** {3:backquadinit Initialization} *)

        (** These functions compute the quadrature equation right-hand side for the
            backward problem.

            @cvodes <node7#ss:ODErhs_quad_b> CVQuadRhsFnB
            @cvodes <node7#ss:ODErhs_quad_sens_B> CVQuadRhsFnBS *)
        type 'a bquadrhsfn =
            QuadBasic of (float         (* t *)
                           -> 'a        (* y *)
                           -> 'a        (* yb *)
                           -> 'a        (* qbdot *)
                           -> unit)
          | QuadWithSens of (float           (* t *)
                              -> 'a          (* y *)
                              -> 'a array    (* ys *)
                              -> 'a          (* yb *)
                              -> 'a          (* qbdot *)
                              -> unit)

        (** This function, [init s fQB yQB0], activates integration of
            quadrature equations, with or without sensitivities, where [fQB]
            computes the right-hand side of the backward quadrature equations,
            and [yQB0] contains the values of the quadrature variables at [tB0].

            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitB
            @cvodes <node6#sss:cvquadinitb> CVodeQuadInitBS *)
        val init : 'a bsession -> 'a bquadrhsfn -> 'a nvector -> unit

        (** This function reinitializes the integration of quadrature equations
            during the backward phase.

            @cvodes <node6#ss:quad_sens_init> CVodeQuadReInitB *)
        val reinit : 'a bsession -> 'a nvector -> unit

        (** {4:adjextraction Extraction function} *)

        (**
          [tret = get s w yqs] fills [yqs] with the quadrature solution vector
          after a successful return from {!backward_normal} or
          {!backward_one_step}, and returns the time reached by the solver.

          @cvodes <node7#sss:quad_get_b> CVodeGetQuadB
         *)
        val get : 'a bsession -> 'a nvector -> float

        (** {4 Optional input functions} *)

        type 'a tolerance =
            NoStepSizeControl
            (** Do not use quadrature variables for step-size control (default). *)
          | SSTolerances of float * float
            (** [(rel, abs)] : scalar relative and absolute tolerances. *)
          | SVTolerances of float * 'a nvector
            (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

        (** Specify whether and how quadrature variables should be used in the step
            size control mechanism.

            @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
            @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances *)
        val set_tolerances : 'a bsession -> 'a tolerance -> unit

        (** {3 Optional output functions} *)

        (** Returns the number of calls to the user's quadrature right-hand side
            function.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output
                                                      functions for backward
                                                      quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals *)
        val get_num_rhs_evals       : 'a bsession -> int

        (** Returns the number of local error test failures due to quadrature
            variables.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output
                                                      functions for backward
                                                      quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails *)
        val get_num_err_test_fails  : 'a bsession -> int

        (** Returns the quadrature error weights at the current time.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output
                                                      functions for backward
                                                      quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights *)
        val get_err_weights : 'a bsession -> 'a nvector -> unit

        (** [nfqevals, nqetfails = get_stats s] returns
            - [fqevals], the number of calls to the user's quadrature function, and,
            - [nqetfails], the number of error test failures due to quadrature variables.

            @cvodes <node7#sss:quad_optional_input_B> Optional input/output
                                                      functions for backward
                                                      quadrature integration
            @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats *)
        val get_stats : 'a bsession -> int * int
      end
  end

