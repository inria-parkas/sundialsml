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

(** Abstract nvector interface to the CVODES Extensions.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

(** {2:quadrature Integration of pure quadrature equations} *)

module Quadrature :
  sig

(* XXX WORKING XXX

   CVodeQuadInit
     fQ   : float -> 'a -> 'a -> unit (CVQuadRhsFn)
     yQ0  : 'a nvector

      @cvodes <node5#ss:quad_malloc> CVodeQuadInit
      @cvodes <node5#ss:user_fct_quad> CVQuadRhsFn
    

   CVodeQuadReInit
     yQ0  : 'a nvector
     (check for CV_NO_QUAD if cvode_mem not QuadInited... )
      @cvodes <node5#ss:quad_malloc> CVodeQuadReInit

  Accept tolerances during initialization?

  Design questions:
    - How to incorporate extra elements into the session type.
      (i.e., qrhs function )
    - How to mark the session type as valid for the quadrature functions?
      1. Fail dynamically (add an exception for CV_NO_QUAD).
      2. Add a phantom type.
      3. create a quad_session from which a normal session can be extracted?
      4. upgrade a normal session to a quad_session?
 *)
    (** {3 Tolerance specification} *)

    (* TODO: necessitates call to one of the *_tolerances functions... *)
    (**
        Set whether quadrature variables should be used in the error control
        mechanism (the default is [true]). 

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
     *)
    val set_err_con : 'a session -> bool -> unit

    (**
        Specify that quadrature variables should be used in the step size
        control mechanism. [ss_tolerances s reltol abstol] sets the relative and
        absolute tolerances using scalar values.

        @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
     *)
    val ss_tolerances : session -> float -> float -> unit

    (**
        Specify that quadrature variables should be used in the step size
        control mechanism. [sv_tolerances s reltol abstol] sets the relative
        tolerance using a scalar value, and the absolute tolerance as a vector.

        @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances
     *)
    val sv_tolerances : session -> float -> 'a vector -> unit

    (** {3:exceptions Exceptions} *)

    (** @cvodes <node5#SECTION00572000000000000000> CV_QRHSFUNC_FAIL *)
    exception QuadRhsFuncFailure

    (** @cvodes <node5#SECTION00572000000000000000> CV_FIRST_QRHSFUNC_ERR *)
    exception FirstQuadRhsFuncErr

    (** @cvodes <node5#SECTION00572000000000000000> CV_REPTD_QRHSFUNC_ERR *)
    exception RepeatedQuadRhsFuncErr

    (** @cvodes <node5#SECTION00572000000000000000> CV_UNREC_QRHSFUNC_ERR *)
    exception UnrecoverableQuadRhsFuncErr

    (** {3:output Output functions} *)

    (**
      [tret = get s yq] fills [yq] with the quadrature solution vector after a
      successful return from {!Cvode.solve_normal} or {!Cvode.solve_one_step},
      and returns the time reached by the solver.

      @cvodes <node5#ss:quad_get> CVodeGetQuad
     *)
    val get : 'a session -> 'a nvector -> float

    (**
      [tret = get_dky s t k dkyq] fills [dkyq] with the derivatives of the
      quadrature solution vector after a successful return from
      {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested, [t],
      must fall within the interval defined by the last successful step
      ({!Cvode.get_last_step}). The requested order, [k], must be less than or
      equal to the value returned by {!Cvode.get_last_order}.

      @cvodes <node5#ss:quad_get> CVodeGetQuadDky
     *)
    val get_dky : 'a session -> float -> int -> 'a nvector -> unit


    (** {3 Optional output functions} *)

    (**
      Returns the number of calls to the user's quadrature right-hand side function.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals
     *)
    val get_num_rhs_evals       : 'a session -> int

    (**
      Returns the number of local error test failures due to quadrature variables.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails
     *)
    val get_num_err_test_fails  : 'a session -> int

    (**
      Returns the quadrature error weights at the current time.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights
     *)
    val get_quad_err_weights : 'a session -> 'a nvector -> unit

    (**
      [nfqevals, nqetfails = get_stats s] returns
      - [fqevals], the number of calls to the user's quadrature function, and,
      - [nqetfails], the number of error test failures due to quadrature variables.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats
     *)
    val get_stats : session -> int * int

  end

(** {2:forward Forward Sensitivity Analysis} *)

module Forward :
  sig

(* XXX WORKING XXX
 
  type method = Simultaneous (* CV_SIMULTANEOUS *)
              | Staggered    (* CV_STAGGERED *)
              | Staggered1   (* CV_STAGGERED1 *)

  (** @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn *)
  type 'a rhsfn =
     float            (* t *)
       -> 'a          (* y *)
       -> 'a          (* ydot *)
       -> 'a array    (* yS *)
       -> 'a array    (* ySdot *)
       -> 'a          (* tmp1 *)
       -> 'a          (* tmp2 *)
       -> unit

  (** @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn *)
  type 'a rhs1fn =
     float            (* t *)
       -> 'a          (* y *)
       -> 'a          (* ydot *)
       -> int         (* iS *)
       -> 'a          (* yS *)
       -> 'a          (* ySdot *)
       -> 'a          (* tmp1 *)
       -> 'a          (* tmp2 *)
       -> unit

   CVodeSensInit
     m    : method      (* Staggered1 not allowed *)
     fS   : 'a rhsfn
     ys0  : 'a nvector array (get ns from here)

      @cvodes <node6#ss:sensi_malloc> CVodeSensInit
      @cvodes <node6#ss:user_fct_fwd> CVSensRhsFn

   CVodeSensInit1
     m    : method
     fS   : 'a rhs1fn
     ys0  : 'a nvector array (get ns from here)

      @cvodes <node6#ss:sensi_malloc> CVodeSensInit1
      @cvodes <node6#ss:user_fct_fwd> CVSensRhs1Fn
    

   CVodeSensReInit
     m    : method      (* may not always pass Staggered1 *)
     ys0  : 'a nvector array
     (check for CV_NO_SENS if cvode_mem not SensInited... )
      @cvodes <node6#ss:sensi_malloc> CVodeSensReInit
  *)

    (**
     Deactivates forward sensitivity calculations without deallocating memory.
     Sensitivities can be reactivated with {!reinit}.

     @cvodes <node6#ss:sensi_malloc> CVodeSensToggleOff
     *)
    val toggle_off : 'a session -> unit

    (** {3 Tolerance specification} *)

    (**
        Specify the integration tolerances for sensitivities. [ss_tolerances s
        reltol abstol] sets the relative and absolute tolerances using scalar
        values.

        @cvodes <node6#sss:cvfwdtolerances> CVodeSensSStolerances
     *)
    val ss_tolerances : session -> float -> float -> unit

    (**
        Specify the integration tolerances for sensitivities. [sv_tolerances s
        reltol abstol] sets the relative tolerance using a scalar value, and the
        absolute tolerance as a vector.

        @cvodes <node6#ss:cvfwdtolerances> CVodeSensSVtolerances
     *)
    val sv_tolerances : session -> float -> 'a vector -> unit

    (**
        Specify the integration tolerances for sensitivities based on those
        supplied for state variables and the scaling factors.

        @cvodes <node6#ss:cvfwdtolerances> CVodeSensEEtolerances
     *)
    val ee_tolerances : unit -> unit

    (** {3:exceptions Exceptions} *)

    (** @cvodes <node6#SECTION00623000000000000000> CV_SRHSFUNC_FAIL *)
    exception SensRhsFuncFailure

    (** @cvodes <node6#SECTION00623000000000000000> CV_FIRST_SRHSFUNC_ERR *)
    exception FirstSensRhsFuncErr

    (** @cvodes <node6#SECTION00623000000000000000> CV_REPTD_SRHSFUNC_ERR *)
    exception RepeatedSensRhsFuncErr

    (** @cvodes <node6#SECTION00623000000000000000> CV_UNREC_SRHSFUNC_ERR *)
    exception UnrecoverableSensRhsFuncErr

    (** {3:output Output functions} *)

    (**
      [tret = get s ys] fills [ys] with the sensitivity solution vectors after a
      successful return from {!Cvode.solve_normal} or {!Cvode.solve_one_step},
      and returns the time reached by the solver.

      @cvodes <node6#ss:sensi_get> CVodeGetSens
     *)
    val get : 'a session -> 'a nvector array -> float

    (**
      [tret = get_dky s t k dkys] fills [dkys] with the derivatives of the
      sensitivity solution vectors after a successful return from
      {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested, [t],
      must fall within the interval defined by the last successful step
      ({!Cvode.get_last_step}). The requested order, [k], must be less than or
      equal to the value returned by {!Cvode.get_last_order}.

      @cvodes <node6#ss:sensi_get> CVodeGetSensDky
     *)
    val get_dky : 'a session -> float -> int -> 'a nvector array -> unit

    (**
      [tret = get s i ys] fills [ys] with the [i]th sensitivity solution vector
      after a successful return from {!Cvode.solve_normal} or
      {!Cvode.solve_one_step}, and returns the time reached by the solver.

      @cvodes <node6#ss:sensi_get> CVodeGetSens1
     *)
    val get' : 'a session -> int -> 'a nvector -> float

    (**
      [tret = get_dky s t k i dkys] fills [dkys] with the derivatives of the
      [i]th sensitivity solution vector after a successful return from
      {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested, [t],
      must fall within the interval defined by the last successful step
      ({!Cvode.get_last_step}). The requested order, [k], must be less than or
      equal to the value returned by {!Cvode.get_last_order}.

      @cvodes <node6#ss:sensi_get> CVodeGetSensDky1
     *)
    val get_dky' : 'a session -> float -> int -> int -> 'a nvector -> unit

    (** {3 Optional input functions} *)

    (* TODO: check that a lint_array is compatabile with int* *)
    (* TODO:  If non-NULL, p must point to a field in the user's data structure
              user_data passed to the right-hand side function. (See §5.1). *)
    (* TODO: check that pbar and plist are ns long. *)
    (**
      [set_params s p pbar plist] specifies problem parameter information for
      sensitivity calculations:
        - [p], the parameters used to evaluate {i f(t, y, p)},
        - [pbar], an array of {i ns} positive scaling factors, and,
        - [plist], an array of non-negative indices to specify which components
        to use in estimating the sensitivity equations.

      @cvodes <node6#ss:sens_optional_input> CVodeSetSensParams
     *)
    val set_params : 'a session -> Sundials.real_array option
                     -> Sundials.real_array option
                     -> Sundials.lint_array option
                     -> unit

    type dq_method = DQCentered (* CV_CENTERED *)
                   | DQForward  (* CV_FORWARD *)

    (**
      [set_dq_method s dqtype dqrhomax] specifies the difference quotient
      strategy in the case in which the right-hand side of the sensitivity
      equations is to be computed by CVODES; [dqrhomax] is used in deciding the
      switching between simultaneous or separate approximations of the two tersm
      in the sensitivity right-hand side.

      @cvodes <node6#ss:sens_optional_input> CVodeSetSensDQMethod
      @cvodes <node3#ss:fwd_sensi> Forward Sensitivity Analysis
     *)
    val set_dq_method : 'a session -> dq_method -> float -> unit

    (**
        Set whether sensitivity variables should be used in the error control
        mechanism (the default is [true]). 

        @cvodes <node5#ss:sens_optional_input> CVodeSetSensErrCon
     *)
    val set_err_con : 'a session -> bool -> unit

    (**
      Specifies the maximum number of nonlinear solver iterations for
      sensitivity variables permitted per step.

      @cvode <node5#ss:sens_optional_input> CVodeSetSensMaxNonlinIters
     *)
    val set_max_nonlin_iters : 'a session -> int -> unit

    (** {3 Optional output functions} *)

    (**
      Returns the number of calls to the sensitivity right-hand side function.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensNumRhsEvals
     *)
    val get_num_sens_evals       : 'a session -> int

    (**
      Returns the number of calls to the user's right-hand side function due to
      the internal finite difference approximation of the sensitivity
      right-hand sides.

      @cvode <node6#ss:sens_optional_output> CVodeGetNumRhsEvalsSens
     *)
    val get_num_rhs_evals      : 'a session -> int

    (**
      Returns the number of local error test failures for the sensitivity
      variables that have occurred.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensNumErrTestFails
     *)
    val get_num_err_test_fails  : 'a session -> int

    (**
      Returns the number of calls made to the linear solver's setup function due
      to forward sensitivity calculations.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensNumLinSolvSetups
     *)
    val get_num_lin_solv_setups : 'a session -> int

    type sensitivity_stats = {
        num_rhs_evals : int;
        num_sens_evals :int
        num_err_test_fails : int;
        num_lin_solv_setups :int;
      }

    (**
      Returns all of the sensitivity-related solver statistics as a group.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensStats
     *)
    val get_stats : 'a session -> sensitivity_stats

    (**
      Returns the sensitivity error weight vectors at the current time.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensErrWeights
      @cvode <node3#e:errwt> Eq. (2.7) IVP solution (W_i)
     *)
    val get_err_weights : 'a session -> 'a nvector array -> unit

    (**
      Returns the number of nonlinear iterations performed for sensitivity
      calculations.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvIters
     *)
    val get_num_nonlin_solv_iters : 'a session -> int

    (**
      Returns the number of nonlinear convergence failures that have occurred
      for sensitivity calculations.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensNumNonlinSolvConvFails
     *)
    val get_num_nonlin_solv_conv_fails : 'a session -> int

    type nonlin_stats = {
        num_nonlin_solv_iters : int;
        num_nonlin_solv_conv_fails : int;
      }

    (**
      Returns the sensitivity-related nonlinear solver statistics as a group.

      @cvode <node6#ss:sens_optional_output> CVodeGetSensNonlinSolvStats
     *)
    val get_nonlin_solv_stats : 'a session -> nonlin_stats

    (**
      Returns the number of nonlinear (functional or Newton) iterations
      performed for each sensitivity equation separately, in the {!Staggered1}
      case.

      @cvode <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvIters
     *)
    val get_num_stgr_nonlin_solv_iters : 'a session
                                         -> Sundials.lint_array -> unit

    (**
      Returns the number of nonlinear convergence failures that have occurred
      for each sensitivity equation separately, in the {!Staggered1} case.

      @cvode <node6#ss:sens_optional_output> CVodeGetStgrSensNumNonlinSolvConvFails
     *)
    val get_num_stgr_nonlin_solv_conv_fails : 'a session
                                              -> Sundials.lint_array -> unit

    (** {3 Integration of quadrature equations depending on forward sensitiviites} *)
    module Quadrature :
      sig
    (* XXX WORKING XXX

       (** @cvodes <node6#ss:user_fct_quad_sens> CVodeQuadSensRhsFn *)
       type 'a quadsensrhsfn =
          float
          -> 'a nvector
          -> 'a nvector
          -> 'a nvector
          -> 'a nvector array
          -> 'a nvector
          -> 'a nvector
          -> unit

       TODO: must call CVodeSensInit or CVodeSensInit1 first ...
       CVodeQuadSensInit
         fQ   : 'a quadsensrhsfn
         yqs0  : 'a nvector array

          @cvodes <node6#ss:quad_sens_init> CVodeQuadSensInit

       CVodeQuadSensReInit
         yqs0  : 'a nvector array
         (check for CV_NO_QUAD if cvode_mem not QuadInited... )
          @cvodes <node6#ss:quad_sens_init> CVodeQuadSensReInit

      Accept tolerances during initialization?

     *)
        (** {4 Tolerance specification} *)

        (* TODO: necessitates call to one of the *_tolerances functions... *)
        (**
            Set whether quadrature variables should be used in the step size
            mechanism (the default is [true]).

            @cvodes <node6#ss:quad_sens_optional_input> CVodeSetQuadSensErrCon
         *)
        val set_err_con : 'a session -> bool -> unit

        (**
            Specify that quadrature variables should be used in the step size
            control mechanism. [ss_tolerances s reltol abstol] sets the relative and
            absolute tolerances using a scalar value and an array of {i ns} scalar
            values.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSStolerances
         *)
        val ss_tolerances : session -> float -> Sundials.real_array -> unit

        (**
            Specify that quadrature variables should be used in the step size
            control mechanism. [sv_tolerances s reltol abstol] sets the relative
            tolerance using a scalar value, and the absolute tolerance as an
            array of {i ns} vectors.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSVtolerances
         *)
        val sv_tolerances : session -> float -> 'a vector array -> unit

        (* TODO: Requires that CVodeQuad*Tolerances be set... *)
        (**
            Specify that the tolerances for sensitivity-dependent quadratures be
            estimated from those provided for pure quadrature variables.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensEEtolerances
         *)
        val ee_tolerances : unit -> unit

        (** {4:exceptions Exceptions} *)

        (** @cvodes <node6#SECTION00642000000000000000> CV_QSRHSFUNC_FAIL *)
        exception QuadSensRhsFuncFailure

        (** @cvodes <node6#SECTION00642000000000000000> CV_FIRST_QSRHSFUNC_ERR *)
        exception FirstQuadSensRhsFuncErr

        (** @cvodes <node6#SECTION00642000000000000000> CV_REPTD_QSRHSFUNC_ERR *)
        exception RepeatedQuadSensRhsFuncErr

        (** @cvodes <node6#SECTION00642000000000000000> CV_UNREC_QSRHSFUNC_ERR *)
        exception UnrecoverableQuadSensRhsFuncErr

        (** {4:extraction Extraction functions} *)

        (**
          [tret = get s yqs] fills [yqs] with quadrature solution vectors after a
          successful return from {!Cvode.solve_normal} or {!Cvode.solve_one_step},
          and returns the time reached by the solver.

          @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSens
         *)
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

        (**
          [tret = get_dky s t k i dkyqs] fills [dkyqs] with the derivatives of
          the [i]th quadrature solution vector after a successful return from
          {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested,
          [t], must fall within the interval defined by the last successful step
          ({!Cvode.get_last_step}). The requested order, [k], must be less than
          or equal to the value returned by {!Cvode.get_last_order}.

          @cvodes <node6#ss:quad_sens_get> CVodeGetQuadSensDky1
         *)
        val get_dky1 : 'a session -> float -> int -> int -> 'a nvector array -> unit

        (** {4 Optional output functions} *)

        (**
          Returns the number of calls to the user's quadrature right-hand side function.

          @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumRhsEvals
         *)
        val get_num_rhs_evals       : 'a session -> int

        (**
          Returns the number of local error test failures due to quadrature variables.

          @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensNumErrTestFails
         *)
        val get_num_err_test_fails  : 'a session -> int

        (**
          Returns the quadrature error weights at the current time.

          @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensErrWeights
         *)
        val get_quad_err_weights : 'a session -> 'a nvector array -> unit

        (**
          [nfqevals, nqetfails = get_stats s] returns
          - [fqevals], the number of calls to the user's quadrature function, and,
          - [nqetfails], the number of error test failures due to quadrature variables.

          @cvodes <node6#ss:quad_sens_optional_output> CVodeGetQuadSensStats
         *)
        val get_stats : session -> int * int

      end
    end

(** {2:adjoint Adjoint Sensitivity Analysis} *)

module Adjoint :
  sig
  end

