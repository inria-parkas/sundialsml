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
    val ss_tolerances : 'a session -> float -> float -> unit

    (**
        Specify that quadrature variables should be used in the step size
        control mechanism. [sv_tolerances s reltol abstol] sets the relative
        tolerance using a scalar value, and the absolute tolerance as a vector.

        @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances
     *)
    val sv_tolerances : 'a session -> float -> 'a vector -> unit

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
    val ss_tolerances : 'a session -> float -> float -> unit

    (**
        Specify the integration tolerances for sensitivities. [sv_tolerances s
        reltol abstol] sets the relative tolerance using a scalar value, and the
        absolute tolerance as a vector.

        @cvodes <node6#ss:cvfwdtolerances> CVodeSensSVtolerances
     *)
    val sv_tolerances : 'a session -> float -> 'a vector -> unit

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
        val ss_tolerances : 'a session -> float -> Sundials.real_array -> unit

        (**
            Specify that quadrature variables should be used in the step size
            control mechanism. [sv_tolerances s reltol abstol] sets the relative
            tolerance using a scalar value, and the absolute tolerance as an
            array of {i ns} vectors.

            @cvodes <node6#ss:quad_sens_optional_input> CVodeQuadSensSVtolerances
         *)
        val sv_tolerances : 'a session -> float -> 'a vector array -> unit

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

    (**
       This function evaluates the right-hand side of the backward ODE system.

       @cvodes <node7#ss:ODErhs_b> CVRhsFnB
       @cvodes <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis
       @cvodes <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis
     *)
    type 'a rhsfnb =
       float            (* t *)
         -> 'a          (* y *)
         -> 'a          (* yb *)
         -> 'a          (* ybdot *)
         -> unit

    (**
       This function evaluates the right-hand side of the backward ODE system
       when it depends on forward sensitivities.

       @cvodes <node7#ss:ODErhs_bs> CVRhsFnBS
       @cvodes <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis
       @cvodes <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis

    *)
    type 'a rhsfnbs =
       float            (* t *)
         -> 'a          (* y *)
         -> 'a array    (* ys *)
         -> 'a          (* yb *)
         -> 'a          (* ybdot *)
         -> unit

    (**
       This function computes the quadrature equation right-hand side for the
       backward problem.

       @cvodes <node7#ss:ODErhs_quad_b> CVQuadRhsFnB *)
    type 'a quadrhsfnb =
       float            (* t *)
         -> 'a          (* y *)
         -> 'a          (* yb *)
         -> 'a          (* qbdot *)
         -> unit

    (**
       This function computes the sensitivity-dependent quadrature equation
       right-hand side for the backward problem.

       @cvodes <node7#ss:ODErhs_quad_sens_B> CVQuadRhsFnBS *)
    type 'a quadrhsfnbs =
       float            (* t *)
         -> 'a          (* y *)
         -> 'a array    (* ys *)
         -> 'a          (* yb *)
         -> 'a          (* qbdot *)
         -> unit

    type 'a single_tmp = 'a
    type 'a triple_tmp = 'a * 'a * 'a

    (**
      Arguments common to all Jacobian callback functions.    
     
      @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB
      @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB 
      @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
      @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
      @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB
    *)
    type ('t, 'a) jacobian_arg =
      {
        jac_t   : float;        (** The independent variable. *)
        jac_y   : 'a;           (** The forward solution vector. *)
        jac_yb  : 'a;           (** The backward dependent variable vector. *)
        jac_fyb : 'a;           (** The backward right-hand side function fB. *)
        jac_tmp : 't            (** Workspace data,
                                    either {!single_tmp} or {!triple_tmp}. *)
      }

    (* TODO: consolidate all bandranges in Sundials ? *)
    (** The range of nonzero entries in a band matrix.  *)
    type bandrange = { mupper : int; (** The upper half-bandwidth.  *)
                       mlower : int; (** The lower half-bandwidth.  *) }

    (**
      This function computes the dense Jacobian of the backward problem (or an
      approximation to it).

      @cvodes <node7#ss:densejac_b> CVDlsDenseJacFnB
    *)
    type 'a dense_jac_fnb =
      ('a triple_tmp, 'a) jacobian_arg -> Dls.DenseMatrix.t -> unit

    (**
      This function computes the banded Jacobian of the backward problem (or an
      approximation to it).

      @cvodes <node7#ss:bandjac_b> CVDlsBandJacFnB 
    *)
    type 'a band_jac_fnb =
      band_range -> ('a triple_tmp, 'a) jacobian_arg -> Dls.DenseMatrix.t -> unit

    (**
      This function computes the action of the Jacobian for the backward problem
      on a given vector.

      @cvodes <node7#ss:jtimesv_b> CVSpilsJacTimesVecFnB
    *)
    type 'a jac_times_vec_fnb =
      ('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> unit

    (** Arguments passed to the preconditioner solve callback function.  See
        [prec_solve_fn] in {!spils_callbacks}.

        @cvode <node7#ss:psolveFn> CVSpilsPrecSolveFnB
     *)
    type 'a prec_solve_arg =
      {
        rvecB   : 'a;       (** The right-hand side vector, {i r}, of the
                                linear system. *)
        gammaB : float;     (** The scalar {i g} appearing in the Newton
                                matrix given by M = I - {i g}J. *)
        deltaB : float;     (** Input tolerance to be used if an
                                iterative method is employed in the
                                solution. *)
      }

    (**
      This function solves the preconditioning system {i Pz = r} for the
      backward problem.

      @cvodes <node7#ss:psolve_b> CVSpilsPrecSolveFnB
    *)
    type 'a prec_solve_fnb =
      ('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg -> 'a -> unit

    (**
      This function preprocesses and/or evaluates Jacobian-related data needed
      by the preconditioner for the backward problem.

      @cvodes <node7#ss:psetup_b> CVSpilsPrecSetupFnB
    *)
    type 'a prec_setup_fnb =
      ('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool

(* XXX WORKING XXX
 
  (**
     Specifies the type of interpolation.

      @cvodes <node3#ss:checkpointing> Checkpointing scheme
   *)
  type interpolation = IPolynomial (* CV_POLYNOMIAL *)
                     | IHermite    (* CV_HERMITE *)

   CVodeAdjInit
     Nd         : the number of integration steps between two consecutive checkpoints
     interptype : interpolation

      @cvodes <node7#ss:cvadjinit> CVodeAdjInit
    
  *)

    (** {3:adjforward Forward integration functions} *)

    (**
        [tret, ncheck = forward_normal s tout yret] integrates the forward
        problem over an interval and saves checkpointing data. The function
        takes as arguments the next time at which a solution is desired
        ([tout]), a vector for storing the computed result ([yret]), and returns
        the time reached by the solver ([tret]) and the number of checkpoints
        stored so far ([ncheck]).

        This call asks the solver to take internal steps until it has reached or
        just passed the [tout] parameter ([CV_NORMAL]). The solver then
        interpolates in order to return an approximate value of [y(tout)].

        @cvodes <node7#sss:cvsolvef> CVodeF
        TODO: list of exceptions raised.
     *)
    val forward_normal :
      'a session
      -> float
      -> 'a nvector
      -> float * int

    (**
        [tret, ncheck = forward_normal s tout yret] integrates the forward
        problem over an interval and saves checkpointing data. The function
        takes as arguments the next time at which a solution is desired
        ([tout]), a vector for storing the computed result ([yret]), and returns
        the time reached by the solver ([tret]) and the number of checkpoints
        stored so far ([ncheck]).

        This call asks the solver to take one internal step and to return the
        solution at the point reached by that step ([CV_ONE_STEP]).

        @cvodes <node7#sss:cvsolvef> CVodeF
        TODO: list of exceptions raised.
     *)
    val forward_one_step :
      'a session
      -> float
      -> 'a nvector
      -> float * int


(* XXX WORKING XXX
 
   CVodeCreateB
     ImmB       : multi-step method: CV_ADAMS or CV_BDF
     iterB      : nonlinear solver iteration: CV_NEWTON or CV_FUNCTIONAL.

    returns
     which      : an indentifier for the newly created backward problem.
     @cvodes <node6#sss:cvinitb> CVodeCreateB

   CVodeInitB
     which
     rhsB : CVRhsFnB
     tB0  : endpoint where final conditions are provided for the backward
            problem, normally equal to the endpoint of the forward integration
     yB0  : final value of the backward problem

   CVodeInitBS
     which
     rhsBS : CVRhsFnBS
     tB0   : endpoint where final conditions are provided for the backward
             problem, normally equal to the endpoint of the forward integration
     yB0   : final value of the backward problem

   CVodeReInitB
     which
     tB0
     yB0

   Other functions:
   CVodeSetUserDataB(cvode_mem, which, user_dataB)
   CVDlsSetDenseJacFnB
   CVDlsSetBandJacFnB
   CVSpilsSetPreconditionerB
   CVSpilsSetJacTimesVecFnB
   CVSpilsSetGSTypeB
   CVSpilsSetMaxlB
   CVSpilsSetEpsLinB
   CVSpilsSetPrecTypeB

   CVBandPrecInitB
   
   (* TODO: linear solver initialization functions...
            adapt from CVODE. *)

  *)

    (* TODO: understand how this works... *)
    (** Identifies a backward problem. *)
    type which

    (** {3 Tolerance specification} *)

    (**
        Specify the integration tolerances for the backward problem.
        [ss_tolerances s reltol abstol] sets the relative and absolute
        tolerances using scalar values.

        @cvodes <node7#sss:cvtolerances_b> CVodeSStolerancesB
     *)
    val ss_tolerances : 'a session -> which -> float -> float -> unit

    (**
        Specify the integration tolerances for the backward problem.
        [sv_tolerances s reltol abstol] sets the relative tolerance using a
        scalar value, and the absolute tolerance as a vector.

        @cvodes <node7#sss:cvtolerances_b> CVodeSVtolerancesB
     *)
    val sv_tolerances : 'a session -> which -> float -> 'a vector -> unit

    (** {3:adjbackward Backward integration functions} *)

    (**
        [backward_normal s tbout] integrates the backward ODE problem. The
        function takes internal steps until it has reached or just passed the
        user-specified value [tbout] ([CV_NORMAL]). The solver then interpolates
        in order to return an approximate value of [y(tbout)].

        @cvodes <node7#sss:cvsolveb> CVodeB
        TODO: list of exceptions raised.
     *)
    val backward_normal : 'a session -> float -> unit

    (**
        [backward_one_step s tbout] integrates the backward ODE problem. The
        function takes one internal step ([CV_ONE_STEP]).

        @cvodes <node7#sss:cvsolveb> CVodeB
        TODO: list of exceptions raised.
     *)
    val backward_one_step : 'a session -> float -> unit

    (**
        [tret = get s which yb] returns the solution of the backward ODE problem
        in [yb] at time [tret].

        @cvodes <node7#sss:cvsolveb> CVodeGetB
     *)
    val get : 'a session -> which -> 'a nvector -> float

    (** {3 Optional input functions} *)

    (**
        Instructs {!forward_normal} and {!forward_one_step} not to save
        checkpointing data for forward sensitivities anymore.

        @cvodes <node7#SECTION00727000000000000000> CVodeAdjSetNoSensi
     *)
    val set_no_sensitivity : 'a session -> unit

    (**
      Specifies the maximum order of the linear multistep method.

      @cvodes <node7#ss:optional_input_b> CVodeSetMaxOrdB
     *)
    val set_max_ord : 'a session -> which -> int -> unit

    (**
      Specifies the maximum number of steps to be taken by the solver in its attempt
      to reach the next output time.

      @cvodes <node7#ss:optional_input_b> CVodeSetMaxNumStepsB
     *)
    val set_max_num_steps : 'a session -> which -> int -> unit

    (**
      Specifies the initial step size.

      @cvodes <node7#ss:optional_input_b> CVodeSetInitStepB
     *)
    val set_init_step : 'a session -> which -> float -> unit

    (**
      Specifies a lower bound on the magnitude of the step size.

      @cvodes <node7#ss:optional_input_b> CVodeSetMinStepB
     *)
    val set_min_step : 'a session -> which -> float -> unit

    (**
      Specifies an upper bound on the magnitude of the step size.

      @cvodes <node7#ss:optional_input_b> CVodeSetMaxStepB
     *)
    val set_max_step : 'a session -> which -> float -> unit

    (**
      Indicates whether the BDF stability limit detection algorithm should be
      used.

      @cvode <node7#ss:optional_input_b> CVodeSetStabLimDet
     *)
    val set_stab_lim_det : 'a session -> which -> bool -> unit

    (** {3 Optional output functions} *)

    (* TODO:
      
        cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, which)

        gets the cvode_mem associated with 'which'. The CVodeGet* and CVode*Get*
        functions can be called on this. "The user should not modify in any way
        cvode_memB"...

        How do we handle this?
     *)

    (** {3 Integration of quadrature equations depending on forward sensitiviites} *)
    module Quadrature :
      sig
    (* XXX WORKING XXX
        
        CVodeQuadInitB:
            which
            rhsQB : CVQuadRhsFnB
            yQB0  : 'a nvector

        CVodeQuadInitBS
            which
            rhsQBS : CVQuadRhsFnBS
            yQBS0  : 'a nvector

        CVodeQuadReInitB
            which
            yQB0   : 'a nvector
     *)
        (** {4:adjextraction Extraction function} *)

        (**
          [tret = get s w yqs] fills [yqs] with the quadrature solution vector
          after a successful return from {!backward_normal} or
          {!backward_one_step}, and returns the time reached by the solver.

          @cvodes <node7#sss:quad_get_b> CVodeGetQuadB
         *)
        val get : 'a session -> 'a nvector array -> float

        (** {4 Tolerance specification} *)

        (**
            Set whether quadrature variables should be used in the step size
            mechanism (the default is [true]).

            @cvodes <node7#sss:quad_optional_input_B> CVodeSetQuadErrConB
         *)
        val set_err_con : 'a session -> bool -> unit

        (**
            Specify that quadrature variables should be used in the step size
            control mechanism. [ss_tolerances s reltol abstol] sets the relative and
            absolute tolerances using scalar values.

            @cvodes <node7#sss:quad_optional_input_B> CVodeQuadSStolerancesB
         *)
        val ss_tolerances : 'a session -> which -> float -> float -> unit

        (**
            Specify that quadrature variables should be used in the step size
            control mechanism. [sv_tolerances s reltol abstol] sets the relative
            tolerance using a scalar value, and the absolute tolerances from a
            vector.

            @cvodes <node7#sss:quad_optional_input_B> CVodeQuadSVtolerancesB
         *)
        val sv_tolerances : 'a session -> which -> float -> 'a vector -> unit

        (* TODO: use CVodeGetAdjCVodeBmem for indirect access via the
           CVodeGetQuad* functions. *)
      end

  end

