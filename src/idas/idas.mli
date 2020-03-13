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
(*               User Documentation for IDAS v1.1.0                    *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Sensitivity analysis (forward and adjoint) and quadrature equations for
    IDA.

    These submodules extend basic state equations,
    {% $F(t, y, \dot{y}) = 0$%},
    with parameters $p$ and quadrature variables $y_Q$. The former exposes
    a new dependency, {% $F(t, y, \dot{y}, p) = 0$%}, relative to which the
    sensitivity of a solution can be approximated. The latter introduces
    an efficient means to calculate a subset of variables defined by ODEs:
    {% $\dot{y}_Q = f_Q(t, y, \dot{y}, p)$%}.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** Alias for Ida sessions. The Idas submodules all work by ‘extending’
    an existing session created with {!Ida.init}. *)
type ('data, 'kind) session = ('data, 'kind) Ida.session

(** Integration of pure quadrature equations.

    Adds a vector $y_Q$ of $N_Q$ quadrature variables defined by
    {% $\frac{\mathrm{d} y_Q}{\mathrm{d}t} = f_Q(t, y, \dot{y}, p)$%}. These
    are treated more efficiently since they are excluded from the nonlinear
    solution stage.

    @idas <node5#SECTION00570000000000000000> Integration of pure quadrature equations
    @idas <node3#s:quad> Pure quadrature integration *)
module Quadrature : sig (* {{{ *)
  (** {2:init Initialization and access} *)

  (** Functions defining quadrature variables. They are passed four
      arguments:
      - [t], the value of the independent variable, i.e., the simulation time,
      - [y], the vector of dependent-variable values, i.e., $y(t)$,
      - [y'], the vector of dependent-variable derivatives,
                i.e., {% $\dot{y}(t)$%} and,
      - [yq'], a vector for storing the computed value of
               {% $\dot{y}_Q = f_Q(t, y, \dot{y}, p)$%}.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning The vectors in the function arguments should not
               be accessed after the function returns.}

      @idas <node5#ss:user_fct_quad> IDAQuadRhsFn *)
  type 'd quadrhsfn = float -> 'd -> 'd -> 'd -> unit

  (** Activates the integration of quadrature equations. The vector
      gives the initial value of $y_Q$.

      @idas <node5#ss:quad_init> IDAQuadInit *)
  val init : ('d, 'k) Ida.session
              -> 'd quadrhsfn -> ('d, 'k) Nvector.t -> unit

  (** Reinitializes the integration of quadrature equations. The vector
      gives a new value for $y_Q$.

      @idas <node5#ss:quad_init> IDAQuadReInit *)
  val reinit : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t -> unit

  (** Returns the quadrature solutions and time reached after a successful
      solver step. The given vector is filled with values calculated during
      either {!Ida.solve_normal} or {!Ida.solve_one_step} and the
      value of the independent variable is returned.

      @idas <node5#ss:quad_get> IdaGetQuad *)
  val get : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t -> float

  (** Returns the interpolated solution or derivatives of quadrature
      variables.

      [get_dky s dkyq t k] computes the [k]th derivative at time [t], i.e.,
      {% $\frac{d^\mathtt{k}y_Q(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
      and stores it in [dkyq]. The arguments must satisfy {% $t_n - h_u \leq
      \mathtt{t} \leq t_n$%}—where $t_n$ denotes {!Ida.get_current_time}
      and $h_u$ denotes {!Ida.get_last_step},— and
      {% $0 \leq \mathtt{k} \leq q_u$%}—where
      $q_u$ denotes {!Ida.get_last_order}.

      This function may only be called after a successful return from either
      {!Ida.solve_normal} or {!Ida.solve_one_step}.

      @idas <node5#ss:quad_get> IdaGetQuadDky
      @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
      @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
  val get_dky : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t
                  -> float -> int -> unit

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

      @idas <node5#ss:quad_optional_input> IDASetQuadErrCon
      @idas <node5#ss:quad_optional_input> IDAQuadSStolerances
      @idas <node5#ss:quad_optional_input> IDAQuadSVtolerances *)
  val set_tolerances : ('d, 'k) Ida.session -> ('d, 'k) tolerance -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the number of calls to the quadrature function.

      @idas <node5#ss:quad_optional_output> IDAGetQuadNumRhsEvals *)
  val get_num_rhs_evals : ('d, 'k) Ida.session -> int

  (** Returns the number of local error test failures that have occurred
      due to quadrature variables.

      @idas <node5#ss:quad_optional_output> IDAGetQuadNumErrTestFails *)
  val get_num_err_test_fails : ('d, 'k) Ida.session -> int

  (** Returns the quadrature error weights at the current time.

      @idas <node5#ss:quad_optional_output> IDAGetQuadErrWeights *)
  val get_err_weights : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t -> unit

  (** Returns quadrature-related statistics. These are the
      number of calls to the quadrature function ([nfqevals]) and the number
      of error test failures due to quadrature variables ([nqetfails]).

      @idas <node5#ss:quad_optional_output> IDAGetQuadStats
      @return ([nfqevals], [nqetfails]) *)
  val get_stats : ('d, 'k) Ida.session -> int * int

  (** {2:exceptions Exceptions} *)

  (** Quadrature integration was not initialized.

      @idas <node5#ss:quad_get> IDA_NO_QUAD *)
  exception QuadNotInitialized

  (** The quadrature function failed in an unrecoverable manner.

      @idas <node5#SECTION00572000000000000000> IDA_QRHS_FAIL *)
  exception QuadRhsFuncFailure

  (** The quadrature function failed at the first call.

      @idas <node5#SECTION00572000000000000000> IDA_FIRST_QRHS_ERR *)
  exception FirstQuadRhsFuncFailure

  (** Convergence test failures occurred too many times due to repeated
      recoverable errors in the quadrature function. Also raised if the
      quadrature function had repeated recoverable errors during the
      estimation of an initial step size (if quadrature variables are
      included in error tests).

      @idas <node5#SECTION00572000000000000000> IDA_REP_QRHS_ERR *)
  exception RepeatedQuadRhsFuncFailure
end (* }}} *)

(** (Forward) Sensitivity analysis of DAEs with respect to their parameters.

    Formalizes the dependency of a set of DAEs on $N_p$ parameters $p$ and
    calculates the sensitivities $s$ and their derivatives {% $\dot{s}$%}
    of the solutions $y$ and {% $\dot{y}$%} relative to a subset of
    {% $N_s \leq N_p$%} of those parameters such that the $i$th sensitivity
    satisfies the {i (forward) sensitivity equations}:
    {% $\frac{\partial F}{\partial y}s_i(t) + \frac{\partial F}{\partial \dot{y}}\dot{s}_i(t) + \frac{\partial F}{\partial p_i} = 0$%},
    {% $s_i(t_0) = \frac{\partial y_0(p)}{\partial p_i}$%}, and,
    {% $\dot{s_i}(t_0) = \frac{\partial \dot{y}_0(p)}{\partial p_i}$%}, where
    $F$, $y$, and {% $\dot{y}$%} are from the $N$ equations of the original
    model.

    This documented interface is structured as follows.
    {ol
      {- {{:#init}Initialization} (including {{!Quadrature}Quadrature equations})}
      {- {{:#sensout}Output functions}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#exceptions}Exceptions}}}

    @idas <node6#SECTION00610000000000000000> Enhanced skeleton for sensitivity analysis
    @idas <node6#s:forward> Using IDAS for Forward Sensitivity Analysis
    @idas <node3#ss:adj_sensi> Adjoint sensitivity analysis *)
module Sensitivity : sig (* {{{ *)
  (** {2:init Initialization} *)

  (** Arguments to {!sensresfn}. *)
  type 'd sensresfn_args =
    {
      t : float;
      (** The value of the independent variable. *)

      y : 'd;
      (** The vector of dependent-variable values $y(t)$. *)

      y' : 'd;
      (** The vector of derivatives of dependent variables
          {% $\dot{y}(t)$%}. *)

      res : 'd;
      (** The vector of current residuals. *)

      s : 'd array;
      (** The array of sensitivity vectors. *)

      s' : 'd array;
      (** The array of sensitivity derivative vectors. *)

      tmp : 'd Ida.triple
      (** Temporary storage vectors. *)
    }

  (** Sensitivity functions that calculate the residuals of all
      sensitivity equations. They are passed the arguments:
      - [args], which summarizes current values of state and sensitivty
                variables, and
      - [rs] an array of vectors for storing the calculated sensitivity
             residual values,
             {% $\mathit{rs}_i = \frac{\partial F}{\partial y}s_i(t) + \frac{\partial F}{\partial \dot{y}}\dot{s}_i(t) + \frac{\partial F}{\partial p_i}$ %}

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning The vectors in the function arguments should not
               be accessed after the function returns.}

      @idas <node6#s:user_fct_fwd> IDASensResFn *)
  type 'd sensresfn = 'd sensresfn_args -> 'd array -> unit

  (** Specifies a sensitivity solution method.
      @idas <node6#ss:sensi_init> IDASensInit *)
  type sens_method =
      Simultaneous
      (** Correct state and sensitivity variables at the same time.
          {cconst IDA_SIMULTANEOUS} *)
    | Staggered
      (** Correct sensitivity variables simultaneously, but only
          after the convergence of state variable corrections and
          when state variables pass the local error test.
          {cconst IDA_STAGGERED} *)

  (** Specifies problem parameter information for sensitivity
      calculations.  Which fields are required varies as follows:
      - [pvals] is mandatory if {!Sensitivity.init} and/or
        {!Sensitivity.Quadrature.init} is not given a
        user-defined callback.  Otherwise, it is ignored.
      - [plist] is optional if [pvals] is given.  Otherwise, it is
        ignored.
      - [pbar] is always meaningful and optional.

      @idas <node6#ss:sens_optional_input> IDASetSensParams *)
  type sens_params = {
    pvals : RealArray.t option;
      (** An array of $N_p$ parameters $p$ in {% $F(t, y, \dot{y},
          p)$%}.  If specified, this array is updated by the solver
          to pass perturbed parameter values to the original
          residual and root functions.  Those functions must
          (re-)read parameter values from this array on every
          invocation. *)
    pbar : RealArray.t option;
      (** An array of $N_s$ positive scaling factors. These are needed
          when estimating tolerances or using the internal difference
          quotient function. *)
    plist : int array option;
      (** An array of $N_s$ ($< N_p$) indices specifying the parameters
          to analyze. If not specified, the first $N_s$ parameters are
          analyzed. *)
  }

  (** Tolerances for calculating sensitivities. *)
  type ('d, 'k) tolerance =
      SStolerances of float * RealArray.t
      (** [(rel, abs)] : scalar relative and absolute tolerances. *)
    | SVtolerances of float * ('d, 'k) Nvector.t array
      (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
    | EEtolerances
      (** Calculate the integration tolerances for sensitivities
          from those of state variables and the scaling factors
          in {{!sens_params}pbar}. *)

  (** Activates the calculation of forward sensitivities. The call
      {[init s tol ism ps ~fs:fs s0 s0']} has as arguments:
      - [s], a session created with {!Ida.init},
      - [tol], the tolerances desired,
      - [sm], the solution method,
      - [~sens_params], the parameter information (see {!sens_params}),
      - [~fs], the sensitivity function,
      - [s0], initial values of the sensitivities for each parameter, and,
      - [s0'], initial values of the sensitivity derivatives for each
               parameter.

      If [~fs] is not given, an internal difference quotient routine
      is used.  In that case, or if the internal difference quotient
      routine will be specified in a subsequent call to
      {!Sensitivity.Quadrature.init}, then [sens_params] must be
      given with [pvals] set to non-[None].

      @idas <node6#ss:sensi_init> IDASensInit
      @idas <node6#sss:idafwdtolerances> IDASetSensParams
      @idas <node6#sss:idafwdtolerances> IDASensSStolerances
      @idas <node6#sss:idafwdtolerances> IDASensSVtolerances
      @idas <node6#sss:idafwdtolerances> IDASensEEtolerances
      @idas <node6> IDASetNonlinearSolverSensSim
      @idas <node6> IDASetNonlinearSolverSensStg *)
  val init : ('d, 'k) Ida.session
             -> ('d, 'k) tolerance
             -> sens_method
             -> ?sens_nlsolver:
                 (('d, 'k) Sundials_NonlinearSolver.Senswrapper.t, 'k,
                  (('d, 'k) Ida.session) Sundials_NonlinearSolver.integrator)
                                Sundials_NonlinearSolver.nonlinear_solver
             -> ?sens_params:sens_params
             -> ?fs:'d sensresfn
             -> ('d, 'k) Nvector.t array
             -> ('d, 'k) Nvector.t array
             -> unit

  (** Reinitializes the forward sensitivity computation.

      @idas <node6#ss:sensi_init> IDASensReInit *)
  val reinit : ('d, 'k) Ida.session
               -> sens_method
               -> ?sens_nlsolver:
                   (('d, 'k) Sundials_NonlinearSolver.Senswrapper.t, 'k,
                    (('d, 'k) Ida.session) Sundials_NonlinearSolver.integrator)
                                  Sundials_NonlinearSolver.nonlinear_solver
               -> ('d, 'k) Nvector.t array
               -> ('d, 'k) Nvector.t array
               -> unit

  (** Deactivates forward sensitivity calculations without deallocating
      memory. Sensitivities can be reactivated with {!reinit}.

      @idas <node6#ss:sensi_init> IDASensToggleOff *)
  val toggle_off : ('d, 'k) Ida.session -> unit

  (** Support for quadrature sensitivity equations.

      Adds a vector {% $s_\mathit{Q}$%} of {% $N_\mathit{Q}$%}
      quadrature sensitivities
      {% $\frac{\mathrm{d} s_\mathit{Q}}{\mathrm{d}t}
              = f_{\mathit{QS}}(t, y, \dot{y}, s, \dot{s}, \dot{y}_Q, p)$%}.
      This mechanism allows, in particular, the calculation of the
      sensitivities of the ‘pure’ {!Idas.Quadrature} variables, $y_Q$.

      @idas <node3#SECTION00354000000000000000> Quadratures depending on forward sensitivities
      @idas <node6#SECTION00640000000000000000> Integration of quadrature equations depending on forward sensitivities *)
  module Quadrature : sig (* {{{ *)
    (** {2:init Initialization} *)

    (** Arguments to {!quadsensrhsfn}.  *)
    type 'd quadsensrhsfn_args =
      {
        t : float;
        (** The value of the independent variable. *)

        y : 'd;
        (** The vector of dependent-variable values $y(t)$. *)

        y' : 'd;
        (** The vector of dependent-variable derivatives
            {% $\dot{y}(t)$ %}. *)

        s : 'd array;
        (** The array of sensitivity vectors. *)

        s' : 'd array;
        (** The array of sensitivity derivative vectors. *)

        q : 'd;
        (** The current value of the quadrature right-hand side $q$. *)

        tmp : 'd Ida.triple;
        (** Temporary storage vectors. *)
      }

    (** Functions defining quadrature sensitivities.
        They are passed the arguments:
        - [args], the current values of state, sensitivity, and quadrature
                  variables, and,
        - [sq'], an array of vectors for storing the computed values of
          {% $\dot{s}_\mathit{Q} =
              f_\mathit{QS}(t, y, \dot{y}, s, \dot{s}, \dot{y}_Q)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure}
        exception indicates a recoverable error. Any other exception is
        treated as an unrecoverable error.

        {warning The vectors in the function arguments should not
                 be accessed after the function returns.}

       @idas <node6#ss:user_fct_quad_sens> IDAQuadSensRhsFn *)
    type 'd quadsensrhsfn =
      'd quadsensrhsfn_args
      -> 'd array
      -> unit

    (** Activate the integration of quadrature sensitivities.  The
        arguments are:
        - [~fqs], a function that computes the right-hand sides
          of the quadrature sensitivities, and
        - [q0], an array of vectors specifying initial values for the
          quadrature sensitivities.

        If [~fqs] is not given, an internal difference quotient
        routine is used.  In that case, {!Sensitivity.init} must
        have been invoked with a [sens_params] whose [pvals] is
        non-[None].

        @idas <node6#ss:quad_sens_init> IDAQuadSensInit
        @raise QuadNotInitialized {!Quadrature.init} has not been called. *)
    val init : ('d, 'k) Ida.session -> ?fqs:'d quadsensrhsfn
                  -> ('d, 'k) Nvector.t array -> unit

    (** Reinitializes the quadrature sensitivity integration.

        @idas <node6#ss:quad_sens_init> IDAQuadSensReInit *)
    val reinit : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t array -> unit

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

        @idas <node6#ss:quad_sens_optional_input> IDASetQuadSensErrCon
        @idas <node6#ss:quad_sens_optional_input> IDAQuadSensSStolerances
        @idas <node6#ss:quad_sens_optional_input> IDAQuadSensSVtolerances
        @idas <node6#ss:quad_sens_optional_input> IDAQuadSensEEtolerances *)
    val set_tolerances : ('d, 'k) Ida.session -> ('d, 'k) tolerance -> unit

    (** {2:quadout Output functions} *)

    (** Returns the quadrature sensitivity solutions and time reached
        after a successful solver step. The given vectors are filled with
        values calculated during either {!Ida.solve_normal} or
        {!Ida.solve_one_step} and the value of the independent variable
        is returned.

        @idas <node6#ss:quad_sens_get> IDAGetQuadSens *)
    val get : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t array -> float

    (** Returns a single quadrature sensitivity vector after a successful
        solver step. Like {!get}, but the argument [i] specifies a specific
        vector to return.

        @idas <node6#ss:quad_sens_get> IDAGetQuadSens1
        @raise BadIS The index is not in the allowed range. *)
    val get1 : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t -> int -> float

    (** Returns the interpolated solution or derivatives of the quadrature
        sensitivity solution.

        [get_dky s dksq t k] computes the [k]th derivative at time [t],
        i.e.,
        {% $\frac{d^\mathtt{k}s_\mathit{Q}(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
        and stores it in [dksq]. The arguments must satisfy {% $t_n - h_u
        \leq \mathtt{t} \leq t_n$%}—where $t_n$ denotes
        {!Ida.get_current_time} and $h_u$ denotes
        {!Ida.get_last_step},—and
        {% $0 \leq \mathtt{k} \leq q_u$%}—where $q_u$ denotes
        {!Ida.get_last_order}.

        @idas <node6#ss:quad_sens_get> IDAGetQuadSensDky
        @raise BadIS The index is not in the allowed range.
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t array
                    -> float -> int -> unit

    (** Returns the interpolated solution or derivatives of a single
        quadrature sensitivity solution vector.
        [get_dky s dksq t k i] is like
        {!get_dky} but restricted to the [i]th sensitivity solution vector.

        @idas <node6#ss:quad_sens_get> IDAGetQuadSensDky1
        @raise BadK [k] is not in the range 0, 1, ..., [qlast].
        @raise BadT [t] is not in the allowed range. *)
    val get_dky1 : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t
                     -> float -> int -> int -> unit

    (** {2:get Querying the solver (optional output functions)} *)

    (** Returns the number of calls to the quadrature right-hand side
        function.

        @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensNumRhsEvals *)
    val get_num_rhs_evals : ('d, 'k) Ida.session -> int

    (** Returns the number of local error test failures due to quadrature
        variables.

        @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensNumErrTestFails *)
    val get_num_err_test_fails : ('d, 'k) Ida.session -> int

    (** Returns the quadrature error weights at the current time.

        @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensErrWeights *)
    val get_err_weights :
      ('d, 'k) Ida.session -> ('d, 'k) Nvector.t array -> unit

    (** Returns quadrature-related statistics. These are the
        number of calls to the quadrature function ([nfqevals]) and the
        number of error test failures due to quadrature variables
        ([nqetfails]).

      @idas <node6#ss:quad_sens_optional_output> IDAGetQuadSensStats
      @return ([nfqevals], [nqetfails]) *)
    val get_stats : ('d, 'k) Ida.session -> int * int

    (** {2:exceptions Exceptions} *)

    (** Quadrature integration was not initialized.

        @idas <node6#ss:quad_sens_get> IDA_NO_QUADSENS *)
    exception QuadSensNotInitialized

    (** The sensitivity quadrature function failed in an unrecoverable
        manner.

        @idas <node5#SECTION00642000000000000000> IDA_QSRHS_FAIL *)
    exception QuadSensRhsFuncFailure

    (** The sensitivity quadrature function failed at the first call.

        @idas <node5#SECTION00642000000000000000> IDA_FIRST_QSRHS_ERR *)
    exception FirstQuadSensRhsFuncFailure

    (** Convergence test failures occurred too many times due to repeated
        recoverable errors in the quadrature function. Also raised if the
        sensitivity quadrature function had repeated recoverable errors
        during the estimation of an initial step size (if quadrature
        variables are included in error tests).

        @idas <node6#SECTION00642000000000000000> IDA_REP_QSRHS_ERR *)
    exception RepeatedQuadSensRhsFuncFailure
  end (* }}} *)

  (** {3:calcic Initial Condition Calculation} *)

  (** Identical to {!Ida.calc_ic_ya_yd'}, but with the possibility of
      filling [s] and [s'] with the corrected sensitivity and sensitivity
      derivative values.

      @ida <node5#ss:idacalcic> IDACalcIC
      @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
      @idas <node6#sss:sens_optout_iccalc> IDAGetSensConsistentIC *)
  val calc_ic_ya_yd' :
    ('d, 'k) Ida.session
    -> ?y :('d, 'k) Nvector.t
    -> ?y':('d, 'k) Nvector.t
    -> ?s :('d, 'k) Nvector.t array
    -> ?s':('d, 'k) Nvector.t array
    -> ?varid:('d, 'k) Nvector.t
    -> float
    -> unit

  (** Identical to {!Ida.calc_ic_y}, but with the possibility of
      filling [s] with the corrected sensitivity values.

      @ida <node5#ss:idacalcic> IDACalcIC
      @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
      @idas <node6#sss:sens_optout_iccalc> IDAGetSensConsistentIC *)
  val calc_ic_y :
    ('d, 'k) Ida.session ->
    ?y:('d, 'k) Nvector.t -> ?s:('d, 'k) Nvector.t array -> float -> unit

  (** {2:sensout Output functions} *)

  (** Returns the sensitivity solution vectors after a successful solver
      step. The given vectors are filled with values calculated during
      either {!Ida.solve_normal} or {!Ida.solve_one_step} and the
      value of the independent variable is returned.

      @idas <node6#ss:sensi_get> IDAGetSens *)
  val get : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t array -> float

  (** Returns the interpolated solution or derivatives of the sensitivity
      solution vectors.

      [get_dky s dks t k] computes the [k]th derivative at time [t], i.e.,
      {% $\frac{d^\mathtt{k}s(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
      and stores it in [dks]. The arguments must satisfy {% $t_n - h_u \leq
      \mathtt{t} \leq t_n$%}—where $t_n$ denotes {!Ida.get_current_time}
      and $h_u$ denotes {!Ida.get_last_step},— and
      {% $0 \leq \mathtt{k} \leq q_u$%}—where
      $q_u$ denotes {!Ida.get_last_order}.

      This function may only be called after a successful return from either
      {!Ida.solve_normal} or {!Ida.solve_one_step}.

      @idas <node6#ss:sensi_get> IDAGetSensDky
      @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
      @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
  val get_dky :
    ('d, 'k) Ida.session -> ('d, 'k) Nvector.t array -> float -> int -> unit

  (** Returns a single sensitivity solution vector after a successful solver
      step. The given vector is filled with values calculated for the [i]th
      sensitivity vector during either {!Ida.solve_normal} or
      {!Ida.solve_one_step} and the value of the independent variable is
      returned.

      @idas <node6#ss:sensi_get> IDAGetSens1
      @raise BadIS The index [i] is not in the allowed range. *)
  val get1 : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t -> int -> float

  (** Returns the interpolated solution or derivatives of a single
      sensitivity solution vector. [get_dky s dks t k i] is like {!get_dky}
      but restricted to the [i]th sensitivity solution vector.

      @idas <node6#ss:sensi_get> IDAGetSensDky1
      @raise BadIS The index [i] is not in the allowed range.
      @raise BadK [k] is not in the range 0, 1, ..., [qlast].
      @raise BadT [t] is not in the allowed range. *)
  val get_dky1 :
    ('d, 'k) Ida.session -> ('d, 'k) Nvector.t -> float -> int -> int -> unit

  (** {2:set Modifying the solver (optional input functions)} *)

  (** Specify the integration tolerances for sensitivities.

      {b NB}: Unlike the other [set_tolerances] functions in [Idas], this
      one does {b not} call {!set_err_con} (which defaults to [false]).

      @idas <node6#ss:quad_sens_optional_input> IDASensSStolerances
      @idas <node6#ss:quad_sens_optional_input> IDASensSVtolerances
      @idas <node6#ss:quad_sens_optional_input> IDASensEEtolerances *)
  val set_tolerances : ('d, 'k) Ida.session -> ('d, 'k) tolerance -> unit

  (** Sets whether sensitivity variables are used in the error control
      mechanism. The default is [false].

      @idas <node6#ss:sens_optional_input> IDASetSensErrCon *)
  val set_err_con : ('d, 'k) Ida.session -> bool -> unit

  (** A difference quotient strategy. See {!set_dq_method}. *)
  type dq_method = DQCentered (** IDA_CENTERED *)
                 | DQForward  (** IDA_FORWARD *)

  (** Sets the difference quotient strategy when sensitivity equations
      are computed internally by the solver rather than via callback. A
      method and [dqrhomax] value must be given. The latter determines when
      to switch between simultaneous or separate approximations of the two
      terms in the sensitivity residual.

      @idas <node6#ss:sens_optional_input> IDASetSensDQMethod
      @idas <node3#ss:fwd_sensi> Forward Sensitivity Analysis *)
  val set_dq_method : ('d, 'k) Ida.session -> dq_method -> float -> unit

  (** Specifies the maximum number of nonlinear solver iterations for
      sensitivity variables permitted per step.

      @idas <node6#ss:sens_optional_input> IDASetSensMaxNonlinIters *)
  val set_max_nonlin_iters : ('d, 'k) Ida.session -> int -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the number of calls to the sensitivity residual function.

      @idas <node6#ss:sens_optional_output> IDAGetSensNumResEvals *)
  val get_num_res_evals : ('d, 'k) Ida.session -> int

  (** Returns the number of calls to the residual function due
      to the internal finite difference approximation of the sensitivity
      residual.

      @idas <node6#ss:sens_optional_output> IDAGetNumResEvalsSens *)
  val get_num_res_evals_sens : ('d, 'k) Ida.session -> int

  (** Returns the number of local error test failures for the sensitivity
      variables that have occurred.

      @idas <node6#ss:sens_optional_output> IDAGetSensNumErrTestFails *)
  val get_num_err_test_fails : ('d, 'k) Ida.session -> int

  (** Returns the number of calls made to the linear solver's setup function
      due to forward sensitivity calculations.

      @idas <node6#ss:sens_optional_output> IDAGetSensNumLinSolvSetups *)
  val get_num_lin_solv_setups : ('d, 'k) Ida.session -> int

  (** Summaries of sensitivity stats. *)
  type sensitivity_stats = {
    num_sens_evals : int;
        (** Number of calls to the sensitivity function. *)
    num_res_evals : int;
        (** Number of calls to the residual function to
            calculate sensitivities. *)
    num_err_test_fails : int;
        (** Number of local error test failures for sensitivity variables. *)
    num_lin_solv_setups : int;
        (** Number of setups calls to the linear solver for sensitivity
            calculations. *)
  }

  (** Returns the sensitivity-related statistics as a group.

      @idas <node6#ss:sens_optional_output> IDAGetSensStats *)
  val get_stats : ('d, 'k) Ida.session -> sensitivity_stats

  (** Returns the sensitivity error weights at the current time.

      @idas <node6#ss:sens_optional_output> IDAGetSensErrWeights
      @idas <node3#e:errwt> Eq. (2.7) IVP solution (W_i) *)
  val get_err_weights : ('d, 'k) Ida.session
                          -> ('d, 'k) Nvector.t array -> unit

  (** Returns the number of nonlinear iterations performed for sensitivity
      calculations.

      @idas <node6#ss:sens_optional_output> IDAGetSensNumNonlinSolvIters *)
  val get_num_nonlin_solv_iters : ('d, 'k) Ida.session -> int


  (** Returns the number of nonlinear convergence failures that have occurred
      during sensitivity calculations.

      @idas <node6#ss:sens_optional_output> IDAGetSensNumNonlinSolvConvFails *)
  val get_num_nonlin_solv_conv_fails : ('d, 'k) Ida.session -> int

  (** Returns both the numbers of nonlinear iterations performed [nniters]
      and nonlinear convergence failures [nncfails] during sensitivity
      calculations.

      @idas <node6#ss:sens_optional_output> IDAGetSensNonlinSolvStats
      @return ([nniters], [nncfails]) *)
  val get_nonlin_solv_stats : ('d, 'k) Ida.session -> int * int

  (** {2:exceptions Exceptions} *)

  (** Sensitivity analysis was not initialized.

      @idas <node6#ss:sensi_get> IDA_NO_SENS *)
  exception SensNotInitialized

  (** The sensitivity residual function failed in an unrecoverable manner.

      @idas <node6#SECTION00624000000000000000> IDA_SRES_FAIL *)
  exception SensResFuncFailure

  (** Too many convergence test failures, or unable to estimate the initial
      step size, due to repeated recoverable errors in the sensitivity
      residual function.

      @idas <node6#SECTION00624000000000000000> IDA_REP_SRES_ERR *)
  exception RepeatedSensResFuncFailure

  (** The index passed to identify a particular sensitivity is invalid.

      @idas <node6#SECTION00625000000000000000> IDA_BAD_IS *)
  exception BadSensIdentifier
end (* }}} *)

(** (Adjoint) Sensitivity analysis of DAEs with respect to their parameters.

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

    @idas <node7#ss:skeleton_adj> Enhanced Skeleton for Adjoint Sensitivity Analysis
    @idas <node7#s:adjoint> Using IDAS for Adjoint Sensitivity Analysis
    @idas <node3#ss:adj_sensi> Adjoint sensitivity analysis *)
module Adjoint : sig (* {{{ *)
  (** A backward session with the IDAS solver. Multiple backward sessions
      may be associated with a single parent session.

      @idas <node7#sss:idainitb> Backward problem initialization functions *)
  type ('d, 'k) bsession = ('d, 'k) Ida_impl.AdjointTypes.bsession

  (** Alias for backward sessions based on serial nvectors. *)
  type 'k serial_bsession = (RealArray.t, 'k) bsession
                            constraint 'k = [>Nvector_serial.kind]

  (** {2:fwd Forward solution} *)

  (** Specifies the type of interpolation to use between checkpoints.

      @idas <node3#ss:checkpointing> Checkpointing scheme *)
  type interpolation = IPolynomial    (** {cconst IDA_POLYNOMIAL} *)
                     | IHermite       (** {cconst IDA_HERMITE} *)

  (** Activates the forward-backward problem. The arguments specify the number
      of integration steps between consecutive checkpoints, and the type of
      variable-degree interpolation.

      @idas <node7#sss:idaadjinit> IDAAdjInit *)
  val init : ('d, 'k) Ida.session -> int -> interpolation -> unit

  (** Integrates the forward problem over an interval and saves
      checkpointing data. The arguments are the next time at which a solution
      is desired ([tout]) and two vectors to receive the computed results
      ([y] and [y']). The function returns a triple [tret, ncheck, sr]: the time
      reached by the solver, the cumulative number of checkpoints stored, and
      whether [tout] was reached. The solver takes internal steps until it
      has reached or just passed the [tout] parameter {cconst IDA_NORMAL},
      it then interpolates to approximate [y(tout)].

      @idas <node7#sss:idasolvef> IDASolveF
      @raise Ida.IllInput           One of the inputs is invalid.
      @raise Ida.TooMuchWork        Could not reach [tout] in [mxstep] steps
      @raise Ida.TooMuchAccuracy    Could not satisfy the demanded accuracy
      @raise Ida.ErrFailure         Too many error test failures.
      @raise Ida.ConvergenceFailure Too many convergence test failures.
      @raise Ida.LinearSetupFailure Unrecoverable failure in linear solver setup function.
      @raise Ida.LinearSolveFailure Unrecoverable failure in linear solver solve function.
      @raise AdjointNotInitialized  The {!init} function has not been called. *)
  val forward_normal :
    ('d, 'k) Ida.session
    -> float
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) Nvector.t
    -> float * int * Ida.solver_result

  (** Integrates the forward problem over an interval and saves
      checkpointing data. The arguments are the next time at which a solution
      is desired ([tout]) and two vectors to receive the computed results
      ([y] and [y']). The function returns a triple [tret, ncheck, sr]: the
      time reached by the solver, the cumulative number of checkpoints
      stored, and whether [tout] was reached. The solver takes one step
      {cconst IDA_ONE_STEP} and returns the solution reached.

      @idas <node7#sss:idasolvef> IDASolveF
      @raise Ida.IllInput           One of the inputs is invalid.
      @raise Ida.TooMuchWork        Could not reach [tout] in [mxstep] steps
      @raise Ida.TooMuchAccuracy    Could not satisfy the demanded accuracy
      @raise Ida.ErrFailure         Too many error test failures.
      @raise Ida.ConvergenceFailure Too many convergence test failures.
      @raise Ida.LinearSetupFailure Unrecoverable failure in linear solver setup function.
      @raise Ida.LinearSolveFailure Unrecoverable failure in linear solver solve function.
      @raise AdjointNotInitialized  The {!init} function has not been called. *)
  val forward_one_step :
    ('d, 'k) Ida.session
    -> float
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) Nvector.t
    -> float * int * Ida.solver_result

  (** {2:linear Linear solvers} *)

  (** Linear solvers used in backward problems.

      @idas <node7#sss:lin_solv_b> Linear Solver Initialization Functions *)
  type ('data, 'kind) session_linear_solver =
    ('data, 'kind) Ida_impl.AdjointTypes.session_linear_solver

  (** Alias for linear solvers that are restricted to serial nvectors. *)
  type 'kind serial_session_linear_solver =
    (RealArray.t, 'kind) session_linear_solver
    constraint 'kind = [>Nvector_serial.kind]

  (** Workspaces with three temporary vectors. *)
  type 'd triple = 'd * 'd * 'd

  (** Arguments common to Jacobian callback functions.

      @idas <node7#ss:densejac_b> IDALsJacFnB
      @idas <node7#ss:jactimesvec_b> IDAJacTimesVecFnB
      @idas <node7#ss:psolve_b> IDALsPrecSolveFnB
      @idas <node7#ss:psetup_b> IDALsPrecSetupFnB *)
  type ('t, 'd) jacobian_arg = ('t, 'd) Ida_impl.AdjointTypes.jacobian_arg =
    {
      jac_t : float;        (** The independent variable. *)
      jac_y : 'd;           (** The forward solution vector. *)
      jac_y' : 'd;          (** The forward derivatives vector. *)
      jac_yb : 'd;          (** The backward solution vector. *)
      jac_yb' : 'd;         (** The forward derivatives vector. *)
      jac_resb : 'd;        (** The current residual for the backward problem. *)
      jac_coef : float;     (** The scalar {% $c_\mathit{jB}$%} in the
                                system Jacobian, proportional to the inverse
                                of the step size. *)
      jac_tmp : 't;         (** Workspace data. *)
    }

  (** Direct Linear Solvers operating on dense, banded, and sparse matrices. *)
  module Dls : sig (* {{{ *)
    include module type of Sundials_LinearSolver.Direct

    (** Callback functions that compute dense approximations to a Jacobian
        matrix without forward sensitivities. In the call
        [jac arg jm], [arg] is a {!jacobian_arg} with three work
        vectors and the computed Jacobian must be stored in [jm].

        The callback should load the [(i,j)]th entry of [jm] with
        {% $\frac{\partial F_i}{\partial y_j} + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%},
        i.e., the partial derivative of the [i]th equation with respect to
        the [j]th variable, evaluated at the values of [t], [y], and [y']
        obtained from [arg]. Only nonzero elements need be loaded into [jm].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jm] should
                 be accessed after the function has returned.}

        @idas <node7#ss:densejac_b> IDALsJacFnB *)
    type 'm jac_fn_no_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg -> 'm -> unit

    (** Callback functions that compute dense approximations to a Jacobian
        matrix with forward sensitivities. In the call
        [jac arg ys yps jm], [arg] is a {!jacobian_arg} with
        three work vectors, [ys] contains the sensitivities of the forward
        solution, [yps] contains the derivatives of the forward solution
        sensitivities, and the computed Jacobian must be stored in [jm].

        The callback should load the [(i,j)]th entry of [jm] with
        {% $\frac{\partial F_i}{\partial y_j} + c_j\frac{\partial F_i}{\partial\dot{y}_j}$%},
        i.e., the partial derivative of the [i]th equation with respect to
        the [j]th variable, evaluated at the values of [t], [y], and [y']
        obtained from [arg]. Only nonzero elements need be loaded into [jm].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jm] should
                 be accessed after the function has returned.}

        @noidas <node7#ss:densejac_bs> IDALsJacFnBS *)
    type 'm jac_fn_with_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg
      -> RealArray.t array
      -> RealArray.t array
      -> 'm
      -> unit

    (** Callback functions that compute dense approximations to a Jacobian
        matrix.

        @noidas <node7#ss:densejac_b> IDALsJacFnB
        @noidas <node7#ss:densejac_bs> IDALsJacFnBS *)
    type 'm jac_fn =
        NoSens of 'm jac_fn_no_sens
        (** Does not depend on forward sensitivities. *)
      | WithSens of 'm jac_fn_with_sens
        (** Depends on forward sensitivities. *)

    (** Create an Idas-specific linear solver from a Jacobian approximation
        function and a generic direct linear solver.
        The Jacobian approximation function is optional for dense and banded
        solvers (if not given an internal difference quotient approximation is
        used), but must be provided for other solvers (or [Invalid_argument]
        is raised).

        @noidas <node> IDASetLinearSolverB
        @noidas <node> IDASetJacFnB
        @noidas <node> IDASetJacFnBS *)
    val solver :
      ?jac:'m jac_fn ->
      ('m, 'kind, 't) LinearSolver.Direct.serial_linear_solver ->
      'kind serial_session_linear_solver

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by a direct
        linear solver.

        @ida <node5#sss:optout_dls> IDAGetWorkSpace
        @idas <node7#SECTION007210100000000000000> IDAGetAdjIDABmem
        @return ([real_size], [integer_size]) *)
    val get_work_space : 'k serial_bsession -> int * int

    (** Returns the number of calls made by a direct linear solver to the
        Jacobian approximation function.

        @ida <node5#sss:optout_dls> IDAGetNumJacEvals
        @idas <node7#SECTION007210100000000000000> IDAGetAdjIDABmem *)
    val get_num_jac_evals : 'k serial_bsession -> int

    (** Returns the number of calls to the residual callback due to
        the finite difference Jacobian approximation.

        @ida <node5#sss:optout_dls> IDAGetNumResEvals
        @idas <node7#SECTION007210100000000000000> IDAGetAdjIDABmem *)
    val get_num_lin_res_evals : 'k serial_bsession -> int

  end (* }}} *)

  (** Scaled Preconditioned Iterative Linear Solvers

      @idas <node7#ss:optional_output_b> Optional output functions for the backward problem.
      @idas <node7#ss:psolve_b> IDALsPrecSolveFnB
      @idas <node7#ss:psetup_b> IDALsPrecSetupFnB *)
  module Spils : sig (* {{{ *)
    include module type of Sundials_LinearSolver.Iterative

    (** {3:precond Preconditioners} *)

    (** Callback functions that solve a linear system involving a
        preconditioner matrix without forward sensitivities.
        In the call [prec_solve_fn jac r z delta],
        - [jac] is a {!jacobian_arg} with no work vectors,
        - [r] is the right-hand side vector,
        - [z] is computed to solve {% $Pz = r$%}, and
        - [delta] is the input tolerance.
        $P$ is a preconditioner matrix, which approximates, however crudely,
        the Jacobian matrix
        {% $\frac{\partial F}{\partial y} + \mathtt{arg.jac\_coef}\frac{\partial F}{\partial\dot{y}}$%}.
        If the solution is found via an iterative method, it must satisfy
        {% $\sqrt{\sum_i (\mathit{Res}_i \cdot \mathit{ewt}_i)^2}
              < \mathtt{delta}$%},
        where {% $\mathit{Res} = r - Pz$%} and {% $\mathit{ewt}$%} comes from
        {!get_err_weights}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac], [r], and [z] should not
                 be accessed after the function has returned.}

        @idas <node7#ss:psolve_b> IDALsPrecSolveFnB *)
    type 'd prec_solve_fn =
      (unit, 'd) jacobian_arg
      -> 'd
      -> 'd
      -> float
      -> unit

    (** Callback functions that solve a linear system involving a
        preconditioner matrix with forward sensitivities.
        In the call [prec_solve_fn jac ys yps r z delta],
        - [jac] is a {!jacobian_arg} with no work vectors,
        - [ys] contains the sensitivities of the forward solution,
        - [yps] contains the derivatives of the forward solution
                sensitivities,
        - [r] is the right-hand side vector,
        - [z] is computed to solve {% $Pz = r$%}, and
        - [delta] is the input tolerance.
        $P$ is a preconditioner matrix, which approximates, however crudely,
        the Jacobian matrix
        {% $\frac{\partial F}{\partial y} + \mathtt{arg.jac\_coef}\frac{\partial F}{\partial\dot{y}}$%}.
        If the solution is found via an iterative method, it must satisfy
        {% $\sqrt{\sum_i (\mathit{Res}_i \cdot \mathit{ewt}_i)^2}
              < \mathtt{delta}$%},
        where {% $\mathit{Res} = r - Pz$%} and {% $\mathit{ewt}$%} comes from
        {!get_err_weights}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac], [ys], [yps], [r], and [z] should not
                 be accessed after the function has returned.}

        @noidas <node7#ss:psolve_bs> IDALsPrecSolveFnBS *)
    type 'd prec_solve_fn_with_sens =
      (unit, 'd) jacobian_arg
      -> 'd array
      -> 'd array
      -> 'd
      -> 'd
      -> float
      -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related data
        need by {!prec_solve_fn} without forward sensitivities.
        The only argument is a {!jacobian_arg} with no work vectors.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of the argument should not be accessed after
                 the function has returned.}

        @idas <node7#ss:psetup_b> IDALsPrecSetupFnB *)
    type 'd prec_setup_fn =
      (unit, 'd) jacobian_arg
      -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related data
        need by {!prec_solve_fn} with forward sensitivities.
        In the call [prec_setup_fn jac ys yps],
        - [jac] is a {!jacobian_arg} with no work vectors,
        - [ys] contains the sensitivities of the forward solution, and
        - [yps] contains the derivatives of the forward solution
                sensitivities.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of the arguments should not be accessed after
                 the function has returned.}

        @noidas <node7#ss:psetup_bs> IDALsPrecSetupFnBS *)
    type 'd prec_setup_fn_with_sens =
      (unit, 'd) jacobian_arg
      -> 'd array
      -> 'd array
      -> unit

    (** Specifies a preconditioner and its callback functions.
        The following functions and those in {!Idas_bbd} construct
        preconditioners.

        The {!prec_solve_fn} is mandatory. The {!prec_setup_fn} can be
        omitted if not needed.

        @idas <node7#SECTION00729400000000000000> IDALsSetPreconditionerB
        @noidas <node7> IDALsSetPreconditionerBS
        @idas <node7#ss:psolve_b> IDALsPrecSolveFnB
        @idas <node7#ss:psetup_b> IDALsPrecSetupFnB
        @noidas <node7#ss:psolve_bs> IDALsPrecSolveFnBS
        @noidas <node7#ss:psetup_bs> IDAlsPrecSetupFnBS *)
    type ('d, 'k) preconditioner =
      ('d, 'k) Ida_impl.AdjointTypes.SpilsTypes.preconditioner

    (** No preconditioning.  *)
    val prec_none : ('d, 'k) preconditioner

    (** Left preconditioning without forward sensitivities.
        {% $Pz = r$%}, where $P$ approximates, perhaps crudely,
        {% $J = \frac{\partial F}{\partial y} + c_j\frac{\partial F}{\partial\dot{y}}$%}. *)
    val prec_left :
      ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

    (** Left preconditioning with forward sensitivities.
        {% $Pz = r$%}, where $P$ approximates, perhaps crudely,
        {% $J = \frac{\partial F}{\partial y} + c_j\frac{\partial F}{\partial\dot{y}}$%}. *)
    val prec_left_with_sens :
      ?setup:'d prec_setup_fn_with_sens
      -> 'd prec_solve_fn_with_sens
      -> ('d, 'k) preconditioner

    (** {3:lsolvers Solvers} *)

    (** Callback functions that preprocess or evaluate Jacobian-related
        data needed by the jac_times_vec_fn. In the call
        [jac_times_setup_fn arg], [arg] is a {!jacobian_arg} with no
        work vectors.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [arg] should not be accessed after the
                 function has returned.}

        @noidas <node> IDASetJacTimesVecFnB
        @noidas <node> IDALsJacTimesSetupFnB *)
    type 'd jac_times_setup_fn_no_sens = (unit, 'd) jacobian_arg -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related
        data needed by the jac_times_vec_fn. In the call
        [jac_times_setup_fn arg s], [arg] is a {!jacobian_arg} with no
        work vectors and [s] is an array of forward sensitivity vectors.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The elements of [arg] should not be accessed after the
                 function has returned.}

        @noidas <node> IDASetJacTimesVecFnBS
        @noidas <node> IDALsJacTimesSetupFnBS *)
    type 'd jac_times_setup_fn_with_sens =
      (unit, 'd) jacobian_arg -> 'd array -> unit

    (** Callback functions that compute the Jacobian times a vector without
        forward sensitivities.
        In the call [jac_times_vec_fn arg v jv],
        - [arg] is a {!jacobian_arg} with two work vectors,
        - [v] is the vector multiplying the Jacobian, and
        - [jv] is the vector in which to store the
               result—{% $\mathtt{jv} = J\mathtt{v}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor [v] or [jv] should be
                 accessed after the function has returned.}

        @idas <node7#ss:jactimesvec_b> IDALsJacTimesVecFnB *)
    type 'd jac_times_vec_fn_no_sens =
      ('d, 'd) jacobian_arg
      -> 'd
      -> 'd
      -> unit

    (** Callback functions that compute the Jacobian times a vector with
        forward sensitivities.
        In the call [jac_times_vec_fn arg ys yps v jv],
        - [arg] is a {!jacobian_arg} with two work vectors,
        - [ys] contains the sensitivities of the forward solution,
        - [yps] contains the derivatives of the forward solution
                sensitivities.
        - [v] is the vector multiplying the Jacobian, and
        - [jv] is the vector in which to store the
               result—{% $\mathtt{jv} = J\mathtt{v}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg], [ys], [yps], [v] nor [jv]
                 should be accessed after the function has returned.}

        @noidas <node7#ss:jactimesvec_bs> IDALsJacTimesVecFnBS *)
    type 'd jac_times_vec_fn_with_sens =
      ('d, 'd) jacobian_arg
      -> 'd array
      -> 'd array
      -> 'd
      -> 'd
      -> unit

    (** Callback functions that compute the Jacobian times a vector.

        @idas <node7#ss:jactimesvec_b> IDALsJacTimesVecFnB
        @noidas <node7#ss:jactimesvec_bs> IDALsJacTimesVecFnBS *)
    type 'd jac_times_vec_fn =
      | NoSens of 'd jac_times_setup_fn_no_sens option
                  * 'd jac_times_vec_fn_no_sens
        (** Does not depend on forward sensitivities. *)
      | WithSens of 'd jac_times_setup_fn_with_sens option
                    * 'd jac_times_vec_fn_with_sens
        (** Depends on forward sensitivities. *)
        (** Depends on forward sensitivities. *)

    (** Create a Idas-specific linear solver from a generic iterative
        linear solver.

        NB: the [jac_times_setup] argument is not supported in
            {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

        @noidas <node> IDASetLinearSolverB
        @noidas <node> IDASetJacTimesVecFnB
        @noidas <node> IDASetJacTimesVecFnBS *)
    val solver :
      ('d, 'k, 'f) LinearSolver.Iterative.linear_solver
      -> ?jac_times_vec:'d jac_times_vec_fn
      -> ('d, 'k) preconditioner
      -> ('d, 'k) session_linear_solver

    (** {3:set Solver parameters} *)

    (** Sets the factor by which the Krylov linear solver's convergence test
        constant is reduced from the Newton iteration test constant.
        This factor must be >= 0; passing 0 specifies the default (0.05).

        @idas <node7#SECTION00729400000000000000> IDASetEpsLinB *)
    val set_eps_lin : ('d, 'k) bsession -> float -> unit

    (** Sets the increment factor ([dqincfac]) to use in the difference-quotient
        approximation for the backward problem.

        @ida <node5> IDASetIncrementFactorB *)
    val set_increment_factor : ('d, 'k) bsession -> float -> unit

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the
        linear solver.

        @idas <node5#sss:optout_spils> IDAGetLinWorkSpace
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
        @return ([real_size], [integer_size]) *)
    val get_work_space : ('d, 'k) bsession -> int * int

    (** Returns the cumulative number of linear iterations.

        @idas <node5#sss:optout_spils> IDAGetNumLinIters
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_lin_iters : ('d, 'k) bsession -> int

    (** Returns the cumulative number of linear convergence failures.

        @idas <node5#sss:optout_spils> IDAGetNumLinConvFails
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_lin_conv_fails : ('d, 'k) bsession -> int

    (** Returns the number of calls to the setup function.

        @idas <node5#sss:optout_spils> IDAGetNumPrecEvals
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_prec_evals : ('d, 'k) bsession -> int

    (** Returns the cumulative number of calls to the preconditioner solve
        function.

        @idas <node5#sss:optout_spils> IDAGetNumPrecSolves
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_prec_solves : ('d, 'k) bsession -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        setup function.

        @since 3.0.0
        @noidas <node> IDAGetNumJTSetupEvals
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_jtsetup_evals : ('d, 'k) bsession -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        function.

        @idas <node5#sss:optout_spils> IDAGetNumJtimesEvals
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_jtimes_evals : ('d, 'k) bsession -> int

    (** Returns the number of calls to the residual callback for
        finite difference Jacobian-vector product approximation. This counter is
        only updated if the default difference quotient function is used.

        @idas <node5#sss:optout_spils> IDAGetNumLinResEvals
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_lin_res_evals : ('d, 'k) bsession -> int

    (** {3:lowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change preconditioner functions. These low-level functions are provided
        for experts who want to avoid resetting internal counters and other
        associated side-effects. *)

    (** Change the preconditioner functions without using forward
        sensitivities.

        @idas <node7#SECTION00729400000000000000> IDASetPreconditionerB
        @idas <node7#ss:psolve_b> IDALsPrecSolveFnB
        @idas <node7#ss:psetup_b> IDALsPrecSetupFnB *)
    val set_preconditioner :
      ('d,'k) bsession
      -> ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> unit

    (** Change the preconditioner functions using forward sensitivities.

        @noidas <node7> IDASetPreconditionerBS
        @noidas <node7#ss:psolve_bs> IDALsPrecSolveFnBS
        @noidas <node7#ss:psetup_bs> IDALsPrecSetupFnBS *)
    val set_preconditioner_with_sens :
      ('d,'k) bsession
      -> ?setup:'d prec_setup_fn_with_sens
      -> 'd prec_solve_fn_with_sens
      -> unit

    (** Change the Jacobian-times-vector function.

        @idas <node7#SECTION00729400000000000000> IDASetJacTimesVecFnB
        @noidas <node7> IDASetJacTimesVecFnBS
        @idas <node7#ss:jactimesvec_b> IDALsJacTimesVecFnB
        @noidas <node7#ss:jactimesvec_bs> IDAJacTimesVecFnBS *)
    val set_jac_times :
      ('d,'k) bsession
      -> 'd jac_times_vec_fn
      -> unit

    (** Remove a Jacobian-times-vector function and use the default
        implementation.

        @idas <node7#SECTION00729400000000000000> IDASetJacTimesVecFnB *)
    val clear_jac_times : ('d, 'k) bsession -> unit
  end (* }}} *)

  (** Alternate Linear Solvers.

      @noidas <node8#s:new_linsolv> Providing Alternate Linear Solver Modules *)
  module Alternate : sig (* {{{ *)

    (** Functions that initialize linear solver data, like counters and
        statistics.

        Raising any exception in this function (including
        {!Sundials.RecoverableFailure}) is treated as an unrecoverable error.

        @idas <node8#SECTION00810000000000000000> linit *)
    type ('data, 'kind) linit = ('data, 'kind) bsession -> unit

  (** Functions that prepare the linear solver for subsequent calls to
      {{!callbacks}lsolve}. The call [lsetup s y y' res tmp] has as
      arguments

      - [s], the solver session,
      - [y],  the predicted $y$ vector for the current internal step,
      - [y'], the predicted {% $\dot{y}$%} vector for the current internal
              step,
      - [res], the value of the residual function at [y] and [y'], i.e.
               {% $F(t_n, y_{\text{pred}}, \dot{y}_{\text{pred}})$%}, and,
      - [tmp], temporary variables for use by the routine.

      This function may raise a {!Sundials.RecoverableFailure} exception to
      indicate that a recoverable error has occurred. Any other exception is
      treated as an unrecoverable error.

      {warning The vectors in {!Ida.Alternate.lsetup_args} should not be
               accessed after the function returns.}

      @idas <node8#SECTION00820000000000000000> lsetup *)
    type ('data, 'kind) lsetup =
      ('data, 'kind) bsession
      -> 'data Ida.Alternate.lsetup_args
      -> unit

    (** Functions that solve the linear equation $Mx = b$.
        $M$ is a preconditioning matrix chosen by the user, and $b$ is the
        right-hand side vector calculated within the function.
        $M$ should approximate {% $J = \frac{\partial F}{\partial y} + c_j\frac{\partial F}{\partial \dot{y}}$%},
        and $c_j$ is available through {!get_cj}.
        The call [lsolve s b weight ycur y'cur rescur] has as arguments:

        - [s], the solver session,
        - [args], the current approximation to the solution, and,
        - [b], for returning the calculated solution,

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        @idas <node8#SECTION00830000000000000000> lsolve
        @idas <node3#e:DAE_Jacobian> IVP solution (Eq. 2.5) *)
    type ('data, 'kind) lsolve =
      ('data, 'kind) bsession
      -> 'data Ida.Alternate.lsolve_args
      -> 'data
      -> unit

    (** The callbacks needed to implement an alternate linear solver. *)
    type ('data, 'kind) callbacks =
      {
        linit  : ('data, 'kind) linit option;
        lsetup : ('data, 'kind) lsetup option;
        lsolve : ('data, 'kind) lsolve;
      }

    (** Creates a linear solver from a function returning a set of
        callbacks. The creation function is passed a session and a vector.
        The latter indicates the problem size and can, for example, be
        cloned. *)
    val solver :
          (('data, 'kind) bsession
            -> ('data, 'kind) Nvector.t
            -> ('data, 'kind) callbacks)
          -> ('data, 'kind) session_linear_solver

    (** {3:internals Solver internals} *)

    (** Returns the current [cj] value. *)
    val get_cj : ('data, 'kind) bsession -> float

    (** Returns the current [cjratio] value. *)
    val get_cjratio : ('data, 'kind) bsession -> float
  end (* }}} *)

  (** {2:bsolve Backward solutions} *)

  type 'd bresfn_args = 'd Ida_impl.AdjointTypes.bresfn_args =
    {
      t : float;
      (** The value of the independent variable. *)

      y : 'd;
      (** The vector of dependent-variable values $y(t)$. *)

      y' : 'd;
      (** The vector of dependent-variable derivatives {% $\dot{y}(t)$ %}. *)

      yb : 'd;
      (** The vector of backward dependent-variable values $y_B(t)$. *)

      yb' : 'd;
      (** The vector of backward dependent-variable derivatives
          {% $\dot{y}_B(t)$ %}. *)
    }

  (** Backward functions without forward sensitivities. They are passed
      the arguments:
      - [args], the current values of forward and backward state variables,
                and,
      - [rb], a vector for storing the residual value
              {% $F_B(t, y, \dot{y}, y_B, \dot{y}_B)$%}.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning Vectors held in this function's arguments should not
               be accessed after the function returns.}

      @idas <node7#ss:ODErhs_b> IDAResFnB
      @idas <node3#e:adj_eqns> Eq 2.19, Adjoint sensitivity analysis *)
  type 'd bresfn_no_sens = 'd bresfn_args -> 'd -> unit

  (** Backward functions with forward sensitivities. They are passed the
      arguments:
      - [args], the current values of forward and backward state variables,
      - [s], the array of forward sensitivity vectors,
      - [s'], the array of forward sensitivity derivative vectors, and,
      - [resb], a vector for storing the residual value
              {% $F_B(t, y, \dot{y}, s, \dot{s}, y_B, \dot{y}_B)$%}.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      {warning Vectors held in this function's arguments should not
               be accessed after the function returns.}

      @idas <node7#ss:ODErhs_bs> IDAResFnBS
      @idas <node3#e:adj1_eqns> Eq 2.21, Adjoint sensitivity analysis *)
  type 'd bresfn_with_sens = 'd bresfn_args -> 'd array -> 'd array
                           -> 'd -> unit

  (** Functions that evaluate the right-hand side of a backward DAE system
      with or without forward sensitivities. *)
  type 'd bresfn =
    | NoSens of 'd bresfn_no_sens
        (** No dependency on forward sensitivities. *)
    | WithSens of 'd bresfn_with_sens
        (** Dependency on forward sensitivities. *)

  (** Tolerance specifications. *)
  type ('d, 'k) tolerance =
      SStolerances of float * float
      (** [(rel, abs)] : scalar relative and absolute tolerances. *)
    | SVtolerances of float * ('d, 'k) Nvector.t
      (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)

  (** Creates and initializes a backward session attached to an existing
      (forward) session. The call
      {[init_backward s linsolv tol fb tb0 yb0 yb0']} has as arguments:
      - [s], the parent (forward) session,
      - [linsolv], the linear solver to use,
      - [tol], the integration tolerances,
      - [fb], the backward residual function,
      - [varid], (optionally) classifies variables as algebraic or
                 differential,
      - [tb0], specifies the endpoint where final conditions are provided
               for the backward problem, which is normally the endpoint of
               forward integration, and,
      - [yb0], the final value of the backward variables.
      - [yb0'], the final value of the backward derivatives.

      This function does everything necessary to initialize a backward
      session, i.e., it makes the calls referenced below. The
      {!backward_normal} and {!backward_one_step} functions may be called
      directly.

      If an [nlsolver] is not specified, then the
      {{!Sundials_NonlinearSolver.Newton}Newton} module is used by default.
      The [nlsolver] must be of type
      {{!Sundials_NonlinearSolver.nonlinear_solver_Type}RootFind}, otherwise an
      {!Ida.IllInput} exception is raised.

      @idas <node7#sss:idainitb> IDACreateB
      @idas <node7#sss:idainitb> IDAInitB
      @idas <node7#sss:idainitb> IDAInitBS
      @idas <node7#sss:idatolerances_b> IDASStolerancesB
      @idas <node7#sss:idatolerances_b> IDASVtolerancesB
      @raise AdjointNotInitialized The {!init} function has not been called.
      @raise BadFinalTime      The final time is outside the interval over which the forward problem was solved. *)
  val init_backward :
       ('d, 'k) Ida.session
    -> ('d, 'k) tolerance
    -> ?nlsolver:('d, 'k,
                  (('d, 'k) session) Sundials_NonlinearSolver.integrator)
                  Sundials_NonlinearSolver.nonlinear_solver
    -> lsolver:('d, 'k) session_linear_solver
    -> 'd bresfn
    -> ?varid:('d, 'k) Nvector.t
    -> float
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) bsession

  (** Support for backward quadrature equations that may or may
      not depend on forward sensitivities.

      @idas <node7#SECTION007211000000000000000> Backward integration of quadrature equations *)
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

        y' : 'd;
        (** The vector of dependent-variable derivatives
            {% $\dot{y}(t)$ %}. *)

        yb : 'd;
        (** The vector of backward dependent-variable values $y_B(t)$. *)

        yb' : 'd;
        (** The vector of backward dependent-variable derivatives
            {% $\dot{y}_B(t)$ %}. *)
      }

    (** Functions defining backward quadrature variables without forward
        sensitivities.  They are passed the arguments:
        - [args], the current values of forward and backward state
                  variables, and,
        - [qb'], a vector for storing the computed value of
                 {% $\dot{y}_\mathit{BQ} =
                     f_\mathit{BQ}(t, y, \dot{y}, y_B, \dot{y}_B)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure}
        exception indicates a recoverable error. Any other exception is
        treated as an unrecoverable error.

        {warning Vectors held in this function's arguments should not
                 be accessed after the function returns.}

        @idas <node7#sss:rhs_quad_B> IDAQuadRhsFnB *)
    type 'd bquadrhsfn_no_sens = 'd bquadrhsfn_args -> 'd -> unit

    (** Functions defining backward quadrature variables that
        depend on forward sensitivities.  They are passed the arguments:

        - [args], current values of forward and backward state variables,
        - [s], the array of forward sensitivity vectors,
        - [s'], the array of forward sensitivity derivative vectors, and,
        - [qb'], a vector for storing the computed value of
             {% $\dot{y}_\mathit{BQ} =
               f_\mathit{BQ}(t, y, \dot{y}, s, \dot{s}, y_B, \dot{y}_B)$%}.

        Within the function, raising a {!Sundials.RecoverableFailure}
        exception indicates a recoverable error. Any other exception is
        treated as an unrecoverable error.

        {warning Vectors held in this function's arguments should not
                 be accessed after the function returns.}

        @idas <node7#sss:rhs_quad_sens_B> IDAQuadRhsFnBS *)
    type 'd bquadrhsfn_with_sens =
      'd bquadrhsfn_args -> 'd array -> 'd array -> 'd -> unit

    (** These functions compute the quadrature equation right-hand side for
        the backward problem. *)
    type 'd bquadrhsfn =
      | NoSens of 'd bquadrhsfn_no_sens
        (** Does not depend on forward sensitivities. *)
      | WithSens of 'd bquadrhsfn_with_sens
        (** Depends on forward sensitivities. *)

    (** This function activates the integration of quadrature equations.
        The arguments are the function that computes the right-hand side of
        the backward quadrature equations, and a vector giving the values
        of the quadrature variables at [tB0].

        @idas <node7#SECTION007211000000000000000> IDAQuadInitB
        @idas <node7#SECTION007211000000000000000> IDAQuadInitBS *)
    val init :
      ('d, 'k) bsession -> 'd bquadrhsfn -> ('d, 'k) Nvector.t -> unit

    (** This function reinitializes the integration of quadrature equations
        during the backward phase.

        @idas <node7#SECTION007211000000000000000> IDAQuadReInitB *)
    val reinit : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

    (** {2:tols Tolerances} *)

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

        @idas <node7#sss:quad_optional_input_B> IDASetQuadErrConB
        @idas <node7#sss:quad_optional_input_B> IDAQuadSStolerancesB
        @idas <node7#sss:quad_optional_input_B> IDAQuadSVtolerancesB *)
    val set_tolerances : ('d, 'k) bsession -> ('d, 'k) tolerance -> unit

    (** {2:quadout Output functions}

        @idas <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration *)

    (** Returns the backward quadrature solutions and time reached
        after a successful solver step. The given vectors are filled with
        values calculated during either {!backward_normal} or
        {!backward_one_step} and the value of the independent variable
        is returned.

        @idas <node7#sss:quad_get_b> IDAGetQuadB *)
    val get : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> float

    (** {2:get Querying the solver (optional output functions)}

        @idas <node7#sss:quad_optional_input_B> Optional input/output functions for backward quadrature integration *)

    (** Returns the number of calls to the backward quadrature right-hand
        side function.

        @idas <node5#ss:quad_optional_output> IDAGetQuadNumRhsEvals
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_rhs_evals : ('d, 'k) bsession -> int

    (** Returns the number of local error test failures due to quadrature
        variables.

        @idas <node5#ss:quad_optional_output> IDAGetQuadNumErrTestFails
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_num_err_test_fails : ('d, 'k) bsession -> int

    (** Returns the quadrature error weights at the current time.

        @idas <node5#ss:quad_optional_output> IDAGetQuadErrWeights
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
    val get_err_weights : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

    (** Returns quadrature-related statistics. These are the
        number of calls to the quadrature function ([nfqevals]) and the
        number of error test failures due to quadrature variables
        ([nqetfails]).

        @idas <node5#ss:quad_optional_output> IDAGetQuadStats
        @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
        @return ([nfqevals], [nqetfails]) *)
    val get_stats : ('d, 'k) bsession -> int * int
  end (* }}} *)

  (** Integrates a backward ODE system over an interval. The solver takes
      internal steps until it has reached or just passed the specified value.

      @idas <node7#sss:idasolveb> IDASolveB (IDA_NORMAL)
      @raise AdjointNotInitialized  The {!init} function has not been called.
      @raise NoBackwardProblem      The {!init_backward} function has not been called.
      @raise NoForwardCall          Neither {!forward_normal} nor {!forward_one_step} has been called.
      @raise Ida.IllInput           One of the inputs is invalid.
      @raise Ida.TooMuchWork        Could not reach [tout] in [mxstep] steps
      @raise Ida.TooMuchAccuracy    Could not satisfy the demanded accuracy
      @raise Ida.ErrFailure         Too many error test failures.
      @raise Ida.ConvergenceFailure Too many convergence test failures.
      @raise Ida.LinearSetupFailure Unrecoverable failure in linear solver setup function.
      @raise Ida.LinearSolveFailure Unrecoverable failure in linear solver solve function.
      @raise BadOutputTime          The requested output time is outside the interval over which the forward problem was solved.
      @raise ForwardReinitializationFailed Reinitialization of the forward problem failed at the first checkpoint (corresponding to the initial time of the forward problem).
      @raise ForwardFail              An error occurred during the integration of the forward problem. *)
  val backward_normal : ('d, 'k) Ida.session -> float -> unit

  (** Like {!backward_normal} but returns after one internal solver step.

      @idas <node7#sss:idasolveb> IDASolveB (IDA_ONE_STEP) *)
  val backward_one_step : ('d, 'k) Ida.session -> float -> unit

  (** Fills the given vectors, [yb] and [yb'], with the solution of the
      backward DAE problem at the returned time, interpolating if necessary.

      @idas <node7#sss:idasolveb> IDAGetB *)
  val get :
    ('d, 'k) bsession -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> float

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

      @idas <node5#ss:optional_dky> IDAGetDky
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
      @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
      @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
  val get_dky :
    ('d, 'k) bsession -> ('d, 'k) Nvector.t -> float -> int -> unit

  (** Fills the vector with the interpolated forward solution and its
      derivative at the given time during a backward simulation.

      @noidas <node5#ss:get_adjy> IdaGetAdjY *)
  val get_y : ('d, 'k) Ida.session -> ('d, 'k) Nvector.t
                -> ('d, 'k) Nvector.t  -> float -> unit

  (** Reinitializes the backward problem with new parameters and state
      values. The values of the independent variable, i.e., the simulation
      time, and the state variables and derivatives must be given. It is also
      possible to change the linear solver.

      @idas <node7#sss:idainitb> IDAReInitB
      @raise AdjointNotInitialized The {!init} function has not been called.
      @raise BadFinalTime      The final time is outside the interval over which the forward problem was solved. *)
  val reinit :
    ('d, 'k) bsession
    -> ?nlsolver:('d, 'k,
                  (('d, 'k) session) Sundials_NonlinearSolver.integrator)
                  Sundials_NonlinearSolver.nonlinear_solver
    -> ?lsolver:('d, 'k) session_linear_solver
    -> float
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) Nvector.t
    -> unit

  (** {3:calcic Initial Condition Calculation} *)

  (** Class components of the state vector as either algebraic or
      differential.
      These classifications are required by {!calc_ic}, {!calc_ic_sens}, and
      {!set_suppress_alg}. See also {!Ida.VarId}.

      @ida <node5#sss:optin_main> IDASetIdB *)
  val set_id : ('d, 'k) bsession -> ('d,'k) Nvector.t -> unit

  (** Indicates whether or not to ignore algebraic variables in the local error
      test. When ignoring algebraic variables ([true]), a [varid] vector must be
      specified either in the call or by a prior call to {!init} or {!set_id}.
      Suppressing local error tests for algebraic variables is {i discouraged}
      for DAE systems of index 1 and {i encouraged} for systems of index 2 or
      more.

      @ida <node5#sss:optin_main> IDASetId
      @idas <node7#ss:optional_input_b> IDASetSuppressAlgB *)
  val set_suppress_alg : ('d, 'k) bsession
                         -> ?varid:('d, 'k) Nvector.t -> bool -> unit

  (** Computes the algebraic components of the initial state and the
      differential components of the derivative vectors for certain
      index-one problems.
      The elements of $y_B$ marked algebraic and of {% $\dot{y}_B$%} marked
      differential are computed from the differential components of $y_B$,
      to satisfy the constraint
      {% $F(t_0, y_0, \dot{y}_0, y_\mathit{B0}, \dot{y}_\mathit{B0}) = 0$%}.
      The variable ids must be given in [~varid] or by a prior call to {!init} or
      {!set_id}.
      The call [calc_ic s ~yb ~yb' tbout0 y0 dy0] gives the first value
      at which a solution will be requested [tbout0], and the vectors of
      forward solutions [y0] and forward derivatives [dy0]. If given,
      the [~yb] and [~yb'] vectors are filled with the corrected backward
      states and derivatives. A {!reinit} is required before calling this
      function after {!forward_normal} or {!forward_one_step}.

      @idas <node7#sss:idacalcicB> IDACalcICB
      @idas <node7#ss:optional_output_b> IDAGetConsistentICB *)
  val calc_ic :
    ('d, 'k) bsession
    -> ?yb:('d, 'k) Nvector.t
    -> ?yb':('d, 'k) Nvector.t
    -> float                    (* tbout0 *)
    -> ('d, 'k) Nvector.t       (* y0 *)
    -> ('d, 'k) Nvector.t       (* dy0 *)
    -> unit

  (** Computes the algebraic components of the initial state and the
      differential components of the derivative vectors for certain
      index-one problems.
      The elements of $y_B$ marked algebraic and of {% $\dot{y}_B$%} marked
      differential are computed from the differential components of $y_B$,
      to satisfy the constraint
      {% $F(t_0, y_0, \dot{y}_0, y_\mathit{B0}, \dot{y}_\mathit{B0},
            s0, \dot{s}0) = 0$%}.
      The variable ids must be given in [~varid] or by a prior call to {!init} or
      {!set_id}.
      The call [calc_ic s ~yb ~yb' tbout0 y0 dy0 s0 ds0] gives the first
      value at which a solution will be requested [tbout0], the vectors of
      forward solutions [y0] and forward derivatives [dy0], and arrays of
      vectors of sensitivities and sensitivity derivatives.
      If given, the [~yb] and [~yb'] vectors are filled with the corrected
      backward states and derivatives. A {!reinit} is required before
      calling this function after {!forward_normal} or {!forward_one_step}.

      @idas <node7#sss:idacalcicB> IDACalcICBS
      @idas <node7#ss:optional_output_b> IDAGetConsistentICB *)
  val calc_ic_sens :
    ('d, 'k) bsession
    -> ?yb:('d, 'k) Nvector.t
    -> ?yb':('d, 'k) Nvector.t
    -> ?varid:('d, 'k) Nvector.t
    -> float                          (* tbout0 *)
    -> ('d, 'k) Nvector.t             (* y0 *)
    -> ('d, 'k) Nvector.t             (* dy0 *)
    -> ('d, 'k) Nvector.t array       (* s0 *)
    -> ('d, 'k) Nvector.t array       (* ds0 *)
    -> unit

  (** {2:set Modifying the solver (optional input functions)} *)

  (** Cancels the storage of sensitivity checkpointing data during forward
      solution (with {!forward_normal} or {!forward_one_step}).

      @idas <node7#SECTION00728000000000000000> IDAAdjSetNoSensi *)
  val set_no_sensitivity : ('d, 'k) Ida.session -> unit

  (** Specifies the maximum order of the linear multistep method.

      @idas <node7#ss:optional_input_b> IDASetMaxOrdB *)
  val set_max_ord : ('d, 'k) bsession -> int -> unit

  (** Specifies the maximum number of steps taken in attempting to reach
      a given output time.

      @idas <node7#ss:optional_input_b> IDASetMaxNumStepsB *)
  val set_max_num_steps : ('d, 'k) bsession -> int -> unit

  (** Specifies the initial step size.

      @idas <node7#ss:optional_input_b> IDASetInitStepB *)
  val set_init_step : ('d, 'k) bsession -> float -> unit

  (** Specifies an upper bound on the magnitude of the step size.

      @idas <node7#ss:optional_input_b> IDASetMaxStepB *)
  val set_max_step : ('d, 'k) bsession -> float -> unit

  (** Specifies a vector defining inequality constraints for each
      component of the solution vector [u].  See {!Sundials.Constraint}.

      @noidas <node> IDASetConstraintsB *)
  val set_constraints : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

  (** Disables constraint checking.

      @noidas <node> IDASetConstraints *)
  val clear_constraints : ('d, 'k) bsession -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the real and integer workspace sizes.

      @idas <node5#sss:optout_main> IDAGetWorkSpace
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
      @return ([real_size], [integer_size]) *)
  val get_work_space : ('d, 'k) bsession -> int * int

  (** Returns the cumulative number of internal steps taken by the solver.

      @idas <node5#sss:optout_main> IDAGetNumSteps
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_num_steps : ('d, 'k) bsession -> int

  (** Returns the number of calls to the backward residual function.

      @idas <node5#sss:optout_main> IDAGetNumRhsEvals
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_num_res_evals : ('d, 'k) bsession -> int

  (** Returns the number of calls made to the linear solver's setup function.

      @idas <node5#sss:optout_main> IDAGetNumLinSolvSetups
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_num_lin_solv_setups : ('d, 'k) bsession -> int

  (** Returns the number of local error test failures that have occurred.

      @idas <node5#sss:optout_main> IDAGetNumErrTestFails
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_num_err_test_fails : ('d, 'k) bsession -> int

  (** Returns the integration method order used during the last internal step.

      @idas <node5#sss:optout_main> IDAGetLastOrder
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_last_order : ('d, 'k) bsession -> int

  (** Returns the integration method order to be used on the next internal
      step.

      @idas <node5#sss:optout_main> IDAGetCurrentOrder
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_current_order : ('d, 'k) bsession -> int

  (** Returns the integration step size taken on the last internal step.

      @idas <node5#sss:optout_main> IDAGetLastStep
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_last_step : ('d, 'k) bsession -> float

  (** Returns the integration step size to be attempted on the next internal
      step.

      @idas <node5#sss:optout_main> IDAGetCurrentStep
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_current_step : ('d, 'k) bsession -> float

  (** Returns the the value of the integration step size used on the first
      step.

      @idas <node5#sss:optout_main> IDAGetActualInitStep
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_actual_init_step : ('d, 'k) bsession -> float

  (** Returns the the current internal time reached by the solver.

      @idas <node5#sss:optout_main> IDAGetCurrentTime
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_current_time : ('d, 'k) bsession -> float

  (** Returns a suggested factor by which the user's tolerances should be
      scaled when too much accuracy has been requested for some internal
      step.

      @idas <node5#sss:optout_main> IDAGetTolScaleFactor
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_tol_scale_factor : ('d, 'k) bsession -> float

  (** Returns the solution error weights at the current time.

      @idas <node5#sss:optout_main> IDAGetErrWeights
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
      @idas <node3#ss:ivp_sol> IVP solution (W_i) *)
  val get_err_weights : ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

  (** Returns the vector of estimated local errors.

      @idas <node5#sss:optout_main> IDAGetEstLocalErrors
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_est_local_errors :
    ('d, 'k) bsession -> ('d, 'k) Nvector.t -> unit

  (** Returns the integrator statistics as a group.

      @idas <node5#sss:optout_main> IDAGetIntegratorStats
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_integrator_stats : ('d, 'k) bsession -> Ida.integrator_stats

  (** Prints the integrator statistics on the given channel.

      @idas <node5#sss:optout_main> IDAGetIntegratorStats
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val print_integrator_stats : ('d, 'k) bsession -> out_channel -> unit

  (** Returns the number of nonlinear (functional or Newton) iterations
      performed.

      @idas <node5#sss:optout_main> IDAGetNumNonlinSolvIters
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_num_nonlin_solv_iters : ('d, 'k) bsession -> int

  (** Returns the number of nonlinear convergence failures that have occurred.

      @idas <node5#sss:optout_main> IDAGetNumNonlinSolvConvFails
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem *)
  val get_num_nonlin_solv_conv_fails : ('d, 'k) bsession -> int

  (** Returns both the numbers of nonlinear iterations performed [nniters] and
      nonlinear convergence failures [nncfails].

      @idas <node5#sss:optout_main> IDAGetNonlinSolvStats
      @idas <node7#ss:optional_output_b> IDAGetAdjIDABmem
      @return ([nniters], [nncfails]) *)
  val get_nonlin_solv_stats : ('d, 'k) bsession -> int * int

  (** {2:exceptions Exceptions} *)

  (** Adjoint sensitivity analysis was not initialized.

      @idas <node7#sss:idasolvef> IDA_NO_ADJ *)
  exception AdjointNotInitialized

  (** Neither {!forward_normal} nor {!forward_one_step} has been called.

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

  (** The final time was outside the interval over which the forward
      problem was solved.

      @idas <node7#sss:idainitb> IDA_BAD_TB0 *)
  exception BadFinalTime

end (* }}} *)

