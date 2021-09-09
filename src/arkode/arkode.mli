(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(***********************************************************************)
(* Parts of the comment text are taken directly from:                  *)
(*                                                                     *)
(*               User Documentation for Arkode v2.6.0                  *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(*              User Documentation for ARKODE v1.0.2                   *)
(*              Daniel R. Reynolds, David J. Gardner,                  *)
(*              Alan C. Hindmarsh, Carol S. Woodward                   *)
(*                     and Jean M. Sexton                              *)
(*                                                                     *)
(***********************************************************************)

(** Adaptive-step time integration for stiff, nonstiff, and mixed
    stiff/nonstiff systems of ODE initial value problems with zero-crossing
    detection.

    The interface is structured as follows.
    {ol
      {- {{:#generic}Generic constants and types}}
      {- {{:#MODULEARKStep}The ARKStep time-stepping module}}
      {- {{:#MODULEERKStep}The ERKStep time-stepping module}}
      {- {{:#MODULEMRIStep}The MRIStep time-stepping module}}
      {- {{:#exceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** {2:generic Generic constants and types}

    Types and values used by all time-stepping modules. *)

(** Values returned by the step functions. Failures are indicated by
    exceptions.

    @noarkode <node> ARKStepEvolve
    @noarkode <node> ERKStepEvolve
    @noarkode <node> MRIStepEvolve *)
type solver_result =
  | Success             (** The solution was advanced. {cconst ARK_SUCCESS} *)
  | RootsFound          (** A root was found. See {!ARKStep.get_root_info},
                            {!ERKStep.get_root_info}, or
                            {!MRIStep.get_root_info}.
                            {cconst ARK_ROOT_RETURN} *)
  | StopTimeReached     (** The stop time was reached. See
                            {!ARKStep.set_stop_time},
                            {!ERKStep.set_stop_time},
                            or {!MRIStep.set_stop_time}.
                            {cconst ARK_TSTOP_RETURN} *)

(** Summaries of integrator statistics. *)
type step_stats = {
    num_steps : int;
      (** Cumulative number of internal solver steps. *)
    actual_init_step : float;
      (** Integration step sized used on the first step. *)
    last_step : float;
      (** Integration step size of the last internal step. *)
    current_step : float;
      (** Integration step size to attempt on the next internal step. *)
    current_time : float
      (** Current internal time reached by the solver. *)
  }

(** Used to specify the nterpolation method used for output values and
    implicit method predictors. See, for example,
    {!MRIStep.set_interpolation_type}. *)
type interpolant_type =
    Hermite     (** Polynomial interpolants of Hermite form.
                    {cconst ARK_INTERP_HERMITE} *)
  | Lagrange    (** Polynomial interpolants of Lagrange form, for stiff
                    problems.
                    {cconst ARK_INTERP_LAGRANGE} *)

(** {3:arkodetol Tolerances} *)

(** Functions that set the multiplicative error weights for use in the weighted
    RMS norm. The call [efun y ewt] takes the dependent variable vector [y] and
    fills the error-weight vector [ewt] with positive values or raises
    {!Sundials.NonPositiveEwt}. Other exceptions are eventually propagated, but
    should be avoided ([efun] is not allowed to abort the solver). *)
type 'data error_weight_fun = 'data -> 'data -> unit

(** Tolerance specifications. *)
type ('data, 'kind) tolerance =
  | SStolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVtolerances of float * ('data, 'kind) Nvector.t
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFtolerances of 'data error_weight_fun
    (** Set the multiplicative error weights for the weighted RMS norm. *)

(** A default relative tolerance of 1.0e-4 and absolute tolerance of 1.0e-9. *)
val default_tolerances : ('data, 'kind) tolerance

(** {3:arkoderoots Roots} *)

(** Called by the solver to calculate the values of root functions. These
    ‘zero-crossings’ are used to detect significant events. The function is
    passed three arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
    - [gout], a vector for storing the value of $g(t, y)$.

    {warning [y] and [gout] should not be accessed after the function has
             returned.}

    @noarkode <node> ARKRootFn *)
type 'd rootsfn = float -> 'd -> RealArray.t -> unit

(** A convenience value for signalling that there are no roots to monitor. *)
val no_roots : (int * 'd rootsfn)

(** {3:arkodeadapt Adaptivity} *)

type adaptivity_args = {
    h1 : float;  (** the current step size, {% $t_m - t_{m-1}$%}. *)
    h2 : float;  (** the previous step size, {% $t_{m-1} - t_{m-2}$%}. *)
    h3 : float;  (** the step size {% $t_{m-2} - t_{m-3}$%}. *)
    e1 : float;  (** the error estimate from the current step, {% $m$%}. *)
    e2 : float;  (** the error estimate from the previous step, {% $m-1$%}. *)
    e3 : float;  (** the error estimate from the step {% $m-2$%}. *)
    q  : int;    (** the global order of accuracy for the integration method. *)
    p  : int;    (** the global order of accuracy for the embedding. *)
  }

(** A function implementing a time step adaptivity algorithm that chooses an
    $h$ that satisfies the error tolerances. The call [hnew = adapt_fn t y args]
    has as arguments
    - [t], the value of the independent variable,
    - [y], the value of the dependent variable vector {% $y(t)$%}, and
    - [args], information on step sizes, error estimates, and accuracies.
    and returns the next step size [hnew]. The function should raise an
    exception if it cannot set the next step size. The step size should be the
    maximum value where the error estimates remain below 1.

    This function should focus on accuracy-based time step estimation; for
    stability based time steps, {!ARKStep.set_stability_fn} and
    {!ERKStep.set_stability_fn} should be used.

    @noarkode <node> ARKAdaptFn *)
type 'd adaptivity_fn = float -> 'd -> adaptivity_args -> float

(** Parameters for the standard adaptivity algorithms.
    There are two:
    - [adaptivity_ks], the [k1], [k2], and [k3] parameters, or [None]
      to use the defaults, and
    - [adaptivity_method_order], [true] specifies the method order of
      accuracy $q$ and [false] specifies the embedding order of
      accuracy $p$. *)
type adaptivity_params = {
    ks : (float * float * float) option;
    method_order : bool;
  }

(** Asymptotic error control algorithms.

    @noarkode <node> Asymptotic error control *)
type 'd adaptivity_method =
  | PIDcontroller of adaptivity_params
        (** The default time adaptivity controller. *)
  | PIcontroller of adaptivity_params
        (** Uses the two most recent step sizes. *)
  | Icontroller of adaptivity_params
        (** Standard time adaptivity control algorithm. *)
  | ExplicitGustafsson of adaptivity_params
        (** Primarily used with explicit RK methods. *)
  | ImplicitGustafsson of adaptivity_params
        (** Primarily used with implicit RK methods. *)
  | ImExGustafsson of adaptivity_params
        (** An ImEx version of the two preceding controllers. *)
  | AdaptivityFn of 'd adaptivity_fn
        (** A custom time-step adaptivity function. *)

(** {3:arkodecallbacks Callback functions} *)

(** Right-hand side functions for calculating ODE derivatives. They are passed
    three arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
    - [y'], a vector for storing the value of $f(t, y)$.

    Within the function, raising a {!Sundials.RecoverableFailure} exception
    indicates a recoverable error. Any other exception is treated as an
    unrecoverable error.

    {warning [y] and [y'] should not be accessed after the function returns.}

    @noarkode <node> ARKRhsFn *)
type 'd rhsfn = float -> 'd -> 'd -> unit

(** A function that predicts the maximum stable step size for the explicit
    portions of an ImEx ODE system. The call [hstab = stab_fn t y]
    has as arguments
    - [t], the value of the independent variable, and
    - [y], the value of the dependent variable vector {% $y(t)$%}.
    and returns the absolute value of the maximum stable step size [hstab].
    Returning {% $\mathtt{hstab}\le 0$%} indicates that there is no explicit
    stability restriction on the time step size. The function should raise
    an exception if it cannot set the next step size.

    @noarkode <node> ARKExpStabFn *)
type 'd stability_fn = float -> 'd -> float

(** Called to resize a vector to match the dimensions of another. The call
    [resizefn y ytemplate] must resize [y] to match the size of [ytemplate].

    {warning [y] and [ytemplate] should not be accessed after the function
             has returned.}

    @noarkode <node> ARKVecResizeFn *)
type 'd resize_fn = 'd -> 'd -> unit

(** A function to process the results of each timestep solution.
    The arguments are
    - [t], the value of the independent variable, and
    - [y], the value of the dependent variable vector {% $y(t)$%}.

    @noarkode <node> ARKPostprocessStepFn *)
type 'd postprocess_step_fn = float -> 'd -> unit

(** Butcher tables *)
module ButcherTable : sig (* {{{ *)

  (** A butcher table.

      {% $\begin{array}{ c|c}
          c & A \\ \hline
          q & b \\
          p & \tilde{b}
         \end{array}$ %}

      Instantiated according to [stages]. For example, when [stages = 3]:
        {% $\begin{array}{ c|ccc}
            c_1 & A_{1,1} & A_{1,2} & A_{1,3} \\
            c_2 & A_{2,1} & A_{2,2} & A_{2,3} \\
            c_3 & A_{3,1} & A_{3,2} & A_{3,3} \\ \hline
            q   & b_1     & b_2     & b_3 \\
            p   & \widetilde{b_1} & \widetilde{b_2} & \widetilde{b_3}
            \end{array}$ %}

      @noarkode <node> ARKodeButcherTable *)
  type t = {
      method_order : int;          (** Method order of accuracy ({% $q$ %}). *)
      embedding_order : int;       (** Embedding order of accuracy ({% $p$ %}). *)
      stages : int;                (** Number of stages ({% $s$ %}). *)
      stage_values : Sundials.RealArray2.t; (** Matrix ([stages * stages])
                                                of coefficients ({% $A$ %}) *)
      stage_times : RealArray.t;   (** Array (of length [stages]) of
                                       stage times ({% $c$ %}). *)
      coefficients : RealArray.t;  (** Array (of length [stages]) of
                                       solution coefficients ({% $b$ %}). *)
      bembed : RealArray.t option; (** Optional array (of length [stages]) of
                                       embedding coefficients
                                       ({% $\tilde{b}$ %}). *)
    }

  (** Explicit Butcher tables

      @noarkode <node> Explicit Butcher tables *)
  type erk_table =
    | HeunEuler_2_1_2       (** Default 2nd order explicit method. *)
    | BogackiShampine_4_2_3 (** Default 3rd order explicit method. *)
    | ARK_4_2_3_Explicit    (** Explicit portion of default 3rd order additive
                                method. *)
    | Zonneveld_5_3_4       (** Default 4th order explicit method. *)
    | ARK_6_3_4_Explicit    (** Explicit portion of default 3rd order additive
                                method. *)
    | SayfyAburub_6_3_4     (** From Sayfy and Aburub 2002. *)
    | CashKarp_6_4_5        (** Default 5th order explicit method. *)
    | Fehlberg_6_4_5        (** From Fehlberg 1969. *)
    | DormandPrince_7_4_5   (** From Dormand Prince 1980. *)
    | ARK_8_4_5_Explicit    (** Explicit portion of default 5th order additive
                                method. *)
    | Verner_8_5_6          (** Default 6th order explicit method. *)
    | Fehlberg_13_7_8       (** Default 8th order explicit method. *)
    | Knoth_Wolke_3_3       (** Default 3rd order slow and fast method
                                (Sundials >= 4.0.0). *)
    | ARK_7_3_4_Explicit    (** Explicit portion of the 4th order additive
                                method from Kennedy and Carpenter 2019.
                                (Sundials >= 5.0.0). *)
    | ARK_8_4_5b_Explicit   (** Explicit portion of the 5th order additive
                                method from Kennedy and Carpenter 2019.
                                (Sundials >= 5.0.0). *)

  (** Implicit Butcher tables

      @noarkode <node> Implicit Butcher tables *)
  type dirk_table =
    | SDIRK_2_1_2           (** Default 2nd order implicit method. *)
    | Billington_3_2_3      (** From Billington 1983. *)
    | TRBDF2_3_2_3          (** From Billington 1985. *)
    | Kvaerno_4_2_3         (** From Kvaerno 2004. *)
    | ARK_4_2_3_Implicit    (** Default 3rd order implicit method and the
                                implicit portion of the default 3rd order
                                additive method. *)
    | Cash_5_2_4            (** From Cash 1979. *)
    | Cash_5_3_4            (** From Cash 1979. *)
    | SDIRK_5_3_4           (** Default 4th order implicit method. *)
    | Kvaerno_5_3_4         (** From Kvaerno 2004. *)
    | ARK_6_3_4_Implicit    (** Implicit portion of the default 4th order
                                additive method. *)
    | Kvaerno_7_4_5         (** From Kvaerno 2004. *)
    | ARK_8_4_5_Implicit    (** Default 5th order method and the implicit
                                portion of the default 5th order additive
                                method. *)
    | ARK_7_3_4_Implicit    (** Implicit portion of the 4th order additive
                                method from Kennedy and Carpenter 2019.
                                (Sundials >= 5.0.0). *)
    | ARK_8_4_5b_Implicit   (** Implicit portion of the 5th order additive
                                method from Kennedy and Carpenter 2019.
                                (Sundials >= 5.0.0). *)

  (** Additive Butcher tables

      @noarkode <node> Additive Butcher tables *)
  type ark_table =
    | ARK_4_2_3             (** 3rd-order pair combining
                                BogackiShampine_4_2_3 and ARK_4_2_3_Implicit. *)
    | ARK_6_3_4             (** 4th-order pair combining
                                ARK_6_3_4_Explicit and ARK_6_3_4_Implicit. *)
    | ARK_8_4_5             (** 5th-order pair combining
                                ARK_8_4_5_Explicit and ARK_8_4_5_Implicit. *)

  (** Retrieves an explicit Butcher table.

      @noarkode <node> ARKodeButcherTable_LoadERK
      @since 4.0.0 *)
  val load_erk  : erk_table -> t

  (** Retrieves a diagonally-implicit Butcher table.

      @noarkode <node> ARKodeButcherTable_LoadDIRK
      @since 4.0.0 *)
  val load_dirk : dirk_table -> t

  (** Writes a Butcher table to a file.

      @noarkode <node> ARKodeButcherTable_Write
      @since 4.0.0 *)
  val write : t -> Logfile.t -> unit

  (** Indicates that a check on the analytic order of accuracy failed. *)
  exception ButcherTableCheckFailed

  (** Determines the analytic order of accuracy for a Butcher table.
      The analytic (necessary) conditions are checked up to order 6. For
      orders greater than 6, the Butcher simplifying (sufficient) assumptions
      are used. In the call, [(q, p, warn) = check_order bt],
      - [q] is the measured order of accuracy,
      - [p] is the measured order of accuracy for the embedding if applicable,
      - [warn] is true if the table values are lower than the measured values,
        or the measured values reach the maximum order possible with this
        function and the values of [q] and [p] in the provided table are higher.

      The logfile, if given, is used to print results.

      @noarkode <node> ARKodeButcherTable_CheckOrder
      @raise ButcherTableCheckFailed Table values are higher than measured ones
      @since 4.0.0 *)
  val check_order : ?outfile:Logfile.t -> t -> int * int option * bool

  (** Determines the analytic order of accuracy for a pair of Butcher tables.
      The analytic (necessary) conditions are checked up to order 6. For
      orders greater than 6, the Butcher simplifying (sufficient) assumptions
      are used. In the call, [(q, p, warn) = check_order b1 b2],
      - [q] is the measured order of accuracy,
      - [p] is the measured order of accuracy for the embedding if applicable,
      - [warn] is true if the table values are lower than the measured values,
        or the measured values reach the maximum order possible with this
        function and the values of [q] and [p] in the provided table are higher.

      The logfile, if given, is used to print results.

      @noarkode <node> ARKodeButcherTable_CheckARKOrder
      @since 4.0.0 *)
  val check_ark_order : ?outfile:Logfile.t -> t -> t -> int * int option * bool

end (* }}} *)

(** ARKStep Time-Stepping Module for ODE systems in split, linearly-implicit
    form.

    This module solves problems of the form
    {% $M\dot{y} = f_E(t, y) + f_I(t, y)$%}, {% $y(t_0) = y_0$%}.

    Its interface is structured as follows.
    {ol
      {- {{:#linear}Linear and mass matrix solvers}}
      {- {{:#tols}Tolerances}}
      {- {{:#solver}Solver initialization and use}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}}

    @noarkode <node> Using ARKStep for C and C++ Applications *)
module ARKStep : sig (* {{{ *)

  (** A session with the ARKStep time-stepping solver.

      An example session with ARKStep ({openfile arkode_ark_skel.ml}): {[
#include "../../examples/ocaml/skeletons/arkode_ark_skel.ml"
      ]}

      @noarkode <node> Skeleton of main program *)
  type ('d, 'k) session = ('d, 'k, Arkode_impl.arkstep) Arkode_impl.session

  (** Alias for sessions based on serial nvectors. *)
  type 'k serial_session = (Nvector_serial.data, 'k) session
                           constraint 'k = [>Nvector_serial.kind]

  (** {2:linear Linear and mass matrix solvers} *)

  (** Linear solvers used by Arkode.

      @noarkode <node> Linear Solver Specification Functions *)
  type ('data, 'kind) linear_solver = ('data, 'kind) Arkode_impl.linear_solver

  (** Alias for linear solvers that are restricted to serial nvectors. *)
  type 'kind serial_linear_solver =
    (Nvector_serial.data, 'kind) linear_solver
    constraint 'kind = [>Nvector_serial.kind]

  (** Workspaces with three temporary vectors. *)
  type 'd triple = 'd * 'd * 'd

  (** Arguments common to Jacobian callback functions.

      @noarkode <node> ARKLsJacFn
      @noarkode <node> ARKLsJacTimesVecFn
      @noarkode <node> ARKLsPrecSolveFn
      @noarkode <node> ARKLsPrecSetupFn *)
  type ('t, 'd) jacobian_arg = ('t, 'd) Arkode_impl.jacobian_arg =
    {
      jac_t   : float;        (** The independent variable. *)
      jac_y   : 'd;           (** The dependent variable vector. *)
      jac_fy  : 'd;           (** The derivative vector {% $f_I(t, y)$%}. *)
      jac_tmp : 't            (** Workspace data. *)
    }

  (** Direct Linear Solvers operating on dense, banded, and sparse matrices.

      @noarkode <node> Linear solver specification functions
      @noarkode <node> Dense/band direct linear solvers optional input functions
      @noarkode <node> Dense/band direct linear solvers optional output functions
  *)
  module Dls : sig (* {{{ *)
    include module type of Sundials_LinearSolver.Direct

    (** Callback functions that compute dense approximations to a Jacobian
        matrix. In the call [jac arg jm], [arg] is a {!jacobian_arg}
        with three work vectors and the computed Jacobian must be stored
        in [jm].

        The callback should load the [(i,j)]th entry of [jm] with
        {% $\partial (f_I)_i/\partial y_j$%}, i.e., the partial derivative
        of the [i]th implicit equation with respect to the [j]th variable,
        evaluated at the values of [t] and [y] obtained from [arg]. Only
        nonzero elements need be loaded into [jm].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jm] should
                 be accessed after the function has returned.}

        @noarkode <node> ARKLsJacFn *)
    type 'm jac_fn = (RealArray.t triple, RealArray.t) jacobian_arg
                     -> 'm -> unit

  (** Function to compute the linear system matrix {% $A = M - \gamma J$ %}
      or an approximation of it. Offers an alternative to evaluating the
      Jacobian of the right-hand-side function.

      In addition to those shared with the Jacobian function, the arguments of
      this function are
      - [a], storage for the computed linear system matrix,
      - [m], the current mass matrix if {% $M \neq I$ %},
      - [jok], indicates whether the Jacobian-related data needs to be
               updated, and
      - [gamma], the scalar in the formula above.

      The function should return true only if the Jacobian data was
      recomputed.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the Jacobian argument elements nor the matrices
               [a] and [m] should be accessed after the function has
               returned.}

      @since 5.0.0
      @nocvode <node> ARKLsLinSysFn *)
  type 'm linsys_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg
    -> 'm
    -> 'm option
    -> bool
    -> float
    -> bool

    (** Create an Arkode-specific linear solver from a Jacobian approximation
        function and a generic direct linear solver.
        The Jacobian approximation function is optional for dense and banded
        solvers (if not given an internal difference quotient approximation
        is used), but must be provided for other solvers (or [Invalid_argument]
        is raised).

        @noarkode <node> ARKStepSetLinearSolver
        @noarkode <node> ARKStepSetJacFn
        @noarkode <node> ARKStepSetLinSysFn *)
    val solver :
      ?jac:'m jac_fn ->
      ?linsys:'m linsys_fn ->
      ('m, RealArray.t, 'kind, [>`Dls]) LinearSolver.t ->
      'kind serial_linear_solver

    (** {3:arkdlsstats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by a direct
        linear solver.

        @noarkode <node> ARKStepGetLinWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : 'k serial_session -> int * int

    (** Returns the number of calls made by a direct linear solver to the
        Jacobian approximation function.

        @noarkode <node> ARKStepGetNumJacEvals *)
    val get_num_jac_evals : 'k serial_session -> int

    (** Returns the number of calls to the right-hand side callback due to
        the finite difference Jacobian approximation.

        @noarkode <node> ARKStepGetNumLinRhsEvals *)
    val get_num_lin_rhs_evals : 'k serial_session -> int

  end (* }}} *)

  (** Scaled Preconditioned Iterative Linear Solvers.

      @noarkode <node> Linear solver specification functions
      @noarkode <node> Iterative linear solvers optional input functions.
      @noarkode <node> Iterative linear solvers optional output functions. *)
  module Spils : sig (* {{{ *)
    include module type of Sundials_LinearSolver.Iterative

    (** {3:arkspilsprecond Preconditioners} *)

    (** Arguments passed to the preconditioner solver function.

        @noarkode <node> ARKLsPrecSolveFn *)
    type 'd prec_solve_arg =
      {
        rhs   : 'd;         (** Right-hand side vector of the linear system. *)
        gamma : float;      (** Scalar $\gamma$ in the Newton
                                matrix given by $A = M - \gamma J$. *)
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
        the Newton matrix {% $A = M - \gamma J$%} where
        {% $J = \frac{\partial f_I}{\partial y}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac], [arg], and [z] should not
                 be accessed after the function has returned.}

        @noarkode <node> ARKLsPrecSolveFn *)
    type 'd prec_solve_fn =
      (unit, 'd) jacobian_arg
      -> 'd prec_solve_arg
      -> 'd
      -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related data
        needed by {!prec_solve_fn}. In the call [prec_setup_fn jac jok gamma],
        [jac] is a {!jacobian_arg} with three work vectors, [jok] indicates
        whether any saved Jacobian-related data can be reused with the current
        value of [gamma], and [gamma] is the scalar $\gamma$ in the Newton
        matrix {% $A = M - \gamma J$%} where $J$ is the Jacobian matrix.
        A function should return [true] if Jacobian-related data was updated
        and [false] if saved data was reused.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac] should not be accessed after the
                 function has returned.}

        @noarkode <node> ARKStepSetPreconditioner
        @noarkode <node> ARKLsPrecSetupFn *)
    type 'd prec_setup_fn =
      (unit, 'd) jacobian_arg
      -> bool
      -> float
      -> bool

    (** Specifies a preconditioner, including the type of preconditioning
        (none, left, right, or both) and callback functions.
        The following functions and those in {!Banded} and {!Arkode_bbd}
        construct preconditioners.

        The {!prec_solve_fn} is usually mandatory. The {!prec_setup_fn} can be
        omitted if not needed.

        @noarkode <node> ARKStepSetPreconditioner
        @noarkode <node> ARKLsPrecSetupFn
        @noarkode <node> ARKLsPrecSolveFn *)
    type ('d,'k) preconditioner = ('d,'k) Arkode_impl.SpilsTypes.preconditioner

    (** No preconditioning.  *)
    val prec_none : ('d, 'k) preconditioner

    (** Left preconditioning. {% $(P^{-1}A)x = P^{-1}b$ %}. *)
    val prec_left :
      ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

    (** Right preconditioning. {% $(AP^{-1})Px = b$ %}. *)
    val prec_right :
      ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

    (** Left and right preconditioning.
        {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)
    val prec_both :
      ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> ('d, 'k) preconditioner

    (** Banded preconditioners.  *)
    module Banded : sig (* {{{ *)

      (** The range of nonzero entries in a band matrix. *)
      type bandrange =
        { mupper : int; (** The upper half-bandwidth. *)
          mlower : int; (** The lower half-bandwidth. *) }

      (** A band matrix {!preconditioner} based on difference quotients.
          The call [prec_left br] instantiates a left preconditioner which
          generates a banded approximation to the Jacobian with [br.mlower]
          sub-diagonals and [br.mupper] super-diagonals.

          NB: Banded preconditioners may not be used for problems involving
          a non-identity mass matrix.

          @noarkode <node> ARKBandPrecInit *)
      val prec_left : bandrange -> (Nvector_serial.data,
                                    [>Nvector_serial.kind]) preconditioner

      (** Like {!prec_left} but preconditions from the right.

          @noarkode <node> ARKBandPrecInit *)
      val prec_right :
           bandrange -> (Nvector_serial.data,
                         [>Nvector_serial.kind]) preconditioner

      (** Like {!prec_left} but preconditions from both sides.

          @noarkode <node> ARKBandPrecInit *)
      val prec_both :
           bandrange -> (Nvector_serial.data,
                         [>Nvector_serial.kind]) preconditioner

      (** {4:stats Banded statistics} *)

      (** Returns the sizes of the real and integer workspaces used by the
          banded preconditioner module.

          @noarkode <node> ARKBandPrecGetWorkSpace
          @return ([real_size], [integer_size]) *)
      val get_work_space : 'k serial_session -> int * int

      (** Returns the number of calls to the right-hand side callback for
          the difference banded Jacobian approximation. This counter is only
          updated if the default difference quotient function is used.

          @noarkode <node> ARKBandPrecGetNumRhsEvals *)
      val get_num_rhs_evals : 'k serial_session -> int
    end (* }}} *)

    (** {3:arkspilslsolvers Solvers} *)

    (** Callback functions that preprocess or evaluate Jacobian-related data
        needed by the jac_times_vec_fn. In the call [jac_times_setup_fn arg],
        [arg] is a {!jacobian_arg} with no work vectors.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The elements of [arg] should not be accessed after the
                 function has returned.}

        @noarkode <node> ARKLsJacTimesSetupFn *)
    type 'd jac_times_setup_fn =
      (unit, 'd) jacobian_arg
      -> unit

    (** Callback functions that compute the Jacobian times a vector. In the
        call [jac_times_vec_fn arg v jv], [arg] is a {!jacobian_arg} with one
        work vector, [v] is the vector multiplying the Jacobian, and [jv] is
        the vector in which to store the
        result—{% $\mathtt{jv} = J\mathtt{v}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor [v] or [jv] should be
                 accessed after the function has returned.}

        @noarkode <node> ARKLsJacTimesVecFn *)
    type 'd jac_times_vec_fn =
      ('d, 'd) jacobian_arg
      -> 'd (* v *)
      -> 'd (* Jv *)
      -> unit

    (** Create an Arkode-specific linear solver from a generic iterative
        linear solver.

        The [jac_times_rhs] argument specifies an alternative right-hand-side
        function for use in the internal Jacobian-vector product difference
        quotient approximation. It is incorrect to specify both this argument
        and [jac_times_vec].

        NB: a [jac_times_setup_fn] is not supported in
            {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

        NB: a [jac_times_rhs] function is not supported in
            {{!Sundials_Config.sundials_version}Config.sundials_version} < 5.3.0.

        @since 3.0.0
        @noarkode <node> ARKStepSetLinearSolver
        @noarkode <node> ARKStepSetJacTimes
        @noarkode <node> ARKStepSetJacTimesRhsFn *)
    val solver :
      ('m, 'd, 'k, [>`Iter]) LinearSolver.t
      -> ?jac_times_vec:'d jac_times_setup_fn option * 'd jac_times_vec_fn
      -> ?jac_times_rhs:'d rhsfn
      -> ('d, 'k) preconditioner
      -> ('d, 'k) linear_solver

    (** {3:arkspilsset Solver parameters} *)

    (** Sets the maximum number of time steps to wait before recomputation of
        the Jacobian or recommendation to update the preconditioner. If the
        integer argument is less than or equal to 0, a default value of 50 is
        used.

        @noarkode <node5> ARKStepSetJacEvalFrequency
        @since 4.0.0 *)
    val set_jac_eval_frequency : ('d, 'k) session -> int -> unit

    (** Enables or disables scaling of the linear system solution to account
        for a change in {% $\gamma$ %} in the linear system.
        Linear solution scaling is enabled by default when a matrix-based
        linear solver is attached.

        @noarkode <node5> ARKStepSetLinearSolutionScaling
        @since 5.2.0 *)
    val set_linear_solution_scaling : ('d, 'k) session -> bool -> unit

    (** Sets the factor by which the Krylov linear solver's convergence test
        constant is reduced from the Newton iteration test constant.
        This factor must be >= 0; passing 0 specifies the default (0.05).

        @noarkode <node> ARKStepSetEpsLin *)
    val set_eps_lin : ('d, 'k) session -> float -> unit

    (** Sets the factor for converting from the integrator tolerance (WRMS
        norm) to the linear solver tolerance (L2 norm). That is,
        {% $\mathit{tol}_{\mathsf{L2}} =
            \mathit{fact}\cdot\mathit{tol}_{\mathsf{WRMS}}$ %}.
        The given value is used directly if it is greater than zero.
        If it is zero (the default), then the square root of the state
        vector length is used.
        If it is less than zero, then the square root of the dot product of a
        state vector full of ones with itself is used.

        @since 5.4.0
        @noarkode <node> ARKStepSetLSNormFactor *)
    val set_ls_norm_factor : ('d, 'k) session -> float -> unit

    (** {3:arkspilsstats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the spils
        linear solver.

        @noarkode <node> ARKStepGetLinWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space       : ('d, 'k) session -> int * int

    (** Returns the cumulative number of linear iterations.

        @noarkode <node> ARKStepGetNumLinIters *)
    val get_num_lin_iters    : ('d, 'k) session -> int

    (** Returns the cumulative number of linear convergence failures.

        @noarkode <node> ARKStepGetNumLinConvFails *)
    val get_num_lin_conv_fails   : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the setup function with
        [jok=false].

        @noarkode <node> ARKStepGetNumPrecEvals *)
    val get_num_prec_evals   : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the preconditioner solve
        function.

        @noarkode <node> ARKStepGetNumPrecSolves *)
    val get_num_prec_solves  : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        setup function.

        @since 3.0.0
        @noarkode <node> ARKStepGetNumJTSetupEvals *)
    val get_num_jtsetup_evals : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        function.

        @noarkode <node> ARKStepGetNumJtimesEvals *)
    val get_num_jtimes_evals : ('d, 'k) session -> int

    (** Returns the number of calls to the right-hand side callback for
        finite difference Jacobian-vector product approximation. This counter is
        only updated if the default difference quotient function is used.

        @noarkode <node> ARKStepGetNumLinRhsEvals *)
    val get_num_lin_rhs_evals    : ('d, 'k) session -> int

    (** {3:arkspilslowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change preconditioner functions. These low-level functions are
        provided for experts who want to avoid resetting internal counters
        and other associated side-effects. *)

    (** Change the preconditioner functions.

        @noarkode <node> ARKStepSetPreconditioner
        @noarkode <node> ARKLsPrecSolveFn
        @noarkode <node> ARKLsPrecSetupFn *)
    val set_preconditioner :
      ('d, 'k) session
      -> ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> unit

    (** Change the Jacobian-times-vector function.

        @noarkode <node> ARKStepSetJacTimes
        @noarkode <node> ARKLsJacTimesVecFn *)
    val set_jac_times :
      ('d, 'k) session
      -> ?jac_times_setup:'d jac_times_setup_fn
      -> 'd jac_times_vec_fn
      -> unit

    (** Remove a Jacobian-times-vector function and use the default
        implementation.

        @noarkode <node> ARKStepSetJacTimes *)
    val clear_jac_times : ('d, 'k) session -> unit

  end (* }}} *)

  (** Mass Matrix Solvers

      @noarkode <node> Mass matrix solver *)
  module Mass : sig (* {{{ *)

    (** Mass matrix solvers used by Arkode.

        @noarkode <node> Mass matrix solver specification functions *)
    type ('data, 'kind) solver = ('data, 'kind) Arkode_impl.MassTypes.solver

    (** Alias for mass matrix solvers that are restricted to serial nvectors. *)
    type 'kind serial_solver = (Nvector_serial.data, 'kind) solver
                               constraint 'kind = [>Nvector_serial.kind]

    (** {3:arkdlsmass Direct mass matrix solvers} *)
    module Dls : sig (* {{{ *)
      include module type of Sundials_LinearSolver.Direct

      (** Functions that compute a mass matrix (or an approximation of one).
          In the call [mass t work m],
          - [t] is the independent variable,
          - [work] is workspace data, and
          - [m] is the output dense mass matrix.

          The callback should load the [N] by [N] matrix [m] with an
          approximation to the mass matrix {% $M(t)$%}. Only nonzero elements
          need be loaded into [m].

          Raising {!Sundials.RecoverableFailure} indicates a recoverable
          error. Any other exception is treated as an unrecoverable error.

          {warning Neither the elements of [work] nor the matrix [m]
                   should be accessed after the function has returned.}

          @noarkode <node> ARKLsMassFn *)
      type 'm mass_fn = float -> RealArray.t triple -> 'm -> unit

      (** Create an Arkode-specific mass linear solver from a mass-matrix
          constructor function and a generic dense linear solver. The
          boolean argument indicates whether the mass matrix depends on the
          independent variable [t], if not it is only computed and factored
          once.

          NB: The boolean argument is ignored in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

          @noarkode <node> ARKStepSetMassLinearSolver
          @noarkode <node> ARKStepSetMassFn *)
      val solver :
        'm mass_fn
        -> bool
        -> ('m, RealArray.t, 'kind, [>`Dls]) LinearSolver.t
        -> 'kind serial_solver

      (** {3:arkdlsmassstats Solver statistics} *)

      (** Returns the sizes of the real and integer workspaces used by a
          direct linear mass matrix solver.

          @noarkode <node> ARKStepGetMassWorkSpace
          @return ([real_size], [integer_size]) *)
      val get_work_space : 'k serial_session -> int * int

      (** Returns the number of calls made to the mass matrix solver setup
          routine.

          @since 3.0.0
          @noarkode <node> ARKStepGetNumMassSetups *)
      val get_num_setups : 'k serial_session -> int

      (** Returns the number of calls made to the mass matrix {em matvec
          setup} routine.

          @since 5.0.0
          @noarkode <node> ARKStepGetNumMassMultSetups *)
      val get_num_mult_setups : 'k serial_session -> int

      (** Returns the number of calls made to the mass matrix solver solve
          routine.

          @noarkode <node> ARKStepGetNumMassSolves *)
      val get_num_solves : 'k serial_session -> int

      (** Returns the number of calls made to the mass matrix-times-vector
          routine.

          @since 3.0.0
          @noarkode <node> ARKStepGetNumMassMult *)
      val get_num_mult : 'k serial_session -> int

    end (* }}} *)

    (** {3:arkspilsmass Iterative mass matrix solvers} *)
    module Spils : sig (* {{{ *)
      include module type of Sundials_LinearSolver.Iterative

      (** {3:arkspilsmassprecond Preconditioners} *)

      (** Arguments passed to the mass matrix preconditioner solver function.

          @noarkode <node> ARKLsMassPrecSolveFn *)
      type 'd prec_solve_arg =
        {
          rhs   : 'd;      (** Right-hand side vector of the linear system. *)
          delta : float;   (** Input tolerance for iterative methods. *)
          left  : bool;    (** [true] for left preconditioning and
                               [false] for right preconditioning. *)
        }

      (** Callback functions that solve a linear mass matrix system involving
          a preconditioner matrix. In the call [prec_solve_fn t arg z], [t]
          is the independent variable, [arg] is a {!prec_solve_arg} that
          specifies the linear system, and [z] is computed to solve
          {% $P\mathtt{z} = \mathtt{arg.rhs}$%}. {% $P$%} is a left or right
          preconditioning matrix, if preconditioning is done on both sides,
          the product of the two preconditioner matrices should approximate
          {% $M$%}.

          Raising {!Sundials.RecoverableFailure} indicates a recoverable
          error. Any other exception is treated as an unrecoverable error.

          {warning The elements of [arg] and [z] should not
                   be accessed after the function has returned.}

          @noarkode <node> ARKLsMassPrecSolveFn *)
      type 'd prec_solve_fn =
           float
        -> 'd prec_solve_arg
        -> 'd
        -> unit

      (** Callback functions that preprocess or evaluate mass matrix-related
          data needed by {!prec_solve_fn}. The argument gives the independent
          variable [t].

          Raising {!Sundials.RecoverableFailure} indicates a recoverable
          error. Any other exception is treated as an unrecoverable error.

          @noarkode <node> ARKLsMassPrecSetupFn *)
      type 'd prec_setup_fn =
           float
        -> unit

      (** Specifies a preconditioner, including the type of preconditioning
          (none, left, right, or both) and callback functions.
          The following functions construct preconditioners.

          @noarkode <node> ARKStepSetMassPreconditioner
          @noarkode <node> ARKLsMassPrecSetupFn
          @noarkode <node> ARKLsMassPrecSolveFn *)
      type ('d, 'k) preconditioner =
        ('d, 'k) Arkode_impl.MassTypes.Iterative'.preconditioner

      (** No preconditioning.  *)
      val prec_none : ('d, 'k) preconditioner

      (** Left preconditioning. {% $(P^{-1}M)x = P^{-1}b$ %}. *)
      val prec_left :
        ?setup:'d prec_setup_fn
        -> 'd prec_solve_fn
        -> ('d, 'k) preconditioner

      (** Right preconditioning. {% $(MP^{-1})Px = b$ %}. *)
      val prec_right :
        ?setup:'d prec_setup_fn
        -> 'd prec_solve_fn
        -> ('d, 'k) preconditioner

      (** Left and right preconditioning.
          {% $(P_L^{-1}MP_R^{-1})P_Rx = P_L^{-1}b$ %} *)
      val prec_both :
        ?setup:'d prec_setup_fn
        -> 'd prec_solve_fn
        -> ('d, 'k) preconditioner

      (** {3:arkspilsmasslsolvers Solvers} *)

      (** Callback functions that preprocess or evaluate Jacobian-related
          data needed by the mass_times_vec_fn. The argument gives the
          independent variable [t].

          Raising {!Sundials.RecoverableFailure} indicates a recoverable
          error. Any other exception is treated as an unrecoverable error.

          @noarkode <node> ARKLsMassTimesSetupFn *)
      type mass_times_setup_fn = float -> unit

      (** Callback functions that compute the mass matrix times a vector. In
          the call [mass_times_vec_fn t v mv],
          - [t] is the independent variable,
          - [v] is the vector to multiply, and
          - [mv] is the computed output
                 vector—{% $\mathtt{mv} = M\mathtt{v}$%}.

          Raising {!Sundials.RecoverableFailure} indicates a recoverable
          error. Any other exception is treated as an unrecoverable error.

          {warning Neither the elements of [v] nor [mv] should be
                   accessed after the function has returned.}

          @noarkode <node> ARKLsMassTimesVecFn *)
      type 'd mass_times_vec_fn =
           float (* t *)
        -> 'd    (* v *)
        -> 'd    (* Mv *)
        -> unit

      (** Create an Arkode-specific mass linear solver from a generic
          iterative linear solver. The boolean argument indicates whether
          the mass matrix depends on the independent variable [t], if not
          it is only computed and factored once.

          NB: a [mass_times_setup_fn] is not supported in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

          NB: The boolean argument is ignored in
          {{!Sundials_Config.sundials_version}Config.sundials_version} < 3.0.0.

          @noarkode <node> ARKStepSetMassLinearSolver
          @noarkode <node> ARKStepSetMassTimes *)
      val solver :
        ('m, 'd, 'k, [>`Iter]) LinearSolver.t
        -> ?mass_times_setup:mass_times_setup_fn
        -> 'd mass_times_vec_fn
        -> bool
        -> ('d, 'k) preconditioner
        -> ('d, 'k) solver

      (** {3:arkspilsmassset Solver parameters} *)

      (** Sets the factor by which the Krylov linear solver's convergence
          test constant is reduced from the Newton iteration test constant.
          This factor must be >= 0; passing 0 specifies the default (0.05).

          @noarkode <node> ARKStepSetMassEpsLin *)
      val set_eps_lin : ('d, 'k) session -> float -> unit

      (** Sets the factor for converting from the integrator tolerance (WRMS
          norm) to the linear solver tolerance (L2 norm). That is,
          {% $\mathit{tol}_{\mathsf{L2}} =
              \mathit{fact}\cdot\mathit{tol}_{\mathsf{WRMS}}$ %}.
          The given value is used directly if it is greater than zero.
          If it is zero (the default), then the square root of the state
          vector length is used.
          If it is less than zero, then the square root of the dot product of a
          state vector full of ones with itself is used.

          @since 5.4.0
          @noarkode <node> ARKStepSetMassLSNormFactor *)
      val set_ls_norm_factor : ('d, 'k) session -> float -> unit

      (** {3:arkspilsmassstats Solver statistics} *)

      (** Returns the sizes of the real and integer workspaces used by the
          spils linear solver.

          @noarkode <node> ARKStepGetMassWorkSpace
          @return ([real_size], [integer_size]) *)
      val get_work_space       : ('d, 'k) session -> int * int

      (** Returns the cumulative number of linear iterations.

          @noarkode <node> ARKStepGetNumMassIters *)
      val get_num_lin_iters    : ('d, 'k) session -> int

      (** Returns the cumulative number of linear convergence failures.

          @noarkode <node> ARKStepGetNumMassConvFails *)
      val get_num_conv_fails   : ('d, 'k) session -> int

      (** Returns the cumulative number of calls to the mass-matrix-vector
          setup function.

          @since 3.0.0
          @noarkode <node> ARKStepGetNumMTSetups *)
      val get_num_mtsetups : ('d, 'k) session -> int

      (** Returns the cumulative number of calls to the mass-matrix-vector
          product function ({!mass_times_vec_fn}).

          @noarkode <node> ARKStepGetNumMassMult
          @since Sundials 2.6.3 *)
      val get_num_mass_mult : ('d, 'k) session -> int

      (** Returns the cumulative number of calls to the setup function with
          [jok=false].

          @noarkode <node> ARKStepGetNumMassPrecEvals *)
      val get_num_prec_evals   : ('d, 'k) session -> int

      (** Returns the cumulative number of calls to the preconditioner solve
          function.

          @noarkode <node> ARKStepGetNumMassPrecSolves *)
      val get_num_prec_solves  : ('d, 'k) session -> int

      (** {3:arkspilsmasslowlevel Low-level solver manipulation}

          The {!init} and {!reinit} functions are the preferred way to set or
          change preconditioner functions. These low-level functions are
          provided for experts who want to avoid resetting internal counters
          and other associated side-effects. *)

      (** Change the preconditioner functions.

          @noarkode <node> ARKStepSetMassPreconditioner
          @noarkode <node> ARKLsMassPrecSolveFn
          @noarkode <node> ARKLsMassPrecSetupFn *)
      val set_preconditioner :
        ('d, 'k) session
        -> ?setup:'d prec_setup_fn
        -> 'd prec_solve_fn
        -> unit

      (** Change the mass matrix-times-vector function.

          @noarkode <node> ARKStepSetMassTimes
          @noarkode <node> ARKLsMassTimesVecFn *)
      val set_times :
        ('d, 'k) session
        -> ?mass_times_setup:mass_times_setup_fn
        -> 'd mass_times_vec_fn
        -> unit

    end (* }}} *)
  end (* }}} *)

  (** {2:tols Tolerances} *)

  (** Functions that compute the weighted RMS residual weights. The call
      [rfun y rwt] takes the dependent variable vector [y] and fills the
      residual-weight vector [rwt] with positive values or raises
      {!Sundials.NonPositiveEwt}. Other exceptions are eventually propagated,
      but should be avoided ([ffun] is not allowed to abort the solver). *)
  type 'data res_weight_fun = 'data -> 'data -> unit

  type ('data, 'kind) res_tolerance =
    | ResStolerance of float
      (** [abs] : scalar absolute residual tolerance. *)
    | ResVtolerance of ('data, 'kind) Nvector.t
      (** [abs] : vector of absolute residual tolerances. *)
    | ResFtolerance of 'data res_weight_fun
      (** Compute the residual weight vector. *)

  (** {2:solver Solver initialization and use} *)

  (** The linearity of the implicit portion of the problem. *)
  type linearity =
    | Linear of bool  (** Implicit portion is linear. Specifies whether
                          $f_I(t, y)$ or the preconditioner is
                          time-dependent. *)
    | Nonlinear       (** Implicit portion is nonlinear. *)

  (** The form of the initial value problem. *)
  type ('d, 'k) problem

  (** Diagonally Implicit Runge-Kutta (DIRK) solution of stiff problem.
      The fields are as follows.
      - [nlsolver], the nonlinear solver used for implicit stage solves,
      - [lsolver], used by [nlsolver]s based on Newton interation.
      - [linearity], specifies whether the implicit portion is linear or
                     nonlinear (the default), and
      - [irhsfn], the implicit portion of the right-hand-side function,

      If an [nlsolver] is not specified, then the
      {{!Sundials_NonlinearSolver.Newton}Newton} module is used by default.
      Specifying an [nlsolver] that requires a linear solver without
      specifying an [lsolver] results in a {!NonlinearInitFailure} (or
      {!IllInput} for Sundials < 4.0.0) exception on the first call to
      {!solve_normal} or {!solve_one_step}. *)
  val implicit :
    ?nlsolver : ('data, 'kind,
                  (('data, 'kind) session) Sundials_NonlinearSolver.integrator)
                Sundials_NonlinearSolver.t
    -> ?lsolver : ('data, 'kind) linear_solver
    -> ?linearity : linearity
    -> 'data rhsfn
    -> ('data, 'kind) problem

  (** Explicit Runge-Kutta (ERK) solution of non-stiff problem. The argument
      specifies the explicit portion of the right-hand-side problem. *)
  val explicit : 'data rhsfn -> ('data, 'kind) problem

  (** Additive Runge-Kutta (ARK) solution of multi-rate problem. The arguments
      are as described under {!implicit} and {!explicit}. *)
  val imex :
    ?nlsolver : ('data, 'kind,
                  (('data, 'kind) session) Sundials_NonlinearSolver.integrator)
                Sundials_NonlinearSolver.t
    -> ?lsolver  : ('data, 'kind) linear_solver
    -> ?linearity : linearity
    -> fi:'data rhsfn
    -> 'data rhsfn
    -> ('data, 'kind) problem

  (** Creates and initializes a session with the solver. The call
      {[init problem tol ~restol ~order ~mass:msolver ~roots:(nroots, g) t0 y0]}
      has as arguments:
      - [problem], specifies the problem to solve (see {!problem}),
      - [tol],     the integration tolerances,
      - [restol],  (optional) mass matrix residual tolerances,
      - [order],   the order of accuracy for the integration method,
      - [msolver], a linear mass matrix solver is required only if the problem
                   involves a non-identity mass matrix,
      - [nroots],  the number of root functions,
      - [g],       the root function ([(nroots, g)] defaults to {!no_roots}),
      - [t0],     the initial value of the independent variable, and
      - [y0],     a vector of initial values that also determines the number
                  of equations.

      The allowed values for [order] are:
      - for explicit methods: {% $2 \le \mathtt{order} \le 6$%},
      - for implicit methods: {% $2 \le \mathtt{order} \le 5$%}, and
      - for imex methods: {% $3 \le \mathtt{order} \le 5$%}.

      This function does everything necessary to initialize a session, i.e.,
      it makes the calls referenced below. The {!solve_normal} and
      {!solve_one_step} functions may be called directly.

      @noarkode <node> ARKStepCreate
      @noarkode <node> ARKStepSetLinear
      @noarkode <node> ARKStepSetNonlinear
      @noarkode <node> ARKStepSetLinearSolver
      @noarkode <node> ARKStepSetMassLinearSolver
      @noarkode <node> ARKStepSetNonlinearSolver
      @noarkode <node> ARKStepRootInit
      @noarkode <node> ARKStepSStolerances
      @noarkode <node> ARKStepSVtolerances
      @noarkode <node> ARKStepWFtolerances
      @noarkode <node> ARKStepResStolerance
      @noarkode <node> ARKStepResVtolerance
      @noarkode <node> ARKStepResFtolerance
      @noarkode <node> ARKStepSetOrder *)
  val init :
      ('data, 'kind) problem
      -> ('data, 'kind) tolerance
      -> ?restol:(('data, 'kind) res_tolerance)
      -> ?order:int
      -> ?mass:('data, 'kind) Mass.solver
      -> ?roots:(int * 'data rootsfn)
      -> float
      -> ('data, 'kind) Nvector.t
      -> ('data, 'kind) session

  (** Integrates an ODE system over an interval. The call
      [tret, r = solve_normal s tout yout] has as arguments
      - [s], a solver session,
      - [tout], the next time at which a solution is desired, and,
      - [yout], a vector to store the computed solution.

      It returns [tret], the time reached by the solver, which will be equal to
      [tout] if no errors occur, and, [r], a {!solver_result}.

      @noarkode <node> ARKStepEvolve (ARK_NORMAL)
      @raise IllInput Missing or illegal solver inputs.
      @raise TooClose The initial and final times are too close to each other and not initial step size was specified.

      @raise TooMuchWork The requested time could not be reached in [mxstep] internal steps.
      @raise TooMuchAccuracy The requested accuracy could not be satisfied.
      @raise ErrFailure Too many error test failures within a step or at the minimum step size.
      @raise ConvergenceFailure Too many convergence test failures within a step or at the minimum step size.
      @raise LinearInitFailure Linear solver initialization failed.
      @raise LinearSetupFailure Linear solver setup failed unrecoverably.
      @raise LinearSolveFailure Linear solver solution failed unrecoverably.
      @raise MassInitFailure Mass matrix solver initialization failed.
      @raise MassSetupFailure Mass matrix solver setup failed unrecoverably.
      @raise MassSolveFailure Mass matrix solver solution failed unrecoverably.
      @raise RhsFuncFailure Unrecoverable failure in one of the RHS functions.
      @raise FirstRhsFuncFailure Initial unrecoverable failure in one of the RHS functions.
      @raise RepeatedRhsFuncFailure Too many convergence test failures, or unable to estimate the initial step size, due to repeated recoverable errors in one of the right-hand side functions.
      @raise UnrecoverableRhsFuncFailure One of the right-hand side functions had a recoverable error, but no recovery was possible. This error can only occur after an error test failure at order one.
      @raise RootFuncFailure Failure in the rootfinding function [g].
      @raise PostprocStepFailure Failure in the Postprocess Step function. *)
  val solve_normal : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                          -> float * solver_result

  (** Like {!solve_normal} but returns after one internal solver step.

      @noarkode <node> ARKStepEvolve (ARK_ONE_STEP) *)
  val solve_one_step : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                          -> float * solver_result

  (** Returns the interpolated solution or derivatives.
      [get_dky s dky t k] computes the [k]th derivative of the function
      at time [t], i.e.,
      {% $\frac{d^\mathtt{k}}{\mathit{dt}^\mathtt{k}}y(\mathtt{t})$%},
      and stores it in [dky]. The arguments must satisfy
      {% $t_n - h_n \leq \mathtt{t} \leq t_n$%}—where $t_n$
      denotes {!get_current_time} and $h_n$ denotes {!get_last_step},—
      and {% $0 \leq \mathtt{k} \leq 3$%}.

      This function may only be called after a successful return from either
      {!solve_normal} or {!solve_one_step}.

      @noarkode <node> ARKStepGetDky
      @raise BadT [t] is not in the interval {% $[t_n - h_n, t_n]$%}.
      @raise BadK [k] is not in the range \{0, 1, ..., dord\}. *)
  val get_dky : ('d, 'k) session -> ('d, 'k) Nvector.t -> float -> int -> unit

  (** Reinitializes the solver with new parameters and state values. The
      values of the independent variable, i.e., the simulation time, and the
      state variables must be given. If given, [problem] specifies new
      callback functions and solvers, [order] changes the order of accuracy,
      [roots] specifies a new root finding function, and [mass] specifies a
      new linear mass matrix solver. The new problem must have the same size
      as the previous one.

      The allowed values for [order] are as for {!init}.

      @noarkode <node> ARKStepReInit
      @noarkode <node> ARKStepRootInit
      @noarkode <node> ARKStepSetOrder
      @noarkode <node> ARKStepSetLinearSolver
      @noarkode <node> ARKStepSetMassLinearSolver
      @noarkode <node> ARKStepSetNonlinearSolver *)
  val reinit :
    ('d, 'k) session
    -> ?problem:('d, 'k) problem
    -> ?order:int
    -> ?mass:('d, 'k) Mass.solver
    -> ?roots:(int * 'd rootsfn)
    -> float
    -> ('d, 'k) Nvector.t
    -> unit

  (** Resets the state to the given independent variable value and dependent
      variable vector. All previously set options, internal counter values,
      and step-size/error histories are retained.

      @since 5.4.0
      @noarkode <node> ARKStepReset *)
  val reset : ('d, 'k) session -> float -> ('d, 'k) Nvector.t

  (** Change the number of equations and unknowns between integrator steps.
      The call
      [resize s ~resize_nvec:rfn ~lsolver ~mass tol ~restol hscale ynew t0]
      has as arguments:
      - [s], the solver session to resize,
      - [rfn], a resize function that transforms nvectors in place-otherwise
               they are simply destroyed and recloned,
      - [lsolver], specify a new linear solver for the changed size,
      - [mass], specify a new mass matrix solver for the changed size,
      - [tol], tolerance values (ensures that any tolerance vectors are resized),
      - [restol], (optional) mass matrix residual tolerances,
      - [hscale], the next step will be of size {% $h \mathtt{hscale}$%},
      - [ynew], the newly-sized solution vector with the value {% $y(t_0)$%}, and
      - [t0], the current value of the independent variable $t_0$.

      If a new linear solver is not specified, any existing linear solver
      in use is destroyed and reinitialized; settings must be reconfigured
      after resizing. The same holds for the mass matrix solver.

      If the mass matrix residual tolerance was previously set to
      {{!res_tolerance}ResVtolerance} and [restol] is not given, then it is
      reset to {{!res_tolerance}ResStolerance} with the default value.

      The [tol] argument is ignored in versions 2.6.1 and 2.6.2 since it may
      cause a segmentation error.

      A call to this function disables inequality constraint checking. It can
      be reenabled by calling {!set_constraints}.

      @noarkode <node> ARKStepResize
      @noarkode <node> ARKStepSetLinearSolver
      @noarkode <node> ARKStepSetMassLinearSolver *)
  val resize :
    ('d, 'k) session
    -> ?resize_nvec:('d resize_fn)
    -> ?lsolver:('d, 'k) linear_solver
    -> ?mass:('d, 'k) Mass.solver
    -> ('d, 'k) tolerance
    -> ?restol:(('d, 'k) res_tolerance)
    -> float
    -> ('d, 'k) Nvector.t
    -> float
    -> unit

  (** {2:set Modifying the solver (optional input functions)} *)

  (** {3:arkset Optional inputs for ARKStep} *)

  (** Sets the integration tolerances.

      @noarkode <node> ARKStepSStolerances
      @noarkode <node> ARKStepSVtolerances
      @noarkode <node> ARKStepWFtolerances *)
  val set_tolerances : ('d, 'k) session -> ('d, 'k) tolerance -> unit

  (** Sets the residual tolerance.

      @noarkode <node> ARKStepResStolerance
      @noarkode <node> ARKStepResVtolerance
      @noarkode <node> ARKStepResFtolerance *)
  val set_res_tolerance : ('d, 'k) session -> ('d, 'k) res_tolerance -> unit

  (** Resets all optional input parameters to their default values. Neither
      the problem-defining functions nor the root-finding functions are
      changed.

      @noarkode <node> ARKStepSetDefaults *)
  val set_defaults : ('d, 'k) session -> unit

  (** Write step adaptivity and solver diagnostics to the given file.

      @noarkode <node> ARKStepSetDiagnostics *)
  val set_diagnostics : ('d, 'k) session -> Logfile.t -> unit

  (** Specifies the interpolation module used for output value interpolation
      and implicit method predictors.

      @since 5.2.0
      @noarkode <node> ARKStepSetInterpolantType *)
  val set_interpolant_type : ('d, 'k) session -> interpolant_type -> unit

  (** Specifies the degree of the polynomial interpolant used for output
      values and implicit method predictors.

      @since 5.2.0
      @noarkode <node> ARKStepSetInterpolantDegree *)
  val set_interpolant_degree : ('d, 'k) session -> int -> unit

  (** Do not write step adaptivity or solver diagnostics of a file.

      @noarkode <node> ARKStepSetDiagnostics *)
  val clear_diagnostics : ('d, 'k) session -> unit

  (** Configure the default error handler to write messages to a file.
      By default it writes to Logfile.stderr.

      @noarkode <node> ARKStepSetErrFile *)
  val set_error_file : ('d, 'k) session -> Logfile.t -> unit

  (** Specifies a custom function for handling error messages.
      The handler must not fail: any exceptions are trapped and discarded.

      @noarkode <node> ARKStepSetErrHandlerFn
      @noarkode <node> ARKErrHandlerFn *)
  val set_err_handler_fn
    : ('d, 'k) session -> (Util.error_details -> unit) -> unit

  (** Restores the default error handling function.

      @noarkode <node> ARKStepSetErrHandlerFn *)
  val clear_err_handler_fn : ('d, 'k) session -> unit

  (** Specifies the initial step size.

      @noarkode <node> ARKStepSetInitStep *)
  val set_init_step : ('d, 'k) session -> float -> unit

  (** Disables time step adaptivity and fix the step size for all internal
      steps. Pass [None] to restore time step adaptivity. This function is not
      recommended since there is no assurance of the validity of the computed
      solutions. It is provided primarily for code-to-code verification
      testing. Use in conjunction with {!set_min_step} or {!set_max_step}.

      @noarkode <node> ARKStepSetFixedStep *)
  val set_fixed_step : ('d, 'k) session -> float option -> unit

  (** Specifies the maximum number of constraint failures in a step before an
      error is signalled.

      @since 5.0.0
      @noarkode <node> ARKStepSetMaxNumConstrFails *)
  val set_max_num_constr_fails : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of messages warning that [t + h = t] on
      the next internal step.

      @noarkode <node> ARKStepSetMaxHnilWarns *)
  val set_max_hnil_warns : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of steps taken in attempting to reach
      a given output time.

      @noarkode <node> ARKStepSetMaxNumSteps *)
  val set_max_num_steps : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of error test failures permitted in
      attempting one step.

      @noarkode <node> ARKStepSetMaxErrTestFails *)
  val set_max_err_test_fails : ('d, 'k) session -> int -> unit

  (** Specifies a lower bound on the magnitude of the step size.

      @noarkode <node> ARKStepSetMinStep *)
  val set_min_step : ('d, 'k) session -> float -> unit

  (** Specifies an upper bound on the magnitude of the step size.

      @noarkode <node> ARKStepSetMaxStep *)
  val set_max_step : ('d, 'k) session -> float -> unit

  (** Sets all adaptivity and solver parameters to ‘best guess’ values. This
      routine takes into account the integration method (ERK, DIRK, or ARK)
      and a given method order; it should only be called after these have
      been set.

      @noarkode <node> ARKStepSetOptimalParams *)
  val set_optimal_params : ('d, 'k) session -> unit

  (** Limits the value of the independent variable [t] when solving.
      By default no stop time is imposed.

      @noarkode <node> ARKStepSetStopTime *)
  val set_stop_time : ('d, 'k) session -> float -> unit

  (** {3:arksetivp Optional inputs for IVP method selection} *)

  (** Enables both the implicit and explicit portions of a problem.

      @raise IllInput If $f_I$ and $f_E$ are not already specified.
      @noarkode <node> ARKStepSetImEx *)
  val set_imex : ('d, 'k) session -> unit

  (** Disables the implicit portion of a problem.

      @raise IllInput If $f_E$ is not already specified.
      @noarkode <node> ARKStepSetExplicit *)
  val set_explicit : ('d, 'k) session -> unit

  (** Disables the explicit portion of a problem.

      @raise IllInput If $f_I$ is not already specified.
      @noarkode <node> ARKStepSetImplicit *)
  val set_implicit : ('d, 'k) session -> unit

  (** Specifies a customized Butcher table or pair for the ERK, DIRK, or ARK
      method.

      To set an explicit table, give the [explicit_table] argument only.
      This also has the effect of calling {!set_explicit}. It may
      be better, though, to use the {!ERKStep} time-stepper in this case.

      To set an implicit table, give the [implicit_table] argument only.
      This also has the effect of calling {!set_implicit}.

      Otherwise, provide all arguments. This also has the effect of calling
      {!set_imex}.

      If neither [implicit_table] or [explicit_table] is given, the solver will
      run in fixed step mode and the step size must be set, or have been set,
      using either {!set_fixed_step} or {!set_init_step}. This feature is not
      available for Sundials versions prior to 2.7.0
      (the
      {{!Sundials_Config.NotImplementedBySundialsVersion}Config.NotImplementedBySundialsVersion}
      exception is raised).

      In versions of Sundials prior to 4.0.0, if both [implicit_table] and
      [explicit_table] are given, they must have the same values for
      [global_order] and [global_embeddded_order].

      In versions of Sundials prior to 2.7.0, the [coefficients] and [bembed]
      fields of the [explicit_table] are ignored when both [implicit_table]
      and [explicit_table] are given.

      @raise IllInput If $f_I$ and $f_E$ are not already specified.
      @noarkode <node> ARKStepSetTables *)
  val set_tables
    : ('d, 'k) session
      -> ?global_method_order:int
      -> ?global_embedding_order:int
      -> ?implicit_table:ButcherTable.t
      -> ?explicit_table:ButcherTable.t
      -> unit
      -> unit

  (** Use specific built-in Butcher tables for an ImEx system.

      @raise IllInput If $f_I$ and $f_E$ are not already specified.
      @noarkode <node> ARKStepSetTables
      @noarkode <node> ARKStepSetImEx *)
  val set_ark_table_num : ('d, 'k) session -> ButcherTable.ark_table -> unit

  (** Use specific built-in Butcher tables for an explicit integration of the
      problem.

      The {{!ButcherTable.erk_table}Fehlberg_13_7_8} method is not available
      prior to Sundials 2.7.0.
      The {{!ButcherTable.erk_table}Knoth_Wolke_3_3} method is not available
      prior to Sundials 4.0.0.

      @raise IllInput If $f_E$ is not already specified.
      @noarkode <node> ARKStepSetTables
      @noarkode <node> ARKStepSetExplicit *)
  val set_erk_table_num : ('d, 'k) session -> ButcherTable.erk_table -> unit

  (** Use specific built-in Butcher tables for an implicit integration of the
      problem.

      @raise IllInput If $f_I$ is not already specified.
      @noarkode <node> ARKStepSetTables
      @noarkode <node> ARKStepSetImplicit *)
  val set_dirk_table_num : ('d, 'k) session -> ButcherTable.dirk_table -> unit

  (** {3:arkadapt Optional inputs for time step adaptivity} *)

  (** Specifies the method and associated parameters used for time step
      adaptivity.

      @noarkode <node> ARKStepSetAdaptivityMethod
      @noarkode <node> ARKStepSetAdaptivityFn *)
  val set_adaptivity_method : ('d, 'k) session -> 'd adaptivity_method -> unit

  (** Specifies the fraction of the estimated explicitly stable step to use.
      Any non-positive argument resets to the default value (0.5).

      @noarkode <node> ARKStepSetCFLFraction *)
  val set_cfl_fraction : ('d, 'k) session -> float -> unit

  (** Specifies the bias to apply to the error estimates within accuracy-based
      adaptivity strategies.

      @noarkode <node> ARKStepSetErrorBias *)
  val set_error_bias : ('d, 'k) session -> float -> unit

  (** Specifies the step growth interval in which the step size will remain
      unchanged. In the call [set_fixed_step_bounds s lb ub], [lb] specifies a
      lower bound on the window to leave the step size fixed and [ub]
      specifies an upper bound. Any interval not containing 1.0 resets to
      default values.

      @noarkode <node> ARKStepSetFixedStepBounds *)
  val set_fixed_step_bounds : ('d, 'k) session -> float -> float -> unit

  (** Specifies the maximum step size growth factor upon a convergence
      failure on a stage solve within a step. Any value outside the
      interval {% $(0, 1]$%} resets to the default value.

      @noarkode <node> ARKStepSetMaxCFailGrowth *)
  val set_max_cfail_growth : ('d, 'k) session -> float -> unit

  (** Specifies the maximum step size growth factor upon multiple successive
      accuracy-based error failures in the solver. Any value outside the
      interval {% $(0, 1]$%} resets to the default value.

      @noarkode <node> ARKStepSetMaxEFailGrowth *)
  val set_max_efail_growth : ('d, 'k) session -> float -> unit

  (** Specifies the maximum allowed step size change following the very first
      integration step. Any value $\le 1$ resets to the default value.

      @noarkode <node> ARKStepSetMaxFirstGrowth *)
  val set_max_first_growth : ('d, 'k) session -> float -> unit

  (** Specifies the maximum growth of the step size between consecutive time
      steps. Any value $\le 1$ resets to the default value.

      @noarkode <node> ARKStepSetMaxGrowth *)
  val set_max_growth : ('d, 'k) session -> float -> unit

  (** Specifies the minimum allowed reduction factor in step size between
      step attempts that result from a temporal error failure in the
      integration process.

      @since 5.3.0
      @noarkode <node> ARKStepSetMinReduction *)
  val set_min_reduction : ('d, 'k) session -> float -> unit

  (** Specifies the safety factor to be applied to the accuracy-based
      estimated step. Any non-positive value resets to the default value.

      @noarkode <node> ARKStepSetSafetyFactor *)
  val set_safety_factor : ('d, 'k) session -> float -> unit

  (** Specifies the threshold for “multiple” successive error failures
      before the factor from {!set_max_efail_growth} is applied. Any
      non-positive value resets to the default value.

      @noarkode <node> ARKStepSetSmallNumEFails *)
  val set_small_num_efails : ('d, 'k) session -> float -> unit

  (** Sets a problem-dependent function to estimate a stable time step size
      for the explicit portion of the ODE system.

      @noarkode <node> ARKStepSetStabilityFn *)
  val set_stability_fn : ('d, 'k) session -> 'd stability_fn -> unit

  (** Clears the problem-dependent function that estimates a stable time step
      size for the explicit portion of the ODE system.

      @noarkode <node> ARKStepSetStabilityFn *)
  val clear_stability_fn : ('d, 'k) session -> unit

  (** {3:arksetimplicit Optional inputs for implicit stage solves} *)

  (** Solve the implicit portion of the problem using the accelerated
      fixed-point solver instead of the modified Newton iteration. The integer
      argument gives the maximum dimension of the accelerated subspace (i.e.,
      the number of vectors to store within the Anderson acceleration
      subspace).

      @noarkode <node> ARKStepSetFixedPoint *)
  val set_fixed_point : ('d, 'k) session -> int -> unit

  (** Solve the implicit portion of the problem using the modified Newton
      solver. Overrides a previous call to {!set_fixed_point}.

      @noarkode <node> ARKStepSetNewton *)
  val set_newton : ('d, 'k) session -> unit

  (** Specifies that the implicit portion of the problem is linear. The flag
      indicates whether the Jacobian of {% $f_I(t,y)$%}, or, when using an
      iterative linear solver, the preconditioner is time-dependent ([true])
      or not ([false]).

      @noarkode <node> ARKStepSetLinear *)
  val set_linear : ('d, 'k) session -> bool -> unit

  (** Specifies that the implicit portion of the problem is nonlinear.

      @noarkode <node> ARKStepSetNonlinear *)
  val set_nonlinear : ('d, 'k) session -> unit

  (** Method choices for predicting implicit solutions.

      @noarkode <node> Implicit predictors *)
  type predictor_method =
    | TrivialPredictor        (** Piece-wise constant interpolant
                                  {% $p_0(\tau) = y_{n-1}$%} %*)
    | MaximumOrderPredictor   (** An interpolant {% $p_q(t)$%} of polynomial
                                  order up to {% $q=3$%}. *)
    | VariableOrderPredictor  (** Decrease the polynomial degree for later RK
                                  stages. *)
    | CutoffOrderPredictor    (** Maximum order for early RK stages and
                                  first-order for later ones. *)
    | BootstrapPredictor      (** Second-order predictor based only on the
                                  current step. *)
    | MinimumCorrectionPredictor
                              (** Uses all preceding stage information within
                                  the current step for prediction. *)

  (** Specifies the method for predicting implicit solutions.

      @noarkode <node> ARKStepSetPredictorMethod *)
  val set_predictor_method : ('d, 'k) session -> predictor_method -> unit

  (** Specifies the maximum number of nonlinear solver iterations permitted
      per RK stage at each step.

      @raise NonlinearOperationError Nonlinear solver not configured
      @noarkode <node> ARKStepSetMaxNonlinIters *)
  val set_max_nonlin_iters : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of nonlinear solver convergence failures
      permitted during one step.

      @noarkode <node> ARKStepSetMaxConvFails *)
  val set_max_conv_fails : ('d, 'k) session -> int -> unit

  (** A function to be called {e after} the predictor algorithm to update the
      predictor.

      The first argument is the current value of the independent variable. The
      second is the predicated stage solution which may be updated directly.

      {warning The array and its contents should not be accessed after the
               function has returned}.

      All exceptions are treated as unrecoverable errors.

      @noarkode <node> ARKStepStagePredictFn *)
  type 'd stage_predict_fn = float -> 'd -> unit

  (** Set the function called after the predictor algorithm and before the
      calculation of an implicit stage solution.

      If a stage prediction function is set and the current predictor method,
      see {!set_predictor_method}, is
      {{!predictor_method}MinimumCorrectionPredictor}, then the
      {{!predictor_method}TrivialPredictor} will be used instead.

      @since 5.0.0
      @noarkode <node> ARKStepSetStagePredictFn *)
  val set_stage_predict_fn : ('d, 'k) session -> 'd stage_predict_fn -> unit

  (** Clear the function called after the predictor algorithm.

      @since 5.0.0
      @noarkode <node> ARKStepSetStagePredictFn *)
  val clear_stage_predict_fn : ('d, 'k) session -> unit

  (** Specifies the safety factor used in the nonlinear convergence test.

      @noarkode <node> ARKStepSetNonlinConvCoef *)
  val set_nonlin_conv_coef : ('d, 'k) session -> float -> unit

  (** Specifies a vector defining inequality constraints for each
      component of the solution vector.  See {!Sundials.Constraint}.

      @since 5.0.0
      @noarkode <node> ARKStepSetConstraints *)
  val set_constraints : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

  (** Specifies the constant used in estimating the nonlinear solver
      convergence rate.

      @noarkode <node> ARKStepSetNonlinCRDown *)
  val set_nonlin_crdown : ('d, 'k) session -> float -> unit

  (** Specifies the nonlinear correction threshold beyond which the iteration
      will be declared divergent.

      @noarkode <node> ARKStepSetNonlinRDiv *)
  val set_nonlin_rdiv : ('d, 'k) session -> float -> unit

  (** Specifies a scaled step size ratio tolerance beyond which the linear
      solver setup routine will be signalled.

      @noarkode <node> ARKStepSetDeltaGammaMax *)
  val set_delta_gamma_max : ('d, 'k) session -> float -> unit

  (** Specifies the frequency of calls to the linear solver setup routine.
      Positive values specify the number of time steps between setup calls,
      negative values force recomputation at each Newton step, and zero values
      reset to the default (20).

      @noarkode <node> ARKStepSetLSetupFrequency *)
  val set_lsetup_frequency : ('d, 'k) session -> int -> unit

  (** Set a post processing step function.

      @since 2.7.0
      @raise Config.NotImplementedBySundialsVersion Post processing not available
      @noarkode <node> ARKSetPostprocessStepFn *)
  val set_postprocess_step_fn : ('d, 'k) session -> 'd postprocess_step_fn -> unit

  (** Clear the post processing step function.

      @since 2.7.0
      @raise Config.NotImplementedBySundialsVersion Post processing not available
      @noarkode <node> ARKSetPostprocessStepFn *)
  val clear_postprocess_step_fn : ('d, 'k) session -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the real and integer workspace sizes.

      @noarkode <node> ARKStepGetWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space          : ('d, 'k) session -> int * int

  (** Returns the cumulative number of internal steps taken by the solver.

      @noarkode <node> ARKStepGetNumSteps *)
  val get_num_steps           : ('d, 'k) session -> int

  (** Returns the cumulative number of stability-limited steps taken by the
      solver.

      @noarkode <node> ARKStepGetNumExpSteps *)
  val get_num_exp_steps       : ('d, 'k) session -> int

  (** Returns the cumulative number of accuracy-limited steps taken by the
      solver.

      @noarkode <node> ARKStepGetNumAccSteps *)
  val get_num_acc_steps       : ('d, 'k) session -> int

  (** Returns the cumulative number of steps attempted by the solver.

      @noarkode <node> ARKStepGetNumStepAttempts *)
  val get_num_step_attempts   : ('d, 'k) session -> int

  (** Returns the number of calls to the right-hand side functions.
      In the call [(nfe_evals, nfi_evals) = get_num_rhs_evals s],
      - [nfe_evals] is the number of calls to {% $f_E$%}, and
      - [nfi_evals] is the number of calls to {% $f_I$%}.

      @noarkode <node> ARKStepGetNumRhsEvals *)
  val get_num_rhs_evals       : ('d, 'k) session -> int * int

  (** Returns the number of local error test failures that have occurred.

      @noarkode <node> ARKStepGetNumErrTestFails *)
  val get_num_err_test_fails  : ('d, 'k) session -> int

  (** Returns the integration step size taken on the last successful internal
      step.

      @noarkode <node> ARKStepGetLastStep *)
  val get_last_step           : ('d, 'k) session -> float

  (** Returns the integration step size to be attempted on the next internal
      step.

      @noarkode <node> ARKStepGetCurrentStep *)
  val get_current_step        : ('d, 'k) session -> float

  (** Returns the the value of the integration step size used on the first
      step.

      @noarkode <node> ARKStepGetActualInitStep *)
  val get_actual_init_step    : ('d, 'k) session -> float

  (** Returns the the current internal time reached by the solver.

      @noarkode <node> ARKStepGetCurrentTime *)
  val get_current_time        : ('d, 'k) session -> float

  (** Returns the current state vector. This vector provides direct access to
      the data within the integrator.

      @since 5.0.0
      @noarkode <node> ARKStepGetCurrentState *)
  val get_current_state : ('d, 'k) session -> 'd

  (** Internal data required to construct the current nonlinear implicit
      system within a nonlinear solver. *)
  type 'd nonlin_system_data = {
    tcur  : float;
      (** Independent variable value for slow stage {% $t^S_{n,i}$ %}. *)
    zpred : 'd;
      (** Predicted nonlinear solution {% $z_{\mathit{pred}}$ %}. This
          data must not be changed. *)
    zi    : 'd;
      (** Stage vector {% $z_i$ %}. This data may not be current and may
          need to be filled. *)
    fi    : 'd;
      (** Memory available for evaluating the slow right-hand side
          {% $f^S(t^S_{n,i}, z_i)$ %}. This data may not be current and may
          need to be filled. *)
    gamma : float;
      (** Current {% $\gamma$ %} for slow-stage calculation. *)
    sdata : 'd;
      (** Accumulated data from previous solution and stages
          {% $\tilde{a}_i$ %}. This data must not be changed. *)
  }

  (** Gives direct access to the internal data required to construct the
      current nonlinear implicit system within a nonlinear solver. This
      function should be called inside the nonlinear system function.
      If the nonlinear solver uses the [lsetup] or [lsolve] functions, then
      the nonlinear solver system function must fill the [zi] and [fi]
      vectors with, respectively, the current state and corresponding
      evaluation of the right-hand-side function:
      {% $z_i = z_{\mathit{pred}} + z_{\mathit{cor}}$ %} and
      {% $F_i = f^S(t^S_{n,i}, z_i)$ %} where {% $z_{\mathit{cor}}$ %} is the
      first argument of the nonlinear solver system function. Within a custom
      linear solver, then the vectors [zi] and [fi] are only current after
      an evaluation of the nonlinear system function.

      @since 5.4.0
      @noarkode <node> ARKStepGetNonlinearSystemData *)
  val get_nonlin_system_data : ('d, 'k) session -> 'd nonlin_system_data

  (** Computes the current stage state vector using the stored prediction and
      the supplied correction from the nonlinear solver. The call
      [compute_state s zcor z] computes {% $z_i(t) = z_{\mathit{pred}}
      + z_{\mathit{cor}}$ %}.

      @since 5.4.0
      @noarkode <node> ARKStepComputeState *)
  val compute_state : ('d, 'k) session
                      -> ('d, 'k) Nvector.t
                      -> ('d, 'k) Nvector.t
                      -> unit

  (** Returns the current value of {% $\gamma$ %}.
      This scalar appears in the internal Newton equation, either
      {% $A = I - \gamma J$ %} or {% $A = M - \gamma J$ %}.

      @since 5.0.0
      @noarkode <node> ARKStepGetCurrentGamma *)
  val get_current_gamma : ('d, 'k) session -> float

  (** Returns the implicit and explicit Butcher tables in use by the solver.
      In the call [bi, be = get_current_butcher_tables s], [bi] is the
      implicit butcher table and [be] is the explicit one.

      For versions of Sundials prior to 2.7.0, the fields of [vbe] are identical
      to those of [vbi] except for [state_values].

      @noarkode <node> ARKStepGetCurrentButcherTables *)
  val get_current_butcher_tables
    : ('d, 'k) session -> ButcherTable.t option * ButcherTable.t option

  (** Returns a suggested factor by which the user's tolerances should be
      scaled when too much accuracy has been requested for some internal step.

      @noarkode <node> ARKStepGetTolScaleFactor *)
  val get_tol_scale_factor : ('d, 'k) session -> float

  (** Returns the solution error weights at the current time.

      @noarkode <node> ARKStepGetErrWeights *)
  val get_err_weights : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

  (** Returns the residual error weights at the current time.

      @since 3.0.0
      @noarkode <node> ARKStepGetResWeights *)
  val get_res_weights : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

  (** Returns the vector of estimated local errors.

      @noarkode <node> ARKStepGetEstLocalErrors *)
  val get_est_local_errors : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

  (** Summaries of time-stepper statistics. *)
  type timestepper_stats = {
      exp_steps : int;
        (** Cumulative number of stability-limited solver steps. *)
      acc_steps : int;
        (** Cumulative number of accuracy-limited solver steps. *)
      step_attempts : int;
        (** Cumulative number of steps attempted by the solver. *)
      num_nfe_evals : int;
        (** Number of calls to the explicit right-hand side function. *)
      num_nfi_evals : int;
        (** Number of calls to the implicit right-hand side function. *)
      num_lin_solv_setups : int;
        (** Number of setups calls to the linear solver. *)
      num_err_test_fails : int;
        (** Number of local error test failures. *)
    }

  (** Returns a grouped set of time-stepper statistics.

      @noarkode <node> ARKStepGetTimestepperStats *)
  val get_timestepper_stats    : ('d, 'k) session -> timestepper_stats

  (** Returns a grouped set of integrator statistics.

      @noarkode <node> ARKStepGetStepStats *)
  val get_step_stats           : ('d, 'k) session -> step_stats

  (** Prints time-stepper statistics on the given channel.

      @noarkode <node> ARKStepGetTimestepperStats *)
  val print_timestepper_stats  : ('d, 'k) session -> out_channel -> unit

  (** Prints integrator statistics on the given channel.

      @noarkode <node> ARKStepGetStepStats *)
  val print_step_stats  : ('d, 'k) session -> out_channel -> unit

  (** {3:arkgetimplicit Implicit solver optional output functions} *)

  (** Returns the number of calls made to the linear solver's setup function.

      @noarkode <node> ARKStepGetNumLinSolvSetups *)
  val get_num_lin_solv_setups : ('d, 'k) session -> int

  (** Returns the cumulative number of nonlinear (functional or Newton)
      iterations.

      @raise NonlinearOperationError Nonlinear solver not configured
      @noarkode <node> ARKStepGetNumNonlinSolvIters *)
  val get_num_nonlin_solv_iters : ('d, 'k) session -> int

  (** Returns the cumulative number of nonlinear convergence failures.

      @noarkode <node> ARKStepGetNumNonlinSolvConvFails *)
  val get_num_nonlin_solv_conv_fails : ('d, 'k) session -> int

  (** Returns both the numbers of nonlinear iterations performed [nniters] and
      nonlinear convergence failures [nncfails].

      @raise NonlinearOperationError Nonlinear solver not configured
      @noarkode <node> ARKStepGetNonlinSolvStats
      @return ([nniters], [nncfails]) *)
  val get_nonlin_solv_stats : ('d, 'k) session -> int * int

  (** {2:roots Additional root-finding functions} *)

  (** [set_root_direction s dir] specifies the direction of zero-crossings to
      be located and returned. [dir] may contain one entry for each root
      function.

      @noarkode <node> ARKStepSetRootDirection *)
  val set_root_direction : ('d, 'k) session -> RootDirs.d array -> unit

  (** Like {!set_root_direction} but specifies a single direction for all root
      functions.

      @noarkode <node> ARKStepSetRootDirection *)
  val set_all_root_directions : ('d, 'k) session -> RootDirs.d -> unit

  (** Disables issuing a warning if some root function appears to be
      identically zero at the beginning of the integration.

      @noarkode <node> ARKStepSetNoInactiveRootWarn *)
  val set_no_inactive_root_warn : ('d, 'k) session -> unit

  (** Returns the number of root functions. *)
  val get_num_roots : ('d, 'k) session -> int

  (** Fills an array showing which functions were found to have a root.

      @noarkode <node> ARKStepGetRootInfo *)
  val get_root_info : ('d, 'k) session -> Roots.t -> unit

  (** Returns the cumulative number of calls made to the user-supplied root
      function g.

      @noarkode <node> ARKStepGetNumGEvals *)
  val get_num_g_evals : ('d, 'k) session -> int

  (** Returns the cumulative number of test failures.

      @since 5.0.0
      @noarkode <node> ARKStepGetNumConstrFails *)
  val get_num_constr_fails : ('d, 'k) session -> int

  (** Outputs all the solver parameters to the given file.

      @since 4.1.0
      @noarkode <node> ARKStepWriteParameters *)
  val write_parameters : ('d, 'k) session -> Logfile.t -> unit

  (** Outputs the current butcher table to the given file.

      @since 4.1.0
      @noarkode <node> ARKStepWriteButcher *)
  val write_butcher : ('d, 'k) session -> Logfile.t -> unit

end (* }}} *)

(** ERKStep Time-Stepping Module for nonstiff initial value problems.

    This module solves problems of the form
    {% $\dot{y} = f(t, y)$%}, {% $y(t_0) = y_0$%}.

    Its interface is structured as follows.
    {ol
      {- {{:#solver}Solver initialization and use}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}}

    @noarkode <node> Using ERKStep for C and C++ Applications *)
module ERKStep : sig (* {{{ *)

  (** A session with the ERKStep time-stepping solver.

      An example session with ERKStep ({openfile arkode_erk_skel.ml}): {[
#include "../../examples/ocaml/skeletons/arkode_erk_skel.ml"
      ]}

      @noarkode <node> Skeleton of main program *)
  type ('d, 'k) session = ('d, 'k, Arkode_impl.erkstep) Arkode_impl.session

  (** {2:solver Solver initialization and use} *)

  (** Creates and initializes a session with the solver. The call
      {[init tol ~order f ~roots:(nroots, g) t0 y0]}
      has as arguments:
      - [tol],    the integration tolerances,
      - [order],  the order of accuracy for the integration method,
      - [f],      the ODE right-hand side function,
      - [nroots], the number of root functions,
      - [g],      the root function ([(nroots, g)] defaults to {!no_roots}),
      - [t0],     the initial value of the independent variable, and
      - [y0],     a vector of initial values that also determines the number
                  of equations.

      This function does everything necessary to initialize a session, i.e.,
      it makes the calls referenced below. The {!solve_normal} and
      {!solve_one_step} functions may be called directly.

      @since 4.0.0
      @noarkode <node> ERKStepCreate
      @noarkode <node> ERKStepSetOrder
      @noarkode <node> ERKStepRootInit
      @noarkode <node> ERKStepSStolerances
      @noarkode <node> ERKStepSVtolerances
      @noarkode <node> ERKStepWFtolerances *)
  val init :
      ('data, 'kind) tolerance
      -> ?order:int
      -> 'data rhsfn
      -> ?roots:(int * 'data rootsfn)
      -> float
      -> ('data, 'kind) Nvector.t
      -> ('data, 'kind) session

  (** Integrates an ODE system over an interval. The call
      [tret, r = solve_normal s tout yout] has as arguments
      - [s], a solver session,
      - [tout], the next time at which a solution is desired, and,
      - [yout], a vector to store the computed solution.

      It returns [tret], the time reached by the solver, which will be equal to
      [tout] if no errors occur, and, [r], a {!solver_result}.

      @noarkode <node> ERKStepEvolve (ARK_NORMAL)
      @raise IllInput Missing or illegal solver inputs.
      @raise TooMuchWork The requested time could not be reached in [mxstep] internal steps.
      @raise TooMuchAccuracy The requested accuracy could not be satisfied.
      @raise ErrFailure Too many error test failures within a step or at the minimum step size.
      @raise RhsFuncFailure Unrecoverable failure in one of the RHS functions.
      @raise FirstRhsFuncFailure Initial unrecoverable failure in one of the RHS functions.
      @raise RepeatedRhsFuncFailure Too many convergence test failures, or unable to estimate the initial step size, due to repeated recoverable errors in one of the right-hand side functions.
      @raise UnrecoverableRhsFuncFailure One of the right-hand side functions had a recoverable error, but no recovery was possible. This error can only occur after an error test failure at order one.
      @raise RootFuncFailure Failure in the rootfinding function [g]. *)
  val solve_normal : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                          -> float * solver_result

  (** Like {!solve_normal} but returns after one internal solver step.

      @noarkode <node> ERKStepEvolve (ARK_ONE_STEP) *)
  val solve_one_step : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                          -> float * solver_result

  (** Returns the interpolated solution or derivatives.
      [get_dky s dky t k] computes the [k]th derivative of the function
      at time [t], i.e.,
      {% $\frac{d^\mathtt{k}}{\mathit{dt}^\mathtt{k}}y(\mathtt{t})$%},
      and stores it in [dky]. The arguments must satisfy
      {% $t_n - h_n \leq \mathtt{t} \leq t_n$%}—where $t_n$
      denotes {!get_current_time} and $h_n$ denotes {!get_last_step},—
      and {% $0 \leq \mathtt{k} \leq 3$%}.

      This function may only be called after a successful return from either
      {!solve_normal} or {!solve_one_step}.

      @noarkode <node> ERKStepGetDky
      @raise BadT [t] is not in the interval {% $[t_n - h_n, t_n]$%}.
      @raise BadK [k] is not in the range \{0, 1, ..., dord\}. *)
  val get_dky : ('d, 'k) session -> ('d, 'k) Nvector.t -> float -> int -> unit

  (** Reinitializes the solver with new parameters and state values. The
      values of the independent variable, i.e., the simulation time, and the
      state variables must be given. If given, [order] changes the order of
      accuracy, and [roots] specifies a new root finding function. The new
      problem must have the same size as the previous one.

      The allowed values for [order] are as for {!init}.

      @noarkode <node> ERKStepReInit
      @noarkode <node> ERKStepRootInit
      @noarkode <node> ERKStepSetOrder *)
  val reinit :
    ('d, 'k) session
    -> ?order:int
    -> ?roots:(int * 'd rootsfn)
    -> float
    -> ('d, 'k) Nvector.t
    -> unit

  (** Resets the state to the given independent variable value and dependent
      variable vector. All previously set options, internal counter values,
      and step-size/error histories are retained.

      @since 5.4.0
      @noarkode <node> ERKStepReset *)
  val reset : ('d, 'k) session -> float -> ('d, 'k) Nvector.t

  (** Change the number of equations and unknowns between integrator steps.
      The call
      [resize s ~resize_nvec:rfn tol hscale ynew t0]
      has as arguments:
      - [s], the solver session to resize,
      - [rfn], a resize function that transforms nvectors in place-otherwise
               they are simply destroyed and recloned,
      - [tol], tolerance values (ensures that any tolerance vectors are resized),
      - [hscale], the next step will be of size {% $h \mathtt{hscale}$%},
      - [ynew], the newly-sized solution vector with the value {% $y(t_0)$%}, and
      - [t0], the current value of the independent variable $t_0$.

      A call to this function disables inequality constraint checking. It can
      be reenabled by calling {!set_constraints}.

      @noarkode <node> ERKStepResize *)
  val resize :
    ('d, 'k) session
    -> ?resize_nvec:('d resize_fn)
    -> ('d, 'k) tolerance
    -> float
    -> ('d, 'k) Nvector.t
    -> float
    -> unit

  (** {2:set Modifying the solver (optional input functions)} *)

  (** {3:erkset Optional inputs for ERKStep} *)

  (** Sets the integration tolerances.

      @noarkode <node> ERKStepSStolerances
      @noarkode <node> ERKStepSVtolerances
      @noarkode <node> ERKStepWFtolerances *)
  val set_tolerances : ('d, 'k) session -> ('d, 'k) tolerance -> unit

  (** Resets all optional input parameters to their default values. Neither
      the problem-defining functions nor the root-finding functions are
      changed.

      @noarkode <node> ERKStepSetDefaults *)
  val set_defaults : ('d, 'k) session -> unit

  (** Write step adaptivity and solver diagnostics to the given file.

      @noarkode <node> ERKStepSetDiagnostics *)
  val set_diagnostics : ('d, 'k) session -> Logfile.t -> unit

  (** Specifies the interpolation module used for output value interpolation
      and implicit method predictors.

      @since 5.2.0
      @noarkode <node> ERKStepSetInterpolantType *)
  val set_interpolant_type : ('d, 'k) session -> interpolant_type -> unit

  (** Specifies the degree of the polynomial interpolant used for output
      values and implicit method predictors.

      @since 5.2.0
      @noarkode <node> ERKStepSetInterpolantDegree *)
  val set_interpolant_degree : ('d, 'k) session -> int -> unit

  (** Do not write step adaptivity or solver diagnostics of a file.

      @noarkode <node> ERKStepSetDiagnostics *)
  val clear_diagnostics : ('d, 'k) session -> unit

  (** Configure the default error handler to write messages to a file.
      By default it writes to Logfile.stderr.

      @noarkode <node> ERKStepSetErrFile *)
  val set_error_file : ('d, 'k) session -> Logfile.t -> unit

  (** Specifies a custom function for handling error messages.
      The handler must not fail: any exceptions are trapped and discarded.

      @noarkode <node> ERKStepSetErrHandlerFn
      @noarkode <node> ARKErrHandlerFn *)
  val set_err_handler_fn
    : ('d, 'k) session -> (Util.error_details -> unit) -> unit

  (** Restores the default error handling function.

      @noarkode <node> ERKStepSetErrHandlerFn *)
  val clear_err_handler_fn : ('d, 'k) session -> unit

  (** Specifies the initial step size.

      @noarkode <node> ERKStepSetInitStep *)
  val set_init_step : ('d, 'k) session -> float -> unit

  (** Disables time step adaptivity and fix the step size for all internal
      steps. Pass [None] to restore time step adaptivity. This function is not
      recommended since there is no assurance of the validity of the computed
      solutions. It is provided primarily for code-to-code verification
      testing. Use in conjunction with {!set_min_step} or {!set_max_step}.

      @noarkode <node> ERKStepSetFixedStep *)
  val set_fixed_step : ('d, 'k) session -> float option -> unit

  (** Specifies the maximum number of constraint failures in a step before an
      error is signalled.

      @since 5.0.0
      @noarkode <node> ERKStepSetMaxNumConstrFails *)
  val set_max_num_constr_fails : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of messages warning that [t + h = t] on
      the next internal step.

      @noarkode <node> ERKStepSetMaxHnilWarns *)
  val set_max_hnil_warns : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of steps taken in attempting to reach
      a given output time.

      @noarkode <node> ERKStepSetMaxNumSteps *)
  val set_max_num_steps : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of error test failures permitted in
      attempting one step.

      @noarkode <node> ERKStepSetMaxErrTestFails *)
  val set_max_err_test_fails : ('d, 'k) session -> int -> unit

  (** Specifies a lower bound on the magnitude of the step size.

      @noarkode <node> ERKStepSetMinStep *)
  val set_min_step : ('d, 'k) session -> float -> unit

  (** Specifies an upper bound on the magnitude of the step size.

      @noarkode <node> ERKStepSetMaxStep *)
  val set_max_step : ('d, 'k) session -> float -> unit

  (** Limits the value of the independent variable [t] when solving.
      By default no stop time is imposed.

      @noarkode <node> ERKStepSetStopTime *)
  val set_stop_time : ('d, 'k) session -> float -> unit

  (** {3:erksetivp Optional inputs for IVP method selection} *)

  (** Specifies a customized Butcher table.

      @noarkode <node> ERKStepSetTable *)
  val set_table
    : ('d, 'k) session
      -> ButcherTable.t
      -> unit

  (** Use a specific built-in Butcher table for integration.

      @noarkode <node> ERKStepSetTableNum *)
  val set_table_num : ('d, 'k) session -> ButcherTable.erk_table -> unit

  (** {3:erksetadap Optional inputs for time step adaptivity} *)

  (** Specifies the method and associated parameters used for time step
      adaptivity.

      @noarkode <node> ERKStepSetAdaptivityMethod
      @noarkode <node> ERKStepSetAdaptivityFn *)
  val set_adaptivity_method : ('d, 'k) session -> 'd adaptivity_method -> unit

  (** Specifies the fraction of the estimated explicitly stable step to use.
      Any non-positive argument resets to the default value (0.5).

      @noarkode <node> ERKStepSetCFLFraction *)
  val set_cfl_fraction : ('d, 'k) session -> float -> unit

  (** Specifies the bias to apply to the error estimates within accuracy-based
      adaptivity strategies.

      @noarkode <node> ERKStepSetErrorBias *)
  val set_error_bias : ('d, 'k) session -> float -> unit

  (** Specifies the step growth interval in which the step size will remain
      unchanged. In the call [set_fixed_step_bounds s lb ub], [lb] specifies a
      lower bound on the window to leave the step size fixed and [ub]
      specifies an upper bound. Any interval not containing 1.0 resets to
      default values.

      @noarkode <node> ERKStepSetFixedStepBounds *)
  val set_fixed_step_bounds : ('d, 'k) session -> float -> float -> unit

  (** Specifies the maximum step size growth factor upon multiple successive
      accuracy-based error failures in the solver. Any value outside the
      interval {% $(0, 1]$%} resets to the default value.

      @noarkode <node> ERKStepSetMaxEFailGrowth *)
  val set_max_efail_growth : ('d, 'k) session -> float -> unit

  (** Specifies the maximum allowed step size change following the very first
      integration step. Any value $\le 1$ resets to the default value.

      @noarkode <node> ERKStepSetMaxFirstGrowth *)
  val set_max_first_growth : ('d, 'k) session -> float -> unit

  (** Specifies the maximum growth of the step size between consecutive time
      steps. Any value $\le 1$ resets to the default value.

      @noarkode <node> ERKStepSetMaxGrowth *)
  val set_max_growth : ('d, 'k) session -> float -> unit

  (** Specifies the minimum allowed reduction factor in step size between
      step attempts that result from a temporal error failure in the
      integration process.

      @since 5.3.0
      @noarkode <node> ERKStepSetMinReduction *)
  val set_min_reduction : ('d, 'k) session -> float -> unit

  (** Specifies the safety factor to be applied to the accuracy-based
      estimated step. Any non-positive value resets to the default value.

      @noarkode <node> ERKStepSetSafetyFactor *)
  val set_safety_factor : ('d, 'k) session -> float -> unit

  (** Specifies the threshold for “multiple” successive error failures
      before the factor from {!set_max_efail_growth} is applied. Any
      non-positive value resets to the default value.

      @noarkode <node> ERKStepSetSmallNumEFails *)
  val set_small_num_efails : ('d, 'k) session -> float -> unit

  (** Sets a problem-dependent function to estimate a stable time step size
      for the explicit portion of the ODE system.

      @noarkode <node> ERKStepSetStabilityFn *)
  val set_stability_fn : ('d, 'k) session -> 'd stability_fn -> unit

  (** Clears the problem-dependent function that estimates a stable time step
      size for the explicit portion of the ODE system.

      @noarkode <node> ERKStepSetStabilityFn *)
  val clear_stability_fn : ('d, 'k) session -> unit

  (** Set a post processing step function.

      @noarkode <node> ERKSetPostprocessStepFn *)
  val set_postprocess_step_fn : ('d, 'k) session -> 'd postprocess_step_fn -> unit

  (** Clear the post processing step function.

      @noarkode <node> ERKSetPostprocessStepFn *)
  val clear_postprocess_step_fn : ('d, 'k) session -> unit

  (** Specifies a vector defining inequality constraints for each
      component of the solution vector.  See {!Sundials.Constraint}.

      @since 5.0.0
      @noarkode <node> ERKStepSetConstraints *)
  val set_constraints : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the real and integer workspace sizes.

      @noarkode <node> ERKStepGetWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space          : ('d, 'k) session -> int * int

  (** Returns the cumulative number of internal steps taken by the solver.

      @noarkode <node> ERKStepGetNumSteps *)
  val get_num_steps           : ('d, 'k) session -> int

  (** Returns the cumulative number of stability-limited steps taken by the
      solver.

      @noarkode <node> ERKStepGetNumExpSteps *)
  val get_num_exp_steps       : ('d, 'k) session -> int

  (** Returns the cumulative number of accuracy-limited steps taken by the
      solver.

      @noarkode <node> ERKStepGetNumAccSteps *)
  val get_num_acc_steps       : ('d, 'k) session -> int

  (** Returns the cumulative number of steps attempted by the solver.

      @noarkode <node> ERKStepGetNumStepAttempts *)
  val get_num_step_attempts   : ('d, 'k) session -> int

  (** Returns the number of calls to the right-hand side function.

      @noarkode <node> ERKStepGetNumRhsEvals *)
  val get_num_rhs_evals       : ('d, 'k) session -> int

  (** Returns the number of local error test failures that have occurred.

      @noarkode <node> ERKStepGetNumErrTestFails *)
  val get_num_err_test_fails  : ('d, 'k) session -> int

  (** Returns the integration step size taken on the last successful internal
      step.

      @noarkode <node> ERKStepGetLastStep *)
  val get_last_step           : ('d, 'k) session -> float

  (** Returns the integration step size to be attempted on the next internal
      step.

      @noarkode <node> ERKStepGetCurrentStep *)
  val get_current_step        : ('d, 'k) session -> float

  (** Returns the the value of the integration step size used on the first
      step.

      @noarkode <node> ERKStepGetActualInitStep *)
  val get_actual_init_step    : ('d, 'k) session -> float

  (** Returns the the current internal time reached by the solver.

      @noarkode <node> ERKStepGetCurrentTime *)
  val get_current_time        : ('d, 'k) session -> float

  (** Returns the Butcher table in use by the solver.

      @noarkode <node> ERKStepGetCurrentButcherTable *)
  val get_current_butcher_table : ('d, 'k) session -> ButcherTable.t

  (** Returns a suggested factor by which the user's tolerances should be
      scaled when too much accuracy has been requested for some internal step.

      @noarkode <node> ERKStepGetTolScaleFactor *)
  val get_tol_scale_factor : ('d, 'k) session -> float

  (** Returns the solution error weights at the current time.

      @noarkode <node> ERKStepGetErrWeights *)
  val get_err_weights : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

  (** Returns the vector of estimated local errors.

      @noarkode <node> ERKStepGetEstLocalErrors *)
  val get_est_local_errors : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

  (** Summaries of time-stepper statistics. *)
  type timestepper_stats = {
      exp_steps : int;
        (** Cumulative number of stability-limited solver steps. *)
      acc_steps : int;
        (** Cumulative number of accuracy-limited solver steps. *)
      step_attempts : int;
        (** Cumulative number of steps attempted by the solver. *)
      num_nf_evals : int;
        (** Number of calls to the right-hand side function. *)
      num_err_test_fails : int;
        (** Number of error test failures. *)
    }

  (** Returns a grouped set of time-stepper statistics.

      @noarkode <node> ERKStepGetTimestepperStats *)
  val get_timestepper_stats    : ('d, 'k) session -> timestepper_stats

  (** Returns a grouped set of integrator statistics.

      @noarkode <node> ERKStepGetStepStats *)
  val get_step_stats           : ('d, 'k) session -> step_stats

  (** Prints time-stepper statistics on the given channel.

      @noarkode <node> ERKStepGetTimestepperStats *)
  val print_timestepper_stats  : ('d, 'k) session -> out_channel -> unit

  (** Prints integrator statistics on the given channel.

      @noarkode <node> ERKStepGetStepStats *)
  val print_step_stats  : ('d, 'k) session -> out_channel -> unit

  (** {2:roots Additional root-finding functions} *)

  (** [set_root_direction s dir] specifies the direction of zero-crossings to
      be located and returned. [dir] may contain one entry for each root
      function.

      @noarkode <node> ERKStepSetRootDirection *)
  val set_root_direction : ('d, 'k) session -> RootDirs.d array -> unit

  (** Like {!set_root_direction} but specifies a single direction for all root
      functions.

      @noarkode <node> ERKStepSetRootDirection *)
  val set_all_root_directions : ('d, 'k) session -> RootDirs.d -> unit

  (** Disables issuing a warning if some root function appears to be
      identically zero at the beginning of the integration.

      @noarkode <node> ERKStepSetNoInactiveRootWarn *)
  val set_no_inactive_root_warn : ('d, 'k) session -> unit

  (** Returns the number of root functions. *)
  val get_num_roots : ('d, 'k) session -> int

  (** Fills an array showing which functions were found to have a root.

      @noarkode <node> ERKStepGetRootInfo *)
  val get_root_info : ('d, 'k) session -> Roots.t -> unit

  (** Returns the cumulative number of calls made to the user-supplied root
      function g.

      @noarkode <node> ERKStepGetNumGEvals *)
  val get_num_g_evals : ('d, 'k) session -> int

  (** Returns the cumulative number of test failures.

      @since 5.0.0
      @noarkode <node> ERKStepGetNumConstrFails *)
  val get_num_constr_fails : ('d, 'k) session -> int

  (** Outputs all the solver parameters to the given file.

      @since 4.1.0
      @noarkode <node> ERKStepWriteParameters *)
  val write_parameters : ('d, 'k) session -> Logfile.t -> unit

  (** Outputs the current butcher table to the given file.

      @since 4.1.0
      @noarkode <node> ERKStepWriteButcher *)
  val write_butcher : ('d, 'k) session -> Logfile.t -> unit

end (* }}} *)

(** MRIStep Time-Stepping Module for two-rate initial value problems.

    This module solves problems of the form
    {% $\dot{y} = f_s(t, y) + f_f(t, y)$%}, {% $y(t_0) = y_0$%}.

    Its interface is structured as follows.
    {ol
      {- {{:#solver}Solver initialization and use}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}}

    @noarkode <node> Using MRIStep for C and C++ Applications *)
module MRIStep : sig (* {{{ *)

  (** A session with the MRIStep time-stepping solver.

      An example session with MRIStep ({openfile arkode_mri_skel.ml}): {[
#include "../../examples/ocaml/skeletons/arkode_mri_skel.ml"
      ]}

      @noarkode <node> Skeleton of main program *)
  type ('d, 'k) session = ('d, 'k, Arkode_impl.mristep) Arkode_impl.session

  (** {2:solver Solver initialization and use} *)

  (** Creates and initializes a session with the solver. The call
      {[init inner f_s h_s ~roots:(nroots, g) t0 y0]}
      has as arguments:
      - [inner],  a session to use for the (fast) inner integrator,
      - [f_s],    the slow portion of the right-hand side function,
      - [h_s],    the slow step size,
      - [nroots], the number of root functions,
      - [g],      the root function ([(nroots, g)] defaults to {!no_roots}),
      - [t0],     the initial value of the independent variable, and
      - [y0],     a vector of initial values that also determines the number
                  of equations.

      This function does everything necessary to initialize a session, i.e.,
      it makes the calls referenced below. The {!solve_normal} and
      {!solve_one_step} functions may be called directly.

      If the inner (fast) stepper uses a fixed step size [h_f] that does not
      evenly divide the time interval between the stages of the outer (slow)
      method, then the actual value used for the fast steps will
      be slightly smaller than [h_f] to ensure that
      {% $(c_i^s - c_{i-1}^s)h_s/h_f$ %} is an integer value. Specifically,
      the fast step for the ith slow stage will be
      {% $h = \frac{(c_i^s - c_{i-1}^s)h_s}
                   {\left\lceil (c_i^s - c_{i-1}^s)h_s/h_f \right\rceil}$ %}.

      If fixed step sizes and a stop time are used, then the fixed step size
      will be used for all steps except the final step, which may be shorter.
      To resume use of the previous fixed step size, another call to
      {!set_fixed_step} is required before calling either {!solve_normal} or
      {!solve_one_step}.

      Neither Root finding nor non-identity mass matrices are supported in the
      inner session.

      @since 5.0.0
      @noarkode <node> MRIStepCreate
      @noarkode <node> MRIStepRootInit
      @noarkode <node> MRIStepSetFixedStep *)
  val init :
        ('data, 'kind) ARKStep.session
      -> 'data rhsfn
      -> float
      -> ?roots:(int * 'data rootsfn)
      -> float
      -> ('data, 'kind) Nvector.t
      -> ('data, 'kind) session

  (** Integrates an ODE system over an interval. The call
      [tret, r = solve_normal s tout yout] has as arguments
      - [s], a solver session,
      - [tout], the next time at which a solution is desired, and,
      - [yout], a vector to store the computed solution.

      It returns [tret], the time reached by the solver, which will be equal to
      [tout] if no errors occur, and, [r], a {!solver_result}.

      @noarkode <node> MRIStepEvolve (ARK_NORMAL)
      @raise IllInput Missing or illegal solver inputs.
      @raise TooMuchWork The requested time could not be reached in [mxstep] internal steps.
      @raise VectorOpErr A vector operation error occurred
      @raise InnerStepFail The inner stepper had an unrecoverable error *)
  val solve_normal : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                          -> float * solver_result

  (** Like {!solve_normal} but returns after one internal solver step.

      @noarkode <node> MRIStepEvolve (ARK_ONE_STEP) *)
  val solve_one_step : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                          -> float * solver_result

  (** Returns the interpolated solution or derivatives.
      [get_dky s dky t k] computes the [k]th derivative of the function
      at time [t], i.e.,
      {% $\frac{d^\mathtt{k}}{\mathit{dt}^\mathtt{k}}y(\mathtt{t})$%},
      and stores it in [dky]. The arguments must satisfy
      {% $t_n - h_n \leq \mathtt{t} \leq t_n$%}—where $t_n$
      denotes {!get_current_time} and $h_n$ denotes {!get_last_step},—
      and {% $0 \leq \mathtt{k} \leq 3$%}.

      This function may only be called after a successful return from either
      {!solve_normal} or {!solve_one_step}.

      @noarkode <node> MRIStepGetDky
      @raise BadT [t] is not in the interval {% $[t_n - h_n, t_n]$%}.
      @raise BadK [k] is not in the range \{0, 1, ..., dord\}. *)
  val get_dky : ('d, 'k) session -> ('d, 'k) Nvector.t -> float -> int -> unit

  (** Reinitializes the solver with new parameters and state values. The
      values of the independent variable, i.e., the simulation time, and the
      state variables must be given. If given, [roots] specifies a new root
      finding function. The new problem must have the same size as the
      previous one. The number of Runge Kutta stages for both the slow and
      fast methods in the new problem must not be larger than that of the
      previous one.

      @noarkode <node> MRIStepReInit
      @noarkode <node> MRIStepRootInit *)
  val reinit :
    ('d, 'k) session
    -> ?roots:(int * 'd rootsfn)
    -> float
    -> ('d, 'k) Nvector.t
    -> unit

  (** Resets the state to the given independent variable value and dependent
      variable vector. All previously set options, internal counter values,
      and step-size/error histories are retained.

      @since 5.4.0
      @noarkode <node> MRIStepReset *)
  val reset : ('d, 'k) session -> float -> ('d, 'k) Nvector.t

  (** Change the number of equations and unknowns between integrator steps.
      The call
      [resize s ~resize_nvec:rfn ynew t0]
      has as arguments:
      - [s], the solver session to resize,
      - [rfn], a resize function that transforms nvectors in place-otherwise
               they are simply destroyed and recloned,
      - [ynew], the newly-sized solution vector with the value {% $y(t_0)$%}, and
      - [t0], the current value of the independent variable $t_0$.

      @noarkode <node> MRIStepResize *)
  val resize :
    ('d, 'k) session
    -> ?resize_nvec:('d resize_fn)
    -> ('d, 'k) Nvector.t
    -> float
    -> unit

  (** {2:set Modifying the solver (optional input functions)} *)

  (** {3:mriset Optional inputs for MRIStep} *)

  (** Resets all optional input parameters to their default values. Neither
      the problem-defining functions nor the root-finding functions are
      changed.

      @noarkode <node> MRIStepSetDefaults *)
  val set_defaults : ('d, 'k) session -> unit

  (** Specifies the interpolation module used for output value interpolation
      and implicit method predictors.

      @since 5.2.0
      @noarkode <node> MRIStepSetInterpolantType *)
  val set_interpolant_type : ('d, 'k) session -> interpolant_type -> unit

  (** Specifies the degree of the polynomial interpolant used for output
      values and implicit method predictors.

      @since 5.2.0
      @noarkode <node> MRIStepSetInterpolantDegree *)
  val set_interpolant_degree : ('d, 'k) session -> int -> unit

  (** Write step adaptivity and solver diagnostics to the given file.

      @noarkode <node> MRIStepSetDiagnostics *)
  val set_diagnostics : ('d, 'k) session -> Logfile.t -> unit

  (** Do not write step adaptivity or solver diagnostics of a file.

      @noarkode <node> MRIStepSetDiagnostics *)
  val clear_diagnostics : ('d, 'k) session -> unit

  (** Configure the default error handler to write messages to a file.
      By default it writes to Logfile.stderr.

      @noarkode <node> MRIStepSetErrFile *)
  val set_error_file : ('d, 'k) session -> Logfile.t -> unit

  (** Specifies a custom function for handling error messages.
      The handler must not fail: any exceptions are trapped and discarded.

      @noarkode <node> MRIStepSetErrHandlerFn
      @noarkode <node> ARKErrHandlerFn *)
  val set_err_handler_fn
    : ('d, 'k) session -> (Util.error_details -> unit) -> unit

  (** Restores the default error handling function.

      @noarkode <node> MRIStepSetErrHandlerFn *)
  val clear_err_handler_fn : ('d, 'k) session -> unit

  (** Disables time step adaptivity and fix the step size for all internal
      steps. See the notes under {!init}.

      @noarkode <node> MRIStepSetFixedStep *)
  val set_fixed_step
    : ('d, 'k) session
      -> float
      -> unit

  (** Specifies the maximum number of messages warning that [t + h = t] on
      the next internal step.

      @noarkode <node> MRIStepSetMaxHnilWarns *)
  val set_max_hnil_warns : ('d, 'k) session -> int -> unit

  (** Specifies the maximum number of steps taken in attempting to reach
      a given output time.

      @noarkode <node> MRIStepSetMaxNumSteps *)
  val set_max_num_steps : ('d, 'k) session -> int -> unit

  (** Limits the value of the independent variable [t] when solving.
      By default no stop time is imposed.

      @noarkode <node> MRIStepSetStopTime *)
  val set_stop_time : ('d, 'k) session -> float -> unit

  (** A function to be called {e before} each inner integration. This function
      may be used, for instance, to perform communications or memory transfers
      of forcing data supplied by the outer integrator to the inner
      integrator.

      The first argument is the current value of the independent variable. The
      second is an array of outer forcing vectors.

      {warning The array and its contents should not be accessed after the
               function has returned}.

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      @noarkode <node> MRIStepPreInnerFn *)
  type 'd pre_inner_fn = float -> 'd array -> unit

  (** Set the function called {e before} each inner integration.

      @noarkode <node> MRIStepSetPreInnerFn *)
  val set_pre_inner_fn : ('d, 'k) session -> 'd pre_inner_fn -> unit

  (** Clear the function called {e before} each inner integration.

      @noarkode <node> MRIStepSetPreInnerFn *)
  val clear_pre_inner_fn : ('d, 'k) session -> unit

  (** A function to be called {e after} each inner integration. This function
      may be used, for instance, to perform communications or memory transfers
      of state data supplied by the inner integrator to the outer
      integrator.

      The first argument is the current value of the independent variable. The
      second is the current value of the dependent variable vector.

      {warning The dependent variable vector should not be accessed after the
               function has returned.}

      Within the function, raising a {!Sundials.RecoverableFailure} exception
      indicates a recoverable error. Any other exception is treated as an
      unrecoverable error.

      @noarkode <node> MRIStepPreInnerFn *)
  type 'd post_inner_fn = float -> 'd -> unit

  (** Set the function called {e after} each inner integration.

      @noarkode <node> MRIStepSetPostInnerFn *)
  val set_post_inner_fn : ('d, 'k) session -> 'd post_inner_fn -> unit

  (** Clear the function called {e after} each inner integration.

      @noarkode <node> MRIStepSetPostInnerFn *)
  val clear_post_inner_fn : ('d, 'k) session -> unit

  (** {3:mrisetivp Optional inputs for IVP method selection} *)

  (** Specifies a customized Butcher table for the outer (slow) method.
      The integer argument sets the global order of accuracy.

      @noarkode <node> MRIStepSetTable *)
  val set_table : ('d, 'k) session -> int -> ButcherTable.t -> unit

  (** Use a specific built-in Butcher table for the outer (slow) method.

      @noarkode <node> MRIStepSetMRITableNum *)
  val set_table_num : ('d, 'k) session -> ButcherTable.erk_table -> unit

  (** {3:mrisetadap Optional inputs for time step adaptivity} *)

  (** Set a post processing step function.

      @noarkode <node> MRISetPostprocessStepFn *)
  val set_postprocess_step_fn : ('d, 'k) session -> 'd postprocess_step_fn -> unit

  (** Clear the post processing step function.

      @noarkode <node> MRISetPostprocessStepFn *)
  val clear_postprocess_step_fn : ('d, 'k) session -> unit

  (** {2:get Querying the solver (optional output functions)} *)

  (** Returns the real and integer workspace sizes.

      @noarkode <node> MRIStepGetWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space          : ('d, 'k) session -> int * int

  (** Returns the cumulative number of internal steps taken by the solver.

      @noarkode <node> MRIStepGetNumSteps *)
  val get_num_steps           : ('d, 'k) session -> int

  (** Returns the integration step size taken on the last successful internal
      step.

      @noarkode <node> MRIStepGetLastStep *)
  val get_last_step           : ('d, 'k) session -> float

  (** Returns the number of calls to the (outer) right-hand side function.

      @noarkode <node> MRIStepGetNumRhsEvals *)
  val get_num_rhs_evals       : ('d, 'k) session -> int

  (** Returns the the current internal time reached by the solver.

      @noarkode <node> MRIStepGetCurrentTime *)
  val get_current_time        : ('d, 'k) session -> float

  (** Returns the current state vector. This vector provides direct access to
      the data within the integrator.

      @noarkode <node> MRIStepGetCurrentState *)
  val get_current_state : ('d, 'k) session -> 'd

  (** Internal data required to construct the current nonlinear implicit
      system within a nonlinear solver. *)
  type 'd nonlin_system_data = {
    tcur  : float;
      (** Independent variable value for slow stage {% $t^S_{n,i}$ %}. *)
    zpred : 'd;
      (** Predicted nonlinear solution {% $z_{\mathit{pred}}$ %}. This
          data must not be changed. *)
    zi    : 'd;
      (** Stage vector {% $z_i$ %}. This data may not be current and may
          need to be filled. *)
    fi    : 'd;
      (** Memory available for evaluating the slow right-hand side
          {% $f^S(t^S_{n,i}, z_i)$ %}. This data may not be current and may
          need to be filled. *)
    gamma : float;
      (** Current {% $\gamma$ %} for slow-stage calculation. *)
    sdata : 'd;
      (** Accumulated data from previous solution and stages
          {% $\tilde{a}_i$ %}. This data must not be changed. *)
  }

  (** Gives direct access to the internal data required to construct the
      current nonlinear implicit system within a nonlinear solver. This
      function should be called inside the nonlinear system function.
      If the nonlinear solver uses the [lsetup] or [lsolve] functions, then
      the nonlinear solver system function must fill the [zi] and [fi]
      vectors with, respectively, the current state and corresponding
      evaluation of the right-hand-side function:
      {% $z_i = z_{\mathit{pred}} + z_{\mathit{cor}}$ %} and
      {% $F_i = f^S(t^S_{n,i}, z_i)$ %} where {% $z_{\mathit{cor}}$ %} is the
      first argument of the nonlinear solver system function. Within a custom
      linear solver, then the vectors [zi] and [fi] are only current after
      an evaluation of the nonlinear system function.

      @since 5.4.0
      @noarkode <node> MRIStepGetNonlinSystemData *)
  val get_nonlin_system_data : ('d, 'k) session -> 'd nonlin_system_data

  (** Computes the current stage state vector using the stored prediction and
      the supplied correction from the nonlinear solver. The call
      [compute_state s zcor z] computes {% $z_i(t) = z_{\mathit{pred}}
      + z_{\mathit{cor}}$ %}.

      @since 5.4.0
      @noarkode <node> MRIStepComputeState *)
  val compute_state : ('d, 'k) session
                      -> ('d, 'k) Nvector.t
                      -> ('d, 'k) Nvector.t
                      -> unit


  (** {2:roots Additional root-finding functions} *)

  (** [set_root_direction s dir] specifies the direction of zero-crossings to
      be located and returned. [dir] may contain one entry for each root
      function.

      @noarkode <node> MRIStepSetRootDirection *)
  val set_root_direction : ('d, 'k) session -> RootDirs.d array -> unit

  (** Like {!set_root_direction} but specifies a single direction for all root
      functions.

      @noarkode <node> MRIStepSetRootDirection *)
  val set_all_root_directions : ('d, 'k) session -> RootDirs.d -> unit

  (** Disables issuing a warning if some root function appears to be
      identically zero at the beginning of the integration.

      @noarkode <node> MRIStepSetNoInactiveRootWarn *)
  val set_no_inactive_root_warn : ('d, 'k) session -> unit

  (** Returns the number of root functions. *)
  val get_num_roots : ('d, 'k) session -> int

  (** Fills an array showing which functions were found to have a root.

      @noarkode <node> MRIStepGetRootInfo *)
  val get_root_info : ('d, 'k) session -> Roots.t -> unit

  (** Returns the cumulative number of calls made to the user-supplied root
      function g.

      @noarkode <node> MRIStepGetNumGEvals *)
  val get_num_g_evals : ('d, 'k) session -> int

  (** Outputs all the solver parameters to the given file.

      @noarkode <node> MRIStepWriteParameters *)
  val write_parameters : ('d, 'k) session -> Logfile.t -> unit

end (* }}} *)

(** {2:exceptions Exceptions} *)

(** Raised on missing or illegal solver inputs. Also raised if an element
    of the error weight vector becomes zero during time stepping, or the
    linear solver initialization function failed, or a root was found both at
    [t] and very near [t].

 @noarkode <node> ARK_ILL_INPUT *)
exception IllInput

(** The initial and final times are too close to each other and an initial step
    size was not specified.

    @noarkode <node> ARK_TOO_CLOSE *)
exception TooClose

(** The requested time could not be reached in [mxstep] internal steps.
    See {!ARKStep.set_max_num_steps}, {!ERKStep.set_max_num_steps},
    or {!MRIStep.set_max_num_steps}.

    @noarkode <node> ARK_TOO_MUCH_WORK *)
exception TooMuchWork

(** The requested accuracy could not be satisfied.

    @noarkode <node> ARK_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** The inner stepper returned with an unrecoverable error.
    If possible, the exception in the inner stepper is specified.

    @noarkode <node> MRIStepGetLastInnerStepFlag
    @noarkode <node> ARK_INNER_STEP_FAILED *)
exception InnerStepFail of exn option

(** Too many error test failures within a step or at the minimum step size.
    See {!ARKStep.set_max_err_test_fails} or {!ERKStep.set_max_err_test_fails},
    and {!ARKStep.set_min_step} or {!ERKStep.set_min_step}.

    @noarkode <node> ARK_ERR_FAILURE *)
exception ErrFailure

(** Too many convergence test failures within a step or at the minimum step
    size. See {!ARKStep.set_max_conv_fails}, and {!ARKStep.set_min_step} or
    {!ERKStep.set_min_step}.

    @noarkode <node> ARK_CONV_FAILURE *)
exception ConvergenceFailure

(** Linear solver initialization failed.

    @noarkode <node> ARK_LINIT_FAIL *)
exception LinearInitFailure

(** Linear solver setup failed in an unrecoverable manner.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.PSetFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @noarkode <node> ARKStepGetLastLinFlag
    @noarkode <node> ARK_LSETUP_FAIL *)
exception LinearSetupFailure of exn option

(** Linear solver solution failed in an unrecoverable manner.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.ATimesFailure},
    {!Sundials_LinearSolver.PSolveFailure},
    {!Sundials_LinearSolver.GSFailure},
    {!Sundials_LinearSolver.QRSolFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @noarkode <node> ARKStepGetLastLinFlag
    @noarkode <node> ARK_LSOLVE_FAIL *)
exception LinearSolveFailure of exn option

(** Mass matrix solver initialization failed.

    @noarkode <node> ARK_MASSINIT_FAIL *)
exception MassInitFailure

(** Mass matrix solver setup failed in an unrecoverable manner.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.PSetFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @noarkode <node> ARKStepGetLastMassFlag
    @noarkode <node> ARK_MASSSETUP_FAIL *)
exception MassSetupFailure of exn option

(** Mass matrix solver solution failed in an unrecoverable manner.
    If possible, the exception in the underlying linear solver is specified.
    It is typically one of
    {!Sundials_LinearSolver.ZeroInDiagonal},
    {!Sundials_LinearSolver.ATimesFailure},
    {!Sundials_LinearSolver.PSolveFailure},
    {!Sundials_LinearSolver.GSFailure},
    {!Sundials_LinearSolver.QRSolFailure},
    or
    {!Sundials_LinearSolver.PackageFailure}.

    @noarkode <node> ARKStepGetLastMassFlag
    @noarkode <node> ARK_MASSSOLVE_FAIL *)
exception MassSolveFailure of exn option

(** Mass matrix-vector multiplication failed.

    @noarkode <node> ARK_MASSMULT_FAIL *)
exception MassMultFailure

(** The right-hand side function failed in an unrecoverable manner.

    @noarkode <node> ARK_RHSFUNC_FAIL *)
exception RhsFuncFailure

(** The right-hand side function had a recoverable error when first called.

    @noarkode <node> ARK_FIRST_RHSFUNC_ERR *)
exception FirstRhsFuncFailure

(** Too many convergence test failures, or unable to estimate the initial step
    size, due to repeated recoverable errors in the right-hand side function.

    @noarkode <node> ARK_REPTD_RHSFUNC_ERR *)
exception RepeatedRhsFuncFailure

(** The right-hand side function had a recoverable error, but no recovery was
    possible. This error can only occur after an error test failure at order
    one.

    @noarkode <node> ARK_UNREC_RHSFUNC_ERR *)
exception UnrecoverableRhsFuncFailure

(** Nonlinear solver initialization failed.

    @noarkode <node> ARK_NLS_INIT_FAIL *)
exception NonlinearInitFailure

(** Nonlinear solver setup failed in an unrecoverable manner.

    @noarkode <node> ARK_NLS_SETUP_FAIL *)
exception NonlinearSetupFailure

(** Nonlinear solver setup failed in a recoverable manner.

    @noarkode <node> ARK_NLS_SETUP_RECVR *)
exception NonlinearSetupRecoverable

(** The function requires a nonlinear solver, but one has not been configured.

    @noarkode <node> ARK_NLS_OP_ERR *)
exception NonlinearOperationError

(** The rootfinding function failed.

    @noarkode <node> ARK_RTFUNC_FAIL *)
exception RootFuncFailure

(** The postprocess step function failed.

    @noarkode <node> ARK_POSTPROCESS_FAIL *)
exception PostprocStepFailure

(** Raised by {!ARKStep.get_dky}, {!ERKStep.get_dky}, and
    {!MRIStep.get_dky} for invalid order values.

    @noarkode <node> ARKStepGetDky (ARK_BAD_K) *)
exception BadK

(** Raised by {!ARKStep.get_dky}, {!ERKStep.get_dky}, and
    {!MRIStep.get_dky} for invalid time values.

    @noarkode <node> ARKStepGetDky (ARK_BAD_T) *)
exception BadT

(** A fused vector operation failed.

    @noarkode <node> ARK_VECTOROP_ERR *)
exception VectorOpErr

