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
(*               User Documentation for Arkode v2.6.0                   *)
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

(** Adaptive-step time integration for stiff, nonstiff, and multi-rate
    systems of ODE initial value problems with zero-crossing
    detection.

    This module solves numerically problems of the form
    {% $M\dot{y} = f_E(t, y) + f_I(t, y)$%}, {% $y(t_0) = y_0$%}.

    This documented interface is structured as follows.
    {ol
      {- {{:#linear}Linear and mass matrix solvers}}
      {- {{:#tols}Tolerances}}
      {- {{:#solver}Solver initialization and use}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}
      {- {{:#exceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

open Sundials

(** A session with the ARKODE solver.

    An example session with Arkode ({openfile arkode_skel.ml}): {[
#include "../../examples/ocaml/skeletons/arkode_skel.ml"
    ]}

    @noarkode <node> Skeleton of main program *)
type ('d, 'k) session = ('d, 'k) Arkode_impl.session

(** Alias for sessions based on serial nvectors. *)
type 'k serial_session = (Nvector_serial.data, 'k) session
                         constraint 'k = [>Nvector_serial.kind]

(** {2:linear Linear and mass matrix solvers} *)

(** Linear solvers used by Arkode.

    @noarkode <node> Linear Solver Specification Functions *)
type ('data, 'kind) session_linear_solver =
  ('data, 'kind) Arkode_impl.session_linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type 'kind serial_session_linear_solver =
  (Nvector_serial.data, 'kind) session_linear_solver
  constraint 'kind = [>Nvector_serial.kind]

(** Workspaces with three temporary vectors. *)
type 'd triple = 'd * 'd * 'd

(** Arguments common to Jacobian callback functions.

    @noarkode <node> ARKDlsDenseJacFn
    @noarkode <node> ARKDlsBandJacFn
    @noarkode <node> ARKSpilsJacTimesVecFn
    @noarkode <node> ARKSpilsPrecSolveFn
    @noarkode <node> ARKSpilsPrecSetupFn *)
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
  include module type of LinearSolver.Direct

  (** Callback functions that compute dense approximations to a Jacobian
      matrix. In the call [jac arg jm], [arg] is a {!jacobian_arg}
      with three work vectors and the computed Jacobian must be stored in [jm].

      The callback should load the [(i,j)]th entry of [jm] with
      {% $\partial (f_I)_i/\partial y_j$%}, i.e., the partial derivative of
      the [i]th implicit equation with respect to the [j]th variable,
      evaluated at the values of [t] and [y] obtained from [arg]. Only
      nonzero elements need be loaded into [jm].

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning Neither the elements of [arg] nor the matrix [jm] should
               be accessed after the function has returned.}

      @noarkode <node> ARKDlsJacFn *)
  type 'm jac_fn = (RealArray.t triple, RealArray.t) jacobian_arg
                   -> 'm -> unit

  (** Create an Arkode-specific linear solver from a Jacobian approximation
      function and a generic direct linear solver.
      The Jacobian approximation function is optional for dense and banded
      solvers (if not given an internal difference quotient approximation is
      used), but must be provided for other solvers (or {Invalid_argument} is
      raised).

      @nocvode <node> ARKDlsSetLinearSolver
      @nocvode <node> ARKDlsSetJacFn *)
  val solver :
    ?jac:'m jac_fn ->
    ('m, 'kind, 't) LinearSolver.Direct.serial_linear_solver ->
    'kind serial_session_linear_solver

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by a direct
      linear solver.

      @noarkode <node> ARKDlsGetWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space : 'k serial_session -> int * int

  (** Returns the number of calls made by a direct linear solver to the
      Jacobian approximation function.

      @noarkode <node> ARKDlsGetNumJacEvals *)
  val get_num_jac_evals : 'k serial_session -> int

  (** Returns the number of calls to the right-hand side callback due to
      the finite difference Jacobian approximation.

      @noarkode <node> ARKDlsGetNumRhsEvals *)
  val get_num_rhs_evals : 'k serial_session -> int

end (* }}} *)

(** Scaled Preconditioned Iterative Linear Solvers.

    @noarkode <node> Linear solver specification functions
    @noarkode <node> Iterative linear solvers optional input functions.
    @noarkode <node> Iterative linear solvers optional output functions. *)
module Spils : sig (* {{{ *)
  include module type of LinearSolver.Iterative

  (** {3:precond Preconditioners} *)

  (** Arguments passed to the preconditioner solver function.

      @noarkode <node> ARKSpilsPrecSolveFn *)
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

      @noarkode <node> ARKSpilsPrecSolveFn *)
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

      @noarkode <node> ARKSpilsSetPreconditioner
      @noarkode <node> ARKSpilsPrecSetupFn *)
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

      @noarkode <node> ARKSpilsSetPreconditioner
      @noarkode <node> ARKSpilsPrecSetupFn
      @noarkode <node> ARKSpilsPrecSolveFn *)
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

    (** Returns the number of calls to the right-hand side callback for the
        difference banded Jacobian approximation. This counter is only updated
        if the default difference quotient function is used.

        @noarkode <node> ARKBandPrecGetNumRhsEvals *)
    val get_num_rhs_evals : 'k serial_session -> int
  end (* }}} *)

  (** {3:lsolvers Solvers} *)

  (** Callback functions that preprocess or evaluate Jacobian-related data
      needed by the jac_times_vec_fn. In the call [jac_times_setup_fn arg],
      [arg] is a {!jacobian_arg} with no work vectors.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning The elements of [arg] should not be accessed after the
               function has returned.}

      @nocvode <node> ARKSpilsJacTimesSetupFn *)
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

      @noarkode <node> ARKSpilsJacTimesVecFn *)
  type 'd jac_times_vec_fn =
    ('d, 'd) jacobian_arg
    -> 'd (* v *)
    -> 'd (* Jv *)
    -> unit

    (** Create an Arkode-specific linear solver from a generic iterative
        linear solver.

        NB: a [jac_times_setup_fn] is not supported in
            {!Sundials.sundials_version} < 3.0.0.

        @nocvode <node> ARKSpilsSetLinearSolver
        @nocvode <node> ARKSpilsSetJacTimes *)
    val solver :
      ('d, 'k, 'f) LinearSolver.Iterative.linear_solver
      -> ?jac_times_vec:'d jac_times_setup_fn option * 'd jac_times_vec_fn
      -> ('d, 'k) preconditioner
      -> ('d, 'k) session_linear_solver

  (** {3:set Solver parameters} *)

  (** Sets the factor by which the Krylov linear solver's convergence test
      constant is reduced from the Newton iteration test constant.
      This factor must be >= 0; passing 0 specifies the default (0.05).

      @noarkode <node> ARKSpilsSetEpsLin *)
  val set_eps_lin : ('d, 'k) session -> float -> unit

  (** {3:stats Solver statistics} *)

  (** Returns the sizes of the real and integer workspaces used by the spils
      linear solver.

      @noarkode <node> ARKSpilsGetWorkSpace
      @return ([real_size], [integer_size]) *)
  val get_work_space       : ('d, 'k) session -> int * int

  (** Returns the cumulative number of linear iterations.

      @noarkode <node> ARKSpilsGetNumLinIters *)
  val get_num_lin_iters    : ('d, 'k) session -> int

  (** Returns the cumulative number of linear convergence failures.

      @noarkode <node> ARKSpilsGetNumConvFails *)
  val get_num_conv_fails   : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the setup function with
      [jok=false].

      @noarkode <node> ARKSpilsGetNumPrecEvals *)
  val get_num_prec_evals   : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the preconditioner solve
      function.

      @noarkode <node> ARKSpilsGetNumPrecSolves *)
  val get_num_prec_solves  : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the Jacobian-vector
      setup function.

      @since 3.0.0
      @nocvode <node> ARKSpilsGetNumJTSetupEvals *)
  val get_num_jtsetup_evals : ('d, 'k) session -> int

  (** Returns the cumulative number of calls to the Jacobian-vector
      function.

      @noarkode <node> ARKSpilsGetNumJtimesEvals *)
  val get_num_jtimes_evals : ('d, 'k) session -> int

  (** Returns the number of calls to the right-hand side callback for
      finite difference Jacobian-vector product approximation. This counter is
      only updated if the default difference quotient function is used.

      @noarkode <node> ARKSpilsGetNumRhsEvals *)
  val get_num_rhs_evals    : ('d, 'k) session -> int

  (** {3:lowlevel Low-level solver manipulation}

      The {!init} and {!reinit} functions are the preferred way to set or
      change preconditioner functions. These low-level functions are provided
      for experts who want to avoid resetting internal counters and other
      associated side-effects. *)

  (** Change the preconditioner functions.

      @noarkode <node> ARKSpilsSetPreconditioner
      @noarkode <node> ARKSpilsPrecSolveFn
      @noarkode <node> ARKSpilsPrecSetupFn *)
  val set_preconditioner :
    ('d, 'k) session
    -> ?setup:'d prec_setup_fn
    -> 'd prec_solve_fn
    -> unit

  (** Change the Jacobian-times-vector function.

      @noarkode <node> ARKSpilsSetJacTimesVecFn
      @noarkode <node> ARKSpilsJacTimesVecFn *)
  val set_jac_times :
    ('d, 'k) session
    -> ?jac_times_setup:'d jac_times_setup_fn
    -> 'd jac_times_vec_fn
    -> unit

  (** Remove a Jacobian-times-vector function and use the default
      implementation.

      @noarkode <node> ARKSpilsSetJacTimesVecFn
      @noarkode <node> ARKSpilsJacTimesVecFn *)
  val clear_jac_times : ('d, 'k) session -> unit

end (* }}} *)

(** Alternate Linear Solvers.

    @noarkode <node> Providing Alternate Linear Solver Modules *)
module Alternate : sig (* {{{ *)

  (** Indicates problems during the solution of nonlinear equations at a
      step. Helps decide whether to update the Jacobian data kept by a
      linear solver. *)
  type conv_fail =
    | NoFailures
        (** Either the first call for a step or the local error test
            failed on the previous attempt at this setup but the Newton
            iteration converged. {cconst ARK_NO_FAILURES} *)

    | FailBadJ
        (** The setup routine indicates that its Jacobian related data is not
            current and either the previous Newton corrector iteration did not
            converge, or, during the previous Newton corrector iteration, the
            linear solver's {{!callbacks}lsolve} routine failed in a
            recoverable manner. {cconst ARK_FAIL_BAD_J} *)

    | FailOther
        (** The previous Newton iteration failed to converge even
            though the linear solver was using current
            Jacobian-related data. {cconst ARK_FAIL_OTHER} *)

  (** Functions that initialize linear solver data, like counters and
      statistics.

      Raising any exception in this function (including
      {!Sundials.RecoverableFailure}) is treated as an unrecoverable error.

      @noarkode <node> linit *)
  type ('data, 'kind) linit = ('data, 'kind) session -> unit

  (** Functions that prepare the linear solver for subsequent calls
      to {{!callbacks}lsolve}.  This function must return [true]
      only if the Jacobian-related data is current after the call.

      This function may raise a {!Sundials.RecoverableFailure} exception to
      indicate that a recoverable error has occurred. Any other exception is
      treated as an unrecoverable error.

      {warning The vectors [ypred], [fpred], and those in [tmp] should not be
               accessed after the function returns.}

      @noarkode <node> lsetup *)
  type ('data, 'kind) lsetup =
    ('data, 'kind) session
    -> 'data lsetup_args
    -> bool

  (** Arguments to {!lsetup}. *)
  and 'data lsetup_args =
    {
      lsetup_conv_fail : conv_fail;
      (** Indicates that a problem occurred during the solution of
          the nonlinear equation at the current time step. *)

      lsetup_y : 'data;
      (** The predicted $y$ vector for the current internal step. *)

      lsetup_rhs : 'data;
      (** The value of the implicit right-hand side at the predicted $y$
          vector. *)

      lsetup_tmp : 'data triple;
      (** Temporary variables for use by the routine. *)
    }

  (** Functions that solve the linear equation $Ax = b$.
      $A$ is arises in the Newton iteration and gives some approximation to
      the Newton matrix {% $M-\gamma J$%} where
      {% $J = \frac{\partial}{\partial y}f_I(t_n, y_{\text{cur}})$%},
      $\gamma$ is available through {!get_gammas}.
      The call [lsolve s args b] has as arguments:
      - [s], the solver session,
      - [args], summarizing current approximations to the solution, and
      - [b], the right-hand side vector, also used for returning the result.

      Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
      Any other exception is treated as an unrecoverable error.

      {warning The vectors in {!lsolve_args} should not be accessed
               after the function returns.}

      @noarkode <node> lsolve *)
  type ('data, 'kind) lsolve =
    ('data, 'kind) session
    -> 'data lsolve_args
    -> 'data
    -> unit

  (** Arguments to {!lsolve}. *)
  and 'data lsolve_args =
    {
      lsolve_y : 'data;
      (** The solver's current approximation to $y(t_n)$. *)

      lsolve_rhs : 'data;
      (** A vector containing {% $f_I(t_n, y_{\text{cur}})$%}. *)
    }

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
        (('data, 'kind) session
          -> ('data, 'kind) Nvector.t
          -> ('data, 'kind) callbacks)
        -> ('data, 'kind) session_linear_solver

  (** {3:internals Solver internals} *)

  (** Internal values used in Newton iteration. *)
  type gammas = {
    gamma : float;    (** The current $\gamma$ value. *)
    gammap : float;   (** The value of $\gamma$ at the last setup call. *)
  }

  (** Returns the current and previous gamma values. *)
  val get_gammas : ('data, 'kind) session -> gammas
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

  (** {3:dlsmass Direct mass matrix solvers} *)
  module Dls : sig (* {{{ *)
    include module type of LinearSolver.Direct

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

        @noarkode <node> ARKDlsMassFn *)
    type 'm mass_fn = float -> RealArray.t triple -> 'm -> unit

    (** Create an Arkode-specific mass linear solver from a mass-matrix
        constructor function and a generic dense linear solver. The boolean
        argument indicates whether the mass matrix depends on the independent
        variable [t], if not it is only computed and factored once.

        NB: The boolean argument is ignored in
        {!Sundials.sundials_version} < 3.0.0.

        @nocvode <node> ARKDlsSetMassLinearSolver
        @nocvode <node> ARKDlsSetMassFn *)
    val solver :
      'm mass_fn
      -> bool
      -> ('m, 'kind, 't) serial_linear_solver
      -> 'kind serial_solver

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by a
        direct linear mass matrix solver.

        @noarkode <node> ARKDlsGetMassWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : 'k serial_session -> int * int

    (** Returns the number of calls made to the mass matrix solver setup
        routine.

        NB: This function is not supported by
        {!Sundials.sundials_version} < 3.0.0.

        @noarkode <node> ARKDlsGetNumMassSetups *)
    val get_num_setups : 'k serial_session -> int

    (** Returns the number of calls made to the mass matrix solver solve
        routine.

        @noarkode <node> ARKDlsGetNumMassSolves *)
    val get_num_solves : 'k serial_session -> int

    (** Returns the number of calls made to the mass matrix-times-vector
        routine.

        NB: This function is not supported by
        {!Sundials.sundials_version} < 3.0.0.

        @noarkode <node> ARKDlsGetNumMassMult *)
    val get_num_mult : 'k serial_session -> int

  end (* }}} *)

  (** {3:spilsmass Iterative mass matrix solvers} *)
  module Spils : sig (* {{{ *)
    include module type of LinearSolver.Iterative

    (** {3:precond Preconditioners} *)

    (** Arguments passed to the mass matrix preconditioner solver function.

        @noarkode <node> ARKSpilsMassPrecSolveFn *)
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

        @noarkode <node> ARKSpilsMassPrecSolveFn *)
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

        @noarkode <node> ARKSpilsMassPrecSetupFn *)
    type 'd prec_setup_fn =
         float
      -> unit

    (** Specifies a preconditioner, including the type of preconditioning
        (none, left, right, or both) and callback functions.
        The following functions construct preconditioners.

        @noarkode <node> ARKSpilsSetMassPreconditioner
        @noarkode <node> ARKSpilsMassPrecSetupFn
        @noarkode <node> ARKSpilsMassPrecSolveFn *)
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

    (** {3:lsolvers Solvers} *)

    (** Callback functions that preprocess or evaluate Jacobian-related data
        needed by the mass_times_vec_fn. The argument gives the independent
        variable [t].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        @nocvode <node> ARKSpilsMassTimesSetupFn *)
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

        @noarkode <node> ARKSpilsMassTimesVecFn *)
    type 'd mass_times_vec_fn =
         float (* t *)
      -> 'd    (* v *)
      -> 'd    (* Mv *)
      -> unit

    (** Create an Arkode-specific mass linear solver from a generic iterative
        linear solver. The boolean argument indicates whether the mass matrix
        depends on the independent variable [t], if not it is only computed
        and factored once.

        NB: a [mass_times_setup_fn] is not supported in
        {!Sundials.sundials_version} < 3.0.0.

        NB: The boolean argument is ignored in
        {!Sundials.sundials_version} < 3.0.0.

        @nocvode <node> ARKSpilsSetMassLinearSolver
        @nocvode <node> ARKSpilsSetMassTimes *)
    val solver :
      ('d, 'k, 'f) linear_solver
      -> ?mass_times_setup:mass_times_setup_fn
      -> 'd mass_times_vec_fn
      -> bool
      -> ('d, 'k) preconditioner
      -> ('d, 'k) solver

    (** {3:set Solver parameters} *)

    (** Sets the factor by which the Krylov linear solver's convergence
        test constant is reduced from the Newton iteration test constant.
        This factor must be >= 0; passing 0 specifies the default (0.05).

        @noarkode <node> ARKSpilsSetMassEpsLin *)
    val set_eps_lin : ('d, 'k) session -> float -> unit

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the
        spils linear solver.

        @noarkode <node> ARKSpilsGetMassWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space       : ('d, 'k) session -> int * int

    (** Returns the cumulative number of linear iterations.

        @noarkode <node> ARKSpilsGetNumMassIters *)
    val get_num_lin_iters    : ('d, 'k) session -> int

    (** Returns the cumulative number of linear convergence failures.

        @noarkode <node> ARKSpilsGetNumMassConvFails *)
    val get_num_conv_fails   : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the mass-matrix-vector
        setup function.

        @since 3.0.0
        @nocvode <node> ARKSpilsGetNumMTSetupEvals *)
    val get_num_mtsetup_evals : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the mass-matrix-vector
        product function ({!mass_times_vec_fn}).

        @noarkode <node> ARKSpilsGetNumMtimesEvals
        @since Sundials 2.6.3 *)
    val get_num_mtimes_evals : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the setup function with
        [jok=false].

        @noarkode <node> ARKSpilsGetNumMassPrecEvals *)
    val get_num_prec_evals   : ('d, 'k) session -> int

    (** Returns the cumulative number of calls to the preconditioner solve
        function.

        @noarkode <node> ARKSpilsGetNumMassPrecSolves *)
    val get_num_prec_solves  : ('d, 'k) session -> int

    (** {3:lowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change preconditioner functions. These low-level functions are
        provided for experts who want to avoid resetting internal counters
        and other associated side-effects. *)

    (** Change the preconditioner functions.

        @noarkode <node> ARKSpilsSetMassPreconditioner
        @noarkode <node> ARKSpilsMassPrecSolveFn
        @noarkode <node> ARKSpilsMassPrecSetupFn *)
    val set_preconditioner :
      ('d, 'k) session
      -> ?setup:'d prec_setup_fn
      -> 'd prec_solve_fn
      -> unit

    (** Change the mass matrix-times-vector function.

        @noarkode <node> ARKSpilsSetMassTimes
        @noarkode <node> ARKSpilsMassTimesVecFn *)
    val set_times :
      ('d, 'k) session
      -> ?mass_times_setup:mass_times_setup_fn
      -> 'd mass_times_vec_fn
      -> unit

  end (* }}} *)

  (** {3:altmass Alternate mass matrix solvers} *)
  module Alternate : sig (* {{{ *)

    (** Functions that initialize mass matrix solver data, like counters and
        statistics.

        Raising any exception in this function (including
        {!Sundials.RecoverableFailure}) is treated as an unrecoverable
        error.

        @noarkode <node> minit *)
    type ('data, 'kind) minit = ('data, 'kind) session -> unit

    (** Functions that prepare the mass matrix solver for subsequent calls
        to {{!callbacks}lsolve}.

        This function may raise a {!Sundials.RecoverableFailure} exception
        to indicate that a recoverable error has occurred. Any other
        exception is treated as an unrecoverable error.

        {warning The vectors [ypred], [fpred], and those in [tmp] should
                 not be accessed after the function returns.}

        @noarkode <node> msetup *)
    type ('data, 'kind) msetup =
      ('data, 'kind) session
      -> 'data triple
      -> unit

    (** Functions that solve the linear equation $Mx = b$. $M$ is the system
        mass matrix. The call [lsolve s b weight] has as arguments:
        - [s], the solver session,
        - [b], the right-hand side vector, also used for returning the
               result, and
        - [weight], a vector containing the error weights.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable
        error. Any other exception is treated as an unrecoverable error.

        {warning The vectors [b] and [weight] should not be accessed after
                 the function returns.}

        @noarkode <node> msolve
        @noarkode <node> Mass matrix solver *)
    type ('data, 'kind) msolve =
      ('data, 'kind) session
      -> 'data (* b *)
      -> unit

    (** The callbacks needed to implement an alternate linear solver. *)
    type ('data, 'kind) callbacks =
      {
        minit  : ('data, 'kind) minit option;
        msetup : ('data, 'kind) msetup option;
        msolve : ('data, 'kind) msolve;
      }

    (** Creates a mass matrix solver from a function returning a set of
        callbacks. The creation function is passed a session and a vector.
        The latter indicates the problem size and can, for example, be
        cloned. *)
    val solver :
          (('data, 'kind) session
            -> ('data, 'kind) Nvector.t
            -> ('data, 'kind) callbacks)
          -> ('data, 'kind) solver
  end (* }}} *)

end (* }}} *)

(** {2:tols Tolerances} *)

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

(** Functions that compute the weighted RMS residual weights. The call
    [rfun y rwt] takes the dependent variable vector [y] and fills the
    residual-weight vector [rwt] with positive values or raises
    {!Sundials.NonPositiveEwt}. Other exceptions are eventually propagated, but
    should be avoided ([ffun] is not allowed to abort the solver). *)
type 'data res_weight_fun = 'data -> 'data -> unit

type ('data, 'kind) res_tolerance =
  | ResStolerance of float
    (** [abs] : scalar absolute residual tolerance. *)
  | ResVtolerance of ('data, 'kind) Nvector.t
    (** [abs] : vector of absolute residual tolerances. *)
  | ResFtolerance of 'data res_weight_fun
    (** Compute the residual weight vector. *)

(** {2:solver Solver initialization and use} *)

(** Choice of method for solving non-linear systems that arise in the solution
    of implicit systems.

    @noarkode <node> Nonlinear solver methods
    @noarkode <node> ARKodeInit *)
type ('d, 'kind) iter =
  | Newton of ('d, 'kind) session_linear_solver
    (** Modified Newton iteration with a given linear solver. *)
  | FixedPoint of int
    (** Accelerated fixed-point solver. Specifies the number of vectors to
        store within the Anderson acceleration subspace. *)

(** The linearity of the implicit portion of the problem. *)
type linearity =
  | Linear of bool  (** Implicit portion is linear. Specifies whether
                        $f_I(t, y)$ or the preconditioner is time-dependent. *)
  | Nonlinear       (** Implicit portion is nonlinear. *)

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

type ('d, 'k) imex = {
    implicit: 'd rhsfn * ('d, 'k) iter * linearity;
    explicit: 'd rhsfn
  }

(** The form of the initial value problem. *)
type ('d, 'k) problem =
  | Implicit of 'd rhsfn * ('d, 'k) iter * linearity
      (** Diagonally Implicit Runge-Kutta (DIRK) solution of stiff problem. *)
  | Explicit of 'd rhsfn
      (** Explicit Runge-Kutta (ERK) solution of non-stiff problem. *)
  | ImEx of ('d, 'k) imex
      (** Additive Runge-Kutta (ARK) solution of multi-rate problem. *)

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

(** Creates and initializes a session with the solver. The call
    {[init problem tol ~restol ~order:ord ~mass:msolver ~roots:(nroots, g) t0 y0]}
    has as arguments:
    - [problem], specifies the problem to solve (see {!problem}),
    - [tol],     the integration tolerances,
    - [restol],  (optional) mass matrix residual tolerances,
    - [ord],     the order of accuracy for the integration method,
    - [msolver], optionally, a linear mass matrix solver,
    - [nroots],  the number of root functions,
    - [g],       the root function ([(nroots, g)] defaults to {!no_roots}),
    - [t0],     the initial value of the independent variable, and
    - [y0],     a vector of initial values that also determines the number
                of equations.

    The allowed values for [ord] are:
    - for explicit methods: {% $2 \le \mathtt{ord} \le 6$%},
    - for implicit methods: {% $2 \le \mathtt{ord} \le 5$%}, and
    - for imex methods: {% $3 \le \mathtt{ord} \le 5$%}.

    This function does everything necessary to initialize a session, i.e.,
    it makes the calls referenced below. The {!solve_normal} and
    {!solve_one_step} functions may be called directly.

    @noarkode <node> ARKodeCreate
    @noarkode <node> ARKodeInit
    @noarkode <node> ARKodeSetFixedPoint
    @noarkode <node> ARKodeSetLinear
    @noarkode <node> ARKodeSetNonlinear
    @noarkode <node> Linear solver specification functions
    @noarkode <node> Mass matrix solver specification functions
    @noarkode <node> ARKodeRootInit
    @noarkode <node> ARKodeSStolerances
    @noarkode <node> ARKodeSVtolerances
    @noarkode <node> ARKodeWFtolerances
    @noarkode <node> ARKodeResStolerance
    @noarkode <node> ARKodeResVtolerance
    @noarkode <node> ARKodeResFtolerance
    @noarkode <node> ARKodeSetOrder *)
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

(** Values returned by the step functions. Failures are indicated by
    exceptions.

    @noarkode <node> ARKode *)
type solver_result =
  | Success             (** The solution was advanced. {cconst ARK_SUCCESS} *)
  | RootsFound          (** A root was found. See {!get_root_info}.
                            {cconst ARK_ROOT_RETURN} *)
  | StopTimeReached     (** The stop time was reached. See {!set_stop_time}.
                            {cconst ARK_TSTOP_RETURN} *)

(** Integrates an ODE system over an interval. The call
    [tret, r = solve_normal s tout yout] has as arguments
    - [s], a solver session,
    - [tout], the next time at which a solution is desired, and,
    - [yout], a vector to store the computed solution.

    It returns [tret], the time reached by the solver, which will be equal to
    [tout] if no errors occur, and, [r], a {!solver_result}.

    @noarkode <node> ARKode (ARK_NORMAL)
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
    @raise PostprocStepFailure Failure in the Postprocess Step function.
    *)
val solve_normal : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                        -> float * solver_result

(** Like {!solve_normal} but returns after one internal solver step.

    @noarkode <node> ARKode (ARK_ONE_STEP) *)
val solve_one_step : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                        -> float * solver_result

(** Returns the interpolated solution or derivatives.
    [get_dky s dky t k] computes the [k]th derivative of the function at time
    [t], i.e., {% $\frac{d^\mathtt{k}}{\mathit{dt}^\mathtt{k}}y(\mathtt{t})$%},
    and stores it in [dky]. The arguments must satisfy
    {% $t_n - h_n \leq \mathtt{t} \leq t_n$%}—where $t_n$
    denotes {!get_current_time} and $h_n$ denotes {!get_last_step},—
    and {% $0 \leq \mathtt{k} \leq 3$%}.

    This function may only be called after a successful return from either
    {!solve_normal} or {!solve_one_step}.

    @noarkode <node> ARKodeGetDky
    @raise BadT [t] is not in the interval {% $[t_n - h_n, t_n]$%}.
    @raise BadK [k] is not in the range \{0, 1, 2, 3\}. *)
val get_dky : ('d, 'k) session -> ('d, 'k) Nvector.t -> float -> int -> unit

(** Reinitializes the solver with new parameters and state values. The
    values of the independent variable, i.e., the simulation time, and the
    state variables must be given. The new problem must have the same size
    as the previous one. Only the functions in the given problem description
    have any effect; the [iter] and [linear] fields are ignored.

    The allowed values for [order] are:
    - for explicit methods: {% $2 \le \mathtt{order} \le 6$%},
    - for implicit methods: {% $2 \le \mathtt{order} \le 5$%}, and
    - for imex methods: {% $3 \le \mathtt{order} \le 5$%}.

    @noarkode <node> ARKodeReInit
    @noarkode <node> ARKodeRootInit
    @noarkode <node> ARKodeSetOrder *)
val reinit :
  ('d, 'kind) session
  -> ?problem:('d, 'kind) problem
  -> ?order:int
  -> ?roots:(int * 'd rootsfn)
  -> float
  -> ('d, 'kind) Nvector.t
  -> unit

(** Called to resize a vector to match the dimensions of another. The call
    [resizefn y ytemplate] must resize [y] to match the size of [ytemplate].

    {warning [y] and [ytemplate] should not be accessed after the function
             has returned.}

    @noarkode <node> ARKVecResizeFn *)
type 'd resize_fn = 'd -> 'd -> unit

(** Change the number of equations and unknowns between integrator steps.
    The call [resize s ~resize_nvec:rfn ~linsolv:ls tol ~restol hscale ynew t0]
    has as arguments:
    - [s], the solver session to resize,
    - [rfn], a resize function that transforms nvectors in place-otherwise
             they are simply destroyed and recloned,
    - [ls], specify a different linear solver,
    - [tol], tolerance values (ensures that any tolerance vectors are resized),
    - [restol], (optional) mass matrix residual tolerances,
    - [hscale], the next step will be of size {% $h \mathtt{hscale}$%},
    - [ynew], the newly-sized solution vector with the value {% $y(t_0)$%}, and
    - [t0], the current value of the independent variable $t_0$.

    If a new linear solver is not specified, any existing linear solver in use
    is destroyed and reinitialized; settings must be reconfigured after
    resizing.

    If the mass matrix residual tolerance was previously set to
    {{!res_tolerance}ResVtolerance} and [restol] is not given, then it is reset
    to {{!res_tolerance}ResStolerance} with the default value.

    The [tol] argument is ignored in versions 2.6.1 and 2.6.2 since it may cause
    a segmentation error.

    @noarkode <node> ARKodeResize *)
val resize :
  ('d, 'kind) session
  -> ?resize_nvec:('d resize_fn)
  -> ?linsolv:(('d, 'kind) session_linear_solver)
  -> ('d, 'kind) tolerance
  -> ?restol:(('d, 'kind) res_tolerance)
  -> float
  -> ('d, 'kind) Nvector.t
  -> float
  -> unit

(** {2:set Modifying the solver (optional input functions)} *)

(** {3:setarkode Optional inputs for ARKode} *)

(** Sets the integration tolerances.

    @noarkode <node> ARKodeSStolerances
    @noarkode <node> ARKodeSVtolerances
    @noarkode <node> ARKodeWFtolerances *)
val set_tolerances : ('d, 'k) session -> ('d, 'k) tolerance -> unit

(** Sets the residual tolerance.

    @noarkode <node> ARKodeResStolerance
    @noarkode <node> ARKodeResVtolerance
    @noarkode <node> ARKodeResFtolerance *)
val set_res_tolerance : ('d, 'k) session -> ('d, 'k) res_tolerance -> unit

(** Resets all optional input parameters to their default values. Neither the
    problem-defining functions nor the root-finding functions are changed.

    @noarkode <node> ARKodeSetDefaults *)
val set_defaults : ('d, 'k) session -> unit

(** Specifies the order of accuracy for the polynomial interpolant used for
    dense output. The interpolant is used both for solution output values and
    implicit method predictors.

    @noarkode <node> ARKodeSetDenseOrder *)
val set_dense_order : ('d, 'k) session -> int -> unit

(** Write step adaptivity and solver diagnostics to the given file.

    @noarkode <node> ARKodeSetDiagnostics *)
val set_diagnostics : ('d, 'k) session -> Sundials.Logfile.t -> unit

(** Do not write step adaptivity or solver diagnostics of a file.

    @noarkode <node> ARKodeSetDiagnostics *)
val clear_diagnostics : ('d, 'k) session -> unit

(** Configure the default error handler to write messages to a file.
    By default it writes to Sundials.Logfile.stderr.

    @noarkode <node> ARKodeSetErrFile *)
val set_error_file : ('d, 'k) session -> Sundials.Logfile.t -> unit

(** Specifies a custom function for handling error messages.
    The handler must not fail: any exceptions are trapped and discarded.

    @noarkode <node> ARKodeSetErrHandlerFn
    @noarkode <node> ARKErrHandlerFn *)
val set_err_handler_fn : ('d, 'k) session -> (error_details -> unit) -> unit

(** Restores the default error handling function.

    @noarkode <node> ARKodeSetErrHandlerFn *)
val clear_err_handler_fn : ('d, 'k) session -> unit

(** Specifies the initial step size.

    @noarkode <node> ARKodeSetInitStep *)
val set_init_step : ('d, 'k) session -> float -> unit

(** Disables time step adaptivity and fix the step size for all internal steps.
    Pass [None] to restore time step adaptivity. This function is not
    recommended since there is no assurance of the validity of the computed
    solutions. It is provided primarily for code-to-code verification testing.
    Use in conjunction with {!set_min_step} or {!set_max_step}.

    @noarkode <node> ARKodeSetFixedStep *)
val set_fixed_step : ('d, 'k) session -> float option -> unit

(** Specifies the maximum number of messages warning that [t + h = t] on
    the next internal step.

    @noarkode <node> ARKodeSetMaxHnilWarns *)
val set_max_hnil_warns : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of steps taken in attempting to reach
    a given output time.

    @noarkode <node> ARKodeSetMaxNumSteps *)
val set_max_num_steps : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of error test failures permitted in attempting
    one step.

    @noarkode <node> ARKodeSetMaxErrTestFails *)
val set_max_err_test_fails : ('d, 'k) session -> int -> unit

(** Specifies a lower bound on the magnitude of the step size.

    @noarkode <node> ARKodeSetMinStep *)
val set_min_step : ('d, 'k) session -> float -> unit

(** Specifies an upper bound on the magnitude of the step size.

    @noarkode <node> ARKodeSetMaxStep *)
val set_max_step : ('d, 'k) session -> float -> unit

(** Sets all adaptivity and solver parameters to ‘best guess’ values. This
    routine takes into account the integration method (ERK, DIRK, or ARK) and a
    given method order; it should only be called after these have been set.

    @noarkode <node> ARKodeSetOptimalParams *)
val set_optimal_params : ('d, 'k) session -> unit

(** Limits the value of the independent variable [t] when solving.
    By default no stop time is imposed.

    @noarkode <node> ARKodeSetStopTime *)
val set_stop_time : ('d, 'k) session -> float -> unit

(** {3:setivp Optional inputs for IVP method selection} *)

(** Enables both the implicit and explicit portions of a problem.

    @raise IllInput If $f_I$ and $f_E$ are not already specified.
    @noarkode <node> ARKodeSetImEx *)
val set_imex : ('d, 'k) session -> unit

(** Disables the implicit portion of a problem.

    @raise IllInput If $f_E$ is not already specified.
    @noarkode <node> ARKodeSetExplicit *)
val set_explicit : ('d, 'k) session -> unit

(** Disables the explicit portion of a problem.

    @raise IllInput If $f_I$ is not already specified.
    @noarkode <node> ARKodeSetImplicit *)
val set_implicit : ('d, 'k) session -> unit

(** Parameters for the RK method. *)
type rk_method = {
    stages : int;                (** Number of stages ($s$). *)
    global_order : int;          (** Global order of accuracy ($q$). *)
    global_embedded_order : int; (** Global order of accuracy for the embedded
                                     RK method ($p$). *)
  }

(** Coefficients for the RK method. *)
type rk_timescoefs = {
    stage_times : RealArray.t;   (** Array (of length [stages]) of
                                     stage times ($c$). *)
    coefficients : RealArray.t;  (** Array (of length [stages]) of
                                     solution coefficients ($b$). *)
    bembed : RealArray.t option; (** Optional array (of length [stages]) of
                                     embedding coefficients ($b2$). *)
  }

(** Specifies a customized Butcher table pair for the additive RK method. The
    call [set_ark_tables rkm ai ae] specifies:
    - [rkm], the RK method parameters,
    - [ai], the coefficients defining the implicit RK stages (of length
            [rkm.stages * rkm.stages] in row-major order),
    - [ae], the coefficients defining the explicit RK stages (of length
            [rkm.stages * rkm.stages] in row-major order),
    - [i], the implicit stage times, solution coefficients, and embedding
          coefficients, and
    - [e], the explicit stage times, solution coefficients, and embedding
          coefficients.

    In versions of Sundials prior to 2.7.0, the [e] parameter is ignored; only
    the [i] parameter is used.

    If the [i.bembed] or the [e.bembed] field is [None] then the solver will
    run in fixed step mode and the step size must be set, or have been set,
    using either {!set_fixed_step} or {!set_init_step}. This feature is not
    available for Sundials versions prior to 2.7.0
    (the {!Sundials.NotImplementedBySundialsVersion} exception is raised).

    @raise IllInput If $f_I$ and $f_E$ are not already specified.
    @noarkode <node> ARKodeSetARKTables
    @noarkode <node> ARKodeSetImEx *)
val set_ark_tables
  : ('d, 'k) session
    -> rk_method
    -> RealArray.t -> RealArray.t
    -> rk_timescoefs -> rk_timescoefs
    -> unit

(** Specifies a customized Butcher table pair for the explicit portion of the
    system. The call [set_erk_table rkm ae] specifies:
    - [rkm], the RK method parameters,
    - [ae], the coefficients defining the explicit RK stages (of length
            [rkm.stages * rkm.stages] in row-major order), and
    - [e], the explicit stage times, solution coefficients, and embedding
          coefficients.

    If the [e.bembed] field is [None] then the solver will
    run in fixed step mode and the step size must be set, or have been set,
    using either {!set_fixed_step} or {!set_init_step}. This feature is not
    available for Sundials versions prior to 2.7.0
    (the {!Sundials.NotImplementedBySundialsVersion} exception is raised).

    @raise IllInput If $f_E$ is not already specified.
    @noarkode <node> ARKodeSetERKTable
    @noarkode <node> ARKodeSetExplicit *)
val set_erk_table
  : ('d, 'k) session -> rk_method -> RealArray.t -> rk_timescoefs -> unit

(** Specifies a customized Butcher table pair for the implicit portion of the
    system. The call [set_irk_table rkm ai] specifies:
    - [rkm], the RK method parameters,
    - [ai], the coefficients defining the implicit RK stages (of length
            [rkm.stages * rkm.stages] in row-major order), and
    - [i], the implicit stage times, solution coefficients, and embedding
          coefficients.

    If the [i.bembed] field is [None] then the solver will
    run in fixed step mode and the step size must be set, or have been set,
    using either {!set_fixed_step} or {!set_init_step}. This feature is not
    available for Sundials versions prior to 2.7.0
    (the {!Sundials.NotImplementedBySundialsVersion} exception is raised).

    @raise IllInput If $f_I$ is not already specified.
    @noarkode <node> ARKodeSetIRKTable
    @noarkode <node> ARKodeSetImplicit *)
val set_irk_table
  : ('d, 'k) session -> rk_method -> RealArray.t -> rk_timescoefs -> unit

(** Explicit Butcher tables

    @noarkode <node> Explicit Butcher tables *)
type erk_table =
  | HeunEuler_2_1_2       (** Default 2nd order explicit method (table 0). *)
  | BogackiShampine_4_2_3 (** Default 3rd order explicit method (table 1). *)
  | ARK_4_2_3_Explicit    (** Explicit portion of default 3rd order additive
                              method (table 2). *)
  | Zonneveld_5_3_4       (** Default 4th order explicit method (table 3). *)
  | ARK_6_3_4_Explicit    (** Explicit portion of default 3rd order additive
                              method (table 4). *)
  | SayfyAburub_6_3_4     (** Butcher table number 5. *)
  | CashKarp_6_4_5        (** Default 5th order explicit method (table 6). *)
  | Fehlberg_6_4_5        (** Butcher table number 7. *)
  | DormandPrince_7_4_5   (** Butcher table number 8. *)
  | ARK_8_4_5_Explicit    (** Explicit portion of default 5th order additive
                              method (table 9). *)
  | Verner_8_5_6          (** Default 6th order explicit method (table 10). *)
  | Fehlberg_13_7_8       (** Default 8th order explicit method (table 11). *)

(** Implicit Butcher tables

    @noarkode <node> Implicit Butcher tables *)
type irk_table =
  | SDIRK_2_1_2           (** Default 2nd order implicit method (table 12). *)
  | Billington_3_2_3      (** Butcher table number 13. *)
  | TRBDF2_3_2_3          (** Butcher table number 14. *)
  | Kvaerno_4_2_3         (** Butcher table number 15. *)
  | ARK_4_2_3_Implicit    (** Default 3rd order implicit method and the implicit
                              portion of the default 3rd order additive method
                              (table 16). *)
  | Cash_5_2_4            (** Butcher table number 17. *)
  | Cash_5_3_4            (** Butcher table number 18. *)
  | SDIRK_5_3_4           (** Default 4th order implicit method (table 19). *)
  | Kvaerno_5_3_4         (** Butcher table number 20. *)
  | ARK_6_3_4_Implicit    (** Implicit portion of the default 4th order additive
                              method (table 21). *)
  | Kvaerno_7_4_5         (** Butcher table number 22. *)
  | ARK_8_4_5_Implicit    (** Default 5th order method and the implicit portion
                              of the default 5th order additive method
                              (table 23). *)

(** Additive Butcher tables

    @noarkode <node> Additive Butcher tables *)
type ark_table =
  | ARK_4_2_3             (** 3rd-order pair combining tables 2 and 15. *)
  | ARK_6_3_4             (** 4th-order pair combining tables 4 and 20. *)
  | ARK_8_4_5             (** 5th-order pair combining tables 9 and 22. *)

(** Use specific built-in Butcher tables for an ImEx system.

    @raise IllInput If $f_I$ and $f_E$ are not already specified.
    @noarkode <node> ARKodeSetARKTableNum
    @noarkode <node> ARKodeSetImEx *)
val set_ark_table_num : ('d, 'k) session -> ark_table -> unit

(** Use specific built-in Butcher tables for an explicit integration of the
    problem.

    The {{!erk_table}Fehlberg_13_7_8} method is not available prior to
    Sundials 2.7.0.

    @raise IllInput If $f_E$ is not already specified.
    @noarkode <node> ARKodeSetERKTableNum
    @noarkode <node> ARKodeSetExplicit *)
val set_erk_table_num : ('d, 'k) session -> erk_table -> unit

(** Use specific built-in Butcher tables for an implicit integration of the
    problem.

    @raise IllInput If $f_I$ is not already specified.
    @noarkode <node> ARKodeSetIRKTableNum
    @noarkode <node> ARKodeSetImplicit *)
val set_irk_table_num : ('d, 'k) session -> irk_table -> unit

(** {3:setadap Optional inputs for time step adaptivity} *)

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
    stability based time steps, {!set_stability_fn} should be used.

    @noarkode <node> ARKAdaptFn *)
type 'd adaptivity_fn = float -> 'd -> adaptivity_args -> float

(** Parameters for the standard adaptivity algorithms.
    There are two:
    - [adaptivity_ks], the [k1], [k2], and [k3] parameters, or [None] to use
      the defaults, and
    - [adaptivity_method_order], [true] specifies the method order of
      accuracy $q$ and [false] specifies the embedding order of accuracy $p$. *)
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
  | AdaptivityFn of 'd adaptivity_fn (** A custom time-step adaptivity function. *)

(** Specifies the method and associated parameters used for time step
    adaptivity.

    @noarkode <node> ARKodeSetAdaptivityMethod
    @noarkode <node> ARKodeSetAdaptivityFn *)
val set_adaptivity_method : ('d, 'k) session -> 'd adaptivity_method -> unit

(** Specifies the fraction of the estimated explicitly stable step to use. Any
    non-positive argument resets to the default value (0.5).

    @noarkode <node> ARKodeSetCFLFraction *)
val set_cfl_fraction : ('d, 'k) session -> float -> unit

(** Specifies the bias to apply to the error estimates within accuracy-based
    adaptivity strategies.

    @noarkode <node> ARKodeSetErrorBias *)
val set_error_bias : ('d, 'k) session -> float -> unit

(** Specifies the step growth interval in which the step size will remain
    unchanged. In the call [set_fixed_step_bounds s lb ub], [lb] specifies a
    lower bound on the window to leave the step size fixed and [ub] specifies
    an upper bound. Any interval not containing 1.0 resets to default
    values.

    @noarkode <node> ARKodeSetFixedStepBounds *)
val set_fixed_step_bounds : ('d, 'k) session -> float -> float -> unit

(** Specifies the maximum step size growth factor upon a convergence failure on
    a stage solve within a step. Any value outside the interval {% $(0, 1]$%}
    resets to the default value.

    @noarkode <node> ARKodeSetMaxCFailGrowth *)
val set_max_cfail_growth : ('d, 'k) session -> float -> unit

(** Specifies the maximum step size growth factor upon multiple successive
    accuracy-based error failures in the solver. Any value outside the interval
    {% $(0, 1]$%} resets to the default value.

    @noarkode <node> ARKodeSetMaxEFailGrowth *)
val set_max_efail_growth : ('d, 'k) session -> float -> unit

(** Specifies the maximum allowed step size change following the very first
    integration step. Any value $\le 1$ resets to the default value.

    @noarkode <node> ARKodeSetMaxFirstGrowth *)
val set_max_first_growth : ('d, 'k) session -> float -> unit

(** Specifies the maximum growth of the step size between consecutive time
    steps. Any value $\le 1$ resets to the default value.

    @noarkode <node> ARKodeSetMaxGrowth *)
val set_max_growth : ('d, 'k) session -> float -> unit

(** Specifies the safety factor to be applied to the accuracy-based estimated
    step. Any non-positive value resets to the default value.

    @noarkode <node> ARKodeSetSafetyFactor *)
val set_safety_factor : ('d, 'k) session -> float -> unit

(** Specifies the threshold for “multiple” successive error failures before the
    factor from {!set_max_efail_growth} is applied. Any non-positive value
    resets to the default value.

    @noarkode <node> ARKodeSetSmallNumEFails *)
val set_small_num_efails : ('d, 'k) session -> float -> unit

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

(** Sets a problem-dependent function to estimate a stable time step size for
    the explicit portion of the ODE system.

    @noarkode <node> ARKodeSetStabilityFn *)
val set_stability_fn : ('d, 'k) session -> 'd stability_fn -> unit

(** Clears the problem-dependent function that estimates a stable time step
    size for the explicit portion of the ODE system.

    @noarkode <node> ARKodeSetStabilityFn *)
val clear_stability_fn : ('d, 'k) session -> unit

(** {3:setimplicit Optional inputs for implicit stage solves} *)

(** Solve the implicit portion of the problem using the accelerated fixed-point
    solver instead of the modified Newton iteration. The integer argument gives
    the maximum dimension of the accelerated subspace (i.e., the number of
    vectors to store within the Anderson acceleration subspace).

    @noarkode <node> ARKodeSetFixedPoint *)
val set_fixed_point : ('d, 'k) session -> int -> unit

(** Solve the implicit portion of the problem using the modified Newton solver.
    Overrides a previous call to {!set_fixed_point}.

    @noarkode <node> ARKodeSetNewton *)
val set_newton : ('d, 'k) session -> unit

(** Specifies that the implicit portion of the problem is linear. The flag
    indicates whether the Jacobian of {% $f_I(t,y)$%}, or, when using an
    iterative linear solver, the preconditioner is time-dependent ([true])
    or not ([false]).

    @noarkode <node> ARKodeSetLinear *)
val set_linear : ('d, 'k) session -> bool -> unit

(** Specifies that the implicit portion of the problem is nonlinear.

    @noarkode <node> ARKodeSetNonlinear *)
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

(** Specifies the method for predicting implicit solutions.

    @noarkode <node> ARKodeSetPredictorMethod *)
val set_predictor_method : ('d, 'k) session -> predictor_method -> unit

(** Specifies the maximum number of nonlinear solver iterations permitted per RK
    stage at each step.

    @noarkode <node> ARKodeSetMaxNonlinIters *)
val set_max_nonlin_iters : ('d, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver convergence failures
    permitted during one step.

    @noarkode <node> ARKodeSetMaxConvFails *)
val set_max_conv_fails : ('d, 'k) session -> int -> unit

(** Specifies the safety factor used in the nonlinear convergence test.

    @noarkode <node> ARKodeSetNonlinConvCoef *)
val set_nonlin_conv_coef : ('d, 'k) session -> float -> unit

(** Specifies the constant used in estimating the nonlinear solver convergence
    rate.

    @noarkode <node> ARKodeSetNonlinCRDown *)
val set_nonlin_crdown : ('d, 'k) session -> float -> unit

(** Specifies the nonlinear correction threshold beyond which the iteration
    will be declared divergent.

    @noarkode <node> ARKodeSetNonlinRDiv *)
val set_nonlin_rdiv : ('d, 'k) session -> float -> unit

(** Specifies a scaled step size ratio tolerance beyond which the linear solver
    setup routine will be signalled.

    @noarkode <node> ARKodeSetDeltaGammaMax *)
val set_delta_gamma_max : ('d, 'k) session -> float -> unit

(** Specifies the frequency of calls to the linear solver setup routine.
    Positive values specify the number of time steps between setup calls,
    negative values force recomputation at each Newton step, and zero values
    reset to the default.

    @noarkode <node> ARKodeSetMaxStepsBetweenLSet *)
val set_max_steps_between_lset : ('d, 'k) session -> int -> unit

(** A function to process the results of each timestep solution.
    The arguments are
    - [t], the value of the independent variable, and
    - [y], the value of the dependent variable vector {% $y(t)$%}.

    @noarkode <node> ARKPostprocessStepFn *)
type 'd postprocess_step_fn = float -> 'd -> unit

(** Set a post processing step function.

    @since 2.7.0
    @raise Sundials.NotImplementedBySundialsVersion Post processing not available
    @noarkode <node> ARKSetPostprocessStepFn *)
val set_postprocess_step_fn : ('d, 'k) session -> 'd postprocess_step_fn -> unit

(** Clear the post processing step function.

    @since 2.7.0
    @raise Sundials.NotImplementedBySundialsVersion Post processing not available
    @noarkode <node> ARKSetPostprocessStepFn *)
val clear_postprocess_step_fn : ('d, 'k) session -> unit

(** {2:get Querying the solver (optional output functions)} *)

(** Returns the real and integer workspace sizes.

    @noarkode <node> ARKodeGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space          : ('d, 'k) session -> int * int

(** Returns the cumulative number of internal steps taken by the solver.

    @noarkode <node> ARKodeGetNumSteps *)
val get_num_steps           : ('d, 'k) session -> int

(** Returns the cumulative number of stability-limited steps taken by the
    solver.

    @noarkode <node> ARKodeGetNumExpSteps *)
val get_num_exp_steps       : ('d, 'k) session -> int

(** Returns the cumulative number of accuracy-limited steps taken by the
    solver.

    @noarkode <node> ARKodeGetNumAccSteps *)
val get_num_acc_steps       : ('d, 'k) session -> int

(** Returns the cumulative number of steps attempted by the solver.

    @noarkode <node> ARKodeGetNumStepAttempts *)
val get_num_step_attempts   : ('d, 'k) session -> int

(** Returns the number of calls to the right-hand side functions.
    In the call [(nfe_evals, nfi_evals) = get_num_rhs_evals s],
    - [nfe_evals] is the number of calls to {% $f_E$%}, and
    - [nfi_evals] is the number of calls to {% $f_I$%}.

    @noarkode <node> ARKodeGetNumRhsEvals *)
val get_num_rhs_evals       : ('d, 'k) session -> int * int

(** Returns the number of local error test failures that have occurred.

    @noarkode <node> ARKodeGetNumErrTestFails *)
val get_num_err_test_fails  : ('d, 'k) session -> int

(** Returns the integration step size taken on the last successful internal
    step.

    @noarkode <node> ARKodeGetLastStep *)
val get_last_step           : ('d, 'k) session -> float

(** Returns the integration step size to be attempted on the next internal step.

    @noarkode <node> ARKodeGetCurrentStep *)
val get_current_step        : ('d, 'k) session -> float

(** Returns the the value of the integration step size used on the first step.

    @noarkode <node> ARKodeGetActualInitStep *)
val get_actual_init_step    : ('d, 'k) session -> float

(** Returns the the current internal time reached by the solver.

    @noarkode <node> ARKodeGetCurrentTime *)
val get_current_time        : ('d, 'k) session -> float

(** Returns the explicit and implicit Butcher tables in use by the solver.
    In the call [(rkm, ai, ae, i, e) = get_current_butcher_tables s],
    - [rkm] are the number of stages and global order accuracies,
    - [ai] are the coefficients of the DIRK method,
    - [ae] are the coefficients of the ERK method,
    - [i] are the implicit stage times, solution coefficients, and embedding
          coefficients, and
    - [e] are the explicit stage times, solution coefficients, and embedding
          coefficients.

    For versions of Sundials prior to 2.7.0, the [i] and [e] results are
    identical.

    @noarkode <node> ARKodeGetCurrentButcherTables *)
val get_current_butcher_tables
  : ('d, 'k) session
    -> rk_method * RealArray.t * RealArray.t * rk_timescoefs * rk_timescoefs

(** Returns a suggested factor by which the user's tolerances should be scaled
    when too much accuracy has been requested for some internal step.

    @noarkode <node> ARKodeGetTolScaleFactor *)
val get_tol_scale_factor : ('d, 'k) session -> float

(** Returns the solution error weights at the current time.

    @noarkode <node> ARKodeGetErrWeights *)
val get_err_weights : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Returns the vector of estimated local errors.

    @noarkode <node> ARKodeGetEstLocalErrors *)
val get_est_local_errors : ('d, 'k) session -> ('d, 'k) Nvector.t -> unit

(** Summaries of integrator statistics. *)
type integrator_stats = {
    num_steps : int;
      (** Cumulative number of internal solver steps. *)
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
    actual_init_step : float;
      (** Integration step sized used on the first step. *)
    last_step : float;
      (** Integration step size of the last internal step. *)
    current_step : float;
      (** Integration step size to attempt on the next internal step. *)
    current_time : float
      (** Current internal time reached by the solver. *)
  }

(** Returns the integrator statistics as a group.

    @noarkode <node> ARKodeGetIntegratorStats *)
val get_integrator_stats    : ('d, 'k) session -> integrator_stats

(** Prints the integrator statistics on the given channel.

    @noarkode <node> ARKodeGetIntegratorStats *)
val print_integrator_stats  : ('d, 'k) session -> out_channel -> unit

(** {3:setimplicit Implicit solver optional output functions} *)

(** Returns the number of calls made to the linear solver's setup function.

    @noarkode <node> ARKodeGetNumLinSolvSetups *)
val get_num_lin_solv_setups : ('d, 'k) session -> int

(** Returns the number of calls made to the mass matrix solver.

    @noarkode <node> ARKodeGetNumMassSolves *)
val get_num_mass_solves     : ('d, 'k) session -> int

(** Returns the number of nonlinear (functional or Newton) iterations performed.

    @noarkode <node> ARKodeGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('d, 'k) session -> int

(** Returns the number of nonlinear convergence failures that have occurred.

    @noarkode <node> ARKodeGetNumNonlinSolvConvFails *)
val get_num_nonlin_solv_conv_fails : ('d, 'k) session -> int

(** Returns both the numbers of nonlinear iterations performed [nniters] and
    nonlinear convergence failures [nncfails].

    @noarkode <node> ARKodeGetNonlinSolvStats
    @return ([nniters], [nncfails]) *)
val get_nonlin_solv_stats : ('d, 'k) session -> int * int

(** {2:roots Additional root-finding functions} *)

(** [set_root_direction s dir] specifies the direction of zero-crossings to be
    located and returned. [dir] may contain one entry for each root function.

    @noarkode <node> ARKodeSetRootDirection *)
val set_root_direction : ('d, 'k) session -> RootDirs.d array -> unit

(** Like {!set_root_direction} but specifies a single direction for all root
    functions.

    @noarkode <node> ARKodeSetRootDirection *)
val set_all_root_directions : ('d, 'k) session -> RootDirs.d -> unit

(** Disables issuing a warning if some root function appears to be identically
    zero at the beginning of the integration.

    @noarkode <node> ARKodeSetNoInactiveRootWarn *)
val set_no_inactive_root_warn : ('d, 'k) session -> unit

(** Returns the number of root functions. *)
val get_num_roots : ('d, 'k) session -> int

(** Fills an array showing which functions were found to have a root.

    @noarkode <node> ARKodeGetRootInfo *)
val get_root_info : ('d, 'k) session -> Roots.t -> unit

(** Returns the cumulative number of calls made to the user-supplied root
    function g.

    @noarkode <node> ARKodeGetNumGEvals *)
val get_num_g_evals : ('d, 'k) session -> int

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
    See {!set_max_num_steps}

    @noarkode <node> ARK_TOO_MUCH_WORK *)
exception TooMuchWork

(** The requested accuracy could not be satisfied.

    @noarkode <node> ARK_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** Too many error test failures within a step or at the minimum step size.
    See {!set_max_err_test_fails} and {!set_min_step}.

    @noarkode <node> ARK_ERR_FAILURE *)
exception ErrFailure

(** Too many convergence test failures within a step or at the minimum step
    size. See {!set_max_conv_fails} and {!set_min_step}.

    @noarkode <node> ARK_CONV_FAILURE *)
exception ConvergenceFailure

(** Linear solver initialization failed.

    @noarkode <node> ARK_LINIT_FAIL *)
exception LinearInitFailure

(** Linear solver setup failed in an unrecoverable manner.

    @noarkode <node> ARK_LSETUP_FAIL *)
exception LinearSetupFailure

(** Linear solver solution failed in an unrecoverable manner.

    @noarkode <node> ARK_LSOLVE_FAIL *)
exception LinearSolveFailure

(** Mass matrix solver initialization failed.

    @noarkode <node> ARK_MASSINIT_FAIL *)
exception MassInitFailure

(** Mass matrix solver setup failed in an unrecoverable manner.

    @noarkode <node> ARK_MASSSETUP_FAIL *)
exception MassSetupFailure

(** Mass matrix solver solution failed in an unrecoverable manner.

    @noarkode <node> ARK_MASSSOLVE_FAIL *)
exception MassSolveFailure

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

(** The rootfinding function failed.

    @noarkode <node> ARK_RTFUNC_FAIL *)
exception RootFuncFailure

(** The postprocess step function failed.

    @noarkode <node> ARK_POSTPROCESS_FAIL *)
exception PostprocStepFailure

(** Raised by {!get_dky} for invalid order values.

    @noarkode <node> ARKodeGetDky (ARK_BAD_K) *)
exception BadK

(** Raised by {!get_dky} for invalid time values.

    @noarkode <node> ARKodeGetDky (ARK_BAD_T) *)
exception BadT

