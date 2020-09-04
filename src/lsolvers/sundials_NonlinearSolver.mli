(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2020 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Generic nonlinear solvers.

    Sundials provides generic nonlinear solvers of two main types:
    {!module:Newton} and {!module:FixedPoint}. An instance of a nonlinear
    solver may only be associated with at most one integrator session at
    a time.

    This module supports calling both Sundials and custom OCaml nonlinear
    solvers from both Sundials integrators and OCaml applications. It is not
    possible, however, to use a single instance in both a Sundials integrator
    and an OCaml application, nor to override the solve and setup functions
    set by a Sundials integrator. Special support is provided for solving
    sensitivity problems.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNNonlinearSolver module
    @since 4.0.0 *)

open Sundials

(** A generic nonlinear solver.
    The type variables specify the {!Nvector.nvector} data (['data]) and
    kind (['kind]), and the type (['s]) of session data to be passed
    into the nonlinear solver and through to callbacks.

    @nocvode <node> Description of the SUNNonlinearSolver module
    @nocvode <node> SUNNonlinearSolver *)
type ('data, 'kind, 's) t
    = ('data, 'kind, 's) Sundials_NonlinearSolver_impl.nonlinear_solver

(** Signifies a nonlinear solver used from an OCaml application.
    Used as a ['s] argument to {!t} to signal an internal
    calling convention for {!sysfn}, {!lsetupfn}, and {!lsolvefn}, which
    are OCaml functions. Only such nonlinear solvers can be invoked
    directly from user code (in OCaml) to initialize, setup, and solve
    a problem. The {!convtestfn} may always be overridden. *)
type user = Sundials_NonlinearSolver_impl.user

(** Signifies a nonlinear solver used from a Sundials integrator.
    Used as a ['s] argument to {!t} to signal that
    direct calls are made to {!sysfn}, {!lsetupfn}, and {!lsolvefn}, which
    are C functions. Such nonlinear solvers cannot be reinvoked from OCaml.
    The ['a] argument specifies the integrator; a nonlinear solver used for
    one integrator cannot be reused with another. *)
type 'a integrator = 'a Sundials_NonlinearSolver_impl.integrator

(** A limited interface to arrays of {!Nvector.nvector}s sometimes required
    to apply nonlinear solvers to sensitivity problems. *)
module Senswrapper : sig (* {{{ *)

  (** A senswrapper is an {!Nvector.nvector} of {!Nvector.nvector}s that
      cannot be created or manipulated from OCaml. *)
  type ('d, 'k) t = ('d, 'k) Sundials_NonlinearSolver_impl.Senswrapper.t

  (** Creates an array to the nvector data within a senswrapper.
      Given [s1, s2 : ('d, 'k) t] where [s1 = s2], then [data s1]
      and [data s2] access the same underlying data. This fact can
      be exploited for caching the array in callbacks.

      @raise IncorrectUse Attempt to access an invalidated senswrapper *)
  val data : ('d, 'k) t -> 'd array

end (* }}} *)

(** {2:nlscore Core functions} *)

(** The problem specification expected by a nonlinear solver. *)
type nonlinear_solver_type =
  | RootFind   (** Solves {% $F(y) = 0$ %} *)
  | FixedPoint (** Solves {% $G(y) = y$ %} *)

(** Returns the type of a nonlinear solver.

    @nocvode <node> SUNNonlinSolGetType *)
val get_type : ('d, 'k, 's) t -> nonlinear_solver_type

(** Initializes a nonlinear solver.

    @nocvode <node> SUNNonlinSolInitialize *)
val init  : ('d, 'k, 's) t -> unit

(** Setup a nonlinear solver with an initial iteration value.

    @nocvode <node> SUNNonlinSolSetup *)
val setup :
  ('d, 'k, user) t
  -> y:('d, 'k) Nvector.t
  -> unit

(** Solves a nonlinear system.
    The call [solve ls ~y0 ~y ~w tol callLSetup] solves the
    nonlinear system {% $F(y) = 0$ %} or {% $G(y) = y$ %}, given
    an initial iterate [y0] (which must not be modified), a solution
    error-weight vector [w] used for computing weighted error norms,
    the requested solution tolerance in the weighted root-mean-squared norm
    [tol], and a flag [callLSetup] indicating whether the integrator
    recommends calling the setup function.

    @nocvode <node> SUNNonlinSolSolve *)
val solve :
  ('d, 'k, user) t
  -> y0:('d, 'k) Nvector.t
  ->  y:('d, 'k) Nvector.t
  ->  w:('d, 'k) Nvector.t
  -> float
  -> bool
  -> unit

(** {2:nlsset Set functions} *)

(** A function [sysfn y fg mem] to evaluate the nonlinear system
    {% $F(y)$ %} (for {{!t}RootFind})
    or {% $G(y)$ %} (for {{!t}FixedPoint}).
    The contents of [y] must not be modified.

    This function raises {!exception:Sundials.RecoverableFailure} to
    indicate a recoverable failure. Other exceptions signal unrecoverable
    failures.

    @nocvode <node> SUNNonlinSolSysFn *)
type ('data, 's) sysfn = 'data -> 'data -> 's -> unit

(** Specify a system function callback.
    The system function specifies the problem, either {% $F(y)$ %} or
    {% $G(y)$ %}.

    @nocvode <node> SUNNonlinSolSetSysFn *)
val set_sys_fn : ('d, 'k, user) t -> ('d, user) sysfn -> unit

(** A function to setup linear solves.
    For direct linear solvers, sets up the system {% $Ax = b$ %}
    where {% $A = \frac{\partial F}{\partial y}$ %} is the linearization
    of the nonlinear residual function {% $F(y) = 0$ %}. For iterative
    linear solvers, calls a preconditioner setup function.

    The call [jcur = lsetupfn y f jbad mem] has as arguments

    - [y], the state vector where the system should be set up,
    - [f], the value of the nonlinear system at [y],
    - [jbad], indicates if the solver believes that {% $A$ %} has gone stale,
       and
    - [mem], a token passed by the function provider.

    A true return value ([jcur]) signals that the Jacobian {% $A$ %} has been
    updated.

    This function raises {!exception:Sundials.RecoverableFailure} to
    indicate a recoverable failure. Other exceptions signal unrecoverable
    failures.

    @nocvode <node> SUNNonlinSolLSetupFn *)
type ('data, 's) lsetupfn = 'data -> 'data -> bool -> 's -> bool

(** Specify a linear solver setup callback.

    @nocvode <node> SUNNonlinSolSetLSetupFn *)
val set_lsetup_fn : ('d, 'k, user) t -> ('d, user) lsetupfn -> unit

(** A function to solve linear systems.
    Solves the system {% $Ax = b$ %} where
    {% $A = \frac{\partial F}{\partial y}$ %} is the linearization of the
    nonlinear residual function {% $F(y)= 0$ %}.

    The call [lsolvefn y b mem] has as arguments

    - [y], the input vector containing the current nonlinear iteration;
    - [b], on input: the right-hand-side vector for the linear solve,
           set on output to the solution {% $x$ %}; and,
    - [mem], a token passed by the function provider.

    This function raises {!exception:Sundials.RecoverableFailure} to
    indicate a recoverable failure. Other exceptions signal unrecoverable
    failures.

    @nocvode <node> SUNNonlinSolLSolveFn *)
type ('data, 's) lsolvefn = 'data -> 'data -> 's -> unit

(** Specify a linear solver callback.

    @nocvode <node> SUNNonlinSolSetLSolveFn *)
val set_lsolve_fn : ('d, 'k, user) t -> ('d, user) lsolvefn -> unit

(** Values returned by convergence tests.
    @nocvode <node> SUNNonlinSolConvTestFn *)
type convtest =
  | Success  (** Converged ([SUN_NLS_SUCCESS]) *)
  | Continue (** Not converged, keep iterating ([SUN_NLS_CONTINUE]) *)
  | Recover  (** Appears to diverge, try to recover ([SUN_NLS_CONV_RECVR]) *)

(** A function providing an integrator-specific convergence test.
    The call [convtestfn y del tol ewt mem] has as arguments

    - [y], the current nonlinear iterate,
    - [del], the difference between current and prior nonlinear iterates,
    - [tol], the nonlinear solver tolerance (in a weighted root-mean-squared
             norm with the given error-weight vector),
    - [ewt], the error-weight vector used in computing weighted norms, and,
    - [mem], a token passed by the function provider.

    @nocvode <node> SUNNonlinSolConvTestFn *)
type ('data, 's) convtestfn = 'data -> 'data -> float -> 'data -> 's -> convtest

(** Specify a convergence test callback for the nonlinear solver iteration.

    @nocvode <node> SUNNonlinSolSetConvTestFn *)
val set_convtest_fn : ('d, 'k, 's) t -> ('d, 's) convtestfn -> unit

(** Sets the maximum number of nonlinear solver iterations.

    @nocvode <node> SUNNonlinSolSetMaxIters *)
val set_max_iters : ('d, 'k, 's) t -> int -> unit

(** {2:nlsget Get functions} *)

(** Returns the total number of nonlinear solver iterations.

    @nocvode <node> SUNNonlinSolGetNumIters *)
val get_num_iters : ('d, 'k, 's) t -> int

(** Returns the iteration index of the current nonlinear solve.

    @nocvode <node> SUNNonlinSolGetCurIter *)
val get_cur_iter : ('d, 'k, 's) t -> int

(** Returns the total number of nonlinear solver convergence failures.

    @nocvode <node> SUNNonlinSolGetNumConvFails *)
val get_num_conv_fails : ('d, 'k, 's) t -> int

(** {2:nlsolvers Nonlinear Solver Implementations} *)

(** Generic nonlinear solver based on Newton's method.

    @nocvode <node> The SUNNonlinearSolver_Newton implementation *)
module Newton : sig (* {{{ *)

  (** Creates a nonlinear solver based on Newton's method.
      Solves nonlinear systems of the form {% $F(y) = 0$ %}.

      @nocvode <node> SUNNonlinSol_Newton *)
  val make : ('d, 'k) Nvector.t -> ('d, 'k, 's) t

  (** Creates a nonlinear solver based on Newton's method for
      sensitivity-enabled integrators.
      Solves nonlinear systems of the form {% $F(y) = 0$ %}.

      In the call [make_sens count y],

      - [count] is the number of vectors in the nonlinear problem,
        if there are {% $N_s$ %} sensitivities, then [count] should be
        {% $N_s + 1$ %} if using a simultaneous corrector
        or {% $N_s$ %} if using a staggered corrector; and,
      - [y] is a template for cloning vectors.

      @nocvode <node> SUNNonlinSol_NewtonSens *)
  val make_sens : int -> ('d, 'k) Nvector.t
    -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) t

  (** Returns the residual function that defines the nonlinear system.

      Raises [Invalid_argument] if called on a nonlinear solver that was not
      created by this module.

      @nocvode <node> SUNNonlinSolGetSysFn_Newton *)
  val get_sys_fn : ('d, 'k, 's) t -> (('d, 'k) Nvector.t, 's) sysfn option

end (* }}} *)

(** Generic nonlinear solver for fixed-point (functional) iteration with
    optional Anderson acceleration.

    @nocvode <node> The SUNNonlinearSolver_FixedPoint implementation *)
module FixedPoint : sig (* {{{ *)

  (** Creates a nonlinear solver using fixed-point (functional) iteration.
      Solves nonlinear systems of the form {% $G(y) = y$ %}.
      The number of [acceleration_vectors] defaults to zero.

      @nocvode <node> SUNNonlinSol_FixedPoint *)
  val make : ?acceleration_vectors:int
             -> ('d, 'k) Nvector.t
             -> ('d, 'k, 's) t

  (** Creates a nonlinear solver using fixed-point (functional) iteration for
      sensitivity-enabled integrators.
      Solves nonlinear systems of the form {% $G(y) = y$ %}.

      In the call [make_sens count y],
      - [count] is the number of vectors in the nonlinear problem,
        if there are {% $N_s$ %} sensitivities, then [count] should be
        {% $N_s + 1$ %} if using a simultaneous corrector
        or {% $N_s$ %} if using a staggered corrector;
      - [y] is a template for cloning vectors; and,

      The number of [acceleration_vectors] defaults to zero.

      @nocvode <node> SUNNonlinSol_FixedPointSens *)
  val make_sens : ?acceleration_vectors:int -> int -> ('d, 'k) Nvector.t
    -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) t

  (** Returns the residual function that defines the nonlinear system.

      Raises [Invalid_argument] if called on a nonlinear solver that was not
      created by this module.

      @nocvode <node> SUNNonlinSolGetSysFn_FixedPoint *)
  val get_sys_fn : ('d, 'k, 's) t -> (('d, 'k) Nvector.t, 's) sysfn option

end (* }}} *)

(** Custom nonlinear solvers.

    @nocvode <node> Implementing a Custom SUNNonlinearSolver Module *)
module Custom : sig (* {{{ *)

  (** Create a nonlinear solver from a set of callback functions.

      The callbacks should indicate failure by raising an exception (preferably
      one of the exceptions in this package). Raising
      {!exception:Sundials.RecoverableFailure} indicates a generic recoverable
      failure.

      The expected operations are:

      - [init]: initializes the nonlinear solver.

      - [setup]: sets up the nonlinear solver with an initial iteration value.

      - [set_lsetup_fn]: receive a linear solver setup callback.

      - [set_lsolve_fn]: receive a linear solver callback.

      - [set_convtest_fn]: receive a convergence test callback.

      - [set_max_iters]: sets the maximum number of iterations.

      - [get_num_iters]: returns the total number of iterations.

      - [get_cur_iter]: returns the iteration index of the current solve. This function is required when using a convergence test provided by Sundials or one of the spils linear solvers.

      - [get_num_conv_fails]: return the total number of convergence failures.

      - [nls_type]: the type of problem solved.

      - [solve]: the call [solve y0 y w tol callLSetup mem] should solve the nonlinear system {% $F(y) = 0$ %} or {% $G(y) = y$ %}, given the initial iterate [y0], which must not be modified, the solution error-weight vector [w] used for computing weighted error norms, the requested solution tolerance in the weighted root-mean-squared norm [tol], a flag [callLSetup] indicating whether the integrator recommends calling the setup function, and a memory value to be passed to the system function.

      - [set_sys_fn]: receive the system callback.
  *)
  val make :
       ?init               : (unit -> unit)
    -> ?setup              : (('d, 'k) Nvector.t -> 's -> unit)
    -> ?set_lsetup_fn      : ((('d, 'k) Nvector.t, 's) lsetupfn -> unit)
    -> ?set_lsolve_fn      : ((('d, 'k) Nvector.t, 's) lsolvefn -> unit)
    -> ?set_convtest_fn    : ((('d, 'k) Nvector.t, 's) convtestfn -> unit)
    -> ?set_max_iters      : (int  -> unit)
    -> ?get_num_iters      : (unit -> int)
    -> ?get_cur_iter       : (unit -> int)
    -> ?get_num_conv_fails : (unit -> int)
    -> nls_type            : nonlinear_solver_type
    -> solve               : (('d, 'k) Nvector.t
                              -> ('d, 'k) Nvector.t
                              -> ('d, 'k) Nvector.t
                              -> float -> bool -> 's -> unit)
    -> set_sys_fn          : ((('d, 'k) Nvector.t, 's) sysfn -> unit)
    -> ('d, 'k, 's) t

  (** Create a nonlinear solver from a set of callback functions for
      sensitivity problems that pass arrays of nvectors. As for the
      {!make} function except that the callbacks receive arrays of
      values.

      Writing custom nonlinear solvers for use with some forward sensitivity
      methods requires the "internal" senswrapper type.
      Any attempt to use {!Senswrapper.t}s outside of the call to
      setup or solve that provides them will result in an {!IncorrectUse}
      exception. They must only be used to extract the underlying data with
      {!Senswrapper.data} or as arguments for lsetup_fn, lsolve_fn,
      convtest_fn, and sys_fn. There are no restrictions on the arrays
      extracted with {!Senswrapper.data}. *)
  val make_sens :
       ?init               : (unit -> unit)
    -> ?setup              : (('d, 'k) Senswrapper.t -> 'a integrator -> unit)
    -> ?set_lsetup_fn      : ((('d, 'k) Senswrapper.t, 'a integrator) lsetupfn -> unit)
    -> ?set_lsolve_fn      : ((('d, 'k) Senswrapper.t, 'a integrator) lsolvefn -> unit)
    -> ?set_convtest_fn    : ((('d, 'k) Senswrapper.t, 'a integrator) convtestfn -> unit)
    -> ?set_max_iters      : (int  -> unit)
    -> ?get_num_iters      : (unit -> int)
    -> ?get_cur_iter       : (unit -> int)
    -> ?get_num_conv_fails : (unit -> int)
    -> nls_type            : nonlinear_solver_type
    -> solve               : (('d, 'k) Senswrapper.t
                              -> ('d, 'k) Senswrapper.t
                              -> ('d, 'k) Senswrapper.t
                              -> float -> bool -> 'a integrator -> unit)
    -> set_sys_fn          : ((('d, 'k) Senswrapper.t, 'a integrator) sysfn -> unit)
    -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) t

end (* }}} *)

(** {2:exceptions Exceptions} *)

(** An error occurred in a vector operation.
    {cconst SUN_NLS_VECTOROP_ERR} *)
exception VectorOpError

(** Raised when a nonlinear solver is used incorrectly.
    For example, calling {!solve} without having first called {!set_sys_fn}
    ([SUN_NLS_MEM_NULL]). *)
exception IncorrectUse

(** Raised on an attempt to associate a nonlinear solver instance with more
    than one session. *)
exception NonlinearSolverInUse

