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
    solvers from both Sundials integrators and OCaml applications.

    This documentation is structured as follows.
    {ol
      {- {{:#nlscore}Core functions}}
      {- {{:#nlsset}Set functions}}
      {- {{:#nlsget}Get functions}}
      {- {{:#nlssolvers}Nonlinear Solver Implementations}}
      {- {{:#nlsexceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nonlinsol <SUNNonlinSol_API_link.html#the-sunnonlinearsolver-api> The SUNNonlinearSolver API
    @since 4.0.0 *)

open Sundials

(** A generic nonlinear solver.
    The type variables specify
    - ['data], the {!Nvector.nvector} data,
    - ['kind], the {!Nvector.nvector} kind,
    - ['s], the type of of session data to be passed
      into the nonlinear solver and through to callbacks, and,
    - ['v], a type indicating that the solver manipulates nvectors ([`Nvec])
      or {{!Senswrapper.t}senswrappers} ([`Sens]).

    @nonlinsol SUNNonlinearSolver *)
type ('data, 'kind, 's, 'v) t
    = ('data, 'kind, 's, 'v) Sundials_NonlinearSolver_impl.nonlinear_solver

(** A limited interface to arrays of {!Nvector.nvector}s required
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

    @nonlinsol SUNNonlinSolGetType *)
val get_type : ('d, 'k, 's, 'v) t -> nonlinear_solver_type

(** Initializes a nonlinear solver.

    @nonlinsol SUNNonlinSolInitialize *)
val init  : ('d, 'k, 's, 'v) t -> unit

(** Setup a nonlinear solver with an initial iteration value.

    @nonlinsol SUNNonlinSolSetup *)
val setup : ('d, 'k, 's, [`Nvec]) t -> y:('d, 'k) Nvector.t -> 's -> unit

(** Solves a nonlinear system.
    The call [solve ls ~y0 ~y ~w tol callLSetup s] solves the
    nonlinear system {% $F(y) = 0$ %} or {% $G(y) = y$ %}, given the following
    arguments.

    - [y0], a predicted value for the new solution state (which must not be
      modified),
    - [ycor], on input, an initial guess for the correction to the predicted
      states, and on output, the final correction to the predicted state,
    - [w], a solution error-weight vector used for computing weighted error
      norms,
    - [tol], the requested solution tolerance in the weighted
      root-mean-squared norm,
    - [callLSetup], a flag indicating whether the integrator recommends
      calling the setup function, and,
    - [s], the state to pass through to callbacks.

    @nonlinsol SUNNonlinSolSolve *)
val solve :
  ('d, 'k, 's, [`Nvec]) t
  ->  y0:('d, 'k) Nvector.t
  -> ycor:('d, 'k) Nvector.t
  ->   w:('d, 'k) Nvector.t
  -> float
  -> bool
  -> 's
  -> unit

(** {2:nlsset Set functions} *)

(** A function [sysfn y fg mem] to evaluate the nonlinear system
    {% $F(y)$ %} (for {{!t}RootFind})
    or {% $G(y)$ %} (for {{!t}FixedPoint}).
    The contents of [y] must not be modified.

    This function raises {!exception:Sundials.RecoverableFailure} to
    indicate a recoverable failure. Other exceptions signal unrecoverable
    failures.

    @nonlinsol SUNNonlinSolSysFn *)
type ('nv, 's) sysfn = 'nv -> 'nv -> 's -> unit

(** Specify a system function callback.
    The system function specifies the problem, either {% $F(y)$ %} or
    {% $G(y)$ %}.

    @nonlinsol SUNNonlinSolSetSysFn *)
val set_sys_fn : ('d, 'k, 's, [`Nvec]) t -> ('d, 's) sysfn -> unit

(** A function to setup linear solves.
    For direct linear solvers, sets up the system {% $Ax = b$ %}
    where {% $A = \frac{\partial F}{\partial y}$ %} is the linearization
    of the nonlinear residual function {% $F(y) = 0$ %}. For iterative
    linear solvers, calls a preconditioner setup function.

    The call [jcur = lsetupfn jbad mem] has as arguments [jbad], which
    indicates if the solver believes that {% $A$ %} has gone stale, and [mem],
    a token passed by the function provider. A true return value ([jcur])
    signals that the Jacobian {% $A$ %} has been updated.

    This function raises {!exception:Sundials.RecoverableFailure} to
    indicate a recoverable failure. Other exceptions signal unrecoverable
    failures.

    @nonlinsol SUNNonlinSolLSetupFn *)
type 's lsetupfn = bool -> 's -> bool

(** Specify a linear solver setup callback.

    @nonlinsol SUNNonlinSolSetLSetupFn *)
val set_lsetup_fn : ('d, 'k, 's, 'v) t -> 's lsetupfn -> unit

(** A function to solve linear systems.
    Solves the system {% $Ax = b$ %} where
    {% $A = \frac{\partial F}{\partial y}$ %} is the linearization of the
    nonlinear residual function {% $F(y)= 0$ %}.

    The call [lsolvefn b mem] has as arguments

    - [b], on input: the right-hand-side vector for the linear solve,
           set on output to the solution {% $x$ %}; and,
    - [mem], a token passed by the function provider.

    This function raises {!exception:Sundials.RecoverableFailure} to
    indicate a recoverable failure. Other exceptions signal unrecoverable
    failures.

    @nonlinsol SUNNonlinSolLSolveFn *)
type ('nv, 's) lsolvefn = 'nv -> 's -> unit

(** Specify a linear solver callback.

    @nonlinsol SUNNonlinSolSetLSolveFn *)
val set_lsolve_fn : ('d, 'k, 's, [`Nvec]) t -> ('d, 's) lsolvefn -> unit

(** Values returned by convergence tests.
    @nonlinsol SUNNonlinSolConvTestFn *)
type convtest =
  | Success  (** Converged ([SUN_NLS_SUCCESS]) *)
  | Continue (** Not converged, keep iterating ([SUN_NLS_CONTINUE]) *)
  | Recover  (** Appears to diverge, try to recover ([SUN_NLS_CONV_RECVR]) *)

(** A function providing a convergence test.
    The call [convtestfn y del tol ewt mem] has as arguments
    - [y],   the current nonlinear iterate,
    - [del], the difference between current and prior nonlinear iterates,
    - [tol], the nonlinear solver tolerance (in a weighted root-mean-squared
             norm with the given error-weight vector),
    - [ewt], the error-weight vector used in computing weighted norms, and,
    - [mem], a token passed by the function provider.

    @nonlinsol SUNNonlinSolConvTestFn *)
type ('nv, 's) convtestfn' = 'nv -> 'nv -> float -> 'nv -> 's -> convtest

(** A convergence test callback provided by an integrator.
    Such callbacks require an additional first argument, the nonlinear solver
    invoking the function, and otherwise expect nvector arguments.
    They access the linear solver and nvector arguments using generic
    functions, which is why the type variables are universally quantified. *)
type 's convtest_callback =
  { f : 'd1 'k1 't2 'd2 'k2. ('d1, 'k1, 't2, [`Nvec]) t
      -> (('d2, 'k2) Nvector.t, 's) convtestfn' }
  [@@unboxed]

(** A convergence test callback provided by an integrator with sensitivities.
    Such callbacks require an additional first argument, the nonlinear solver
    invoking the function, and otherwise expect senswrapper arguments.
    They access the linear solver and senswrapper arguments using generic
    functions, which is why the type variables are universally quantified. *)
type 's convtest_callback_sens =
  { f : 'd1 'k1 't2 'd2 'k2. ('d1, 'k1, 't2, [`Sens]) t
      -> (('d2, 'k2) Senswrapper.t, 's) convtestfn' }
  [@@unboxed]

(** A convergence test provided either by an integrator or a user program.
    The OCaml interface distinguishes callback functions set by the
    underlying library ([CConvTest]) from those supplied by user programs
    ([OConvTest]). This reflects the different underlying mechanisms used
    to create and invoke such functions. Callback functions provied by the
    underlying library can be invoked with any kind of linear solver and
    (homogeneous) nvectors since they manipulate these values generically. *)
type ('nv, 's, 'v) convtestfn =
  | CConvTest
    : 's convtest_callback cfun -> ('nv, 's, [`Nvec]) convtestfn
  | CSensConvTest
    : 's convtest_callback_sens cfun -> ('nv, 's, [`Sens]) convtestfn
  | OConvTest of ('nv, 's) convtestfn'

(** Ignore the nvector type argument in a convtestfn.

    @raise Invalid_argument if the value was constructed with [OConvTest] *)
val assert_not_oconvtestfn
  : ('nv1, 's, [`Nvec]) convtestfn -> ('nv2, 's, [`Nvec]) convtestfn

(** Specify a convergence test callback for the nonlinear solver iteration.

    @nonlinsol SUNNonlinSolSetConvTestFn *)
val set_convtest_fn :
  ('d, 'k, 's, [`Nvec]) t -> ('d, 's, [`Nvec]) convtestfn -> unit

(** Support for nonlinear solvers with sensitivities. *)
module Sens : sig (* {{{ *)

  (** Setup a nonlinear solver for sensitivities with an initial iteration
      value. See {!setup}.

      @nonlinsol SUNNonlinSolSetup *)
  val setup :
    ('d, 'k, 's, [`Sens]) t -> y:('d, 'k) Senswrapper.t -> 's -> unit

  (** Solves a nonlinear system with sensitivities. See {!solve}.

      @nonlinsol SUNNonlinSolSolve *)
  val solve :
    ('d, 'k, 's, [`Sens]) t
    ->  y0:('d, 'k) Senswrapper.t
    -> ycor:('d, 'k) Senswrapper.t
    ->   w:('d, 'k) Senswrapper.t
    -> float
    -> bool
    -> 's
    -> unit

  (** Specify a system function callback with sensitivities.

      @nonlinsol SUNNonlinSolSetSysFn *)
  val set_sys_fn :
    ('d, 'k, 's, [`Sens]) t -> (('d, 'k) Senswrapper.t, 's) sysfn -> unit

  (** Specify a linear solver callback with sensitivities.
      See {!set_lsolve_fn}.

      @nonlinsol SUNNonlinSolSetLSolveFn *)
  val set_lsolve_fn :
    ('d, 'k, 's, [`Sens]) t -> (('d, 'k) Senswrapper.t, 's) lsolvefn -> unit

  (** Ignore the nvector type argument in a convtestfn.

      @raise Invalid_argument if the value was constructed with [OConvTest] *)
  val assert_not_oconvtestfn
    : ('nv1, 's, [`Sens]) convtestfn -> ('nv2, 's, [`Sens]) convtestfn

  (** Specify a convergence test callback for the nonlinear solver iteration
      when using sensitivities. See {!set_convtest_fn}.

      @nonlinsol SUNNonlinSolSetConvTestFn *)
  val set_convtest_fn :
       ('d, 'k, 's, [`Sens]) t
    -> (('d, 'k) Senswrapper.t, 's, [`Sens]) convtestfn
    -> unit

end (* }}} *)

(** Sets the maximum number of nonlinear solver iterations.

    @nonlinsol SUNNonlinSolSetMaxIters *)
val set_max_iters : ('d, 'k, 's, 'v) t -> int -> unit

(** {2:nlsget Get functions} *)

(** Returns the number of nonlinear solver iterations in the most recent solve.

    @nonlinsol SUNNonlinSolGetNumIters *)
val get_num_iters : ('d, 'k, 's, 'v) t -> int

(** Returns the iteration index of the current nonlinear solve.

    @nonlinsol SUNNonlinSolGetCurIter *)
val get_cur_iter : ('d, 'k, 's, 'v) t -> int

(** Returns the number of nonlinear solver convergence failures in the most
    recent solve.

    @nonlinsol SUNNonlinSolGetNumConvFails *)
val get_num_conv_fails : ('d, 'k, 's, 'v) t -> int

(** {2:nlssolvers Nonlinear Solver Implementations} *)

(** Generic nonlinear solver based on Newton's method.

    @nonlinsol <SUNNonlinSol_links.html#the-sunnonlinsol-newton-implementation> The SUNNonlinearSolver_Newton implementation *)
module Newton : sig (* {{{ *)

  (** Creates a nonlinear solver based on Newton's method.
      Solves nonlinear systems of the form {% $F(y) = 0$ %}.

      @nonlinsol_module SUNNonlinSol_Newton *)
  val make :
       ?context:Context.t
    -> ('d, 'k) Nvector.t
    -> ('d, 'k, 's, [`Nvec]) t

  (** Creates a nonlinear solver based on Newton's method for
      sensitivity-enabled integrators.
      Solves nonlinear systems of the form {% $F(y) = 0$ %}.

      In the call [make_sens count y],

      - [count] is the number of vectors in the nonlinear problem,
        if there are {% $N_s$ %} sensitivities, then [count] should be
        {% $N_s + 1$ %} if using a simultaneous corrector
        or {% $N_s$ %} if using a staggered corrector; and,
      - [y] is a template for cloning vectors.

      @nonlinsol_module SUNNonlinSol_Newton *)
  val make_sens :
       ?context:Context.t
    -> int
    -> ('d, 'k) Nvector.t
    -> ('d, 'k, 's, [`Sens]) t

  (** Returns the residual function that defines the nonlinear system.

      Raises [Invalid_argument] if called on a nonlinear solver that was not
      created by this module.

      @nonlinsol_module SUNNonlinSolGetSysFn_Newton *)
  val get_sys_fn
    : ('d, 'k, 's, [`Nvec]) t -> (('d, 'k) Nvector.t, 's) sysfn option

end (* }}} *)

(** Generic nonlinear solver for fixed-point (functional) iteration with
    optional Anderson acceleration.

    @nonlinsol <SUNNonlinSol_links.html#the-sunnonlinsol-fixedpoint-implementation> The SUNNonlinearSolver_FixedPoint implementation *)
module FixedPoint : sig (* {{{ *)

  (** Creates a nonlinear solver using fixed-point (functional) iteration.
      Solves nonlinear systems of the form {% $G(y) = y$ %}.
      The number of [acceleration_vectors] defaults to zero.

      @nonlinsol_module SUNNonlinSol_FixedPoint *)
  val make :
       ?context:Context.t
    -> ?acceleration_vectors:int
    -> ('d, 'k) Nvector.t
    -> ('d, 'k, 's, [`Nvec]) t

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

      @nonlinsol_module SUNNonlinSol_FixedPoint *)
  val make_sens :
       ?context:Context.t
    -> ?acceleration_vectors:int
    -> int
    -> ('d, 'k) Nvector.t
    -> ('d, 'k, 's, [`Sens]) t

  (** Returns the residual function that defines the nonlinear system.

      Raises [Invalid_argument] if called on a nonlinear solver that was not
      created by this module.

      @nonlinsol_module SUNNonlinSolGetSysFn_FixedPoint *)
  val get_sys_fn
    : ('d, 'k, 's, [`Nvec]) t -> (('d, 'k) Nvector.t, 's) sysfn option

  (** Sets the damping parameter {% $\beta$ %} to use with Anderson
      acceleration. Damping is disabled by default {% $\beta = 1.0$ %}.

      @nonlinsol_module SUNNonlinSolSetDamping_FixedPoint
      @since 5.1.0 *)
  val set_damping : ('d, 'k, 's, 'v) t -> float -> unit

end (* }}} *)

(** Custom nonlinear solvers.

    @nonlinsol <SUNNonlinSol_API_link.html#implementing-a-custom-sunnonlinearsolver-module> Implementing a Custom SUNNonlinearSolver Module *)
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

      - [get_num_iters]: returns the number of iterations in the most recent
                         solve.

      - [get_cur_iter]: returns the iteration index of the current solve. This function is required when using a convergence test provided by Sundials or one of the spils linear solvers.

      - [get_num_conv_fails]: return the number of convergence failures in the
                              most recent solve.

      - [nls_type]: the type of problem solved.

      - [solve]: the call [solve y0 y w tol callLSetup mem] should solve the nonlinear system {% $F(y) = 0$ %} or {% $G(y) = y$ %}, given the initial iterate [y0], which must not be modified, the solution error-weight vector [w] used for computing weighted error norms, the requested solution tolerance in the weighted root-mean-squared norm [tol], a flag [callLSetup] indicating whether the integrator recommends calling the setup function, and a memory value to be passed to the system function.

      - [set_sys_fn]: receive the system callback.

      Note that the [setup] and [solve] functions are passed the payload data
      directly, whereas the [lsolvefn] and [sysfn]s require
      the data to be wrapped in an nvector. This asymmetry is awkward but,
      unfortunately, unavoidable given the implementation of nvectors and the
      different constraints for C-to-OCaml calls and OCaml-to-C calls. *)
  val make :
       ?init               : (unit -> unit)
    -> ?setup              : ('d -> 's -> unit)
    -> ?set_lsetup_fn      : ('s lsetupfn -> unit)
    -> ?set_lsolve_fn      : ((('d, 'k) Nvector.t, 's) lsolvefn -> unit)
    -> ?set_convtest_fn    : (('d, 's, [`Nvec]) convtestfn -> unit)
    -> ?set_max_iters      : (int  -> unit)
    -> ?get_num_iters      : (unit -> int)
    -> ?get_cur_iter       : (unit -> int)
    -> ?get_num_conv_fails : (unit -> int)
    -> nls_type            : nonlinear_solver_type
    -> solve               : ('d -> 'd -> 'd -> float -> bool -> 's -> unit)
    -> set_sys_fn          : ((('d, 'k) Nvector.t, 's) sysfn -> unit)
    -> ?context:Context.t
    -> unit
    -> ('d, 'k, 's, [`Nvec]) t

  (** Create a nonlinear solver from a set of callback functions for
      sensitivity problems that pass arrays of nvectors. As for the
      {!make} function except that the callbacks receive arrays of
      values.

      Writing custom nonlinear solvers for use with some forward sensitivity
      methods requires the "internal" senswrapper type.
      Any attempt to use {!Senswrapper.t}s outside of the call to
      setup or solve that provides them will result in an {!IncorrectUse}
      exception. They must only be used to extract the underlying data with
      {!Senswrapper.data} or as arguments for lsolve_fn, convtest_fn, or sys_fn.
      There are no restrictions on the arrays extracted with
      {!Senswrapper.data}. *)
  val make_sens :
       ?init               : (unit -> unit)
    -> ?setup              : (('d, 'k) Senswrapper.t -> 's -> unit)
    -> ?set_lsetup_fn      : ('s lsetupfn -> unit)
    -> ?set_lsolve_fn      : ((('d, 'k) Senswrapper.t, 's) lsolvefn -> unit)
    -> ?set_convtest_fn    : ((('d, 'k) Senswrapper.t, 's, [`Sens]) convtestfn -> unit)
    -> ?set_max_iters      : (int  -> unit)
    -> ?get_num_iters      : (unit -> int)
    -> ?get_cur_iter       : (unit -> int)
    -> ?get_num_conv_fails : (unit -> int)
    -> nls_type            : nonlinear_solver_type
    -> solve               : (('d, 'k) Senswrapper.t
                              -> ('d, 'k) Senswrapper.t
                              -> ('d, 'k) Senswrapper.t
                              -> float -> bool -> 's -> unit)
    -> set_sys_fn          : ((('d, 'k) Senswrapper.t, 's) sysfn -> unit)
    -> ?context:Context.t
    -> unit
    -> ('d, 'k, 's, [`Sens]) t

end (* }}} *)

(** {2:nlsexceptions Exceptions} *)

(** An error occurred in a vector operation.

    @nodoc SUN_NLS_VECTOROP_ERR *)
exception VectorOpError

(** Raised when a nonlinear solver is used incorrectly.
    For example, calling {!solve} without having first called {!set_sys_fn}
    ([SUN_NLS_MEM_NULL]). *)
exception IncorrectUse

(** Raised if an external library call fails. *)
exception ExtFail

(** Raised on an attempt to associate a nonlinear solver instance with more
    than one session. *)
exception NonlinearSolverInUse

