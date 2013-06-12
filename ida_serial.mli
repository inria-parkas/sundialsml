(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*           Timothy Bourke (INRIA) and Marc Pouzet (LIENS)            *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

(***********************************************************************)
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for IDA v2.7.0                     *)
(*         Alan C. Hindmarsh, Radu Serban, and Aaron Collier           *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

include module type of Ida
  with type Roots.t = Ida.Roots.t
  and type RootDirs.t = Ida.RootDirs.t
  and type linear_solver = Ida.linear_solver
  and type Bandmatrix.t = Dls.Bandmatrix.t
  and type Directbandmatrix.t = Dls.Directbandmatrix.t
  and type Densematrix.t = Dls.Densematrix.t
  and type Directdensematrix.t = Dls.Directdensematrix.t

(*STARTINTRO*)
(** Serial nvector interface to the IDA solver.

  Serial vectors are passed between Sundials and Ocaml programs as
  Bigarrays.
  These vectors are manipulated within the solver using the original low-level
  vector operations (cloning, linear sums, adding constants, and etcetera).
  While direct interfaces to these operations are not provided, there are
  equivalent implementations written in Ocaml for arrays of floats
  ({! Nvector_array}) and bigarrays ({! Nvector_array.Bigarray}) of floats.

  @version VERSION()
  @author Timothy Bourke (INRIA)
  @author Marc Pouzet (LIENS)
 *)

(**
    This type represents a session with the IDA solver using serial nvectors
    accessed as {{:OCAML_DOC_ROOT(Bigarray.Array1)} Bigarray.Array1}s.

    A skeleton of the main program:
    + {b Set vector of initial values}
    {[let y = Ida.Carray.of_array [| 0.0; 0.0; 0.0 |] ]}
    The length of this vector determines the problem size.
    + {b Create and initialize a solver session}
    {[let s = Ida.init Ida.Dense f (2, g) y]}
    This will initialize a specific linear solver and the root-finding
    mechanism, if necessary.
    + {b Specify integration tolerances (optional)}, e.g.
    {[ss_tolerances s reltol abstol]}
    + {b Set optional inputs}, e.g.
    {[set_stop_time s 10.0; ...]}
    Call any of the [set_*] functions to change solver parameters from their
    defaults.
    + {b Advance solution in time}, e.g.
    {[let (t', result) = Ida.solve_normal s !t y in
...
t := t' + 0.1]}
    Repeatedly call either [solve_normal] or [solve_one_step] to advance the
    simulation.
    + {b Get optional outputs}
    {[let stats = get_integrator_stats s in ...]}
    Call any of the [get_*] functions to examine solver statistics.

    @ida <node5#ss:skeleton_sim> Skeleton of main program
 *)
(*ENDINTRO*)
type session

(** The type of vectors passed to the solver. *)
type nvec = Sundials.Carray.t

(** The type of vectors containing dependent variable values, passed from the
   solver to callback functions. *)
type val_array = Sundials.Carray.t

(** The type of vectors containing derivative values, passed from the
   solver to callback functions. *)
type der_array = Sundials.Carray.t

(** The type of vectors containing detected roots (zero-crossings). *)
type root_array = Sundials.Roots.t

(** The type of vectors containing the values of root functions
   (zero-crossings). *)
type root_val_array = Sundials.Roots.val_array

(** {2 Initialization} *)

(**
    [init linsolv f (nroots, g) y0 y'0] initializes the IDA solver to solve
    the DAE f t y y' = 0 and returns a {!session}.
    - [linsolv] is the linear solver to attach to this solver.
    - [f]       is the residual function (see below).
    - [nroots]  specifies the number of root functions (zero-crossings).
    - [g]       calculates the values of the root functions.
    - [y0]      is a vector of initial values for the dependent-variable vector
                [y].  This vector's size determines the number of equations
                in the session, see {!Sundials.Carray.t}.
    - [y'0]     is a vector of initial values for [y'], i.e. the derivative
                of [y] with respect to t.  This vector's size must match the
                size of [y0].

    The start time defaults to 0. It can be set manually by instead using
    {!init_at_time}.

    This function calls IDACreate, IDAInit, IDARootInit, an appropriate
    linear solver function, and IDASStolerances (with default values for
    relative tolerance of 1.0e-4 and absolute tolerance as 1.0e-8; these can be
    changed with {!ss_tolerances}, {!sv_tolerances}, or {!wf_tolerances}).
    It does everything necessary to initialize an IDA session; the {!normal} or
    {!one_step} functions can be called directly afterward.

    The residual function [f] is called by the solver like [f t y y' r] to
    compute the problem residual, where:
    - [t] is the current value of the independent variable,
          i.e., the simulation time.
    - [y] is a vector of dependent-variable values, i.e. y(t).
    - [y'] is the derivative of [y] with respect to [t], i.e. dy/dt.
    - [r] is the output vector to fill in with the value of the residual
          function for the given values of t, y, and y'.
    The residual function should return normally if successful, raise
    RecoverableFailure if a recoverable error occurred (e.g. yy has an
    illegal value), or raise some other exception if a nonrecoverable error
    occurred.  If a recoverable error occurred, the integrator will attempt
    to correct and retry.  If a nonrecoverable error occurred, the integrator
    will halt and propagate the exception to the caller.

    {b NB:} [y], [y'], and [r] must no longer be accessed after [f] has
            returned a result, i.e. if their values are needed outside of
            the function call, then they must be copied to separate physical
            structures.

    The roots function [g] is called by the solver to calculate the values of
    root functions (zero-crossing expressions) which are used to detect
    significant events.  It is passed four arguments [t], [y], [y'], and
    [gout]:
    - [t], [y], [y'] are as for [f].
    - [gout] is a vector for storing the values of g(t, y, y').
    The {!Ida.no_roots} value can be passed for the [(nroots, g)] argument if
    root functions are not required.  If the root function raises an exception,
    the integrator will halt immediately and propagate the exception to the
    caller.

    {b NB:} [y] and [gout] must no longer be accessed after [g] has returned
            a result, i.e. if their values are needed outside of the function
            call, then they must be copied to separate physical structures.

    @ida <node5#sss:idamalloc>     IDACreate/IDAInit
    @ida <node5#ss:resFn>          ODE right-hand side function
    @ida <node5#ss:idarootinit>    IDARootInit
    @ida <node5#ss:rootFn>         Rootfinding function
    @ida <node5#sss:lin_solv_init> Linear solvers
    @ida <node5#sss:idatolerances> IDASStolerances
 *)
val init :
    linear_solver
    -> (float -> val_array -> der_array -> val_array -> unit)
    -> (int * (float -> val_array -> der_array -> root_val_array -> unit))
    -> nvec
    -> nvec
    -> session

(**
  [init_at_time linsolv roots t0 y0 y'0] is the same as init except that a
  start time [t0], can be given explicitly.
 *)
val init_at_time :
    linear_solver
    -> (float -> val_array -> der_array -> val_array -> unit)
    -> (int * (float -> val_array -> der_array -> root_val_array -> unit))
    -> float (* start time *)
    -> nvec
    -> nvec
    -> session

(** Return the number of root functions. *)
val nroots : session -> int

(** Return the number of equations. *)
val neqs : session -> int

(** {2 Tolerance specification} *)

(**
    [ss_tolerances s reltol abstol] sets the relative and absolute
    tolerances using scalar values.

    @ida <node5#sss:cvtolerances> IDASStolerances
 *)
val ss_tolerances : session -> float -> float -> unit

(**
    [sv_tolerances s reltol abstol] sets the relative tolerance using a scalar
    value, and the absolute tolerance as a vector.

    @ida <node5#sss:cvtolerances> IDASVtolerances
 *)
val sv_tolerances : session -> float -> nvec -> unit

(**
    [wf_tolerances s efun] specifies a function [efun] that sets the multiplicative
    error weights Wi for use in the weighted RMS norm.

    [efun y ewt] is passed the dependent variable vector [y] and is expected to
    set the values inside the error-weight vector [ewt].

    @ida <node5#sss:cvtolerances> IDAWFtolerances
    @ida <node5#ss:ewtsetFn> Error weight function
 *)
val wf_tolerances : session -> (val_array -> val_array -> unit) -> unit

(** {2 Solver functions } *)

(**
   [(tret, r) = normal s tout yout y'out] integrates the DAE over an interval
   in t.

   The arguments are:
   - [s] a session with the solver.
   - [tout] the next time at which a computed solution is desired.
   - [yout] a vector to store the computed solution. The same vector as was
   - [y'out] a vector to store the computed solution's derivative. The same
   vector as was passed to {!init} can (but does not have to) be used again
   for this argument.

   Two values are returned:
    - [tret] the time reached by the solver, which will be equal to [tout] if
      no errors occur.
    - [r] indicates whether roots were found, or whether an optional stop time,
   set by {!set_stop_time}, was reached; see {!Ida.solver_result}.

   This routine will throw one of the solver {!Ida.exceptions} if an error
   occurs.

   @ida <node5#sss:ida> IDASolve
   @ida <node5#sss:ida> IDA (IDA_NORMAL)
 *)
val solve_normal :
  session -> float -> val_array -> der_array -> float * solver_result

(**
   This function is identical to {!normal}, except that it returns after one
   internal solver step.

   @ida <node5#sss:ida> IDASolve
   @ida <node5#sss:ida> IDA (IDA_ONE_STEP)
 *)
val solve_one_step :
  session -> float -> val_array -> der_array -> float * solver_result

(** {2 Main optional functions} *)

(** {3 Input} *)

(**
  [set_error_file s fname trunc] opens the file named [fname] and to which all
  messages from the default error handler are then directed.
  If the file already exists it is either trunctated ([trunc] = [true]) or
  appended to ([trunc] = [false]).

  The error file is closed if set_error_file is called again, or otherwise when
  the session is garbage collected.
   
  @ida <node5#sss:optin_main> IDASetErrFile
 *)
val set_error_file : session -> string -> bool -> unit

(**
  [set_err_handler_fn s efun] specifies a custom function [efun] for handling
  error messages.  The error handler function must not fail -- any exceptions
  raised from it will be captured and silently discarded.

  @ida <node5#sss:optin_main> IDASetErrHandlerFn
  @ida <node5#ss:ehFn> Error message handler function
 *)
val set_err_handler_fn : session -> (error_details -> unit) -> unit

(**
  This function restores the default error handling function. It is equivalent
  to calling IDASetErrHandlerFn with an argument of [NULL].

  @ida <node5#sss:optin_main> IDASetErrHandlerFn
 *)
val clear_err_handler_fn : session -> unit

(**
  Specifies the maximum order of the linear multistep method.

  @ida <node5#sss:optin_main> IDASetMaxOrd
 *)
val set_max_ord : session -> int -> unit

(**
  Specifies the maximum number of steps to be taken by the solver in its attempt
  to reach the next output time.

  @ida <node5#sss:optin_main> IDASetMaxNumSteps
 *)
val set_max_num_steps : session -> int -> unit

(**
  Specifies the initial step size.

  @ida <node5#sss:optin_main> IDASetInitStep
 *)
val set_init_step : session -> float -> unit

(**
  Specifies an upper bound on the magnitude of the step size.

  @ida <node5#sss:optin_main> IDASetMaxStep
 *)
val set_max_step : session -> float -> unit

(**
  Specifies the value of the independent variable t past which the solution is
  not to proceed.
  The default, if this routine is not called, is that no stop time is imposed.

  @ida <node5#sss:optin_main> IDASetStopTime
 *)
val set_stop_time : session -> float -> unit

(**
  Specifies the maximum number of error test failures permitted in attempting
  one step.

  @ida <node5#sss:optin_main> IDASetMaxErrTestFails
 *)
val set_max_err_test_fails : session -> int -> unit

(**
  Specifies the maximum number of nonlinear solver iterations permitted per
  step.

  @ida <node5#sss:optin_main> IDASetMaxNonlinIters
 *)
val set_max_nonlin_iters : session -> int -> unit

(**
  Specifies the maximum number of nonlinear solver convergence failures
  permitted during one step.

  @ida <node5#sss:optin_main> IDASetMaxConvFails
 *)
val set_max_conv_fails : session -> int -> unit

(**
  Specifies the safety factor used in the nonlinear convergence test.

  @ida <node5#sss:optin_main> IDASetNonlinConvCoef
  @ida <node3#ss:ivp_sol> IVP Solution
 *)
val set_nonlin_conv_coef : session -> float -> unit

(** {3 Output } *)

(**
  Returns the real and integer workspace sizes.

  @ida <node5#sss:optout_main> IDAGetWorkSpace
  @return ([real_size], [integer_size])
 *)
val get_work_space          : session -> int * int

(**
  Returns the cumulative number of internal steps taken by the solver.

  @ida <node5#sss:optout_main> IDAGetNumSteps
 *)
val get_num_steps           : session -> int

(**
  Returns the number of calls to the user's right-hand side function.

  @ida <node5#sss:optout_main> IDAGetNumResEvals
 *)
val get_num_res_evals       : session -> int

(**
  Returns the number of calls made to the linear solver's setup function.

  @ida <node5#sss:optout_main> IDAGetNumLinSolvSetups
 *)
val get_num_lin_solv_setups : session -> int

(**
  Returns the number of local error test failures that have occurred.

  @ida <node5#sss:optout_main> IDAGetNumErrTestFails
 *)
val get_num_err_test_fails  : session -> int

(**
  Returns the integration method order used during the last internal step.

  @ida <node5#sss:optout_main> IDAGetLastOrder
 *)
val get_last_order          : session -> int

(**
  Returns the integration method order to be used on the next internal step.

  @ida <node5#sss:optout_main> IDAGetCurrentOrder
 *)
val get_current_order       : session -> int

(**
  Returns the integration step size taken on the last internal step.

  @ida <node5#sss:optout_main> IDAGetLastStep
 *)
val get_last_step           : session -> float

(**
  Returns the integration step size to be attempted on the next internal step.

  @ida <node5#sss:optout_main> IDAGetCurrentStep
 *)
val get_current_step        : session -> float

(**
  Returns the the value of the integration step size used on the first step.

  @ida <node5#sss:optout_main> IDAGetActualInitStep
 *)
val get_actual_init_step    : session -> float

(**
  Returns the the current internal time reached by the solver.

  @ida <node5#sss:optout_main> IDAGetCurrentTime
 *)
val get_current_time        : session -> float

(* IDAGetNumStabLimOrderReds appears in the sundials manual on p.52 but there's
   no such function in the implementation.  It's probably a typo or a leftover
   from earlier versions.

(**
   Returns the number of order reductions dictated by the BDF stability limit
   detection algorithm.

   @ida <node5#sss:optout_main> IDAGetNumStabLimOrderReds
   @ida <node3#s:bdf_stab> BDF stability limit detection
 *)
val get_num_stab_lim_order_reds : session -> int
 *)

(**
  Returns a suggested factor by which the user's tolerances should be scaled
  when too much accuracy has been requested for some internal step.

  @ida <node5#sss:optout_main> IDAGetTolScaleFactor
 *)
val get_tol_scale_factor : session -> float

(**
  Returns the solution error weights at the current time.

  @ida <node5#sss:optout_main> IDAGetErrWeights
  @ida <node3#ss:ivp_sol> IVP solution (W_i)
 *)
val get_err_weights : session -> nvec -> unit

(**
  Returns the vector of estimated local errors.

  @ida <node5#sss:optout_main> IDAGetEstLocalErrors
 *)
val get_est_local_errors : session -> nvec -> unit

(**
  Returns the integrator statistics as a group.

  @ida <node5#sss:optout_main> IDAGetIntegratorStats
 *)
val get_integrator_stats    : session -> Ida.integrator_stats

(**
  Convenience function that calls get_integrator_stats and prints the results to
  stdout.

  @ida <node5#sss:optout_main> IDAGetIntegratorStats
 *)
val print_integrator_stats  : session -> unit


(**
  Returns the number of nonlinear (functional or Newton) iterations performed.

  @ida <node5#sss:optout_main> IDAGetNumNonlinSolvIters
 *)
val get_num_nonlin_solv_iters : session -> int

(**
  Returns the number of nonlinear convergence failures that have occurred.

  @ida <node5#sss:optout_main> IDAGetNumNonlinSolvConvFails
 *)
val get_num_nonlin_solv_conv_fails : session -> int

(**
  Changes the linear solver.

  @ida <node5#sss:optout_main> IDADense
  @ida <node5#sss:optout_main> IDABand
  @ida <node5#sss:optout_main> IDASpgmr
  @ida <node5#sss:optout_main> IDASpbcg
  @ida <node5#sss:optout_main> IDASptfqmr
 *)
val set_linear_solver : session -> linear_solver -> unit

(** {2 Root finding optional functions} *)

(** {3 Input} *)

(**
  [set_root_direction s dir] specifies the direction of zero-crossings to be
  located and returned. [dir] may contain one entry of type
  {!Ida.root_direction} for each root function.

  @ida <node5#sss:optin_root> IDASetRootDirection
 *)
val set_root_direction : session -> root_direction array -> unit

(**
  Like {!set_root_direction} but specifies a single direction of type
  {!Ida.root_direction} for all root
  functions.

  @ida <node5#sss:optin_root> IDASetRootDirection
 *)
val set_all_root_directions : session -> root_direction -> unit

(**
  Disables issuing a warning if some root function appears to be identically
  zero at the beginning of the integration.

  @ida <node5#sss:optin_root> IDASetNoInactiveRootWarn
 *)
val set_no_inactive_root_warn : session -> unit

(** {3 Output} *)

(**
  Fills an array showing which functions were found to have a root.

  @ida <node5#sss:optout_root> IDAGetRootInfo
 *)
val get_root_info : session -> root_array -> unit

(**
  Returns the cumulative number of calls made to the user-supplied root function g.

  @ida <node5#sss:optout_root> IDAGetNumGEvals
 *)
val get_num_g_evals : session -> int

(** {2 Interpolated output function } *)

(**
  [get_dky s t k dky] computes the [k]th derivative of the function y at time
  [t], i.e. d(k)y/dt(k)(t). The function requires that tn - hu <= [t] <=
  tn, where tn denotes the current internal time reached, and hu is the last
  internal step size successfully used by the solver.
  The user may request [k] = 0, 1,..., qu, where qu is the current order.

  This function may only be called after a successful return from either
  {!normal} or {!one_step}.

  Values for the limits may be obtained:
    - tn = {!get_current_time}
    - qu = {!get_last_order}
    - hu = {!get_last_step}

  @ida <node5#sss:optin_root> IDAGetDky
 *)
val get_dky : session -> float -> int -> nvec -> unit

(** {2 Reinitialization} *)

(**
  [reinit s t0 y0 y'0] reinitializes the solver session [s] with a new time
  [t0] and new values for the variables [y0].  There are two optional arguments
  to change the linear solver and the set of root functions.

  [linsolv] sets the linear solver.  If omitted, the current linear solver will
  be kept.

  [roots] sets the root functions.  {!no_roots} may be passed in to turn off
  root finding.  If omitted, the current root functions will be kept.

  @ida <node5#sss:cvreinit> IDAReInit
 *)
val reinit :
  session
  -> ?linsolv:linear_solver
  -> ?roots:(int * (float -> val_array -> der_array -> root_val_array -> unit))
  -> float
  -> val_array
  -> der_array
  -> unit

(** {2 Linear Solvers} *)

type single_tmp = val_array
type double_tmp = val_array * val_array
type triple_tmp = val_array * val_array * val_array

(**
  Arguments common to all Jacobian callback functions.    
 
  @ida <node5#ss:djacFn> Dense Jacobian function
  @ida <node5#ss:bjacFn> Banded Jacobian function
  @ida <node5#ss:jtimesFn> Product Jacobian function
  @ida <node5#ss:psolveFn> Linear preconditioning function
  @ida <node5#ss:precondFn> Jacobian preconditioning function
 *)
type 't jacobian_arg =
  {
    jac_t    : float;        (** The independent variable. *)
    jac_y    : val_array;    (** The dependent variable vector. *)
    jac_y'   : der_array;    (** The derivative vector (i.e. dy/dt). *)
    jac_res  : val_array;    (** The current value of the residual vector. *)
    jac_coef : float;        (** The coefficient [a] in the system Jacobian
                                 to compute,
                                   [J = dF/dy + a*dF/d(y')]
                                 where [F(t,y,y')] is the residual vector.
                                 See Eq (2.5) of IDA's user documentation.  *)
    jac_tmp  : 't            (** Workspace data,
                                either {!single_tmp} or {!triple_tmp}. *)
  }

(** {3 Direct Linear Solvers (DLS)} *)

(** Control callbacks and get optional outputs for the Direct Linear Solvers
    that operate on dense and banded matrices.
    
    @ida <node5#sss:optin_dls> Direct linear solvers optional input functions
    @ida <node5#sss:optout_dls> Direct linear solvers optional output functions
    @ida <node5#ss:djacFn> Dense Jacobian function
  *)
module Dls :
  sig
    (** {4 Callback functions} *)
    (**
     Specify a callback function that computes an approximation to the Jacobian
     matrix J(t, y) for the Dense and Lapackdense {!Ida.linear_solver}s.

     The callback function takes the {!jacobian_arg} as an input and must store
     the computed Jacobian as a {!Ida.Densematrix.t}.  The Jacobian has the
     form

       [dF0/dy0 + c dF0/dy'0, dF0/dy1 + c dF0/dy'1, ..., dF0/dyn + c dF0/dy'n]
       [dF1/dy0 + c dF1/dy'0, dF1/dy1 + c dF1/dy'1, ..., dF1/dyn + c dF1/dy'n]
           :        :        :             :
       [dFn/dy0 + c dFn/dy'0, dFn/dy1 + c dFn/dy'1, ..., dFn/dyn + c dFn/dy'n]

     i.e. each row should be a gradient, or put differently, the row index
     matches the equation index while the column index matches the variable
     index.

     {b NB:} the elements of the Jacobian argument and the output matrix must no
     longer be accessed after callback function has returned a result, i.e. if
     their values are needed outside of the function call, then they must be
     copied to separate physical structures.

     @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
     @ida <node5#ss:djacFn> Dense Jacobian function
     *)
    val set_dense_jac_fn :
         session
      -> (triple_tmp jacobian_arg -> Densematrix.t -> unit)
      -> unit

    (**
      This function disables the user-supplied dense Jacobian function, and
      switches back to the default internal difference quotient approximation
      that comes with the Dense and Lapackdense {!Ida.linear_solver}s. It is
      equivalent to calling IDASetDenseJacFn with an argument of [NULL].

      @ida <node5#ss:djacFn> Dense Jacobian function
    *)
    val clear_dense_jac_fn : session -> unit

    (**
     Specify a callback function that computes an approximation to the Jacobian
     matrix J(t, y) for the Band and Lapackband {!Ida.linear_solver}s.

     The callback function takes three input arguments:
     - [jac] the standard {!jacobian_arg} with three work vectors.
     - [mupper] the upper half-bandwidth of the Jacobian.
     - [mlower] the lower half-bandwidth of the Jacobian.
     and it must store the computed Jacobian as a {!Ida.Bandmatrix.t}.

    {b NB:} [jac] and the computed Jacobian must no longer be accessed after the
            calback function has returned a result, i.e. if their values are
            needed outside of the function call, then they must be copied to
            separate physical structures.

     @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
     @ida <node5#ss:bjacFn> Banded Jacobian function
     *)
    val set_band_jac_fn :
         session
      -> (triple_tmp jacobian_arg -> int -> int -> Bandmatrix.t -> unit)
      -> unit

    (**
      This function disables the user-supplied band Jacobian function, and
      switches back to the default internal difference quotient approximation
      that comes with the Band and Lapackband {!Ida.linear_solver}s. It is
      equivalent to calling IDASetBandJacFn with an argument of [NULL].

      @ida <node5#ss:bjacFn> Banded Jacobian function
    *)
    val clear_band_jac_fn : session -> unit

    (** {4 Optional input functions} *)

    (**
      Returns the sizes of the real and integer workspaces used by the Dense and
      Band direct linear solvers .

      @ida <node5#sss:optout_dls> IDADlsGetWorkSpace
      @return ([real_size], [integer_size])
     *)
    val get_work_space : session -> int * int


    (**
      Returns the number of calls made to the Dense and Band direct linear
      solvers Jacobian approximation function.

      @ida <node5#sss:optout_dls> IDADlsGetNumJacEvals
    *)
    val get_num_jac_evals : session -> int

    (**
      Returns the number of calls made to the user-supplied right-hand side
      function due to the finite difference (Dense or Band) Jacobian
      approximation.

      @ida <node5#sss:optout_dls> IDADlsGetNumResEvals
    *)
    val get_num_res_evals : session -> int
  end

(** {3 Scaled Preconditioned Iterative Linear Solvers (SPILS)} *)

(** Set callback functions, set optional outputs, and get optional inputs for
    the Scaled Preconditioned Iterative Linear Solvers: SPGMR, SPBCG, SPTFQMR.
    @ida <node5#sss:optin_spils> Iterative linear solvers optional input functions.
    @ida <node5#sss:optout_spils> Iterative linear solvers optional output functions.
    @ida <node5#ss:psolveFn> Linear preconditioning function
    @ida <node5#ss:precondFn> Jacobian preconditioning function
 *)
module Spils :
  sig
    (** {4 Callback functions} *)

    (**
      Setup preconditioning for any of the SPILS linear solvers. Two functions
      are required: [psetup] and [psolve].

      [psetup jac] preprocesses and/or evaluates any Jacobian-related
       data needed by the preconditioner.  There is one argument:
        - [jac] supplies the basic problem data as a {!jacobian_arg}.

       Note that unlike in CVODE, whatever data this function computes has to
       be recomputed every time it is called.

       It should raise RecoverableError to instruct the solver to retry with a
       different step size, or simply raise any other exception to abort the
       solver.  In the latter case, the exception will be propagated out of the
       solver.

      {b NB:} The elements of [jac] must no longer be accessed after [psetup]
              has returned a result, i.e. if their values are needed outside
              of the function call, then they must be copied to a separate
              physical structure.

      [psolve jac r z delta] is called to solve the linear system
      {i P}[z] = [r], where {i P} is the (left) preconditioner matrix.
      {i P} should approximate, at least crudely, the system Jacobian matrix
      J = dF/dy + {jac.coef} * dF/d(y') where F is the residual function.
      - [jac] supplies the basic problem data as a {!jacobian_arg}.
      - [r] is the right-hand side vector.
      - [z] is the vector in which the result must be stored.
      - [delta] is an input tolerance.

      [delta] is be used if an iterative method is employed in the solution.
      In that than case, the residual vector Res = r - {i P} {i z} of the
      system should be made less than [delta] in weighted l2 norm, i.e.
        sqrt (sum_i ((Res_i * ewt_i)^2)) < [delta],
      where the nvector ewt can be obtained through
      {!Ida_serial.get_err_weights}.

      {b NB:} The elements of [jac], [r], and [z] must no longer be accessed
              after [psolve] has returned a result, i.e. if their values are
              needed outside of the function call, then they must be copied
              to separate physical structures.

      @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
      @ida <node5#ss:psolveFn> Linear preconditioning function
      @ida <node5#ss:precondFn> Jacobian preconditioning function
    *)
    val set_preconditioner :
      session
      -> (triple_tmp jacobian_arg -> unit)
      -> (single_tmp jacobian_arg -> val_array -> val_array -> float -> unit)
      -> unit

    (**
      Specifies a Jacobian-times-vector function.

      The function given, [jactimes jac v Jv], computes the matrix-vector
      product {i J}[v].
      - [jac] provides the data necessary to compute the Jacobian.
      - [v] is the vector by which the Jacobian must be multiplied.
      - [Jv] is the vector in which the result must be stored.

      {b NB:} The elements of [jac], [v], and [Jv] must no longer be accessed
              after [psolve] has returned a result, i.e. if their values are
              needed outside of the function call, then they must be copied
              to separate physical structures.

      @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
      @ida <node5#ss:jtimesFn> Product Jacobian function
    *)
    val set_jac_times_vec_fn :
      session
      -> (double_tmp jacobian_arg
          -> val_array (* v *)
          -> val_array (* Jv *)
          -> unit)
      -> unit

    (**
      This function disables the user-supplied Jacobian-vector function, and
      switches back to the default internal difference quotient approximation.
      It is equivalent to calling IDASpilsSetJacTimesVecFn with an argument of
      [NULL].

      @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
      @ida <node5#ss:jtimesFn> Product Jacobian function
    *)
    val clear_jac_times_vec_fn : session -> unit

    (** {4 Optional output functions} *)

    (** Constants representing the types of Gram-Schmidt orthogonalization
        possible for the Spgmr {Ida.linear_solver}. *)
    type gramschmidt_type =
      | ModifiedGS
            (** Modified Gram-Schmidt orthogonalization (MODIFIED_GS) *)
      | ClassicalGS
            (** Classical Gram Schmidt orthogonalization (CLASSICAL_GS) *)

    (**
      Sets the Gram-Schmidt orthogonalization to be used with the
      Spgmr {!Ida.linear_solver}.

      @ida <node5#sss:optin_spils> IDASpilsSetGSType
    *)
    val set_gs_type : session -> gramschmidt_type -> unit

    (**
      [set_eps_lin eplifac] sets the factor by which the Krylov linear solver's
      convergence test constant is reduced from the Newton iteration test
      constant. [eplifac]  must be >= 0. Passing a value of 0 specifies the
      default (which is 0.05).

      @ida <node5#sss:optin_spils> IDASpilsSetEpsLin
    *)
    val set_eps_lin : session -> float -> unit

    (**
      [set_maxl maxl] resets the maximum Krylov subspace dimension for the
      Bi-CGStab or TFQMR methods. [maxl] is the maximum dimension of the Krylov
      subspace, a value of [maxl] <= 0 specifies the default (which is 5.0).

      @ida <node5#sss:optin_spils> IDASpilsSetMaxl
    *)
    val set_maxl : session -> int -> unit

    (** {4 Optional input functions} *)

    (**
      Returns the sizes of the real and integer workspaces used by the SPGMR
      linear solver.

      @ida <node5#sss:optout_spils> IDASpilsGetWorkSpace
      @return ([real_size], [integer_size])
    *)
    val get_work_space       : session -> int * int

    (**
      Returns the cumulative number of linear iterations.

      @ida <node5#sss:optout_spils> IDASpilsGetNumLinIters
    *)
    val get_num_lin_iters    : session -> int

    (**
      Returns the cumulative number of linear convergence failures.

      @ida <node5#sss:optout_spils> IDASpilsGetNumConvFails
    *)
    val get_num_conv_fails   : session -> int

    (**
      Returns the number of preconditioner evaluations, i.e., the number of
      calls made to psetup with jok = [false] (see {!set_preconditioner}).

      @ida <node5#sss:optout_spils> IDASpilsGetNumPrecEvals
    *)
    val get_num_prec_evals   : session -> int

    (**
      Returns the cumulative number of calls made to the preconditioner solve
      function, psolve (see {!set_preconditioner}).

      @ida <node5#sss:optout_spils> IDASpilsGetNumPrecSolves
    *)
    val get_num_prec_solves  : session -> int

    (**
      Returns the cumulative number of calls made to the Jacobian-vector
      function, jtimes (see {! set_jac_times_vec_fn}).

      @ida <node5#sss:optout_spils> IDASpilsGetNumJtimesEvals
    *)
    val get_num_jtimes_evals : session -> int

    (**
      Returns the number of calls to the user right-hand side function for
      finite difference Jacobian-vector product approximation. This counter is
      only updated if the default difference quotient function is used.

      @ida <node5#sss:optout_spils> IDASpilsGetNumResEvals
    *)
    val get_num_res_evals    : session -> int
  end

(** Inequality constraints on variables.

 @ida <node5#sss:idasetconstraints> IDASetConstraints
 *)
module Constraints :
  sig
    (** An abstract array type, whose i-th component specifies that the i-th
        component of the dependent variable vector y should be:

        NonNegative  i.e. >= 0, or
        NonPositive  i.e. <= 0, or
        Positive     i.e. > 0, or
        Negative     i.e. < 0, or
        Unconstrained
     *)
    type t
    type constraint_type =
    | Unconstrained
    | NonNegative
    | NonPositive
    | Positive
    | Negative

    (** [create n] returns an array with [n] elements, each set to
        Unconstrained.  *)
    val create : int -> t

    (** [init n x] returns an array with [n] elements, each set to [x]. *)
    val init : int -> constraint_type -> t

    (** [of_array a] converts an OCaml array [a] of {!constraint_type}s into an
        abstract array suitable for passing into IDA.  *)
    val of_array : constraint_type array -> t

    (** Returns the length of an array *)
    val length : t -> int

    (** [get c i] returns the constraint on the i-th variable in the DAE.  *)
    val get : t -> int -> constraint_type

    (** [set c i x] sets the constraint on the i-th variable in the DAE to
        [x].  *)
    val set : t -> int -> constraint_type -> unit

    (** [fill c x] fills the array so that all variables will have constraint
        [x].  *)
    val fill : t -> constraint_type -> unit

    (** [blit a b] copies the contents of [a] to [b].  *)
    val blit : t -> t -> unit
  end

(** Variable classification that needs to be specified for computing consistent
 initial values and for suppressing local error tests on some variables.

 @ida <node5#sss:idasetid> IDASetId
 *)
module VarTypes :
  sig
    (** An abstract array type, whose i-th component specifies whether the i-th
        component of the dependent variable vector y is an algebraic or
        differential variable, for each i.  *)
    type t
    type var_type =
    | Algebraic    (** Algebraic variable; residual function must not depend
                       on this component's derivative.  *)
    | Differential (** Differential variable; residual function can depend on
                       this component's derivative.  *)

    (** [create n] returns an array with [n] elements, each set to
        Algebraic.  *)
    val create : int -> t

    (** [init n x] returns an array with [n] elements, each set to [x]. *)
    val init : int -> var_type -> t

    (** [of_array a] converts an OCaml array [a] of {!var_type}s into an
        abstract array suitable for passing into IDA.  *)
    val of_array : var_type array -> t

    (** Returns the length of an array *)
    val length : t -> int

    (** [get c i] returns the component type of the i-th variable in the
        DAE.  *)
    val get : t -> int -> var_type

    (** [set c i x] sets the component type of the i-th variable in the DAE to
        [x].  *)
    val set : t -> int -> var_type -> unit

    (** [set_algebraic c i] sets the component type of the i-th variable in
        the DAE to algebraic.  *)
    val set_algebraic : t -> int -> unit

    (** [set_differential c i] sets the component type of the i-th variable
        in the DAE to differential.  *)
    val set_differential : t -> int -> unit

    (** [fill c x] fills the array so that all variables will have component
        type [x].  *)
    val fill : t -> var_type -> unit

    (** [blit a b] copies the contents of [a] to [b].  *)
    val blit : t -> t -> unit
  end

(** An unpreferred alias for {!VarTypes}.  SUNDIALS calls variable types by the
    cryptic name "Id", and this OCaml binding preserves this alternative naming
    to help users transition from the C interface.  *)
module Id :
  sig
    (** An abstract array type, whose i-th component specifies whether the i-th
        component of the dependent variable vector y is an algebraic or
        differential variable, for each i.  *)
    type t = VarTypes.t
    type var_type = VarTypes.var_type = Algebraic | Differential

    (** [create n] returns an array with [n] elements, each set to
        Algebraic.  *)
    val create : int -> t

    (** [init n x] returns an array with [n] elements, each set to [x]. *)
    val init : int -> var_type -> t

    (** [of_array a] converts an OCaml array [a] of {!var_type}s into an
        abstract array suitable for passing into IDA.  *)
    val of_array : var_type array -> t

    (** Returns the length of an array *)
    val length : t -> int

    (** [get c i] returns the component type of the i-th variable in the
        DAE.  *)
    val get : t -> int -> var_type

    (** [set c i x] sets the component type of the i-th variable in the DAE to
        [x].  *)
    val set : t -> int -> var_type -> unit

    (** [set_algebraic c i] sets the component type of the i-th variable in
        the DAE to algebraic.  *)
    val set_algebraic : t -> int -> unit

    (** [set_differential c i] sets the component type of the i-th variable
        in the DAE to differential.  *)
    val set_differential : t -> int -> unit

    (** [fill c x] fills the array so that all variables will have component
        type [x].  *)
    val fill : t -> var_type -> unit

    (** [blit a b] copies the contents of [a] to [b].  *)
    val blit : t -> t -> unit
  end

val set_constraints : session -> Constraints.t -> unit

(**
   Specify whether or not to ignore algebraic variables in local error tests.
   If you set this option to [true], you must also specify which variables are
   algebraic through {!calc_ic_ya_yd'} or {!set_var_types} -- exactly one of
   these functions should be called, exactly once, before the first call to
   {!solve_normal}, {!solve_one_step}, or {!calc_ic}.  Forgetting to do so will
   cause an {!Ida.IllInput} exception.

   Note: {!set_var_types} is the preferred alias to {!set_id}, which
   corresponds to [IDASetId] in the C interface.

   In general, suppressing local error tests for algebraic variables is {i
   discouraged} when solving DAE systems of index 1, whereas it is generally {i
   encouraged} for systems of index 2 or more.  See pp. 146-147 of the
   following reference for more on this issue:

   K. E. Brenan, S. L. Campbell, and L. R. Petzold.  Numerical Solution of
   Initial-Value Problems in Differential-Algebraic Equations.  SIAM,
   Philadelphia, Pa, 1996.

   @ida <node#sss:idacalcic> IDASetSuppressAlg
 *)
val set_suppress_alg : session -> bool -> unit

(** An unpreferred alias for {!set_var_types}.  SUNDIALS calls variable types
    by the cryptic name "Id", and this OCaml binding preserves this alternative
    naming to help users transition from the C interface.  *)
val set_id : session -> Id.t -> unit

(** Specify which variables are algebraic and which variables are differential,
    needed for {!set_suppress_alg}.  The function {!calc_ic_ya_yd'} also sets
    the same information, so you only need to call one of these functions.
    FIXME: it's not clear from the SUNDIALS manual if it's supposed to be safe
    to change the variable types after you've already set it.

    [set_var_types] corresponds to [IDASetId] in the C interface, and an alias
    {!set_id} is also available in this binding.  We prefer the more
    descriptive name {!set_var_types}, however.

    @ida <node#sss:idasetid> IDASetId
 *)
val set_var_types : session -> VarTypes.t -> unit

(** [calc_ic_y ida ~y:yvar tout1] corrects the initial values y0 at time t0,
    using the initial values of the derivatives y'0 that are stored in the
    session.  That is, if the t0,y0,y'0 that were given to {!init_at_time},
    {!init}, or {!reinit} does not satisfy F(t0,y0,y'0) = 0, where F is the
    residual function, then [calc_ic_y] will modify y'0 so that this equation
    holds.  If F(t0,y0,y'0) = 0 is already true, a call to [calc_ic_y] is
    unnecessary.

    The optional parameter [~y], if given, will receive the corrected y vector.
    [tout1] is the first value of t at which a solution will be requested
    (using {!solve_normal} or {!solve_one_step}). This value is needed here
    only to determine the direction of integration and rough scale in the
    independent variable t.

    IDA's initial value correction works for certain index-one problems
    including a class of systems of semi-implicit form, and uses Newton
    iteration combined with a linesearch algorithm.  See Section 2.1 of the IDA
    User Guide and the following reference for more information:

    P. N. Brown, A. C. Hindmarsh, and L. R. Petzold. Consistent Initial Condition Calculation for Differential-Algebraic Systems. SIAM J. Sci. Comput., 19:1495–1512, 1998.

    @ida <node#sss:idacalcic> IDACalcIC
    @ida <node#sss:idagetconsistentic> IDAGetConsistentIC
 *)
val calc_ic_y : session -> ?y:val_array -> float -> unit

(** [calc_ic_ya_yd' ida ~y:yvar ~y':y'var vartypes tout1] corrects the initial
    values y0 and y0' at time t0.  The optional parameters [~y] and [~y'], if
    given, will be overwritten by the corrected vectors.  [tout1] is the first
    value of t at which a solution will be requested, and is needed here only
    to determine the direction of integration and rough scale in the
    independent variable t.  See {!calc_ic_y} for more general information
    about initial value correction.

    [calc_ic_ya_yd'] differs from {!calc_ic_y} in that, [calc_ic_ya_yd'] can
    compute the initial values of some derivatives whereas {!calc_ic_y} can
    only compute the initial values of non-derivatives.

    The [vartypes] argument specifies which components of y are algebraic --
    meaning their derivatives do not appear in the DAE -- and which components
    are differential -- meaning their derivatives appear in the DAE.  Note that
    y here means the vector formed by collecting scalar variables that appear
    in the mathematical description of your DAE system; it does not mean the
    array passed in as the labeled argument whose name also happens to be [y].
    Likewise for [y'].

    [calc_ic_ya_yd'] modifies the algebraic components of y and differential
    components of y', using the differential components of y as input.  So if
    Ia is the set of indices at which [vartypes] is [Algebraic] and Id is the
    set of indices at which [vartypes] is [Differential], then y and y' are
    partitioned into two subsequences each:

      y  splits into A  = { y.{i}  | i in Ia } and D  = { y.{i}  | i in Id }
      y' splits into A' = { y'.{i} | i in Ia } and D' = { y'.{i} | i in Id }

    The residual function must be such that it ignores all values in A'.
    [calc_ic_ya_yd'] then computes (i.e. modifies) A and D' while treating D as
    read-only and ignoring A'.

      input:   D
      output:  A, D'
      ignored: A'

    Note: [vartypes] is called "id" in the C interface, e.g. [IDASetId].

    [calc_ic_ya_yd'] sets the variable types that {!set_suppress_alg} uses, so
    you do not need to set it again with {!set_var_types} (or its alias
    {!set_id}) before calling {!set_suppress_alg}.

    @ida <node#sss:idacalcic> IDACalcIC
    @ida <node#sss:idasetid> IDASetId
    @ida <node#sss:idasetsuppressalg> IDASetSuppressAlg
    @ida <node#sss:idagetconsistentic> IDAGetConsistentIC
 *)
val calc_ic_ya_yd' :
  session
  -> ?y:val_array
  -> ?y':der_array
  -> VarTypes.t
  -> float
  -> unit

(** [get_num_backtrack_ops ida] gets the number of backtrack operations done in
    the linesearch algorithm in {!calc_ic_ya_yd'} or {!calc_ic_y}.
    @ida <ndoe#sss:idagetnumbcktrackops> IDAGetNumBcktrackOps
 *)
val get_num_backtrack_ops : session -> int
