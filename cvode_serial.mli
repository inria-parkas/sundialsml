(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*       Timothy Bourke (INRIA Rennes) and Marc Pouzet (LIENS)         *)
(*                                                                     *)
(*  Copyright 2011 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

(***********************************************************************)
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(*STARTINTRO*)
(** Serial nvector interface to the CVODE solver.

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 *)
(*ENDINTRO*)

include module type of Cvode

(**
    This type represents a session with the CVODE solver using serial nvectors
    accessed as {{:OCAML_DOC_ROOT(Bigarray.Array1)} Bigarray.Array1}s.

    A skeleton of the main program:
    + {b Set vector of initial values}
    {[let y = Cvode.Carray.of_array [| 0.0; 0.0; 0.0 |] ]}
    The length of this vector determines the problem size.    
    + {b Create and initialize a solver session}
    {[let s = Cvode.init Cvode.Adams Cvode.Functional f (2, g) y]}
    This will initialize a specific linear solver and the root-finding
    mechanism, if necessary.
    + {b Specify integration tolerances (optional)}, e.g.
    {[ss_tolerances s reltol abstol]}
    + {b Set optional inputs}, e.g.
    {[set_stop_time s 10.0; ...]}
    Call any of the [set_*] functions to change solver parameters from their
    defaults.
    + {b Advance solution in time}, e.g.
    {[let (t', result) = Cvode.normal s !t y in
...
t := t' + 0.1]}
    Repeatedly call either [normal] or [one_step] to advance the simulation.
    + {b Get optional outputs}
    {[let stats = get_integrator_stats s in ...]}
    Call any of the [get_*] functions to examine solver statistics.

    @cvode <node5#ss:skeleton_sim> Skeleton of main program
 *)
type session

(** The type of vectors passed to the solver, see {!Sundials.Carray.t}. *)
type nvec = Carray.t

(** The type of vectors containing dependent variable values, passed from the
   solver to callback functions, see {!Sundials.Carray.t} *)
type val_array = Carray.t

(** The type of vectors containing derivative values, passed from the
   solver to callback functions, see {!Sundials.Carray.t} *)
type der_array = Carray.t

(** The type of vectors containing detected roots (zero-crossings), see
   {!Sundials.Roots.t} *)
type root_array = Roots.t

(** The type of vectors containing the values of root functions
   (zero-crossings), see {!Sundials.Roots.val_array} *)
type root_val_array = Roots.val_array

(** {2 Initialization} *)

(**
    [init lmm iter f (nroots, g) y0] initializes the CVODE solver and returns a
    {!session}.
    - [lmm]     specifies the linear multistep method, see {!Cvode.lmm}.
    - [iter]    specifies either functional iteration or Newton iteration
                with a specific linear solver, see {!Cvode.iter}.
    - [f]       is the ODE right-hand side function.
    - [nroots]  specifies the number of root functions (zero-crossings).
    - [g]       calculates the values of the root functions.
    - [y0]      is a vector of initial values, the size of this vector
                determines the number of equations in  the session, see
                {!Sundials.Carray.t}.

    The start time defaults to 0. It can be set manually by instead using
    {!init'}.

    This function calls CVodeCreate, CVodeInit, CVodeRootInit, an appropriate
    linear solver function, and CVodeSStolerances (with default values for
    relative tolerance of 1.0e-4 and absolute tolerance as 1.0e-8; these can be
    changed with {!ss_tolerances}, {!sv_tolerances}, or {!wf_tolerances}).
    It does everything necessary to initialize a CVODE session; the {!normal} or
    {!one_step} functions can be called directly afterward.

    The right-hand side function [f] is called by the solver to calculate the
    instantaneous derivative values, it is passed three arguments: [t], [y], and
    [dy].
    - [t] is the current value of the independent variable,
          i.e., the simulation time.
    - [y] is a vector of dependent-variable values, i.e. y(t).
    - [dy] is a vector for storing the value of f(t, y).

    The roots function [g] is called by the solver to calculate the values of
    root functions (zero-crossing expressions) which are used to detect
    significant events, it is passed three arguments: [t], [y], and [gout].
    - [t] and [y] are as for [f].
    - [gout] is a vector for storing the values of g(t, y).
    The {!Cvode.no_roots} value can be passed for the [(nroots, g)] argument if
    root functions are not required.

    @cvode <node5#sss:cvodemalloc>   CVodeCreate/CVodeInit
    @cvode <node5#ss:rhsFn>          ODE right-hand side function
    @cvode <node5#ss:cvrootinit>     CVodeRootInit
    @cvode <node5#ss:rootFn>         Rootfinding function
    @cvode <node5#sss:lin_solv_init> Linear solvers
    @cvode <node5#sss:cvtolerances> CVodeSStolerances
 *)
val init :
    lmm
    -> iter
    -> (float -> val_array -> der_array -> unit)
    -> (int * (float -> val_array -> root_val_array -> unit))
    -> nvec
    -> session

(**
  [init lmm iter f roots y0 t0] is the same as init' except that the start time,
  [t0], must be given explicitly.
 *)
val init' :
    lmm
    -> iter
    -> (float -> val_array -> der_array -> unit)
    -> (int * (float -> val_array -> root_val_array -> unit))
    -> nvec
    -> float (* start time *)
    -> session

(** Return the number of root functions. *)
val nroots : session -> int

(** Return the number of equations. *)
val neqs : session -> int

(** {2 Tolerance specification} *)

(**
    [ss_tolerances s reltol abstol] sets the relative and absolute
    tolerances using scalar values.

    @cvode <node5#sss:cvtolerances> CVodeSStolerances
 *)
val ss_tolerances : session -> float -> float -> unit

(**
    [sv_tolerances s reltol abstol] sets the relative tolerance using a scalar
    value, and the absolute tolerance as a vector.

    @cvode <node5#sss:cvtolerances> CVodeSVtolerances
 *)
val sv_tolerances : session -> float -> nvec -> unit

(**
    [wf_tolerances s efun] specifies a function [efun] that sets the multiplicative
    error weights Wi for use in the weighted RMS norm.

    [efun y ewt] is passed the dependent variable vector [y] and is expected to
    set the values inside the error-weight vector [ewt].

    @cvode <node5#sss:cvtolerances> CVodeWFtolerances
    @cvode <node5#ss:ewtsetFn> Error weight function
 *)
val wf_tolerances : session -> (val_array -> val_array -> unit) -> unit

(** {2 Solver functions } *)

(**
    {b TODO}: write this description.

    @cvode <node5#sss:cvode> CVode (CV_NORMAL)
 *)
val normal : session -> float -> nvec -> float * solver_result

(**
    {b TODO}: write this description.

    @cvode <node5#sss:cvode> CVode (CV_ONE_STEP)
 *)
val one_step : session -> float -> nvec -> float * solver_result

(** {2 Main optional functions} *)

(** {3 Input} *)

(**
  [set_error_file s fname trunc] opens the file named [fname] and to which all
  messages from the default error handler are then directed.
  If the file already exists it is either trunctated ([trunc] = [true]) or
  appended to ([trunc] = [false]).

  The error file is closed if set_error_file is called again, or otherwise when
  the session is garbage collected.
   
  @cvode <node5#sss:optin_main> CVodeSetErrFile
 *)
val set_error_file : session -> string -> bool -> unit

(**
  [set_err_handler_fn s efun] specifies a custom function [efun] for handling
  error messages.

  @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
  @cvode <node5#ss:ehFn> Error message handler function
 *)
val set_err_handler_fn : session -> (error_details -> unit) -> unit

(**
  This function restores the default error handling function. It is equivalent
  to calling CVodeSetErrHandlerFn with an argument of [NULL].

  @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
 *)
val clear_err_handler_fn : session -> unit

(**
  Specifies the maximum order of the linear multistep method.

  @cvode <node5#sss:optin_main> CVodeSetMaxOrd
 *)
val set_max_ord : session -> int -> unit

(**
  Specifies the maximum number of steps to be taken by the solver in its attempt
  to reach the next output time.

  @cvode <node5#sss:optin_main> CVodeSetMaxNumSteps
 *)
val set_max_num_steps : session -> int -> unit

(**
  Specifies the maximum number of messages issued by the solver warning that t +
  h = t on the next internal step.

  @cvode <node5#sss:optin_main> CVodeSetMaxHnilWarns
 *)
val set_max_hnil_warns : session -> int -> unit

(**
  Indicates whether the BDF stability limit detection algorithm should be used.

  @cvode <node5#sss:optin_main> CVodeSetStabLimDet
  @cvode <node3#s:bdf_stab> BDF Stability Limit Detection
 *)
val set_stab_lim_det : session -> bool -> unit

(**
  Specifies the initial step size.

  @cvode <node5#sss:optin_main> CVodeSetInitStep
 *)
val set_init_step : session -> float -> unit

(**
  Specifies a lower bound on the magnitude of the step size.

  @cvode <node5#sss:optin_main> CVodeSetMinStep
 *)
val set_min_step : session -> float -> unit

(**
  Specifies an upper bound on the magnitude of the step size.

  @cvode <node5#sss:optin_main> CVodeSetMaxStep
 *)
val set_max_step : session -> float -> unit

(**
  Specifies the value of the independent variable t past which the solution is
  not to proceed.
  The default, if this routine is not called, is that no stop time is imposed.

  @cvode <node5#sss:optin_main> CVodeSetStopTime
 *)
val set_stop_time : session -> float -> unit

(**
  Specifies the maximum number of error test failures permitted in attempting
  one step.

  @cvode <node5#sss:optin_main> CVodeSetMaxErrTestFails
 *)
val set_max_err_test_fails : session -> int -> unit

(**
  Specifies the maximum number of nonlinear solver iterations permitted per
  step.

  @cvode <node5#sss:optin_main> CVodeSetMaxNonlinIters
 *)
val set_max_nonlin_iters : session -> int -> unit

(**
  Specifies the maximum number of nonlinear solver convergence failures
  permitted during one step.

  @cvode <node5#sss:optin_main> CVodeSetMaxConvFails
 *)
val set_max_conv_fails : session -> int -> unit

(**
  Specifies the safety factor used in the nonlinear convergence test.

  @cvode <node5#sss:optin_main> CVodeSetNonlinConvCoef
  @cvode <node3#ss:ivp_sol> IVP Solution
 *)
val set_nonlin_conv_coef : session -> float -> unit

(**
  [set_iter_type s iter] resets the nonlinear solver iteration type to [iter]
  ({!Cvode.iter}).

  {b TODO}: describe what happens internally.

  @cvode <node5#sss:optin_main> CVodeSetIterType
 *)
val set_iter_type : session -> iter -> unit

(** {3 Output } *)

(**
  Returns the real and integer workspace sizes.

  @cvode <node5#sss:optout_main> CVodeGetWorkSpace
  @return ([lenrw], [leniw])
 *)
val get_work_space          : session -> int * int

(**
  Returns the cumulative number of internal steps taken by the solver.

  @cvode <node5#sss:optout_main> CVodeGetNumSteps
 *)
val get_num_steps           : session -> int

(**
  Returns the number of calls to the user's right-hand side function.

  @cvode <node5#sss:optout_main> CVodeGetNumRhsEvals
 *)
val get_num_rhs_evals       : session -> int

(**
  Returns the number of calls made to the linear solver's setup function.

  @cvode <node5#sss:optout_main> CVodeGetNumLinSolvSetups
 *)
val get_num_lin_solv_setups : session -> int

(**
  Returns the number of local error test failures that have occurred.

  @cvode <node5#sss:optout_main> CVodeGetNumErrTestFails
 *)
val get_num_err_test_fails  : session -> int

(**
  Returns the integration method order used during the last internal step.

  @cvode <node5#sss:optout_main> CVodeGetLastOrder
 *)
val get_last_order          : session -> int

(**
  Returns the integration method order to be used on the next internal step.

  @cvode <node5#sss:optout_main> CVodeGetCurrentOrder
 *)
val get_current_order       : session -> int

(**
  Returns the integration step size taken on the last internal step.

  @cvode <node5#sss:optout_main> CVodeGetLastStep
 *)
val get_last_step           : session -> float

(**
  Returns the integration step size to be attempted on the next internal step.

  @cvode <node5#sss:optout_main> CVodeGetCurrentStep
 *)
val get_current_step        : session -> float

(**
  Returns the the value of the integration step size used on the first step.

  @cvode <node5#sss:optout_main> CVodeGetActualInitStep
 *)
val get_actual_init_step    : session -> float

(**
  Returns the the current internal time reached by the solver.

  @cvode <node5#sss:optout_main> CVodeGetCurrentTime
 *)
val get_current_time        : session -> float

(**
  Returns the number of order reductions dictated by the BDF stability limit
  detection algorithm.

  @cvode <node5#sss:optout_main> CVodeGetNumStabLimOrderReds
  @cvode <node3#s:bdf_stab> BDF stability limit detection
 *)
val get_num_stab_lim_order_reds : session -> int

(**
  Returns a suggested factor by which the user's tolerances should be scaled
  when too much accuracy has been requested for some internal step.

  @cvode <node5#sss:optout_main> CVodeGetTolScaleFactor
 *)
val get_tol_scale_factor : session -> float

(**
  Returns the solution error weights at the current time.

  @cvode <node5#sss:optout_main> CVodeGetErrWeights
  @cvode <node3#ss:ivp_sol> IVP solution (W_i)
 *)
val get_err_weights : session -> nvec -> unit

(**
  Returns the vector of estimated local errors.

  @cvode <node5#sss:optout_main> CVodeGetEstLocalErrors
 *)
val get_est_local_errors : session -> nvec -> unit

(**
  Returns the integrator statistics as a group.

  @cvode <node5#sss:optout_main> CVodeGetIntegratorStats
 *)
val get_integrator_stats    : session -> Cvode.integrator_stats

(**
  Convenience function that calls get_integrator_stats and prints the results to
  stdout.

  @cvode <node5#sss:optout_main> CVodeGetIntegratorStats
 *)
val print_integrator_stats  : session -> unit


(**
  Returns the number of nonlinear (functional or Newton) iterations performed.

  @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvIters
 *)
val get_num_nonlin_solv_iters : session -> int

(**
  Returns the number of nonlinear convergence failures that have occurred.

  @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails
 *)
val get_num_nonlin_solv_conv_fails : session -> int

(** {2 Root finding optional functions} *)

(** {3 Input} *)

(**
  [set_root_direction s dir] specifies the direction of zero-crossings to be
  located and returned. [dir] may contain one entry for each root function.

  @cvode <node5#sss:optin_root> CVodeSetRootDirection
 *)
val set_root_direction : session -> root_direction array -> unit

(**
  Like {!set_root_direction} but specifies a single direction for all root
  functions.

  @cvode <node5#sss:optin_root> CVodeSetRootDirection
 *)
val set_all_root_directions : session -> root_direction -> unit

(**
  Disables issuing a warning if some root function appears to be identically
  zero at the beginning of the integration.

  @cvode <node5#sss:optin_root> CVodeSetNoInactiveRootWarn
 *)
val set_no_inactive_root_warn : session -> unit

(** {3 Output} *)

(**
  Fills an array showing which functions were found to have a root.

  @cvode <node5#sss:optout_root> CVodeGetRootInfo
 *)
val get_root_info : session -> root_array -> unit

(**
  Returns the cumulative number of calls made to the user-supplied root function g.

  @cvode <node5#sss:optout_root> CVodeGetNumGEvals
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

  @cvode <node5#sss:optin_root> CVodeGetDky
 *)
val get_dky : session -> float -> int -> nvec -> unit

(** {2 Reinitialization} *)

(**
  [reinit s t0 y0] reinitializes the solver session [s] with a new time [t0] and
  new values for the variables [y0].

  @cvode <node5#sss:cvreinit> CVodeReInit
 *)
val reinit : session -> float -> nvec -> unit


(** {2 Linear Solvers} *)

(** {b TODO} *)
type 't jacobian_arg =
  {
    jac_t   : float;
    jac_y   : val_array;
    jac_fy  : val_array;
    jac_tmp : 't
  }

type triple_tmp = val_array * val_array * val_array

(** {3 Direct Linear Solvers (DLS)} *)

module Dls :
  sig
    val set_dense_jac_fn :
         session
      -> (triple_tmp jacobian_arg -> Densematrix.t -> unit)
      -> unit

    val clear_dense_jac_fn : session -> unit

    val set_band_jac_fn :
         session
      -> (triple_tmp jacobian_arg -> int -> int -> Bandmatrix.t -> unit)
      -> unit

    val clear_band_jac_fn : session -> unit

    val get_work_space : session -> int * int

    (* No. of Jacobian evaluations *)
    val get_num_jac_evals : session -> int

    (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
    val get_num_rhs_evals : session -> int
  end

(** {3 Diagonal approximation} *)

module Diag :
  sig
    val get_work_space : session -> int * int

    (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
    val get_num_rhs_evals : session -> int
  end

(** {3 Banded preconditioning} *)

module BandPrec :
  sig
    val get_work_space : session -> int * int

    (* No. of r.h.s. calls for finite diff. banded Jacobian evals. *)
    val get_num_rhs_evals : session -> int
  end

(** {3 Scaled Preconditioned Iterative Linear Solvers (SPILS)} *)

module Spils :
  sig
    type solve_arg =
      {
        rhs   : val_array;
        gamma : float;
        delta : float;
        left  : bool; (* true: left, false: right *)
      }

    type single_tmp = nvec

    type gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    val set_preconditioner :
      session
      -> (triple_tmp jacobian_arg -> bool -> float -> bool)
      -> (single_tmp jacobian_arg -> solve_arg -> nvec -> unit)
      -> unit

    val set_jac_times_vec_fn :
      session
      -> (single_tmp jacobian_arg
          -> val_array (* v *)
          -> val_array (* Jv *)
          -> unit)
      -> unit
    val clear_jac_times_vec_fn : session -> unit

    val set_prec_type : session -> preconditioning_type -> unit

    val set_gs_type :
      session -> gramschmidt_type -> unit

    val set_eps_lin : session -> float -> unit

    val set_maxl : session -> int -> unit

    val get_work_space       : session -> int * int
    val get_num_prec_evals   : session -> int
    val get_num_prec_solves  : session -> int
    val get_num_lin_iters    : session -> int
    val get_num_conv_fails   : session -> int
    val get_num_jtimes_evals : session -> int
    val get_num_rhs_evals    : session -> int
  end

