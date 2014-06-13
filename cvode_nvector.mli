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
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

include module type of Cvode
  with type Roots.t = Cvode.Roots.t
  and type RootDirs.t = Cvode.RootDirs.t
  and type lmm = Cvode.lmm
  and type solver_result = Cvode.solver_result
  and type integrator_stats = Cvode.integrator_stats
  and type bandrange = Cvode.bandrange
  and type spils_params = Cvode.spils_params

(** Abstract nvector interface to the CVODE Solver.

  Nvectors can be created for an {!Nvector.Mutable} or {!Nvector.Immutable}
  data type by providing a set of low-level operations which are called by
  the solver.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

(**
    This type represents a session with the CVODE solver using
    {!Nvector.nvector}s.

    A skeleton of the main program:
    + {b Set vector of initial values}
    {[let y = Nvector_array.wrap [| 0.0; 0.0; 0.0 |] ]}
    The length of this vector determines the problem size.    
    + {b Create and initialize a solver session}
    {[let s = Cvode.init Cvode.Adams Cvode.Functional tols f ~roots:(2, g) y]}
    This will initialize a specific linear solver and the root-finding
    mechanism, if necessary.
    + {b Specify integration tolerances (optional)}, e.g.
    {[set_tolerances s SStolerances (reltol, abstol)]}
    + {b Set optional inputs}, e.g.
    {[set_stop_time s 10.0; ...]}
    Call any of the [set_*] functions to change solver parameters from their
    defaults.
    + {b Advance solution in time}, e.g.
    {[let (t', result) = Cvode.solve_normal s !t y in
...
t := t' + 0.1]}
    Repeatedly call either [solve_normal] or [solve_one_step] to advance the
    simulation.
    + {b Get optional outputs}
    {[let stats = get_integrator_stats s in ...]}
    Call any of the [get_*] functions to examine solver statistics.

    @cvode <node5#ss:skeleton_sim> Skeleton of main program
 *)
type 'a session = 'a Cvode_session_nvector.session

(** The type of vectors passed to the solver. *)
type 'a nvector = 'a Nvector.nvector

(** The type of vectors containing detected roots (zero-crossings). *)
type root_array = Sundials.Roots.t

(** The type of vectors containing the values of root functions
   (zero-crossings). *)
type root_val_array = Sundials.Roots.val_array
(** {2 Linear Solvers} *)

type 'a single_tmp = 'a
type 'a triple_tmp = 'a * 'a * 'a

(**
  Arguments common to all Jacobian callback functions.    
 
  @cvode <node5#ss:jtimesFn> Jacobian-times-vector function
  @cvode <node5#ss:psolveFn> Linear preconditioning function
  @cvode <node5#ss:precondFn> Jacobian preconditioning function
 *)
type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;        (** The independent variable. *)
    jac_y   : 'a;           (** The dependent variable vector. *)
    jac_fy  : 'a;           (** The derivative vector (i.e., f(t, y)). *)
    jac_tmp : 't            (** Workspace data,
                                either {!single_tmp} or {!triple_tmp}. *)
  }

(** {3 Setting linear solvers to sessions} *)

(**
 Specify a solution method.

 @cvode <node3#ss:ivp_sol> IVP Solution
 @cvode <node5#sss:cvodemalloc> CVodeCreate
 *)
type 'a iter =
  | Newton of 'a linear_solver (** Newton iteration with a given linear solver *)
  | Functional              (** Functional iteration (non-stiff systems only) *)
(**
   Specify a linear solver.

   The Banded Krylov solvers imply an additional call to
   {{:CVODE_DOC_ROOT(node5#sss:cvbandpre)} CVBandPrecInit}.

   @cvode <node5#sss:lin_solv_init> Linear Solver Specification Functions
*)
and 'a linear_solver =
  | Diag
  (** Diagonal approximation of the Jacobian by difference quotients.

      @cvode <node5#sss:lin_solve_init> CVDiag
    *)
  | Spgmr of Cvode.spils_params * 'a spils_callbacks
  (** Krylov iterative solver with the scaled preconditioned GMRES method.  The
      arguments specify the maximum dimension of the Krylov subspace and
      preconditioning type ({!Cvode.spils_params}) and the preconditioner
      callback functions ({!spils_callbacks}).  See also {!Spils}.

      @cvode <node5#sss:lin_solve_init> CVSpgmr
      @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
      @cvode <node5#ss:psolveFn> Linear preconditioning function
      @cvode <node5#ss:precondFn> Jacobian preconditioning function
    *)
  | Spbcg of Cvode.spils_params * 'a spils_callbacks
  (** Krylov iterative solver with the scaled preconditioned Bi-CGStab method.
      The arguments are the same as [Spgmr].  See also {!Spils}.

      @cvode <node5#sss:lin_solve_init> CVSpbcg
      @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
      @cvode <node5#ss:psolveFn> Linear preconditioning function
      @cvode <node5#ss:precondFn> Jacobian preconditioning function
      (* TODO: change all these titles into function names? *)
    *)
  | Sptfqmr of Cvode.spils_params * 'a spils_callbacks
  (** Krylov iterative with the scaled preconditioned TFQMR method.  The
      arguments are the same as [Spgmr].  See also {!Spils}.

      @cvode <node5#sss:lin_solve_init> CVSptfqmr
      @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
      @cvode <node5#ss:psolveFn> Linear preconditioning function
      @cvode <node5#ss:precondFn> Jacobian preconditioning function
    *)

(** Callbacks for Krylov subspace linear solvers.  Ignored if the
    {!Spils.preconditioning_type} is set to [PrecNone].  In that case, you
    should use {!spils_no_precond} as [spils_callbacks].  *)
and 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                     -> 'a -> unit) option;
    (** Called like [prec_solve_fn jac_arg solve_arg z] to solve the linear
        system {i P}[z] = [solve_arg.rhs], where {i P} may be either a left or
        right preconditioner matrix.  {i P} should approximate, however
        crudely, the Newton matrix [M = I - jac.gamma * J], where {i J} =
        delr({i f}) / delr({i y}).

        - [jac_arg] supplies the basic problem data as a {!jacobian_arg}.
        - [solve_arg] specifies the linear system as a {!prec_solve_arg}.
        - [z] is the vector in which the result must be stored.

        The function may raise a {!Sundials.RecoverableFailure} exception to
        indicate that a recoverable error has occurred. Any other exception is
        treated as an unrecoverable error.

        {b NB:} The elements of [jac], [arg], and [z] must no longer be
        accessed after [psolve] has returned a result, i.e. if their values are
        needed outside of the function call, then they must be copied to
        separate physical structures.

        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
    *)

    prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool)
                    option;
    (** A function that preprocesses and/or evaluates any Jacobian-related data
        needed by [prec_solve_fn] above.  When [prec_solve_fn] doesn't need any
        such data, this field can be [None].

        This callback is invoked as [prec_setup_fn arg jok gamma], where the
        arguments have the following meaning:
        - [arg] supplies the basic problem data as a {!jacobian_arg}.
        - [jok] indicates whether any saved Jacobian-related data can be
          reused.  If [false] any such data must be recomputed from scratch,
          otherwise, if [true], any such data saved from a previous call to the
          function can be reused, with the current value of [gamma].  A call
          with [jok = true] can only happen after an earlier call with [jok =
          false].
        - [gamma] is the scalar {i g} appearing in the Newton matrix
          {i M = I - gJ} where {i I} is the identity and {i J} is the Jacobian.
        This function must return [true] if the Jacobian-related data was
        updated, or [false] otherwise, i.e. if the saved data was reused.

        The function may raise a {!Sundials.RecoverableFailure} exception to
        indicate that a recoverable error has occurred. Any other exception is
        treated as an unrecoverable error.

        {b NB:} The fields of [arg] must no longer be accessed after this
        callback has returned, i.e. if their values are needed outside of the
        function call, then they must be copied to a separate physical
        structure.

        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn
     *)

    jac_times_vec_fn :
      (('a single_tmp, 'a) jacobian_arg
       -> 'a (* v *)
       -> 'a (* Jv *)
       -> unit) option;
    (** Specifies a function that computes the Jacobian times a vector.

        The function is called like [jac_times_vec_fn arg v jv] and should
        compute the matrix-vector product {i J}[v], where {i J} is the
        system Jacobian (not explicitly constructed).
        - [arg] is the standard {!jacobian_arg} with one work vector.
        - [v] is the vector to multiply to the Jacobian.
        - [jv] is the vector in which to store the result.
        The ({i i,j}) entry of the Jacobian {i J} is {i dFi/dyj}, i.e. the
        partial derivative of the {i i}-th equation with respect to the
        {i j}-th variable.

        The function may raise a {!Sundials.RecoverableFailure} exception to
        indicate that a recoverable error has occurred. Any other exception is
        treated as an unrecoverable error.

        {b NB:} The elements of [arg], [v], and [jv] must no longer be
        accessed after this callback has returned, i.e. if their values are
        needed outside of the function call, then they must be copied to a
        separate physical structure.

        @cvode <node5#ss:jtimesfn> CVSpilsJacTimesVecFn
        @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
     *)
  }

(** Arguments passed to the preconditioner solve callback function.  See
    [prec_solve_fn] in {!spils_callbacks}.

    @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
 *)
and 'a prec_solve_arg =
  {
    rhs   : 'a;         (** The right-hand side vector, {i r}, of the
                            linear system. *)
    gamma : float;      (** The scalar {i g} appearing in the Newton
                            matrix given by M = I - {i g}J. *)
    delta : float;      (** Input tolerance to be used if an
                            iterative method is employed in the
                            solution. *)
    left  : bool;       (** [true] if the left preconditioner
                            is to be used and [false] if the
                            right preconditioner is to be used. *)
  }

(** See {!spils_callbacks}.  *)
val spils_no_precond : 'a spils_callbacks

(** {3 Diagonal approximation} *)

(** Get optional inputs for the linear solver that gives diagonal approximations
    of the Jacobian matrix.
    @cvode <node5#sss:optout_diag> Diagonal linear solver optional output functions
  *)
module Diag :
  sig
    (** {4 Optional input functions} *)

    (**
      Returns the sizes of the real and integer workspaces used by the Diagonal
      linear solver.

      @cvode <node5#sss:optout_diag> CVDiagGetWorkSpace
      @return ([real_size], [integer_size])
     *)
    val get_work_space : 'a session -> int * int

    (**
      Returns the number of calls made to the user-supplied right-hand side
      function due to finite difference Jacobian approximation in the Diagonal
      linear solver.

      @cvode <node5#sss:optout_diag> CVDiagGetNumRhsEvals
    *)
    val get_num_rhs_evals : 'a session -> int
  end

(** {3 Scaled Preconditioned Iterative Linear Solvers (SPILS)} *)

(** Set callback functions, set optional outputs, and get optional inputs for
    the Scaled Preconditioned Iterative Linear Solvers: SPGMR, SPBCG, SPTFQMR.
    @cvode <node5#sss:optin_spils> Iterative linear solvers optional input functions.
    @cvode <node5#sss:optout_spils> Iterative linear solvers optional output functions.
    @cvode <node5#ss:psolveFn> Linear preconditioning function
    @cvode <node5#ss:precondFn> Jacobian preconditioning function
 *)
module Spils :
  sig
    (** {4 Low-level solver manipulation} *)

    (** The type of one of the parameters passed to the preconditioner
        function.  See {!spils_callbacks}.  *)
    type 'a solve_arg = 'a prec_solve_arg =
      {
        rhs   : 'a;
        gamma : float;
        delta : float;
        left  : bool;
      }

    (**
      Set preconditioning functions (see {!spils_callbacks}).  It may be unsafe
      to use this function without a {!reinit}.  Users are encouraged to use
      the [iter_type] parameter of {!reinit} instead, unless they are desperate
      for performance.

      @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
      @cvode <node5#ss:psolveFn> Linear preconditioning function
      @cvode <node5#ss:precondFn> Jacobian preconditioning function
    *)
    val set_preconditioner :
      'a session
      -> (('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool)
      -> (('a single_tmp, 'a) jacobian_arg -> 'a solve_arg -> 'a -> unit)
      -> unit

    (**
      Set the Jacobian-times-vector function (see {!Cvode.spils_params}).  It
      may be unsafe to use this function without a {!reinit}.  Users are
      encouraged to use the [iter_type] parameter of {!reinit} instead, unless
      they are desperate for performance.

      @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
      @cvode <node5#ss:jtimesFn> Jacobian-times-vector function
    *)
    val set_jac_times_vec_fn :
      'a session
      -> (('a single_tmp, 'a) jacobian_arg
          -> 'a (* v *)
          -> 'a (* Jv *)
          -> unit)
      -> unit

    (**
      This function restores the default Jacobian-times-vector function. It is
      equivalent to calling CVodeSetJacTimesVecFn with an argument of [NULL]. It
      may be unsafe to use this function without a {!reinit}.  Users are
      encouraged to use the [iter_type] parameter of {!reinit} instead, unless
      they are desperate for performance.

      @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
      @cvode <node5#ss:jtimesFn> Jacobian-times-vector function
    *)
    val clear_jac_times_vec_fn : 'a session -> unit

    (** {4 Optional input functions} *)

    (**
      This function resets the type of preconditioning to be used using a value
      of type {!Spils.preconditioning_type}.

      @cvode <node5#sss:optin_spils> CVSpilsSetPrecType
    *)
    val set_prec_type : 'a session -> Spils.preconditioning_type -> unit

    (**
      Sets the Gram-Schmidt orthogonalization to be used with the
      Spgmr {!linear_solver}.

      @cvode <node5#sss:optin_spils> CVSpilsSetGSType
    *)
    val set_gs_type : 'a session -> Spils.gramschmidt_type -> unit

    (**
      [set_eps_lin eplifac] sets the factor by which the Krylov linear solver's
      convergence test constant is reduced from the Newton iteration test
      constant. [eplifac]  must be >= 0. Passing a value of 0 specifies the
      default (which is 0.05).

      @cvode <node5#sss:optin_spils> CVSpilsSetEpsLin
    *)
    val set_eps_lin : 'a session -> float -> unit

    (**
      [set_maxl maxl] resets the maximum Krylov subspace dimension for the
      Bi-CGStab or TFQMR methods. [maxl] is the maximum dimension of the Krylov
      subspace, a value of [maxl] <= 0 specifies the default (which is 5.0).

      @cvode <node5#sss:optin_spils> CVSpilsSetMaxl
    *)
    val set_maxl : 'a session -> int -> unit

    (** {4 Optional output functions} *)

    (**
      Returns the sizes of the real and integer workspaces used by the SPGMR
      linear solver.

      @cvode <node5#sss:optout_spils> CVSpilsGetWorkSpace
      @return ([real_size], [integer_size])
    *)
    val get_work_space       : 'a session -> int * int

    (**
      Returns the cumulative number of linear iterations.

      @cvode <node5#sss:optout_spils> CVSpilsGetNumLinIters
    *)
    val get_num_lin_iters    : 'a session -> int

    (**
      Returns the cumulative number of linear convergence failures.

      @cvode <node5#sss:optout_spils> CVSpilsGetNumConvFails
    *)
    val get_num_conv_fails   : 'a session -> int

    (**
      Returns the number of preconditioner evaluations, i.e., the number of
      calls made to psetup with jok = [false] (see {!set_preconditioner}).

      @cvode <node5#sss:optout_spils> CVSpilsGetNumPrecEvals
    *)
    val get_num_prec_evals   : 'a session -> int

    (**
      Returns the cumulative number of calls made to the preconditioner solve
      function, psolve (see {!set_preconditioner}).

      @cvode <node5#sss:optout_spils> CVSpilsGetNumPrecSolves
    *)
    val get_num_prec_solves  : 'a session -> int

    (**
      Returns the cumulative number of calls made to the Jacobian-vector
      function, jtimes (see {! set_jac_times_vec_fn}).

      @cvode <node5#sss:optout_spils> CVSpilsGetNumJtimesEvals
    *)
    val get_num_jtimes_evals : 'a session -> int

    (**
      Returns the number of calls to the user right-hand side function for
      finite difference Jacobian-vector product approximation. This counter is
      only updated if the default difference quotient function is used.

      @cvode <node5#sss:optout_spils> CVSpilsGetNumRhsEvals
    *)
    val get_num_rhs_evals    : 'a session -> int
  end

(** {3 Banded preconditioner} *)

(** Get optional outputs for the banded preconditioner module of the
    Scaled Preconditioned Iterative Linear Solvers:
      SPGMR, SPBCG, SPTFQMR.
    @cvode <node5#sss:cvbandpre> Serial banded preconditioner module
  *)
module BandPrec :
  sig
    (** {4 Optional output functions} *)

    (**
      Returns the sizes of the real and integer workspaces used by the serial
      banded preconditioner module.

      @cvode <node5#sss:cvbandpre> CVBandPrecGetWorkSpace
      @return ([real_size], [integer_size])
     *)
    val get_work_space : 'a session -> int * int

    (**
      Returns the number of calls made to the user-supplied right-hand side
      function due to finite difference banded Jacobian approximation in the
      banded preconditioner setup function.

      @cvode <node5#sss:cvbandpre> CVBandPrecGetNumRhsEvals
    *)
    val get_num_rhs_evals : 'a session -> int
  end

type 'a tolerance =
  | SStolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVtolerances of float * 'a nvector
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFtolerances of ('a -> 'a -> unit)
    (** Specifies a function [efun y ewt] that sets the multiplicative
        error weights Wi for use in the weighted RMS norm. The function is
        passed the dependent variable vector [y] and is expected to set the
        values inside the error-weight vector [ewt]. *)

(** A default relative tolerance of 1.0e-4 and absolute tolerance of 1.0e-8. *)
val default_tolerances : 'a tolerance

(** {2 Initialization} *)

(**
    [init lmm iter tol f ~roots:(nroots, g) ~t0:t0 (neqs, y0)] initializes the
    CVODE solver and returns a {!session}.
    - [lmm]     specifies the linear multistep method, see {!Cvode.lmm},
    - [iter]    specifies either functional iteration or Newton iteration
                with a specific linear solver, see {!iter},
    - [tol]     specifies the integration tolerances,
    - [f]       is the ODE right-hand side function,
    - [nroots]  specifies the number of root functions (zero-crossings).
    - [g]       calculates the values of the root functions,
    - [t0]      is the initial value of the independent variable,
    - [neqs]    specifies the number of equations (continuous state
                variables). There is no operation for demanding the length
                of an 'a nvector, and,
    - [y0]      is a vector of initial values, the size of this vector
                determines the number of equations in the session, see
                {!Sundials.Carray.t}.

    The labeled arguments [roots] and [t0] are both optional and default to
    {!Cvode.no_roots} (i.e. no root finding is done) and [0.0], respectively.

    This function calls CVodeCreate, CVodeInit, CVodeRootInit, an appropriate
    linear solver function, and one of CVodeSStolerances, CVodeSVtolerances, or
    CVodeWFtolerances. It does everything necessary to initialize a CVODE
    session; the {!solve_normal} or {!solve_one_step} functions can be called
    directly afterward.

    The right-hand side function [f] is called by the solver to calculate the
    instantaneous derivative values, and is passed three arguments: [t], [y],
    and [dy].
    - [t] is the current value of the independent variable,
          i.e., the simulation time.
    - [y] is a vector of dependent-variable values, i.e. y(t).
    - [dy] is a vector for storing the value of f(t, y).
    The function may raise a {!Sundials.RecoverableFailure} exception to
    indicate that a recoverable error has occurred. Any other exception is
    treated as an unrecoverable error.

    {b NB:} [y] and [dy] must no longer be accessed after [f] has returned a
            result, i.e. if their values are needed outside of the function
            call, then they must be copied to separate physical structures.

    The roots function [g] is called by the solver to calculate the values of
    root functions (zero-crossing expressions) which are used to detect
    significant events, it is passed three arguments: [t], [y], and [gout].
    - [t] and [y] are as for [f].
    - [gout] is a vector for storing the values of g(t, y).
    The {!Cvode.no_roots} value can be passed for the [(nroots, g)] argument if
    root functions are not required.

    {b NB:} [y] and [gout] must no longer be accessed after [g] has returned
            a result, i.e. if their values are needed outside of the function
            call, then they must be copied to separate physical structures.

    @cvode <node5#sss:cvodemalloc>   CVodeCreate/CVodeInit
    @cvode <node5#ss:rhsFn>          ODE right-hand side function
    @cvode <node5#ss:cvrootinit>     CVodeRootInit
    @cvode <node5#ss:rootFn>         Rootfinding function
    @cvode <node5#sss:lin_solv_init> Linear solvers
    @cvode <node5#sss:cvtolerances>  CVodeSStolerances
    @cvode <node5#sss:cvtolerances>  CVodeSVtolerances
    @cvode <node5#sss:cvtolerances>  CVodeWFtolerances
    @cvode <node5#ss:ewtsetFn>       Error weight function
 *)
val init :
    lmm
    -> 'a iter
    -> 'a tolerance
    -> (float -> 'a -> 'a -> unit)
    -> ?roots:(int * (float -> 'a -> root_val_array -> unit))
    -> ?t0:float
    -> 'a nvector
    -> 'a session

(** Return the number of root functions. *)
val nroots : 'a session -> int

(** {2 Solver functions } *)

(**
   [(tret, r) = solv_normal s tout yout] integrates the ODE over an interval
   in t.

   The arguments are:
   - [s] a session with the solver.
   - [tout] the next time at which a computed solution is desired.
   - [yout] a vector to store the computed solution. The same vector as was
   passed to {!init} can be used again for this argument.

   Two values are returned:
    - [tret] the time reached by the solver, which will be equal to [tout] if
      no errors occur.
    - [r] indicates whether roots were found, or whether an optional stop time, set by
    {!set_stop_time}, was reached; see {!Sundials.solver_result}.

   This routine will throw one of the solver {!Cvode.exceptions} if an error
   occurs.

   @cvode <node5#sss:cvode> CVode (CV_NORMAL)
 *)
val solve_normal : 'a session -> float -> 'a nvector -> float * solver_result

(**
   This function is identical to {!solve_normal}, except that it returns after
   one internal solver step.

   @cvode <node5#sss:cvode> CVode (CV_ONE_STEP)
 *)
val solve_one_step : 'a session -> float -> 'a nvector -> float * solver_result

(** {2 Main optional functions} *)

(** {3 Input} *)

(** Set the integration tolerances.

    @cvode <node5#sss:cvtolerances> CVodeSStolerances
    @cvode <node5#sss:cvtolerances> CVodeSVtolerances
    @cvode <node5#sss:cvtolerances> CVodeWFtolerances
    @cvode <node5#ss:ewtsetFn> Error weight function
 *)
val set_tolerances : 'a session -> 'a tolerance -> unit

(**
  [set_error_file s fname trunc] opens the file named [fname] and to which all
  messages from the default error handler are then directed.
  If the file already exists it is either trunctated ([trunc] = [true]) or
  appended to ([trunc] = [false]).

  The error file is closed if set_error_file is called again, or otherwise when
  the session is garbage collected.
   
  @cvode <node5#sss:optin_main> CVodeSetErrFile
 *)
val set_error_file : 'a session -> string -> bool -> unit

(**
  [set_err_handler_fn s efun] specifies a custom function [efun] for handling
  error messages.

  @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
  @cvode <node5#ss:ehFn> Error message handler function
 *)
val set_err_handler_fn : 'a session -> (Sundials.error_details -> unit) -> unit

(**
  This function restores the default error handling function. It is equivalent
  to calling CVodeSetErrHandlerFn with an argument of [NULL].

  @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
 *)
val clear_err_handler_fn : 'a session -> unit

(**
  Specifies the maximum order of the linear multistep method.

  @cvode <node5#sss:optin_main> CVodeSetMaxOrd
 *)
val set_max_ord : 'a session -> int -> unit

(**
  Specifies the maximum number of steps to be taken by the solver in its attempt
  to reach the next output time.

  @cvode <node5#sss:optin_main> CVodeSetMaxNumSteps
 *)
val set_max_num_steps : 'a session -> int -> unit

(**
  Specifies the maximum number of messages issued by the solver warning that t +
  h = t on the next internal step.

  @cvode <node5#sss:optin_main> CVodeSetMaxHnilWarns
 *)
val set_max_hnil_warns : 'a session -> int -> unit

(**
  Indicates whether the BDF stability limit detection algorithm should be used.

  @cvode <node5#sss:optin_main> CVodeSetStabLimDet
  @cvode <node3#s:bdf_stab> BDF Stability Limit Detection
 *)
val set_stab_lim_det : 'a session -> bool -> unit

(**
  Specifies the initial step size.

  @cvode <node5#sss:optin_main> CVodeSetInitStep
 *)
val set_init_step : 'a session -> float -> unit

(**
  Specifies a lower bound on the magnitude of the step size.

  @cvode <node5#sss:optin_main> CVodeSetMinStep
 *)
val set_min_step : 'a session -> float -> unit

(**
  Specifies an upper bound on the magnitude of the step size.

  @cvode <node5#sss:optin_main> CVodeSetMaxStep
 *)
val set_max_step : 'a session -> float -> unit

(**
  Specifies the value of the independent variable t past which the solution is
  not to proceed.
  The default, if this routine is not called, is that no stop time is imposed.

  @cvode <node5#sss:optin_main> CVodeSetStopTime
 *)
val set_stop_time : 'a session -> float -> unit

(**
  Specifies the maximum number of error test failures permitted in attempting
  one step.

  @cvode <node5#sss:optin_main> CVodeSetMaxErrTestFails
 *)
val set_max_err_test_fails : 'a session -> int -> unit

(**
  Specifies the maximum number of nonlinear solver iterations permitted per
  step.

  @cvode <node5#sss:optin_main> CVodeSetMaxNonlinIters
 *)
val set_max_nonlin_iters : 'a session -> int -> unit

(**
  Specifies the maximum number of nonlinear solver convergence failures
  permitted during one step.

  @cvode <node5#sss:optin_main> CVodeSetMaxConvFails
 *)
val set_max_conv_fails : 'a session -> int -> unit

(**
  Specifies the safety factor used in the nonlinear convergence test.

  @cvode <node5#sss:optin_main> CVodeSetNonlinConvCoef
  @cvode <node3#ss:ivp_sol> IVP Solution
 *)
val set_nonlin_conv_coef : 'a session -> float -> unit

(** {3 Output } *)

(**
  Returns the real and integer workspace sizes.

  @cvode <node5#sss:optout_main> CVodeGetWorkSpace
  @return ([real_size], [integer_size])
 *)
val get_work_space          : 'a session -> int * int

(**
  Returns the cumulative number of internal steps taken by the solver.

  @cvode <node5#sss:optout_main> CVodeGetNumSteps
 *)
val get_num_steps           : 'a session -> int

(**
  Returns the number of calls to the user's right-hand side function.

  @cvode <node5#sss:optout_main> CVodeGetNumRhsEvals
 *)
val get_num_rhs_evals       : 'a session -> int

(**
  Returns the number of calls made to the linear solver's setup function.

  @cvode <node5#sss:optout_main> CVodeGetNumLinSolvSetups
 *)
val get_num_lin_solv_setups : 'a session -> int

(**
  Returns the number of local error test failures that have occurred.

  @cvode <node5#sss:optout_main> CVodeGetNumErrTestFails
 *)
val get_num_err_test_fails  : 'a session -> int

(**
  Returns the integration method order used during the last internal step.

  @cvode <node5#sss:optout_main> CVodeGetLastOrder
 *)
val get_last_order          : 'a session -> int

(**
  Returns the integration method order to be used on the next internal step.

  @cvode <node5#sss:optout_main> CVodeGetCurrentOrder
 *)
val get_current_order       : 'a session -> int

(**
  Returns the integration step size taken on the last internal step.

  @cvode <node5#sss:optout_main> CVodeGetLastStep
 *)
val get_last_step           : 'a session -> float

(**
  Returns the integration step size to be attempted on the next internal step.

  @cvode <node5#sss:optout_main> CVodeGetCurrentStep
 *)
val get_current_step        : 'a session -> float

(**
  Returns the the value of the integration step size used on the first step.

  @cvode <node5#sss:optout_main> CVodeGetActualInitStep
 *)
val get_actual_init_step    : 'a session -> float

(**
  Returns the the current internal time reached by the solver.

  @cvode <node5#sss:optout_main> CVodeGetCurrentTime
 *)
val get_current_time        : 'a session -> float

(**
  Returns the number of order reductions dictated by the BDF stability limit
  detection algorithm.

  @cvode <node5#sss:optout_main> CVodeGetNumStabLimOrderReds
  @cvode <node3#s:bdf_stab> BDF stability limit detection
 *)
val get_num_stab_lim_order_reds : 'a session -> int

(**
  Returns a suggested factor by which the user's tolerances should be scaled
  when too much accuracy has been requested for some internal step.

  @cvode <node5#sss:optout_main> CVodeGetTolScaleFactor
 *)
val get_tol_scale_factor : 'a session -> float

(**
  Returns the solution error weights at the current time.

  @cvode <node5#sss:optout_main> CVodeGetErrWeights
  @cvode <node3#ss:ivp_sol> IVP solution (W_i)
 *)
val get_err_weights : 'a session -> 'a nvector -> unit

(**
  Returns the vector of estimated local errors.

  @cvode <node5#sss:optout_main> CVodeGetEstLocalErrors
 *)
val get_est_local_errors : 'a session -> 'a nvector -> unit

(**
  Returns the integrator statistics as a group.

  @cvode <node5#sss:optout_main> CVodeGetIntegratorStats
 *)
val get_integrator_stats    : 'a session -> Cvode.integrator_stats

(**
  Convenience function that calls get_integrator_stats and prints the results to
  stdout.

  @cvode <node5#sss:optout_main> CVodeGetIntegratorStats
 *)
val print_integrator_stats  : 'a session -> unit


(**
  Returns the number of nonlinear (functional or Newton) iterations performed.

  @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvIters
 *)
val get_num_nonlin_solv_iters : 'a session -> int

(**
  Returns the number of nonlinear convergence failures that have occurred.

  @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails
 *)
val get_num_nonlin_solv_conv_fails : 'a session -> int

(**
  [nniters, nncfails = get_nonlin_solv_stats s] returns both the numbers of
  nonlinear iterations performed [nniters] and of nonlinear convergence
  failures that have occurred [nncfails].

  @cvode <node5#sss:optout_main> CVodeGetNonlinSolvStats
 *)
val get_nonlin_solv_stats : 'a session -> int *int

(** {2 Root finding optional functions} *)

(** {3 Input} *)

(**
  [set_root_direction s dir] specifies the direction of zero-crossings to be
  located and returned. [dir] may contain one entry of type
  {!Cvode.root_direction} for each root function.

  @cvode <node5#sss:optin_root> CVodeSetRootDirection
 *)
val set_root_direction : 'a session -> root_direction array -> unit

(**
  Like {!set_root_direction} but specifies a single direction of type
  {!Cvode.root_direction} for all root
  functions.

  @cvode <node5#sss:optin_root> CVodeSetRootDirection
 *)
val set_all_root_directions : 'a session -> root_direction -> unit

(**
  Disables issuing a warning if some root function appears to be identically
  zero at the beginning of the integration.

  @cvode <node5#sss:optin_root> CVodeSetNoInactiveRootWarn
 *)
val set_no_inactive_root_warn : 'a session -> unit

(** {3 Output} *)

(**
  Fills an array showing which functions were found to have a root.

  @cvode <node5#sss:optout_root> CVodeGetRootInfo
 *)
val get_root_info : 'a session -> root_array -> unit

(**
  Returns the cumulative number of calls made to the user-supplied root function g.

  @cvode <node5#sss:optout_root> CVodeGetNumGEvals
 *)
val get_num_g_evals : 'a session -> int

(** {2 Interpolated output function } *)

(**
  [get_dky s t k dky] computes the [k]th derivative of the function y at time
  [t], i.e. d(k)y/dt(k)(t). The function requires that tn - hu <= [t] <=
  tn, where tn denotes the current internal time reached, and hu is the last
  internal step size successfully used by the solver.
  The user may request [k] = 0, 1,..., qu, where qu is the current order.

  This function may only be called after a successful return from either
  {!solve_normal} or {!solve_one_step}.

  Values for the limits may be obtained:
    - tn = {!get_current_time}
    - qu = {!get_last_order}
    - hu = {!get_last_step}

  @cvode <node5#sss:optin_root> CVodeGetDky
 *)
val get_dky : 'a session -> float -> int -> 'a nvector -> unit

(** {2 Reinitialization} *)

(**

  [reinit s ~iter_type:iter_type ~roots:roots t0 y0] reinitializes the
  solver session [s] with new values for the variables [y0].

  The labeled arguments are all optional, and if omitted, their current
  settings are kept.

  [iter_type] sets the nonlinear solver iteration type.  If omitted, the
  current iteration type will be kept.

  [roots] sets the root functions; see {!init} for what each component does.
  {!Cvode.no_roots} may be specified to turn off root finding (if the session
  has been doing root finding until now).  If omitted, the current root
  functions or lack thereof will be kept.

  [t0] sets the value of the independent variable.  If omitted, the current
  time value will be kept.

  @cvode <node5#sss:cvreinit> CVodeReInit
 *)
val reinit :
  'a session
  -> ?iter_type:'a iter
  -> ?roots:(int * (float -> 'a -> root_val_array -> unit))
  -> float
  -> 'a nvector
  -> unit


