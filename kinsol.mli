(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
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

(** KINSOL solves nonlinear systems using Newton-Krylov techniques.

 @version VERSION()
 @author Timothy Bourke (Inria)
 @author Jun Inoue (Inria)
 @author Marc Pouzet (LIENS)
 *)
open Kinsol_impl

(** {3:exceptions Exceptions} *)

(** @kinsol <node5#sss:kinsol> KIN_ILL_INPUT *)
exception IllInput

(** @kinsol <node5#sss:kinsol> KIN_LINESEARCH_NONCONV *)
exception LineSearchNonConvergence

(** @kinsol <node5#sss:kinsol> KIN_MAXITER_REACHED *)
exception MaxIterationsReached

(** @kinsol <node5#sss:kinsol> KIN_MXNEWT_5X_EXCEEDED *)
exception MaxNewtonStepExceeded

(** @kinsol <node5#sss:kinsol> KIN_LINESEARCH_BCFAIL *)
exception LineSearchBetaConditionFailure

(** @kinsol <node5#sss:kinsol> KIN_LINSOLV_NO_RECOVERY *)
exception LinearSolverNoRecovery

(** @kinsol <node5#sss:kinsol> KIN_LINIT_FAIL *)
exception LinearSolverInitFailure

(** @kinsol <node5#sss:kinsol> KIN_LSETUP_FAIL *)
exception LinearSetupFailure

(** @kinsol <node5#sss:kinsol> KIN_LSOLVE_FAIL *)
exception LinearSolverFailure

(** @kinsol <node5#sss:kinsol> KIN_SYSFUNC_FAIL *)
exception SystemFunctionFailure

(** @kinsol <node5#sss:kinsol> KIN_FIRST_SYSFUNC_FAIL *)
exception FirstSystemFunctionFailure

(** @kinsol <node5#sss:kinsol> KIN_REPTD_SYSFUNC_ERR *)
exception RepeatedSystemFunctionFailure

(** This type represents a session with the KINSOL solver.

    A skeleton of the main program:
    + {b Set vector with initial guess}
    {[let u = Nvector_array.wrap [| 0.0; 0.0; 0.0 |] ]}
    The length of this vector determines the problem size.    
    + {b Create and initialize a solver session}
    {[let s = Kinsol.init (Kinsol.Spgmr callbacks) f u]}
    This will also initialize the specified linear solver.
    + {b Set optional inputs}, e.g.
    {[set_num_max_iters s 500; ...]}
    Call any of the [set_*] functions to change solver and linear solver
    parameters from their defaults.
    + {b Solve problem}
    {[let result = Kinsol.solve s u u_scale f_scale]}
    + {b Get optional outputs}
    {[let fnorm = Kinsol.get_func_norm s in ...]}
    Call any of the [get_*] functions to examine solver and linear solver
    statistics.

    @kinsol <node5#ss:skeleton_sol> Skeleton of main program *)
type ('data, 'kind) session = ('data, 'kind) Kinsol_impl.session

type real_array = Sundials.RealArray.t
type serial_session = (real_array, Nvector_serial.kind) session

(** The type of vectors passed to the solver. *)
type ('data, 'kind) nvector = ('data, 'kind) Sundials.nvector

(** {2 Linear Solvers} *)

(**
   Specify a linear solver.

   The Lapack solvers require that both Sundials and the OCaml interface were
   built to link with a LAPACK library.

   @kinsol <node5#sss:lin_solv_init> Linear Solver Specification Functions *)
type ('data, 'kind) linear_solver = ('data, 'kind) Kinsol_impl.linear_solver
type serial_linear_solver = (real_array, Nvector_serial.kind) linear_solver

type 'a single_tmp = 'a
type 'a double_tmp = 'a * 'a

(**
  Arguments common to all Jacobian callback functions.    
 
  @kinsol <node5#ss:djacFn> Dense Jacobian function
  @kinsol <node5#ss:bjacFn> Banded Jacobian function
  @kinsol <node5#ss:psolveFn> Linear preconditioning function
  @kinsol <node5#ss:precondFn> Jacobian preconditioning function *)
type ('t, 'a) jacobian_arg =
  {
    jac_u   : 'a;   (** The current (unscaled) iterate. *)
    jac_fu  : 'a;   (** The current value of the vector F([u]). *)
    jac_tmp : 't    (** Workspace data: {!single_tmp} or {!double_tmp}. *)
  }

(** The range of nonzero entries in a band matrix.  *)
type bandrange = Kinsol_impl.bandrange =
  { mupper : int; (** The upper half-bandwidth.  *)
    mlower : int; (** The lower half-bandwidth.  *) }

module Dls :
  sig
    (** Direct Linear Solvers (DLS) operating on dense and banded matrices.
        
        @kinsol <node5#sss:optin_dls> Direct linear solvers optional input functions
        @kinsol <node5#sss:optout_dls> Direct linear solvers optional output functions
        @kinsol <node5#ss:djacFn> Dense Jacobian function *)

    (** The type of a user-supplied callback function that computes an
        approximation to the Jacobian matrix for the [Dense] and [LapackDense]
        {!linear_solver}s.

        This callback is invoked as [dense_jac_fn arg jac], where
        - [arg] is the {!jacobian_arg} with two work vectors.
        - [jac] is the matrix in which to store the computed Jacobian.

        Only nonzero elements need to be loaded into [jac] because [jac] is set
        to the zero matrix before the call to the Jacobian function.

        {b NB:} The elements of both arguments to this function must no longer
        be accessed after the callback has returned, i.e. if their values are
        needed outside of the function call, then they must be copied to
        separate physical structures.

        @kinsol <node5#sss:lin_solv_init> KINDense/KINLapackDense
        @kinsol <node5#sss:optin_dls> KINDlsSetDenseJacFn
        @kinsol <node5#ss:djacFn> KINDlsDenseJacFn *)
    type dense_jac_fn = (real_array double_tmp, real_array) jacobian_arg
                                                -> Dls.DenseMatrix.t -> unit

    (** Direct linear solver with dense matrix.  The optional argument specifies
        a callback function that computes an approximation to the Jacobian
        matrix (see {!dense_jac_fn} for details).  If this argument is [None],
        then KINSOL uses a default implementation based on difference quotients.

        @kinsol <node5#sss:lin_solv_init> KINDense
        @kinsol <node5#sss:optin_dls> KINDlsSetDenseJacFn
        @kinsol <node5#ss:djacFn> Dense Jacobian function *)
    val dense : dense_jac_fn option -> serial_linear_solver

    (** Direct linear solver with dense matrix, using LAPACK.  The argument is
        the same as [Dense].

        @kinsol <node5#sss:lin_solv_init> KINLapackDense
        @kinsol <node5#sss:optin_dls> KINDlsSetDenseJacFn
        @kinsol <node5#ss:djacFn> Dense Jacobian function *)
    val lapack_dense : dense_jac_fn option -> serial_linear_solver

    (** A user-supplied callback function that computes an approximation to
        the Jacobian matrix for the Band and Lapackband {!linear_solver}s.

        The callback is invoked as [band_jac_fn arg mupper mlower jac] where:
        - [arg] is the {!jacobian_arg} with two work vectors.
        - [mupper] is the upper half-bandwidth of the Jacobian.
        - [mlower] is the lower half-bandwidth of the Jacobian.
        - [jac] is the matrix in which to store the computed Jacobian.

        Only nonzero elements need to be loaded into [jac] because [jac] is set to
        the zero matrix before the call to the Jacobian function.

        {b NB:} The elements of [arg] and [jac] must no longer be accessed after
        this function has returned.  If their values are needed outside of the
        function call, then they must be copied to separate physical structures.

        @kinsol <node5#sss:lin_solv_init> KINBand/KINLapackBand
        @kinsol <node5#sss:optin_dls> KINDlsSetBandJacFn
        @kinsol <node5#ss:bjacFn> KINDlsBandJacFn *)
    type band_jac_fn = bandrange
                        -> (real_array double_tmp, real_array) jacobian_arg
                        -> Dls.BandMatrix.t -> unit

    (** Direct linear solver with banded matrix.  The arguments specify the
        width of the band ({!bandrange}) and an optional Jacobian function
        ({!band_jac_fn}).  If the Jacobian function is [None], KINSOL uses an
        internal implementation based on difference quotients.

        @kinsol <node5#sss:lin_solv_init> KINBand
        @kinsol <node5#sss:optin_dls> KINDlsSetBandJacFn
        @kinsol <node5#ss:bjacFn> Banded Jacobian function *)
    val band : bandrange -> band_jac_fn option -> serial_linear_solver

    (** Direct linear solver with banded matrix using LAPACK.  The arguments
        are the same as [Band].

        @kinsol <node5#sss:lin_solv_init> KINLapackBand
        @kinsol <node5#sss:optin_dls> KINDlsSetBandJacFn
        @kinsol <node5#ss:bjacFn> Banded Jacobian function *)
    val lapack_band : bandrange -> band_jac_fn option -> serial_linear_solver

    (** {4 Low-level solver manipulation} *)

    (** Change the dense Jacobian function (see [Dense] in {!linear_solver}). *)
    val set_dense_jac_fn : serial_session -> dense_jac_fn -> unit

    (** Remove the user-supplied dense Jacobian function, if any, and fall back
        to KINSOLS's internal different quotient approximation (see [Dense] in
        {!linear_solver}). *)
    val clear_dense_jac_fn : serial_session -> unit

    (** Change the band Jacobian function (see [Band] in {!linear_solver}). *)
    val set_band_jac_fn : serial_session -> band_jac_fn -> unit

    (** Remove the user-supplied band Jacobian function, if any, and fall back
        to KINSOL's internal different qutoient implementation (see [Band] in
        {!linear_solver}). *)
    val clear_band_jac_fn : serial_session -> unit

    (** {4 Optional output functions} *)

    (** Returns the sizes of the real and integer workspaces used by the Dense
        and Band direct linear solvers .

        @kinsol <node5#sss:optout_dense> KINDlsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : serial_session -> int * int

    (** Returns the number of calls made to the Dense and Band direct linear
        solvers Jacobian approximation function.

        @kinsol <node5#sss:optout_dense> KINDlsGetNumJacEvals *)
    val get_num_jac_evals : serial_session -> int

    (** Returns the number of calls to the user system function used to compute
        the difference quotient approximation to the dense or banded Jacobian.

        @kinsol <node5#sss:optout_dense> KINDlsGetNumFuncEvals *)
    val get_num_func_evals : serial_session -> int
  end

module Spils :
  sig
    (** Scaled Preconditioned Iterative Linear Solvers (SPILS)

        @kinsol <node5#sss:optin_spils> Iterative linear solvers optional input functions.
        @kinsol <node5#sss:optout_spils> Iterative linear solvers optional output functions.
        @kinsol <node5#ss:psolveFn> Linear preconditioning function
        @kinsol <node5#ss:precondFn> Jacobian preconditioning function *)

    (** Arguments passed to the preconditioner solve callback function.  See
        {!prec_solve_fn}. *)
    type 'a solve_arg =
      {
        uscale : 'a;  (** A vector containing diagonal elements of the
                          scaling matrix for [u] *)
        fscale : 'a;  (** A vector containing diagonal elements of the
                          scaling matrix for [fval]. *)
      }

    (** Called like [prec_solve_fn jarg sarg v] to solve the
        preconditioning system {i P}[z] = [r] where {i P} is the
        preconditioning matrix chosen by the user.  {i P} should
        approximate the system Jacobian [J] (partial [F] / partial [u]).
        - [jarg] supplies the basic problem data as a {!jacobian_arg}.
        - [sarg] specifies the linear system as a {!solve_arg}.
        - [v] is initially the right-hand side vector [r]. On completion, it
        must contain the solution [z].

        The {!Sundials.RecoverableFailure} exception can be raised to
        indicate a recoverable failure. Any other exception indicates
        an unrecoverable failure.

        {b NB:} The fields of [jarg], [sarg], and [z] must no longer
        be accessed after [psolve] has returned a result, i.e. if
        their values are needed outside of the function call, then
        they must be copied to separate physical structures.

        See also {!preconditioner}.

        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
        @kinsol <node5#ss:psolveFn> Linear preconditioning function
      *)
    and 'a prec_solve_fn =
      ('a single_tmp, 'a) jacobian_arg
      -> 'a solve_arg
      -> 'a
      -> unit

    (** A function that preprocesses and/or evaluates any Jacobian-related data
        needed by [prec_solve_fn] above.  When [prec_solve_fn] doesn't need any
        such data, this field can be [None].

        This callback is invoked as [prec_setup_fn jarg sarg], where
        - [jarg] supplies the basic problem data as a {!jacobian_arg}.
        - [sarg] specifies the linear system as a {!solve_arg}.
        It should return whether it was successful or not.

        {b NB:} The fields of [arg] must no longer be accessed after this
        callback has returned, i.e. if their values are needed outside of the
        function call, then they must be copied to separate physical
        structures.

        See also {!preconditioner}.

        @kinsol <node5#ss:precondFn> Jacobian preconditioning function
        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
     *)
    and 'a prec_setup_fn =
      ('a double_tmp, 'a) jacobian_arg
      -> 'a solve_arg
      -> unit

    (** This callback is invoked as [new_u' = jac_times_vec_fn v jv u new_u],
        to compute the Jacobian at [u] times a vector (or an approximation
        thereof), i.e., [jv = J(u) v]. That is,
        - [v] is the vector to multiply from the left,
        - [jv] is where the result is to be stored,
        - [u] is the current value of the dependent variable vector, and,

        The flag [new_u] indicates whether [u] has been updated since the last
        callback, and thus whether and saved Jacobian data should be recomputed.
        The callback routine can return false updated value for this flag, [new_u'],
        allowing it to reuse computed information at the next call
        ([new_u' = false]).

        {b NB:} The elements of [v], [jv], and [u] must no longer be accessed
        after this callback has returned, i.e. if their values are needed
        outside of the function call, then they must be copied to separate
        physical structures.

        See also {!preconditioner}.

        @kinsol <node5#ss:jtimesFn> KINSpilsJacTimesVecFn
        @kinsol <node5#sss:optin_spils> KINSpilsSetJacTimesVecFn
     *)
    and 'a jac_times_vec_fn =
      'a      (* v *)
      -> 'a   (* jv *)
      -> 'a   (* u *)
      -> bool (* new_u *)
      -> bool

    (** Specifies a preconditioner, including the type of
        preconditioning to be done (none or right), and a set of
        three callbacks if applicable:

        - [solve], the main function that solves the preconditioning
          system $Pz = r$, where $P$ is a preconditioning matrix
          chosen by the user.  See {!prec_solve_fn} for details.
        - [setup], which preprocesses and/or evaluates
          Jacobian-related data needed by [solve].  It can be omitted
          if there are no such data.  See {!prec_setup_fn} for
          details.
        - [jac_times_vec], which multiplies the system Jacobian to a
          given vector.  See {!jac_times_vec_fn} for details.  If the
          user doesn't give such a function, KINSOL uses a default
          implementation based on difference quotients.

        The following convenience functions are provided for
        constructing values of this type concisely:

        - {!prec_none} is just [PrecNone].
        - {!prec_right} creates [PrecRight] but takes optional fields
          as optional arguments: e.g. [prec_right ~setup:setup solve]
          returns [PrecRight (solve, setup, None)].

        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
        @kinsol <node5#sss:optin_spils> KINSpilsSetJacTimesVecFn
        @kinsol <node5#ss:psolveFn> KINSpilsPrecSolveFn
        @kinsol <node5#ss:precondFn> KINSpilsPrecSetupFn
        @kinsol <node5#ss:jtimesFn> KINSpilsJacTimesVecFn
    *)
    type 'a preconditioner =
      | PrecNone
      | PrecRight of 'a prec_solve_fn
                     * 'a prec_setup_fn option
                     * 'a jac_times_vec_fn option

    (** See {!preconditioner}.  *)
    val prec_none : 'a preconditioner

    (** See {!preconditioner}.  *)
    val prec_right :
      ?setup:'a prec_setup_fn
      -> ?jac_times_vec:'a jac_times_vec_fn
      -> 'a prec_solve_fn
      -> 'a preconditioner

    (** Krylov iterative solver with the scaled preconditioned GMRES method.
        Called like [spgmr ~maxl:maxl ~max_restarts:maxr prec], where:

        - [~maxl] is the maximum dimension of the Krylov subspace.
          Defaults to [5].
        - [~max_restarts] is the maximum number of restarts.  Defaults
          to [5].  Passing [0] disables restarts.
        - [prec] is a preconditioner.  See {!preconditioner}.

        @kinsol <node5#sss:lin_solv_init> KINSpgmr
        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
        @kinsol <node5#sss:optin_spils> KINSpilsSetMaxRestarts
        @kinsol <node5#ss:psolveFn> Linear preconditioning function
        @kinsol <node5#ss:precondFn> Jacobian preconditioning function *)
    val spgmr : ?maxl:int -> ?max_restarts:int -> 'a preconditioner
                    -> ('a, 'k) linear_solver

    (** Krylov iterative linear solver with the scaled preconditioned
        Bi-CGStab method.  The arguments are the same as {!spgmr},
        except the maximum number of restarts ([~max_restarts]) cannot
        be specified.

        @kinsol <node5#sss:lin_solv_init> KINSpbcg
        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
        @kinsol <node5#ss:psolveFn> Linear preconditioning function
        @kinsol <node5#ss:precondFn> Jacobian preconditioning function *)
    val spbcg : ?maxl:int -> 'a preconditioner -> ('a, 'k) linear_solver

    (** Krylov iterative linear solver with the scaled preconditioned
        TFQMR method.  The arguments are the same as {!spgmr}, except
        the maximum number of restarts ([~max_restarts]) cannot be
        specified.

        @kinsol <node5#sss:lin_solv_init> KINSptfqmr
        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
        @kinsol <node5#ss:psolveFn> Linear preconditioning function
        @kinsol <node5#ss:precondFn> Jacobian preconditioning function *)
    val sptfqmr : ?maxl:int -> 'a preconditioner -> ('a, 'k) linear_solver

    (** {4 Low-level solver manipulation} *)

    (** [set_preconditioner s psetup psolve] sets the preconditioning functions
        (see {!callbacks}).

        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
        @kinsol <node5#ss:precondFn> KINSpilsPrecSetupFn
        @kinsol <node5#ss:psolveFn> KINSpilsPrecSolveFn *)
    val set_preconditioner :
      ('a, 'k) session
      -> ?setup:'a prec_setup_fn
      -> 'a prec_solve_fn
      -> unit

    (** Set the Jacobian-times-vector function.

        @kinsol <node5#sss:optin_spils> KINSpilsSetJacTimesVecFn
        @kinsol <node5#ss:jtimesFn> KINSpilsJacTimesVecFn *)
    val set_jac_times_vec_fn :
      ('a, 'k) session
      -> 'a jac_times_vec_fn
      -> unit

    (** Use the default Jacobian-times-vector function.

        @kinsol <node5#sss:optin_spils> KINSpilsSetJacTimesVecFn *)
    val clear_jac_times_vec_fn : ('a, 'k) session -> unit

    (** {4 Optional input functions} *)

    (** Returns the sizes of the real and integer workspaces used by the linear
        solver.

        @kinsol <node5#sss:optout_spils> KINSpilsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space       : ('a, 'k) session -> int * int

    (** Returns the cumulative number of linear iterations.

        @kinsol <node5#sss:optout_spils> KINSpilsGetNumLinIters *)
    val get_num_lin_iters    : ('a, 'k) session -> int

    (** Returns the cumulative number of linear convergence failures.

        @kinsol <node5#sss:optout_spils> KINSpilsGetNumConvFails *)
    val get_num_conv_fails   : ('a, 'k) session -> int

    (** Returns the number of preconditioner evaluations, i.e., the number of
        calls made to psetup (see {!set_preconditioner}).

        @kinsol <node5#sss:optout_spils> KINSpilsGetNumPrecEvals *)
    val get_num_prec_evals   : ('a, 'k) session -> int

    (** Returns the cumulative number of calls made to the preconditioner solve
        function, psolve (see {!set_preconditioner}).

        @kinsol <node5#sss:optout_spils> KINSpilsGetNumPrecSolves *)
    val get_num_prec_solves  : ('a, 'k) session -> int

    (** Returns the cumulative number of calls made to the Jacobian-vector
        function, jtimes (see {! set_jac_times_vec_fn}).

        @kinsol <node5#sss:optout_spils> KINSpilsGetNumJtimesEvals *)
    val get_num_jtimes_evals : ('a, 'k) session -> int

    (** Returns the number of calls to the user system function for difference
        quotient Jacobian-vector product approximations.

        @kinsol <node5#sss:optout_spils> KINSpilsGetNumFuncEvals *)
    val get_num_func_evals    : ('a, 'k) session -> int
  end

module Alternate :
  sig
    (** Alternate Linear Solvers

        @kinsol <node8#s:new_linsolv> Providing Alternate Linear Solver Modules *)

    type ('data, 'kind) callbacks =
      {
        linit  : ('data, 'kind) linit option;
        lsetup : ('data, 'kind) lsetup option;
        lsolve : ('data, 'kind) lsolve;
      }

    (** Complete initializations for a specific linear solver, such as
        counters and statistics.

        Raising any exception in this function (including
        {!Sundials.RecoverableFailure}) is treated as an unrecoverable
        error.

        See also {!callbacks}.

        @cvode <node8#SECTION00810000000000000000> linit *)
    and ('data, 'kind) linit = ('data, 'kind) session -> unit

    (** The job of lsetup is to prepare the linear solver for subsequent
        calls to {!lsolve}. It may recompute Jacobian-related data if it
        deems necessary.

        Raising any exception in this function (including
        {!Sundials.RecoverableFailure}) is treated as an unrecoverable
        error.  Note this behavior is different from CVODE's lsetup.

        See also {!callbacks}.

        @kinsol <node8#SECTION00820000000000000000> lsetup *)
    and ('data, 'kind) lsetup = ('data, 'kind) session -> unit

    (** [res_norm = lsolve x b] must solve the linear equation given:
        - [x], on entry: an initial guess, on return: it should contain
          the solution to [Jx = b].
        - [b] is the right-hand side vector, set to [-F(u)], evaluated at
          the current iterate.

        This function should return the L2 norm of the residual
        vector. If an error occurs and recovery could be possible by
        calling again the {!lsetup} function, this function may raise
        a {!Sundials.RecoverableFailure} exception to indicate this
        fact.  Any other exception is treated as an unrecoverable
        error.

        See also {!callbacks}.

        @kinsol <node8#SECTION00830000000000000000> lsolve *)
    and ('data, 'kind) lsolve =
      ('data, 'kind) session
      -> 'data
      -> 'data
      -> float option

    (** Returns the internal [u] and [uscale] values. *)
    val get_u_uscale : ('data, 'kind) session -> 'data * 'data

    (** Returns the internal [f] and [fscale] values. *)
    val get_f_fscale  : ('data, 'kind) session -> 'data * 'data

    (** Sets the internal [sJpnorm] value. *)
    val set_sjpnorm   : ('data, 'kind) session -> float -> unit

    (** Sets the internal [sfdotJp] value. *)
    val set_sfdotjp   : ('data, 'kind) session -> float -> unit

    (** Create a linear solver from a function returning a set of callback
        functions *)
    val make_solver :
          (('data, 'kind) session -> ('data, 'kind) nvector option
                                                  -> ('data, 'kind) callbacks)
          -> ('data, 'kind) linear_solver
  end

(** Increasing levels of verbosity for informational messages. *)
type print_level =
  | NoInformation     (** [0] no information displayed. *)
  | ShowScaledNorms   (** [1] at each nonlinear iteration, display the scaled
                          Euclidean norm of the system function evaluated at the
                          current iterate, the scaled norm of the Newton step
                          (if no globalization strategy is used), and the number
                          of function evaluations performed so far. *)
  | ShowScaledDFNorm  (** [2] additionally display [||F(u)||_DF] (if no
                          globalization strategy is used) and
                          [||F(u)||_DF,infty]. *)
  | ShowGlobalValues  (** [3] additionally display the values used by the global
                          strategy and statistical information for the linear
                          solver. *)

(** The parameters {i gamma} and {i alpha} in the formula for the Eisenstat and
    Walker Choice 2 for {i eta}. Set either to [None] to specify its default
    value. The legal values are 0 < [egamma] <= 1.0 and 1 < [ealpha] <= 2.0.

    @kinsol <node3#SECTION00300900000000000000>   Stopping criteria for iterative linear solvers *)
type eta_params = {
  egamma : float option; (** default = 0.9 *)
  ealpha : float option; (** default = 2.0 *)
}

(** The choice of {i eta} used in the stopping criteria for the linear system
    solver.

    @kinsol <node3#SECTION00300900000000000000>   Stopping criteria for iterative linear solvers *)
type eta_choice =
  | EtaChoice1                   (** Eisenstat and Walker Choice 1 *)
  | EtaChoice2 of eta_params     (** Eisenstat and Walker Choice 2 *)
  | EtaConstant of float option  (** Constant (default = 0.1) *)

(** {2 Initialization} *)

(** The system function of a nonlinear problem.

    See also {!init}.

    @kinsol <node5#ss:sysFn>           Problem-defining function
  *)
type 'a sysfn = 'a -> 'a -> unit

(** [init linsolv f tmpl] creates a KINSOL session, where
     - [linsolv] specify the linear solver to attach.
     - [f]       the system function of the nonlinear problem,
     - [tmpl]    used as a template to initialize the session (e.g., the
                 initial guess vector), and,

     The {!Sundials.RecoverableFailure} exception can be raised to indicate a
     recoverable failure. Any other exception indicates an unrecoverable
     failure.

     @kinsol <node5#sss:kinmalloc>     KINCreate/KINInit
     @kinsol <node5#sss:lin_solv_init> Linear solver specification functions *)
val init : ('a, 'kind) linear_solver -> 'a sysfn
           -> ('a, 'kind) nvector -> ('a, 'kind) session

(** {2 Solver functions } *)

type result =
  | Success           (** KIN_SUCCESS *)
  | InitialGuessOK    (** KIN_INITIAL_GUESS_OK *)
  | StoppedOnStepTol  (** KIN_STEP_LT_STPTOL *)

(** [r = solve s u linesearch u_scale f_scale] computes an approximate solution
    to the nonlinear system, where
    - [s] is a session created with {!init},
    - [u] is the initial guess, and, on return, the approximate solution
    for [F(u) = 0],
    - [linesearch] specifies whether to usea  globalization strategy (line
    search) or not,
    - [u_scale] contains the diagonal elements of the scaling matrix {i D_u} for
    vector [u] chosen so that the components of {i D_u}[u] all have roughly the
    same magnitude when [u] is close to a root of {i F(}[u]{i )},
    - [f_scale] contains the diagonal elements of the scaling matrix {i D_f} for
    {i F(}[u]{i )} chosen so that the components of {i D_f}{i F(}[u]{i )} all
    have roughtly the same magnitude when [u] is not too near a root of {i
    F(}[u]{i )}.

    On success, the result [r] indicates
    - [Success], that the scaled norm of {i F(}[u]{i )} is less than {i fnormtol}
    (see {!set_func_norm_tol}),
    - [InitialGuessOK], that the initial guess already satisfied the system
    within the tolerances given, or,
    - [StoppedOnStepTol], that the current iterate may be an approximate
    solution of the given nonlinear system, but it is also quite possible that
    the algorithm stalled near an invalid solution or that scsteptol is too
    large (see {!set_scaled_step_tol}).

    On failure, an exception is thrown.
 
    @kinsol <node5#sss:kinsol> KINSol
    @raise IllInput                     An input parameter was invalid.
    @raise LineSearchNonConvergence     The line search algorithm was unable to find an iterate sufficiently distinct from the current iterate, or could not find an iterate satisfying the sufficient decrease condition.
    @raise MaxIterationsReached         The maximum number of nonlinear iterations was reached.
    @raise MaxNewtonStepExceeded        Five consecutive steps have been taken that satisfy a scaled step length test.
    @raise LineSearchBetaConditionFailure   The line search algorithm was unable to satisfy the beta-condition for [nbcfails] iterations.
    @raise LinearSolverNoRecovery       The linear solver's solve function failed recoverably, but the Jacobian data is already current.
    @raise LinearSolverInitFailure      The linear solver's init routine failed.
    @raise LinearSetupFailure           The linear solver's setup function failed in an unrecoverable manner.
    @raise LinearSolverFailure          The linear solver's solve function failed in an unrecoverable manner.
    @raise SystemFunctionFailure        The system function failed in an unrecoverable manner.
    @raise FirstSystemFunctionFailure   The system function failed at the first call.
    @raise RepeatedSystemFunctionFailure   Unable to correct repeated recoverable system function errors. *)
val solve :
    ('a, 'k) session
    -> ('a, 'k) nvector
    -> bool
    -> ('a, 'k) nvector
    -> ('a, 'k) nvector
    -> result

(** {2 Main optional functions} *)

(** {3 Input} *)

(** [set_error_file s fname trunc] opens the file named [fname] and to which all
    messages from the default error handler are then directed.
    If the file already exists it is either trunctated ([trunc] = [true]) or
    appended to ([trunc] = [false]).

    The error file is closed if {!set_error_file} is called again, or otherwise
    when the session is garbage collected.

    @kinsol <node5#ss:optin_main> KINSetErrFile *)
val set_error_file : ('a, 'k) session -> string -> bool -> unit

(** [set_err_handler_fn s efun] specifies a custom function [efun] for
    handling error messages.  The error handler function must not fail
    -- any exceptions raised from it will be captured and discarded.

    @kinsol <node5#ss:optin_main> KINSetErrHandlerFn
    @kinsol <node5#ss:ehFn> KINErrHandlerFn *)
val set_err_handler_fn : ('a, 'k) session -> (Sundials.error_details -> unit) -> unit

(** This function restores the default error handling function. It is equivalent
    to calling KINSetErrHandlerFn with an argument of [NULL].

    @kinsol <node5#ss:optin_main> KINSetErrHandlerFn *)
val clear_err_handler_fn : ('a, 'k) session -> unit

(** [set_info_file s fname trunc] opens the file named [fname] to which all
    informative (non-error) messages are then written. If the file already
    exists it is either trunctated ([trunc] = [true]) or appended to ([trunc] =
    [false]).

    The info file is closed if {!set_info_file} is called again, or otherwise
    when the session is garbage collected.
   
    @kinsol <node5#ss:optin_main> KINSetInfoFile *)
val set_info_file : ('a, 'k) session -> string -> bool -> unit

(** [set_info_handler_fn s efun] specifies a custom function [efun] for handling
    informative (non-error) messages. The [error_code] field of the
    {!Sundials.error_details} value is set to zero for informative messages.

    @kinsol <node5#ss:optin_main> KINSetInfoHandlerFn
    @kinsol <node5#ss:ihFn> KINInfoHandlerFn *)
val set_info_handler_fn : ('a, 'k) session -> (Sundials.error_details -> unit) -> unit

(** This function restores the default info handling function. It is equivalent
    to calling KINSetInfoHandlerFn with an argument of [NULL].

    @kinsol <node5#ss:optin_main> KINSetErrHandlerFn *)
val clear_info_handler_fn : ('a, 'k) session -> unit

(** Sets the level of verbosity of the output (see {!print_level}).

    @kinsol <node5#ss:optin_main> KINSetPrintLevel *)
val set_print_level : ('a, 'k) session -> print_level -> unit

(** Specify the maximum number of nonlinear iterations allowed (defaults to
    200).

    @kinsol <node5#ss:optin_main> KINSetNumMaxIters *)
val set_num_max_iters : ('a, 'k) session -> int -> unit

(** Specifies wether an initial call to the preconditioner setup function is
    made ([false], the default) or not made ([true]).

    This function is useful when solving a sequence of problems where the final
    preconditioner values of a problem becomes the initial value for the next
    problem.

    @kinsol <node5#ss:optin_main> KINSetNoInitSetup *)
val set_no_init_setup : ('a, 'k) session -> bool -> unit

(** Controls whether the nonlinear residual monitoring scheme is used ([false],
    the default), or not ([true]) to control Jacobian updating. It only has an
    effect for the Dense and Band solvers.

    @kinsol <node5#ss:optin_main> KINSetNoResMon *)
val set_no_res_mon : serial_session -> bool -> unit

(** Specifies the maximum number of nonlinear iterations between calls to the
    preconditioner setup function. Pass [None] to set the default (10).

    @kinsol <node5#ss:optin_main> KINSetMaxSetupCalls *)
val set_max_setup_calls : ('a, 'k) session -> int option -> unit

(** Specifies the maximum number of nonlinear iterations between checks by the
    residual monitoring algorithm. Pass [None] to set the default (5). It only
    has an effect for the Dense and Band solvers.

    @kinsol <node5#ss:optin_main> KINSetMaxSubSetupCalls *)
val set_max_sub_setup_calls : serial_session -> int option -> unit

(** Specifies the method for computing the value of the {i eta} coefficient used
    in the calculation of the linear solver convergence tolerance. (See
    {!eta_choice}; the default is [EtaChoice1])

    @kinsol <node5#ss:optin_main> KINSetEtaForm
    @kinsol <node5#ss:optin_main> KINSetEtaConstValue
    @kinsol <node5#ss:optin_main> KINSetEtaParams *)
val set_eta_choice : ('a, 'k) session -> eta_choice -> unit

(** Specifies the constant value of {i omega} when using residual monitoring.
    Pass [None] to specify the default value (0.9). The legal values are 0 <
    [omegaconst] < 1.0.

    @kinsol <node5#ss:optin_main> KINSetResMonConstValue *)
val set_res_mon_const_value : ('a, 'k) session -> float option -> unit

(** [set_res_mon_params omegamin omegamax] specifies the parameters in the
    formula for {i omega}. Pass [None] to specify the default values (0.00001
    and 0.9, respectively).
    The legal values are 0 < [omegamin] < [omegamax] < 1.0.

    @kinsol <node5#ss:optin_main> KINSetResMonParams
    @kinsol <node3#SECTION00300800000000000000> Residual monitoring for Modified Newton method (2.3) *)
val set_res_mon_params : ('a, 'k) session -> float option -> float option -> unit

(** Specifies a flag that controls whether or not the value of {i epsilon}, the
    scaled linear residual tolerance, is bounded from below.

    The default value for this flag is [false] meaning that a positive minimum
    value equal to [0.01*fnormtol] is applied to {i epsilon}.

    @kinsol <node5#ss:optin_main> KINSetNoMinEps
    @kinsol <node5#ss:optin_main> KINSetFuncNormTol *)
val set_no_min_eps : ('a, 'k) session -> bool -> unit

(** Specifies the maximum allowable scaled length of the Newton step. Pass
    [None] to specify the default value (1000*||u_0||_Du), otherwise the given
    value must be greater than zero.

    @kinsol <node5#ss:optin_main> KINSetMaxNewtonStep *)
val set_max_newton_step : ('a, 'k) session -> float option -> unit

(** Specifies the maximum number of {i beta}-condition failures in the
   linesearch algorithm. Pass [None] to specify the default (10).

    @kinsol <node5#ss:optin_main> KINSetMaxBetaFails *)
val set_max_beta_fails : ('a, 'k) session -> float option -> unit

(** Specifies the relative error in computing {i F(u)}, which is used in the
    difference quotient approximation of the Jacobian-vector product. Pass
    [None] to specify the default value (the square root of unit roundoff).

    @kinsol <node5#ss:optin_main> KINSetRelErrFunc *)
val set_rel_err_func : ('a, 'k) session -> float option -> unit

(** Specifies the stopping tolerance on the scaled maximum norm of the system
    function {i F(u)}, which must be greater than zero ([fnormtol]). Pass [None]
    to specify the default value (unit roundoff^1/3).

    @kinsol <node5#ss:optin_main> KINSetFuncNormTol *)
val set_func_norm_tol : ('a, 'k) session -> float option -> unit

(** Specifies the stopping tolerance on the minimum scaled step length, which
    must be greater than zeor ([scsteptol]). Pass [None] to specify the default
    value (unit roundoff^1/3).

    @kinsol <node5#ss:optin_main> KINSetScaledStepTol *)
val set_scaled_step_tol : ('a, 'k) session -> float option -> unit

(** Specifies a vecotr that defines inequality constraints for each component of
    the solution vector [u]. For each value:
    -  0.0: no constraint is imposed on u_i,
    -  1.0: constrain u_i >= 0.0,
    - -1.0: constrain u_i <= 0.0,
    -  2.0: constrain u_i > 0.0, and,
    - -2.0: constrain u_i < 0.0.

    @kinsol <node5#ss:optin_main> KINSetConstraints *)
val set_constraints : ('a, 'k) session -> ('a, 'k) nvector -> unit

(** Allows one to change the linear solver so as to try solve a problem using
    different tools or parameters.

    @kinsol <node5#sss:lin_solv_init> Linear solver specification functions *)
val set_linear_solver : ('a, 'k) session -> ('a, 'k) linear_solver -> unit

(** Allows one to change the system function so as to solve several problems of
    the same size but with different functions.

    @kinsol <node5#ss:optin_main> KINSetSysFunc
    @kinsol <node5#ss:sysFn> Problem-defining function *)
val set_sys_func : ('a, 'k) session -> ('a -> 'a -> unit) -> unit

(** {3 Output} *)

(** [lenrw, leniw = get_work_space s] returns the number of realtype ([lenrw])
    and integer ([leniw]) values in the workspace.

    @kinsol <node5#sss:output_main> KINGetWorkSpace *)
val get_work_space : ('a, 'k) session -> int * int

(** Returns the number of evaluations of the system function.

    @kinsol <node5#ss:optout_main> KINGetNumFuncEvals *)
val get_num_func_evals : ('a, 'k) session -> int

(** Returns the number of nonlinear iterations.

    @kinsol <node5#ss:optout_main> KINGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('a, 'k) session -> int

(** Returns the number of {i beta}-condition failures.

    @kinsol <node5#ss:optout_main> KINGetNumBetaCondFails *)
val get_num_beta_cond_fails : ('a, 'k) session -> int

(** Returns the number of backtrack operations (step length adjustments)
    performed by the linesearch algorithm.

    @kinsol <node5#ss:optout_main> KINGetNumBacktrackOps *)
val get_num_backtrack_ops : ('a, 'k) session -> int

(** Returns the scaled Euclidiean {i l2} norm of the nonlinear system function
    {i F(u)} evaluated at the current iterate.

    @kinsol <node5#ss:optout_main> KINGetFuncNorm *)
val get_func_norm : ('a, 'k) session -> float

(** Returns the scaled Euclidiean {i l2} norm of the step used during the
    previous iteration.

    @kinsol <node5#ss:optout_main> KINGetStepLength *)
val get_step_length : ('a, 'k) session -> float

