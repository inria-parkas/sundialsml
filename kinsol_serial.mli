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
(*               User Documentation for KINSOL v2.7.0                  *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(* TODO: cvode = CVSpilsSetMaxl, only for CVSPBCG and CVSPTFQMR. *)
(* TODO: cvode = CVSpilsSetGSType, only for CVSPGMR. *)

(** Serial nvector interface to the KINSOL Solver.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

(**
    This type represents a session with the KINSOL solver using
    {!Nvector.nvector}s.

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
    {[let result = Kinsol.solve s u u_scale f_scale}
    + {b Get optional outputs}
    {[let fnorm = Kinsol.get_func_norm s in ...]}
    Call any of the [get_*] functions to examine solver and linear solver
    statistics.

    @kinsol <node5#ss:skeleton_sol> Skeleton of main program
 *)
type session

(** The type of vectors passed to the solver. *)
type nvec = Sundials.Carray.t

(** {2 Linear Solvers} *)

type single_tmp = nvec
type double_tmp = nvec * nvec

(**
  Arguments common to all Jacobian callback functions.    
 
  @kinsol <node5#ss:djacFn> Dense Jacobian function
  @kinsol <node5#ss:bjacFn> Banded Jacobian function
  @kinsol <node5#ss:psolveFn> Linear preconditioning function
  @kinsol <node5#ss:precondFn> Jacobian preconditioning function
 *)
type 't jacobian_arg =
  {
    jac_u   : nvec;   (** The current (unscaled) iterate. *)
    jac_fu  : nvec;   (** The current value of the vector F([u]). *)
    jac_tmp : 't      (** Workspace data,
                          either {!single_tmp} or {!double_tmp}. *)
  }

(** {3 Setting linear solvers to sessions} *)

(**
   Specify a linear solver.

   The Lapack solvers require that both Sundials and the OCaml interface were
   built to link with a LAPACK library.

   The Banded Krylov solvers imply an additional call to
   {{:KINSOL_DOC_ROOT(node5#sss:cvbandpre)} CVBandPrecInit}.

   @kinsol <node5#sss:lin_solv_init> Linear Solver Specification Functions
*)
type linear_solver =
  | Dense of dense_jac_fn option
  (** Direct linear solver with dense matrix.  The optional argument specifies
      a callback function that computes an approximation to the Jacobian matrix
      (see {!dense_jac_fn} for details).  If this argument is [None], then
      KINSOL uses a default implementation based on difference quotients.  See
      also {!Dls}.

      @kinsol <node5#sss:lin_solve_init> CVDense
      @kinsol <node5#sss:optin_dls> CVDlsSetDenseJacFn
      @kinsol <node5#ss:djacFn> Dense Jacobian function
  *)
  | LapackDense of dense_jac_fn option
  (** Direct linear solver with dense matrix, using LAPACK.  The argument is
      the same as [Dense].  See also {!Dls}.

      @kinsol <node5#sss:lin_solve_init> CVLapackDense
      @kinsol <node5#sss:optin_dls> CVDlsSetDenseJacFn
      @kinsol <node5#ss:djacFn> Dense Jacobian function
   *)
  | Band of bandrange * band_jac_fn option
  (** Direct linear solver with banded matrix.  The arguments specify the width
      of the band ({!bandrange}) and an optional Jacobian function
      ({!band_jac_fn}).  If the Jacobian function is [None], KINSOL uses an
      internal implementation based on difference quotients.  See also {!Dls}.

      @kinsol <node5#sss:lin_solve_init> CVBand
      @kinsol <node5#sss:optin_dls> CVDlsSetBandJacFn
      @kinsol <node5#ss:bjacFn> Banded Jacobian function
   *)
  | LapackBand of bandrange * band_jac_fn option
  (** Direct linear solver with banded matrix using LAPACK.  The arguments
      are the same as [Band].

      @kinsol <node5#sss:lin_solve_init> CVLapackBand
      @kinsol <node5#sss:optin_dls> CVDlsSetBandJacFn
      @kinsol <node5#ss:bjacFn> Banded Jacobian function
   *)
  | Spgmr of int option * spils_callbacks
  (** Krylov iterative solver with the scaled preconditioned GMRES method.
      The arguments specify the maximum dimension of the Krylov subspace (Pass
      None to use the default value 5) and the preconditioner callback functions
      ({!spils_callbacks}).  See also {!Spils}.

      @kinsol <node5#sss:lin_solve_init> KINSpgmr
      @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
      @kinsol <node5#sss:optin_spils> KINSpilsSetMaxRestarts
      @kinsol <node5#ss:psolveFn> Linear preconditioning function
      @kinsol <node5#ss:precondFn> Jacobian preconditioning function
    *)
  | Spbcg of int option * spils_callbacks
  (** Krylov iterative solver with the scaled preconditioned Bi-CGStab method.
      The arguments specify the maximum dimension of the Krylov subspace (Pass
      None to use the default value 5) and the preconditioner callback functions
      ({!spils_callbacks}).  See also {!Spils}.

      @kinsol <node5#sss:lin_solve_init> KINSpbcg
      @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
      @kinsol <node5#ss:psolveFn> Linear preconditioning function
      @kinsol <node5#ss:precondFn> Jacobian preconditioning function
    *)
  | Sptfqmr of int option * spils_callbacks
  (** Krylov iterative with the scaled preconditioned TFQMR method.
      The arguments specify the maximum dimension of the Krylov subspace (Pass
      None to use the default value 5) and the preconditioner callback functions
      ({!spils_callbacks}).  See also {!Spils}.

      @kinsol <node5#sss:lin_solve_init> KINSptfqmr
      @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
      @kinsol <node5#ss:psolveFn> Linear preconditioning function
      @kinsol <node5#ss:precondFn> Jacobian preconditioning function
    *)

(** The type of a user-supplied callback function that computes an
    approximation to the Jacobian matrix for the [Dense] and [LapackDense]
    {!linear_solver}s.

    This callback is invoked as [dense_jac_fn arg jac], where
    - [arg] is the {!jacobian_arg} with two work vectors.
    - [jac] is the matrix in which to store the computed Jacobian.

    Only nonzero elements need to be loaded into [jac] because [jac] is set to
    the zero matrix before the call to the Jacobian function.

    {b NB:} The elements of both arguments to this function must no longer be
    accessed after the callback has returned, i.e. if their values are needed
    outside of the function call, then they must be copied to separate physical
    structures.

    @kinsol <node5#sss:lin_solve_init> KINDense/KINLapackDense
    @kinsol <node5#sss:optin_dls> KINDlsSetDenseJacFn
    @kinsol <node5#ss:djacFn> KINDlsDenseJacFn
 *)
and dense_jac_fn = double_tmp jacobian_arg -> Dls.DenseMatrix.t -> unit

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

    @kinsol <node5#sss:lin_solve_init> KINBand/KINLapackBand
    @kinsol <node5#sss:optin_dls> KINDlsSetBandJacFn
    @kinsol <node5#ss:bjacFn> KINDlsBandJacFn
 *)
and band_jac_fn = double_tmp jacobian_arg
                      -> int -> int -> Dls.BandMatrix.t -> unit

(** The range of nonzero entries in a band matrix.  *)
and bandrange = { mupper : int; (** The upper half-bandwidth.  *)
                  mlower : int; (** The lower half-bandwidth.  *) }

(** Callbacks for Krylov subspace linear solvers.  Ignored if the
    {!Spils.preconditioning_type} is set to [PrecNone].  In that case, you
    should use {!spils_no_precond} as [spils_callbacks].  *)
and spils_callbacks =
  {
    prec_solve_fn : (single_tmp jacobian_arg -> prec_solve_arg
                      -> nvec-> unit) option;
    (**
        This callback is invoked as [prec_solve_fn jarg sarg v] to solve
        the preconditioning system [P z = r]. The preconditioning matrix, [P],
        approximates the system Jacobian [J] (partial [F] / partial [u]).
        - [jarg] supplies the basic problem data as a {!jacobian_arg}.
        - [sarg] specifies the linear system as a {!prec_solve_arg}.
        - [v] is initially the right-hand side vector [r]. On completion, it
        must contain the solution [z].

        The {!Sundials.RecoverableFailure} exception can be raised to indicate a
        recoverable failure. Any other exception indicates an unrecoverable
        failure.

        {b NB:} The fields of [jarg], [sarg], and [z] must no longer be accessed
        after [psolve] has returned a result, i.e. if their values are needed
        outside of the function call, then they must be copied to separate
        physical structures.

        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
        @kinsol <node5#ss:psolveFn> Linear preconditioning function
    *)

    prec_setup_fn : (double_tmp jacobian_arg -> prec_solve_arg
                     -> unit) option;
    (** A function that preprocesses and/or evaluates any Jacobian-related data
        needed by [prec_solve_fn] above.  When [prec_solve_fn] doesn't need any
        such data, this field can be [None].

        This callback is invoked as [prec_setup_fn jarg sarg], where
        - [jarg] supplies the basic problem data as a {!jacobian_arg}.
        - [sarg] specifies the linear system as a {!prec_solve_arg}.
        It should return whether it was successful or not.

        {b NB:} The fields of [arg] must no longer be accessed after this
        callback has returned, i.e. if their values are needed outside of the
        function call, then they must be copied to separate physical
        structures.

        @kinsol <node5#ss:precondFn> Jacobian preconditioning function
        @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
     *)

    jac_times_vec_fn :
      (nvec      (* v *)
       -> nvec   (* jv *)
       -> nvec   (* u *)
       -> bool   (* new_u *)
       -> bool) option;
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

        @kinsol <node5#ss:jtimesFn> KINSpilsJacTimesVecFn
        @kinsol <node5#sss:optin_spils> KINSpilsSetJacTimesVecFn
     *)
  }

(** Arguments passed to the preconditioner solve callback function.  See
    [prec_solve_fn] in {!spils_callbacks}. 
 *)
and prec_solve_arg =
  {
    uscale : nvec;  (** A vector containing diagonal elements of the
                        scaling matrix for [u] *)
    fscale : nvec;  (** A vector containing diagonal elements of the
                        scaling matrix for [fval]. *)
  }

(** See {!spils_callbacks}. *)
val spils_no_precond : spils_callbacks

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

    @kinsol <node3#SECTION00300900000000000000>   Stopping criteria for
                                                  iterative linear solvers
 *)
type eta_params = {
  egamma : float option; (** default = 0.9 *)
  ealpha : float option; (** default = 2.0 *)
}

(** The choice of {i eta} used in the stopping criteria for the linear system
    solver.

    @kinsol <node3#SECTION00300900000000000000>   Stopping criteria for
                                                  iterative linear solvers
  *)
type eta_choice =
  | EtaChoice1                   (** Eisenstat and Walker Choice 1 *)
  | EtaChoice2 of eta_params     (** Eisenstat and Walker Choice 2 *)
  | EtaConstant of float option  (** Constant (default = 0.1) *)

(** {3 Direct Linear Solvers (DLS)} *)

(** Control callbacks and get optional outputs for the Direct Linear Solvers
    that operate on dense and banded matrices.
    
    @kinsol <node5#sss:optin_dls> Direct linear solvers optional input functions
    @kinsol <node5#sss:optout_dls> Direct linear solvers optional output functions
    @kinsol <node5#ss:djacFn> Dense Jacobian function
  *)
module Dls :
  sig
    (** {4 Low-level solver manipulation} *)

    (** Change the dense Jacobian function (see [Dense] in {!linear_solver}). *)
    val set_dense_jac_fn : session -> dense_jac_fn -> unit

    (** Remove the user-supplied dense Jacobian function, if any, and fall back
        to KINSOLS's internal different quotient approximation (see [Dense] in
        {!linear_solver}). *)
    val clear_dense_jac_fn : session -> unit

    (** Change the band Jacobian function (see [Band] in {!linear_solver}). *)
    val set_band_jac_fn : session -> band_jac_fn -> unit

    (** Remove the user-supplied band Jacobian function, if any, and fall back
        to KINSOL's internal different qutoient implementation (see [Band] in
        {!linear_solver}). *)
    val clear_band_jac_fn : session -> unit

    (** {4 Optional output functions} *)

    (**
      Returns the sizes of the real and integer workspaces used by the Dense and
      Band direct linear solvers .

      @kinsol <node5#sss:optout_dense> KINDlsGetWorkSpace
      @return ([real_size], [integer_size])
     *)
    val get_work_space : session -> int * int

    (**
      Returns the number of calls made to the Dense and Band direct linear
      solvers Jacobian approximation function.

      @kinsol <node5#sss:optout_dense> KINDlsGetNumJacEvals
    *)
    val get_num_jac_evals : session -> int

    (**
      Returns the number of calls to the user system function used to compute
      the difference quotient approximation to the dense or banded Jacobian.

      @kinsol <node5#sss:optout_dense> KINDlsGetNumFuncEvals
    *)
    val get_num_func_evals : session -> int
  end

(** {3 Scaled Preconditioned Iterative Linear Solvers (SPILS)} *)

(** Set callback functions, set optional outputs, and get optional inputs for
    the Scaled Preconditioned Iterative Linear Solvers: SPGMR, SPBCG, SPTFQMR.
    @kinsol <node5#sss:optin_spils> Iterative linear solvers optional input functions.
    @kinsol <node5#sss:optout_spils> Iterative linear solvers optional output functions.
    @kinsol <node5#ss:psolveFn> Linear preconditioning function
    @kinsol <node5#ss:precondFn> Jacobian preconditioning function
 *)
module Spils :
  sig
    (** {4 Low-level solver manipulation} *)

    (** Alias for {!prec_solve_arg}.  *)
    type solve_arg = prec_solve_arg =
      {
        uscale : nvec;
        fscale : nvec;
      }

    (**
      Set preconditioning functions (see {!spils_callbacks}).

      @kinsol <node5#sss:optin_spils> KINSpilsSetPreconditioner
      @kinsol <node5#ss:precondFn> Jacobian preconditioning function
      @kinsol <node5#ss:psolveFn> Linear preconditioning function
    *)
    val set_preconditioner :
      session
      -> (double_tmp jacobian_arg -> prec_solve_arg -> unit) option
      -> (single_tmp jacobian_arg -> prec_solve_arg -> nvec -> unit)
      -> unit

    (* TODO: add a clear_preconditioner function? *)

    (**
      Set the Jacobian-times-vector function.

      @kinsol <node5#sss:optin_spils> KINSpilsSetJacTimesVecFn
      @kinsol <node5#ss:jtimesFn> KINSpilsJacTimesVecFn
    *)
    val set_jac_times_vec_fn :
      session
      -> (nvec     (* v *)
          -> nvec  (* jv *)
          -> nvec  (* u *)
          -> bool  (* new_u *)
          -> bool)
      -> unit

    (**
      Use the default Jacobian-times-vector function.

      @kinsol <node5#sss:optin_spils> KINSpilsSetJacTimesVecFn
    *)
    val clear_jac_times_vec_fn : session -> unit

    (** {4 Optional input functions} *)

    (**
      Returns the sizes of the real and integer workspaces used by the linear
      solver.

      @kinsol <node5#sss:optout_spils> KINSpilsGetWorkSpace
      @return ([real_size], [integer_size])
    *)
    val get_work_space       : session -> int * int

    (**
      Returns the cumulative number of linear iterations.

      @kinsol <node5#sss:optout_spils> KINSpilsGetNumLinIters
    *)
    val get_num_lin_iters    : session -> int

    (**
      Returns the cumulative number of linear convergence failures.

      @kinsol <node5#sss:optout_spils> KINSpilsGetNumConvFails
    *)
    val get_num_conv_fails   : session -> int

    (**
      Returns the number of preconditioner evaluations, i.e., the number of
      calls made to psetup (see {!set_preconditioner}).

      @kinsol <node5#sss:optout_spils> KINSpilsGetNumPrecEvals
    *)
    val get_num_prec_evals   : session -> int

    (**
      Returns the cumulative number of calls made to the preconditioner solve
      function, psolve (see {!set_preconditioner}).

      @kinsol <node5#sss:optout_spils> KINSpilsGetNumPrecSolves
    *)
    val get_num_prec_solves  : session -> int

    (**
      Returns the cumulative number of calls made to the Jacobian-vector
      function, jtimes (see {! set_jac_times_vec_fn}).

      @kinsol <node5#sss:optout_spils> KINSpilsGetNumJtimesEvals
    *)
    val get_num_jtimes_evals : session -> int

    (**
      Returns the number of calls to the user system function for difference
      quotient Jacobian-vector product approximations.

      @kinsol <node5#sss:optout_spils> KINSpilsGetNumFuncEvals
    *)
    val get_num_func_evals    : session -> int
  end

(** {2 Initialization} *)

(**
   [init linsolv f tmpl] creates a KINSOL session, where
   - [linsolv] specify the linear solver to attach.
   - [f]       the system function of the nonlinear problem,
   - [tmpl]    used as a template to initialize the session (e.g., the
               initial guess vector), and,

   The {!Sundials.RecoverableFailure} exception can be raised to indicate a
   recoverable failure. Any other exception indicates an unrecoverable failure.

   @kinsol <node5#sss:kinmalloc>      KINCreate/KINInit
   @kinsol <node5#ss:sysFn>           Problem-defining function
   @kinsol <node5#sss:lin_solve_init> Linear solver specification functions
 *)
val init : linear_solver -> (nvec -> nvec -> unit) -> nvec -> session

(** {2 Solver functions } *)

type result =
  | Success           (** KIN_SUCCESS *)
  | InitialGuessOK    (** KIN_INITIAL_GUESS_OK *)
  | StoppedOnStepTol  (** KIN_STEP_LT_STPTOL *)

(**
    [r = solve s u linesearch u_scale f_scale] computes an approximate solution
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
    @raise LineSearchNonConvergence     The line search algorithm was unable
                                        to find an iterate sufficiently distinct
                                        from the current iterate, or could not
                                        find an iterate satisfying the
                                        sufficient decrease condition.
    @raise MaxIterationsReached         The maximum number of nonlinear
                                        iterations was reached.
    @raise MaxNewtonStepExceeded        Five consecutive steps have been taken
                                        that satisfy a scaled step length test.
    @raise LineSearchBetaConditionFailure   The line search algorithm was
                                        unable to satisfy the beta-condition
                                        for [nbcfails] iterations.
    @raise LinearSolverNoRecovery       The linear solver's solve function
                                        failed recoverably, but the Jacobian
                                        data is already current.
    @raise LinearSolverInitFailure      The linear solver's init routine failed.
    @raise LinearSetupFailure           The linear solver's setup function failed
                                        in an unrecoverable manner.
    @raise LinearSolverFailure          The linear solver's solve function failed
                                        in an unrecoverable manner.
    @raise SystemFunctionFailure        The system function failed in an
                                        unrecoverable manner.
    @raise FirstSystemFunctionFailure   The system function failed at the first
                                        call.
    @raise RepeatedSystemFunctionFailure   Unable to correct repeated
                                        recoverable system function errors.
 *)
val solve : session -> nvec -> bool -> nvec -> nvec -> result

(** {2 Main optional functions} *)

(** {3 Input} *)

(**
  [set_error_file s fname trunc] opens the file named [fname] and to which all
  messages from the default error handler are then directed.
  If the file already exists it is either trunctated ([trunc] = [true]) or
  appended to ([trunc] = [false]).

  The error file is closed if {!set_error_file} is called again, or otherwise
  when the session is garbage collected.
   
  @kinsol <node5#ss:optin_main> KINSetErrFile
 *)
val set_error_file : session -> string -> bool -> unit

(**
  [set_err_handler_fn s efun] specifies a custom function [efun] for handling
  error messages.

  @kinsol <node5#ss:optin_main> KINSetErrHandlerFn
  @kinsol <node5#ss:ehFn> KINErrHandlerFn
 *)
val set_err_handler_fn : session -> (Sundials.error_details -> unit) -> unit

(**
  This function restores the default error handling function. It is equivalent
  to calling KINSetErrHandlerFn with an argument of [NULL].

  @kinsol <node5#ss:optin_main> KINSetErrHandlerFn
 *)
val clear_err_handler_fn : session -> unit

(**
  [set_info_file s fname trunc] opens the file named [fname] to which all
  informative (non-error) messages are then written.
  If the file already exists it is either trunctated ([trunc] = [true]) or
  appended to ([trunc] = [false]).

  The info file is closed if {!set_info_file} is called again, or otherwise when
  the session is garbage collected.
   
  @kinsol <node5#ss:optin_main> KINSetInfoFile
 *)
val set_info_file : session -> string -> bool -> unit

(**
  [set_info_handler_fn s efun] specifies a custom function [efun] for handling
  informative (non-error) messages. The [error_code] field of the
  {!Sundials.error_details} value is set to zero for informative messages.

  @kinsol <node5#ss:optin_main> KINSetInfoHandlerFn
  @kinsol <node5#ss:ihFn> KINInfoHandlerFn
 *)
val set_info_handler_fn : session -> (Sundials.error_details -> unit) -> unit

(**
  This function restores the default info handling function. It is equivalent
  to calling KINSetInfoHandlerFn with an argument of [NULL].

  @kinsol <node5#ss:optin_main> KINSetErrHandlerFn
 *)
val clear_info_handler_fn : session -> unit

(**
  Sets the level of verbosity of the output (see {!print_level}).

  @kinsol <node5#ss:optin_main> KINSetPrintLevel
 *)
val set_print_level : session -> print_level -> unit

(**
  Specify the maximum number of nonlinear iterations allowed (defaults to 200).

  @kinsol <node5#ss:optin_main> KINSetNumMaxIters
 *)
val set_num_max_iters : session -> int -> unit

(**
  Specifies wether an initial call to the preconditioner setup function is made
  ([false], the default) or not made ([true]).

  This function is useful when solving a sequence of problems where the final
  preconditioner values of a problem becomes the initial value for the next
  problem.

  @kinsol <node5#ss:optin_main> KINSetNoInitSetup
 *)
val set_no_init_setup : session -> bool -> unit

(**
  Controls whether the nonlinear residual monitoring scheme is used ([false],
  the default), or not ([true]) to control Jacobian updating. It only has an
  effect for the Dense and Band solvers.

  @kinsol <node5#ss:optin_main> KINSetNoResMon
 *)
val set_no_res_mon : session -> bool -> unit

(**
  Specifies the maximum number of nonlinear iterations between calls to the
  preconditioner setup function. Pass [None] to set the default (10).

  @kinsol <node5#ss:optin_main> KINSetMaxSetupCalls
 *)
val set_max_setup_calls : session -> int option -> unit

(**
  Specifies the maximum number of nonlinear iterations between checks by the
  residual monitoring algorithm. Pass [None] to set the default (5). It only has
  an effect for the Dense and Band solvers.

  @kinsol <node5#ss:optin_main> KINSetMaxSubSetupCalls
 *)
val set_max_sub_setup_calls : session -> int option -> unit

(**
  Specifies the method for computing the value of the {i eta} coefficient used
  in the calculation of the linear solver convergence tolerance. (See
  {!eta_choice}; the default is [EtaChoice1])

  @kinsol <node5#ss:optin_main> KINSetEtaForm
  @kinsol <node5#ss:optin_main> KINSetEtaConstValue
  @kinsol <node5#ss:optin_main> KINSetEtaParams
 *)
val set_eta_choice : session -> eta_choice -> unit

(**
  Specifies the constant value of {i omega} when using residual monitoring.
  Pass [None] to specify the default value (0.9). The legal values are
  0 < [omegaconst] < 1.0.

  @kinsol <node5#ss:optin_main> KINSetResMonConstValue
 *)
val set_res_mon_const_value : session -> float option -> unit

(**
  [set_res_mon_params omegamin omegamax] specifies the parameters in the formula
  for {i omega}. Pass [None] to specify the default values (0.00001 and 0.9,
  respectively). The legal values are 0 < [omegamin] < [omegamax] < 1.0.

  @kinsol <node5#ss:optin_main> KINSetResMonParams
  @kinsol <node3#SECTION00300800000000000000> Residual monitoring for Modified
                                              Newton method (2.3)
 *)
val set_res_mon_params : session -> float option -> float option -> unit

(**
  Specifies a flag that controls whether or not the value of {i epsilon}, the
  scaled linear residual tolerance, is bounded from below.

  The default value for this flag is [false] meaning that a positive minimum
  value equal to [0.01*fnormtol] is applied to {i epsilon}.

  @kinsol <node5#ss:optin_main> KINSetNoMinEps
  @kinsol <node5#ss:optin_main> KINSetFuncNormTol
 *)
val set_no_min_eps : session -> bool -> unit

(**
  Specifies the maximum allowable scaled length of the Newton step. Pass [None]
  to specify the default value (1000*||u_0||_Du), otherwise the given value must
  be greater than zero.

  @kinsol <node5#ss:optin_main> KINSetMaxNewtonStep
 *)
val set_max_newton_step : session -> float option -> unit

(**
  Specifies the maximum number of {i beta}-condition failures in the linesearch
  algorithm. Pass [None] to specify the default (10).

  @kinsol <node5#ss:optin_main> KINSetMaxBetaFails
 *)
val set_max_beta_fails : session -> float option -> unit

(**
  Specifies the relative error in computing {i F(u)}, which is used in the
  difference quotient approximation of the Jacobian-vector product. Pass [None]
  to specify the default value (the square root of unit roundoff).

  @kinsol <node5#ss:optin_main> KINSetRelErrFunc
 *)
val set_rel_err_func : session -> float option -> unit

(**
  Specifies the stopping tolerance on the scaled maximum norm of the system
  function {i F(u)}, which must be greater than zero ([fnormtol]). Pass [None]
  to specify the default value (unit roundoff^1/3).

  @kinsol <node5#ss:optin_main> KINSetFuncNormTol
 *)
val set_func_norm_tol : session -> float option -> unit

(**
  Specifies the stopping tolerance on the minimum scaled step length, which must
  be greater than zeor ([scsteptol]). Pass [None] to specify the default value
  (unit roundoff^1/3).

  @kinsol <node5#ss:optin_main> KINSetScaledStepTol
 *)
val set_scaled_step_tol : session -> float option -> unit

(**
  Specifies a vecotr that defines inequality constraints for each component of
  the solution vector [u]. For each value:
  -  0.0: no constraint is imposed on u_i,
  -  1.0: constrain u_i >= 0.0,
  - -1.0: constrain u_i <= 0.0,
  -  2.0: constrain u_i > 0.0, and,
  - -2.0: constrain u_i < 0.0.

  @kinsol <node5#ss:optin_main> KINSetConstraints
 *)
val set_constraints : session -> nvec -> unit

(**
  Allows one to change the system function so as to solve several problems of
  the same size but with different functions.

  @kinsol <node5#ss:optin_main> KINSetSysFunc
  @kinsol <node5#ss:sysFn> Problem-defining function
 *)
val set_sys_func : session -> (nvec -> nvec -> unit) -> unit

(** {3 Output} *)

(**
  [lenrw, leniw = get_work_space s] returns the number of realtype ([lenrw]) and
  integer ([leniw]) values in the workspace.

  @kinsol <node5#sss:output_main> KINGetWorkSpace
 *)
val get_work_space : session -> int * int

(**
  Returns the number of evaluations of the system function.

  @kinsol <node5#ss:optout_main> KINGetNumFuncEvals
 *)
val get_num_func_evals : session -> int

(**
  Returns the number of nonlinear iterations.

  @kinsol <node5#ss:optout_main> KINGetNumNonlinSolvIters
 *)
val get_num_nonlin_solv_iters : session -> int

(**
  Returns the number of {i beta}-condition failures.

  @kinsol <node5#ss:optout_main> KINGetNumBetaCondFails
 *)
val get_num_beta_cond_fails : session -> int

(**
  Returns the number of backtrack operations (step length adjustments) performed
  by the linesearch algorithm.

  @kinsol <node5#ss:optout_main> KINGetNumBacktrackOps
 *)
val get_num_backtrack_ops : session -> int

(**
  Returns the scaled Euclidiean {i l2} norm of the nonlinear system function
  {i F(u)} evaluated at the current iterate.

  @kinsol <node5#ss:optout_main> KINGetFuncNorm
 *)
val get_func_norm : session -> float

(**
  Returns the scaled Euclidiean {i l2} norm of the step used during the previous
  iteration.

  @kinsol <node5#ss:optout_main> KINGetStepLength
 *)
val get_step_length : session -> float

