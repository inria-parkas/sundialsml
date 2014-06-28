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

(** Interface to the CVODE solver.

 @version VERSION()
 @author Timothy Bourke (Inria)
 @author Jun Inoue (Inria)
 @author Marc Pouzet (LIENS)
 *)

(** {2:exceptions Exceptions} *)

(** @cvode <node5#sss:cvode> CV_ILL_INPUT *)
exception IllInput

(** @cvode <node5#sss:cvode> CV_TOO_CLOSE *)
exception TooClose

(** @cvode <node5#sss:cvode> CV_TOO_MUCH_WORK *)
exception TooMuchWork

(** @cvode <node5#sss:cvode> CV_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** @cvode <node5#sss:cvode> CV_ERR_FAIL *)
exception ErrFailure                

(** @cvode <node5#sss:cvode> CV_CONV_FAIL *)
exception ConvergenceFailure        

(** @cvode <node5#sss:cvode> CV_LINIT_FAIL *)
exception LinearInitFailure         

(** @cvode <node5#sss:cvode> CV_LSETUP_FAIL *)
exception LinearSetupFailure        

(** @cvode <node5#sss:cvode> CV_LSOLVE_FAIL *)
exception LinearSolveFailure        

(** @cvode <node5#sss:cvode> CV_RHSFUNC_FAIL *)
exception RhsFuncFailure

(** @cvode <node5#sss:cvode> CV_FIRST_RHSFUNC_ERR *)
exception FirstRhsFuncErr

(** @cvode <node5#sss:cvode> CV_REPTD_RHSFUNC_ERR *)
exception RepeatedRhsFuncErr        

(** @cvode <node5#sss:cvode> CV_UNREC_RHSFUNC_ERR *)
exception UnrecoverableRhsFuncErr   

(** @cvode <node5#sss:cvode> CV_RTFUNC_FAIL *)
exception RootFuncFailure           

exception BadK      (** k is not in the range 0, 1, ..., q_u (CV_BAD_K)
                        @cvode <node5#ss:optional_dky> CVodeGetDky *)

exception BadT      (** t is not in the interval
                        \[t_n - h_u, t_n\] (CV_BAD_T)
                        @cvode <node5#ss:optional_dky> CVodeGetDky *)
exception BadDky    (** invalid dky argument (CV_BAD_DKY)
                        @cvode <node5#ss:optional_dky> CVodeGetDky *)

(**
    This type represents a session with the CVODE solver using
    {!Nvector.nvector}s.

    A skeleton of the main program:
    + {b Set vector of initial values}
    {[let y = Nvector_array.wrap [| 0.0; 0.0; 0.0 |] ]}
    The length of this vector determines the problem size.    
    + {b Create and initialize a solver session}
    {[let s = init Adams Functional tols f ~roots:(2, g) y]}
    This will initialize a specific linear solver and the root-finding
    mechanism, if necessary.
    + {b Specify integration tolerances (optional)}, e.g.
    {[set_tolerances s SStolerances (reltol, abstol)]}
    + {b Set optional inputs}, e.g.
    {[set_stop_time s 10.0; ...]}
    Call any of the [set_*] functions to change solver parameters from their
    defaults.
    + {b Advance solution in time}, e.g.
    {[let (t', result) = solve_normal s !t y in
...
t := t' + 0.1]}
    Repeatedly call either [solve_normal] or [solve_one_step] to advance the
    simulation.
    + {b Get optional outputs}
    {[let stats = get_integrator_stats s in ...]}
    Call any of the [get_*] functions to examine solver statistics.

    @cvode <node5#ss:skeleton_sim> Skeleton of main program
 *)
type ('a, 'k) session = ('a, 'k) Cvode_session.session

type real_array = Sundials.RealArray.t
type serial_session = (real_array, Nvector_serial.kind) session

(** The type of vectors passed to the solver. *)
type ('data, 'kind) nvector = ('data, 'kind) Sundials.nvector

(** {2 Linear Solvers} *)

(** Specify a linear solver.

    @cvode <node5#sss:lin_solv_init> Linear Solver Specification Functions *)
type ('data, 'kind) linear_solver
type serial_linear_solver = (real_array, Nvector_serial.kind) linear_solver

module Diag :
  sig
    (** Diagonal approximation of the Jacobian by difference quotients.

        @cvode <node5#sss:lin_solve_init> CVDiag *)
    val solver : ('data, 'kind) linear_solver

    (** {4 Optional output functions} *)

    (** Get optional inputs for the linear solver that gives diagonal
        approximations of the Jacobian matrix.
        
        @cvode <node5#sss:optout_diag> Diagonal linear solver optional output functions *)

    (** Returns the sizes of the real and integer workspaces used by the Diagonal
        linear solver.

        @cvode <node5#sss:optout_diag> CVDiagGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : ('a, 'k) session -> int * int

    (** Returns the number of calls made to the user-supplied right-hand side
        function due to finite difference Jacobian approximation in the Diagonal
        linear solver.

        @cvode <node5#sss:optout_diag> CVDiagGetNumRhsEvals *)
    val get_num_rhs_evals : ('a, 'k) session -> int
  end

type 'a single_tmp = 'a
type 'a triple_tmp = 'a * 'a * 'a

(** Arguments common to all Jacobian callback functions.    
 
    @cvode <node5#ss:jtimesFn> Jacobian-times-vector function
    @cvode <node5#ss:psolveFn> Linear preconditioning function
    @cvode <node5#ss:precondFn> Jacobian preconditioning function *)
type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;        (** The independent variable. *)
    jac_y   : 'a;           (** The dependent variable vector. *)
    jac_fy  : 'a;           (** The derivative vector (i.e., f(t, y)). *)
    jac_tmp : 't            (** Workspace data,
                                either {!single_tmp} or {!triple_tmp}. *)
  }

(** The range of nonzero entries in a band matrix.  *)
type bandrange = { mupper : int; (** The upper half-bandwidth.  *)
                   mlower : int; (** The lower half-bandwidth.  *) }

module Dls :
  sig
    (** Direct Linear Solvers operating on dense and banded matrices.
        
        @cvode <node5#sss:optin_dls> Direct linear solvers optional input functions
        @cvode <node5#sss:optout_dls> Direct linear solvers optional output functions
        @cvode <node5#ss:djacFn> Dense Jacobian function
      *)

    (** The type of a user-supplied callback function that computes an
        approximation to the Jacobian matrix for the [Dense] and [LapackDense]
        {!linear_solver}s.

        The function is called like [dense_jac_fn arg jac] where:
        - [arg] is the standard {!jacobian_arg} with three work vectors.
        - [jac] is the matrix in which to store the computed Jacobian.
        The function should load the ({i i,j}) entry of the Jacobian with {i
        dFi/dyj}, i.e. the partial derivative of the right-hand side of the {i
        i}-th equation with respect to the {i j}-th variable, evaluated at the
        values of ({i t,y}) that can be obtained from [arg].

        Only nonzero elements need to be loaded into [jac] because [jac] is set
        to the zero matrix before the call to the Jacobian function.

        The function may raise a {!Sundials.RecoverableFailure} exception to
        indicate that a recoverable error has occurred. Any other exception is
        treated as an unrecoverable error.

        {b NB:} The elements of both arguments to this function must no longer
        be accessed after the callback has returned, i.e. if their values are
        needed outside of the function call, then they must be copied to
        separate physical structures.

        @cvode <node5#sss:lin_solve_init> CVDense
        @cvode <node5#sss:optin_dls> CVDlsSetDenseJacFn
        @cvode <node5#ss:djacFn> Dense Jacobian function *)
    type dense_jac_fn = (real_array triple_tmp, real_array) jacobian_arg
                                                  -> Dls.DenseMatrix.t -> unit

    (** Direct linear solver with dense matrix.  The optional argument specifies
        a callback function that computes an approximation to the Jacobian
        matrix (see {!dense_jac_fn} for details).  If this argument is [None],
        then CVODE uses a default implementation based on difference quotients.
        See also {!Dls}.

        @cvode <node5#sss:lin_solve_init> CVDense
        @cvode <node5#sss:optin_dls> CVDlsSetDenseJacFn
        @cvode <node5#ss:djacFn> Dense Jacobian function *)
    val dense : dense_jac_fn option -> serial_linear_solver

    (** Direct linear solver with dense matrix, using LAPACK.  The argument is
        the same as [Dense].  See also {!Dls}.

        Both Sundials and the OCaml interface must have been built to link with
        a LAPACK library.

        @cvode <node5#sss:lin_solve_init> CVLapackDense
        @cvode <node5#sss:optin_dls> CVDlsSetDenseJacFn
        @cvode <node5#ss:djacFn> Dense Jacobian function *)
    val lapack_dense : dense_jac_fn option -> serial_linear_solver

    (** A user-supplied callback function that computes an approximation to
        the Jacobian matrix for the Band and Lapackband {!linear_solver}s.  If
        this field is [None], CVODE uses a default implementation based on
        difference quotients.

        The function is called as [band_jac_fn {mupper; mlower} arg jac] where:
        - [mupper] is the upper half-bandwidth of the Jacobian.
        - [mlower] is the lower half-bandwidth of the Jacobian.
        - [arg] is the standard {!jacobian_arg} with three work vectors.
        - [jac] is the matrix in which to store the computed Jacobian.
        The function should load the ({i i,j}) entry of the Jacobian with {i
        dFi/dyj}, i.e. the partial derivative of the right-hand side of the {i
        i}-th equation with respect to the {i j}-th variable, evaluated at the
        values of ({i t,y}) that can be obtained from [arg].

        Only nonzero elements need to be loaded into [jac] because [jac] is set
        to the zero matrix before the call to the Jacobian function.

        The function may raise a {!Sundials.RecoverableFailure} exception to
        indicate that a recoverable error has occurred. Any other exception is
        treated as an unrecoverable error.

        {b NB:} The elements of [arg] and [jac] must no longer be accessed after
        this function has returned.  If their values are needed outside of the
        function call, then they must be copied to separate physical
        structures. *)
    type band_jac_fn = bandrange
                        -> (real_array triple_tmp, real_array) jacobian_arg
                        -> Dls.BandMatrix.t -> unit

    (** Direct linear solver with banded matrix.  The arguments specify the
        width of the band ({!bandrange}) and an optional Jacobian function
        ({!band_jac_fn}).  If the Jacobian function is [None], CVODE uses an
        internal implementation based on difference quotients.

        @cvode <node5#sss:lin_solve_init> CVBand
        @cvode <node5#sss:optin_dls> CVDlsSetBandJacFn
        @cvode <node5#ss:bjacFn> Banded Jacobian function *)
    val band : bandrange -> band_jac_fn option -> serial_linear_solver

    (** Direct linear solver with banded matrix using LAPACK.  The arguments
        are the same as [Band].

        Both Sundials and the OCaml interface must have been built to link with
        a LAPACK library.

        @cvode <node5#sss:lin_solve_init> CVLapackBand
        @cvode <node5#sss:optin_dls> CVDlsSetBandJacFn
        @cvode <node5#ss:bjacFn> Banded Jacobian function *)
    val lapack_band : bandrange -> band_jac_fn option -> serial_linear_solver

    (** {4 Low-level solver manipulation} *)

    (** Change the dense Jacobian function (see [Dense] in {!linear_solver}).
        It may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.  *)
    val set_dense_jac_fn : serial_session -> dense_jac_fn -> unit

    (** Remove the user-supplied dense Jacobian function, if any, and fall back
        to IDA's internal implementation (see [Dense] in {!linear_solver}).  It
        may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.  *)
    val clear_dense_jac_fn : serial_session -> unit

    (** Change the band Jacobian function (see [Band] in {!linear_solver}).
        It may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.  *)
    val set_band_jac_fn : serial_session
                                                        -> band_jac_fn -> unit

    (** Remove the user-supplied band Jacobian function, if any, and fall back
        to IDA's internal implementation (see [Band] in {!linear_solver}).  It
        may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.  *)
    val clear_band_jac_fn : serial_session -> unit

    (** {4 Optional output functions} *)

    (** Returns the sizes of the real and integer workspaces used by the Dense
        and Band direct linear solvers .

        @cvode <node5#sss:optout_dls> CVDlsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : serial_session -> int * int


    (** Returns the number of calls made to the Dense and Band direct linear
        solvers Jacobian approximation function.

        @cvode <node5#sss:optout_dls> CVDlsGetNumJacEvals *)
    val get_num_jac_evals : serial_session -> int

    (** Returns the number of calls made to the user-supplied right-hand side
        function due to the finite difference (Dense or Band) Jacobian
        approximation.

        @cvode <node5#sss:optout_dls> CVDlsGetNumRhsEvals *)
    val get_num_rhs_evals : serial_session -> int
  end

module Spils :
  sig
    (** Scaled Preconditioned Iterative Linear Solvers (SPILS)

        @cvode <node5#sss:optin_spils> Iterative linear solvers optional input functions.
        @cvode <node5#sss:optout_spils> Iterative linear solvers optional output functions.
        @cvode <node5#ss:psolveFn> Linear preconditioning function
        @cvode <node5#ss:precondFn> Jacobian preconditioning function *)

    type gramschmidt_type = Spils.gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    type preconditioning_type = Spils.preconditioning_type =
      | PrecNone
      | PrecLeft
      | PrecRight
      | PrecBoth

    (** Arguments passed to the preconditioner solve callback function.  See
        [prec_solve_fn] in {!callbacks}.

        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn *)
    type 'a solve_arg =
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

    (** Callbacks for Krylov subspace linear solvers.  Ignored if the
        {!preconditioning_type} is set to [PrecNone].  In that case, you
        should use {!no_precond} as [callbacks].  *)
    type 'a callbacks =
      {
        prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a solve_arg
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
            @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn *)

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
            @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)

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
            @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn *)
      }

    val no_precond : 'a callbacks

    (** Krylov iterative solver with the scaled preconditioned GMRES method.
        The arguments specify the maximum dimension of the Krylov subspace (Pass
        [None] to use the default value [5].), preconditioning type, and the
        preconditioner callback functions ({!callbacks}). See also {!Spils}.

        @cvode <node5#sss:lin_solve_init> CVSpgmr
        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
    val spgmr : int option -> preconditioning_type -> 'data callbacks
                    -> ('data, 'kind) linear_solver

    (** Krylov iterative solver with the scaled preconditioned Bi-CGStab method.
        The arguments are the same as [Spgmr].  See also {!Spils}.

        @cvode <node5#sss:lin_solve_init> CVSpbcg
        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
    val spbcg : int option -> preconditioning_type -> 'data callbacks
                    -> ('data, 'kind) linear_solver

    (** Krylov iterative with the scaled preconditioned TFQMR method.  The
        arguments are the same as [Spgmr].  See also {!Spils}.

        @cvode <node5#sss:lin_solve_init> CVSptfqmr
        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
    val sptfqmr : int option -> preconditioning_type -> 'data callbacks
                    -> ('data, 'kind) linear_solver

    (** {4 Low-level solver manipulation} *)

    (** Set preconditioning functions (see {!callbacks}).  It may be
        unsafe to use this function without a {!reinit}.  Users are encouraged
        to use the [iter_type] parameter of {!reinit} instead, unless they are
        desperate for performance.

        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
    val set_preconditioner :
      ('a, 'k) session
      -> (('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool)
      -> (('a single_tmp, 'a) jacobian_arg -> 'a solve_arg -> 'a -> unit)
      -> unit

    (** Set the Jacobian-times-vector function.  It may be unsafe to use this
        function without a {!reinit}.  Users are encouraged to use the
        [iter_type] parameter of {!reinit} instead, unless they are desperate
        for performance.

        @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
        @cvode <node5#ss:jtimesFn> Jacobian-times-vector function *)
    val set_jac_times_vec_fn :
      ('a, 'k) session
      -> (('a single_tmp, 'a) jacobian_arg
          -> 'a (* v *)
          -> 'a (* Jv *)
          -> unit)
      -> unit

    (** This function restores the default Jacobian-times-vector function. It is
        equivalent to calling CVodeSetJacTimesVecFn with an argument of [NULL].
        It may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead, unless
        they are desperate for performance.

        @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
        @cvode <node5#ss:jtimesFn> Jacobian-times-vector function *)
    val clear_jac_times_vec_fn : ('a, 'k) session -> unit

    (** {4 Optional input functions} *)

    (** This function resets the type of preconditioning to be used using a
        value of type {!preconditioning_type}.

        @cvode <node5#sss:optin_spils> CVSpilsPrecSetupFn *)
    val set_prec_type : ('a, 'k) session -> preconditioning_type -> unit

    (** Sets the Gram-Schmidt orthogonalization to be used with the
        Spgmr {!linear_solver}.

        @cvode <node5#sss:optin_spils> CVSpilsSetGSType *)
    val set_gs_type : ('a, 'k) session -> gramschmidt_type -> unit

    (** [set_eps_lin eplifac] sets the factor by which the Krylov linear
        solver's convergence test constant is reduced from the Newton iteration
        test constant. [eplifac]  must be >= 0. Passing a value of 0 specifies
        the default (which is 0.05).

        @cvode <node5#sss:optin_spils> CVSpilsSetEpsLin *)
    val set_eps_lin : ('a, 'k) session -> float -> unit

    (** [set_maxl maxl] resets the maximum Krylov subspace dimension for the
        Bi-CGStab or TFQMR methods. [maxl] is the maximum dimension of the Krylov
        subspace, a value of [maxl] <= 0 specifies the default (which is 5.0).

        @cvode <node5#sss:optin_spils> CVSpilsSetMaxl *)
    val set_maxl : ('a, 'k) session -> int -> unit

    (** {4 Optional output functions} *)

    (** Returns the sizes of the real and integer workspaces used by the SPGMR
        linear solver.

        @cvode <node5#sss:optout_spils> CVSpilsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space       : ('a, 'k) session -> int * int

    (** Returns the cumulative number of linear iterations.

        @cvode <node5#sss:optout_spils> CVSpilsGetNumLinIters *)
    val get_num_lin_iters    : ('a, 'k) session -> int

    (** Returns the cumulative number of linear convergence failures.

        @cvode <node5#sss:optout_spils> CVSpilsGetNumConvFails *)
    val get_num_conv_fails   : ('a, 'k) session -> int

    (** Returns the number of preconditioner evaluations, i.e., the number of
        calls made to psetup with jok = [false] (see {!set_preconditioner}).

        @cvode <node5#sss:optout_spils> CVSpilsGetNumPrecEvals *)
    val get_num_prec_evals   : ('a, 'k) session -> int

    (** Returns the cumulative number of calls made to the preconditioner solve
        function, psolve (see {!set_preconditioner}).

        @cvode <node5#sss:optout_spils> CVSpilsGetNumPrecSolves *)
    val get_num_prec_solves  : ('a, 'k) session -> int

    (** Returns the cumulative number of calls made to the Jacobian-vector
        function, jtimes (see {! set_jac_times_vec_fn}).

        @cvode <node5#sss:optout_spils> CVSpilsGetNumJtimesEvals *)
    val get_num_jtimes_evals : ('a, 'k) session -> int

    (** Returns the number of calls to the user right-hand side function for
        finite difference Jacobian-vector product approximation. This counter is
        only updated if the default difference quotient function is used.

        @cvode <node5#sss:optout_spils> CVSpilsGetNumRhsEvals *)
    val get_num_rhs_evals    : ('a, 'k) session -> int

    module Banded :
      sig
        (** Banded preconditioners

            @cvode <node5#sss:cvbandpre> Serial banded preconditioner module *)

        (** Same as Spgmr (the Krylov iterative solver with scaled preconditioned
            GMRES), but the preconditioner is set to CVODE's internal implementation
            using a banded matrix of difference quotients.  The arguments specify
            the maximum dimension of the Krylov subspace (Pass [None] to use the
            default value [5].), preconditioning type, and the width of the band
            matrix ({!bandrange}).

            @cvode <node5#sss:lin_solve_init> CVSpgmr
            @cvode <node5#sss:cvbandpre> CVBandPrecInit *)
        val spgmr : int option -> preconditioning_type -> bandrange
                            -> serial_linear_solver

        (** Same as Spbcg (the Krylov iterative solver with scaled preconditioned
            Bi-CGStab), but the preconditioner is set to CVODE's internal
            implementation using a banded matrix of difference quotients.  The
            arguments are the same as [BandedSpgmr].

            @cvode <node5#sss:lin_solve_init> CVSpbcg
            @cvode <node5#sss:cvbandpre> CVBandPrecInit *)
        val spbcg : int option -> preconditioning_type -> bandrange
                            -> serial_linear_solver

        (** Same as Spbcg (the Krylov iterative solver with scaled preconditioned
            Bi-CGStab), but the preconditioner is set to CVODE's internal
            implementation using a banded matrix of difference quotients.  The
            arguments are the same as [BandedSpgmr].

            @cvode <node5#sss:lin_solve_init> CVSptfqmr
            @cvode <node5#sss:cvbandpre> CVBandPrecInit *)
        val sptfqmr : int option -> preconditioning_type -> bandrange
                               -> serial_linear_solver

        (** {4 Optional output functions} *)

        (** Returns the sizes of the real and integer workspaces used by the serial
            banded preconditioner module.

            @cvode <node5#sss:cvbandpre> CVBandPrecGetWorkSpace
            @return ([real_size], [integer_size]) *)
        val get_work_space : serial_session -> int * int

        (** Returns the number of calls made to the user-supplied right-hand side
            function due to finite difference banded Jacobian approximation in the
            banded preconditioner setup function.

            @cvode <node5#sss:cvbandpre> CVBandPrecGetNumRhsEvals *)
        val get_num_rhs_evals : serial_session -> int
      end
  end

(** {2 Tolerances} *)

type ('data, 'kind) tolerance =
  | SStolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVtolerances of float * ('data, 'kind) nvector
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFtolerances of ('data -> 'data -> unit)
    (** Specifies a function [efun y ewt] that sets the multiplicative
        error weights Wi for use in the weighted RMS norm. The function is
        passed the dependent variable vector [y] and is expected to set the
        values inside the error-weight vector [ewt]. *)

(** A default relative tolerance of 1.0e-4 and absolute tolerance of 1.0e-8. *)
val default_tolerances : ('data, 'kind) tolerance

(** {2 Initialization} *)

(** Specify a solution method.

    @cvode <node3#ss:ivp_sol> IVP Solution
    @cvode <node5#sss:cvodemalloc> CVodeCreate *)
type ('a, 'kind) iter =
  | Newton of ('a, 'kind) linear_solver
    (** Newton iteration with a given linear solver *)
  | Functional
    (** Functional iteration (non-stiff systems only) *)

(** Specify a linear multistep method.

    @cvode <node3#ss:ivp_sol> IVP Solution
    @cvode <node5#sss:cvodemalloc> CVodeCreate *)
type lmm =
  | Adams   (** Non-stiff systems; Adams-Moulton formulas *)
  | BDF     (** Stiff systems;     Backward Differentiation Formulas *)

(**
    [init lmm iter tol f ~roots:(nroots, g) ~t0:t0 (neqs, y0)] initializes the
    CVODE solver and returns a {!session}.
    - [lmm]     specifies the linear multistep method, see {!lmm},
    - [iter]    specifies either functional iteration or Newton iteration
                with a specific linear solver, see {!iter},
    - [tol]     specifies the integration tolerances,
    - [f]       is the ODE right-hand side function,
    - [nroots]  specifies the number of root functions (zero-crossings),
    - [g]       calculates the values of the root functions,
    - [t0]      is the initial value of the independent variable, and,
    - [y0]      is a vector of initial values, the size of this vector
                determines the number of equations in the session, see
                {!Sundials.RealArray.t}.

    The labeled arguments [roots] and [t0] are both optional and default to
    {!no_roots} (i.e. no root finding is done) and [0.0], respectively.

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
    The {!no_roots} value can be passed for the [(nroots, g)] argument if
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
    @cvode <node5#ss:ewtsetFn>       Error weight function *)
val init :
    lmm
    -> ('a, 'kind) iter
    -> ('a, 'kind) tolerance
    -> (float -> 'a -> 'a -> unit)
    -> ?roots:(int * (float -> 'a -> Sundials.Roots.val_array -> unit))
    -> ?t0:float
    -> ('a, 'kind) nvector
    -> ('a, 'kind) session

(** This is a convenience value for signalling that there are no
    roots (zero-crossings) to monitor. *)
val no_roots : (int * ('a -> 'b -> 'c -> unit))

(** Return the number of root functions. *)
val nroots : ('a, 'k) session -> int

(** {2 Solver functions } *)

(** [(tret, r) = solv_normal s tout yout] integrates the ODE over an interval
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

    This routine will throw one of the solver {!exceptions} if an error
    occurs.

    @cvode <node5#sss:cvode> CVode (CV_NORMAL) *)
val solve_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                        -> float * Sundials.solver_result

(** This function is identical to {!solve_normal}, except that it returns after
    one internal solver step.

    @cvode <node5#sss:cvode> CVode (CV_ONE_STEP) *)
val solve_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                        -> float * Sundials.solver_result

(** {2 Main optional functions} *)

(** {3 Input} *)

(** Set the integration tolerances.

    @cvode <node5#sss:cvtolerances> CVodeSStolerances
    @cvode <node5#sss:cvtolerances> CVodeSVtolerances
    @cvode <node5#sss:cvtolerances> CVodeWFtolerances
    @cvode <node5#ss:ewtsetFn> Error weight function *)
val set_tolerances : ('a, 'k) session -> ('a, 'k) tolerance -> unit

(** [set_error_file s fname trunc] opens the file named [fname] and to which all
    messages from the default error handler are then directed. If the file
    already exists it is either trunctated ([trunc] = [true]) or appended to
    ([trunc] = [false]).

    The error file is closed if set_error_file is called again, or otherwise
    when the session is garbage collected.
   
    @cvode <node5#sss:optin_main> CVodeSetErrFile *)
val set_error_file : ('a, 'k) session -> string -> bool -> unit

(** [set_err_handler_fn s efun] specifies a custom function [efun] for handling
    error messages.

    @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
    @cvode <node5#ss:ehFn> Error message handler function *)
val set_err_handler_fn : ('a, 'k) session -> (Sundials.error_details -> unit) -> unit

(** This function restores the default error handling function. It is equivalent
    to calling CVodeSetErrHandlerFn with an argument of [NULL].

    @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn *)
val clear_err_handler_fn : ('a, 'k) session -> unit

(** Specifies the maximum order of the linear multistep method.

    @cvode <node5#sss:optin_main> CVodeSetMaxOrd *)
val set_max_ord : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of steps to be taken by the solver in its
    attempt to reach the next output time.

    @cvode <node5#sss:optin_main> CVodeSetMaxNumSteps *)
val set_max_num_steps : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of messages issued by the solver warning that
    t + h = t on the next internal step.

    @cvode <node5#sss:optin_main> CVodeSetMaxHnilWarns *)
val set_max_hnil_warns : ('a, 'k) session -> int -> unit

(** Indicates whether the BDF stability limit detection algorithm should be
    used.

    @cvode <node5#sss:optin_main> CVodeSetStabLimDet
    @cvode <node3#s:bdf_stab> BDF Stability Limit Detection *)
val set_stab_lim_det : ('a, 'k) session -> bool -> unit

(** Specifies the initial step size.

    @cvode <node5#sss:optin_main> CVodeSetInitStep *)
val set_init_step : ('a, 'k) session -> float -> unit

(** Specifies a lower bound on the magnitude of the step size.

    @cvode <node5#sss:optin_main> CVodeSetMinStep *)
val set_min_step : ('a, 'k) session -> float -> unit

(** Specifies an upper bound on the magnitude of the step size.

    @cvode <node5#sss:optin_main> CVodeSetMaxStep *)
val set_max_step : ('a, 'k) session -> float -> unit

(** Specifies the value of the independent variable t past which the solution is
    not to proceed. The default, if this routine is not called, is that no stop
    time is imposed.

    @cvode <node5#sss:optin_main> CVodeSetStopTime *)
val set_stop_time : ('a, 'k) session -> float -> unit

(** Specifies the maximum number of error test failures permitted in attempting
    one step.

    @cvode <node5#sss:optin_main> CVodeSetMaxErrTestFails *)
val set_max_err_test_fails : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver iterations permitted per
    step.

    @cvode <node5#sss:optin_main> CVodeSetMaxNonlinIters *)
val set_max_nonlin_iters : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of nonlinear solver convergence failures
    permitted during one step.

    @cvode <node5#sss:optin_main> CVodeSetMaxConvFails *)
val set_max_conv_fails : ('a, 'k) session -> int -> unit

(** Specifies the safety factor used in the nonlinear convergence test.

    @cvode <node5#sss:optin_main> CVodeSetNonlinConvCoef
    @cvode <node3#ss:ivp_sol> IVP Solution *)
val set_nonlin_conv_coef : ('a, 'k) session -> float -> unit

(** {3 Output } *)

(** Returns the real and integer workspace sizes.

    @cvode <node5#sss:optout_main> CVodeGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space          : ('a, 'k) session -> int * int

(** Returns the cumulative number of internal steps taken by the solver.

    @cvode <node5#sss:optout_main> CVodeGetNumSteps *)
val get_num_steps           : ('a, 'k) session -> int

(** Returns the number of calls to the user's right-hand side function.

    @cvode <node5#sss:optout_main> CVodeGetNumRhsEvals *)
val get_num_rhs_evals       : ('a, 'k) session -> int

(** Returns the number of calls made to the linear solver's setup function.

    @cvode <node5#sss:optout_main> CVodeGetNumLinSolvSetups *)
val get_num_lin_solv_setups : ('a, 'k) session -> int

(** Returns the number of local error test failures that have occurred.

    @cvode <node5#sss:optout_main> CVodeGetNumErrTestFails *)
val get_num_err_test_fails  : ('a, 'k) session -> int

(** Returns the integration method order used during the last internal step.

    @cvode <node5#sss:optout_main> CVodeGetLastOrder *)
val get_last_order          : ('a, 'k) session -> int

(** Returns the integration method order to be used on the next internal step.

    @cvode <node5#sss:optout_main> CVodeGetCurrentOrder *)
val get_current_order       : ('a, 'k) session -> int

(** Returns the integration step size taken on the last internal step.

    @cvode <node5#sss:optout_main> CVodeGetLastStep *)
val get_last_step           : ('a, 'k) session -> float

(** Returns the integration step size to be attempted on the next internal step.

    @cvode <node5#sss:optout_main> CVodeGetCurrentStep *)
val get_current_step        : ('a, 'k) session -> float

(** Returns the the value of the integration step size used on the first step.

    @cvode <node5#sss:optout_main> CVodeGetActualInitStep *)
val get_actual_init_step    : ('a, 'k) session -> float

(** Returns the the current internal time reached by the solver.

    @cvode <node5#sss:optout_main> CVodeGetCurrentTime *)
val get_current_time        : ('a, 'k) session -> float

(** Returns the number of order reductions dictated by the BDF stability limit
    detection algorithm.

    @cvode <node5#sss:optout_main> CVodeGetNumStabLimOrderReds
    @cvode <node3#s:bdf_stab> BDF stability limit detection *)
val get_num_stab_lim_order_reds : ('a, 'k) session -> int

(** Returns a suggested factor by which the user's tolerances should be scaled
    when too much accuracy has been requested for some internal step.

    @cvode <node5#sss:optout_main> CVodeGetTolScaleFactor *)
val get_tol_scale_factor : ('a, 'k) session -> float

(** Returns the solution error weights at the current time.

    @cvode <node5#sss:optout_main> CVodeGetErrWeights
    @cvode <node3#ss:ivp_sol> IVP solution (W_i) *)
val get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit

(** Returns the vector of estimated local errors.

    @cvode <node5#sss:optout_main> CVodeGetEstLocalErrors *)
val get_est_local_errors : ('a, 'k) session -> ('a, 'k) nvector -> unit

type integrator_stats = {
    num_steps : int;
    num_rhs_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

(** Returns the integrator statistics as a group.

    @cvode <node5#sss:optout_main> CVodeGetIntegratorStats *)
val get_integrator_stats    : ('a, 'k) session -> integrator_stats

(** Convenience function that calls get_integrator_stats and prints the results
    to stdout.

    @cvode <node5#sss:optout_main> CVodeGetIntegratorStats *)
val print_integrator_stats  : ('a, 'k) session -> unit


(** Returns the number of nonlinear (functional or Newton) iterations performed.

    @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('a, 'k) session -> int

(** Returns the number of nonlinear convergence failures that have occurred.

    @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails *)
val get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int

(** [nniters, nncfails = get_nonlin_solv_stats s] returns both the numbers of
    nonlinear iterations performed [nniters] and of nonlinear convergence
    failures that have occurred [nncfails].

    @cvode <node5#sss:optout_main> CVodeGetNonlinSolvStats *)
val get_nonlin_solv_stats : ('a, 'k) session -> int *int

(** {2 Root finding optional functions} *)

(** {3 Input} *)

(** [set_root_direction s dir] specifies the direction of zero-crossings to be
    located and returned. [dir] may contain one entry for each root function.

    @cvode <node5#sss:optin_root> CVodeSetRootDirection *)
val set_root_direction : ('a, 'k) session -> Sundials.RootDirs.root_direction array
                            -> unit

(** Like {!set_root_direction} but specifies a single direction for all root
    functions.

    @cvode <node5#sss:optin_root> CVodeSetRootDirection *)
val set_all_root_directions : ('a, 'k) session -> Sundials.RootDirs.root_direction
                                -> unit

(** Disables issuing a warning if some root function appears to be identically
    zero at the beginning of the integration.

    @cvode <node5#sss:optin_root> CVodeSetNoInactiveRootWarn *)
val set_no_inactive_root_warn : ('a, 'k) session -> unit

(** {3 Output} *)

(** Fills an array showing which functions were found to have a root.

    @cvode <node5#sss:optout_root> CVodeGetRootInfo *)
val get_root_info : ('a, 'k) session -> Sundials.Roots.t -> unit

(** Returns the cumulative number of calls made to the user-supplied root
    function g.

    @cvode <node5#sss:optout_root> CVodeGetNumGEvals *)
val get_num_g_evals : ('a, 'k) session -> int

(** {2 Interpolated output function } *)

(** [get_dky s t k dky] computes the [k]th derivative of the function y at time
    [t], i.e. d(k)y/dt(k)(t). The function requires that tn - hu <= [t] <= tn,
    where tn denotes the current internal time reached, and hu is the last
    internal step size successfully used by the solver. The user may request [k]
    = 0, 1,..., qu, where qu is the current order.

    This function may only be called after a successful return from either
    {!solve_normal} or {!solve_one_step}.

    Values for the limits may be obtained:
      - tn = {!get_current_time}
      - qu = {!get_last_order}
      - hu = {!get_last_step}

    @cvode <node5#sss:optin_root> CVodeGetDky *)
val get_dky : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit

(** {2 Reinitialization} *)

(** [reinit s ~iter_type:iter_type ~roots:roots t0 y0] reinitializes the
    solver session [s] with new values for the variables [y0].

    The labeled arguments are all optional, and if omitted, their current
    settings are kept.

    [iter_type] sets the nonlinear solver iteration type.  If omitted, the
    current iteration type will be kept.

    [roots] sets the root functions; see {!init} for what each component does.
    {!no_roots} may be specified to turn off root finding (if the session has
    been doing root finding until now).  If omitted, the current root functions
    or lack thereof will be kept.

    [t0] sets the value of the independent variable.  If omitted, the current
    time value will be kept.

    @cvode <node5#sss:cvreinit> CVodeReInit *)
val reinit :
  ('a, 'kind) session
  -> ?iter_type:('a, 'kind) iter
  -> ?roots:(int * (float -> 'a -> Sundials.Roots.val_array -> unit))
  -> float
  -> ('a, 'kind) nvector
  -> unit

