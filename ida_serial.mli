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
(*               User Documentation for IDA v2.7.0                     *)
(*         Alan C. Hindmarsh, Radu Serban, and Aaron Collier           *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

include module type of Ida
  with type Roots.t = Ida.Roots.t
  and type RootDirs.t = Ida.RootDirs.t
  and type solver_result = Ida.solver_result
  and type bandrange = Ida.bandrange

(*STARTINTRO*)
(** Serial nvector interface to the IDA solver.

  Serial vectors are passed between Sundials and OCaml programs as
  Bigarrays.
  These vectors are manipulated within the solver using the original low-level
  vector operations (cloning, linear sums, adding constants, and etcetera).
  While direct interfaces to these operations are not provided, there are
  equivalent implementations written in OCaml for arrays of floats
  ({! Nvector_array}) and bigarrays ({! Nvector_array.Bigarray}) of floats.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

(**
    This type represents a session with the IDA solver using serial nvectors
    accessed as {{:OCAML_DOC_ROOT(Bigarray.Array1)} Bigarray.Array1}s.

    A skeleton of the main program:
    + {b Set vector of initial values}
    {[let y = Ida.RealArray.of_array [| 0.0; 0.0; 0.0 |] ]}
    The length of this vector determines the problem size.
    + {b Create and initialize a solver session}
    {[let s = Ida.init Ida.Dense tols f ~roots:(2, g) y]}
    This will initialize a specific linear solver and the root-finding
    mechanism, if necessary.
    + {b Specify integration tolerances (optional)}, e.g.
    {[set_tolerances s SStolerances (reltol, abstol)]}
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
type nvec = Sundials.RealArray.t

(** The type of vectors containing dependent variable values, passed from the
   solver to callback functions. *)
type val_array = Sundials.RealArray.t

(** The type of vectors containing derivative values, passed from the
   solver to callback functions. *)
type der_array = Sundials.RealArray.t

(** The type of vectors containing detected roots (zero-crossings). *)
type root_array = Sundials.Roots.t

(** The type of vectors containing the values of root functions
   (zero-crossings). *)
type root_val_array = Sundials.Roots.val_array

(** {2 Linear Solvers} *)

type single_tmp = val_array
type double_tmp = val_array * val_array
type triple_tmp = val_array * val_array * val_array

(**
  Arguments common to all Jacobian callback functions.    
 
  @ida <node5#ss:djacFn> Dense Jacobian function
  @ida <node5#ss:bjacFn> Banded Jacobian function
  @ida <node5#ss:jtimesFn> Jacobian-times-vector function
  @ida <node5#ss:psolveFn> Linear preconditioning function
  @ida <node5#ss:precondFn> Jacobian preconditioning function
  @ida <node3#ss:ivp_sol> IVP solution
 *)
type 't jacobian_arg =
  {
    jac_t    : float;        (** The independent variable. *)
    jac_y    : val_array;    (** The dependent variable vector. *)
    jac_y'   : der_array;    (** The derivative vector (i.e. dy/dt). *)
    jac_res  : val_array;    (** The current value of the residual vector. *)
    jac_coef : float;        (** The coefficient [a] in the system Jacobian
                                   [J = dF/dy + a*dF/d(y')],
                                 where [F] is the residual function and
                                 d denotes partial differentiation.
                                 See the IVP solution section linked below.  *)
    jac_tmp  : 't            (** Workspace data, either {!single_tmp},
                                 {!double_tmp}, or {!triple_tmp}. *)
  }

(**
 Specify a linear solver.

 The Lapack solvers require that both Sundials and the OCaml interface were
 built to link with a LAPACK library.

 IDA supports direct linear solvers and Krylov solvers, but for the latter
 does not support banded preconditioning.

 @ida <node5#sss:lin_solv_init> Linear Solver Specification
                                                 Functions
 *)
type linear_solver =
  | Dense of dense_jac_fn option
  (** Direct linear solver with dense matrix.  The optional argument specifies
      a callback function that computes an approximation to the Jacobian matrix
      J(t, y).  If the argument is [None], then IDA uses a default
      implementation based on difference quotients.  [Dense None] can be passed
      into {!reinit} to disable any user-supplied Jacobian function that was
      previously active.

      @ida <node5#sss:lin_solve_init> IDADense
      @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
      @ida <node5#ss:djacFn> Dense Jacobian function
      @ida <node3#ss:ivp_soln> IVP solution  *)
  | LapackDense of dense_jac_fn option
  (** Direct linear solver with dense matrix, using LAPACK.  The argument is
      the same as [Dense].

      @ida <node5#sss:lin_solve_init> IDALapackDense
      @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
      @ida <node5#ss:djacFn> Dense Jacobian function
      @ida <node3#ss:ivp_soln> IVP solution  *)
  | Band of Ida.bandrange * band_jac_fn option
  (** Direct linear solver with banded matrix.  The arguments specify the width
      of the band ({!Ida.bandrange}) and an optional Jacobian function
      ({!band_jac_fn}).  If the Jacobian function is [None], IDA uses an
      internal implementation based on difference quotients.

      @ida <node5#sss:lin_solve_init> IDABand
      @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
      @ida <node5#ss:bjacFn> Banded Jacobian function
      @ida <node3#ss:ivp_soln> IVP solution *)
  | LapackBand of Ida.bandrange * band_jac_fn option
  (** Direct linear solver with banded matrix using LAPACK.  The arguments are
      the same as [Band].

      @ida <node5#sss:lin_solve_init> IDALapackBand
      @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
      @ida <node5#ss:bjacFn> Banded Jacobian function
      @ida <node3#ss:ivp_soln> IVP solution *)
  | Spgmr of spils_params
  (** Krylov iterative linear solver with the scaled preconditioned GMRES
      method.  See {!spils_params} for what the argument should contain.

      @ida <node5#sss:lin_solve_init> IDASpgmr
      @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
      @ida <node5#ss:psolveFn> Linear preconditioning function
      @ida <node5#ss:precondFn> Jacobian preconditioning function
    *)
  | Spbcg of spils_params
  (** Krylov iterative linear solver with the scaled preconditioned Bi-CGStab
      method. See {!spils_params} for what the argument should contain.

      @ida <node5#sss:lin_solve_init> IDASpbcg
      @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
      @ida <node5#ss:psolveFn> Linear preconditioning function
      @ida <node5#ss:precondFn> Jacobian preconditioning function
    *)
  | Sptfqmr of spils_params
  (** Krylov iterative linear solver with the scaled preconditioned TFQMR
      method.  See {!spils_params} for what the argument should contain.

      @ida <node5#sss:lin_solve_init> IDASptfqmr
      @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
      @ida <node5#ss:psolveFn> Linear preconditioning function
      @ida <node5#ss:precondFn> Jacobian preconditioning function
    *)

(** The type of a user-supplied callback function that computes an
    approximation to the Jacobian matrix for the [Dense] and [LapackDense]
    {!linear_solver}s.

    The function is called like [dense_jac_fn arg jac] where:
    - [arg] is the standard {!jacobian_arg} with three work vectors.
    - [jac] is the matrix in which to store the computed Jacobian.

    The function should load the ({i i,j}) entry of the Jacobian with {i
    dFi/dyj + c*dFi/dy'j}, i.e. the partial derivative of the {i i}-th equation
    with respect to the {i j}-th variable, evaluated at the values of ({i
    t,y,y'}) that can be obtained from [arg].  Note that in IDA, we have two
    terms due to the chain rule: one differentiated by the {i j}-th component
    of the non-derivative vector and the other by the {i j}-th component of the
    derivative vector.  The coefficient {i c} is the [jac_coef] field of the
    record [arg] (see {!jacobian_arg}).

    Only nonzero elements need to be loaded into [jac] because [jac] is set to
    the zero matrix before the call to the Jacobian function.

    If the user-supplied Jacobian function uses difference quotient
    approximations, then it may need to access quantities not in the argument
    list.  These include the current step size, the error weights, etc. To
    obtain these values, use the [get_*] functions defined in this module. The
    unit roundoff can be accessed as {!Sundials.unit_roundoff}.

    {b NB:} The elements of [arg] and [jac] must no longer be accessed after
    this function has returned.  If their values are needed outside of the
    function call, then they must be copied to separate physical
    structures.

    @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
    @ida <node5#ss:djacFn> Dense Jacobian function
    @ida <node3#ss:ivp_soln> IVP solution
*)
and dense_jac_fn = triple_tmp jacobian_arg -> Dls.DenseMatrix.t -> unit

(** The type of a user-supplied callback function that computes an
    approximation to the Jacobian matrix for the [Band] and [LapackBand]
    {!linear_solver}s.

    A user-supplied Jacobian function takes four arguments, in this order:
    - [arg] a {!jacobian_arg} with three work vectors.
    - [mupper] the upper half-bandwidth of the Jacobian.
    - [mlower] the lower half-bandwidth of the Jacobian.
    - [jac] the matrix to fill in with the values of the Jacobian.

    The function should load the ({i i,j}) entry of the Jacobian with {i
    dFi/dyj + c*dFi/dy'j}, i.e. the partial derivative of the {i i}-th equation
    with respect to the {i j}-th variable, evaluated at the values of ({i
    t,y,y'}) that can be obtained from [arg].  Note that in IDA, we have two
    terms due to the chain rule: one differentiated by the {i j}-th component
    of the non-derivative vector and the other by the {i j}-th component of the
    derivative vector.  The coefficient {i c} is the [jac_coef] field of the
    record [arg] (see {!jacobian_arg}).

    Only nonzero elements need to be loaded into [jac] because [jac] is set to
    the zero matrix before the call to the Jacobian function.

    If the user-supplied Jacobian function uses difference quotient
    approximations, then it may need to access quantities not in the argument
    list.  These include the current step size, the error weights, etc. To
    obtain these values, use the [get_*] functions defined in this module. The
    unit roundoff can be accessed as {!Sundials.unit_roundoff}.

    {b NB:} The elements of [arg] and [jac] must no longer be accessed after
    this function has returned.  If their values are needed outside of the
    function call, then they must be copied to separate physical
    structures.

    @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
    @ida <node5#ss:bjacFn> Banded Jacobian function
    @ida <node3#ss:ivp_soln> IVP solution
 *)
and band_jac_fn = Ida.bandrange -> triple_tmp jacobian_arg
                                -> Dls.BandMatrix.t -> unit

(** Initialization parameters and callbacks for Krylov iterative
    {!linear_solver}s.  Used with the {!linear_solver}s: [Spgmr], [Spbcg], and
    [Sptfqmr].  If you don't want any preconditioning, you should use
    {!spils_no_precond}, optionally with the [maxl] field overridden like
    [{ spils_no_precond with maxl = ... }].  *)
and spils_params =
  {
    maxl : int option;   (** Maximum dimension of the Krylov subspace to be
                             used. Pass [None] to use the default value [5]. *)

    prec_solve_fn : (single_tmp jacobian_arg -> val_array -> val_array -> float
                     -> unit) option;
    (** Called like [prec_solve_fn arg r z delta] to solve the linear system
        {i P}[z] = [r], where {i P} is the (left) preconditioner matrix.

        If set to [None] then no preconditioning is performed, and
        [prec_setup_fn] and [jac_times_vec_fn] are ignored.

        {i P} should approximate, at least crudely, the system Jacobian matrix
        {i J = dF/dy + {jac.coef} * dF/dy'} where {i F} is the residual
        function.
        - [arg] supplies the basic problem data as a {!jacobian_arg}.
        - [r] is the right-hand side vector.
        - [z] is the vector in which the result must be stored.
        - [delta] is an input tolerance.

        [delta] is an input tolerance to be used if an iterative method is
        employed in the solution.  In that case, the residual vector res = [r]
        - {i P} [z] of the system should be made less than [delta] in weighted
        l2 norm, i.e. [sqrt (sum over i ((res.{i} * ewt.{i})^2)) < delta],
        where the vector ewt can be obtained through {!get_err_weights}.

        This function can raise {!Sundials.RecoverableFailure} to instruct the
        integrator to retry with a different step size.  Raising any other kind
        of exception aborts the integrator.

        {b NB:} [r], [z], and the elements of [arg] must no longer be accessed
                after [prec_solve_fn] has returned, i.e. if their values are
                needed outside of the function call, then they must be copied
                to separate physical structures.

        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#ss:psolveFn> Linear preconditioning function
        @ida <node5#ss:precondFn> Jacobian preconditioning function
    *)

    prec_setup_fn : (triple_tmp jacobian_arg -> unit) option;
    (** A function that preprocesses and/or evaluates any Jacobian-related data
        needed by [prec_solve_fn] above.  When [prec_solve_fn] doesn't need any
        such data, this field can be [None].

        The sole argument to this function specifies the basic problem data as
        a {!jacobian_arg}.

        Note that unlike in CVODE, whatever data this function computes has to
        be recomputed every time it is called.

        This function can raise {!Sundials.RecoverableFailure} to instruct the
        integrator to retry with a different step size.  Raising any other kind
        of exception aborts the integrator.

        {b NB:} The elements of [jac] must no longer be accessed after [psetup]
                has returned a result, i.e. if their values are needed outside
                of the function call, then they must be copied to a separate
                physical structure.

        The operations performed by this function might include forming a crude
        approximate Jacobian, and performing an LU factorization on the
        resulting approximation.

        Each call to the preconditioner setup function is preceded by a call to
        the IDAResFn user function with the same (tt, yy, yp) arguments. Thus
        the preconditioner setup function can use any auxiliary data that is
        computed and saved during the evaluation of the DAE residual.

        This function is not called in advance of every call to the
        preconditioner solve function, but rather is called only as often as
        needed to achieve convergence in the Newton iteration.

        If this function uses difference quotient approximations, it may need
        to access quantities not in the argument.  These include the current
        step size, the error weights, etc.  To obtain these, use the [get_*]
        functions defined in this module.

        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#ss:psolveFn> Linear preconditioning function
        @ida <node5#ss:precondFn> Jacobian preconditioning function
    *)

    jac_times_vec_fn : (double_tmp jacobian_arg -> val_array -> val_array
                        -> unit) option;
    (**

       Specifies a Jacobian-times-vector function.  When this field is [None],
       IDA uses a default implementation based on difference quotients.

       [jac_times_vec_fn arg v jv] should compute the matrix-vector product {i
       J}[v], where {i J} is the system Jacobian.
       - [arg] provides the data necessary to compute the Jacobian.
       - [v] is the vector by which the Jacobian must be multiplied.
       - [jv] is the vector in which the result must be stored.

       The Jacobian {i J} (which is not explicitly constructed) has ({i i,j})
       entry {i dFi/dyj + c*dFi/dy'j} where {i F} is the residual function,
       i.e. the partial derivative of the [i]-th equation with respect to the
       [j]-th component of the non-derivative vector.  [c] is the [jac_coef]
       field of [arg] (see {!jacobian_arg}).  See the [Dense] {!linear_solver}
       for a more detailed explanation.

       {b NB:} The elements of [jac], [v], and [Jv] must no longer be accessed
               after [psolve] has returned a result, i.e. if their values are
               needed outside of the function call, then they must be copied to
               separate physical structures.

       Raising any kind of exception (including {!Sundials.RecoverableFailure})
       from this function results in the integrator being aborted.

       @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
       @ida <node5#ss:jtimesFn> Jacobian-times-vector function
     *)
  }

(** Specifies no preconditioning, with default dimension size.  See
    {!spils_params}.  *)
val spils_no_precond : spils_params

(** {3 Direct Linear Solvers (DLS)} *)

(** Get optional outputs for the Direct Linear Solvers that operate on dense
    and banded matrices.

    @ida <node5#sss:optin_dls> Direct linear solvers optional input functions
    @ida <node5#sss:optout_dls> Direct linear solvers optional output functions
  *)
module Dls :
  sig
    (** {4 Low-level solver manipulation} *)

    (** Change the dense Jacobian function (see [Dense] in {!linear_solver}).
        It may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.

        @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
        @ida <node5#ss:djacFn> Dense Jacobian function
        @ida <node3#ss:ivp_soln> IVP solution
      *)
    val set_dense_jac_fn : session -> dense_jac_fn -> unit

    (** Remove the user-supplied dense Jacobian function, if any, and fall back
        to IDA's internal implementation (see [Dense] in {!linear_solver}).
        This is the same as calling IDADlsSetDenseJacFn with an argument of
        [NULL].

        It may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.

        @ida <node5#sss:optin_dls> IDADlsSetDenseJacFn
        @ida <node5#ss:djacFn> Dense Jacobian function
     *)
    val clear_dense_jac_fn : session -> unit

    (** Change the band Jacobian function (see [Band] in {!linear_solver}).
        It may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.

        @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
        @ida <node5#ss:bjacFn> Banded Jacobian function
      *)
    val set_band_jac_fn : session -> band_jac_fn -> unit

    (** Remove the user-supplied band Jacobian function, if any, and fall back
        to IDA's internal implementation (see [Band] in {!linear_solver}).  It
        may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.

        @ida <node5#sss:optin_dls> IDADlsSetBandJacFn
        @ida <node5#ss:bjacFn> Banded Jacobian function
      *)
    val clear_band_jac_fn : session -> unit

    (** {4 Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the Dense
      and Band direct linear solvers .

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

(** Set optional inputs, and get optional outputs for the Scaled Preconditioned
    Iterative Linear Solvers: SPGMR, SPBCG, SPTFQMR.

    @ida <node5#sss:optin_spils> Iterative linear solvers optional input
                                 functions.
    @ida <node5#sss:optout_spils> Iterative linear solvers optional output
                                  functions.
    @ida <node5#ss:psolveFn> Linear preconditioning function
    @ida <node5#ss:precondFn> Jacobian preconditioning function
 *)
module Spils :
  sig
    (** {4 Low-level solver manipulation} *)

    (** Set preconditioning functions (see {!spils_params}).  It may be
        unsafe to use this function without a {!reinit}.  Users are encouraged
        to use the [iter_type] parameter of {!reinit} instead, unless they are
        desperate for performance.

        @ida <node5#sss:optin_spils> IDASpilsSetPreconditioner
        @ida <node5#ss:psolveFn> Linear preconditioning function
      *)
     val set_preconditioner :
       session
       -> (triple_tmp jacobian_arg -> unit)
       -> (single_tmp jacobian_arg -> val_array -> val_array -> float -> unit)
       -> unit

    (** Set the Jacobian-times-vector function (see {!spils_params}).  It
        may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead,
        unless they are desperate for performance.

        @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
        @ida <node5#ss:jtimesFn> Jacobian-times-vector function
      *)
    val set_jac_times_vec_fn :
      session
      -> (double_tmp jacobian_arg
          -> val_array (* v *)
          -> val_array (* Jv *)
          -> unit)
      -> unit

    (** This function disables the user-supplied Jacobian-vector function, and
        switches back to the default internal difference quotient approximation
        (see {!spils_params}).  It is equivalent to calling
        IDASpilsSetJacTimesVecFn with an argument of [NULL].

        It may be unsafe to use this function without a {!reinit}.  Users are
        encouraged to use the [iter_type] parameter of {!reinit} instead, unless
        they are desperate for performance.

        @ida <node5#sss:optin_spils> IDASpilsSetJacTimesVecFn
        @ida <node5#ss:jtimesFn> Jacobian-times-vector function
    *)
    val clear_jac_times_vec_fn : session -> unit

    (** {4 Optional output functions} *)

    (** Sets the Gram-Schmidt orthogonalization to be used with the [Spgmr]
        {!linear_solver}.

      @ida <node5#sss:optin_spils> IDASpilsSetGSType
    *)
    val set_gs_type : session -> Spils.gramschmidt_type -> unit

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
      calls made to [prec_setup_fn] (see {!spils_params}).

      @ida <node5#sss:optout_spils> IDASpilsGetNumPrecEvals
    *)
    val get_num_prec_evals   : session -> int

    (**
      Returns the cumulative number of calls made to the preconditioner solve
      function, [prec_solve_fn] (see {!spils_params}).

      @ida <node5#sss:optout_spils> IDASpilsGetNumPrecSolves
    *)
    val get_num_prec_solves  : session -> int

    (**
      Returns the cumulative number of calls made to [jac_times_vec_fn]
      (see {!spils_params}).

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

type tolerance =
  | SStolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVtolerances of float * nvec
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFtolerances of (val_array -> val_array -> unit)
    (** Specifies a function [efun y ewt] that sets the multiplicative
        error weights Wi for use in the weighted RMS norm. The function is
        passed the dependent variable vector [y] and is expected to set the
        values inside the error-weight vector [ewt]. *)

(** A default relative tolerance of 1.0e-4 and absolute tolerance of 1.0e-8. *)
val default_tolerances : tolerance

(** {2 Initialization} *)

(**
    [init linsolv tol f ~roots:(nroots, g) ~t0:t0 y0 y'0] initializes the IDA
    solver to solve the DAE f t y y' = 0 and returns a {!session}.

    - [linsolv] is the linear solver to attach to this solver,
    - [tol]     specifies the integration tolerances,
    - [f]       is the residual function (see below),
    - [nroots]  specifies the number of root functions (zero-crossings),
    - [g]       calculates the values of the root functions,
    - [t0]      is the initial value of the independent variable t, which
                defaults to 0,
    - [y0]      is a vector of initial values for the dependent-variable vector
                [y].  This vector's size determines the number of equations
                in the session, see {!Sundials.RealArray.t}, and,
    - [y'0]     is a vector of initial values for [y'], i.e. the derivative
                of [y] with respect to t.  This vector's size must match the
                size of [y0].

    This function calls IDACreate, IDAInit, IDARootInit, an appropriate linear
    solver function, and one of IDASStolerances, IDASVtolerances, or
    IDAWFtolerances. It does everything necessary to initialize an IDA session;
    the {!solve_normal} or {!solve_one_step} functions can be called directly
    afterward.

    The residual function [f] is called by the solver like [f t y y' r] to
    compute the problem residual, where:
    - [t] is the current value of the independent variable,
          i.e., the simulation time.
    - [y] is a vector of dependent-variable values, i.e. y(t).
    - [y'] is the derivative of [y] with respect to [t], i.e. dy/dt.
    - [r] is the output vector to fill in with the value of the residual
          function for the given values of t, y, and y'.
    The residual function should return normally if successful, raise
    {!Sundials.RecoverableFailure} if a recoverable error occurred (e.g. [y] has
    an illegal value), or raise some other exception if a nonrecoverable error
    occurred.  If a recoverable error occurred, the integrator will attempt to
    correct and retry.  If a nonrecoverable error occurred, the integrator will
    halt and propagate the exception to the caller.

    {b NB:} [y], [y'], and [r] must no longer be accessed after [f] has
            returned a result, i.e. if their values are needed outside of
            the function call, then they must be copied to separate physical
            structures.

    The roots function [g], if supplied, is called by the solver to calculate
    the values of root functions (zero-crossing expressions) which are used to
    detect significant events.  It is passed four arguments [t], [y], [y'], and
    [gout]:
    - [t], [y], [y'] are as for [f].
    - [gout] is a vector for storing the values of g(t, y, y').
    If the labeled argument ~roots is omitted, then no root finding is
    performed.  If the root function raises an exception, the integrator will
    halt immediately and propagate the exception to the caller.

    {b NB:} [y] and [gout] must no longer be accessed after [g] has returned
            a result, i.e. if their values are needed outside of the function
            call, then they must be copied to separate physical structures.

    @ida <node5#sss:idainit>       IDACreate/IDAInit
    @ida <node5#ss:resFn>          DAE residual function
    @ida <node5#ss:idarootinit>    IDARootInit
    @ida <node5#ss:rootFn>         Rootfinding function
    @ida <node5#sss:lin_solv_init> Linear solvers
    @ida <node5#sss:idatolerances> IDASStolerances
 *)
val init :
    linear_solver
    -> tolerance
    -> (float -> val_array -> der_array -> val_array -> unit)
    -> ?roots:(int *
               (float -> val_array -> der_array -> root_val_array -> unit))
    -> ?t0:float
    -> nvec
    -> nvec
    -> session

(** Return the number of root functions. *)
val nroots : session -> int

(** Return the number of equations. *)
val neqs : session -> int

(** {2 Solver functions } *)

(**
   [(tret, r) = solve_normal s tout yout y'out] integrates the DAE over an
   interval in t.

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
   set by {!set_stop_time}, was reached; see {!Sundials.solver_result}.

   This routine will throw one of the solver {!Ida.exceptions} if an error
   occurs.

   @ida <node5#sss:idasolve> IDASolve (IDA_NORMAL)
 *)
val solve_normal :
  session -> float -> val_array -> der_array -> float * solver_result

(**
   This function is identical to {!solve_normal}, except that it returns after one
   internal solver step.

   @ida <node5#sss:idasolve> IDASolve (IDA_ONE_STEP)
 *)
val solve_one_step :
  session -> float -> val_array -> der_array -> float * solver_result

(** {2 Main optional functions} *)

(** {3 Input} *)

(** Set the integration tolerances.

    @ida <node5#sss:idatolerances> IDASStolerances
    @ida <node5#sss:idatolerances> IDASVtolerances
    @ida <node5#sss:idatolerances> IDAWFtolerances
    @ida <node5#ss:ewtsetFn> Error weight function
 *)
val set_tolerances : session -> tolerance -> unit

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
val set_err_handler_fn : session -> (Sundials.error_details -> unit) -> unit

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
  {!solve_normal} or {!solve_one_step}.

  Values for the limits may be obtained:
    - tn = {!get_current_time}
    - qu = {!get_last_order}
    - hu = {!get_last_step}

  @ida <node5#sss:optin_root> IDAGetDky
 *)
val get_dky : session -> float -> int -> nvec -> unit

(** {2 Reinitialization} *)

(**

  [reinit s ~linsolv:linsolv ~roots:(nroots, g) t0 y0 y'0] reinitializes the
  solver session [s] with a new time [t0] and new values for the variables
  [y0].  There are two optional arguments to change the linear solver and the
  set of root functions.

  The optional argument [linsolv] sets the linear solver.  If omitted, the
  current linear solver will be kept.  If a session is created with, say,
  [Dense (Some f)], and then reinitialized with [Dense None], then the linear
  solver is reset and [f] is removed from the session.  The same goes for all
  other optional callbacks.

  The optional argument [roots] sets the root functions; see {!init} for what
  each component does.  {!Ida.no_roots} may be passed in to turn off root
  finding.  If omitted, the current root functions will be kept.

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

(** Inequality constraints on variables.

 @ida <node5#sss:optin_main> IDASetConstraints
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

val set_constraints : session -> Constraints.t -> unit

(**
   Specify whether or not to ignore algebraic variables in local error tests.
   If you set this option to [true], you must also specify which variables are
   algebraic through {!calc_ic_ya_yd'} or {!set_var_types}.  Exactly one of
   these functions should be called, exactly once, before the first call to
   {!solve_normal}, {!solve_one_step}, or {!calc_ic_ya_yd'}.  Forgetting to do
   so will cause an {!Ida.IllInput} exception.

   Note: {!set_var_types} is the preferred alias to {!set_id}, which
   corresponds to [IDASetId] in the C interface.

   In general, suppressing local error tests for algebraic variables is {i
   discouraged} when solving DAE systems of index 1, whereas it is generally {i
   encouraged} for systems of index 2 or more.  See pp. 146-147 of the
   following reference for more on this issue:

   K. E. Brenan, S. L. Campbell, and L. R. Petzold.  Numerical Solution of
   Initial-Value Problems in Differential-Algebraic Equations.  SIAM,
   Philadelphia, Pa, 1996.

   @ida <node5#sss:optin_main> IDASetSuppressAlg
 *)
val set_suppress_alg : session -> bool -> unit

(** An unpreferred alias for {!set_var_types}.  SUNDIALS calls variable types
    by the cryptic name "Id", and this OCaml binding preserves this alternative
    naming to help users transition from the C interface.  *)
val set_id : session -> Id.t -> unit

(** Specify which variables are algebraic and which variables are differential,
    needed for {!set_suppress_alg}.  The function {!calc_ic_ya_yd'} also sets
    the same information, so you only need to call one of these functions.

    The SUNDIALS manual is not clear about whether it's safe to change the
    variable types after you've already set it.

    [set_var_types] corresponds to [IDASetId] in the C interface, and an alias
    {!set_id} is also available in this binding.  We prefer the more
    descriptive name {!set_var_types}, however.

    @ida <node5#sss:optin_main> IDASetId
 *)
val set_var_types : session -> VarTypes.t -> unit

(** [calc_ic_y ida ~y:yvar tout1] corrects the initial values y0 at time t0,
    using the initial values of the derivatives y'0 that are stored in the
    session.  That is, if the {i t0,y0,y'0} that were given to {!init} or
    {!reinit} does not satisfy {i F(t0,y0,y'0) = 0}, where {i F} is the
    residual function, then [calc_ic_y] will modify {i y'0} so that this
    equation holds.  If {i F(t0,y0,y'0) = 0} is already true, a call to
    [calc_ic_y] is unnecessary.

    The optional parameter [~y], if given, will receive the corrected {i y}
    vector.  [tout1] is the first value of {i t} at which a solution will be
    requested (using {!solve_normal} or {!solve_one_step}). This value is
    needed here only to determine the direction of integration and rough scale
    in the independent variable {i t}.

    IDA's initial value correction works for certain index-one problems
    including a class of systems of semi-implicit form, and uses Newton
    iteration combined with a linesearch algorithm.  See Section 2.1 of the IDA
    User Guide and the following reference for more information:

    P. N. Brown, A. C. Hindmarsh, and L. R. Petzold. Consistent Initial Condition Calculation for Differential-Algebraic Systems. SIAM J. Sci. Comput., 19:1495-1512, 1998.

    @ida <node5#ss:idacalcic> IDACalcIC
    @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
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

      - [y]  splits into [A]  = \{ [y.{i}]  | i in Ia \} and D  = \{ [y.{i}]  | [i] in [Id] \}
      - [y'] splits into [A'] = \{ [y'.{i}] | i in Ia \} and D' = \{ [y'.{i}] | [i] in [Id] \}

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

    @ida <node5#ss:idacalcic> IDACalcIC
    @ida <node5#sss:optin_main> IDASetId
    @ida <node5#sss:optin_main> IDASetSuppressAlg
    @ida <node5#sss:optout_iccalc> IDAGetConsistentIC
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
    @ida <node5#sss:optout_iccalc> IDAGetNumBcktrackOps
 *)
val get_num_backtrack_ops : session -> int
