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
(* Parts of the comment text are taken directly from:                  *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Variable-step solution of ODE initial value problems with
    zero-crossing detection.

    This module solves numerically problems of the form
    {% $\dot{y} = f(t, y)$%}, {% $y(t_0) = y_0$%}.

    This documented interface is structured as follows.
    {ol
      {- {{:#linear}Linear solvers}}
      {- {{:#tols}Tolerances}}
      {- {{:#solver}Solver initialization and use}}
      {- {{:#set}Modifying the solver}}
      {- {{:#get}Querying the solver}}
      {- {{:#roots}Additional root finding functions}}
      {- {{:#exceptions}Exceptions}}}

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS) *)

open Sundials

(** A session with the CVODE solver.

    An example session with Cvode ({openfile cvode_skel.ml}): {[
#include "examples/ocaml/skeletons/cvode_skel.ml"
    ]}

    @cvode <node5#ss:skeleton_sim> Skeleton of main program *)
type ('a, 'k) session = ('a, 'k) Cvode_impl.session

(** Alias for sessions based on serial nvectors. *)
type serial_session = (Nvector_serial.data, Nvector_serial.kind) session

(** {2:linear Linear solvers} *)

(** Linear solvers used by Cvode.

    @cvode <node5#sss:lin_solv_init> Linear Solver Specification Functions *)
type ('data, 'kind) linear_solver = ('data, 'kind) Cvode_impl.linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type serial_linear_solver =
      (Nvector_serial.data, Nvector_serial.kind) linear_solver

(** Workspaces with three temporary vectors. *)
type 'a triple = 'a * 'a * 'a

(** Arguments common to Jacobian callback functions.    
 
    @cvode <node5#ss:djacFn> CVDlsDenseJacFn
    @cvode <node5#ss:bjacFn> CVDlsBandJacFn
    @cvode <node5#ss:jtimesfn> CVSpilsJacTimesVecFn
    @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
    @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;        (** The independent variable. *)
    jac_y   : 'a;           (** The dependent variable vector. *)
    jac_fy  : 'a;           (** The derivative vector (i.e., dy/dt). *)
    jac_tmp : 't            (** Workspace data. *)
  }

(** The range of nonzero entries in a band matrix. *)
type bandrange = Cvode_impl.bandrange =
  { mupper : int; (** The upper half-bandwidth. *)
    mlower : int; (** The lower half-bandwidth. *) }

(** Diagonal approximation of Jacobians by difference quotients. *)
module Diag :
  sig
    (** A linear solver based on Jacobian approximation by difference
        quotients.

        @cvode <node5#sss:lin_solv_init> CVDiag *)
    val solver : ('data, 'kind) linear_solver

    (** Returns the sizes of the real and integer workspaces used by the
        Diagonal linear solver.

        @cvode <node5#sss:optout_diag> CVDiagGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : ('a, 'k) session -> int * int

    (** Returns the number of calls made to the right-hand side
        function due to finite difference Jacobian approximation in the
        Diagonal linear solver.

        @cvode <node5#sss:optout_diag> CVDiagGetNumRhsEvals *)
    val get_num_rhs_evals : ('a, 'k) session -> int
  end

(** Direct Linear Solvers operating on dense and banded matrices.

    @cvode <node5#sss:optin_dls> Direct linear solvers optional input functions
    @cvode <node5#sss:optout_dls> Direct linear solvers optional output functions *)
module Dls :
  sig

    (** Callback functions that compute dense approximations to a Jacobian
        matrix. In the call [dense_jac_fn arg jac], [arg] is a {!jacobian_arg}
        with three work vectors and the computed Jacobian must be stored
        in [jac].

        The callback should load the [(i,j)]th entry of [jac] with
        {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of the
        [i]th equation with respect to the [j]th variable, evaluated at the
        values of [t] and [y] obtained from [arg]. Only nonzero elements need
        be loaded into [jac].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jac] should
                 be accessed after the function has returned.}

        @cvode <node5#ss:djacFn> CVDlsDenseJacFn *)
    type dense_jac_fn = (RealArray.t triple, RealArray.t) jacobian_arg
                                                  -> Dls.DenseMatrix.t -> unit

    (** A direct linear solver on dense matrices. The optional argument
        specifies a callback function for computing an approximation to the
        Jacobian matrix. If this argument is omitted, then a default
        implementation based on difference quotients is used.

        @cvode <node5#sss:lin_solv_init> CVDense
        @cvode <node5#sss:optin_dls> CVDlsSetDenseJacFn
        @cvode <node5#ss:djacFn> CVDlsDenseJacFn *)
    val dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

    (** A direct linear solver on dense matrices using LAPACK. See {!dense}.
        Only available if {!Sundials.lapack_enabled}.

        @cvode <node5#sss:lin_solv_init> CVLapackDense
        @cvode <node5#sss:optin_dls> CVDlsSetDenseJacFn
        @cvode <node5#ss:djacFn> CVDlsDenseJacFn *)
    val lapack_dense : ?jac:dense_jac_fn -> unit -> serial_linear_solver

    (** Callback functions that compute banded approximations to
        a Jacobian matrix. In the call [band_jac_fn {mupper; mlower} arg jac],
        - [mupper] is the upper half-bandwidth of the Jacobian,
        - [mlower] is the lower half-bandwidth of the Jacobian,
        - [arg] is a {!jacobian_arg} with three work vectors, and,
        - [jac] is storage for the computed Jacobian.

        The callback should load the [(i,j)]th entry of [jac] with
        {% $\partial y_i/\partial y_j$%}, i.e., the partial derivative of
        the [i]th equation with respect to the [j]th variable, evaluated
        at the values of [t] and [y] obtained from [arg]. Only nonzero
        elements need be loaded into [jac].

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor the matrix [jac] should
                 be accessed after the function has returned.}

        @cvode <node5#ss:bjacFn> CVDlsBandJacFn *)
    type band_jac_fn = bandrange
                        -> (RealArray.t triple, RealArray.t) jacobian_arg
                        -> Dls.BandMatrix.t -> unit

    (** A direct linear solver on banded matrices. The optional argument
        specifies a callback function for computing an approximation to the
        Jacobian matrix. If this argument is omitted, then a default
        implementation based on difference quotients is used. The other
        argument gives the width of the bandrange.

        @cvode <node5#sss:lin_solv_init> CVBand
        @cvode <node5#sss:optin_dls> CVDlsSetBandJacFn
        @cvode <node5#ss:bjacFn> CVDlsBandJacFn *)
    val band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

    (** A direct linear solver on banded matrices using LAPACK. See {!band}.
        Only available if {!Sundials.lapack_enabled}.

        @cvode <node5#sss:lin_solv_init> CVLapackBand
        @cvode <node5#sss:optin_dls> CVDlsSetBandJacFn
        @cvode <node5#ss:bjacFn> CVDlsBandJacFn *)
    val lapack_band : ?jac:band_jac_fn -> bandrange -> serial_linear_solver

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by a direct
        linear solver.

        @cvode <node5#sss:optout_dls> CVDlsGetWorkSpace
        @return ([real_size], [integer_size]) *)
    val get_work_space : serial_session -> int * int

    (** Returns the number of calls made by a direct linear solver to the
        Jacobian approximation function.

        @cvode <node5#sss:optout_dls> CVDlsGetNumJacEvals *)
    val get_num_jac_evals : serial_session -> int

    (** Returns the number of calls to the right-hand side callback due to
        the finite difference Jacobian approximation.

        @cvode <node5#sss:optout_dls> CVDlsGetNumRhsEvals *)
    val get_num_rhs_evals : serial_session -> int

    (** {3:lowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change a Jacobian function. These low-level functions are provided for
        experts who want to avoid resetting internal counters and other
        associated side-effects. *)

    (** Change the dense Jacobian function.
   
        @cvode <node5#sss:optin_dls> CVDlsSetDenseJacFn *)
    val set_dense_jac_fn : serial_session -> dense_jac_fn -> unit

    (** Remove a dense Jacobian function and use the default
        implementation.

        @cvode <node5#sss:optin_dls> CVDlsSetDenseJacFn *)
    val clear_dense_jac_fn : serial_session -> unit

    (** Change the band Jacobian function.

        @cvode <node5#sss:optin_dls> CVDlsSetBandJacFn *)
    val set_band_jac_fn : serial_session -> band_jac_fn -> unit

    (** Remove a banded Jacobian function and use the default
        implementation.

        @cvode <node5#sss:optin_dls> CVDlsSetBandJacFn *)
    val clear_band_jac_fn : serial_session -> unit
  end

(** Scaled Preconditioned Iterative Linear Solvers.

    @cvode <node5#sss:optin_spils> Iterative linear solvers optional input functions.
    @cvode <node5#sss:optout_spils> Iterative linear solvers optional output functions.
    @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
    @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
module Spils :
  sig
    (** {3:precond Preconditioners} *)

    (** Arguments passed to the preconditioner solver function.

        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn *)
    type 'a prec_solve_arg =
      {
        rhs   : 'a;         (** Right-hand side vector of the linear system. *)
        gamma : float;      (** Scalar $\gamma$ in the Newton
                                matrix given by $M = I - \gamma J$. *)
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
        the Newton matrix {% $M = I - \gamma J$%} where
        {% $J = \frac{\partial f}{\partial y}$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac], [arg], and [z] should not
                 be accessed after the function has returned.}

        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn *)
    type 'a prec_solve_fn =
      ('a, 'a) jacobian_arg
      -> 'a prec_solve_arg
      -> 'a
      -> unit

    (** Callback functions that preprocess or evaluate Jacobian-related data
        need by {!prec_solve_fn}. In the call [prec_setup_fn jac jok gamma],
        [jac] is a {!jacobian_arg} with one work vector, [jok] indicates
        whether any saved Jacobian-related data can be reused with the current
        value of [gamma], and [gamma] is the scalar $\gamma$ in the Newton
        matrix {% $M = I - \gamma J$%} where $J$ is the Jacobian matrix.
        A function should return [true] if Jacobian-related data was updated
        and [false] if saved data was reused.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The elements of [jac] should not be accessed after the
                 function has returned.}

        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
    type 'a prec_setup_fn =
      ('a triple, 'a) jacobian_arg
      -> bool
      -> float
      -> bool

    (** Callback functions that compute the Jacobian times a vector. In the
        call [jac_times_vec_fn arg v jv], [arg] is a {!jacobian_arg} with one
        work vector, [v] is the vector multiplying the Jacobian, and [jv] is
        the vector in which to store the
        result—{% $\mathtt{jv} = J\mathtt{v}$%}.
      
        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning Neither the elements of [arg] nor [v] or [jv] should be
                 accessed after the function has returned.}

        @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
        @cvode <node5#ss:jtimesfn> CVSpilsJacTimesVecFn *)
    type 'a jac_times_vec_fn =
      ('a, 'a) jacobian_arg
      -> 'a (* v *)
      -> 'a (* Jv *)
      -> unit

    (** Specifies a preconditioner, including the type of preconditioning
        (none, left, right, or both) and callback functions.
        The following functions and those in {!Banded} and {!Cvode_bbd}
        construct preconditioners.

        The {!prec_solve_fn} is usually mandatory. The {!prec_setup_fn} can be
        omitted if not needed. If the {!jac_times_vec_fn} is omitted, a
        default implementation based on difference quotients is used.

        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn
        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
        @cvode <node5#ss:jtimesfn> CVSpilsJacTimesVecFn *)
    type ('a,'k) preconditioner = ('a,'k) Cvode_impl.SpilsTypes.preconditioner

    (** No preconditioning.  *)
    val prec_none : ('a, 'k) preconditioner

    (** Left preconditioning. {% $(P^{-1}A)x = P^{-1}b$ %}. *)
    val prec_left :
      ?setup:'a prec_setup_fn
      -> ?jac_times_vec:'a jac_times_vec_fn
      -> 'a prec_solve_fn
      -> ('a, 'k) preconditioner

    (** Right preconditioning. {% $(AP^{-1})Px = b$ %}. *)
    val prec_right :
      ?setup:'a prec_setup_fn
      -> ?jac_times_vec:'a jac_times_vec_fn
      -> 'a prec_solve_fn
      -> ('a, 'k) preconditioner

    (** Left and right preconditioning.
        {% $(P_L^{-1}AP_R^{-1})P_Rx = P_L^{-1}b$ %} *)
    val prec_both :
      ?setup:'a prec_setup_fn
      -> ?jac_times_vec:'a jac_times_vec_fn
      -> 'a prec_solve_fn
      -> ('a, 'k) preconditioner

    (** Banded preconditioners.  *)
    module Banded : sig

      (** A band matrix {!preconditioner} based on difference quotients.
          The call [prec_left br] instantiates a left preconditioner which
          generates a banded approximation to the Jacobian with [br.mlower]
          sub-diagonals and [br.mupper] super-diagonals.

          @cvode <node5#sss:cvbandpre> CVBandPrecInit *)
      val prec_left :
        ?jac_times_vec:(RealArray.t jac_times_vec_fn)
        -> bandrange
        -> (Nvector_serial.data, Nvector_serial.kind) preconditioner

      (** Like {!prec_left} but preconditions from the right.

          @cvode <node5#sss:cvbandpre> CVBandPrecInit *)
      val prec_right :
        ?jac_times_vec:(RealArray.t jac_times_vec_fn)
        -> bandrange
        -> (Nvector_serial.data, Nvector_serial.kind) preconditioner

      (** Like {!prec_left} but preconditions from both sides.

          @cvode <node5#sss:cvbandpre> CVBandPrecInit *)
      val prec_both :
        ?jac_times_vec:(RealArray.t jac_times_vec_fn)
        -> bandrange
        -> (Nvector_serial.data, Nvector_serial.kind) preconditioner

      (** {4:stats Banded statistics} *)

      (** Returns the sizes of the real and integer workspaces used by the
          banded preconditioner module.

          @cvode <node5#sss:cvbandpre> CVBandPrecGetWorkSpace
          @return ([real_size], [integer_size]) *)
      val get_work_space : serial_session -> int * int

      (** Returns the number of calls to the right-hand side callback for the
          difference banded Jacobian approximation. This counter is only updated
          if the default difference quotient function is used.

          @cvode <node5#sss:cvbandpre> CVBandPrecGetNumRhsEvals *)
      val get_num_rhs_evals : serial_session -> int
    end

    (** {3:lsolvers Solvers} *)

    (** Krylov iterative solver using the scaled preconditioned generalized
        minimum residual (GMRES) method.
        In the call [spgmr ~maxl:maxl prec], [maxl] is the maximum dimension of
        the Krylov subspace (defaults to 5), and [prec] is a {!preconditioner}.

        @cvode <node5#sss:lin_solv_init> CVSpgmr
        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#sss:optin_spils> CVSpilsSetMaxl *)
    val spgmr : ?maxl:int -> ('a, 'k) preconditioner -> ('a, 'k) linear_solver

    (** Krylov iterative solver using the scaled preconditioned biconjugate
        stabilized (Bi-CGStab) method.
        In the call [spbcg ~maxl:maxl prec], [maxl] is the maximum dimension of
        the Krylov subspace (defaults to 5), and [prec] is a {!preconditioner}.

        @cvode <node5#sss:lin_solv_init> CVSpbcg
        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#sss:optin_spils> CVSpilsSetMaxl *)
    val spbcg : ?maxl:int -> ('a, 'k) preconditioner -> ('a, 'k) linear_solver

    (** Krylov iterative with the scaled preconditioned transpose-free
        quasi-minimal residual (SPTFQMR) method.
        In the call [sptfqmr ~maxl:maxl prec], [maxl] is the maximum dimension
        of the Krylov subspace (defaults to 5), and [prec] is a
        {!preconditioner}.

        @cvode <node5#sss:lin_solv_init> CVSptfqmr
        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#sss:optin_spils> CVSpilsSetMaxl *)
    val sptfqmr : ?maxl:int -> ('a, 'k) preconditioner -> ('a, 'k) linear_solver

    (** {3:set Solver parameters} *)

    (** The type of Gram-Schmidt orthogonalization.

        @cvode <node9#ss:spgmr> ModifiedGS/ClassicalGS *)
    type gramschmidt_type = Spils.gramschmidt_type =
      | ModifiedGS   (** Modified Gram-Schmidt orthogonalization
                         {cconst MODIFIED_GS} *)
      | ClassicalGS  (** Classical Gram Schmidt orthogonalization
                         {cconst CLASSICAL_GS} *)

    (** Sets the Gram-Schmidt orthogonalization to be used with the
        Spgmr {!linear_solver}.

        @cvode <node5#sss:optin_spils> CVSpilsSetGSType *)
    val set_gs_type : ('a, 'k) session -> gramschmidt_type -> unit

    (** Sets the factor by which the Krylov linear solver's convergence test
        constant is reduced from the Newton iteration test constant.
        This factor must be >= 0; passing 0 specifies the default (0.05).

        @cvode <node5#sss:optin_spils> CVSpilsSetEpsLin *)
    val set_eps_lin : ('a, 'k) session -> float -> unit

    (** Resets the maximum Krylov subspace dimension for the Bi-CGStab and
        TFQMR methods. A value <= 0 specifies the default (5.0).

        @cvode <node5#sss:optin_spils> CVSpilsSetMaxl *)
    val set_maxl : ('a, 'k) session -> int -> unit

    (** {3:stats Solver statistics} *)

    (** Returns the sizes of the real and integer workspaces used by the spils
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

    (** Returns the cumulative number of calls to the setup function with
        [jok=false].

        @cvode <node5#sss:optout_spils> CVSpilsGetNumPrecEvals *)
    val get_num_prec_evals   : ('a, 'k) session -> int

    (** Returns the cumulative number of calls to the preconditioner solve
        function.

        @cvode <node5#sss:optout_spils> CVSpilsGetNumPrecSolves *)
    val get_num_prec_solves  : ('a, 'k) session -> int

    (** Returns the cumulative number of calls to the Jacobian-vector
        function.

        @cvode <node5#sss:optout_spils> CVSpilsGetNumJtimesEvals *)
    val get_num_jtimes_evals : ('a, 'k) session -> int

    (** Returns the number of calls to the right-hand side callback for
        finite difference Jacobian-vector product approximation. This counter is
        only updated if the default difference quotient function is used.

        @cvode <node5#sss:optout_spils> CVSpilsGetNumRhsEvals *)
    val get_num_rhs_evals    : ('a, 'k) session -> int

    (** {3:lowlevel Low-level solver manipulation}

        The {!init} and {!reinit} functions are the preferred way to set or
        change preconditioner functions. These low-level functions are provided
        for experts who want to avoid resetting internal counters and other
        associated side-effects. *)

    (** Change the preconditioner functions.

        @cvode <node5#sss:optin_spils> CVSpilsSetPreconditioner
        @cvode <node5#ss:psolveFn> CVSpilsPrecSolveFn
        @cvode <node5#ss:precondFn> CVSpilsPrecSetupFn *)
    val set_preconditioner :
      ('a, 'k) session
      -> ?setup:'a prec_setup_fn
      -> 'a prec_solve_fn
      -> unit

    (** Change the Jacobian-times-vector function.

        @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
        @cvode <node5#ss:jtimesfn> CVSpilsJacTimesVecFn *)
    val set_jac_times_vec_fn :
      ('a, 'k) session
      -> 'a jac_times_vec_fn
      -> unit

    (** Remove a Jacobian-times-vector function and use the default
        implementation.

        @cvode <node5#sss:optin_spils> CVSpilsSetJacTimesVecFn
        @cvode <node5#ss:jtimesfn> CVSpilsJacTimesVecFn *)
    val clear_jac_times_vec_fn : ('a, 'k) session -> unit

    (** Change the preconditioning direction without modifying
        callback functions. If the preconditioning type is changed from
        {{!Spils.preconditioning_type}Spils.PrecNone}
        then {!set_preconditioner} must be called to install the necessary
        callbacks.

        @cvode <node5#sss:optin_spils> CVSpilsSetPrecType *)
    val set_prec_type : ('a, 'k) session -> Spils.preconditioning_type -> unit
  end

(** Alternate Linear Solvers.

    @cvode <node8#s:new_linsolv> Providing Alternate Linear Solver Modules *)
module Alternate :
  sig
    (** Indicates problems during the solution of nonlinear equations at a
        step. Helps decide whether to update the Jacobian data kept by a
        linear solver. *)
    type conv_fail =
      | NoFailures
          (** Either the first call for a step or the local error test
              failed on the previous attempt at this setup but the Newton
              iteration converged. *)

      | FailBadJ
          (** The setup routine indicates that its Jacobian related data is not
              current and either the previous Newton corrector iteration did not
              converge, or, during the previous Newton corrector iteration, the
              linear solver's {!callbacks.lsolve} routine failed in a
              recoverable manner. *)

      | FailOther
          (** The previous Newton iteration failed to converge even
              though the linear solver was using current
              Jacobian-related data. *)

    (** Functions that initialize linear solver data, like counters and
        statistics.

        Raising any exception in this function (including
        {!Sundials.RecoverableFailure}) is treated as an unrecoverable error.

        @cvode <node8#SECTION00810000000000000000> linit *)
    type ('data, 'kind) linit = ('data, 'kind) session -> unit

    (** Functions that prepare the linear solver for subsequent calls to
        {!callbacks.lsolve}. The call
        [jcur = lsetup s convfail ypred fpred tmp] has as arguments

        - [s], the solver session,
        - [convfail], indicating any problem that occurred during the
           solution of the nonlinear equation on the current time step,
        - [ypred], the predicted [y] vector for the current internal
          step,
        - [fpred], the value of the right-hand side at [ypred], and,
        - [tmp], temporary variables for use by the routine.

        This function must return [true] only if the Jacobian-related data is
        current after the call. It may raise a {!Sundials.RecoverableFailure}
        exception to indicate that a recoverable error has occurred. Any other
        exception is treated as an unrecoverable error.

        {warning The vectors [ypred], [fpred], and those in [tmp] should not be
                 accessed after the function returns.}

        @cvode <node8#SECTION00820000000000000000> lsetup *)
    type ('data, 'kind) lsetup =
      ('data, 'kind) session
      -> conv_fail
      -> 'data
      -> 'data
      -> 'data triple
      -> bool

    (** Functions that solve the linear equation $Mx = b$.
        $M$ is a preconditioning matrix chosen by the user, and $b$ is the
        right-hand side vector calculated within the function.
        $M$ should approximate $I - \gamma J$ where
        {% $J = (\frac{\partial f}{\partial y})(t_n, y_{\text{cur}})$%},
        and $\gamma$ is available through {!get_gammas}.
        The call [lsolve s b weight ycur fcur] has as arguments:

        - [s], the solver session,
        - [b], for returning the calculated solution,
        - [weight], the error weights,
        - [ycur], the solver's current approximation to $y(t_n)$, and,
        - [fcur], a vector containing {% $f(t_n, y_{\text{cur}})$%}.

        Raising {!Sundials.RecoverableFailure} indicates a recoverable error.
        Any other exception is treated as an unrecoverable error.

        {warning The vectors [b], [weight], [ycur], and [fcur] should not be
                 accessed after the function returns.}

        @cvode <node8#SECTION00830000000000000000> lsolve
        @cvode <node3#e:Newtonmat> IVP solution (Eq. 2.5) *)
    type ('data, 'kind) lsolve =
      ('data, 'kind) session
      -> 'data
      -> 'data
      -> 'data
      -> 'data
      -> unit

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
    val make :
          (('data, 'kind) session
            -> ('data, 'kind) Nvector.t
            -> ('data, 'kind) callbacks)
          -> ('data, 'kind) linear_solver

    (** {3:internals Solver internals} *)

    (** Internal values used in Newton iteration.

        @cvode <node3#ss:ivp_sol> IVP Solution *)
    type gammas = {
      gamma : float;    (** The current $\gamma$ value. *)
      gammap : float;   (** The value of $\gamma$ at the last setup call. *)
    }

    (** Returns the current and previous gamma values. *)
    val get_gammas : ('data, 'kind) session -> gammas
  end

(** {2:tols Tolerances} *)

(** Functions that set the multiplicative error weights for use in the weighted
    RMS norm. The call [efun y ewt] takes the dependent variable vector [y] and
    fills the error-weight vector [ewt] with positive values or raises
    {!Sundials.NonPositiveEwt}. Other exceptions are eventually propagated, but
    should be avoided ([efun] is not allowed to abort the solver). *)
type 'data error_fun = 'data -> 'data -> unit

(** Tolerance specifications. *)
type ('data, 'kind) tolerance =
  | SStolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVtolerances of float * ('data, 'kind) Nvector.t
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFtolerances of 'data error_fun
    (** Set the multiplicative error weights for the weighted RMS norm. *)

(** A default relative tolerance of 1.0e-4 and absolute tolerance of 1.0e-8. *)
val default_tolerances : ('data, 'kind) tolerance

(** {2:solver Solver initialization and use} *)

(** Choice of method for solving non-linear systems that arise in solver
    formulas.

    @cvode <node3#ss:ivp_sol> IVP Solution
    @cvode <node5#sss:cvodemalloc> CVodeCreate *)
type ('a, 'kind) iter =
  | Newton of ('a, 'kind) linear_solver
    (** Newton iteration with a given linear solver. *)
  | Functional
    (** Functional iteration (non-stiff systems only). *)

(** Choice of linear multistep method.

    @cvode <node3#ss:ivp_sol> IVP Solution
    @cvode <node5#sss:cvodemalloc> CVodeCreate *)
type lmm =
  | Adams   (** Adams-Moulton formulas (non-stiff systems). *)
  | BDF     (** Backward Differentiation Formulas (stiff systems). *)

(** Right-hand side functions for calculating ODE derivatives. They are passed
    three arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
    - [dy], a vector for storing the value of $f(t, y)$.

    Within the function, raising a {!Sundials.RecoverableFailure} exception
    indicates a recoverable error. Any other exception is treated as an
    unrecoverable error.

    {warning [y] and [dy] should not be accessed after the function
             returns.}

    @cvode <node5#ss:rhsFn> CVRhsFn *)
type 'a rhsfn = float -> 'a -> 'a -> unit

(** Called by the solver to calculate the values of root functions. These
    ‘zero-crossings’ are used to detect significant events. The function is
    passed three arguments:
    - [t], the value of the independent variable, i.e., the simulation time,
    - [y], the vector of dependent-variable values, i.e., $y(t)$, and,
    - [gout], a vector for storing the value of $g(t, y)$.

    {warning [y] and [gout] should not be accessed after the function has
             returned.}

    @cvode <node5#ss:rootFn> cvRootFn *)
type 'a rootsfn = (float -> 'a -> RealArray.t -> unit)

(** Creates and initializes a session with the solver. The call
    {[init lmm iter tol f ~roots:(nroots, g) ~t0:t0 (neqs, y0)]} has
    as arguments:
    - [lmm],    the linear multistep method (see {!lmm}),
    - [iter],   either functional or Newton iteration (see {!iter}),
    - [tol],    the integration tolerances,
    - [f],      the ODE right-hand side function,
    - [nroots], the number of root functions,
    - [g],      the root function ([(nroots, g)] defaults to {!no_roots}),
    - [t0],     the initial value of the independent variable, and,
    - [y0],     a vector of initial values that also determines the number
                of equations.

    This function does everything necessary to initialize a session, i.e.,
    it makes the calls referenced below. The {!solve_normal} and
    {!solve_one_step} functions may be called directly.

    @cvode <node5#sss:cvodemalloc>   CVodeCreate/CVodeInit
    @cvode <node5#ss:cvrootinit>     CVodeRootInit
    @cvode <node5#sss:lin_solv_init> Linear solvers
    @cvode <node5#sss:cvtolerances>  CVodeSStolerances
    @cvode <node5#sss:cvtolerances>  CVodeSVtolerances
    @cvode <node5#sss:cvtolerances>  CVodeWFtolerances
    @cvode <node5#ss:ewtsetFn>       CVEwtFn *)
val init :
    lmm
    -> ('data, 'kind) iter
    -> ('data, 'kind) tolerance
    -> 'data rhsfn
    -> ?roots:(int * 'data rootsfn)
    -> float
    -> ('data, 'kind) Nvector.t
    -> ('data, 'kind) session

(** A convenience value for signalling that there are no roots to monitor. *)
val no_roots : (int * 'a rootsfn)

(** Values returned by the step functions. Failures are indicated by
    exceptions.

    @cvode <node5#sss:cvode> CVode *)
type solver_result =
  | Success             (** The solution was advanced. {cconst CV_SUCCESS} *)
  | RootsFound          (** A root was found. See {!get_root_info}.
                            {cconst CV_ROOT_RETURN} *)
  | StopTimeReached     (** The stop time was reached. See {!set_stop_time}.
                            {cconst CV_TSTOP_RETURN} *)

(** Integrates an ODE system over an interval. The call
    [tret, r = solve_normal s tout yout] has as arguments
    - [s], a solver session,
    - [tout], the next time at which a solution is desired, and,
    - [yout], a vector to store the computed solution.

    It returns [tret], the time reached by the solver, which will be equal to
    [tout] if no errors occur, and, [r], a {!solver_result}.

    @cvode <node5#sss:cvode> CVode (CV_NORMAL)
    @raise IllInput Missing or illegal solver inputs.
    @raise TooClose The initial and final times are too close to each other and not initial step size was specified.
    
    @raise TooMuchWork The requested time could not be reached in [mxstep] internal steps.
    @raise TooMuchAccuracy The requested accuracy could not be satisfied.
    @raise ErrFailure Too many error test failures within a step or at the minimum step size.
    @raise ConvergenceFailure Too many convergence test failures within a step or at the minimum step size.
    @raise LinearInitFailure Linear solver initialization failed.
    @raise LinearSetupFailure Linear solver setup failed unrecoverably.
    @raise LinearSolveFailure Linear solver solution failed unrecoverably.
    @raise RhsFuncFailure Unrecoverable failure in the RHS function [f].
    @raise FirstRhsFuncFailure Initial unrecoverable failure in the RHS function [f].
    @raise RepeatedRhsFuncFailure Too many convergence test failures, or unable to estimate the initial step size, due to repeated recoverable errors in the right-hand side function.
    @raise UnrecoverableRhsFuncFailure The right-hand side function had a recoverable error, but no recovery was possible. This error can only occur after an error test failure at order one.
    @raise RootFuncFailure Failure in the rootfinding function [g]. *)
val solve_normal : ('a, 'k) session -> float -> ('a, 'k) Nvector.t
                        -> float * solver_result

(** Like {!solve_normal} but returns after one internal solver step.

    @cvode <node5#sss:cvode> CVode (CV_ONE_STEP) *)
val solve_one_step : ('a, 'k) session -> float -> ('a, 'k) Nvector.t
                        -> float * solver_result

(** Returns the interpolated solution or derivatives.
    [get_dky s dky t k] computes the [k]th derivative of the function at time
    [t], i.e., {% $\frac{d^\mathtt{k}y(\mathtt{t})}{\mathit{dt}^\mathtt{k}}$%},
    and stores it in [dky]. The arguments must satisfy
    {% $t_n - h_u \leq \mathtt{t} \leq t_n$%}—where $t_n$
    denotes {!get_current_time} and $h_u$ denotes {!get_last_step},—
    and {% $0 \leq \mathtt{k} \leq q_u$%}—where $q_u$ denotes
    {!get_last_order}.

    This function may only be called after a successful return from either
    {!solve_normal} or {!solve_one_step}.

    @cvode <node5#sss:optin_root> CVodeGetDky
    @raise BadT [t] is not in the interval {% $[t_n - h_u, t_n]$%}.
    @raise BadK [k] is not in the range 0, 1, ..., $q_u$. *)
val get_dky : ('a, 'k) session -> ('a, 'k) Nvector.t -> float -> int -> unit

(** Reinitializes the solver with new parameters and state values. The
    values of the independent variable, i.e., the simulation time, and the
    state variables must be given.

    @param iter_type the iteration method (defaults to unchanged).
    @param roots new numbers of roots and the function that computes them (defaults to unchanged).
    @cvode <node5#sss:cvreinit> CVodeReInit *)
val reinit :
  ('a, 'kind) session
  -> ?iter_type:('a, 'kind) iter
  -> ?roots:(int * 'a rootsfn)
  -> float
  -> ('a, 'kind) Nvector.t
  -> unit

(** {2:set Modifying the solver (optional input functions)} *)

(** Sets the integration tolerances.

    @cvode <node5#sss:cvtolerances> CVodeSStolerances
    @cvode <node5#sss:cvtolerances> CVodeSVtolerances
    @cvode <node5#sss:cvtolerances> CVodeWFtolerances
    @cvode <node5#ss:ewtsetFn>       CVEwtFn *)
val set_tolerances : ('a, 'k) session -> ('a, 'k) tolerance -> unit

(** Opens the named file to receive messages from the default error handler.
    If the file already exists it is either truncated ([true]) or extended
    ([false]).
    The file is closed if the function is called again or when the session is
    garbage collected.

    @cvode <node5#sss:optin_main> CVodeSetErrFile *)
val set_error_file : ('a, 'k) session -> string -> bool -> unit

(** Specifies a custom function for handling error messages.
    This function must not fail: any exceptions are trapped and discarded.

    @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn
    @cvode <node5#ss:ehFn> CVErrHandlerFn *)
val set_err_handler_fn : ('a, 'k) session -> (error_details -> unit) -> unit

(** Restores the default error handling function.

    @cvode <node5#sss:optin_main> CVodeSetErrHandlerFn *)
val clear_err_handler_fn : ('a, 'k) session -> unit

(** Specifies the maximum order of the linear multistep method.

    @cvode <node5#sss:optin_main> CVodeSetMaxOrd *)
val set_max_ord : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of steps taken in attempting to reach
    a given output time.

    @cvode <node5#sss:optin_main> CVodeSetMaxNumSteps *)
val set_max_num_steps : ('a, 'k) session -> int -> unit

(** Specifies the maximum number of messages warning that [t + h = t] on
    the next internal step.

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

(** Limits the value of the independent variable [t] when solving.
    By default no stop time is imposed.

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

(** {2:get Querying the solver (optional output functions)} *)

(** Returns the real and integer workspace sizes.

    @cvode <node5#sss:optout_main> CVodeGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space          : ('a, 'k) session -> int * int

(** Returns the cumulative number of internal steps taken by the solver.

    @cvode <node5#sss:optout_main> CVodeGetNumSteps *)
val get_num_steps           : ('a, 'k) session -> int

(** Returns the number of calls to the right-hand side function.

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
val get_err_weights : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit

(** Returns the vector of estimated local errors.

    @cvode <node5#sss:optout_main> CVodeGetEstLocalErrors *)
val get_est_local_errors : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit

(** Summaries of integrator statistics. *)
type integrator_stats = {
    num_steps : int;
      (** Cumulative number of internal solver steps. *)
    num_rhs_evals : int;
      (** Number of calls to the right-hand side function. *)
    num_lin_solv_setups : int;
      (** Number of setups calls to the linear solver. *)
    num_err_test_fails : int;
      (** Number of local error test failures. *)
    last_order : int;
      (** Integration method order used in the last internal step. *)
    current_order : int;
      (** Integration method order to be used in the next internal step. *)
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

    @cvode <node5#sss:optout_main> CVodeGetIntegratorStats *)
val get_integrator_stats    : ('a, 'k) session -> integrator_stats

(** Prints the integrator statistics on the given channel.

    @cvode <node5#sss:optout_main> CVodeGetIntegratorStats *)
val print_integrator_stats  : ('a, 'k) session -> out_channel -> unit

(** Returns the number of nonlinear (functional or Newton) iterations performed.

    @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvIters *)
val get_num_nonlin_solv_iters : ('a, 'k) session -> int

(** Returns the number of nonlinear convergence failures that have occurred.

    @cvode <node5#sss:optout_main> CVodeGetNumNonlinSolvConvFails *)
val get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int

(** Returns both the numbers of nonlinear iterations performed [nniters] and
    nonlinear convergence failures [nncfails].

    @cvode <node5#sss:optout_main> CVodeGetNonlinSolvStats
    @return ([nniters], [nncfails]) *)
val get_nonlin_solv_stats : ('a, 'k) session -> int *int

(** {2:roots Additional root-finding functions} *)

(** [set_root_direction s dir] specifies the direction of zero-crossings to be
    located and returned. [dir] may contain one entry for each root function.

    @cvode <node5#sss:optin_root> CVodeSetRootDirection *)
val set_root_direction : ('a, 'k) session -> RootDirs.d array -> unit

(** Like {!set_root_direction} but specifies a single direction for all root
    functions.

    @cvode <node5#sss:optin_root> CVodeSetRootDirection *)
val set_all_root_directions : ('a, 'k) session -> RootDirs.d -> unit

(** Disables issuing a warning if some root function appears to be identically
    zero at the beginning of the integration.

    @cvode <node5#sss:optin_root> CVodeSetNoInactiveRootWarn *)
val set_no_inactive_root_warn : ('a, 'k) session -> unit

(** Returns the number of root functions. *)
val get_num_roots : ('a, 'k) session -> int

(** Fills an array showing which functions were found to have a root.

    @cvode <node5#sss:optout_root> CVodeGetRootInfo *)
val get_root_info : ('a, 'k) session -> Roots.t -> unit

(** Returns the cumulative number of calls made to the user-supplied root
    function g.

    @cvode <node5#sss:optout_root> CVodeGetNumGEvals *)
val get_num_g_evals : ('a, 'k) session -> int

(** {2:exceptions Exceptions} *)

(** Raised on missing or illegal solver inputs. Also raised if an element
    of the error weight vector becomes zero during time stepping, or the
    linear solver initialization function failed, or a root was found both at
    [t] and very near [t].
 
 @cvode <node5#sss:cvode> CV_ILL_INPUT *)
exception IllInput

(** The initial and final times are too close to each other and an initial step
    size was not specified.

    @cvode <node5#sss:cvode> CV_TOO_CLOSE *)
exception TooClose

(** The requested time could not be reached in [mxstep] internal steps.
    See {!set_max_num_steps}

    @cvode <node5#sss:cvode> CV_TOO_MUCH_WORK *)
exception TooMuchWork

(** The requested accuracy could not be satisfied.
 
    @cvode <node5#sss:cvode> CV_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** Too many error test failures within a step or at the minimum step size.
    See {!set_max_err_test_fails} and {!set_min_step}.

    @cvode <node5#sss:cvode> CV_ERR_FAILURE *)
exception ErrFailure                

(** Too many convergence test failures within a step or at the minimum step
    size. See {!set_max_conv_fails} and {!set_min_step}.
 
    @cvode <node5#sss:cvode> CV_CONV_FAILURE *)
exception ConvergenceFailure

(** Linear solver initialization failed.
 
    @cvode <node5#sss:cvode> CV_LINIT_FAIL *)
exception LinearInitFailure         

(** Linear solver setup failed in an unrecoverable manner.

    @cvode <node5#sss:cvode> CV_LSETUP_FAIL *)
exception LinearSetupFailure        

(** Linear solver solution failed in an unrecoverable manner.

    @cvode <node5#sss:cvode> CV_LSOLVE_FAIL *)
exception LinearSolveFailure        

(** The right-hand side function failed in an unrecoverable manner.
  
    @cvode <node5#sss:cvode> CV_RHSFUNC_FAIL *)
exception RhsFuncFailure

(** The right-hand side function had a recoverable error when first called.

    @cvode <node5#sss:cvode> CV_FIRST_RHSFUNC_ERR *)
exception FirstRhsFuncFailure

(** Too many convergence test failures, or unable to estimate the initial step
    size, due to repeated recoverable errors in the right-hand side function.

    @cvode <node5#sss:cvode> CV_REPTD_RHSFUNC_ERR *)
exception RepeatedRhsFuncFailure

(** The right-hand side function had a recoverable error, but no recovery was
    possible. This error can only occur after an error test failure at order
    one.
  
    @cvode <node5#sss:cvode> CV_UNREC_RHSFUNC_ERR *)
exception UnrecoverableRhsFuncFailure

(** The rootfinding function failed.
  
    @cvode <node5#sss:cvode> CV_RTFUNC_FAIL *)
exception RootFuncFailure           

(** Raised by {!get_dky} for invalid order values.

    @cvode <node5#ss:optional_dky> CVodeGetDky (CV_BAD_K) *)
exception BadK

(** Raised by {!get_dky} for invalid time values.
 
    @cvode <node5#ss:optional_dky> CVodeGetDky (CV_BAD_T) *)
exception BadT

