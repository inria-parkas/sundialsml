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
(*               User Documentation for CVODES v2.7.0                  *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Abstract nvector interface to the CVODES Extensions.

  @version VERSION()
  @author Timothy Bourke (Inria)
  @author Jun Inoue (Inria)
  @author Marc Pouzet (LIENS)
 *)

(** {2:quadrature Integration of pure quadrature equations} *)

module Quadrature :
  sig

(* XXX WORKING XXX

   CVodeQuadInit
     fQ   : float -> 'a -> 'a -> unit (CVQuadRhsFn)
     yQ0  : 'a nvector

      @cvodes <node5#ss:quad_malloc> CVodeQuadInit
      @cvodes <node5#ss:user_fct_quad> CVQuadRhsFn
    

   CVodeQuadReInit
     yQ0  : 'a nvector
     (check for CV_NO_QUAD if cvode_mem not QuadInited... )
      @cvodes <node5#ss:quad_malloc> CVodeQuadReInit

  Design questions:
    - How to incorporate extra elements into the session type.
      (i.e., qrhs function )
    - How to mark the session type as valid for the quadrature functions?
      1. Fail dynamically (add an exception for CV_NO_QUAD).
      2. Add a phantom type.
      3. create a quad_session from which a normal session can be extracted?
      4. upgrade a normal session to a quad_session?
 *)

    (** {3:exceptions Exceptions} *)

    (** @cvodes <node5#SECTION00572000000000000000> CV_QRHSFUNC_FAIL *)
    exception QuadRhsFuncFailure

    (** @cvodes <node5#SECTION00572000000000000000> CV_FIRST_QRHSFUNC_ERR *)
    exception FirstQuadRhsFuncErr

    (** @cvodes <node5#SECTION00572000000000000000> CV_REPTD_QRHSFUNC_ERR *)
    exception RepeatedQuadRhsFuncErr

    (** @cvodes <node5#SECTION00572000000000000000> CV_UNREC_QRHSFUNC_ERR *)
    exception UnrecoverableQuadRhsFuncErr

    (**
      [tret = get_quad s yq] fills [yq] with the quadrature solution vector
      after a successful return from {!Cvode.solve_normal} or
      {!Cvode.solve_one_step}, and returns the time reached by the solver.

      @cvodes <node5#ss:quad_get> CVodeGetQuad
     *)
    val get : 'a session -> 'a nvector -> float

    (**
      [tret = get_quad_dky s t k dkyq] fills [dkyq] with the derivatives of the
      quadrature solution vector after a successful return from
      {!Cvode.solve_normal} or {!Cvode.solve_one_step}. The time requested, [t],
      must fall within the interval defined by the last successful step
      ({!Cvode.get_last_step}). The requested order, [k], must be less than or
      equal to the value returned by {!Cvode.get_last_order}.

      @cvodes <node5#ss:quad_get> CVodeGetQuadDky
     *)
    val get_dky : 'a session -> float -> int -> 'a nvector -> unit


    (** {3 Tolerance specification} *)

    (**
        Specify that quadrature variables should be used in the step size
        control mechanism. [ss_tolerances s reltol abstol] sets the relative and
        absolute tolerances using scalar values.

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSStolerances
     *)
    val ss_tolerances : session -> float -> float -> unit

    (**
        Specify that quadrature variables should be used in the step size
        control mechanism. [sv_tolerances s reltol abstol] sets the relative
        tolerance using a scalar value, and the absolute tolerance as a vector.

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
        @cvodes <node5#ss:quad_optional_input> CVodeQuadSVtolerances
     *)
    val sv_tolerances : session -> float -> 'a vector -> unit

    (**
        Specify that quadrature variables should not be used in the step size
        control mechanism (this is the default). 

        @cvodes <node5#ss:quad_optional_input> CVodeSetQuadErrCon
     *)
    val no_tolerances : 'a session -> unit

    (** {3 Optional input functions} *)

    (**
      Returns the number of calls to the user's quadrature right-hand side function.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumRhsEvals
     *)
    val get_num_rhs_evals       : 'a session -> int

    (**
      Returns the number of local error test failures due to quadrature variables.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadNumErrTestFails
     *)
    val get_num_err_test_fails  : 'a session -> int

    (**
      Returns the quadrature error weights at the current time.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadErrWeights
     *)
    val get_quad_err_weights : 'a session -> 'a nvector -> unit

    (**
      [nfqevals, nqetfails = get_stats s] returns
      - [fqevals], the number of calls to the user's quadrature function, and,
      - [nqetfails], the number of error test failures due to quadrature variables.

      @cvodes <node5#ss:quad_optional_output> CVodeGetQuadStats
     *)
    val get_stats : session -> int * int

  end

