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
(*               User Documentation for IDA v2.7.0                     *)
(*         Alan C. Hindmarsh, Radu Serban, and Aaron Collier           *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Vector-independent types and values for the IDA solver.

 @version VERSION()
 @author Timothy Bourke (Inria)
 @author Jun Inoue (Inria)
 @author Marc Pouzet (LIENS)
 *)

include module type of Sundials
  with type Roots.t = Sundials.Roots.t
  and type Roots.root_event = Sundials.Roots.root_event
  and type RootDirs.t = Sundials.RootDirs.t
  and type RootDirs.root_direction = Sundials.RootDirs.root_direction
  and type solver_result = Sundials.solver_result
  and type error_details = Sundials.error_details

(** {2 General} *)

(** {3 Solver initialisation} *)

(**
 Values for root directions.
 @ida <node5#sss:optin_root> IDASetRootDirection
 *)
type root_direction = RootDirs.root_direction

(**
 This is a convenience value for signalling to {!Ida_serial.init} and
 {!Ida_nvector.init} that there are no roots (zero-crossings) to
 monitor.
 *)
val no_roots : (int * ('a -> 'b -> 'c -> 'd -> unit))

(** {3 Statistics} *)

(**
 Aggregated integrator statistics.
 @ida <node5#sss:optout_main> IDAGetIntegratorStats
 *)
type integrator_stats = {
    num_steps : int;
    num_res_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

(** {3:exceptions Exceptions} *)

(** @ida <node5#sss:ida> IDA_ILL_INPUT *)
exception IllInput

(** @ida <node5#sss:ida> IDA_TOO_CLOSE *)
exception TooClose

(** @ida <node5#sss:ida> IDA_TOO_MUCH_WORK *)
exception TooMuchWork

(** @ida <node5#sss:ida> IDA_TOO_MUCH_ACC *)
exception TooMuchAccuracy

(** @ida <node5#sss:ida> IDA_ERR_FAIL *)
exception ErrFailure                

(** @ida <node5#sss:ida> IDA_CONV_FAIL *)
exception ConvergenceFailure        

(** @ida <node5#sss:ida> IDA_LINIT_FAIL *)
exception LinearInitFailure         

(** @ida <node5#sss:ida> IDA_LSETUP_FAIL *)
exception LinearSetupFailure        

(** @ida <node5#sss:ida> IDA_LSOLVE_FAIL *)
exception LinearSolveFailure        

(** @ida <node5#sss:ida> IDA_RES_FAIL *)
exception ResFuncFailure

(** @ida <node5#sss:ida> IDA_FIRST_RES_FAIL *)
exception FirstResFuncFailure       

(** @ida <node5#sss:ida> IDA_REP_RES_ERR *)
exception RepeatedResFuncErr        

(** @ida <node5#sss:ida> IDA_UNREC_RESFUNC_ERR *)
exception UnrecoverableResFuncErr   

(** @ida <node5#sss:ida> IDA_RTFUNC_FAIL *)
exception RootFuncFailure           

exception BadK      (** k is not in the range 0, 1, ..., q_u (IDA_BAD_K)
                        @ida <node5#ss:optional_dky> IDAGetDky *)

exception BadT      (** t is not in the interval
                        \[t_n - h_u, t_n\] (IDA_BAD_T)
                        @ida <node5#ss:optional_dky> IDAGetDky *)
exception BadDky    (** invalid dky argument (IDA_BAD_DKY)
                        @ida <node5#ss:optional_dky> IDAGetDky *)

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

