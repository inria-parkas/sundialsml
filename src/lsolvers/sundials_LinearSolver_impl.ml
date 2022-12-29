(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

(* Hack to ensure that Sundials.c_init_module is executed so that the global
   exceptions are properly registered. *)
let e = Sundials.RecoverableFailure

(* Types underlying LinearSolver.  *)

(* This module defines the lsolver types which are manipulated (abstractly)
   via the LinearSolver module. They are declared here so that they can also be
   used (concretely) from the <solver> modules. To ensure that this type will
   be opaque outside of Sundials/ML, we simply do not install the
   lsolver_impl.cmi file. *)

(* Must correspond with sundials_linearSolver_ml.h:lsolver_solver_type *)
type linear_solver_type =
  | Direct
  | Iterative
  | MatrixIterative
  | MatrixEmbedded

(* Must correspond with sundials_linearSolver_ml.h:lsolver_solver_id_tag *)
type linear_solver_id =
  | Band
  | Dense
  | Klu
  | LapackBand
  | LapackDense
  | Pcg
  | Spbcgs
  | Spfgmr
  | Spgmr
  | Sptfqmr
  | Superludist
  | Superlumt
  | CuSolverSp_batchQR
  | MagmaDense
  | OneMKLDense
  | Gingko
  | KokkosDense
  | Custom

module Klu = struct (* {{{ *)

  (* Must correspond with linearSolver_ml.h:lsolver_klu_ordering_tag *)
  type ordering =
    Amd
  | ColAmd
  | Natural

  type info = {
    mutable ordering : ordering option;

    mutable reinit : int -> int option -> unit;
    mutable set_ordering : ordering -> unit;
  }

  let info () = {
    ordering     = None;
    reinit       = (fun _ _ -> ());
    set_ordering = (fun _ -> ());
  }

end (* }}} *)

module Superlumt = struct (* {{{ *)

  (* Must correspond with linearSolver_ml.h:lsolver_superlumt_ordering_tag *)
  type ordering =
    Natural
  | MinDegreeProd
  | MinDegreeSum
  | ColAmd

  type info = {
    mutable ordering     : ordering option;
    mutable set_ordering : ordering -> unit;
    num_threads          : int;
  }

  let info num_threads = {
    ordering     = None;
    set_ordering = (fun _ -> ());
    num_threads  = num_threads;
  }

end (* }}} *)

module Iterative = struct (* {{{ *)

  (* Must correspond with linearSolver_ml.h:lsolver_gramschmidt_type_tag *)
  type gramschmidt_type =
    | ModifiedGS
    | ClassicalGS

  (* Must correspond with linearSolver_ml.h:lsolver_preconditioning_type_tag *)
  type preconditioning_type =
    | PrecNone
    | PrecLeft
    | PrecRight
    | PrecBoth

  type info = {
    mutable maxl             : int;
    mutable gs_type          : gramschmidt_type option;
    mutable max_restarts     : int option;

    mutable set_maxl         : int -> unit;
    mutable set_gs_type      : gramschmidt_type -> unit;
    mutable set_max_restarts : int -> unit;
    mutable set_prec_type    : preconditioning_type -> unit;
  }

  let info = {
    maxl             = 0;
    gs_type          = None;
    max_restarts     = None;

    set_maxl         = (fun _ -> ());
    set_gs_type      = (fun _ -> ());
    set_max_restarts = (fun _ -> ());
    set_prec_type    = (fun _ -> ());
  }

end (* }}} *)

module Custom = struct (* {{{ *)

  type ('data, 'kind) atimes_with_data

  external call_atimes
    : ('data, 'kind) atimes_with_data
      -> ('data, 'kind) Nvector.t
      -> ('data, 'kind) Nvector.t
      -> unit
    = "sunml_lsolver_call_atimes"

  type ('data, 'kind) precond_with_data

  external call_psetup : ('data, 'kind) precond_with_data -> unit
    = "sunml_lsolver_call_psetup"

  external call_psolve
    : ('data, 'kind) precond_with_data
      -> ('data, 'kind) Nvector.t
      -> ('data, 'kind) Nvector.t
      -> float
      -> bool
      -> unit
    = "sunml_lsolver_call_psolve"

  (* The fields and their order must match linearSolver_ml.h:lsolver_ops_index *)
  type ('matrix, 'data, 'kind, 't) ops = {
    init : 't -> unit;

    setup : 't -> 'matrix -> unit;

    solve : 't -> 'matrix -> 'data -> 'data -> float -> unit;

    set_atimes : 't -> ('data, 'kind) atimes_with_data -> unit;

    set_preconditioner :
      't -> ('data, 'kind) precond_with_data -> bool -> bool -> unit;

    set_scaling_vectors : 't -> 'data option -> 'data option -> unit;

    set_zero_guess : 't -> bool -> unit;

    get_id : 't -> linear_solver_id;

    get_num_iters : 't -> int;

    get_res_norm : 't -> float;

    get_res_id : 't -> ('data, 'kind) Nvector.t;

    get_last_flag : 't -> int;

    get_work_space : 't -> int * int;

    set_prec_type : 't -> Iterative.preconditioning_type -> unit;
  }

  (* The fields and their order must match linearSolver_ml.h:lsolver_hasops_index *)
  type has_ops = {
    has_init                : bool;
    has_setup               : bool;
    has_set_atimes          : bool;
    has_set_preconditioner  : bool;
    has_set_scaling_vectors : bool;
    has_set_zero_guess      : bool;
    has_get_num_iters       : bool;
    has_get_res_norm        : bool;
    has_get_res_id          : bool;
    has_get_last_flag       : bool;
    has_get_work_space      : bool;
  }

end (* }}} *)

exception LinearSolverInUse

type ('m, 'nd, 'nk) cptr

(* Must correspond with linearSolver_ml.h:lsolver_solver_data_tag *)
type (_, 'nd, 'nk, _) solver_data =
  (* Iterative Linear Solvers *)
  | Spbcgs      : ('m, 'nd, 'nk, [>`Spbcgs])  solver_data
  | Spfgmr      : ('m, 'nd, 'nk, [>`Spfgmr])  solver_data
  | Spgmr       : ('m, 'nd, 'nk, [>`Spgmr])   solver_data
  | Sptfqmr     : ('m, 'nd, 'nk, [>`Sptfqmr]) solver_data
  | Pcg         : ('m, 'nd, 'nk, [>`Pcg])     solver_data
  (* Direct Linear Solvers *)
  | Dense       : (Matrix.Dense.t, 'nd, 'nk, [>`Dls]) solver_data
  | LapackDense : (Matrix.Dense.t, 'nd, 'nk, [>`Dls]) solver_data
  | Band        : (Matrix.Band.t,  'nd, 'nk, [>`Dls]) solver_data
  | LapackBand  : (Matrix.Band.t,  'nd, 'nk, [>`Dls]) solver_data
  | Klu         : Klu.info
                  -> ('s Matrix.Sparse.t, 'nd, 'nk, [>`Klu]) solver_data
  | Superlumt   : Superlumt.info
                  -> ('s Matrix.Sparse.t, 'nd, 'nk, [>`Slu]) solver_data
  (* Custom Linear Solvers *)
  | Custom      : ('t * ('m, 'nd, 'nk, 't) Custom.ops)
                   -> ('m, 'nd, 'nk, [>`Custom of 't]) solver_data

type 'd atimesfn = 'd -> 'd -> unit
type psetupfn = unit -> unit
type 'd psolvefn = 'd -> 'd -> float -> bool -> unit

(* Must correspond with linearSolver_ml.h:lsolver_ocaml_callbacks_index *)
type ('d, 'k) ocaml_callbacks = {
  mutable ocaml_atimes : 'd atimesfn;
  mutable ocaml_psetup : psetupfn;
  mutable ocaml_psolve : 'd psolvefn;
  mutable scaling_vector1 : ('d, 'k) Nvector.t option;
  mutable scaling_vector2 : ('d, 'k) Nvector.t option;
}

(* The type t is defined in two parts, record and constructor, for
   compatiblity with older versions of OCaml.

   When using a Direct solver in Sundials, a Jacobian matrix is passed
   twice: once to create the generic solver, where the matrix is only
   used for compatibility checks, and once when creating the
   session-specific solver, where the matrix is kept and used internally
   for calculations and cloning.

   We store the Jacobian matrix in direct solver (i.e., here) for
   two reasons:
   1. To save users the trouble of having to pass it twice.
   2. As a convenient way to ensure that the underlying SUNMatrix is not
      garbage collected while being used by a linear solver.

   Garbage collection could also be avoided by holding a reference to the
   Jacobian matrix in the session object. This requires, however, treating
   the matrix kind type argument ('mk) in some way (since we must stop the
   SUNMatrix and not just its content from being GCed). Normally, the
   session interface is unconcerned by matrix kind type arguments since
   the compatibility with matrix implementations is checked for the generic
   linear solver (which may directly manipulates the matrix data).
 *)
type ('m, 'mk, 'nd, 'nk, 't) linear_solver_data = {
  rawptr : ('m, 'nd, 'nk) cptr;
  solver : ('m, 'nd, 'nk, 't) solver_data;
  matrix : ('mk, 'm, 'nd, 'nk) Matrix.t option;
  compat : Iterative.info;
  context : Sundials.Context.t;
  mutable check_prec_type : Iterative.preconditioning_type -> bool;
  ocaml_callbacks : ('nd, 'nk) ocaml_callbacks Sundials_impl.Vptr.vptr;
  mutable info_file : Logfile.t option;
  mutable attached : bool;
}

let empty_ocaml_callbacks () =
  Sundials_impl.Vptr.make ({
    ocaml_atimes = (fun _ _     -> failwith "internal error: ocaml_atimes");
    ocaml_psetup = (fun _       -> failwith "internal error: ocaml_psetup");
    ocaml_psolve = (fun _ _ _ _ -> failwith "internal error: ocaml_psolve");
    scaling_vector1 = None;
    scaling_vector2 = None;
  })

(* Slight complication to hide matrix kind type argument as existential.
   This type argument is enforced by the module interface -- some
   solver creation functions require it to be Matrix.standard since the
   underlying Sundials code access the matrix data directly and not just
   through the SUNMatrix interface -- but we prefer to hide this detail
   from users (Direct.t already has four type arguments!) *)
type ('m, 'nd, 'nk, 't) linear_solver =
  LS : ('m, 'mk, 'nd, 'nk, 't) linear_solver_data
       -> ('m, 'nd, 'nk, 't) linear_solver

(* Wrapper for a linear solver held within a session to stop its elements
   from being garbage collected.  *)
type held_linear_solver =
  | HLS : ('m, 'mk, 'nd, 'nk, 't) linear_solver_data
          -> held_linear_solver
  | NoHLS : held_linear_solver

type ('nd, 'nk) solver =
  S : ('m, 'nd, 'nk, 't) solver_data -> ('nd, 'nk) solver

let attach (LS ({ attached } as s)) =
  if attached then raise LinearSolverInUse;
  s.attached <- true

external c_set_prec_type
  : ('m, 'nd, 'nk) cptr
    -> ('m, 'nd, 'nk, 't) solver_data
    -> Iterative.preconditioning_type
    -> bool
    -> unit
  = "sunml_lsolver_set_prec_type"

let impl_set_prec_type (type t) rawptr solver prec_type docheck =
  match (solver : ('m, 'nd, 'nk, t) solver_data) with
  | Custom (ldata, { Custom.set_prec_type = f }) -> f ldata prec_type
  | _ -> c_set_prec_type rawptr solver prec_type docheck

external c_make_custom
  : linear_solver_type
    -> ('t * ('m, 'nd, 'nk, 't) Custom.ops) Weak.t
    -> Custom.has_ops
    -> Sundials.Context.t
    -> ('m, 'nd, 'nk) cptr
  = "sunml_lsolver_make_custom"

