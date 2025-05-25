(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

(* Hack to ensure that Sundials.c_init_module is executed so that the global
   exceptions are properly registered. *)
let e = Sundials.RecoverableFailure

(* Types shared between Arkode, Arkode_bbd, Arkode_klu, and Arkode_superlumt. *)

(* This module's purpose is just to define types shared between Arkode,
   Arkode_bbd, Arkode_klu, and Arkode_superlumt. The session and
   session_linear_solver types must have their implementation exposed to all
   modules yet be abstract to code that lives outside of the
   sundialsml library.

   To satisfy this requirement, we define the types in this separate
   module that exports everything, then omit Arkode_impl.cmi at installation.
   This way, the compiled library doesn't divulge implementation details.

   Unfortunately, session depends on many of the other data types used in
   those modules, so we are forced to declare them here.
   To avoid repeating the types' definitions in Arkode, Arkode_bbd, Arkode_klu,
   and Arkode_superlumt, we "include Arkode_impl" from each module. *)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in arkode_ml.h (and code in arkode_ml.c) must also be updated.
 *)

module LSI = Sundials_LinearSolver_impl
module NLSI = Sundials_NonlinearSolver_impl

type 'a triple = 'a * 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;
    jac_y   : 'a;
    jac_fy  : 'a;
    jac_tmp : 't
  }

module DirectTypes = struct
  type 'm jac_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg
    -> 'm
    -> unit

  (* These fields are accessed from arkode_ml.c *)
  type 'm jac_callback =
    {
      jacfn: 'm jac_fn;
      mutable jmat : 'm option (* Not used in Sundials >= 3.0.0 *)
    }

  type 'm linsys_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg
    -> 'm
    -> 'm option
    -> bool
    -> float
    -> bool

  let no_callback = fun _ _ -> Sundials_impl.crash "no direct callback"

  let no_linsysfn = fun _ _ _ _ _ ->
    Sundials_impl.crash "no linear system function callback"

end

module SpilsCommonTypes = struct
  (* Types that don't depend on jacobian_arg.  *)

  type 'a prec_solve_arg =
    {
      rhs   : 'a;
      gamma : float;
      delta : float;
      left  : bool;
    }

end

module SpilsTypes' = struct
  include SpilsCommonTypes

  type 'a prec_solve_fn =
    (unit, 'a) jacobian_arg
    -> 'a prec_solve_arg
    -> 'a
    -> unit

  type 'a prec_setup_fn =
    (unit, 'a) jacobian_arg
    -> bool
    -> float
    -> bool

  type 'd jac_times_setup_fn =
    (unit, 'd) jacobian_arg
    -> unit

  type 'a jac_times_vec_fn =
    ('a, 'a) jacobian_arg
    -> 'a (* v *)
    -> 'a (* Jv *)
    -> unit

  type 'a precfns =
    {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }

end

module MassTypes' = struct

  module Direct' = struct
    type 'm mass_fn =
      float
      -> RealArray.t triple
      -> 'm
      -> unit

    (* These fields are accessed from arkode_ml.c *)
    type 'm mass_callback =
      {
        massfn: 'm mass_fn;
        mutable mmat : 'm option (* Not used in Sundials >= 3.0.0 *)
      }

    let no_mass_callback = fun _ _ -> Sundials_impl.crash "no mass callback"
  end

  module Iterative' = struct

    type 'd prec_solve_arg =
      {
        rhs   : 'd;
        delta : float;
        left  : bool;
      }

    type 'd prec_solve_fn =
         float
      -> 'd prec_solve_arg
      -> 'd
      -> unit

    type 'd prec_setup_fn = float -> unit

    type mass_times_setup_fn = float -> unit

    type 'd mass_times_vec_fn =
         float
      -> 'd
      -> 'd
      -> unit

    type 'a precfns =
      {
        prec_solve_fn : 'a prec_solve_fn;
        prec_setup_fn : 'a prec_setup_fn option;
      }
  end
end

module ArkodeBbdParamTypes = struct
  type 'a local_fn = float -> 'a -> 'a -> unit
  type 'a comm_fn = float -> 'a -> unit
  type 'a precfns =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end

module ArkodeBbdTypes = struct
  type bandwidths =
    {
      mudq    : int;
      mldq    : int;
      mukeep  : int;
      mlkeep  : int;
    }
end

type 'step arkode_mem
type c_weak_ref

type arkstep  = [`ARKStep]
type erkstep  = [`ERKStep]
type sprkstep = [`SPRKStep]
type mristep  = [`MRIStep]

module Global = struct

  type 'a rhsfn = float -> 'a -> 'a -> unit
  type 'a rootsfn = float -> 'a -> RealArray.t -> unit
  type error_handler = Util.error_details -> unit
  type 'a error_weight_fun = 'a -> 'a -> unit

  type adaptivity_args = {
      h1 : float;
      h2 : float;
      h3 : float;
      e1 : float;
      e2 : float;
      e3 : float;
      q  : int;
      p  : int;
    }

  type 'd adaptivity_fn = float -> 'd -> adaptivity_args -> float
  type 'd stability_fn = float -> 'd -> float
  type 'd resize_fn = 'd -> 'd -> unit
  type 'd postprocess_fn = float -> 'd -> unit
  type 'd stage_predict_fn = float -> 'd -> unit
  type 'd pre_inner_fn = float -> 'd array -> unit
  type 'd post_inner_fn = float -> 'd -> unit
end

(* Inner steppers *)

type ('d, 'k) inner_stepper_cptr

(* synchronized with arkode_ml.h: arkode_mri_istepper_full_rhsfn_mode_tag *)
type fullrhs_mode =
  | Start
  | End
  | Other

(* synchronized with arkode_ml.h: arkode_mri_istepper_callbacks_index *)
type 'd inner_stepper_callbacks = {
  evolve_fn   : float -> float -> 'd -> unit;
  full_rhs_fn : float -> 'd -> 'd -> fullrhs_mode -> unit;
  reset_fn    : float -> 'd -> unit;
}

open Global

type 'a res_weight_fun = 'a -> 'a -> unit

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.

   The ls_solver and mass_solver fields only exist to ensure that the linear
   solver is not garbage collected while still being used by a session.
*)

type ('a, 'kind, 'step) session = {
  arkode     : 'step arkode_mem;
  backref    : c_weak_ref;
  nroots     : int;
  mutable checkvec     : (('a, 'kind) Nvector.t -> unit);
  mutable uses_resv    : bool;
  context    : Context.t;

  mutable exn_temp     : exn option;

  mutable problem      : problem_type; (* ARK only *)
  mutable rhsfn1       : 'a rhsfn;  (* ARK: implicit; ERK: f; MRI: slow *)
  mutable rhsfn2       : 'a rhsfn;  (* ARK: explicit; ERK: unused; MRI: unused *)

  mutable rootsfn      : 'a rootsfn;
  mutable errh         : error_handler;
  mutable errw         : 'a error_weight_fun;
  mutable resw         : 'a res_weight_fun;  (* ARK only *)

  mutable error_file   : Logfile.t option;
  mutable diag_file    : Logfile.t option;

  mutable adaptfn      : 'a adaptivity_fn;
  mutable stabfn       : 'a stability_fn;
  mutable resizefn     : 'a resize_fn;
  mutable poststepfn   : 'a postprocess_fn;
  mutable poststagefn  : 'a postprocess_fn;
  mutable stagepredictfn : 'a stage_predict_fn (* ARK only *);
  mutable preinnerfn   : 'a pre_inner_fn;  (* MRI only *)
  mutable postinnerfn  : 'a post_inner_fn; (* MRI only *)
  mutable preinnerarray : 'a array;        (* MRI only *)

  (* ARK only *)
  mutable linsolver      : ('a, 'kind, 'step) lin_solver option;
  mutable ls_solver      : LSI.held_linear_solver;
  mutable ls_callbacks   : ('a, 'kind) linsolv_callbacks;
  mutable ls_precfns     : 'a linsolv_precfns;

  (* ARK only *)
  mutable mass_solver    : LSI.held_linear_solver;
  mutable mass_callbacks : ('a, 'kind) mass_callbacks;
  mutable mass_precfns   : 'a mass_precfns;

  (* ARK and MRI only *)
  mutable nls_solver     : ('a, 'kind, ('a, 'kind, 'step) session, [`Nvec])
                           NLSI.nonlinear_solver option;
  mutable nls_rhsfn      : 'a rhsfn;

  (* MRI only *)
  mutable inner_session  : ('a, 'kind) inner_stepper option;
}

and problem_type =
  | ImplicitOnly
  | ExplicitOnly
  | ImplicitAndExplicit

and ('data, 'kind, 'step) lin_solver =
  ('data, 'kind, 'step) session
  -> ('data, 'kind) Nvector.t
  -> unit

(* Note: When compatibility with Sundials < 3.0.0 is no longer required,
         this type can be greatly simplified since we would no longer
         need to distinguish between different "direct" linear solvers.

   Note: The first field must always hold the callback closure
         (it is accessed as Field(cb, 0) from arkode_ml.c.
         The second argument holds a reference to the Jacobian matrix,
         whose underlying data is used within the solver (Sundials >= 3.0.0),
         to prevent its garbage collection.

   Normally, for the matrix arguments of DlsDenseCallback, etc., we should
   impose that 'kind satisfy [>Nvector_serial.kind], but this constraint
   then propagates to the session type where it is unwanted (otherwise we
   could not create sessions that don't care about the kind of nvector
   expected in matvec operations, i.e., those without a linear solver or
   that do not require and explicit jacobian matrix). Instead, the constraint
   is imposed indirectly by the external interface (i.e., creating a linear
   solver requires complete compatibility between nvectors, matrices, and
   sessions.)
*)
and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  (* Dls *)
  | DlsDenseCallback
      of Matrix.Dense.t DirectTypes.jac_callback
         * Matrix.Dense.t DirectTypes.linsys_fn
  | DlsBandCallback
      of Matrix.Band.t  DirectTypes.jac_callback
         * Matrix.Band.t DirectTypes.linsys_fn

  (* Sls *)
  | SlsKluCallback
      : ('s Matrix.Sparse.t) DirectTypes.jac_callback
        * ('s Matrix.Sparse.t) DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks
  | SlsSuperlumtCallback
      : ('s Matrix.Sparse.t) DirectTypes.jac_callback
        * ('s Matrix.Sparse.t) DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks

  (* Custom *)
  | DirectCustomCallback
      : 'm DirectTypes.jac_callback
        * 'm DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks

  (* Spils *)
  | SpilsCallback1 of 'a SpilsTypes'.jac_times_vec_fn option
                      * 'a SpilsTypes'.jac_times_setup_fn option
  | SpilsCallback2 of 'a rhsfn

and 'a linsolv_precfns =
  | NoPrecFns
  | PrecFns of 'a SpilsTypes'.precfns
  | BandedPrecFns
  | BBDPrecFns of 'a ArkodeBbdParamTypes.precfns

(* Note: When compatibility with Sundials < 3.0.0 is no longer required,
         this type can be greatly simplified since we would no longer
         need to distinguish between different "direct" linear solvers.

   Note: The first field must always hold the callback closure
         (it is accessed as Field(cb, 0) from arkode_ml.c.
         The second argument holds a reference to the Jacobian matrix,
         whose underlying data is used within the solver (Sundials >= 3.0.0),
         to prevent its garbage collection.
*)
and ('a, 'kind) mass_callbacks =
  | NoMassCallbacks

  (* Dls *)
  | DlsDenseMassCallback
      of Matrix.Dense.t MassTypes'.Direct'.mass_callback * Matrix.Dense.t
  | DlsBandMassCallback
      of Matrix.Band.t  MassTypes'.Direct'.mass_callback * Matrix.Band.t

  (* Sls *)
  | SlsKluMassCallback
      : ('s Matrix.Sparse.t) MassTypes'.Direct'.mass_callback * 's Matrix.Sparse.t
        -> ('a, 'kind) mass_callbacks
  | SlsSuperlumtMassCallback
      : ('s Matrix.Sparse.t) MassTypes'.Direct'.mass_callback * 's Matrix.Sparse.t
        -> ('a, 'kind) mass_callbacks

  (* Custom *)
  | DirectCustomMassCallback :
      'm MassTypes'.Direct'.mass_callback * 'm -> ('a, 'kind) mass_callbacks

  (* Spils *)
  | SpilsMassCallback of 'a MassTypes'.Iterative'.mass_times_vec_fn
                       * MassTypes'.Iterative'.mass_times_setup_fn option

and 'a mass_precfns =
  | NoMassPrecFns
  | MassPrecFns of 'a MassTypes'.Iterative'.precfns

(* synchronized with arkode_ml.h: arkode_mri_istepper_tag *)
and ('d, 'k) istepper =
  | ARKStepInnerStepper of ('d, 'k, arkstep) session (* Sundials < 580 *)
  | SundialsInnerStepper
  | CustomInnerStepper of 'd inner_stepper_callbacks

(* synchronized with arkode_ml.h: arkode_mri_istepper_index *)
and ('d, 'k) inner_stepper = {
  rawptr   : ('d, 'k) inner_stepper_cptr;
  istepper : ('d, 'k) istepper;
  mutable icheckvec : (('d, 'k) Nvector.t -> unit) option
}

(* Linear solver check functions *)

let ls_check_direct session =
  if Sundials_configuration.safe then
    match session.ls_callbacks with
    | DlsDenseCallback _ | DlsBandCallback _
    | SlsKluCallback _ | SlsSuperlumtCallback _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils session =
  if Sundials_configuration.safe then
    match session.ls_callbacks with
    | SpilsCallback1 _ | SpilsCallback2 _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils_band session =
  if Sundials_configuration.safe then
    match session.ls_precfns with
    | BandedPrecFns -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils_bbd session =
  if Sundials_configuration.safe then
    match session.ls_precfns with
    | BBDPrecFns _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

(* Mass solver check functions *)

let mass_check_direct session =
  if Sundials_configuration.safe then
    match session.mass_callbacks with
    | DlsDenseMassCallback _ | DlsBandMassCallback _
    | SlsKluMassCallback _ | SlsSuperlumtMassCallback _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let mass_check_spils session =
  if Sundials_configuration.safe then
    match session.mass_callbacks with
    | SpilsMassCallback _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

(* Types that depend on session *)

type ('k, 'step) serial_session = (Nvector_serial.data, 'k, 'step) session
                                  constraint 'k = [>Nvector_serial.kind]

module SpilsTypes = struct
  include SpilsTypes'

  type ('a, 'k, 's) set_preconditioner =
    ('a, 'k, 's) session -> ('a, 'k) Nvector.t -> unit

  type ('a, 'k, 's) preconditioner =
    LSI.Iterative.preconditioning_type * ('a, 'k, 's) set_preconditioner

  type ('k, 's) serial_preconditioner =
    (Nvector_serial.data, 'k, 's) preconditioner
    constraint 'k = [>Nvector_serial.kind]

end

module MassTypes = struct
  type ('data, 'kind) solver =
    ('data, 'kind, arkstep) session
    -> ('data, 'kind) Nvector.t
    -> unit

  type 'k serial_solver =
    (Nvector_serial.data, 'k) solver
    constraint 'k = [>Nvector_serial.kind]

  module Direct' = struct
    include MassTypes'.Direct'
  end

  module Iterative' = struct
    include MassTypes'.Iterative'

    type ('a, 'k) set_preconditioner =
      ('a, 'k, arkstep) session -> ('a, 'k) Nvector.t -> unit

    type ('a, 'k) preconditioner =
      LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner

    type 'k serial_preconditioner =
      (Nvector_serial.data, 'k) preconditioner
      constraint 'k = [>Nvector_serial.kind]
  end
end

let read_weak_ref x : ('a, 'kind, 'step) session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

let empty_preinnerarray () = Array.init 0 (fun _ -> assert false)

(* Dummy callbacks.  These dummies getting called indicates a fatal
   bug.  Rather than raise an exception (which may or may not get
   propagated properly depending on the context), we immediately abort
   the program. *)
let dummy_rhsfn1 _ _ _ =
  Sundials_impl.crash "Internal error: dummy_rhsfn1 called\n"
let dummy_rhsfn2 _ _ _ =
  Sundials_impl.crash "Internal error: dummy_rhsfn2 called\n"
let dummy_nlsrhsfn _ _ _ =
  Sundials_impl.crash "Internal error: dummy_nlsrhsfn called\n"
let dummy_rootsfn _ _ _ =
  Sundials_impl.crash "Internal error: dummy_rootsfn called\n"
let dummy_errh _ =
  Sundials_impl.crash "Internal error: dummy_errh called\n"
let dummy_errw _ _ =
  Sundials_impl.crash "Internal error: dummy_errw called\n"
let dummy_resw _ _ =
  Sundials_impl.crash "Internal error: dummy_resw called\n"
let dummy_adaptfn _ _ _ =
  (Sundials_impl.crash "Internal error: dummy_adaptfn called\n")
let dummy_stabfn _ _ =
  (Sundials_impl.crash "Internal error: dummy_stabfn called\n")
let dummy_resizefn _ _ =
  Sundials_impl.crash "Internal error: dummy_resizefn called\n"
let dummy_poststepfn _ _ =
  Sundials_impl.crash "Internal error: dummy_poststepfn called\n"
let dummy_poststagefn _ _ =
  Sundials_impl.crash "Internal error: dummy_poststagefn called\n"
let dummy_stagepredictfn _ _ =
  Sundials_impl.crash "Internal error: dummy_stagepredictfn called\n"
let dummy_preinnerfn _ _ =
  Sundials_impl.crash "Internal error: dummy_preinnerfn called\n"
let dummy_postinnerfn _ _ =
  Sundials_impl.crash "Internal error: dummy_postinnerfn called\n"

