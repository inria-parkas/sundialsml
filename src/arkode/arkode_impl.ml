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

(* Types shared between Arkode, Arkode_bbd, Arkode_klu, and Arkode_superlumt. *)

(* This module's purpose is just to define types shared between Arkode,
   Arkode_bbd, Arkode_klu, and Arkode_superlumt. The session and
   linear_solver types must have their implementation exposed to all
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
external crash : string -> unit = "sundials_crash"

type ('data, 'kind) nvector = ('data, 'kind) Nvector.t
module RealArray = Sundials.RealArray

type 'a triple = 'a * 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;
    jac_y   : 'a;
    jac_fy  : 'a;
    jac_tmp : 't
  }

type bandrange = { mupper : int; mlower : int; }

module DlsTypes = struct
  type dense_jac_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg
    -> Dls.DenseMatrix.t
    -> unit

  (* These fields are accessed from arkode_ml.c *)
  type dense_jac_callback =
    {
      jacfn: dense_jac_fn;
      mutable dmat : Dls.DenseMatrix.t option
    }

  let no_dense_callback = {
      jacfn = (fun _ _ -> crash "no dense callback");
      dmat = None;
    }

  type band_jac_fn =
    bandrange
    -> (RealArray.t triple, RealArray.t) jacobian_arg
    -> Dls.BandMatrix.t
    -> unit

  (* These fields are accessed from arkode_ml.c *)
  type band_jac_callback =
    {
      bjacfn: band_jac_fn;
      mutable bmat : Dls.BandMatrix.t option
    }

  let no_band_callback = {
      bjacfn = (fun _ _ _ -> crash "no band callback");
      bmat = None;
    }

  module MassTypes = struct
    type dense_fn =
      float
      -> RealArray.t triple
      -> Dls.DenseMatrix.t
      -> unit

    (* These fields are accessed from arkode_ml.c *)
    type dense_callback =
      {
        massfn: dense_fn;
        mutable dmat : Dls.DenseMatrix.t option
      }

    let no_dense_callback = {
        massfn = (fun _ _ _ -> crash "no dense mass callback");
        dmat = None;
      }

    type band_fn =
      bandrange
      -> float
      -> RealArray.t triple
      -> Dls.BandMatrix.t
      -> unit

    (* These fields are accessed from arkode_ml.c *)
    type band_callback =
      {
        bmassfn: band_fn;
        mutable bmat : Dls.BandMatrix.t option
      }

    let no_band_callback = {
        bmassfn = (fun _ _ _ _ -> crash "no band mass callback");
        bmat = None;
      }
  end
end

module SlsTypes = struct

  type sparse_jac_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg
    -> Sls.SparseMatrix.t
    -> unit

  (* These fields are accessed from arkode_ml.c *)
  type sparse_jac_callback =
    {
      jacfn: sparse_jac_fn;
      mutable smat : Sls_impl.t option
    }

  type sparse_mass_fn =
    float
    -> RealArray.t triple
    -> Sls.SparseMatrix.t
    -> unit

  (* These fields are accessed from arkode_ml.c *)
  type sparse_mass_callback =
    {
      massfn: sparse_mass_fn;
      mutable smat : Sls_impl.t option
    }
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

  type gramschmidt_type = Spils.gramschmidt_type =
    | ModifiedGS
    | ClassicalGS

  type preconditioning_type =
    | PrecNone
    | PrecLeft
    | PrecRight
    | PrecBoth
end

module SpilsTypes' = struct
  include SpilsCommonTypes

  type 'a prec_solve_fn =
    ('a, 'a) jacobian_arg
    -> 'a prec_solve_arg
    -> 'a
    -> unit

  type 'a prec_setup_fn =
    ('a triple, 'a) jacobian_arg
    -> bool
    -> float
    -> bool

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

  module MassTypes' = struct
    type 'd prec_solve_arg =
      {
        rhs   : 'd;
        delta : float;
        left  : bool;
        tmp   : 'd;
      }

    type 'd prec_solve_fn =
         float
      -> 'd prec_solve_arg
      -> 'd
      -> unit

    type 'd prec_setup_fn =
         float
      -> 'd triple
      -> unit

    type 'd times_vec_fn =
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

module AlternateTypes' = struct
  type conv_fail =
    | NoFailures
    | FailBadJ
    | FailOther
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

type arkode_mem
type c_weak_ref

type 'a rhsfn = float -> 'a -> 'a -> unit
type 'a rootsfn = float -> 'a -> Sundials.RealArray.t -> unit
type error_handler = Sundials.error_details -> unit
type 'a error_weight_fun = 'a -> 'a -> unit
type 'a res_weight_fun = 'a -> 'a -> unit

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

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.  *)

type ('a, 'kind) session = {
  arkode     : arkode_mem;
  backref    : c_weak_ref;
  nroots     : int;
  mutable checkvec     : (('a, 'kind) Nvector.t -> unit);
  mutable uses_resv    : bool;

  mutable exn_temp     : exn option;

  mutable problem      : problem_type;
  mutable irhsfn       : 'a rhsfn;
  mutable erhsfn       : 'a rhsfn;

  mutable rootsfn      : 'a rootsfn;
  mutable errh         : error_handler;
  mutable errw         : 'a error_weight_fun;
  mutable resw         : 'a res_weight_fun;

  mutable adaptfn      : 'a adaptivity_fn;
  mutable stabfn       : 'a stability_fn;
  mutable resizefn     : 'a resize_fn;

  mutable linsolver      : ('a, 'kind) linear_solver option;
  mutable ls_callbacks   : ('a, 'kind) linsolv_callbacks;
  mutable ls_precfns     : 'a linsolv_precfns;
  mutable mass_callbacks : ('a, 'kind) mass_callbacks;
  mutable mass_precfns   : 'a mass_precfns;
}

and problem_type =
  | ImplicitOnly
  | ExplicitOnly
  | ImplicitAndExplicit

and ('data, 'kind) linear_solver =
  ('data, 'kind) session
  -> ('data, 'kind) nvector
  -> unit

and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  (* Dls *)
  | DlsDenseCallback of DlsTypes.dense_jac_callback
  | DlsBandCallback  of DlsTypes.band_jac_callback

  (* Sls *)
  | SlsKluCallback of SlsTypes.sparse_jac_callback
  | SlsSuperlumtCallback of SlsTypes.sparse_jac_callback

  (* Spils *)
  | SpilsCallback of 'a SpilsTypes'.jac_times_vec_fn option

  (* Alternate *)
  | AlternateCallback of ('a, 'kind) alternate_linsolv

and 'a linsolv_precfns =
  | NoPrecFns
  | PrecFns of 'a SpilsTypes'.precfns
  | BandedPrecFns
  | BBDPrecFns of 'a ArkodeBbdParamTypes.precfns

and ('data, 'kind) alternate_linsolv =
  {
    linit  : ('data, 'kind) linit' option;
    lsetup : ('data, 'kind) lsetup' option;
    lsolve : ('data, 'kind) lsolve';
  }
and 'data alternate_lsetup_args =
  {
    lsetup_conv_fail : AlternateTypes'.conv_fail;
    lsetup_y : 'data;
    lsetup_rhs : 'data;
    lsetup_tmp : 'data triple;
  }
and 'data alternate_lsolve_args =
  {
    lsolve_ewt : 'data;
    lsolve_y : 'data;
    lsolve_rhs : 'data;
  }
and ('data, 'kind) linit' = ('data, 'kind) session -> unit
and ('data, 'kind) lsetup' =
  ('data, 'kind) session
  -> 'data alternate_lsetup_args
  -> bool
and ('data, 'kind) lsolve' =
  ('data, 'kind) session
  -> 'data alternate_lsolve_args
  -> 'data
  -> unit

and ('a, 'kind) mass_callbacks =
  | NoMassCallbacks

  (* Dls *)
  | DlsDenseMassCallback of DlsTypes.MassTypes.dense_callback
  | DlsBandMassCallback  of DlsTypes.MassTypes.band_callback

  (* Sls *)
  | SlsKluMassCallback of SlsTypes.sparse_mass_callback
  | SlsSuperlumtMassCallback of SlsTypes.sparse_mass_callback

  (* Spils *)
  | SpilsMassCallback of 'a SpilsTypes'.MassTypes'.times_vec_fn

  (* Alternate *)
  | AlternateMassCallback of ('a, 'kind) alternate_mass

and 'a mass_precfns =
  | NoMassPrecFns
  | MassPrecFns of 'a SpilsTypes'.MassTypes'.precfns
  | BandedMassPrecFns
  | BBDMassPrecFns of 'a ArkodeBbdParamTypes.precfns

and ('data, 'kind) alternate_mass =
  {
    minit  : ('data, 'kind) minit' option;
    msetup : ('data, 'kind) msetup' option;
    msolve : ('data, 'kind) msolve';
  }
and ('data, 'kind) minit' = ('data, 'kind) session -> unit
and ('data, 'kind) msetup' =
  ('data, 'kind) session
  -> 'data triple
  -> unit
and ('data, 'kind) msolve' =
  ('data, 'kind) session
  -> 'data
  -> 'data
  -> unit

(* Linear solver check functions *)

let ls_check_dls session =
  if Sundials_config.safe then
    match session.ls_callbacks with
    | DlsDenseCallback _ | DlsBandCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let ls_check_klu session =
  if Sundials_config.safe then
    match session.ls_callbacks with
    | SlsKluCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let ls_check_superlumt session =
  if Sundials_config.safe then
    match session.ls_callbacks with
    | SlsSuperlumtCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let ls_check_spils session =
  if Sundials_config.safe then
    match session.ls_callbacks with
    | SpilsCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let ls_check_spils_band session =
  if Sundials_config.safe then
    match session.ls_precfns with
    | BandedPrecFns -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let ls_check_spils_bbd session =
  if Sundials_config.safe then
    match session.ls_precfns with
    | BBDPrecFns _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

(* Mass solver check functions *)

let mass_check_dls session =
  if Sundials_config.safe then
    match session.mass_callbacks with
    | DlsDenseMassCallback _ | DlsBandMassCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let mass_check_klu session =
  if Sundials_config.safe then
    match session.mass_callbacks with
    | SlsKluMassCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let mass_check_superlumt session =
  if Sundials_config.safe then
    match session.mass_callbacks with
    | SlsSuperlumtMassCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let mass_check_spils session =
  if Sundials_config.safe then
    match session.mass_callbacks with
    | SpilsMassCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let mass_check_spils_band session =
  if Sundials_config.safe then
    match session.mass_precfns with
    | BandedMassPrecFns -> ()
    | _ -> raise Sundials.InvalidLinearSolver

(* Types that depend on session *)

type 'k serial_session = (Nvector_serial.data, 'k) session
                         constraint 'k = [>Nvector_serial.kind]

type 'k serial_linear_solver = (Nvector_serial.data, 'k) linear_solver
                               constraint 'k = [>Nvector_serial.kind]

type ('data, 'kind) mass_solver =
  ('data, 'kind) session
  -> ('data, 'kind) nvector
  -> unit

type 'k serial_mass_solver = (Nvector_serial.data, 'k) mass_solver
                             constraint 'k = [>Nvector_serial.kind]

module SpilsTypes = struct
  include SpilsTypes'

  type ('a, 'k) set_preconditioner =
    ('a, 'k) session -> ('a, 'k) nvector -> unit

  type ('a, 'k) preconditioner =
    | InternalPrecNone  of ('a, 'k) set_preconditioner
    | InternalPrecLeft  of ('a, 'k) set_preconditioner
    | InternalPrecRight of ('a, 'k) set_preconditioner
    | InternalPrecBoth  of ('a, 'k) set_preconditioner

  type 'k serial_preconditioner = (Nvector_serial.data, 'k) preconditioner
                                  constraint 'k = [>Nvector_serial.kind]

  module MassTypes = struct
    include SpilsTypes'.MassTypes'

    type ('a, 'k) set_preconditioner =
      ('a, 'k) session -> ('a, 'k) nvector -> unit

    type ('a, 'k) preconditioner =
      | InternalPrecNone  of ('a, 'k) set_preconditioner
      | InternalPrecLeft  of ('a, 'k) set_preconditioner
      | InternalPrecRight of ('a, 'k) set_preconditioner
      | InternalPrecBoth  of ('a, 'k) set_preconditioner

    type 'k serial_preconditioner = (Nvector_serial.data, 'k) preconditioner
                                    constraint 'k = [>Nvector_serial.kind]
  end
end

module AlternateTypes = struct
  include AlternateTypes'
  type ('data, 'kind) callbacks = ('data, 'kind) alternate_linsolv =
    {
      linit  : ('data, 'kind) linit option;
      lsetup : ('data, 'kind) lsetup option;
      lsolve : ('data, 'kind) lsolve;
    }
  and ('data, 'kind) linit = ('data, 'kind) linit'
  and ('data, 'kind) lsetup = ('data, 'kind) lsetup'
  and ('data, 'kind) lsolve = ('data, 'kind) lsolve'
  and 'data lsetup_args = 'data alternate_lsetup_args = {
    lsetup_conv_fail : conv_fail;
    lsetup_y : 'data;
    lsetup_rhs : 'data;
    lsetup_tmp : 'data triple;
  }
  and 'data lsolve_args = 'data alternate_lsolve_args = {
    lsolve_ewt : 'data;
    lsolve_y : 'data;
    lsolve_rhs : 'data;
  }

  module MassTypes = struct
    type ('data, 'kind) callbacks = ('data, 'kind) alternate_mass =
      {
        minit  : ('data, 'kind) minit option;
        msetup : ('data, 'kind) msetup option;
        msolve : ('data, 'kind) msolve;
      }
    and ('data, 'kind) minit = ('data, 'kind) minit'
    and ('data, 'kind) msetup = ('data, 'kind) msetup'
    and ('data, 'kind) msolve = ('data, 'kind) msolve'
  end
end

let read_weak_ref x : ('a, 'kind) session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

(* Dummy callbacks.  These dummies getting called indicates a fatal
   bug.  Rather than raise an exception (which may or may not get
   propagated properly depending on the context), we immediately abort
   the program. *)
let dummy_erhsfn _ _ _ =
  crash "Internal error: dummy_eresfn called\n"
let dummy_irhsfn _ _ _ =
  crash "Internal error: dummy_iresfn called\n"
let dummy_rootsfn _ _ _ =
  crash "Internal error: dummy_rootsfn called\n"
let dummy_errh _ =
  crash "Internal error: dummy_errh called\n"
let dummy_errw _ _ =
  crash "Internal error: dummy_errw called\n"
let dummy_resw _ _ =
  crash "Internal error: dummy_resw called\n"
let dummy_adaptfn _ _ _ =
  (crash "Internal error: dummy_adaptfn called\n"; 0.0)
let dummy_stabfn _ _ =
  (crash "Internal error: dummy_stabfn called\n"; 0.0)
let dummy_resizefn _ _ =
  crash "Internal error: dummy_resizefn called\n"

